/*	Companion script of
	Golse, F., Hecht, F., Pironneau, O., Tournier, P. H., & Smets, D. (2022). Radiative Transfer For Variable 3D Atmospheres. arXiv preprint arXiv:2208.06410.

	This script also illustrates the use of the Htool Hierarchical Matrix library to compress a dense matrix stemming from a custom user-defined operator defined in a user C++ FreeFEM plugin.

	This can be done by defining and interfacing a Generator class deriving from the virtual class VirtualGenerator of Htool. This is done here in the C++ plugin "RadiativeTransfer_htool.cpp". The Generator object can then be created in the FreeFEM script through the user-defined constructor, for example here taking a mesh3, a meshS, one vector and a new user-defined type KappaGrid (also defined in the plugin) as parameters:

	Generator GenSE(Th3t,ThS,seeface,kappag);
	
	The Generator can then be passed to Htool through the Build() function, which returns the compressed H-Matrix built by Htool from the user-defined copy_submatrix (or get_coef) routine of the Generator computing the matrix coefficients:

	HMatrix HSE = Build(GenSE,VhS,Vht,eta=100,eps=1e-2,minclustersize=10

	To launch:  ff-mpirun -np 8 chamonix.edp -ns -wg -Kmax 9

	In order to run the test cases from the paper, you can download the original mesh at https://freefem.org/misc/chamonix_100K_surf.meshb , change the file in the readmesh command, and run with argument '-fullmesh':
	wget https://freefem.org/misc/chamonix_100K_surf.meshb
	ff-mpirun -np 8 chamonix.edp -ns -wg -Kmax 9 -fullmesh -n 15
*/

load "bem"
load "RadiativeTransfer_htool"
load "qf11to25"
load "medit"
load "iovtk"
load "mmg"
load "shell"

include "getARGV.idp"

verbosity=0;

/********************* SCRIPT PARAMETERS *********************/
int dc = 1, 							// samples every dc points in the nu integrals (speedup with dc>1, danger)
	Newton = 30,						// max Newton iter for the tempe eq
	n = getARGV("-n", 10),				// controls the number of vertices in the mesh
	Niter = 7,  						// outerloop for convergence
	Kmax = getARGV("-Kmax", 9), 		// 1 if kappa constant // in [2,10] otherwise
	gamma = getARGV("-gamma", 0), 		// incidence of sunlight (0 is vertical, 1 is morning and -1 evening)
	heat = 0, 							// if 1, temperature eq is solved 
	nosnow = getARGV("-nosnow", 0), 	// 0 or 1 if no snow
	cloud = getARGV("-cloud", 0), 		// 0 or 1 if cloud
	altkappa = getARGV("-altkappa", 0);	// if 1 (resp-1) Kappa is changed in a range for CO2 (resp CH4)

real sal = 10*nosnow+0.25, beta=0.3, 	// snow altitude limit and snow albedo
	Tsun=1.209, Q0=2.0e-5, sigma=sqr(sqr(pi))/15; // the Sunlight

/********************* OUTPUT FILES *********************/
string basedir("chamonix-output/");
if(mpirank==0) mkdir(basedir);
string tempestr("noon");
if (gamma == 1) tempestr = "morning";
if (gamma == -1) tempestr = "evening";
if (nosnow) tempestr = tempestr + "nosnow";
if (cloud) tempestr = tempestr + "cloud";
if (heat) tempestr = tempestr + "heat";
if (Kmax > 1) tempestr = tempestr + "K";  // res file non grey have a xx.K
if (Kmax > 1 && altkappa>0) tempestr = tempestr + "K";  // res file have a xx.KK
if (Kmax > 1 && altkappa<0) tempestr = tempestr + "KK"; // res file have a xx.KKK

int mem = storageused(); // memory used at this stage
/********************* BUILD THE MESH *********************/
mesh3 Th3;

if (mpirank == 0)
{
	// read the surface mesh of the Chamonix valley:
	meshS Thf = readmeshS("chamonix_surf.meshb");
	//meshS Thf = readmeshS("chamonix_100K_surf.meshb"); // uncomment to use the original fine mesh ; you can download it at https://freefem.org/misc/chamonix_100K_surf.meshb
	// build the volume mesh above the valley:
	real[int] bb(6);
	boundingbox(Thf, bb);
	int nn = 5*n;

	mesh Th2 = movemesh(Thf,transfo=[x,y]); // flatten the mesh, will be the input of buildlayers
	fespace Vhf(Thf, P1);
	fespace Vh2(Th2, P1);
	Vhf zzf = z;
	Vh2 z2;
	z2[] = zzf[]; // get the altitude as a P1 function on the flattened 2D mesh
	
	bool fullmesh = usedARGV("-fullmesh") != -1; // use the original fine adapted surface mesh ; beware, expensive case !
	if (!fullmesh) Th2 = square(nn, nn, [bb[0]+x*(bb[1]-bb[0]),bb[2]+y*(bb[3]-bb[2])]);
	z2 = z2;

	real zmin = z2[].min, zmax = 1;
	plot(z2);

	int[int] l23 = [1,1,2,2,3,3,4,4], lup = [0,6], lbot = [0,5];
	real hx = (bb[1]-bb[0])/nn;
	int nnz = fullmesh ? 5 : (zmax-zmin)/hx;

	func cc = (zmax-z2)/(zmax-zmin);
	Th3 = buildlayers(Th2, coef=cc, nnz, zbound=[z2,zmax],
		labelmid=l23, labelup=lup, labeldown=lbot);
	
	int[int] reqt=[5];
	
	if (fullmesh) Th3 = mmg3d(Th3,requiredTriangle=reqt,hmin=hx,hmax=hx,hgrad=1.5,mem=8192);
}
broadcast(processor(0),Th3);

real[int] bounds(6); // Domain bounds
boundingbox(Th3, bounds);
if (mpirank == 0) cout<< "x,y,z bounds: " << bounds << endl;

int[int] labRT = [5]; // emitting boundary
meshS ThS = extract(Th3, label=labRT, angle=pi); // extract as a surface mesh
if (mpirank == 0) cout << "Number of Vertices = " << Th3.nv << endl;
fespace Vh(Th3, P1);
fespace VhS(ThS, P1);

/********************* READ KAPPA[nu] *********************/
int Nnu = 483; // number of frequencies, will be augmented by 200
real[int] nu(200+Nnu), kappa(200+Nnu);
ifstream kappafile("geminitransmittance.txt");
int j =- 1;
while (j++ < Nnu-1) {
	real auxx, kappaux, nuj;
	kappafile >> nuj >> kappaux >> auxx >> auxx >> auxx;
	kappa[200+j] = 1 - min(1.,max(0.,kappaux)); // necessary
	nu[200+j] = nuj;
	if (altkappa<0 && nuj > 3/4. && nuj < 3/3.) kappa[200+j]=1;
	if (altkappa>0 && nuj > 3/18. && nuj < 3/14.) kappa[200+j]=1;
}
for (int j=0; j<200; j++) { nu[j] = 0.01+(nu[200]-0.01)*j/200.; kappa[j] = kappa[200]; }
Nnu += 200;
real kappa0=0;
for (int j=0; j<Nnu; j++) kappa0 += kappa[j]/Nnu;  // mean kappa is used in Newton steps
// kappa is kappa[nu]*rhof(x,y,z), but integrals are sampled by kappaApprox[K]*rhof(x,y,z)
real[int] kappaApprox(Kmax); // kappa(nu,x) = kappaApprox(K(nu))*rhof(x)  cf Lebesgue integral
for (int K=0; K<Kmax; K++) kappaApprox[K] = 1-real(K)/Kmax; // grey case: kappaApprox[0]=1

int[int] kappaInK(Kmax); // detect if kappa in (kappaApprox[K],kappaApprox[K-1])
kappaInK[0] = (Kmax==1); // used for the grey case
for (int K=1; K<Kmax; K++) {
	kappaInK[K]=0;
	for (int k=0; k<Nnu-dc-1; k+=dc) // every dc points
	if ( (kappa[k]>kappaApprox[K]) && (kappa[k]<=kappaApprox[K-1]) )
		kappaInK[K] = 1;
	if (mpirank == 0)
		cout << K << "K kappaApprox " << kappaApprox[K] << " " << kappaInK[K] << " " << kappaApprox[K-1] <<endl;
}

/********************* DEFINE the x-dependence of kappa, called rho *********************/

func rhof = (0.5-z*0.25) * (1 + (Kmax > 1)) * (1 + cloud * (1+0.5*(z>0.2)*(z<0.8)*(sqr(x-1.5)+sqr(y+1.5)<0.5)) );

Vh rhox = rhof; // rhox is the P1 function equal to rhof
real surf = int3d(Th3)(1.);
if (mpirank == 0) cout << "arithmetic mean of kappa(nu)*rhof = " << kappa0*int3d(Th3)(rhox)/surf << endl << endl;
 
/********************* DEFINE WIND *********************/
Vh u1, u2, u3, dT, dTh;
/*
{
ifstream velocity1("maillageCham/velocity1.txt");
ifstream velocity2("maillageCham/velocity2.txt");
ifstream velocity3("maillageCham/velocity3.txt");
velocity1 >> u1[];
velocity2 >> u2[];
velocity3 >> u3[]; 
if (mpirank==0) plot(u2, wait=1);
}
*/
u2 = 1*z; u1 = 0; u3 = 0;

/********************* Light is x cos(angle with normal)^-. Stored in seeface *********************/
VhS Q0S = Q0*(beta+(1-beta)*(z<sal)); // snow
varf vseeface(u,v) = int2d(ThS,qft=qf1pTlump)(Q0S*(Ns.z+gamma*(Ns.x-0.2*Ns.y))*v/sqrt(1+0.04*gamma*gamma));
real[int] seeface = vseeface(0,VhS);
varf vlump(u,v) = int2d(ThS,qft=qf1pTlump)(v);
real[int] blump = vlump(0,VhS);
seeface ./= blump;
for [i, bi : seeface] bi = -min(bi,0.);

/********************* Interpolation operator from Th3 to ThS *********************/
Vh onGamma;
{
       onGamma[] = -1;
       matrix R = interpolate(VhS,Vh);
       int[int] I(1),J(1);
       real[int] K(1);
       [I,J,K] = R;
       for (int i=0; i<I.n; i++)
               if (K[i] > 0.5) onGamma[][J[i]] = I[i];
}

/********************* Define truncated volume mesh Th3t (without the surface) *********************/
// Th3t will be used to remove surface-surface interactions from the boundary operator,
// as they will be treated separately.

fespace Ph0(Th3,P0);
varf vGamma(u,v) = int3d(Th3)(u*v);
matrix RGamma = vGamma(Vh, Ph0);
Ph0 tGamma;
Vh onG = (onGamma == -1 ? 0 : 1);
tGamma[] = RGamma*onG[];
mesh3 Th3t = trunc(Th3, tGamma == 0);
fespace Vht(Th3t,P1);

matrix Rt = interpolate(Vht,Vh); // restriction matrix from Th3 to Th3t

/********************* BUILD H-Matrices *********************/
verbosity = 1;
HMatrix[int] HVolume(Kmax); // The Volume integral matrices
Vh[int] SE(Kmax); // The surface integral vector operator

Vht SEt;
VhS ones = 1;

for (int K=min(Kmax-1,1); K<Kmax; K++)
if (kappaInK[K]) {
	KappaGrid kappag(bounds,0.01,kappaApprox[K-(K>0)]*rhof);

	// We use the truncated volume mesh Th3t (without the surface) for the boundary operator,
	// as we use an analytic formula for surface-surface interactions
	Generator GenSE(Th3t,ThS,seeface,kappag);
	HMatrix HSE = Build(GenSE,VhS,Vht,eta=100,eps=1e-2,minclustersize=10);
	if (mpirank == 0) cout << HSE.infos << endl; 
	display(HSE);
	SEt[] = HSE*ones[]; // SE = \int_∑ ((x-x').n) exp(-\int_(x,x')kappa(nu,x'))/|x-x'|^3
	SE[K][] = Rt'*SEt[]; // extend the result from truncated Vht to Vh
	// use an analytic formula for surface-surface interactions:
	for [i,bi:SE[K][]] if (onGamma[][i] > -0.5) bi = 0.25*seeface[int(onGamma[][i]+0.5)];

	Generator GenVolume(Th3,kappag);
	HVolume[K] = Build(GenVolume,Vh,Vh,eta=100,eps=1e-2,minclustersize=10);
	if (mpirank == 0) cout << HVolume[K].infos << endl; 
	display(HVolume[K]); // = \int_Omega kappa(nu,x') exp(-\int_(x,x')kappa(nu,x'))/|x-x'|^2
}

/********************* OUTER LOOP *********************/ 
Vh T = 0.01; // initial value of the temperature for iterations
Vh res, sT4=sqr(sqr(T)), J=sT4;

problem heatedp(dT, dTh) = int3d(Th3) (0.02*(dx(dT)*dx(dTh)+dy(dT)*dy(dTh)+dz(dT)*dz(dTh))
		+ (u1*dx(dT)+u2*dy(dT)+u3*dz(dT))*dTh + 15*dT*dTh)
	+ int3d(Th3) (0.02*(dx(T)*dx(dTh)+dy(T)*dy(dTh)+dz(T)*dz(dTh)));

if (Kmax == 1) { /********************* GREY CASE kappaApprox[0] = 1 *********************/
	for (int niter=0; niter<Niter; niter++) {
		sT4 = sigma*sqr(sqr(T));
		res[] = HVolume[0]*sT4[];
		J[] = SE[0][];
		J = J*sigma*sqr(sqr(Tsun));
		J[] += res[];
		T = sqrt(sqrt(abs(J/sigma)));
		if (heat) {
			heatedp;
			T = T + dT;
		}
		sT4 = sT4-J;
		if (mpirank == 0) {
			//plot(J, cmm="J at iter "+ niter, value=1, fill=1);
			plot(sT4, cmm="sT4 at iter "+ niter , value=1, fill=1);
			//plot(T, cmm="T at iter "+ niter, value=1, fill=1);
		}
	}
}
else { /********************* General Nonlinear loop *********************/
	// The method:
	// Let 4 pi J(nu,x) = int_om I dw = \int_∑(x') exp(-kappa(nu,x')|x'-x|)(angle(N,x'-x))^2Bnu(nu,Tsun)d∑ 
	// 							+ \int_O(x') exp(-kappa(nu,x')|x'-x|)/|x'-x|^2 Bnu(nu,T(x'))dx'
	// get T(x) from : \int_0^∞ kappa(nu,x) Bnu(nu,T(x)) dnu = \int_0^∞ kappa(nu,x)J(nu,x)dnu =: J by Newton
	// nu-integrals as double sum: \int_0^∞ f(kappa(nu),nu)=\sum_K=1..Kmax \sum_k f(nu_k,kappa=K/Kmax)dnu_k
	// = \sum_K=1..Kmax H(kappa)*\sum_k F(K/Kmax,nu)dnu_k   when f=H(kappa)*F(kappa,nu)

	// Note: function Bnu is defined in the C++ plugin instead in order to speed up computations. The function corresponds to:
	// func real Bnu(real nu, real T) {return nu*nu*nu/(exp(nu/T)-1);}

	real[int] BnuTsun(Nnu);
	for (int k=0; k<Nnu-1 ;k++) BnuTsun[k] = Bnu(nu[k],Tsun); // aim is to speed up below
	Vh aux, bux, pres, paux = 0;

	real[int] compint(Vh.ndof);
	varf vint(u,v) = int3d(Th3)(v);
	compint = vint(0,Vh);

	// cut expensive loops into MPI slices for parallel computation
	int[int] sizes(mpisize);
	int[int] offsets(mpisize);
	offsets = 0;
	int chunksize = Vh.ndof%mpisize == 0 ? Vh.ndof/mpisize : Vh.ndof/mpisize + 1;

	for (int i = 0; i < mpisize; i++)
	if (i != mpisize - 1) {
		sizes[i] = chunksize;
		offsets[i] = i*chunksize;
	}
	else {
		sizes[i] = Vh.ndof - i*chunksize;
		offsets[i] = i*chunksize;
	}

	int ii = offsets[mpirank];
	real[int] chunkbuff(sizes[mpirank]), chunkbuff2(sizes[mpirank]);

	for (int niter=0; niter<Niter; niter++){ // Nonlinear loop
		bux=0; J=0;
		for (int K=1; K<Kmax; K++)
		if (kappaInK[K]) {
			real auxs = 0;
			chunkbuff2 = 0;       
			for (int k=0; k<Nnu-dc-1; k+=dc)  // every dc points, can be changed
			if ( (kappa[k]>kappaApprox[K]) && (kappa[k]<=kappaApprox[K-1]) ) {    
				real dnu = nu[k+dc] - nu[k];
				for [i,bi:chunkbuff] bi = rhox[](i+ii)*Bnu(nu[k],T[][i+ii]);
				chunkbuff2 += dnu*kappaApprox[K-1]*chunkbuff; // aux(x) = int_{kappa[k]=K} kappa Bnu(nu,T(x)) dnu
				auxs += BnuTsun[k]*dnu; // auxs = Q0*int_{kappa(nu)=K} Bnu(nu,Tsun) dnu
			}
			mpiAllgatherv(chunkbuff2, aux[], mpiCommWorld, sizes, offsets);
			res[] = HVolume[K]*aux[];
			for [i,bi:chunkbuff] bi = res[][i+ii] + SE[K][][i+ii]*rhox[][i+ii]*kappaApprox[K-1]*auxs;
			mpiAllgatherv(chunkbuff, res[], mpiCommWorld, sizes, offsets);
			bux[] += aux[]; // bux = \sum_K kappa_K \int_{kappa(nu)=K} Bnu(nu[k],T(x))dnu 
			J[] += res[]; // J(x) = int_{kappa(nu)=K} kappa exp(-kappa|x'-x|)/|x'-x|^2 Bnu(nu,T(x)) dnu
		}
		bux = bux - J;
		if (mpirank == 0) {
			//plot(J, cmm="J at n="+ niter ,value=1,fill=1);
			plot(bux, cmm="Bnu-J at n = "+ niter ,value=1, fill=1); // J = bux => converged
		}

		int newton = 0; // Newton iterations for T

		for [i,bi:chunkbuff] bi = sqrt(sqrt(abs(J[][i+ii]/sigma/rhox[][i+ii]/kappa0))); // guess to start Newton
		mpiAllgatherv(chunkbuff, T[], mpiCommWorld, sizes, offsets);

		pres = abs(bux/J-1);
		while (pres[]'*compint > 0.001*surf) {
			if (newton>0) { // compute bux with new T
				chunkbuff2 = 0; 
				for (int k=0; k<Nnu-dc-1; k+=dc) { 
					real dnu = nu[k+dc]-nu[k];
					for [i,bi:chunkbuff] bi = rhox[](i+ii)*kappa[k]*Bnu(nu[k],T[][i+ii]);
					chunkbuff2 += dnu*chunkbuff;
				}
				mpiAllgatherv(chunkbuff2, bux[], mpiCommWorld, sizes, offsets);
			}
			pres = bux/J-1;
			if (mpirank == 0) {
				if (Kmax == 1) {
					paux = rhox*kappaApprox[0]*sigma*T^4; 
					cout << niter << " <- niter, error ->" << int3d(Th3)(bux/paux-1)/surf << endl;
				} 
				else 
					cout << niter << " <- niter, newton ->" << newton << " error = " << pres[]'*compint <<endl;
			}
			pres = abs(bux/J-1);
			// bux=paux; T=sqrt(sqrt(J/(sigma*rhox*kappaApprox[0]))); // error check in the grey case
			// T = T - 0.5*(bux-J)/(4*sigma*rhox*kappa0*T*T*T);
			// Newton relaxation parameter 0.5
			for [i,bi:chunkbuff] bi = T[][i+ii] - 0.5*(bux[][i+ii]-J[][i+ii])/(4*sigma*rhox[][i+ii]*kappa0*T[][i+ii]*T[][i+ii]*T[][i+ii]);
			mpiAllgatherv(chunkbuff, T[], mpiCommWorld, sizes, offsets);
			if (newton++ > Newton) bux=J;
		}

		if (heat) {
			heatedp; T=T + dT;
		}

		plot(T, cmm="T at iter = "+ niter, value=1, fill=1);
	}
	// End of non linear loop

	/********************* Display J_nu(x,y,z) *********************/
	if (Kmax == 9) {
		string lightstr(basedir+tempestr+"light.txt");
		ofstream myfile2(lightstr);
		for (int k=0; k<Nnu-dc-1; k+=dc) { // every dc points, can be changed
			int K = -1; 
			for (int K1=1; K1<Kmax; K1++)
			if ( (kappa[k]>kappaApprox[K1]) && (kappa[k]<=kappaApprox[K1-1]) ) K = K1;
			if (K > 0) {
				for [i,bi:chunkbuff] bi = Bnu(nu[k],T[][i+ii]);
				mpiAllgatherv(chunkbuff, res[], mpiCommWorld, sizes, offsets);
				bux[] = HVolume[K]*res[];
				for [i,bi:chunkbuff] bi = bux[][i+ii] + BnuTsun[k]*SE[K][][i+ii];
				mpiAllgatherv(chunkbuff, J[], mpiCommWorld, sizes, offsets);
				if (mpirank == 0) myfile2 << 3/nu[k] << " " << 1e5*J(0.3,1.5,-1.5) << " " << kappa[k]
					<< " " << int(kappa[k]*10+0.5)/10. << " " << kappaApprox[K-1] << " " << K-1 << " " << k << endl;
			}
		}
	}
}
int tmem = storagetotal(); // max memory used at this stage


/********************* GRAPHICS *********************/
real hx=bounds[1]-bounds[0], hy=bounds[3]-bounds[2], hz=bounds[5]-bounds[4]; // domain size
mesh Tp = square(10*n, 5*n, [bounds[2]+hx*x, bounds[4]+hz*y]);
fespace Vp(Tp,P1);
Vp Temp = 4798*T(1.5,x,y)-273;
plot(Temp, value=1, fill=1);

if (mpirank == 0) {
	ofstream myfile(basedir+tempestr+"tempe.txt");
	real[int] XX(2000), YY(2000);
	int kk = 0;
	for (real Z=bounds[4]; Z<bounds[5];Z+=0.01) {
		real Tx = 4798*T(1.5,-1.5,Z)-273;
		XX[kk] = Z; YY[kk++] = Tx;
		cout << Z << "  " << Tx << endl;
		myfile << Z << "  " << Tx << endl;
	}
	XX.resize(kk); YY.resize(kk);
//	plot([XX,YY], cmm = "T in ["+YY.min+","+YY.max+"]");
}
VhS Ts = 4798*T-273, dTs=4798*dT;
//	plot(u1); plot(u2); plot(u3);
meshS ThS2 = movemesh(ThS,[x,y,2*z]); // to enhance graphics
fespace VhS2(ThS2,P1);
VhS2 Ts2;
Ts2[]= dTs[];
Ts2[] = Ts[];
plot(Ts2,fill=1,value=1);
//if (mpirank == 0) medit("Temperature",ThS2,Ts2);
//	T= 4798*T-273;
//	if (mpirank == 0) medit("Temperature",Th3,T);
int[int] fforder = [1];
if (mpirank == 0) savevtk(basedir+tempestr+".vtk",ThS2,Ts2,order=fforder);
if (mpirank == 0) cout << "file name: " << basedir+tempestr << endl;
if (mpirank == 0) cout << "Memory used (Go) = " << (((tmem-mem)/1000.)/1000.)/1000. << endl;
