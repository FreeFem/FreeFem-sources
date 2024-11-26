int nn=10;;
mesh Th2=square(nn,nn,region=0);
fespace Vh2(Th2,P2);
// // label  face  numbering 
//      1 :  ( x == xmin)        2 :  ( x == xmax) 
//      3 :  ( y == ymin)        4 :  ( y == ymax) 
//      5 :  ( z == zmin)        6 :  ( z == zmax) 
// ---
int[int] rup=[0,5],  rdown=[0,6], rmid=[4,1,2,2, 1,3 ,3,4];
real zmin=0,zmax=1;

mesh3 Th=buildlayers(Th2,nn,
  zbound=[zmin,zmax],
  // region=r1, 
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);
cout << "Th :  nv = " << Th.nv << " nt =" << Th.nt << endl;
savemesh(Th,"Th.mesh");
verbosity=10;

//  the Finite element space with full periodic condition in 3 axes
fespace Vh(Th,P2,periodic=[[1,y,z],[2,y,z],[3,x,z],[4,x,z],[5,x,y],[6,x,y]]);
verbosity=2;

// a  code to build some verification ....
fespace Vhh(Th,P2);

int[int] num(Vhh.ndof);
num=-1;
int er=0;
for(int k=0;k<Th.nt;++k)
  {
    int err=0;
    for(int i=0;i<4;i++) 
      {
	if(num[Vhh(k,i)]== -1)
	  num[Vhh(k,i)] = Vh(k,i);
	else if(num[Vhh(k,i)] != Vh(k,i))
	  {
	    ++err;
	    cout << " bug " << k <<  " : " << num[Vh(k,i)]  << " !=  " << Vhh(k,i) << endl; 
	  }
      }
    if(err)
      {
	for(int i=0;i<4;i++) cout << Vh(k,i) << " ";     cout << endl;
	for(int i=0;i<4;i++) cout << Vhh(k,i) << " ";    cout << endl << endl;;
      }
    er+=err;
  }

// ++++++
int  n1 = nn+nn+1; //   P2 =>  
int  n2 = n1-1; //
int  nnn=n2*n2*n2;
int nnn1=n1*n1*n1;
cout << " ndf pare= " << Vh.ndof << " " << nnn << endl;
cout << " ndf  = " << Vhh.ndof << " " << nnn1 << endl;
assert(er==0 && nnn == Vh.ndof && nnn1 == Vhh.ndof); // some verification ...
  Vh u,v,uu;
  real x0=2*pi/3,y0=2*pi/4,z0=2*pi*2/3;
  func ue= sin(2*pi*x+x0)*sin(2*pi*y+y0)*sin(2*pi*z+z0);
  real cc= -3*(2*pi)^2 ;
  func f = -cc*ue;
  uu=ue;
macro Grad(u) [dx(u),dy(u),dz(u)] //;
  solve P(u,v,solver=CG)= int3d(Th)(Grad(u)'*Grad(v)) - int3d(Th)(f*v); //') ;
cout << "Err L2 = " << sqrt(int3d(Th)( square(u-uu)) ) << endl;

// FFCS: add 3D view

///Vh2 u0=u(x,y,0);
///Vh2 u1=u(x,y,1);
///plot(u0,u1,wait=1);

plot(u,nbiso=10);
