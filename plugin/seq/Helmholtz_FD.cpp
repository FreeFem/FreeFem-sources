/******************************************************************************/
/*  This plugin implements the 3D 27-point finite-difference stencil with     */
/*  optimized weights for the discretization of the Helmholtz equation from   */
/*                                                                            */
/*  Operto, S., Virieux, J., Amestoy, P., L’Excellent, J. Y., Giraud, L.,     */
/*  & Ali, H. B. H. (2007). 3D finite-difference frequency-domain modeling of */
/*  visco-acoustic wave propagation using a massively parallel direct solver: */
/*  A feasibility study. Geophysics, 72(5), SM195-SM211.                      */
/******************************************************************************/

#include "ff++.hpp"

class HelmholtzFD_Op : public E_F0mps {
 public:
  Expression expTh;
  Expression expomega, expmu;
  static const int n_name_param = 3;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  HelmholtzFD_Op(const basicAC_F0 &args, Expression Th) : expTh(Th) {
    args.SetNameParam(n_name_param, name_param, nargs);
    expTh = to< const Mesh3 * >(args[0]);
    expomega = to< Complex >(args[1]);
    expmu = to< double >(args[2]);
  }

  long arg(int i, Stack stack, long a) const { return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;}
  string* arg(int i, Stack stack, string* a) const { return nargs[i] ? GetAny< string* >((*nargs[i])(stack)) : a;}
  KN_<long>  arg(int i,Stack stack, KN_<long> a ) const { return nargs[i] ? GetAny<KN_<long> >((*nargs[i])(stack)): a;}

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type HelmholtzFD_Op::name_param[] = {
  {"npml", &typeid(long)},
  {"pmlsides", &typeid(KN_<long>)},
  {"scheme", &typeid(string*)}
};

class  HelmholtzFD : public OneOperator {
  public:
    HelmholtzFD() : OneOperator(atype< newpMatrice_Creuse< Complex > >( ), atype< const Mesh3 * >( ), atype< Complex >( ), atype< double >( )) {};
    E_F0 *code(const basicAC_F0 &args) const { return new HelmholtzFD_Op(args, t[0]->CastTo(args[0])); }
};

long kkindex(long i1, long i2, long i3, long n1, long n2, long n3) {
  return i2*n3*n1+i1*n3+i3;
}

void subdamp(long n,int npml,int apml, double h, double omega, std::vector<std::complex<double>>& damp, std::vector<std::complex<double>>& dampb, int sides) {
  double half_pi = 1.570796327;
  std::complex<double> ci(0,1);
  
  for (int i=0; i<=n+1; i++) {
    damp[i] = 1.;
    dampb[i] = 1.;
  }

  double xpml = npml*h;
  double xmax = (n-1)*h;
  double x, xb, eps, epsb;

  for (int i=1; i<=npml; i++) {
    x = (i-1)*h;
    xb = (i-1)*h+0.5*h;
    eps  = apml*(1.-cos((xpml-x)*half_pi/xpml));
    epsb = apml*(1.-cos((xpml-xb)*half_pi/xpml));
    if (sides & 1) {
      damp[i] = 1./(1.+ci*eps/omega);
      dampb[i] = 1./(1.+ci*epsb/omega);
    }
    if (sides & 2)
      damp[n-i+1] = 1./(1.+ci*eps/omega); 
  }

  damp[0] = damp[1];
  damp[n+1] = damp[n];
  
  for (int i=1; i<=npml+1; i++) {
    xb = xmax+0.5*h-(i-1)*h;
    epsb = apml*(1.-cos((xb-(xmax-xpml))*half_pi/xpml));
    if (sides & 2)
      dampb[n-i+1] = 1./(1.+ci*epsb/omega);
  }

  dampb[0] = dampb[1];
  dampb[n+1] = dampb[n];
}

template<class FESpaceT>
MatriceMorse<double> *  buildInterpolationMatrixT1(const FESpaceT & Uh,const KN_<double> & xx,const KN_<double> & yy ,const KN_<double> & zz)
{
  typedef typename FESpaceT::Mesh MeshT;
  typedef typename FESpaceT::FElement FElementT;
  typedef typename MeshT::Element ElementT;
  typedef typename FESpaceT::Rd RdT;
  typedef typename ElementT::RdHat RdHatT;

  int op=op_id; //  value of the function
  int icomp=0;
  bool inside=false;

  int n=Uh.NbOfDF;
  int mm=xx.N();
  int nbxx= mm;
  const MeshT & ThU = Uh.Th; // line
  FElementT Uh0 = Uh[0];
  int nbdfUK= Uh0.NbDoF();
  int NUh= Uh0.N;

  const int sfb1=Uh0.N*last_operatortype*Uh0.NbDoF();
  KN<double> kv(sfb1);
  R * v = kv;
  const R eps = 1e-10;

  What_d whatd= 1 << op;
  MatriceMorse<double> * m = new MatriceMorse<double>(n,mm,0,0);
  RdHatT Phat;
  bool outside;

  for(int ii=0;ii<nbxx;ii++){
    if(verbosity>9) cout << " Find ThU " <<ii << ":" <<  RdT(xx[ii],yy[ii],zz[ii]) << endl;
    const ElementT *ts=ThU.Find(RdT(xx[ii],yy[ii],zz[ii]),Phat,outside);
    if(outside && !inside) continue;
    int it = ThU(ts);
    FElementT KU(Uh[it]);
    KNMK_<double> fb(v,nbdfUK,NUh,last_operatortype);
    Uh0.tfe->FB(whatd,ThU,ThU[it],Phat,fb);
    KN_<double> Fwi(fb('.',icomp,op));
    for (int idfu=0;idfu<nbdfUK;idfu++){
      int  j = ii;
      int  i = KU(idfu);
      R c = Fwi(idfu);
      if(Abs(c)>eps)
        (*m)(i,j) += c;
    }
  }
  return m;
}

template<class T> bool cmp(const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; }

AnyType HelmholtzFD_Op::operator( )(Stack stack) const {

  typedef typename Mesh3::Element Element;
  typedef typename Mesh3::RdHat RdHat;
  static const int nvedgeTet[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  static const int nvedgeTria[3][2] = {{1, 2}, {2, 0}, {0, 1}};
  static const int nvedgeSeg[1][2] = {{0, 1}};

  const int d = RdHat::d;

  MatriceMorse< Complex > *amorse = 0;
  typedef const Mesh3 *pmesh3;
  const Mesh3 *pTh = GetAny< pmesh3 >((*expTh)(stack));
  ffassert(pTh);
  const Mesh3 &Th(*pTh);

  double lx = 1e+30, ly=1e+30, lz=1e+30, ux=-1e+30, uy=-1e+30, uz=-1e+30, h = 1e+30;

  for (int i=0; i<Th.nbe; i++) {
    for (int j=0; j<3; j++) {
      const R3& p = Th.be(i)[j];
      lx = min(p.x,lx); ly = min(p.y,ly); lz = min(p.z,lz);
      ux = max(p.x,ux); uy = max(p.y,uy); uz = max(p.z,uz);
    }
    R3 p1 = Th.be(i)[1] - Th.be(i)[0];
    R3 p2 = Th.be(i)[2] - Th.be(i)[0];
    h = min(h,min(p1.norme(),p2.norme()));
  }

  double dx = ux-lx, dy = uy-ly, dz = uz-lz;

  //std::cout << "h = " << h << std::endl;

  long n1 = dx/h + 1, n2 = dy/h + 1, n3 = dz/h + 1;

  std::complex<double> omega = GetAny<std::complex<double>>((*expomega)(stack));

  //std::cout << n1 << " " << n2 << " " << n3 << " " << Th.nv << std::endl;

  KN<double> mu(n1*n2*n3);
  KN<double> xx(n1*n2*n3), yy(n1*n2*n3), zz(n1*n2*n3);
  KN<long> boundarynodesFD(n1*n2*n3,0L);

  int cpt = 0;
  for (long i2 = 0; i2 < n2; i2++)
  for (long i1 = 0; i1 < n1; i1++)
  for (long i3 = 0; i3 < n3; i3++) {
    xx[cpt] = lx + i1*h; yy[cpt] = ly + i2*h; zz[cpt] = lz + i3*h;
    MeshPoint* mp(Fem2D::MeshPointStack(stack));
    mp->set(xx[cpt], yy[cpt], zz[cpt]);
    mu[cpt] = GetAny<double>( (*expmu)(stack) );
    if ((i1 == 0) || (i1 == n1-1) || (i2 == 0) || (i2 == n2-1) || (i3 == 0) || (i3 == n3-1))
      boundarynodesFD[kkindex(i1,i2,i3,n1,n2,n3)] = 1;
    cpt++;
  }

  {
    std::vector<std::vector<std::complex<double>>> dd(3), db(3);
    
    dd[0].resize(n1+2); dd[1].resize(n2+2); dd[2].resize(n3+2);
    db[0].resize(n1+2); db[1].resize(n2+2); db[2].resize(n3+2);
       
    int npml(arg(0, stack, 8));
 
    KN<long> def(6,1);
    KN<long> pmlsides(arg(1,stack,def));

    subdamp(n1,npml,90,h,omega.real(),dd[0],db[0],pmlsides[0]+2*pmlsides[1]);
    subdamp(n2,npml,90,h,omega.real(),dd[1],db[1],pmlsides[2]+2*pmlsides[3]);
    subdamp(n3,npml,90,h,omega.real(),dd[2],db[2],pmlsides[4]+2*pmlsides[5]);

    MatriceMorse< Complex > *pAij = new MatriceMorse< Complex >(n1*n2*n3, n1*n2*n3, 0, 0), &Aij = *pAij;
    
    std::complex<double> omega2 = omega*omega;

    double c, d, e, f, w1 ,w2, w3; // weights of the mixed grid stencil
    
    string stringdef = string("G4810");
    string *scheme(arg(2,stack,&stringdef));

    if (*scheme == "G4810") {
      c=0.4966390;
      d=7.5123318E-02;
      e=4.3846383E-03;
      f=6.7614019E-07;
      w1=5.0247996E-05;
      w2=0.8900359;
      w3=0.1099138;
    }
    else if (*scheme == "G4") {
      c=0.5915900;
      d=4.9653493E-02;
      e=5.1085097E-03;
      f=6.1483691E-03;
      w1=8.8075437E-02;
      w2=0.8266806;
      w3=8.5243940E-02;
    }
    else {
      cerr << "ERROR HelmholtzFD: available schemes are \"G4\" and \"G4810\"" << endl;
      ffassert(0);
    }

    double w2u=w2/3.;
    double w3u=w3/4.;
    double h2 = 1./(h*h);

    long k, l, l1, l2, l3, ll[3], li, lj, lk, ii, ij, ik;
    long ind6[6][3] = {{1,0,0},{0,1,0},{0,0,1},
                      {-1,0,0},{0,-1,0},{0,0,-1}};
    long ind12[12][3] = {{1,1,0},{0,1,1},{1,0,1},
                        {-1,-1,0},{0,-1,-1},{-1,0,-1},
                        {-1,1,0},{1,-1,0},
                        {0,-1,1},{0,1,-1},
                        {1,0,-1},{-1,0,1}};
    long ik12[12] = {2,0,1,2,0,1,2,2,0,0,1,1};
    long ind8[8][3] = {{1,1,1},{-1,-1,-1},
                      {-1,1,1},{1,-1,-1},
                      {1,-1,1},{-1,1,-1},
                      {1,1,-1},{-1,-1,1}};

    for (long i3 = 0; i3 < n3; i3++)
    for (long i2 = 0; i2 < n2; i2++)
    for (long i1 = 0; i1 < n1; i1++) {    
      /* Node 000 */
      l1 = i1;
      l2 = i2;
      l3 = i3;
      k = kkindex(i1,i2,i3,n1,n2,n3);
      l = k;

      Aij(k, l) = c*omega2/mu[k]
      -w1*h2*(                                         				
               dd[1][l2+1]*(db[1][l2+1]+db[1][l2])               				
              +dd[2][l3+1]*(db[2][l3+1]+db[2][l3])               			
              +dd[0][l1+1]*(db[0][l1+1]+db[0][l1]) // (6)
             )  //R1
      -w2u*h2*2*(
               dd[0][l1+1]*(db[0][l1+1]+db[0][l1])
              +dd[1][l2+1]*(db[1][l2+1]+db[1][l2])
              +dd[2][l3+1]*(db[2][l3+1]+db[2][l3]) // (0.75*8+6) 
                )//R2-R4
      -w3u*h2*0.5*4*(                                            				
               dd[1][l2+1]*(db[1][l2+1]+db[1][l2])               				
              +dd[2][l3+1]*(db[2][l3+1]+db[2][l3])               			
              +dd[0][l1+1]*(db[0][l1+1]+db[0][l1]) // (0.5*24) 
                    ); //B1-B4

      /* 6 nodes */
      for (int q=0; q < 6; q++) {
        l1 = i1+ind6[q][0];
        l2 = i2+ind6[q][1];
        l3 = i3+ind6[q][2];

        ii = q%3;
        ij = (ii+1)%3;
        ik = (ij+1)%3;

        ll[0] = i1; ll[1] = i2; ll[2] = i3;
        li = ll[ii]; lj = ll[ij]; lk = ll[ik];

        if (l1 >= 0 && l1 < n1 && l2 >= 0 && l2 < n2 && l3 >= 0 && l3 < n3) {
          l = kkindex(l1,l2,l3,n1,n2,n3);
          Aij(k, l) = d*omega2/mu[l]

          +w1*h2*dd[ii][li+1]*db[ii][li+1*(q<3)] //R1

          +w2u*0.25*h2*(
             dd[ii][li+1]*db[ii][li+1*(q<3)]*4.
            -dd[ij][lj+1]*(db[ij][lj+1]+db[ij][lj])
            -dd[ik][lk+1]*(db[ik][lk+1]+db[ik][lk])
          )
          +w2u*h2*dd[ii][li+1]*db[ii][li+1*(q<3)] //R2-R4

          +w3u*h2*0.5*4*dd[ii][li+1]*db[ii][li+1*(q<3)]; //B2-B4
        }
      }

      /* 12 nodes */
      for (int q=0; q < 12; q++) {
        l1 = i1+ind12[q][0];
        l2 = i2+ind12[q][1];
        l3 = i3+ind12[q][2];
        
        ik = ik12[q];
        ii = (ik+1)%3;
        ij = (ii+1)%3;
        
        ll[0] = i1; ll[1] = i2; ll[2] = i3;
        li = ll[ii]; lj = ll[ij]; lk = ll[ik];
        
        long si = ind12[q][ii] == 1;
        long sj = ind12[q][ij] == 1;
        
        if (l1 >= 0 && l1 < n1 && l2 >= 0 && l2 < n2 && l3 >= 0 && l3 < n3) {
          l = kkindex(l1,l2,l3,n1,n2,n3);
          Aij(k, l) = e*omega2/mu[l]
            +w2u*h2*0.25*(dd[ii][li+1]*db[ii][li+si]+dd[ij][lj+1]*db[ij][lj+sj]) //R2-R4 0.25*2
            -w3u*h2*0.5*dd[ik][lk+1]*(db[ik][lk]+db[ik][lk+1]); //B1-B4 0.5*2
        }
      }

      /* 8 nodes */
      for (int q=0; q < 8; q++) {
        l1 = i1+ind8[q][0];
        l2 = i2+ind8[q][1];
        l3 = i3+ind8[q][2];
        if (l1 >= 0 && l1 < n1 && l2 >= 0 && l2 < n2 && l3 >= 0 && l3 < n3) {
          l = kkindex(l1,l2,l3,n1,n2,n3);
          Aij(k, l) = f*omega2/mu[l]
          +w3u*h2*0.5*(dd[0][i1+1]*db[0][i1+(ind8[q][0] > 0)]
                      +dd[1][i2+1]*db[1][i2+(ind8[q][1] > 0)]
                      +dd[2][i3+1]*db[2][i3+(ind8[q][2] > 0)]); //B1-B4
        }
      }    
    }
    amorse = pAij;
  }

  FESpace3 Uh(Th);

  MatriceMorse<double>* R = buildInterpolationMatrixT1<FESpace3>(Uh, xx, yy, zz);

  //std::cout << n1*n2*n3 << " " << Th.nv << " " << R->n << " " << R->m << " " << R->nnz << " " << amorse->n << " " << amorse->m <<  std::endl;

  std::vector<char> boundarynodesTh(Th.nv,0);

  R->COO();
  amorse->CSR();
  unsigned int m = R->nnz;
  KN<int> lg(m+1,0);
  std::vector<signed int> tmpVec;
  tmpVec.resize(amorse->m);
  for(long i = 0; i < m; ++i)
    tmpVec[R->j[i]] = i + 1;

  std::vector<std::pair<int, std::complex<double>> > tmp;
  tmp.reserve(amorse->nnz);

  lg[0] = 0;
  for(long i = 0; i < m; ++i) {
    for(long j = amorse->p[R->j[i]]; j < amorse->p[R->j[i] + 1]; ++j) {
      long col = tmpVec[amorse->j[j]];
      if(col != 0)
        tmp.push_back(std::make_pair(col - 1, amorse->aij[j]));
    }
    std::sort(tmp.begin() + lg[i], tmp.end(),cmp<std::complex<double>>);
    lg[i + 1] = tmp.size();

    boundarynodesTh[i] = boundarynodesFD[R->j[i]];
  }

  delete R;

  amorse->clear();
  amorse->resize(m,m);
  MatriceMorse<std::complex<double>> &MA = *amorse;
  MA.half = 0;
  for(int i=0; i<m; ++i)
  for(int k= lg[i]; k < lg[i+1]; ++k) {
    int j= tmp[k].first;
    std::complex<double> aij = tmp[k].second;
    MA(i,j) = aij;
  }

  MA.SetBC(&boundarynodesTh[0], -1);

  return SetAny<newpMatrice_Creuse<std::complex<double>>>(newpMatrice_Creuse<std::complex<double>>(stack, &MA));
}

static void Load_Init( ) {
  Global.Add("HelmholtzFD", "(", new HelmholtzFD);
}

LOADFUNC(Load_Init)
