typedef int integer;
class BijanMO  
  {
  public: 
    typedef double R;
    typedef KN<R> Vect;
    typedef KNM<R> Mat;
    int debug; //  niveau de print 
    bool diagrand;
    int ndim;
    int n100;
    int nbsol;
    Vect cstr, cstropt;
    double finit,fseul, fseulopt;
    integer ncstr;
    double costsaveming, gnormsave, costsavemin;
    integer nbeval, nbevalp;
    Vect  feval,xoptg,xopt1;
    Mat xfeval;
    KN<double> xmin,xmax;
    // data file 
    int nbrestart,nbext1,nbbvp,nbgrad,ifd;
    double epsfd,rho000,epsloc,epsij;
      int typealgo; //   1 CG 
      string * histpath;// 0 => no file
      string * histcpath;// 0 => no file
      BijanMO(int nndim,
	    int wnbrestart=1,
	    int wnbext1=1,
	    int wnbbvp=5,
	    int wnbgrad=5,
	    double wepsfd=1e-5,
	    double wrho000=100,
	    double wepsloc=1e-4,
	    double wepsij=1e-6,
	    int nn100=100)
      : debug(1),
	diagrand(1), //  choix of diag rand vector or not. 
	ndim(nndim),
	n100(nn100),
	nbsol(1000),
	cstr(n100),cstropt(n100),
	feval(nbsol),
	xoptg(ndim),xopt1(ndim),xfeval(nbsol,ndim),
	xmin(ndim),xmax(ndim),
	nbrestart(wnbrestart),nbext1(wnbext1),
	nbbvp(wnbbvp),nbgrad(wnbgrad),//ifd(wifd),
	epsfd(wepsfd),rho000(wrho000),epsloc(wepsloc),epsij(wepsij),
      typealgo(1),
      histpath(0),
      histcpath(0)
    {
      cout << nbrestart << " == " << wnbrestart << endl;
      
      /*   1 1 1 10 rho =1 epsloc=1e4 (unused), epsglo=1, */
      /*      open(1,file='DATA_BMO.txt',status='unknown') */
      /*      read(1,*) nbrestart */
      /*      read(1,*) nbext1 */
      /*      read(1,*) nbbvp */
      /*      read(1,*) nbgrad */
      /*      read(1,*) ifd,epsfd */
      /*      read(1,*) rho000,epsloc */
      /*      read(1,*) epsij */
      /*      read(1,*) ndim */
      /*      close(1) */  
      /*<       call init(ndim,vinit,xmin,xmax)>*/


}
    

    double main(Vect & vinit,Vect & xxmim,Vect & xxmax);
    double  funcapp(Vect & x, Vect &fpx);
    void    funcp(Vect &x, Vect &fpx,double f);
    double  fun(Vect & x,  Vect& temp, Vect& g, double ro);
    double  ropt_dicho(Vect x, Vect temp, double & ro, Vect g, double ccout);
    int     gradopt(Vect & x1, Vect & fpx, 
		      Vect & temp, double &rho, double  &f,
		      double &gnorm,  
		      Vect & fpx0, Vect & hgc);
    void     tir( Vect &v, Vect &fpx);
    void     rand( Vect &v);
      
    double func(Vect & x) 
    {
      double f=J(x);
      if(nbeval>=0) 
	{
	  int ieval=nbeval++%nbsol;
	  xfeval(ieval,':')=x;
	  feval(ieval)=f;
	}
      return f;
    }

    virtual ~BijanMO(){}
    // les 4 functions utilisateur 
   // virtual void init(Vect & xinit) = 0;    
    virtual double J( Vect & x)=0;
    virtual R * DJ(Vect & x, Vect &fpx){return 0;} // do not existe 
    virtual void result(Vect & xoptg,Vect &vinit){} 
 
};
