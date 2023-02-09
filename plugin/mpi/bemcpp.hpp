struct Data_Funct {
    KN_<double> inti;
    pmesh3 pTh;
    double kappa0,eps;
    Data_Funct(KN_<double> ii, pmesh3 pth3,double k,double ee=0.001) : inti(ii),pTh(pth3),kappa0(k),eps(ee) {}
};
Data_Funct **init_datafunc(Data_Funct ** const &  ppf,KN<double>* const &  pinti,pmesh3  const &  pTh,double const &  kappa0){
    *ppf = new Data_Funct(*pinti,pTh,kappa0);
    return ppf;
}
double RTFunc(Data_Funct ** const &  ppf,long const & i,long const & j)
{
    using Fem2D::R3;
    const double fourpi = Pi*4.;
    const Data_Funct & df = **ppf;
    const Mesh3 & Th3f= *df.pTh;
    const double * inti= df.inti;
    R3 O = Th3f(i);
    R3 PP = Th3f(j);
    R3 OP = R3(PP,O);
    double lOP2 = OP.norme2();
    double lOP = sqrt(lOP2);
    double sc = inti[j]*exp(-(lOP)*df.kappa0)/(lOP2+df.eps);
    return sc * df.kappa0/fourpi;

}

long AdrRTFunc()
{ return (long) & RTFunc;}

class BemKernel;
class BemPotential;
class BemKernel : public RefCounter {
public:
    int typeKernel[2]={0,0};  // Laplace, Helmholtz, Yukawa
    // typeKernel ={SL, DL, HS, TDL} and determine equation Laplace, Helmholtz if k==0 or not
    std::complex<double> wavenum[2]={0,0}; // parameter to Helmholtz
    std::complex<double> coeffcombi[2]={0,0};
    BemKernel(){}
    BemKernel(string *tkernel, Complex alpha , Complex k) {
        
        coeffcombi[0]=alpha;
        wavenum[0]=k;
        
        if(!tkernel->compare("SL"))
            typeKernel[0] = 1;
        else if(!tkernel->compare("DL"))
            typeKernel[0] = 2;
        else if(!tkernel->compare("HS"))
            typeKernel[0] = 3;
        else if(!tkernel->compare("TDL"))
            typeKernel[0] = 4;
        else if(!tkernel->compare("CST"))
            typeKernel[0] = 5;
        else
            ExecError("unknown BEM kernel type ");
        
        if(mpirank==0 && verbosity>5)
            cout << "type BEM kernel " << *tkernel <<": " << typeKernel[0] << " coeff combi " << coeffcombi[0] << " wave number "<< wavenum[0] << endl;
    }
    ~BemKernel() {}
    
    BemKernel(const BemKernel &Bk) {
        for(int i=0;i<2;i++) {typeKernel[i]=Bk.typeKernel[i]; wavenum[i]=Bk.wavenum[i]; coeffcombi[i]=Bk.coeffcombi[i]; } } ;
    // alpha * ker
    BemKernel(Stack s,const BemKernel &Bk, Complex alpha) {
           typeKernel[0]=Bk.typeKernel[0]; wavenum[0]=Bk.wavenum[0]; coeffcombi[0]=alpha; } ;
       
    
private:
    //BemKernel(const BemKernel &);
    void operator=(const BemKernel &);
};
//class FoperatorKBEM;
// new type for bem
typedef const BemKernel *pBemKernel;
typedef const BemPotential *pBemPotential;

typedef  const BemKernel fkernel;
typedef  const BemPotential fpotential;

template<class ffmesh, class bemtoolmesh>
void Mesh2Bemtool(const ffmesh &Th, bemtool::Geometry &node, bemtoolmesh &mesh ) {
    
    typedef typename ffmesh::RdHat RdHat;
    typedef typename ffmesh::Element E;
    const int dHat =  RdHat::d;
    
    // create the geometry;
    
    bemtool::R3 p;
    for(int iv=0 ; iv<Th.nv ; iv++){
        p[0]=Th.vertices[iv].x;p[1]=Th.vertices[iv].y;p[2]=Th.vertices[iv].z;
        node.setnodes(p);
    }
    
    node.initEltData();
    
    if(mpirank==0 && verbosity>10) std::cout << "Creating mesh domain (nodes)" << std::endl;
    
    mesh.set_elt(node);
    bemtool::array<dHat+1,int> I;
    if(mpirank==0 && verbosity>10) std::cout << "End creating mesh domain mesh" << std::endl;
    
    if(mpirank==0 && verbosity>10) std::cout << "Creating geometry domain (elements)" << std::endl;
    for(int it=0; it<Th.nt; it++){
        const E &K(Th[it]);
        for(int j=0;j<dHat+1;j++)
            I[j]=Th.operator () (K[j]);
        mesh.setOneElt(node,I);
    }
    
    //mesh = unbounded;
    //Orienting(mesh);
    bemtool::Normal<dHat> N(mesh);
    for(int it=0; it<Th.nt; it++){
        const E &K(Th[it]);
        Fem2D::R3 nn = K.NormalTUnitaire();
        bemtool::R3 mm; mm[0]=nn.x; mm[1]=nn.y; mm[2]=nn.z;
        N.set(it, mm);
    }
    mesh.Orienting(N);
    
    if(mpirank==0 && verbosity>10) std::cout << "end creating geometry domain" << std::endl;
}

template<class ffmesh>
void Mesh2Bemtool(const ffmesh &Th, bemtool::Geometry &node) {
    if(mpirank==0 && verbosity>10) std::cout << "Creating mesh output" << std::endl;
    bemtool::R3 p;
    Fem2D::R3 pp;
    for(int iv=0 ; iv<Th.nv ; iv++){
        pp = Th.vertices[iv];
        p[0]=pp.x; p[1]=pp.y; p[2]=pp.z;
        node.setnodes(p);
    }
}



/// Domain integration - operator



class CPartBemDI: public E_F0mps {
public:
    
    static const int n_name_param =3; //12;
    static basicAC_F0::name_and_type name_param[] ;
    Expression nargs [n_name_param];
    enum typeofkind  { int1dx1d=0, int2dx2d=1, int2dx1d=2, int1dx2d=3  } ;
    typeofkind  kind; //  0
    int d, dHat; // 3d
    typedef const CPartBemDI* Result;
    Expression Th;
    vector<Expression> what;
    vector<int> whatis; // 0 -> long , 1 -> array ???
    
    
    CPartBemDI(const basicAC_F0 & args,typeofkind b=int1dx1d,int ddim=3, int ddimHat=1) // always ddim=3d
    :kind(b),d(ddim),dHat(ddimHat),
    Th(0), what(args.size()-1), whatis(args.size()-1)
    
    {
        args.SetNameParam(n_name_param,name_param,nargs);
        
        if(d==3 && dHat==1)
            Th=CastTo<pmeshL>(args[0]);
        else if(d==3 && dHat==2)
            Th=CastTo<pmeshS>(args[0]);
        else ffassert(0); // a faire
        
        int n=args.size();
        for (int i=1;i<n;i++)
            if(!BCastTo<KN_<long> >(args[i]) ) {
                whatis[i-1]=0;
                what[i-1]=CastTo<long>(args[i]);
            }
            else {
                whatis[i-1]=1;
                what[i-1]=CastTo<KN_<long> >(args[i]);
            }
        
    }
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh>(), true);} // all type
    AnyType operator()(Stack ) const  { return SetAny<const CPartBemDI *>(this);}
    
    operator aType () const { return atype<const CPartBemDI *>();}
    
    static  E_F0 * f(const basicAC_F0 & args) { return new CPartBemDI(args);}
    
    const Fem2D::QuadratureFormular & FIT(Stack) const ;
    const Fem2D::QuadratureFormular1d & FIE(Stack) const ;
    const Fem2D::GQuadratureFormular<Fem2D::R3> & FIV(Stack) const ;  // 3d
    
};



class CBemDomainOfIntegration: public E_F0mps {
public:
    
    static const int n_name_param =3; //12;
    static basicAC_F0::name_and_type name_param[] ;
    Expression nargs_t [n_name_param];
    // typeofkind  kind; //  0
    int d_s, dHat_s, d_t, dHat_t; // 3d
    typedef const CBemDomainOfIntegration* Result;
    Expression Th_s, Th_t, BemPartDI;
    vector<Expression> what_s, what_t;
    vector<int> whatis_s, whatis_t;
    //CBemDomainOfIntegration();
    CBemDomainOfIntegration(const basicAC_F0 & args_t, int ddim=3, int ddimHat=1) // always ddim=3d
    :d_s(0),dHat_s(0),d_t(ddim),dHat_t(ddimHat),
    Th_s(0), what_s(0), whatis_s(0),
    Th_t(0), what_t(args_t.size()-1), whatis_t(args_t.size()-1)
    
    {
        args_t.SetNameParam(n_name_param,name_param,nargs_t);
        // acces to the value of the fist di Th_s
        const CPartBemDI * sourceDI(dynamic_cast<const CPartBemDI*>((Expression) args_t[0]));
        
        CPartBemDI::typeofkind kind_s= sourceDI->kind; // int1dx1d=0, int2dx2d=1, int2dx1d=2, int1dx2d=3
        d_s=sourceDI->d;
        dHat_s=sourceDI->dHat;
        
        // check the integral operator integral
        if(kind_s==0 || kind_s==2)
            Th_t=CastTo<pmeshL>(args_t[1]);
        else if(kind_s==1 || kind_s==3)
            Th_t=CastTo<pmeshS>(args_t[1]);
        else if(kind_s==0 || kind_s==3)
            Th_s=sourceDI->Th; //Th_s=CastTo<pmeshL>(sourceDI->Th);
        else if(kind_s==1 || kind_s==2)
            Th_s=sourceDI->Th;
        else ffassert(0); // a faire
        
        if (mpirank==0 && verbosity >5)
            cout << " CBemDomainOfIntegration " << kind_s << " Th_s: " << &Th_s << " d_s= " << d_s << " dHat_s= " << dHat_s << " " <<
            " Th_t: " << &Th_t << " d_t= " << d_t << " dHat_t= " << dHat_t << " " << endl;
        
        // read the argument (Th_t,.....)
        int n_t=args_t.size();
        for (int i=2;i<n_t;i++)
            if(!BCastTo<KN_<long> >(args_t[i]) ) {
                whatis_t[i-1]=0;
                what_t[i-1]=CastTo<long>(args_t[i]);
            }
            else {
                whatis_t[i-1]=1;
                what_t[i-1]=CastTo<KN_<long> >(args_t[i]);
            }
        int n_s=(sourceDI->what).size();
        for (int i=1;i<n_s;i++) {
            
            whatis_s[i-1]=sourceDI->whatis[i-1];
            what_s[i-1]=sourceDI->what[i-1];
        }
        
    }
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<const CPartBemDI *>(), true);} // all type
    AnyType operator()(Stack ) const  { return SetAny<const CBemDomainOfIntegration *>(this);}
    
    operator aType () const { return atype<const CBemDomainOfIntegration *>();}
    
    //static  E_F0 * f(const basicAC_F0 & args_s,const basicAC_F0 & args_t) { return new CBemDomainOfIntegration(args_s,args_t);}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CBemDomainOfIntegration(args_s);}
    
    const Fem2D::QuadratureFormular & FIT(Stack) const ;
    const Fem2D::QuadratureFormular1d & FIE(Stack) const ;
    const Fem2D::GQuadratureFormular<Fem2D::R3> & FIV(Stack) const ;  // 3d
    
};



basicAC_F0::name_and_type  CPartBemDI::name_param[]= {
    { "qft", &typeid(const Fem2D::QuadratureFormular *)},    //0
    { "qfe", &typeid(const Fem2D::QuadratureFormular1d *)},
    { "qforder",&typeid(long)},     // 2
    //{ "qfnbpT",&typeid(long)},
    //{ "qfnbpE",&typeid(long)},
    //{ "optimize",&typeid(long)},
    //{ "binside",&typeid(double)},
    //{ "mortar",&typeid(bool)},
    { "qfV", &typeid(const Fem2D::GQuadratureFormular<Fem2D::R3> *)},    // 8 ->  3
    //{ "levelset",&typeid(double)},
    //{ "mapt",&typeid(E_Array)},
    //{ "mapu",&typeid(E_Array)}
    
    
};
basicAC_F0::name_and_type  CBemDomainOfIntegration::name_param[]= {
    { "qft", &typeid(const Fem2D::QuadratureFormular *)},    //0
    { "qfe", &typeid(const Fem2D::QuadratureFormular1d *)},
    { "qforder",&typeid(long)},     // 2
    //{ "qfnbpT",&typeid(long)},
    //{ "qfnbpE",&typeid(long)},
    //{ "optimize",&typeid(long)},
    //{ "binside",&typeid(double)},
    //{ "mortar",&typeid(bool)},
    { "qfV", &typeid(const Fem2D::GQuadratureFormular<Fem2D::R3> *)},    // 8 ->  3
    //{ "levelset",&typeid(double)},
    //{ "mapt",&typeid(E_Array)},
    //{ "mapu",&typeid(E_Array)}
    
    
};



class CPartBemDI1d1d: public CPartBemDI {
public:
    CPartBemDI1d1d( const basicAC_F0 & args_s) :CPartBemDI(args_s,int1dx1d,3,1) {}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CPartBemDI(args_s,int1dx1d,3,1);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshL>(), true);} // all type
};

class CPartBemDI2d2d: public CPartBemDI {
public:
    CPartBemDI2d2d(const basicAC_F0 & args_s) :CPartBemDI(args_s,int2dx2d,3,2) {}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CPartBemDI(args_s,int2dx2d,3,2);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshS>(), true);} // all type
};

class CPartBemDI1d2d: public CPartBemDI {
public:
    CPartBemDI1d2d(const basicAC_F0 & args_s) :CPartBemDI(args_s,int1dx2d,3,1) {}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CPartBemDI(args_s,int2dx2d,3,1);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshS>(), true);} // all type
};
class CPartBemDI2d1d: public CPartBemDI {
public:
    CPartBemDI2d1d(const basicAC_F0 & args_s) :CPartBemDI(args_s,int2dx1d,3,2) {}
    static  E_F0 * f(const basicAC_F0 & args_s) { return new CPartBemDI(args_s,int2dx2d,3,2);}
    static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmeshS>(), true);} // all type
};

//// begin type BEM kernel / potential



class listBemKernel {
public:
    list<BemKernel const  *> *lbk;
    void init()  { lbk=new list<BemKernel const  *>;}
    void destroy() { delete lbk;}
    listBemKernel(Stack s,BemKernel const*bk) : lbk(Add2StackOfPtr2Free(s,new list<BemKernel const*>)) { lbk->push_back(bk);}
    listBemKernel(Stack s,BemKernel const*const bka,BemKernel const* const bkb) : lbk(Add2StackOfPtr2Free(s,new list<BemKernel const*>)) { lbk->push_back(bka);lbk->push_back(bkb);}
    listBemKernel(Stack s,const listBemKernel &l,BemKernel const*const bk) : lbk(Add2StackOfPtr2Free(s,new list<BemKernel const*>(*l.lbk))) { lbk->push_back(bk);}
    listBemKernel(){};
};


template<class RR,class AA=RR,class BB=AA>
struct Op_addBemKernel: public binary_function<AA,BB,RR> {
    static RR f(Stack s,const AA & a,const BB & b) {
        if (mpirank==0 && verbosity>10) cout << "test " <<typeid(RR).name() << " " << typeid(AA).name() << " " << typeid(BB).name() <<endl;
        return RR(s,a,b);}
};


template<bool INIT,class RR,class AA=RR,class BB=AA>
struct Op_setBemKernel: public binary_function<AA,BB,RR> {
    static RR f(Stack stack, const AA & a,const BB & b)
    {
        ffassert(a);
        const pBemKernel p=combKernel(b);
        
        if (!INIT && *a)
            (**a).destroy( );
        *a = p;
        return a;
    }
};


template<class RR,class AA=RR,class BB=AA>
struct Op_coeffBemKernel1: public binary_function<AA,BB,RR> {
    static RR f(Stack s,const AA & a,const BB & b) {
        if (mpirank==0 && verbosity>10) cout << "test " <<typeid(RR).name() << " " << typeid(AA).name() << " " << typeid(BB).name() <<endl;
        
        RR ker=new BemKernel(s,*b,a);
        Add2StackOfPtr2Free(s,ker);
        return ker;}
};
// version ok
class OP_MakeBemKernel {
 public:
  
  class Op : public E_F0mps {
   public:
    static const int n_name_param = 1;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    Complex arg(int i, Stack stack, Complex a) const { return nargs[i] ? GetAny< Complex >((*nargs[i])(stack)) : a; }
    typedef pBemKernel *R;
    typedef pBemKernel *A;
    typedef string *B;
    Expression a, b;
      
    Op(const basicAC_F0 &args);

    AnyType operator( )(Stack s) const {
      A bemker = GetAny< A >((*a)(s));
      B type = GetAny< B >((*b)(s));
      Complex alpha=1.;//(arg(0, s, 1));
      Complex k(arg(0, s, 0));
      *bemker = new BemKernel(type,alpha,k);
      return SetAny< R >(bemker);
    }
  };

  typedef Op::R Result;
  static E_F0 *f(const basicAC_F0 &args) { return new Op(args); }
  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< Op::A >( ), atype< Op::B >( ), false);
  }
};



OP_MakeBemKernel::Op::Op(const basicAC_F0 &args)
  : a(to< A >(args[0])), b(to< B >(args[1])) {
  args.SetNameParam(n_name_param, name_param, nargs);
}

basicAC_F0::name_and_type OP_MakeBemKernel::Op::name_param[] = { {"k", &typeid(Complex)} };

class OP_MakeBemKernelFunc {
 public:
  
  class Op : public E_F0mps {
   public:
    static const int n_name_param = 1;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    Complex arg(int i, Stack stack, Complex a) const { return nargs[i] ? GetAny< Complex >((*nargs[i])(stack)) : a; }
    typedef pBemKernel A;
    typedef string *B;
    Expression b;
      
    Op(const basicAC_F0 &args);

    AnyType operator( )(Stack s) const {
      B type = GetAny< B >((*b)(s));
      Complex alpha=1.;//(arg(0, s, 1));
      Complex k(arg(0, s, Complex(0.,0.)));
      A bemker = new BemKernel(type,alpha,k);
      Add2StackOfPtr2Free(s,bemker);
      return SetAny< A >(bemker);
    }
  };

  typedef Op::A Result;
  static E_F0 *f(const basicAC_F0 &args) { return new Op(args); }
  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< Op::B >( ), false);
  }
};

OP_MakeBemKernelFunc::Op::Op(const basicAC_F0 &args)
  : b(to< B >(args[0])) {
   //   cout << "\n\n ****  OP_MakeBemKernelFunc::Op::Op() \n\n" << endl;
  args.SetNameParam(n_name_param, name_param, nargs);
}

basicAC_F0::name_and_type OP_MakeBemKernelFunc::Op::name_param[] = { {"k", &typeid(Complex)} };


class FormalBemKernel : public OneOperator{
public:
    
    FormalBemKernel( ): OneOperator(atype<C_F0>(),atype<string*>()) {}
    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const {
        Expression e=new OP_MakeBemKernelFunc::Op(args);
        aType r=atype<const BemKernel *>();
        return C_F0(e,r) ;}
    
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}

};

class BemPotential : public RefCounter {
public:
    Stack stack;
    int typePotential;  // Laplace, Helmholtz, Yukawa
    // typePotential ={SL=0, DL=1, HS=2, TDL=3} and determine equation Laplace, Helmholtz if k==0 or not
    std::complex<double> wavenum; // parameter to Helmholtz
    const OneOperator *codePot;
    mutable long s;
    C_F0 c_s;
    mutable long t;
    C_F0 c_t;
    Expression f;
    Data_Funct ** data;
    BemPotential(){}
    BemPotential(string *tpotential, Complex k) : typePotential(-1), wavenum(k) {
        
        if(!tpotential->compare("SL"))
            typePotential = 1;
        else if(!tpotential->compare("DL"))
            typePotential = 2;
        else if(!tpotential->compare("CST"))
            typePotential = 3;
        else
            ExecError("unknown BEM Potential type ");
        
        if(mpirank==0 && verbosity>5)
            cout << "type BEM Potential " << *tpotential <<": "<< tpotential <<": " << typePotential << " wave number "<< wavenum << endl;
    }

    BemPotential(Stack stk, const OneOperator* op, Data_Funct **d) : stack(stk),  s(0), c_s(CPValue(s)), t(0), c_t(CPValue(t)),
        f(op ? CastTo< double >(C_F0(op->code(basicAC_F0_wa({c_s,c_t})), (aType)*op)) : 0), data(d) {}

    ~BemPotential() {}
    
private:
    BemPotential(const BemPotential &);
    void operator=(const BemPotential &);
};


class OP_MakeBemPotential {
 public:
  class Op : public E_F0mps {
   public:
    static const int n_name_param = 1;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    Complex arg(int i, Stack stack, Complex a) const { return nargs[i] ? GetAny< Complex >((*nargs[i])(stack)) : a; }
    typedef pBemPotential *R;
    typedef pBemPotential *A;
    typedef Polymorphic *B;
    typedef Data_Funct **C;
    const OneOperator *codePot;
    Expression a, b, c;
    Op(const basicAC_F0 &args);

    AnyType operator( )(Stack s) const {
      A bempot = GetAny< A >((*a)(s));
      B type = GetAny< B >((*b)(s));
      C data = GetAny< C >((*c)(s));
      Complex k(arg(0, s, 0));
      *bempot = new BemPotential(s, codePot, data);
      return SetAny< R >(bempot);
    }
  };    // end Op class

  typedef Op::R Result;
  static E_F0 *f(const basicAC_F0 &args) { return new Op(args); }
  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< Op::A >( ), atype< Op::B >( ), atype< Op::C >( ), false);
  }
};

OP_MakeBemPotential::Op::Op(const basicAC_F0 &args)
  : a(to< A >(args[0])), b(to< B >(args[1])), c(to< C >(args[2])) {
  args.SetNameParam(n_name_param, name_param, nargs);
  const Polymorphic* op = dynamic_cast< const Polymorphic* >(args[1].LeftValue());
  codePot = op->Find("(", ArrayOfaType(atype<long>(),atype<long>(), false));
}

basicAC_F0::name_and_type OP_MakeBemPotential::Op::name_param[] = { {"k", &typeid(Complex)} };

class OP_MakeBemPotentialFunc {
 public:
  
  class Op : public E_F0mps {
   public:
    static const int n_name_param = 1;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    Complex arg(int i, Stack stack, Complex a) const { return nargs[i] ? GetAny< Complex >((*nargs[i])(stack)) : a; }
    typedef pBemPotential A;
    typedef string *B;
    Expression b;
      
    Op(const basicAC_F0 &args);

    AnyType operator( )(Stack s) const {
      B type = GetAny< B >((*b)(s));
      Complex k(arg(0, s, 0));
      A bempot = new BemPotential(type,k);
      return SetAny< A >(bempot);
    }
  };

  typedef Op::A Result;
  static E_F0 *f(const basicAC_F0 &args) { return new Op(args); }
  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< Op::B >( ), false);
  }
};

OP_MakeBemPotentialFunc::Op::Op(const basicAC_F0 &args)
  : b(to< B >(args[0])) {
  args.SetNameParam(n_name_param, name_param, nargs);
}

basicAC_F0::name_and_type OP_MakeBemPotentialFunc::Op::name_param[] = { {"k", &typeid(Complex)} };


class FormalBemPotential : public OneOperator{
public:
    
    FormalBemPotential( ): OneOperator(atype<C_F0>(),atype<Polymorphic*>()) {}
    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const {
        Expression e=new OP_MakeBemPotentialFunc::Op(args);
        aType r=atype<const BemPotential *>();
        return C_F0(e,r) ;}
    
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}

};


// fusion of Bem Kernel in case combined kernels
BemKernel *combKernel (listBemKernel const &lbemker){
    int kk=0;
    bool LaplaceK=false, HelmholtzK=false;
    const list< const BemKernel * > lbk(*lbemker.lbk);
    BemKernel *combBemKernel=new BemKernel();
    for (list< const BemKernel * >::const_iterator i = lbk.begin( ); i != lbk.end( ); i++) {
        if (!*i) continue;
        const BemKernel &bkb(**i);
        combBemKernel->typeKernel[kk] = bkb.typeKernel[0];
        // test same equation kernel
        (combBemKernel->wavenum[kk]!=0.) ? HelmholtzK=true : LaplaceK=true;
        
        //  ExecError(" combined kernel have to be the same type equation Laplace or Helmholtz");
        combBemKernel->coeffcombi[kk] = bkb.coeffcombi[0];
        combBemKernel->wavenum[kk] = bkb.wavenum[0];
        if(kk>4) ExecError(" combined kernel: 4 max kernels  ");
        kk++;
    }
    // check the same wave number ?
    if( HelmholtzK== LaplaceK) ExecError(" combined kernel: must be same equation kernels Laplace or Helmholtz");
    if (mpirank==0 && verbosity>5)
        for (int i=0;i<kk;i++) cout << "combined type BEM kernel " << combBemKernel->typeKernel[i] << " coeff combi " <<
            combBemKernel->coeffcombi[i] << " wave number "<< combBemKernel->wavenum[i] << endl;
    
    return combBemKernel;
}


// BEM variational form

// define a bilinear form for BEM
class FoperatorKBEM : public E_F0mps { public:
    typedef const FoperatorKBEM* Result;
    typedef finconnue * Fi;
    typedef ftest * Ft;
    typedef fkernel * KBem;
    
    Fi fi;
    Ft ft;
    Expression kbem;
    
    FoperatorKBEM(const basicAC_F0 & args) :fi(0),ft(0),kbem(0) {
        ffassert(args.size()==3);
 
        kbem= CastTo<KBem>(args[0]);
        fi= dynamic_cast<Fi>(CastTo<Fi>(args[1]));
        ft= dynamic_cast<Ft>(CastTo<Ft>(args[2]));
        ffassert(kbem && fi && ft);
    }
    
    AnyType operator()(Stack ) const { return SetAny<Result>(this);}
    operator aType () const { return atype<Result>();}
    FoperatorKBEM(const FoperatorKBEM & fk) : kbem(fk.kbem),fi(fk.fi),ft(fk.ft) {}
    
};

// define a bilinear form for BEM



class BemFormBilinear : virtual public E_F0mps { public:
    int type=-1;
    static  E_F0 * f(const basicAC_F0 & args);
};


class BemKFormBilinear : public BemFormBilinear { public:
    typedef const BemFormBilinear* Result;
    typedef const CBemDomainOfIntegration * A;
    typedef const FoperatorKBEM * B;
    A  di;
    FoperatorKBEM * b;
    
    BemKFormBilinear(const basicAC_F0 & args) {
        di= dynamic_cast<A>(CastTo<A>(args[0]));
        B Kb= dynamic_cast<B>(CastTo<B>(args[1]));
        b= new FoperatorKBEM(*Kb);
        ffassert(di && Kb);
        type=1;
        
        
    };
    
    static  E_F0 * f(const basicAC_F0 & args) { return new BemKFormBilinear(args);}
    static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}
    AnyType operator()(Stack ) const { return SetAny<Result>(this);}
    operator aType () const { return atype<Result>();}
    
    
    BemKFormBilinear(A a,Expression bb) : di(a),b(new FoperatorKBEM(*dynamic_cast<FoperatorKBEM *>(bb)))
    {ffassert(b);}
    BemKFormBilinear operator-() const { return  BemKFormBilinear(di,C_F0(TheOperators,"-",C_F0(b,atype<FoperatorKBEM>())));}
    
    BemKFormBilinear(const BemKFormBilinear & fb) : di(fb.di),b(new FoperatorKBEM(*fb.b) ) {}
    
};


class TypeFormBEM: public ForEachType<const BemFormBilinear*> {
public:
    TypeFormBEM() : ForEachType<const BemFormBilinear*>(0,0) {}
    void SetArgs(const ListOfId *lid) const {
        SetArgsFormLinear(lid,2);    }
    
    Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const
    { return Type_Expr(this,CastTo(c));}
    
    
    C_F0 Initialization(const Type_Expr & e) const
    {
        return C_F0(); }
    Type_Expr construct(const Type_Expr & e) const
    {
        return e; }
    
};

// define the function BEM(k,u,v)
class FormalKBEMcode : public OneOperator{
public:
    
    FormalKBEMcode( ): OneOperator(atype<C_F0>(),atype<pBemKernel>(), atype<finconnue*>(), atype<ftest*>()) {}
    FormalKBEMcode(int  ): OneOperator(atype<C_F0>(),atype<pBemKernel>()) {}
    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const {
        Expression e=new FoperatorKBEM(args);
        aType r=atype<const FoperatorKBEM *>();
        return C_F0(e,r) ;}
    
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    
};



// define the function POT(k,u,v)
class FoperatorPBEM : public E_F0mps { public:
    typedef const FoperatorPBEM* Result;
    typedef finconnue * Fi;
    typedef ftest * Ft;
    typedef fpotential * Pot;
    
    Fi fi;
    Ft ft;
    Expression pot;
    
    FoperatorPBEM(const basicAC_F0 & args) :fi(0),ft(0),pot(0) {
        ffassert(args.size()==3);
        
        pot= CastTo<Pot>(args[0]);
        fi= dynamic_cast<Fi>(CastTo<Fi>(args[1]));
        ft= dynamic_cast<Ft>(CastTo<Ft>(args[2]));
        
        ffassert(pot && fi && ft);
    }
    
    AnyType operator()(Stack ) const { return SetAny<Result>(this);}
    operator aType () const { return atype<Result>();}
    FoperatorPBEM(const FoperatorPBEM & fk) : pot(fk.pot),fi(fk.fi),ft(fk.ft){}
    
};


// define a bilinear form for BEM
class BemPFormBilinear : public BemFormBilinear { public:
    typedef const BemFormBilinear* Result;
    typedef const CDomainOfIntegration * A;
    typedef const FoperatorPBEM * B;
    A  di;
    FoperatorPBEM * b;
    BemPFormBilinear(const basicAC_F0 & args) {
        di= dynamic_cast<A>(CastTo<A>(args[0]));
        B Kb= dynamic_cast<B>(CastTo<B>(args[1]));
        b= new FoperatorPBEM(*Kb);
        ffassert(di && Kb);
        type=2;
    }
    
    static  E_F0 * f(const basicAC_F0 & args) { return new BemPFormBilinear(args);}
    static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
    AnyType operator()(Stack ) const { return SetAny<Result>(this);}
    operator aType () const { return atype<Result>();}
    
    
    BemPFormBilinear(A a,Expression bb) : di(a),b(new FoperatorPBEM(*dynamic_cast<FoperatorPBEM *>(bb))/*->Optimize(currentblock) FH1004 */)
    {ffassert(b);}
    BemPFormBilinear operator-() const { return  BemPFormBilinear(di,C_F0(TheOperators,"-",C_F0(b,atype<FoperatorPBEM>())));}
    
    BemPFormBilinear(const BemPFormBilinear & fb) : di(fb.di),b(new FoperatorPBEM(*fb.b) ) {}
    
};

class FormalPBEMcode : public OneOperator{
public:
    
    FormalPBEMcode( ): OneOperator(atype<C_F0>(),atype<pBemPotential>(), atype<finconnue*>(), atype<ftest*>()) {}
    FormalPBEMcode(int  ): OneOperator(atype<C_F0>(),atype<pBemPotential>()) {}
    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const {
        Expression e=new FoperatorPBEM(args);
        aType r=atype<const FoperatorPBEM *>();
        return C_F0(e,r) ;}
    
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    
};


//// end type BEM kernel / potential
int typeVFBEM(const list<C_F0> & largs, Stack stack)
{
    list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end();
    
    int VVFBEM =-1, ik=-1;;
    for (ii=ib;ii != ie;ii++) {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        
        if (r==atype<const  BemFormBilinear *>()) {
            
            BemFormBilinear * bb= GetAny<BemFormBilinear *>((*e)(0));
            VVFBEM = bb->type;
        }
        ffassert(ik);
    }
    return VVFBEM;
}
bool iscombinedKernel(BemKernel *K ) {
    bool iscombined=false;
    for (int i=0;i<2;i++)
        iscombined= K->coeffcombi[i]!=0. ? 1 : 0;
    return iscombined;
}
bemtool::BIOpKernelEnum whatTypeEnum(BemKernel *K,int i) {
    bemtool::BIOpKernelEnum pKernel;
    switch(K->typeKernel[i]) {
        case 1: pKernel=bemtool::SL_OP ; break;
        case 2: pKernel=bemtool::DL_OP ; break;
        case 3: pKernel=bemtool::HS_OP ; break;
        case 4: pKernel=bemtool::TDL_OP ; break;
    }
    const bemtool::BIOpKernelEnum cpKernel=pKernel;
    return cpKernel;
}

pair<BemKernel*,double> getBemKernel(Stack stack, const list<C_F0> & largs)  {
    list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end();
    
    BemKernel* K;
    double alpha=0.;
    
    bool haveBemBilinearOperator=false, haveBilinearOperator=false;
    
    for (ii=ib;ii != ie;ii++) {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        
        if (r==atype<const  BemFormBilinear *>() && !haveBemBilinearOperator) {
            BemKFormBilinear * bb=new BemKFormBilinear(*dynamic_cast<const BemKFormBilinear *>(e));
            FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
            if (b == NULL) {
                if(mpirank == 0) cout << "dynamic_cast error" << endl; }
            else
                K=GetAny<BemKernel*>((*b->kbem)(stack));
            haveBemBilinearOperator=true;
        }
        
        else if (r==atype<const  FormBilinear *>() && !haveBilinearOperator) {
            
            const  FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
            const CDomainOfIntegration & di= *bb->di;
            // check the integration (keyword)
            ffassert( (di.kind == CDomainOfIntegration::int1d && di.dHat==1) || (di.kind == CDomainOfIntegration::int2d && di.dHat==2) );  //check only necessary in surface case
            
            BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
            if (Op == NULL) {
                if(mpirank == 0) cout << "dynamic_cast error" << endl; }
            
            
            
            BilinearOperator::const_iterator l=Op->v.begin();
            
            BilinearOperator::K ll(*l);   //  LinearComb<pair<MGauche,MDroit>,C_F0> BilinearOperator;
            pair<int,int> finc(ll.first.first),ftest(ll.first.second);
                                   
            alpha = GetAny<double>(ll.second.eval(stack));
            if(mpirank == 0 && verbosity>5) cout << " test coeff mass matrix " << alpha << endl;
          
            if(mpirank == 0 && verbosity>5) {
                cout << "FormBilinear: number of unknown finc=" << finc.first << " ,ftest= " << ftest.first << endl;
                cout << "FormBilinear: operator order finc=" << finc.second << " ,ftest= " << ftest.second << endl;      // ordre   only op_id=0
            }
            ffassert(finc.first==0 && ftest.first==0);
            ffassert(finc.second==0 && ftest.second==0);
            haveBilinearOperator=true; // or ffassert(Op->v.size()==1); // check size 1
        }
        else
            ffassert(0);
    }
    return std::make_pair(K, alpha); //K;
}
BemPotential* getBemPotential(Stack stack, const list<C_F0> & largs)  {
    list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end();
    
    BemPotential* P;
    
    for (ii=ib;ii != ie;ii++) {
        Expression e=ii->LeftValue();
        aType r = ii->left();
        ffassert (r==atype<const  BemFormBilinear *>());
        
        BemPFormBilinear * bb=new BemPFormBilinear(*dynamic_cast<const BemPFormBilinear *>(e));
        FoperatorPBEM * b=const_cast<  FoperatorPBEM *>(bb->b);
        if (b == NULL) {
            if(mpirank == 0) cout << "dynamic_cast error" << endl; }
        else
            P=GetAny<BemPotential*>((*b->pot)(stack));
    }
    return P;
}

bemtool::PotKernelEnum whatTypeEnum(BemPotential *P) {
    bemtool::PotKernelEnum pPotential;
    switch(P->typePotential) {
        case 1: pPotential=bemtool::SL_POT ; break;
        case 2: pPotential=bemtool::DL_POT ; break;
        case 3: pPotential=bemtool::CST_POT ; break;
            //case 4: pPotential=TDL_POT ; break;
    }
    const bemtool::PotKernelEnum cpPotential=pPotential;
    return cpPotential;
}

template <class K, typename P, class MMesh>
void ff_BIO_Generator(htool::VirtualGenerator<K>*& generator, BemKernel *typeKernel, bemtool::Dof<P>& dof, double alpha) {
    
    bemtool::BIOpKernelEnum ker1 = whatTypeEnum(typeKernel,0), ker2 = whatTypeEnum(typeKernel,1);;
    double kappaRe1 = typeKernel->wavenum[0].real(), kappaRe2 = typeKernel->wavenum[1].real();
    double kappaIm1 = typeKernel->wavenum[0].imag(), kappaIm2 = typeKernel->wavenum[1].imag();
    
    bool iscombined = iscombinedKernel(typeKernel);
    if(iscombined) ffassert( (kappaRe1==kappaRe2) && (kappaIm1==kappaIm2) );
    std::complex<double> coeff1=typeKernel->coeffcombi[0], coeff2=typeKernel->coeffcombi[1];
    
    // BIO_Generator -> single kernel
    // Equ Helmholtz kappa1.real() > 0 et kappa1.imag() == 0
    if ( (kappaRe1 && !kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && !alpha ) {
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1);
                if(mpirank == 0 && verbosity>5) cout << " call bemtool func BIOpKernel<HE,SL_OP ..." << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<HE,DL_OP ..." << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<HE,HS_OP ..." << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<HE,TDL_OP ..." << endl; break;
            case bemtool::CST_OP :  if(verbosity>5) cout << " no tested BIOpKernel<HE,CST_OP" << endl; ffassert(0); break;
        }
    }
    // Eq Laplace  kappa1 == 0  kappaRe1=0
    else if ( (!kappaRe1 && !kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && !alpha ) {
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,SL_OP ..." << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,TDL_OP ..." << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,HS_OP ..." << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,TDL_OP ..." << endl; break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<LA,CST_OP" << endl; ffassert(0); break;
        }
    }
    // Eq Yukawa  kappa1.real() == 0 et kappa1.imag() > 0
    else if ( (!kappaRe1 && kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && !alpha ) {
        //(!kappa1 && !iscombined && !kappa2 && !alpha ){
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<YU,SL_OP ..." << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<YU,TDL_OP ..." << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<YU,HS_OP ..." << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<YU,TDL_OP ..." << endl; break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<YU,CST_OP" << endl; ffassert(0); break;
        }
    }

    //BIO_Generator + alpha * mass_matrix
    // Equ HELMHOLTZ kappa1.real() > 0 et kappa1.imag() == 0
    else if ( (kappaRe1 && !kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha ) {
        //if (kappa1 && !iscombined && !kappa2 && alpha) {
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha);
                if(mpirank == 0 && verbosity>5) cout << " call bemtool func BIO_Generator_w_mass<HE,SL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<HE,DL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<HE,HS_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<HE,TDL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case bemtool::CST_OP :  if(mpirank == 0 && verbosity>5) cout << " no tested BIO_Generator_w_mass<HE,CST_OP" << endl; ffassert(0); break;
        }
    }
    // Eq LAPLACE  kappa1 == 0   kappaRe1=0
    else if ( (!kappaRe1 && !kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha ) { //if (!kappa1 && !iscombined && !kappa2 && alpha){
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,SL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,TDL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,HS_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,TDL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<LA,CST_OP" << endl; ffassert(0); break;
        }
    }
    // Eq YUKAWA  kappa1.real() == 0 et kappa1.imag() > 0
    else if ( (!kappaRe1 && kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha ) { //if (!kappa1 && !iscombined && !kappa2 && alpha){
       switch (ker1) {
           case bemtool::SL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,alpha);
               if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<YU,SL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
           case bemtool::DL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,alpha);
               if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<YU,TDL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
           case bemtool::HS_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,alpha);
               if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<YU,HS_OP ...alpha coeff mass matrix=" << alpha << endl; break;
           case bemtool::TDL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,alpha);
               if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<YU,TDL_OP ...alpha coeff mass matrix=" << alpha << endl; break;
           case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<YU,CST_OP" << endl; ffassert(0); break;
       }
    }

    // combined kernel ... coeff1 ker1 + coeff2 ker2 + alpha * mass_matrix

    // eq HELMHOLTZ
    else if ( (kappaRe1 && !kappaIm1) && (kappaRe2 &&!kappaIm2) && iscombined && alpha) {
        
        switch (ker1) {
            case bemtool::SL_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,SL_OP ... HE,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,SL_OP ... HE,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,SL_OP ... HE,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,SL_OP ... HE,TDL_OP" << endl; break;
                    case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<HE,SL_OP ... HE,CST_OP" << endl; ffassert(0); break;
                }
                break;
            case bemtool::DL_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... HE,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... HE,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... HE,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... HE,TDL_OP" << endl; break;
                    case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<HE,DL_OP ... HE,CST_OP" << endl; ffassert(0); break;
                }
                    break;
            case bemtool::HS_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,HS_OP ... HE,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,HS_OP ... HE,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,HS_OP ... HE,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,HS_OP ... HE,TDL_OP" << endl; break;
                    case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<HE,HS_OP ... HE,CST_OP" << endl; ffassert(0); break;
                }
                    break;
            case bemtool::TDL_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,TDL_OP ... HE,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,TDL_OP ... HE,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,TDL_OP ... HE,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,TDL_OP ... HE,TDL_OP" << endl; break;
                    case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<HE,TDL_OP ... CST_OP" << endl; ffassert(0); break;
                }
                    break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<HE,CST_OP" << endl; ffassert(0); break;
                
        }
    }
    
    // Eq LAPLACE
    else if ( (!kappaRe1 && !kappaIm1) && (!kappaRe2 && !kappaIm2) && iscombined && alpha) {
        
        switch (ker1) {
            case bemtool::SL_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,SL_OP ... LA,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,SL_OP ... LA,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,SL_OP ... LA,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,SL_OP ... LA,TDL_OP" << endl; break;
                    case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,SL_OP ... LA,CST_OP" << endl; ffassert(0); break;
                }
            case bemtool::DL_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,DL_OP ... LA,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,DL_OP ... LA,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,DL_OP ... LA,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,DL_OP ... LA,TDL_OP" << endl; break;
                    case bemtool::CST_OP :  if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,DL_OP ... LA,CST_OP" << endl; ffassert(0); break;
                }
            case bemtool::HS_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,HS_OP ... LA,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,HS_OP ... LA,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,HS_OP ... LA,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,HS_OP ... LA,TDL_OP" << endl; break;
                    case bemtool::CST_OP :  if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,HS_OP ... LA,CST_OP" << endl; ffassert(0); break;
                }
            case bemtool::TDL_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,TDL_OP ... LA,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,TDL_OP ... LA,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,TDL_OP ... LA,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<LA,TDL_OP ... LA,TDL_OP" << endl; break;
                    case bemtool::CST_OP :  if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,TDL_OP ... LA,CST_OP" << endl; ffassert(0); break;
                }
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,CST_OP ... LA,CST_OP" << endl; ffassert(0); break;
        }
    }
    
    
    // Eq YUKAWA
    else if ( (!kappaRe1 && kappaIm1) && (!kappaRe2 && kappaIm2) && iscombined && alpha) {
        
        switch (ker1) {
            case bemtool::SL_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,SL_OP ... YU,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,SL_OP ... YU,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,SL_OP ... YU,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,SL_OP ... YU,TDL_OP" << endl; break;
                    case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<YU,SL_OP ... YU,CST_OP" << endl; ffassert(0); break;
                }
                break;
            case bemtool::DL_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... YU,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... YU,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<HE,DL_OP ... YU,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,DL_OP ... YU,TDL_OP" << endl; break;
                    case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<YU,DL_OP ... YU,CST_OP" << endl; ffassert(0); break;
                }
                    break;
            case bemtool::HS_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,HS_OP ... YU,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,HS_OP ... YU,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,HS_OP ... YU,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,HS_OP ... YU,TDL_OP" << endl; break;
                    case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<YU,HS_OP ... YU,CST_OP" << endl; ffassert(0); break;
                }
                    break;
            case bemtool::TDL_OP :
                switch (ker2) {
                    case bemtool::SL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,TDL_OP ... YU,SL_OP" << endl; break;
                    case bemtool::DL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,TDL_OP ... YU,DL_OP" << endl; break;
                    case bemtool::HS_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,TDL_OP ... YU,HS_OP" << endl; break;
                    case bemtool::TDL_OP : generator= new bemtool::Combined_BIO_Generator<bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1,coeff2,alpha);
                        if(mpirank == 0 && verbosity>5) cout << "call bemtool func Combined_BIO_Generator<YU,TDL_OP ... YU,TDL_OP" << endl; break;
                    case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested Combined_BIO_Generator<YU,TDL_OP ... CST_OP" << endl; ffassert(0); break;
                }
                    break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<YU,CST_OP" << endl; ffassert(0); break;
                
        }
    }
    
    else {
        if(mpirank == 0) cout << "kernel definition error" << endl; ffassert(0);}
    
    
}

template <class R, typename P, typename MeshBemtool, class MMesh>
void ff_POT_Generator(htool::VirtualGenerator<R>*& generator,BemPotential *typePot, bemtool::Dof<P> &dof, MeshBemtool &mesh, bemtool::Geometry &node_output) {
    
    bemtool::PotKernelEnum pot = whatTypeEnum(typePot);
    double kappaRe = typePot->wavenum.real(),kappaIm = typePot->wavenum.imag();
    
    switch (pot) {
        case bemtool::SL_POT :
            if (kappaRe && !kappaIm) {
                generator = new bemtool::POT_Generator<bemtool::PotKernel<HE,bemtool::SL_POT,MMesh::RdHat::d+1,P>,P>(dof,node_output,kappaRe);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<HE,SL_POT ..." << endl;
            }
            else if (!kappaRe && kappaIm) {
                generator = new bemtool::POT_Generator<bemtool::PotKernel<YU,bemtool::SL_POT,MMesh::RdHat::d+1,P>,P>(dof,node_output,kappaIm);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<YU,SL_POT ..." << endl;
            }
            else if (!kappaRe && !kappaIm){  //kappaRe=0
                generator = new bemtool::POT_Generator<bemtool::PotKernel<LA,bemtool::SL_POT,MMesh::RdHat::d+1,P>,P>(dof,node_output,kappaRe);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<LA,SL_POT ..." << endl;
            }
            break;
            
        case bemtool::DL_POT :
            if (kappaRe && !kappaIm) {
                generator = new bemtool::POT_Generator<bemtool::PotKernel<HE,bemtool::DL_POT,MMesh::RdHat::d+1,P>,P>(dof,node_output,kappaRe);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<HE,DL_POT ..." << endl;
            }
            if (!kappaRe && kappaIm) {
                           generator = new bemtool::POT_Generator<bemtool::PotKernel<YU,bemtool::DL_POT,MMesh::RdHat::d+1,P>,P>(dof,node_output,kappaIm);
                           if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<YU,DL_POT ..." << endl;
                       }
            else if (!kappaRe && !kappaIm){ //kappaRe=0
                generator = new bemtool::POT_Generator<bemtool::PotKernel<LA,bemtool::DL_POT,MMesh::RdHat::d+1,P>,P>(dof,node_output,kappaRe);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<LA,DL_POT ..." << endl;
            }
            break;
            
        case bemtool::CST_POT :
            if (kappaRe && !kappaIm) {
                if(mpirank == 0 && verbosity>5) cout << " no POT_Generator<HE,CST_POT ... " << endl; ffassert(0); break; }
            else if (!kappaRe && kappaIm) {
                if(mpirank == 0 && verbosity>5) cout << " no POT_Generator<YU,CST_POT ... " << endl; ffassert(0); break; }
            else if (!kappaRe && !kappaIm) {
                if(mpirank == 0 && verbosity>5) cout << " no POT_Generator<LA,CST_POT ... " << endl; ffassert(0); break; }
            
    }
}
