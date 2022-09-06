#ifndef BEM_HPP_
#define BEM_HPP_

class BemKernel;
class BemPotential;
class BemKernel : public RefCounter {
public:
    int typeKernel[2]={0,0};  // Laplace, Helmholtz, Yukawa
    // typeKernel ={SL, DL, HS, TDL} and determine equation Laplace, Helmholtz if k==0 or not
    std::complex<double> wavenum[2]={0.0+1i*0.0,0.0+1i*0.0}; // parameter to Helmholtz
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
struct Op_addBemKernel {
    using first_argument_type  = AA;
    using second_argument_type = BB;
    using result_type          = RR;
    static RR f(Stack s,const AA & a,const BB & b) {
        if (mpirank==0 && verbosity>10) cout << "test " <<typeid(RR).name() << " " << typeid(AA).name() << " " << typeid(BB).name() <<endl;
        return RR(s,a,b);}
};


template<bool INIT,class RR,class AA=RR,class BB=AA>
struct Op_setBemKernel {
    using first_argument_type  = AA;
    using second_argument_type = BB;
    using result_type          = RR;
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
struct Op_coeffBemKernel1 {
    using first_argument_type  = AA;
    using second_argument_type = BB;
    using result_type          = RR;
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
      Complex alpha(1.0,0.0);//(arg(0, s, 1));
      Complex k(arg(0, s, Complex(0.,0.)));
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
    
    int typePotential;  // Laplace, Helmholtz, Yukawa
    // typePotential ={SL=0, DL=1, HS=2, TDL=3} and determine equation Laplace, Helmholtz if k==0 or not
    std::complex<double> wavenum; // parameter to Helmholtz

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
    typedef string *B;
    Expression a, b;
    Op(const basicAC_F0 &args);

    AnyType operator( )(Stack s) const {
      A bempot = GetAny< A >((*a)(s));
      B type = GetAny< B >((*b)(s));
      Complex k(arg(0, s, Complex(0.0,0.0)));
      *bempot = new BemPotential(type,k); 
      return SetAny< R >(bempot);
    }
  };    // end Op class

  typedef Op::R Result;
  static E_F0 *f(const basicAC_F0 &args) { return new Op(args); }
  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< Op::A >( ), atype< Op::B >( ), false);
  }
};

OP_MakeBemPotential::Op::Op(const basicAC_F0 &args)
  : a(to< A >(args[0])), b(to< B >(args[1])) {
  args.SetNameParam(n_name_param, name_param, nargs);
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
      Complex k(arg(0, s, Complex(0.0,0.0)));
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
    
    FormalBemPotential( ): OneOperator(atype<C_F0>(),atype<string*>()) {}
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
    FoperatorKBEM( Expression o_kbem, const Finconnue &o_fi, const Ftest & o_ft)  : kbem(o_kbem){
        fi = new Finconnue(o_fi);
        ft = new Ftest(o_ft);
    }
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
    
    BemKFormBilinear(const BemKFormBilinear & fb) : di(fb.di),b(new FoperatorKBEM(*fb.b) ){ type=1;}
    
    BemKFormBilinear(A a, const FoperatorKBEM &fb) : di(a),b(new FoperatorKBEM(fb))
    { type=1;}  //exit(0); }
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
    
    BemPFormBilinear(const BemPFormBilinear & fb) : di(fb.di),b(new FoperatorPBEM(*fb.b) ) {type=2;}
    
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

pair<BemKernel*,Complex> getBemKernel(Stack stack, const list<C_F0> & largs)  {
    list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end();
    
    BemKernel* K;
    Complex alpha=0.;
    
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
            //if(verbosity >5){
            cout << "KKKK combi=" << K->coeffcombi[0] << " " << K->coeffcombi[1] << endl;
            //}
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
                                   
            alpha = GetAny<Complex>(ll.second.eval(stack));
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
void ff_BIO_Generator(htool::VirtualGenerator<K>*& generator, BemKernel *typeKernel, bemtool::Dof<P>& dof, Complex alpha) {
    
    if(mpirank == 0) cout << "ff_BIO::LaplaceHelmoltz:: " << endl;

    bemtool::BIOpKernelEnum ker1 = whatTypeEnum(typeKernel,0), ker2 = whatTypeEnum(typeKernel,1);;
    double kappaRe1 = typeKernel->wavenum[0].real(), kappaRe2 = typeKernel->wavenum[1].real();
    double kappaIm1 = typeKernel->wavenum[0].imag(), kappaIm2 = typeKernel->wavenum[1].imag();
    
    bool iscombined = iscombinedKernel(typeKernel);
    if(iscombined) ffassert( (kappaRe1==kappaRe2) && (kappaIm1==kappaIm2) );
    std::complex<double> coeff1=typeKernel->coeffcombi[0], coeff2=typeKernel->coeffcombi[1];
    
    if(mpirank == 0) cout << "ff_BIO::LaplaceHelmoltz:: alpha=" << alpha << ", coeff1=" << coeff1 << endl;
    
    // BIO_Generator -> single kernel
    // Equ Helmholtz kappa1.real() > 0 et kappa1.imag() == 0
    if ( (kappaRe1 && !kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha == 0. ) {
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << " call bemtool func BIOpKernel<HE,SL_OP ..." << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<HE,DL_OP ..." << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<HE,HS_OP ..." << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<HE,TDL_OP ..." << endl; break;
            case bemtool::CST_OP :  if(verbosity>5) cout << " no tested BIOpKernel<HE,CST_OP" << endl; ffassert(0); break;
        }
    }
    // Eq Laplace  kappa1 == 0  kappaRe1=0
    else if ( (!kappaRe1 && !kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha == 0. ) {
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,SL_OP ..." << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,TDL_OP ..." << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,HS_OP ..." << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<LA,TDL_OP ..." << endl; break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<LA,CST_OP" << endl; ffassert(0); break;
        }
    }
    // Eq Yukawa  kappa1.real() == 0 et kappa1.imag() > 0
    else if ( (!kappaRe1 && kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha == 0. ) {
        //(!kappa1 && !iscombined && !kappa2 && alpha == 0. ){
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<YU,SL_OP ..." << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<YU,TDL_OP ..." << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<YU,HS_OP ..." << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIOpKernel<YU,TDL_OP ..." << endl; break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<YU,CST_OP" << endl; ffassert(0); break;
        }
    }

    //BIO_Generator + alpha * mass_matrix
    // Equ HELMHOLTZ kappa1.real() > 0 et kappa1.imag() == 0
    else if ( (kappaRe1 && !kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha != 0. ) {
        //if (kappa1 && !iscombined && !kappa2 && alpha != 0.) {
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<HE,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << " call bemtool func BIO_Generator_w_mass<HE,SL_OP ...alpha coeff mass matrix=" << alpha << ", coeff1="<< coeff1 <<endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<HE,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<HE,DL_OP ...alpha coeff mass matrix=" << alpha << ", coeff1="<< coeff1 << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<HE,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<HE,HS_OP ...alpha coeff mass matrix=" << alpha << ", coeff1="<< coeff1 << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<HE,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<HE,TDL_OP ...alpha coeff mass matrix=" << alpha << ", coeff1="<< coeff1 << endl; break;
            case bemtool::CST_OP :  if(mpirank == 0 && verbosity>5) cout << " no tested BIO_Generator_w_mass<HE,CST_OP" << endl; ffassert(0); break;

        }
    }
    // Eq LAPLACE  kappa1 == 0   kappaRe1=0
    else if ( (!kappaRe1 && !kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha != 0. ) { //if (!kappa1 && !iscombined && !kappa2 && alpha != 0.){
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<LA,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,SL_OP ...alpha coeff mass matrix=" << alpha << ", coeff1="<< coeff1 << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<LA,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,TDL_OP ...alpha coeff mass matrix=" << alpha << ", coeff1="<< coeff1 << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<LA,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,HS_OP ...alpha coeff mass matrix=" << alpha << ", coeff1="<< coeff1 << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<LA,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaRe1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<LA,TDL_OP ...alpha coeff mass matrix=" << alpha << ", coeff1="<< coeff1 << endl; break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<LA,CST_OP" << endl; ffassert(0); break;
        }
    }
    // Eq YUKAWA  kappa1.real() == 0 et kappa1.imag() > 0
    else if ( (!kappaRe1 && kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha != 0. ) { //if (!kappa1 && !iscombined && !kappa2 && alpha != 0.){
       switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<YU,bemtool::SL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<YU,SL_OP ...alpha coeff mass matrix=" << alpha <<", coeff1="<< coeff1 << endl; break;
            case bemtool::DL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<YU,bemtool::DL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<YU,TDL_OP ...alpha coeff mass matrix=" << alpha <<", coeff1="<< coeff1 << endl; break;
            case bemtool::HS_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<YU,bemtool::HS_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<YU,HS_OP ...alpha coeff mass matrix=" << alpha <<", coeff1="<< coeff1 << endl; break;
            case bemtool::TDL_OP : generator=new bemtool::BIO_Generator_w_mass<bemtool::BIOpKernel<YU,bemtool::TDL_OP,MMesh::RdHat::d+1,P,P>,P>(dof,kappaIm1,alpha,coeff1);
                if(mpirank == 0 && verbosity>5) cout << "call bemtool func BIO_Generator_w_mass<YU,TDL_OP ...alpha coeff mass matrix=" << alpha <<", coeff1="<< coeff1 << endl; break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no tested BIOpKernel<YU,CST_OP" << endl; ffassert(0); break;
       }
    }

    // combined kernel ... coeff1 ker1 + coeff2 ker2 + alpha * mass_matrix

    // eq HELMHOLTZ
    else if ( (kappaRe1 && !kappaIm1) && (kappaRe2 &&!kappaIm2) && iscombined && alpha != 0.) {
        
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
    else if ( (!kappaRe1 && !kappaIm1) && (!kappaRe2 && !kappaIm2) && iscombined && alpha != 0.) {
        
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
                break;
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
                break;
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
                break;
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
                break;
            case bemtool::CST_OP : if(mpirank == 0 && verbosity>5) cout << " no testedCombined_BIO_Generator<LA,CST_OP ... LA,CST_OP" << endl; ffassert(0); break;
        }
    }
    
    
    // Eq YUKAWA
    else if ( (!kappaRe1 && kappaIm1) && (!kappaRe2 && kappaIm2) && iscombined && alpha != 0.) {
        
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
        if(mpirank == 0) cout << "Laplace Helmoltz:: kernel definition error" << endl; ffassert(0);}
    
    
}

/*
template<int d>
struct is_dim_2{static constexpr bool value = false;};
template<>
struct is_dim_2<2>{static constexpr bool value = true;};

template<int d>
struct is_dim_1{static constexpr bool value = false;};
template<>
struct is_dim_1<1>{static constexpr bool value = true;};

template <class K, class mesh, typename std::enable_if_t<is_dim_1<mesh::dim>::value>* = nullptr>
void ff_BIO_Generator(htool::VirtualGenerator<K>*& generator, BemKernel *typeKernel, mesh& m, double alpha) {
    return;
}
*/

template <class K, class mesh>
void ff_BIO_Generator_Maxwell(htool::VirtualGenerator<K>*& generator, BemKernel *typeKernel, bemtool::Dof<bemtool::RT0_2D>& dof, Complex alpha) {
    cout << " mettre un msg d erreur pour dire que cette combi n existe pas" << endl;
    return;
}

template <class K>
void ff_BIO_Generator_Maxwell(htool::VirtualGenerator<K>*& generator, BemKernel *typeKernel, bemtool::Dof<bemtool::RT0_2D>& dof, Complex alpha) {

    if(mpirank == 0) cout << "ff_BIO::Maxwell:: " << endl;

    bemtool::BIOpKernelEnum ker1 = whatTypeEnum(typeKernel,0), ker2 = whatTypeEnum(typeKernel,1);;
    double kappaRe1 = typeKernel->wavenum[0].real(), kappaRe2 = typeKernel->wavenum[1].real();
    double kappaIm1 = typeKernel->wavenum[0].imag(), kappaIm2 = typeKernel->wavenum[1].imag();

    bool iscombined = iscombinedKernel(typeKernel);
    if(iscombined) ffassert( (kappaRe1==kappaRe2) && (kappaIm1==kappaIm2) );
    std::complex<double> coeff1=typeKernel->coeffcombi[0], coeff2=typeKernel->coeffcombi[1];


    // BIO_Generator -> single kernel
    // Equ Helmholtz kappa1.real() > 0 et kappa1.imag() == 0
    if ( (kappaRe1 && !kappaIm1) && !iscombined && (!kappaRe2 && !kappaIm2) && alpha == 0. ) {
        switch (ker1) {
            case bemtool::SL_OP : generator=new bemtool::BIO_Generator<bemtool::BIOpKernel<MA,bemtool::SL_OP,3,bemtool::RT0_2D,bemtool::RT0_2D>,bemtool::RT0_2D>(dof,kappaRe1);
                if(mpirank == 0 && verbosity>5) cout << " call bemtool func BIOpKernel<HE,SL_OP ..." << endl; break;
            default: ffassert(0);
        }
    }
    else {
        if(mpirank == 0) cout << "Maxwell:: kernel definition error" << endl; ffassert(0);}
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

template <class R, typename P, typename MeshBemtool >
void ff_POT_Generator_Maxwell(htool::VirtualGenerator<R>*& generator,BemPotential *typePot, bemtool::Dof<P> &dof, MeshBemtool &mesh, bemtool::Geometry &node_output ){
    cout << " mettre un msg d erreur pour dire que cette combi n existe pas" << endl;
    return;
}

template <class R, typename P>
void ff_POT_Generator_Maxwell(htool::VirtualGenerator<R>*& generator,BemPotential *typePot, bemtool::Dof<P> &dof, bemtool::Mesh2D &mesh, bemtool::Geometry &node_output ) {
    
    bemtool::PotKernelEnum pot = whatTypeEnum(typePot);
    double kappaRe = typePot->wavenum.real(),kappaIm = typePot->wavenum.imag();
    if(mpirank == 0 && verbosity > 5) cout << "typePot->wavenum=" << typePot->wavenum << endl;
    
    switch (pot) {
        case bemtool::SL_POT :
            if (kappaRe && !kappaIm) {
                
                generator = new bemtool::POT_Generator<bemtool::PotKernel<MA,bemtool::SL_POT,3,bemtool::RT0_2D>,bemtool::RT0_2D>(dof,node_output,kappaRe);

                if(mpirank == 0 && verbosity>5) cout << "call bemtool func POT_Generator<MA,SL_POT ..." << endl;
            }
            break;    
        default: ffassert(0);
        
    }
}

double ff_htoolEta=10., ff_htoolEpsilon=1e-3;
long ff_htoolMinclustersize=10, ff_htoolMaxblocksize=1000000, ff_htoolMintargetdepth=0, ff_htoolMinsourcedepth=0;


template<class K>
class HMatrixVirt {
public:
    virtual const std::map<std::string, std::string>& get_infos() const = 0;
    virtual void mvprod_global(const K* const in, K* const out,const int& mu=1) const = 0;
    virtual int nb_rows() const = 0;
    virtual int nb_cols() const = 0;
    virtual void cluster_to_target_permutation(const K* const in, K* const out) const = 0;
    virtual void source_to_cluster_permutation(const K* const in, K* const out) const = 0;
    virtual MPI_Comm get_comm() const = 0;
    virtual int get_rankworld() const = 0;
    virtual int get_sizeworld() const = 0;
    virtual int get_local_size() const = 0;
    virtual int get_local_offset() const = 0;
    virtual const std::vector<htool::SubMatrix<K>*>& get_MyNearFieldMats() const = 0;
    virtual const htool::LowRankMatrix<K>& get_MyFarFieldMats(int i) const = 0;
    virtual int get_MyFarFieldMats_size() const = 0;
    virtual const std::vector<htool::SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const = 0;
    virtual htool::Matrix<K> get_local_dense() const = 0;
    virtual int get_permt(int) const = 0;
    virtual int get_perms(int) const = 0;
    virtual char get_symmetry_type() const = 0;
    virtual void set_compression(std::shared_ptr<htool::VirtualLowRankGenerator<K>> compressor) = 0;

    // Getters/setters for parameters
    virtual double get_epsilon() const                            = 0;
    virtual double get_eta() const                                = 0;
    virtual int get_minsourcedepth() const                        = 0;
    virtual int get_mintargetdepth() const                        = 0;
    virtual int get_maxblocksize() const                          = 0;
    virtual void set_epsilon(double epsilon0)                     = 0;
    virtual void set_eta(double eta0)                             = 0;
    virtual void set_minsourcedepth(unsigned int minsourcedepth0) = 0;
    virtual void set_mintargetdepth(unsigned int mintargetdepth0) = 0;
    virtual void set_maxblocksize(unsigned int maxblocksize)      = 0;

    // Build
    virtual void build(htool::VirtualGenerator<K> &mat, double* xt, double* xs) = 0;
    virtual void build(htool::VirtualGenerator<K> &mat, double* xt)             = 0;

    virtual ~HMatrixVirt() {};
};


template<class K>
class HMatrixImpl : public HMatrixVirt<K> {
public:
    htool::HMatrix<K> H;
public:
    HMatrixImpl(std::shared_ptr<htool::VirtualCluster> t, std::shared_ptr<htool::VirtualCluster> s, double epsilon=1e-6 ,double eta=10, char symmetry='N', char UPLO='N',const int& reqrank=-1, MPI_Comm comm=MPI_COMM_WORLD) : H(t,s,epsilon,eta,symmetry,UPLO,reqrank,comm){}
    const std::map<std::string, std::string>& get_infos() const {return H.get_infos();}
    void mvprod_global(const K* const in, K* const out,const int& mu=1) const {return H.mvprod_global_to_global(in,out,mu);}
    int nb_rows() const { return H.nb_rows();}
    int nb_cols() const { return H.nb_cols();}
    void cluster_to_target_permutation(const K* const in, K* const out) const {return H.cluster_to_target_permutation(in,out);}
    void source_to_cluster_permutation(const K* const in, K* const out) const {return H.source_to_cluster_permutation(in,out);}
    MPI_Comm get_comm() const {return H.get_comm();}
    int get_rankworld() const {return H.get_rankworld();}
    int get_sizeworld() const {return H.get_sizeworld();}
    int get_local_size() const {return H.get_local_size();}
    int get_local_offset() const {return H.get_local_offset();}
    const std::vector<htool::SubMatrix<K>*>& get_MyNearFieldMats() const {return H.get_MyNearFieldMats();}
    const htool::LowRankMatrix<K>& get_MyFarFieldMats(int i) const {return *(H.get_MyFarFieldMats()[i]);}
    int get_MyFarFieldMats_size() const {return H.get_MyFarFieldMats().size();}
    const std::vector<htool::SubMatrix<K>*>& get_MyStrictlyDiagNearFieldMats() const {return H.get_MyStrictlyDiagNearFieldMats();}
    htool::Matrix<K> get_local_dense() const {return H.get_local_dense();}
    int get_permt(int i) const {return H.get_permt(i);}
    int get_perms(int i) const {return H.get_perms(i);}
    char get_symmetry_type() const {return H.get_symmetry_type();}
    void set_compression(std::shared_ptr<htool::VirtualLowRankGenerator<K>> compressor) {H.set_compression(compressor);};

    // Getters/setters for parameters
    double get_epsilon() const  {return H.get_epsilon();}                          
    double get_eta() const    {return H.get_epsilon();}                            
    int get_minsourcedepth() const {return H.get_minsourcedepth();}              
    int get_mintargetdepth() const {return H.get_mintargetdepth();}            
    int get_maxblocksize() const   {return H.get_maxblocksize();}              
    void set_epsilon(double epsilon0) {H.set_epsilon(epsilon0);}             
    void set_eta(double eta0)        {H.set_eta(eta0);}                        
    void set_minsourcedepth(unsigned int minsourcedepth0) {H.set_minsourcedepth(minsourcedepth0);}   
    void set_mintargetdepth(unsigned int mintargetdepth0) {H.set_mintargetdepth(mintargetdepth0);}   
    void set_maxblocksize(unsigned int maxblocksize0)   {H.set_maxblocksize(maxblocksize0);}      

    // Build
    void build(htool::VirtualGenerator<K> &mat, const std::vector<htool::R3> &xt, const std::vector<double> &rt, const std::vector<int> &tabt, const std::vector<double> &gt, const std::vector<htool::R3> &xs, const std::vector<double> &rs, const std::vector<int> &tabs, const std::vector<double> &gs) {H.build(mat,xt,tabt,xs,tabs);}


    void build(htool::VirtualGenerator<K> &mat, double* xt, double* xs) {H.build(mat,xt,xs);}
    void build(htool::VirtualGenerator<K> &mat, double* xt) {H.build(mat,xt);}                       
};


struct Data_Bem_Solver
: public Data_Sparse_Solver {
    double eta;
    int minclustersize,maxblocksize,mintargetdepth,minsourcedepth;
    string compressor;
       
    Data_Bem_Solver()
       : Data_Sparse_Solver(),
       eta(ff_htoolEta),
       minclustersize(ff_htoolMinclustersize),
       maxblocksize(ff_htoolMaxblocksize),
       mintargetdepth(ff_htoolMintargetdepth),
       minsourcedepth(ff_htoolMinsourcedepth),
       compressor("partialACA")
    
        {epsilon=ff_htoolEpsilon;}
     
    template<class R>
       void Init_sym_positive_var();
    
    private:
        Data_Bem_Solver(const Data_Bem_Solver& ); // pas de copie
        
};

template<class R>
void Data_Bem_Solver::Init_sym_positive_var()
{
    
    Data_Sparse_Solver::Init_sym_positive_var<R>(-1);
    
}

template<class R>
inline void SetEnd_Data_Bem_Solver(Stack stack,Data_Bem_Solver & ds,Expression const *nargs ,int n_name_param)
{
    
    {
        bool unset_eps=true;
        ds.initmat=true;
        ds.factorize=0;
        int kk = n_name_param-(NB_NAME_PARM_MAT+NB_NAME_PARM_HMAT)-1;
        if (nargs[++kk]) ds.initmat= ! GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.solver= * GetAny<string*>((*nargs[kk])(stack));
        ds.Init_sym_positive_var<R>();//  set def value of sym and posi
        if (nargs[++kk]) ds.epsilon= GetAny<double>((*nargs[kk])(stack)),unset_eps=false;
        if (nargs[++kk])
        {// modif FH fev 2010 ...
            const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[kk]);
            if(op)
            {
                ds.precon = op->Find("(",ArrayOfaType(atype<KN<R>* >(),false)); // strange bug in g++ is R become a double
                ffassert(ds.precon);
            } // add miss
        }
        
        if (nargs[++kk]) ds.NbSpace= GetAny<long>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.tgv= GetAny<double>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.factorize= GetAny<long>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.strategy = GetAny<long>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.tol_pivot = GetAny<double>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.tol_pivot_sym = GetAny<double>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.itmax = GetAny<long>((*nargs[kk])(stack)); //  frev 2007 OK
        if (nargs[++kk]) ds.data_filename = *GetAny<string*>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.lparams = GetAny<KN_<long> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.dparams = GetAny<KN_<double> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.smap = GetAny<MyMap<String,String> *>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.perm_r = GetAny<KN_<long> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.perm_c = GetAny<KN_<long> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.scale_r = GetAny<KN_<double> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.scale_c = GetAny<KN_<double> >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.sparams = *GetAny<string*>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.commworld = GetAny<pcommworld>((*nargs[kk])(stack));
#ifdef VDATASPARSESOLVER
        if (nargs[++kk]) ds.master = GetAny<long>((*nargs[kk])(stack));
#else
        ++kk;
#endif
        if (nargs[++kk]) ds.rinfo = GetAny<KN<double>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.info = GetAny<KN<long>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.kerneln = GetAny< KNM<double>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.kernelt = GetAny< KNM<double>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.kerneldim = GetAny<long * >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.verb = GetAny<long  >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.x0 = GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.veps= GetAny<double*>((*nargs[kk])(stack));
        if( unset_eps && ds.veps) ds.epsilon = *ds.veps;//  if veps  and no def value  => veps def value of epsilon.
        if (nargs[++kk]) ds.rightprecon= GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.sym= GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.positive= GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk])  { ds.getnbiter= GetAny<long*>((*nargs[kk])(stack));
                   if( ds.getnbiter) *ds.getnbiter=-1; //undef
               }
        if(ds.solver == "") { // SET DEFAULT SOLVER TO HRE ...
            if( ds.sym && ds.positive ) ds.solver=*def_solver_sym_dp;
            else if( ds.sym ) ds.solver=*def_solver_sym;
            else  ds.solver=*def_solver;
            if(mpirank==0 && verbosity>4) cout << "  **Warning: set default solver to " << ds.solver << endl;
        }
        if (nargs[++kk]) ds.eta = GetAny<double>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.minclustersize = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.maxblocksize = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.mintargetdepth = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.minsourcedepth = GetAny<int>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.compressor = *GetAny<string*>((*nargs[kk])(stack));
        ffassert(++kk == n_name_param);
    }   
}

template<class K>
class HMatrixVirt;

template<class Type, class K>
KNM<K>* HMatrixVirtToDense(HMatrixVirt<K> **Hmat)
{ 
    ffassert(Hmat && *Hmat);
    HMatrixVirt<K>& H = **Hmat;
    if(std::is_same<Type, KNM<K>>::value) {
        htool::Matrix<K> mdense = H.get_local_dense();
        KNM<K>* M = new KNM<K>( (long) H.nb_cols(), (long) H.nb_rows() );
        for (int i=0; i< H.nb_rows(); i++)
            for (int j=0; j< H.nb_cols(); j++)
                (*M)(i,j) = 0;
        for (int i=0; i< mdense.nb_rows(); i++)
            for (int j=0; j< mdense.nb_cols(); j++)
                (*M)(H.get_permt(i+H.get_local_offset()),H.get_perms(j)) = mdense(i,j);
        return M;
    }
    else {
        std::cerr << " Error of value type HMatrixVirt<K> and KNM<K> " << std::endl;
        ffassert(0);
    }
}


template <class R>
void builHmat(HMatrixVirt<R>** Hmat, htool::VirtualGenerator<R>* generatorP,const Data_Bem_Solver& data,vector<double> &p1,vector<double> &p2,MPI_Comm comm)  {
    std::shared_ptr<htool::Cluster<htool::PCARegularClustering>> t =std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
    std::shared_ptr<htool::Cluster<htool::PCARegularClustering>> s =std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
    t->set_minclustersize(data.minclustersize);
    s->set_minclustersize(data.minclustersize);
    t->build(generatorP->nb_rows(),p2.data(),2,comm);
    s->build(generatorP->nb_cols(),p1.data(),2,comm);

    *Hmat = new HMatrixImpl<R>(t,s,data.epsilon,data.eta,data.sym?'S':'N',data.sym?'U':'N',-1,comm);
    std::shared_ptr<htool::VirtualLowRankGenerator<R>> LowRankGenerator = nullptr;
    if (data.compressor=="" || data.compressor == "partialACA")
        LowRankGenerator = std::make_shared<htool::partialACA<R>>();
    else if (data.compressor == "fullACA")
        LowRankGenerator = std::make_shared<htool::fullACA<R>>();
    else if (data.compressor == "SVD")
        LowRankGenerator = std::make_shared<htool::SVD<R>>();
    else {
        cerr << "Error: unknown htool compressor \""+data.compressor+"\"" << endl;
        ffassert(0);
    }

    (*Hmat)->set_compression(LowRankGenerator);
    (*Hmat)->set_maxblocksize(data.maxblocksize);
    (*Hmat)->set_mintargetdepth(data.mintargetdepth);
    (*Hmat)->set_minsourcedepth(data.minsourcedepth);
    // TODO set options
    (*Hmat)->build(*generatorP,p2.data(),p1.data());

}

template<class Matrix, class R, typename std::enable_if< std::is_same< Matrix, HMatrixVirt<R>* >::value >::type* = nullptr>
void Assembly(Matrix* A, htool::VirtualGenerator<R>* generator, const Data_Bem_Solver& data,vector<double> &p1,vector<double> &p2,MPI_Comm comm,int dim,bool = false) {
    builHmat<R>(A,generator,data,p1,p2,comm);
}


template<class R, class MMesh, class FESpace1, class FESpace2> 
void creationHMatrixtoBEMForm(const FESpace1 * Uh, const FESpace2 * Vh, const int & VFBEM, 
                             const list<C_F0> & largs, Stack stack, const Data_Bem_Solver &ds, HMatrixVirt<R>** Hmat){


    typedef typename FESpace1::Mesh SMesh;
    typedef typename FESpace2::Mesh TMesh;
    typedef typename SMesh::RdHat SRdHat;
    typedef typename TMesh::RdHat TRdHat;

    typedef typename std::conditional<SMesh::RdHat::d==1,bemtool::Mesh1D,bemtool::Mesh2D>::type MeshBemtool;
    typedef typename std::conditional<SMesh::RdHat::d==1,bemtool::P0_1D,bemtool::P0_2D>::type P0;
    typedef typename std::conditional<SMesh::RdHat::d==1,bemtool::P1_1D,bemtool::P1_2D>::type P1;
    typedef typename std::conditional<SMesh::RdHat::d==1,bemtool::P2_1D,bemtool::P2_2D>::type P2;

    // size of the matrix
    int n=Uh->NbOfDF;
    int m=Vh->NbOfDF;

    /*
    int VFBEM = typeVFBEM(largs,stack);
    if (mpirank == 0 && verbosity>5)
        cout << "test VFBEM type (1 kernel / 2 potential) "  << VFBEM << endl;
    */

    // compression infogfg
    
    MPI_Comm comm = ds.commworld ? *(MPI_Comm*)ds.commworld : MPI_COMM_WORLD;
    
    // source/target meshes
    const SMesh & ThU =Uh->Th; // line
    const TMesh & ThV =Vh->Th; // colunm
    bool samemesh = (void*)&Uh->Th == (void*)&Vh->Th;  // same Fem2D::Mesh     +++ pot or kernel
 
    /*
    if (VFBEM==1)
        ffassert (samemesh);
        if(init)
            *Hmat =0;
        *Hmat =0;
    if(*Hmat)
        delete *Hmat;
    *Hmat =0;
    */

    bemtool::Geometry node; MeshBemtool mesh;
    Mesh2Bemtool(ThU, node, mesh);


    vector<double> p1(3*n);
    vector<double> p2(3*m);
    Fem2D::R3 pp;
    bemtool::R3 p;
    SRdHat pbs;
    TRdHat pbt;
    pbs[0] = 1./(SRdHat::d+1);
    pbs[1] = 1./(SRdHat::d+1);
    if (SRdHat::d == 2) pbs[2] = 1./(SRdHat::d+1);
    pbt[0] = 1./(TRdHat::d+1);
    pbt[1] = 1./(TRdHat::d+1);
    if (TRdHat::d == 2) pbt[2] = 1./(TRdHat::d+1);

    int Snbv = Uh->TFE[0]->ndfonVertex;
    int Snbe = Uh->TFE[0]->ndfonEdge;
    int Snbt = Uh->TFE[0]->ndfonFace;
    bool SP0 = SRdHat::d == 1 ? (Snbv == 0) && (Snbe == 1) && (Snbt == 0) : (Snbv == 0) && (Snbe == 0) && (Snbt == 1);
    bool SP1 = (Snbv == 1) && (Snbe == 0) && (Snbt == 0);
    bool SP2 = (Snbv == 1) && (Snbe == 1) && (Snbt == 0);
    bool SRT0 = (SRdHat::d == 2) && (Snbv == 0) && (Snbe == 1) && (Snbt == 0);

    if(mpirank == 0){
        cout << "Vh->TFE[0]->N=" << Vh->TFE[0]->N << endl;
        //cout << "Vh->TFE[0].N=" << Vh->TFE[0].N << endl;
        cout << "Vh->TFE.N()=" << Vh->TFE.N() << endl;
        cout << "Vh->MaxNbNodePerElement=" << Vh->MaxNbNodePerElement << endl;
        cout << "SRT0=" << SRT0 << endl;
    }

    if (SP2) {
        bemtool::Dof<P2> dof(mesh,true);
        for (int i=0; i<n; i++) {
            const std::vector<bemtool::N2>& jj = dof.ToElt(i);
            p = dof(jj[0][0])[jj[0][1]];
            p1[3*i+0] = p[0];
            p1[3*i+1] = p[1];
            p1[3*i+2] = p[2];
        }
    }
    else if (SRT0) {
        bemtool::Dof<bemtool::RT0_2D> dof(mesh);
        for (int i=0; i<n; i++) {
            const std::vector<bemtool::N2>& jj = dof.ToElt(i);
            p = dof(jj[0][0])[jj[0][1]];
            p1[3*i+0] = p[0];
            p1[3*i+1] = p[1];
            p1[3*i+2] = p[2];
        }
    }
    else
    for (int i=0; i<n; i++) {
        if (SP1)
            pp = ThU.vertices[i];
        else if (SP0)
            pp = ThU[i](pbs);
        else {
            if (mpirank == 0) std::cerr << "ff-BemTool error: only P0, P1 and P2 discretizations are available for now." << std::endl;
            ffassert(0);
        }
        p1[3*i+0] = pp.x;
        p1[3*i+1] = pp.y;
        p1[3*i+2] = pp.z;
    }

    std::shared_ptr<htool::Cluster<htool::PCARegularClustering>> t =std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
    std::shared_ptr<htool::Cluster<htool::PCARegularClustering>> s =std::make_shared<htool::Cluster<htool::PCARegularClustering>>();
    t->set_minclustersize(ds.minclustersize);
    s->set_minclustersize(ds.minclustersize);
    t->build(n,p1.data(),2,comm);

    if(!samemesh) {
        if( Vh->TFE[0]->N == 1){
            // case the targer FE is scalar
            for (int i=0; i<m; i++) {
                if (Vh->MaxNbNodePerElement == TRdHat::d + 1)
                    pp = ThV.vertices[i];
                else if (Vh->MaxNbNodePerElement == 1)
                    pp = ThV[i](pbt);
                else {
                    if (mpirank == 0) std::cerr << "ff-BemTool error: only P0 and P1 FEspaces are available for reconstructions." << std::endl;
                    ffassert(0);
                }
                p2[3*i+0] = pp.x;
                p2[3*i+1] = pp.y;
                p2[3*i+2] = pp.z;
            }
        }
        else{
            // hack for Maxwell case to have one Hmatrix to avoid one Hmatrix by direction
            ffassert(SRT0 && SRdHat::d == 2 && VFBEM==2);

            // Dans un espace verctoriel, [P1,P1,P1] pour les targets, on a:
            //        m correspond au nombre de dof du FEM space
            // Or dans ce cas, on veut que m = mesh_Target.nv 
            // 
            // ==> on n'a pas besoin de resize les points p2
            int nnn= Vh->TFE[0]->N; // the size of the vector FESpace. For [P1,P1,P1], nnn=3;
            if(verbosity){
                cout << "nnn=" << nnn << endl;
                cout << "p2.size= " << p2.size() << endl; 
                cout << "p2.resize :"<< nnn*3*m << " "<< 3*m << endl;
                cout << "m=" << m << endl; 
                cout << "n=" << n << endl;
            }
            int mDofScalar = m/nnn; // computation of the dof of one component 

            for (int i=0; i<mDofScalar; i++) {
                if (Vh->MaxNbNodePerElement == TRdHat::d + 1)
                    pp = ThV.vertices[i];
                else if (Vh->MaxNbNodePerElement == 1)
                    pp = ThV[i](pbt);
                else {
                    if (mpirank == 0) std::cerr << "ff-BemTool error: only P0 and P1 FEspaces are available for reconstructions." << std::endl;
                    ffassert(0);
                }
                //if( i%2 ==0) cout << i << ",  pp.x =" << pp.x << endl;
                for(int iii=0; iii<nnn; iii++){
                    ffassert( nnn*3*i+3*iii+2 < nnn*3*m );
                    //if(i== m-1) cout << "pp.x =" << pp.x << endl;
                    p2[nnn*3*i+3*iii+0] = pp.x;
                    p2[nnn*3*i+3*iii+1] = pp.y;
                    p2[nnn*3*i+3*iii+2] = pp.z;
                }
            }
        }
        //cout << "call s->build( (Vh->TFE[0]->N)*m,p2.data(),2,comm);" << endl;
        s->build( m,p2.data(),2,comm);  
    }
    else{
        p2=p1;
        s=t;
    }


    // creation of the generator for htool and creation the matrix

    if (VFBEM==1) {
        // info kernel
        cout << "call getBemKernel(stack,largs);" << endl;
        pair<BemKernel*,Complex> kernel = getBemKernel(stack,largs);
        BemKernel *Ker = kernel.first;
        Complex alpha = kernel.second;
        
        if(mpirank == 0 && verbosity >5) {
            int nk=-1;
            iscombinedKernel(Ker) ? nk=2 : nk=1;
            for(int i=0;i<nk;i++)
                cout << " kernel info... i: " << i << " typeKernel: " << Ker->typeKernel[i] << " wave number: " << Ker->wavenum[i]  << " coeffcombi: " << Ker->coeffcombi[i] <<endl;
        }
        htool::VirtualGenerator<R>* generator;

        if(mpirank == 0 && verbosity>5)
            std::cout << "Creating dof" << std::endl;

        if (SP0) {
            bemtool::Dof<P0> dof(mesh);
            ff_BIO_Generator<R,P0,SMesh>(generator,Ker,dof,alpha);
        }
        else if (SP1) {
            bemtool::Dof<P1> dof(mesh,true);
            ff_BIO_Generator<R,P1,SMesh>(generator,Ker,dof,alpha);
        }
        else if (SP2) {
            bemtool::Dof<P2> dof(mesh,true);
            ff_BIO_Generator<R,P2,SMesh>(generator,Ker,dof,alpha);
        }
        else if (SRT0 && SRdHat::d == 2) {
            bemtool::Dof<bemtool::RT0_2D> dof(mesh);
            ff_BIO_Generator_Maxwell<R>(generator,Ker,dof,alpha);
        }
        else
            ffassert(0);

        // build the Hmat
        *Hmat = new HMatrixImpl<R>(t,s,ds.epsilon,ds.eta,ds.sym?'S':'N',ds.sym?'U':'N',-1,comm);
        std::shared_ptr<htool::VirtualLowRankGenerator<R>> compressor = nullptr;
        if ( ds.compressor == "" || ds.compressor == "partialACA")
            compressor = std::make_shared<htool::partialACA<R>>();
        else if (ds.compressor == "fullACA")
            compressor = std::make_shared<htool::fullACA<R>>();
        else if (ds.compressor == "SVD")
            compressor = std::make_shared<htool::SVD<R>>();
        else {
            cerr << "Error: unknown htool compressor \""+ds.compressor+"\"" << endl;
            ffassert(0);
        }

        (*Hmat)->set_compression(compressor);
        (*Hmat)->set_maxblocksize(ds.maxblocksize);
        (*Hmat)->set_mintargetdepth(ds.mintargetdepth);
        (*Hmat)->set_minsourcedepth(ds.minsourcedepth);
        (*Hmat)->build(*generator,p1.data());


        delete generator;
    }
    else if (VFBEM==2) {
        cout << "VFBEM==2"<< endl;
        BemPotential *Pot = getBemPotential(stack,largs);
        bemtool::Geometry node_output;
        if (Vh->MaxNbNodePerElement == TRdHat::d + 1)
            Mesh2Bemtool(ThV,node_output);
        else if (Vh->MaxNbNodePerElement == 1) {
            for (int i=0; i<m; i++) {
                pp = ThV[i](pbt);
                p[0]=pp.x; p[1]=pp.y; p[2]=pp.z;
                node_output.setnodes(p);
            }
        }
        else
            ffassert(0);

        htool::VirtualGenerator<R>* generator;
        if (SP0) {
            bemtool::Dof<P0> dof(mesh);
            ff_POT_Generator<R,P0,MeshBemtool,SMesh>(generator,Pot,dof,mesh,node_output);
        }
        else if (SP1) {
            bemtool::Dof<P1> dof(mesh,true);
            ff_POT_Generator<R,P1,MeshBemtool,SMesh>(generator,Pot,dof,mesh,node_output);
        }
        else if (SP2) {
            bemtool::Dof<P2> dof(mesh,true);
            ff_POT_Generator<R,P2,MeshBemtool,SMesh>(generator,Pot,dof,mesh,node_output);
        }
        else if (SRT0 && SRdHat::d == 2) {
            bemtool::Dof<bemtool::RT0_2D> dof(mesh);
            ff_POT_Generator_Maxwell<R,bemtool::RT0_2D>(generator,Pot,dof,mesh,node_output);
        }
        else
            ffassert(0);
        Assembly(Hmat,generator,ds,p1,p2,comm,MMesh::RdHat::d+1);
        delete generator;
    }
    //return Hmat;
}



/*
template void creationHMatrixtoBEMForm<Complex, MeshL, FESpaceL, FESpaceL>(const FESpaceL * Uh, const FESpaceL * Vh, const int & VFBEM, 
                             const std::list<C_F0> & largs, Stack stack, const Data_Bem_Solver &ds, HMatrixVirt<Complex> **Hmat);

template void creationHMatrixtoBEMForm<Complex, MeshS, FESpaceS, FESpaceS>(const FESpaceS * Uh, const FESpaceS * Vh, const int & VFBEM, 
                             const std::list<C_F0> & largs, Stack stack, const Data_Bem_Solver &ds, HMatrixVirt<Complex> **Hmat);
*/

template<> void creationHMatrixtoBEMForm<double, MeshS, FESpaceS, FESpaceS>(const FESpaceS * Uh, const FESpaceS * Vh, const int & VFBEM, 
                             const std::list<C_F0> & largs, Stack stack, const Data_Bem_Solver &ds, HMatrixVirt<double> **Hmat){
                                 cerr << "we can't use bemtool with Real type." << endl;
                                 ffassert(0);
                             }


template<> void creationHMatrixtoBEMForm<double, MeshL, FESpaceL, FESpaceL>(const FESpaceL * Uh, const FESpaceL * Vh, const int & VFBEM, 
                             const std::list<C_F0> & largs, Stack stack, const Data_Bem_Solver &ds, HMatrixVirt<double> **Hmat){
                                 cerr << "we can't use bemtool with Real type." << endl;
                                 ffassert(0);
                             }


#endif