#include <algorithm>


template<class TypeIndex=int,class TypeScalar=double>
struct CGMatVirt {
public:
    typedef TypeIndex I;
    typedef TypeScalar R;
    
    I n,m;
    virtual  R * addmatmul(R *x,R *Ax) const =0;
    virtual ~CGMatVirt() {}

    virtual R * matmul(R *x,R *Ax) const
    {
        std::fill(Ax,Ax+n, 0.);
        return addmatmul(x,Ax);
    }
    virtual void  SetInitWithBC(R*rhs,R *x) const {}// do nothing by default ..
    CGMatVirt(int nn,int mm=-1) : n(nn),m(mm<0 ?nn:mm) {}
    virtual int * pwcl() const {return  0;} // array know if  node with BC (TGV) 
};

template<class TypeIndex=int,class TypeScalar=double>
int ConjugueGradient(CGMatVirt<TypeIndex,TypeScalar> &A, // fonction et pointeur data pour A
		     CGMatVirt<TypeIndex,TypeScalar>  &C, // fonction et pointeur data pour C
		     TypeScalar * b, // second membre
		     TypeScalar * x, // solution qui contient une initialisation
		     int nbitermax,
                     double eps,
                     int niveauimpression)
;

template<typename K,typename Z>
bool fgmres(CGMatVirt<Z,K> &A, // fonction et pointeur data pour A
            CGMatVirt<Z,K> &C,int leftC,
            K *y,
            K *x,
            double tol,
            int maxits,
            int restart=50,
            int verbo=3,
            int *perm=0 )
;

template<class I,class K>
K * myscopy(I n,const K *x,K *y);
template<class I,class K>
K * myscal(I n,K a,K *x);
template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const CGMatVirt<TypeIndex,TypeScalar> *A,TypeScalar *x, TypeScalar *Ax) { return A->matmul(x,Ax);}
template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const CGMatVirt<TypeIndex,TypeScalar> &A,TypeScalar *x, TypeScalar *Ax) { return A.matmul(x,Ax);}
template<class I,class K> K * mysaxpy(I n,K a,const K *x,K *y);

template<class Z=int,class R=double>
struct CGMatVirtId : public CGMatVirt<Z,R> {
    CGMatVirtId(Z nn): CGMatVirt<Z,R> (nn,nn) {}
    R * matmul(R *x,R *Ax) const { myscopy(this->n,x,Ax);  return Ax;}
    R * addmatmul(R *x,R *Ax) const { mysaxpy(this->n,R(1.),x,Ax); return Ax;}
    
};

