/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <complex>

using namespace std;
//#include "ff++.hpp"

#include "error.hpp"
#include "arpackff.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"

#include "HashMatrix.hpp"

#include "SparseLinearSolver.hpp"

#include "Mesh3dn.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
//#include "ffstack.hpp"

extern Block *currentblock;

typedef double R;
static bool dddd = false;
template<class K>
void Show(int ido, KN_<K> w, const char *cmm) {
  cout << cmm << ido << " max = " << w.max() << " - min = " << w.min() << " - sum = " << w.sum() << endl;
}

template<class R>
class FuncMat : public RNM_VirtualMatrix<R> {
public:
  // typedef double K;
  typedef KN<R> Kn;
  typedef KN_<R> Kn_;

  Stack stack;

  mutable Kn x;
  C_F0 c_x;
  Kn *b;
  Expression mat1, mat;
  typedef typename RNM_VirtualMatrix<R>::plusAx plusAx;

  FuncMat (int n, Stack stk, const OneOperator *op, const OneOperator *op1, Kn *bb=0)
    : RNM_VirtualMatrix<R>(n),
    stack(stk),
    x(n), c_x(CPValue(x)), b(bb),
    mat1(op1 ? CastTo<Kn_>(C_F0(op1->code(basicAC_F0_wa(c_x)), (aType)*op1)) : 0),
    mat(op ? CastTo<Kn_>(C_F0(op->code(basicAC_F0_wa(c_x)), (aType)*op)) : 0)
    {}

  FuncMat () { if (mat1 != mat) delete mat; delete mat1; delete c_x.LeftValue(); }

  void addMatMul (const Kn_ &xx, Kn_ &Ax) const {
    // cout << " addMatMul " << xx.N() << " " << Ax.N() << " " << cas << endl;
    ffassert(xx.N() == Ax.N() && mat);
    x =xx;
    Ax += GetAny<Kn_>((*mat)(stack));
    if (b && &Ax != b) Ax += *b; // Ax - b => add b (not in cas of init. b c.a.d &Ax == b)
    WhereStackOfPtr2Free(stack)->clean(); // add dec 2008
  }

  void Solve (KN_<R> &Ax, const KN_<R> & xx) const {
    ffassert(xx.N() == Ax.N() && mat1);
    x = xx;
    Ax = GetAny<Kn_>((*mat1)(stack));
    if (b && &Ax != b) Ax += *b; // Ax -b => add b (not in cas of init. b c.a.d &Ax == b)
    WhereStackOfPtr2Free(stack)->clean(); // add dec 2008
  }

  plusAx operator * (const Kn &x) const { return plusAx(this,x); }
  virtual bool ChecknbLine (int n) const { return true; }
  virtual bool ChecknbColumn (int m) const { return true; }
  bool WithSolver () const { return mat1; }
};

void ToWhich (Stack stack, char *which, Expression ew) {
  if (!ew) strcpy(which, "LM");
  else {
    string *sw = GetAny<string*>((*ew)(stack));
    if (sw && sw->length() == 2)
      strcpy(which, sw->c_str());
    int ok = 0;
    if (strcmp(which, "LM") == 0) ok = 1;
    else if (strcmp(which, "LA") == 0) ok = 1;
    else if (strcmp(which, "SM") == 0) ok = 1;
    else if (strcmp(which, "SA") == 0) ok = 1;
    else if (strcmp(which, "SR") == 0) ok = 1;
    else if (strcmp(which, "LR") == 0) ok = 1;
    else if (strcmp(which, "BE") == 0) ok = 1;
    ffassert(ok);
  }
}

void CheckMode (int mode, bool A, bool A1, bool B, bool B1) {
  if (mode == 1) { ffassert(A); }
  if (mode == 2) { ffassert(A && B && B1 ); }
  if (mode == 3) { ffassert(A1 && B); }
}

template<class K,class Mat>
void DoIdoAction (int ido, char bmat, int mode, KN_<K> &xx, KN_<K> &yy, KN_<K> &zz, KN_<K> &work, Mat &OP1, Mat &B) {
  if ((mode == 3) && (bmat == 'G')) // shift-invert mode, generalized problem
    switch (ido) {
      case -1: //y <--- OP*x = inv[A-SIGMA*M]*M*x
        if(dddd) Show(ido, xx, " <- (xx) ");
        OP1.Solve(yy, work=B*xx);
        if(dddd) Show(ido, yy, " -> (yy) ");
        break;
      case 1:// y <-- OP*x = inv[A-sigma*M]*z
        if(dddd) Show(ido, zz, " <- (zz) ");
        OP1.Solve(yy, zz);
        if(dddd) Show(ido, yy, " -> (yy) ");
        break;
      case 2: // y <--- M*x
        if(dddd) Show(ido, xx, " <- (xx) ");
        yy = B*xx;
        if(dddd) Show(ido, yy, " -> (yy) ");
        break;
      default:
        ffassert(0);
    }
  else if ((mode == 2) && (bmat == 'G')) { // shift-invert mode, generalized problem
    switch (ido) {
      case -1: // y <--- OP*x = inv[M]*A*x
      case 1:
        if (mode == 1) // M = Id //TODO there is a non-sens here: mode == 2 && mode == 1
          yy = OP1*xx;
        else
          B.Solve(yy, work=OP1*xx);
        break;
      case 2: // y <--- M*x. // not use mode = 1
        yy = B*xx;
        break;
      default :
        ffassert(0);
    }
  }
  else if ((mode == 1) && (bmat == 'I')) // direct mode, simple problem
    switch (ido) {
      case 1:
        B.Solve(yy, work=OP1*xx);
        break;
      default :
        ffassert(0);
    }
  else if ((mode == 3) && (bmat == 'I')) // shift-invert mode, simple problem
    switch (ido) {
      case -1:
      case 1:
        OP1.Solve(yy, work=B*xx);
        break;
      default :
        ffassert(0);
    }
}

/*void DoIdoAction(int ido,int mode,KN_<K> &xx,KN_<K> &yy,KN_<K> &zz,KN_<K> &work,Mat &A,Mat &B)
{
    if(mode>2) // inverse mode
        switch (ido) {
            case -1: //y <--- OP*x = inv[A-SIGMA*M]*M*x
                if(dddd) Show(ido,xx," <- (xx) ");
                A.Solve(yy,work=B*xx);
                if(dddd) Show(ido,yy," -> (yy) ");
                break;
            case 1:// y <-- OP*x = inv[A-sigma*M]*z
                if(dddd) Show(ido,zz," <- (zz) ");
                A.Solve(yy,zz);
                if(dddd) Show(ido,yy," -> (yy) ");
                break;
            case 2: // y <--- M*x
                if(dddd) Show(ido,xx," <- (xx) ");
                yy = B*xx;
                if(dddd) Show(ido,yy," -> (yy) ");
                break;
            default :
                ffassert(0);
        }
    else // direct mode
        switch (ido)
    {
        case -1: // y <--- OP*x = inv[M]*A*x
        case 1:
            if(mode== 1) // M = Id
                yy=A*xx;
            else
                B.Solve(yy,work=A*xx);
            break;
        case 2: // y <--- M*x. // not use mode = 1
            yy = B*xx;
            break;
        default :
            ffassert(0);
    }
}*/

template<class K>
const OneOperator *ToCode(Expression e) {
  const OneOperator *code = 0;
  if (e != 0) {
    const Polymorphic *op = dynamic_cast<const Polymorphic*>(e);
    ffassert(op);
    code = op->Find("(", ArrayOfaType(atype<KN<K>* >(), false));
      if(code == 0)
      {
          cerr << " Erreur we are in Eigen "<< typeid(K).name() << " "<< endl;
          cout << " bad type of function : Maybe we need to use  complexEigenValue or EigenValue" << endl;
          ffassert(code);
      }
  }
  return code;
}

class EigenValue : public OneOperator {
public:
  typedef R K;
  typedef KN<K> Kn;
  typedef KN_<K> Kn_;
  const int cas;
  class E_EV : public E_F0mps {
  public:
    const int cas;

    static basicAC_F0::name_and_type name_param[];
    static const int n_name_param = 18;
    Expression nargs[n_name_param];
    Expression expOP1, expB, expn;
    const OneOperator *codeOP1, *codeB, *codeOP, *codeB1;
    template<class T>
    T arg (int i, Stack stack, const T &a) const { return nargs[i] ? GetAny<T>((*nargs[i])(stack)) : a; }
    E_EV (const basicAC_F0 &args, int cc) :
      cas(cc), expOP1(0), expB(0), expn(0), codeOP1(0), codeB(0), codeOP(0), codeB1(0) {
      // OP1 = (A-sigma*B)
      // int nbj= args.size()-1;
      args.SetNameParam(n_name_param, name_param, nargs);
      if (cas == 1 || cas == 2) {
        expn = to< long>(args[0]);
        const Polymorphic *op = dynamic_cast<const Polymorphic*>(args[1].LeftValue());
        ffassert(op);
        codeOP1 = op->Find("(", ArrayOfaType(atype<Kn* >(), false));
        ffassert(codeOP1);
        if (cas == 2) {
          const Polymorphic *op = dynamic_cast<const Polymorphic*>(args[2].LeftValue());
          ffassert(op);
          codeB = op->Find("(", ArrayOfaType(atype<Kn* >(), false));
          ffassert(codeB);
        }
      } else if (cas == 3) {
        expn = to< long>(args[0]);
      } else {
        expOP1 = to<Matrice_Creuse<K> *>(args[0]);
        expB = to<Matrice_Creuse<K> *>(args[1]);
      }
      ffassert((codeOP1 == 0) || (nargs[15] == 0)); // Double def
      ffassert((codeB == 0) || (nargs[14] == 0)); // Double def
      codeOP = ToCode<double>(nargs[13]);
      if (!codeB) codeB = ToCode<double>(nargs[14]);
      if (!codeOP1) codeOP1 = ToCode<double>(nargs[15]);
      codeB1 = ToCode<double>(nargs[16]);
    }

    AnyType operator () (Stack stack) const;
    operator aType () const { return atype<long>(); }
  };

  E_F0 * code (const basicAC_F0 & args) const { return new E_EV(args,cas); }

  EigenValue () :
    OneOperator(atype<long>(),
    atype<Matrice_Creuse<K> *>(),
    atype<Matrice_Creuse<K> *>()),
    cas(0) {}

  EigenValue (int, int) :
    OneOperator(atype<long>(), atype<long>(),
    atype<Polymorphic *>(),
    atype<Polymorphic *>()),
    cas(2) {}

  EigenValue (int) :
    OneOperator(atype<long>(), atype<long>(),
    atype<Polymorphic *>()),
    cas(1) {}

  EigenValue(int, int, int) : // All data are given through A=, B= A1=, B1=
    OneOperator(atype<long>(), atype<long>()),
    cas(3) {}
};

class EigenValueC : public OneOperator {
public:
  typedef Complex K;
  typedef double R;
  typedef KN<K> Kn;
  typedef KN_<K> Kn_;
  const int cas;
  class E_EV : public E_F0mps {
  public:
    const int cas;

    static basicAC_F0::name_and_type name_param[];
    static const int n_name_param = 16;
    Expression nargs[n_name_param];
    Expression expOP1, expB, expn;
    const OneOperator *codeOP1, *codeB, *codeOP, *codeB1;
    template<class T>
    T arg (int i, Stack stack, const T &a) const { return nargs[i] ? GetAny<T>((*nargs[i])(stack)) : a; }
    E_EV (const basicAC_F0 &args, int cc) :
      cas(cc), expOP1(0), expB(0), expn(0), codeOP1(0), codeB(0), codeOP(0), codeB1(0) {
      // OP1 = (A-sigma*B)
      // int nbj = args.size()-1;
      args.SetNameParam(n_name_param, name_param, nargs);
      if (cas == 1 || cas == 2) {
        expn = to<long>(args[0]);
        const Polymorphic *op = dynamic_cast<const Polymorphic*>(args[1].LeftValue());
        ffassert(op);
        codeOP1 = op->Find("(", ArrayOfaType(atype<Kn* >(), false));
        ffassert(codeOP1);
        if (cas == 2) {
          const Polymorphic *op = dynamic_cast<const Polymorphic*>(args[2].LeftValue());
          ffassert(op);
          codeB = op->Find("(", ArrayOfaType(atype<Kn* >(), false));
          ffassert(codeB);
        }
      } else if (cas == 3) {
        expn = to<long>(args[0]);
      } else {
        expOP1 = to<Matrice_Creuse<K> *>(args[0]);
        expB = to<Matrice_Creuse<K> *>(args[1]);
      }
      ffassert((codeOP1 == 0) || (nargs[13] == 0)); // Double def
      ffassert((codeB == 0) || (nargs[12] == 0)); // Double def
      codeOP = ToCode<K>(nargs[11]);
      if (!codeB) codeB = ToCode<K>(nargs[12]);
      if (!codeOP1) codeOP1 = ToCode<K>(nargs[13]);
      codeB1 = ToCode<K>(nargs[14]);
    }

    AnyType operator () (Stack stack) const;
    operator aType () const { return atype<long>(); }
  };

  E_F0 *code (const basicAC_F0 &args) const { return new E_EV(args,cas); }

  EigenValueC () :
    OneOperator(atype<long>(),
    atype<Matrice_Creuse<K> *>(),
    atype<Matrice_Creuse<K> *>()),
    cas(0) {}

  EigenValueC (int, int ) :
    OneOperator(atype<long>(), atype<long>(),
    atype<Polymorphic *>(),
    atype<Polymorphic *>()),
    cas(2) {}

  EigenValueC (int) :
    OneOperator(atype<long>(), atype<long>(),
    atype<Polymorphic *>()),
    cas(1) {}

  EigenValueC (int, int, int) : // All data are given through A=, B= A1=, B1=
    OneOperator(atype<long>(), atype<long>()),
    cas(3) {}
};



basicAC_F0::name_and_type EigenValue::E_EV::name_param[] = {
  {"tol", &typeid(double)},
  {"nev", &typeid(long)},
  {"sym", &typeid(bool)},
  {"sigma", &typeid(double)},
  {"value", &typeid(KN<double> *)},
  {"vector", &typeid(FEbaseArrayKn<double> *)}, // pferarray
  {"ncv", &typeid(long)}, // the number of Arnoldi vectors generated
  {"maxit", &typeid(long)}, // the maximum number of Arnoldi iterations
  {"ivalue", &typeid(KN<double> *)},
  {"rawvector", &typeid(KNM<double> *)},
  {"resid", &typeid(KN<double> *)},
  {"driver", &typeid(long)}, // 11th
  {"which", &typeid(string*)}, // 12th
  {"A", &typeid(Polymorphic*)},
  {"B", &typeid(Polymorphic*)},
  {"A1", &typeid(Polymorphic*)},
  {"B1", &typeid(Polymorphic*)},
  {"mode", &typeid(long)} // 17th
};

basicAC_F0::name_and_type EigenValueC::E_EV::name_param[] = {
  {"tol", &typeid(double)},
  {"nev", &typeid(long)},
  {"sigma", &typeid(Complex)},
  {"value", &typeid(KN<Complex> *)},
  {"vector", &typeid(FEbaseArrayKn<Complex> *)}, // pfecarray
  {"ncv", &typeid(long)}, // the number of Arnoldi vectors generated
  {"maxit", &typeid(long)}, // the maximum number of Arnoldi iterations
  {"rawvector", &typeid(KNM<Complex> *)},
  {"resid", &typeid(KN<Complex> *)},
  {"driver", &typeid(long)},// 9th
  {"which", &typeid(string*)}, // 10th
  {"A", &typeid(Polymorphic*)},
  {"B", &typeid(Polymorphic*)},
  {"A1", &typeid(Polymorphic*)},
  {"B1", &typeid(Polymorphic*)},
  {"mode", &typeid(long)} // 15th
};

AnyType EigenValue::E_EV::operator () (Stack stack) const {
  dddd = (verbosity >= 200);
  double tol = 1e-6;
  long nconv = 0;
  long nbev = 1;
  bool sym = false;
  long ncv = 0; // the number of Arnoldi vectors generated
  long maxit = 0; // the maximum number of Arnoldi iterations
  double sigma = 0;
  int driver = 4; //
  int mymode = -1; // unset
  KN<double> *evalue = 0;
  KN<double> *evaluei = 0;
  KN<double> *resid = 0;
  KNM<double> *rawvector = 0;
  //double ws, vs; // for debugging FH ++++++++++
  // pferarray evector2;
  FEbaseArrayKn<double> *evector = 0; // change may 2009
  tol = arg<double>(0, stack, 0);
  nbev = arg<long>(1, stack, 10);
  sym = arg<bool>(2, stack, false);
  sigma = arg<double>(3, stack, 0.0);
  evalue = arg<KN<double>*>(4, stack, 0);
  evector = arg<FEbaseArrayKn<double>*>(5, stack, 0);
  ncv = arg<long>(6, stack, 0);
  maxit = arg<long>(7, stack, 0);
  evaluei = arg<KN<double>*>(8, stack, 0);
  rawvector = arg<KNM<double>*>(9, stack, 0);
  resid = arg<KN<double>*>(10, stack, 0);
  mymode = arg<long>(17, stack, -1);
  driver = arg<long>(11, stack, mymode==2 ? 3 : 4);

  char which[3];
  ToWhich (stack, which, nargs[12]);

  // evector=evector2.first;

  Matrice_Creuse<K> *pOP1 = 0, *pB = 0;
  if (expOP1)
    pOP1 = GetAny<Matrice_Creuse<K> *>((*expOP1)(stack));
  if (expB)
    pB = GetAny<Matrice_Creuse<K> *>((*expB)(stack));
  double *residptr = resid ? (double*) *resid : 0;
  long n = 0;
  if (pOP1)
    n = pOP1->A->n;
  else {
    ffassert(expn);
    n = GetAny<long>((*expn)(stack));
  }

  if (evalue) nbev = Max((long)evalue->N(), nbev);

  const RNM_VirtualMatrix<K> *ptOP1 = 0, *ptB = 0, *ptOP = 0, *ptB1 = 0;

  if (pOP1)
    ptOP1 = &(const RNM_VirtualMatrix<K>&)(pOP1->A);
  if (pB)
    ptB = &(const RNM_VirtualMatrix<K>&)(pB->A);
  FuncMat<K> *pcOP = 0, *pcB = 0;
  if (codeOP1 || codeOP)
    ptOP1 = pcOP = new FuncMat<K>(n, stack, codeOP, codeOP1);
  if (codeB || codeB1)
    ptB = pcB = new FuncMat<K>(n, stack, codeB, codeB1);

  MatriceIdentite<K> Id(n);

  if (!ptB) ptB = &Id;
  const RNM_VirtualMatrix<K> &OP1 = *ptOP1, &B = *ptB;
  if (sym) {
    nbev = min(n-1, nbev);
    if (!ncv) ncv = min(nbev*2+1, n);
  } else {
    nbev = min(nbev, n-2);
    if (!ncv) ncv = nbev*2 + 1;
  }
  ncv = max(nbev+2, ncv);
  ncv = min(ncv, n);

  if (!maxit) maxit = 100*nbev;

  const char *serr[10];
  int err = 0;
  if (!(nbev < n))
    serr[err++] = "  Number of eigenvalues of OP to be computed nev <= n ";

  if ((driver < 1 || driver > 5) && sym)
    serr[err++] = "  the driver = 1, 2, 3, 4, 5 ";
  if ((driver < 1 || driver > 4) && !sym)
  serr[err++] = "  the driver = 1, 2, 3, 4 ";
  // 2 <= NCV-NEV and NCV <= N

  if (!(ncv <= n) && sym)
    serr[err++] = "   ( ncv <= n) (symmetric case) ";
  if (!((ncv <= n) && 2 <= (ncv-nbev )) && !sym )
    serr[err++] = "   ( ncv <= n) 2 <= (ncv-nev ) (non-symmetric case) ";

  if (n != OP1.M)
    serr[err++] = " the first matrix in EigneValue is not square.";
  if (n != B.N)
    serr[err++] = "Sorry the row's number of the secand matrix in EigneValue is wrong.";
  if (n != B.M)
    serr[err++] = "Sorry the column's number of the secand matrix in EigneValue is wrong.";

  if (verbosity && mpirank == 0) {
    if(sym)
      cout << "Real symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
    else
      cout << "Real non symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
  }

  if (verbosity > 9 || err)
    cout << "    n " << n << ", nev "<< nbev << ", tol =" << tol << ", maxit =" << maxit
      << ", ncv = " <<ncv << ", driver = " << driver << " which = " << which << endl;
  if (err) {
    cerr << " list of the error " << endl;
    for (int i = 0; i < err; ++i)
      cerr << "\n-  " << serr[i] << endl;
    ExecError("Fatal Error in EigenValue (real)");
  }

  KN<K> work(n);
  int dmode[] = {-1, 1, 3, 2, 3};

  if (sym) {
    int ido = 0;
    //char bmat= mode == 1 ? 'I' : 'G';
    char bmat = (driver<3)? 'I' : 'G';
    int mode = dmode[driver];
    if (mymode > 0 && (mymode != mode)) {
      cerr << " Error mode == " << mymode << " == " << mode << " driver =" << driver << endl;
      ExecError("Wrong mod in Arpack sym");
    }

    // char which[3]= "LM";	// larger value
    // if(mode >2) which[0] ='S'; // smaller value
    int ishift = 1; // Auto Shift true by default
    int iparam[12] = {0, ishift, 0, (int)maxit, 1, (int)nconv, 0, mode, 0, 0, 0, 0};
    int ipntr[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    KN<double> workd(3*n+1);
    int lworkl = ncv*(ncv+9);
    KN<double> workl(lworkl+1);
    KN<double> vp(ncv*n+1);

    int info = (residptr !=0);
    KN<double> vresid(residptr ? 1 : n);
    if (!residptr) residptr = &vresid[0];

    while (1) {
      saupp(ido, bmat, n, which, nbev, tol, residptr, ncv, vp, n,
        iparam, ipntr, workd, workl, lworkl, info);
      if (verbosity > 99)
        cout << "    saupp ido: " << ido << " info : " << info << endl;
      if (info < 0) { cerr << "  -- err arpack info = " << info << endl; }
      sauppError(info);

      if (ido == 99) break;

      KN_<double> xx(&workd[ipntr[1]], n);
      KN_<double> yy(&workd[ipntr[2]], n);
      KN_<double> zz(&workd[ipntr[3]], n);
      DoIdoAction(ido, bmat, mode, xx, yy, zz, work, OP1, B);
    }

    nconv = iparam[5];
    if (nconv == 0) cerr << "  -- Strange: no eigens values ??? " << endl;
    // Finding eigenvalues and eigenvectors.
    if (nconv) {
      int newdim = nbev;
      if (nconv > nbev) {
        cerr << "WARNING: nconv(saupp) > nbev: " << nconv << " > " << nbev << endl;
        newdim = nconv;
      }
      KN<double> evr(newdim);
      KNM<double> Z(n, newdim);
      int ldz = n;
      char HowMny ='A';
      int rvec = 1;

      seupp(rvec, HowMny, evr, Z, ldz, sigma, bmat, n,
        which, nbev, tol, residptr, ncv, vp, n, iparam,ipntr, workd, workl, lworkl, info);

      if (verbosity > 5) {
        cout << "Dimension of the system            : " << n << endl;
        cout << "Number of 'requested' eigenvalues  : " << nbev << endl;
        cout << "Number of 'converged' eigenvalues  : " << nconv << endl;
        cout << "Number of Arnoldi vectors generated: " << ncv << endl;
        cout << "Number of iterations taken         : " << iparam[3] << endl;
        cout << endl;
        nconv = min(nconv, nbev);
        //if (prob.EigenvaluesFound()) {
        cout << "Eigenvalues:" << endl;// Add FH feb 2016
        for (int i = 0; i < nconv; i++) {
          cout << "  lambda[" << (i+1) << "]: " << evr(i) << endl;
          KN_<K> vi(Z(':', i)) ;
          if (verbosity > 99)
            cout <<" Eigenvector: :" << vi << endl;
        }
        cout << endl;
      }

      if (evalue) {
        KN<double> & ev(*evalue);
        int m = Min(nconv, ev.N());
        for (int i = 0; i < m; i++)
          ev[i] = evr(i);
      }
      if (rawvector) {
        int m = Min(nconv, rawvector->M());
        ffassert(rawvector->N() == n);
        for (int i = 0; i < m; i++) {
          KN_<K> vi(Z(':', i)) ;
          // cout << " ------ EV--raw " << vi.min() << " " << vi.max() << endl;
          (*rawvector)(':', i) = vi;
        }
      }

      if (evector) {
        FEbaseArrayKn<K> & ev(*evector);
        int m = Min(nconv, (long)ev.N);
        for (int i = 0; i < m; i++) {
          //KN<K> & xx= *(ev(i));
          //if(xx.pVh->NbDoF != n)
          //ExecError("Wrong Type size of FEspace to store the eigen vector ");
          // if (xx.pVh != pOP1->pUh)
          // ExecError("Wrong Type of FEspace to store the eigen vector ");
          //xx.Vh = pOP1->Uh;
          KN_<K> vi(Z(':', i));
          ev.set(i, vi);
        }
      }
    }
  } else {
    // non symmetric case ,
    // Finding an Arnoldi basis.
    //int mode = 3; // Shift invert
    int ido = 0;
    //char bmat='G';
    char bmat = (driver<3) ? 'I' : 'G';
    int mode = dmode[driver];

    if (mymode > 0 && (mymode != mode)) {
      cerr << " Error nonsym mode == " << mymode << " == " << mode << " driver =" << driver << endl;
      ExecError("Wrong mod in Arpack non-sym");
    }

    // char which[] = "LM";
    int ishift = 1; // Auto Shift true by default
    int iparam[12] = {0, ishift, 0, (int)maxit, 1, (int)nconv, 0, mode, 0, 0, 0, 0};
    int ipntr[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    KN<double> workd(3*n+1);
    int lworkl = 3*ncv*(ncv+2);
    KN<double> workl(lworkl+1);
    KN<double> vp(ncv*n+1);
    int info = (residptr != 0);
    KN<double> vresid(residptr ? 1 : n);
    if (!residptr) residptr = &vresid[0];
    if (verbosity > 9)
      cout << " n " << n << " nbev "<< nbev << " tol =" << tol << " maxit =" << maxit << " ncv = " << ncv << endl;

    while (1) {
      naupp(ido, bmat, n, which, nbev, tol, residptr, ncv, vp, n,
        iparam, ipntr, workd, workl, lworkl, info);
      if (verbosity > 99) { cout << "    naupp ido: " << ido << " info : " << info << endl; }
      if (info < 0) { cerr << "  -- err arpack info = " << info << endl; }
      sauppError(info);
      if (ido == 99) break;
      KN_<double> xx(&workd[ipntr[1]], n);
      KN_<double> yy(&workd[ipntr[2]], n);
      KN_<double> zz(&workd[ipntr[3]], n);
      DoIdoAction(ido, bmat, mode, xx, yy, zz, work, OP1, B);
    }

    nconv = iparam[5];
    if (nconv) {
      int newdim = nbev;
      if (nconv > nbev) {
        cerr << "WARNING: nconv(naupp) > nbev: " << nconv << " > " << nbev << endl;
        newdim = nconv;
      }
      KN<double> evr(newdim+1), evi(newdim+1);
      KNM<double> Z(n, newdim+1);
      KN<double> workev(3*ncv);
      int ldz = n;
      char HowMny = 'A';
      int rvec = 1;
      double sigmai = 0;
      neupp(rvec, HowMny, evr,evi, Z , ldz, sigma, sigmai, workev,
        bmat, n, which, nbev, tol, residptr, ncv,
        vp, n, iparam, ipntr, workd, workl, lworkl, info);

      if (verbosity) {
        cout << "Real non-symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
        cout << "mode " << mode << " sigma=" << sigma << endl << endl;

        cout << "Dimension of the system            : " << n << endl;
        cout << "Number of 'requested' eigenvalues  : " << nbev << endl;
        cout << "Number of 'converged' eigenvalues  : " << nconv << endl;
        cout << "Number of Arnoldi vectors generated: " << ncv << endl;
        cout << "Number of iterations taken         : " << iparam[3] << endl;
        cout << endl;
        nconv = min(nconv, nbev);// Add FH feb 2016

        cout << "Eigenvalues:" << endl;
        for (int i = 0; i < nconv; i++) {
          cout << "  lambda[" << (i+1) << "]: ";
          ios::fmtflags oldflag = cout.flags();
          cout.setf(ios::showpos);
          cout << evr(i) << evi(i) << "i" << endl;
          cout.flags(oldflag); // restore flags
          KN_<double> vi(Z(':', i));
          if (verbosity > 99) {
            cout << " Eigenvector: " << vi << endl;
            cout << endl;
          }
        }
      }

      if (evalue) {
        KN<double> &ev(*evalue);
        int m = Min(nconv, ev.N());
        for (int i = 0; i < m; i++)
          ev[i] = evr(i);
      }

      if (evaluei) {
        KN<double> &ev(*evaluei);
        int m = Min(nconv, ev.N());
        for (int i = 0; i < m; i++)
          ev[i] = evi(i);
      }

      if (rawvector) {
        //K* rawev(prob.RawEigenvectors());
        int m = Min(nconv, rawvector->M());
        ffassert(rawvector->N() == n);
        for (int i = 0; i < m; i++) {
          //int k = i;
          KN_<K> vi(Z(':', i));
          if (verbosity > 2 && mpirank ==0)
            cout << "   ------ EV--raw " << vi.min() << " " << vi.max() << endl;
          (*rawvector)(':', i) = vi;
        }
      }

      if (evector) {
        // K*  rawev(prob.RawEigenvectors());
        // rawev + n*k is
        //  iev = prob.EigenvalueImag(k)
        //  iev==0 => the eigen vector
        //  iev> 0 => real
        //      start real :  rawev + n*k
        //      start imag :  ramev +n*(k+1)
        //  iev < 0 =>  complex
        //      start real :  rawev + n*(k-1)
        //      -start imag :  ramev +n*(k)
        FEbaseArrayKn<K> &ev(*evector);
        int m = Min(nconv, (long)ev.N);
        for (int i = 0; i < m; i++) {
          // K ev_i=
          //prob.EigenvalueImag(i);
          //KN<K> & xx= *(ev[i]);
          // if (xx.pVh != pOP1->pUh)
          //    ExecError("Wrong Type of FEspace to store the eigen vector ");
          // xx.Vh = pOP1->Uh;
          // int  k=(ev_i < 0) ? i-1 : i;
          //int k=i;
          KN_<K> vi(Z(':', i));//rawev+n*k,n) ;
          //xx= new KN<K>(vi);
          ev.set(i, vi);//new KN<K>(vi));
        }
      }
    }
  }
  if(pcOP) delete pcOP;
  if(pcB) delete pcB;

  return (long)nconv;
}

AnyType EigenValueC::E_EV::operator () (Stack stack) const {
  dddd = (verbosity >= 200);
  double tol = 1e-6;
  long nconv = 0;
  long nbev = 1;
  long ncv = 0; // the number of Arnoldi vectors generated
  long maxit = 0; // the maximum number of Arnoldi iterations
  long driver = 4;
  K sigma = 0;
  KN<K> *evalue = 0;
  KN<K> *resid = 0;
  KNM<K> *rawvector = 0;

  // pfecarray evector2;
  FEbaseArrayKn<Complex> *evector = 0;
  tol = arg<double>(0, stack, 0);
  nbev = arg<long>(1, stack, 2);
  sigma = arg<K>(2, stack, 0.0);
  evalue = arg<KN<K>*>(3, stack, 0);
  // evector2 =arg<pfecarray>(4,stack,make_pair<pfecbasearray,int>(0,0));
  evector = arg<FEbaseArrayKn<Complex>*>(4, stack, 0);
  ncv = arg<long>(5, stack, 0);
  maxit = arg<long>(6, stack, 0);
  rawvector = arg<KNM<K>*>(7, stack, 0);
  resid = arg<KN<K>*>(8, stack, 0);
  int mymode = arg<long>(15, stack, -1);
  driver = arg<long>(9, stack, mymode==2 ? 3 : 4);
  K *residptr = resid ? (K*)*resid : 0;
  char which[3];
  ToWhich(stack, which, nargs[10]);

  //evector=evector2.first;
  ffassert(driver > 0 && driver < 5);
  // Matrice_Creuse<K> *pOP1 = GetAny<Matrice_Creuse<K> *>((*expOP1)(stack));
  // Matrice_Creuse<K> *pB = GetAny<Matrice_Creuse<K> *>((*expB)(stack));

  if (evalue) nbev = Max((long)evalue->N()-2, nbev);
  if (!maxit) maxit = 100*nbev;

  //const MatriceCreuse<K> & OP1 = pOP1->A;
  //const MatriceCreuse<K> & B = pB->A;

  //long n=OP1.n;

  Matrice_Creuse<K> *pOP1 = 0, *pB = 0;
  if (expOP1)
    pOP1 = GetAny<Matrice_Creuse<K>*>((*expOP1)(stack));
  if (expB)
    pB = GetAny<Matrice_Creuse<K>*>((*expB)(stack));
  // double * residptr=resid? (double*) *resid : 0;
  // cout << " residptr = " << residptr <<endl;
  long n = 0;
  if (pOP1)
    n = pOP1->A->n;
  else {
    ffassert(expn);
    n = GetAny<long>((*expn)(stack));
  }


  const RNM_VirtualMatrix<K> *ptOP1, *ptB;
  if (pOP1)
    ptOP1 = &(const RNM_VirtualMatrix<K>&)(pOP1->A);
  if (pB)
    ptB = &(const RNM_VirtualMatrix<K>&)(pB->A);
  FuncMat<K> *pcOP1 = 0, *pcB = 0;
  if (codeOP1 || codeOP)
    ptOP1 = pcOP1 = new FuncMat<K>(n, stack, codeOP, codeOP1);
  if (codeB || codeB1)
    ptB = pcB = new FuncMat<K>(n, stack, codeB, codeB1);
  MatriceIdentite<K> Id(n);

  if (!ptB) ptB = &Id;
  const RNM_VirtualMatrix<K> &OP1 = *ptOP1, &B = *ptB;

  nbev = min(nbev, n-2);
  if (!ncv) ncv = nbev*2+1;
  ncv = max(nbev+2, ncv);
  ncv = min(ncv, n);
  const char *serr[10];
  int err = 0;
  if (nbev >= n-1)
    serr[err++] = "  Number of eigenvalues of OP to be computed <= n-2 ";
  if (driver < 1 || driver > 4)
    serr[err++] = "  the driver = 1 ,2 ,4 ";
  // 2 <= NCV-NEV and NCV <= N
  if (!(2 <= nbev && ncv <= n))
    serr[err++] = "   ( 2 <= nbev && nvc <= n) ";
  if (n != OP1.N)
    serr[err++] = " the first matrix in EigneValue is not Hermitien.";
  if (n != B.N)
    serr[err++] = "Sorry the row's number of the secand matrix in EigneValue is wrong.";
  if (n != B.M)
    serr[err++] = "Sorry the column's number of the secand matrix in EigneValue is wrong.";
  if (verbosity)
    cout << "Complex eigenvalue problem: A*x - B*x*lambda" << endl;
  if (verbosity > 9 || err)
    cout << "   n " << n << " nev "<< nbev << " tol = " << tol << " maxit = " << maxit << " ncv = " << ncv
      << " driver = " << driver << " which = " << which << endl;
  if (err) {
    cerr << " list of the error " << endl;
    for (int i = 0; i < err; ++i)
      cerr << "\n-  " << serr[i] << endl;
    ExecError("Fatal Error in EigenValue (complex)");
  }

  KN<K> work(n);
  //  ARrcSymGenEig is a class that requires the user to provide a
  //   way to perform the matrix-vector products w = OP*Bv =
  //   inv(A-sigma*B)*B*v and w = B*v, where sigma is the adopted shift.
  // OP1 = (A-sigma*B)
  // OP = inv(OP)
  /*
  ffassert(0);


  ARrcCompGenEig<R> prob( n, nbev, sigma,"LM",ncv,tol,maxit,residptr);


  // OP = inv[A - sigma*I]


  // Finding an Arnoldi basis.

  while (!prob.ArnoldiBasisFound()) {

  // Calling ARPACK FORTRAN code. Almost all work needed to

  prob.TakeStep();
  */

  // cas non symetrique ,
  // Finding an Arnoldi basis.\

  // int mode=3; // Shift invert			\

  int ido = 0;
  //char bmat='G';
  int dmode[] = {-1, 1, 3, 2, 3};
  char bmat = (driver<3) ? 'I' : 'G';
  int mode = dmode[driver];
  // char which[]="LM";
  int ishift = 1; // Auto Shift true by default
  if (mymode > 0 && (mymode != mode)) {
    cerr << " Error complex case: mode == " << mymode << " == " << mode << " driver = " << driver << endl;
    ExecError("Wrong mod in Arpack complex case:");
  }

  int iparam[12] = {0, ishift, 0, (int)maxit, 1, (int)nconv, 0, (int)mode, 0, 0, 0, 0};
  int ipntr[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  KN<K> workd(3*n+1);
  int lworkl = 3*ncv*ncv + 5*ncv;
  KN<K> workl(lworkl+1);
  KN<K> vp(ncv*n+1);
  int info = (residptr != 0);
  KN<double> rwork(ncv+1);
  KN<K> vresid(residptr ? 1 : n);
  if (!residptr) residptr = &vresid[0];

  while(1) {
    caupp(ido, bmat, n, which, nbev,
      tol, residptr, ncv,
      vp, n, iparam, ipntr, workd, workl,
      lworkl, rwork, info);

    if (verbosity > 99) { cout << "    caupp ido: " << ido << " info : " << info << endl; }
    if (info < 0) { cerr << " -- err arpack info = " << info << endl; }
    sauppError(info);

    if (ido == 99) break;

    KN_<K> xx(&workd[ipntr[1]], n);
    KN_<K> yy(&workd[ipntr[2]], n);
    KN_<K> zz(&workd[ipntr[3]], n);

    DoIdoAction(ido, bmat, mode, xx, yy, zz, work, OP1, B);
  }

  nconv = iparam[5];
  if (nconv) {
    int newdim = nbev;
    if (nconv > nbev) {
      cerr << "WARNING: nconv(caupp) > nbev: " << nconv << " > " << nbev << endl;
      newdim = nconv;
    }
    KN<K> evc(newdim);
    KNM<K> Z(n, newdim);
    KN<K> workev(3*ncv);
    //int ldz=n;
    char HowMny ='A';
    int rvec = 1;
    ceupp(rvec, HowMny, evc, Z, n, sigma, workev, bmat, n, which,
      nbev, tol, residptr, ncv,
      vp, n, iparam, ipntr, workd, workl, lworkl, rwork, info);

    if (verbosity) {
      cout << "Complex eigenvalue problem: A*x - B*x*lambda" << endl;
      cout << "driver = " << driver << " - sigma= " << sigma << endl << endl;

      cout << "Dimension of the system            : " << n << endl;
      cout << "Number of 'requested' eigenvalues  : " << nbev << endl;
      cout << "Number of 'converged' eigenvalues  : " << nconv << endl;
      cout << "Number of Arnoldi vectors generated: " << ncv << endl;
      cout << "Number of iterations taken         : " << iparam[3] << endl;
      cout << endl;
      nconv = min(nconv, nbev);
      cout << "Eigenvalues:" << endl;// Add FH feb 2016
      for (int i = 0; i < nconv; i++) {
        cout << "  lambda[" << (i+1) << "]: " << evc(i) << endl;
        KN_<K> vi(Z(':', i)) ;
        if (verbosity > 99)
          cout <<" Eigenvector: :" << vi << endl;
      }
      cout << endl;
    }

    if (evalue) {
      KN<K> & ev(*evalue);
      int m = Min(nconv, ev.N());
      for (int i = 0; i < m; i++)
        ev[i] = evc(i);
    }
    if (evector) {
      //FEbaseArray<K,v_fes> & ev(*evector);
      FEbaseArrayKn<K> &ev(*evector);
      int m = Min(nconv, (long)ev.N);
      for (int i = 0; i < m; i++) {
        //FEbase<K,v_fes> & xx= **(ev[i]);
        //KN_<K> vi(Z(':',i)) ;
        //xx= new KN<K>(vi);
        ev.set(i, Z(':', i));
      }
    }
    if(rawvector) {
      int m = Min(nconv, rawvector->M());
      ffassert(rawvector->N() == n);
      for (int i = 0; i < m; i++) {
        KN_<K> vi(Z(':', i)) ;
        (*rawvector)(':', i) = vi;
        if (verbosity > 2 && mpirank == 0)
          cout << " raw " << vi.l2() << " == " << (*rawvector)(':', i).l2()<< endl;
      }
    }
  }

  return (long)nconv;
}

#ifndef DYNM_LOAD
void init_eigenvalue () {
  if (verbosity && (mpirank == 0)) cout << "eigenvalue ";
  Global.Add("EigenValue", "(", new EigenValue()); // j + dJ
  Global.Add("EigenValue", "(", new EigenValueC()); // j + dJ
  Global.Add("EigenValue", "(", new EigenValue(1)); // j + dJ
  Global.Add("EigenValue", "(", new EigenValue(1, 1)); // j + dJ
  Global.Add("EigenValue", "(", new EigenValue(1, 1, 1)); // A=, B=
  // sorry no to solve ambiguity => change of name
  Global.Add("complexEigenValue", "(", new EigenValueC(1)); // j + dJ
  Global.Add("complexEigenValue", "(", new EigenValueC(1, 1)); // j + dJ
  Global.Add("complexEigenValue", "(", new EigenValueC(1, 1, 1)); // A=, B= ...
}
#else
class Init {
public:
  Init() {
    if (verbosity && (mpirank == 0)) cout << "eigenvalue ";
    Global.Add("EigenValue2", "(", new EigenValue()); // j + dJ
    Global.Add("EigenValue2", "(", new EigenValueC(1)); // j + dJ
  }
};
Init init;
#endif
