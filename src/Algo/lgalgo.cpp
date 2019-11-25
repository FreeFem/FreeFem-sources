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

#include <iostream>
#include <cfloat>
using namespace std;
#include "ff++.hpp"
#include "NRJ.hpp"
#include "RosenBrock.hpp"
#include "LineSearch.hpp"
#include "CubicLS.hpp"
#include "BrentLS.hpp"
#include "Optima.hpp"
#include "BFGS.hpp"
//#include  "CG.hpp"
#include "NewtonRaphson.hpp"

// template<class R>
extern Block *currentblock;

typedef double R;

class OptimAlgo : public OneOperator {
 public:
  typedef KN< R > Kn;
  typedef KN_< R > Kn_;
  typedef R REAL;
  typedef Param< REAL > PARAM;
  typedef KN< REAL > VECT;
  typedef KNM< REAL > MAT;
  typedef Matrice_Creuse< R >::VMat VMAT;

  const int cas;

  class E_LCG : public E_F0mps {
   public:
    const int cas;
    static basicAC_F0::name_and_type name_param[];
    static const int n_name_param = 3;
    Expression nargs[n_name_param];
    Expression X;
    C_F0 inittheparam, theparam, closetheparam;
    Expression J, dJ, hJ;
    long arg(int i, Stack stack, long a) const {
      return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
    }
    R arg(int i, Stack stack, R a) const { return nargs[i] ? GetAny< R >((*nargs[i])(stack)) : a; }

    class lgNRJ : public tNRJ< PARAM, VECT, VMAT, REAL > {
     private:
      Stack stack;
      Expression J, dJ, hJ, theparame;
      VECT *gg;

     protected:
      void setparam(const Param< R > &x) {
        KN< double > *p = GetAny< KN< double > * >((*theparame)(stack));
        ffassert(p->N( ) == x.N( ));
        *p = x;
      }

     public:
      lgNRJ(Stack s, int n, Expression t, Expression JJ, Expression dJJ, Expression hJJ)
        : tNRJ< PARAM, VECT, VMAT, REAL >(n), stack(s), J(JJ), dJ(dJJ), hJ(hJJ), theparame(t),
          gg(0) {
        if (dJ) gg = new VECT(n);
      }

      ~lgNRJ( ) {
        if (gg) delete gg;
      }

      REAL Val(const Param< R > &x) {
        setparam(x);
        assert(J);
        R ret = GetAny< R >((*J)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        return ret;
      }

      VECT *Gradient(const Param< R > &x) {
        setparam(x);
        if (dJ) {
          *gg = GetAny< Kn_ >((*dJ)(stack));
          WhereStackOfPtr2Free(stack)->clean( );
        }
        return gg;    // dJ ? GetAny<Kn*>( (*dJ)(stack)) :0;}
      }

      VMAT *Hessian(const Param< R > &x) {
        setparam(x);
        if (!hJ) return 0;
        Matrice_Creuse< R > *M = GetAny< Matrice_Creuse< R > * >((*hJ)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        assert(M && M->A);
        return (VMAT *)M->A;
      }
    };

    E_LCG(const basicAC_F0 &args, int cc) : cas(cc) {
      int nbj = args.size( ) - 1;
      Block::open(currentblock);    // make a new block to
      X = to< Kn * >(args[nbj]);
      C_F0 X_n(args[nbj], "n");
      // the expression to init the theparam of all
      inittheparam =
        currentblock->NewVar< LocalVariable >("the parameter", atype< KN< R > * >( ), X_n);
      theparam = currentblock->Find("the parameter");    // the expression for the parameter
      args.SetNameParam(n_name_param, name_param, nargs);
      const Polymorphic *opJ = 0;
      const Polymorphic *opdJ = 0;
      const Polymorphic *ophJ = 0;
      if (nbj > 0) {
        opJ = dynamic_cast< const Polymorphic * >(args[0].LeftValue( ));
        assert(opJ);
      }
      if (nbj > 1) {
        opdJ = dynamic_cast< const Polymorphic * >(args[1].LeftValue( ));
        assert(opdJ);
      }
      if (nbj > 2) {
        ophJ = dynamic_cast< const Polymorphic * >(args[2].LeftValue( ));
        assert(ophJ);
      }
      J = dJ = hJ = 0;

      J = to< R >(C_F0(opJ, "(", theparam));
      if (opdJ)
        dJ = to< Kn_ >(C_F0(
          opdJ, "(", theparam));    // May2019 (bug?) Modif FH 17102005 (to verify) to<Kn*> ->to<Kn>
      if (ophJ) hJ = to< Matrice_Creuse< R > * >(C_F0(ophJ, "(", theparam));

      closetheparam = C_F0((Expression)Block::snewclose(currentblock),
                           atype< void >( ));    // The cleanning block expression
      /*
      if (nargs[2]) {
        const Polymorphic *op = dynamic_cast<const Polymorphic *>(nargs[2]);
        assert(op);
        C = op->Find("(", ArrayOfaType(atype<Kn* >(), false));
      } else C = 0;
      */
    }

    virtual AnyType operator( )(Stack stack) const {
      WhereStackOfPtr2Free(stack) = new StackOfPtr2Free(stack);
      typedef LineSearch< PARAM, VECT, MAT, VMAT, REAL > LS;
      typedef CubicLineSearch< LS > CLS;
      KN< double > *delta = 0;
      R tol = arg(0, stack, 1E-6);
      int itMax = arg(1, stack, 100L);
      int itMaxLine = arg(2, stack, 100L);
      bool verbose = verbosity > 3;

      try {
        Kn &x = *GetAny< Kn * >((*X)(stack));

        const int n = x.N( );
        GetAny< KN< double > * >(inittheparam.eval(stack));    // do allocation

        Param< R > param(x);
        lgNRJ nrj1(stack, n, theparam, J, dJ, hJ);
        if (!dJ) delta = new KN< double >(n);
        CLS ls(&nrj1, itMaxLine, delta);
        REAL initialNrj = nrj1.getVal(param);
        Optima< LS > *opt = 0;
        if (cas == 1)
          opt = new BFGS< LS >((CLS *)&ls, itMax, tol, verbose);
        else if (cas == 2)
          opt = new Newt< LS >((CLS *)&ls, itMax, tol, verbose);
        else
          ErrorExec("lgalgo: Case no available. Sorry, internal bug", cas);

        param = opt->optimizer(param);

        REAL finalNrj = nrj1.getVal(param);

        if (verbosity) cout << endl << "*** RESULTS SUMMARY ***" << endl;

        if (verbosity > 1) {
          cout << "  The number of call to NRJ : " << nrj1.Appel_Val( ) << endl;
          cout << "  The number of call to gradient : " << nrj1.Appel_Grad( ) << endl;
          cout << "  The number of call to hessian : " << nrj1.Appel_Hess( ) << endl;
        }
        if (verbosity > 2) {
          cout << "Normalized residue : ";
          display(opt->allResidue( ));
        }
        if (verbosity) {
          cout << "  Initial NRJ value : " << initialNrj << endl;
          cout << "  Final NRJ value : " << finalNrj << endl;
        }
        delete opt;

        // cout << "The final optimized parameters : " << param;
        x = param;
      } catch (...) {
        closetheparam.eval(stack);                // clean memory
        WhereStackOfPtr2Free(stack)->clean( );    // FH mars 2005
        if (delta) delete delta;

        throw;
      }

      if (delta) delete delta;
      closetheparam.eval(stack);                // clean memory
      WhereStackOfPtr2Free(stack)->clean( );    // FH mars 2005

      return 0L;    // SetAny<long>(0);  Modif FH  july 2005
    }

    operator aType( ) const { return atype< long >( ); }
  };

  E_F0 *code(const basicAC_F0 &args) const { return new E_LCG(args, cas); }

  OptimAlgo(int c)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< KN< R > * >( )), cas(c) {}

  OptimAlgo(int c, int cc)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< KN< R > * >( )),
      cas(c) {}

  OptimAlgo(int c, int cc, int ccc)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< Polymorphic * >( ), atype< KN< R > * >( )),
      cas(c) {}
};

basicAC_F0::name_and_type OptimAlgo::E_LCG::name_param[] = {
  {"eps", &typeid(double)}, {"nbiter", &typeid(long)}, {"nbiterline", &typeid(long)}};
void init_algo( );

void init_algo( ) {
  Global.Add("BFGS", "(", new OptimAlgo(1, 1));         // j + dJ
  Global.Add("Newton", "(", new OptimAlgo(2, 2, 2));    // j + dJ

  // Start msh3 debug
  // extern void init_msh3();
  // init_msh3();
  // End mesh3 debug
}

//#include "InitFunct.hpp"
// static addingInitFunct TheaddingInitFunct(-9, init_algo);
