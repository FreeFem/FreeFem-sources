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

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep: bmo.cpp
// *INDENT-ON* //

#include <iostream>
#include <cfloat>
using namespace std;
#include "ff++.hpp"
#include "bmo.hpp"

// template<class R>
extern Block *currentblock;

typedef double R;

class OptimBMO : public OneOperator {
 public:
  typedef KN< R > Kn;
  typedef KN_< R > Kn_;
  typedef R REAL;
  typedef KN< REAL > VECT;
  typedef KNM< REAL > MAT;

  const int cas;

  class E_BMO : public E_F0mps {
   public:
    const int cas;
    static basicAC_F0::name_and_type name_param[];
    static const int n_name_param = 16;
    Expression nargs[n_name_param];
    Expression X;
    C_F0 inittheparam, theparam, closetheparam;
    Expression JJ, dJJ;
    long arg(int i, Stack stack, long a) const {
      return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
    }

    R arg(int i, Stack stack, R a) const { return nargs[i] ? GetAny< R >((*nargs[i])(stack)) : a; }

    string *arg(int i, Stack stack, string *a) const {
      return nargs[i] ? GetAny< string * >((*nargs[i])(stack)) : a;
    }

    void Set_arg(int i, Stack stack, Kn_ v) const {
      if (nargs[i]) {
        v = GetAny< Kn_ >((*nargs[i])(stack));
      }
    }

    class lgBMO : public BijanMO {
     private:
      Stack stack;
      Expression JJ, dJJ, theparame;

     protected:
      void setparam(const KN_< R > &x) {
        KN_< double > *p = GetAny< KN_< double > * >((*theparame)(stack));
        ffassert(p->N( ) == x.N( ));
        *p = x;
      }

     public:
      lgBMO(Stack s, int n, Expression t, Expression J, Expression dJ, int wnbrestart = 1,
            int wnbext1 = 1, int wnbbvp = 5, int wnbgrad = 5, double wepsfd = 1e-5,
            double wrho000 = 100, double wepsloc = 1e-4, double wepsij = 1e-6, int nn100 = 100)

        : BijanMO(n, wnbrestart, wnbext1, wnbbvp, wnbgrad, wepsfd, wrho000, wepsloc, wepsij, nn100),
          stack(s), JJ(J), dJJ(dJ), theparame(t) {}

      ~lgBMO( ) {}

      /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
      /* functional definition */
      /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

      double J(Vect &x) {
        setparam(x);

        double ret = GetAny< R >((*JJ)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        return ret;
      }

      /* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
      /* gradient exact,   no defini => DF */
      /* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

      double *DJ(Vect &x, Vect &fpx) {
        if (!dJJ) {
          return 0;
        }

        setparam(x);
        fpx = GetAny< Kn_ >((*dJJ)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        return fpx;
      }

      void result(Vect &xoptg, Vect &vinit) {}
    };

    E_BMO(const basicAC_F0 &args, int cc) : cas(cc) {
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
      if (nbj > 0) {
        opJ = dynamic_cast< const Polymorphic * >(args[0].LeftValue( ));
        assert(opJ);
      }

      if (nbj > 1) {
        opdJ = dynamic_cast< const Polymorphic * >(args[1].LeftValue( ));
        assert(opdJ);
      }

      JJ = dJJ = 0;

      JJ = to< R >(C_F0(opJ, "(", theparam));
      if (opdJ) {
        dJJ = to< Kn_ >(
          C_F0(opdJ, "(", theparam));    // Modif FH 17102005 (a verifier) to<Kn*> ->to<Kn>
      }

      closetheparam = C_F0((Expression)Block::snewclose(currentblock), atype< void >( ));
    }

    virtual AnyType operator( )(Stack stack) const {
      WhereStackOfPtr2Free(stack) = new StackOfPtr2Free(stack);    // FH mars 2005

      /*
       * basicAC_F0::name_and_type  OptimBMO::E_BMO::name_param[]= {
       * {  "eps", &typeid(double)  },
       * { "nbrestart",&typeid(long) },
       * { "nbbvp",&typeid(long)},
       * { "nbgrad",&typeid(long)},
       * { "epsfd",&typeid(double)},
       * { "epsloc",&typeid(double)},
       * { "epsij",&typeid(double)},
       * { "n100",&typeid(long)} // 7
       * };
       *
       */

      int nbrestart = arg(1, stack, 5L);
      int nbext1 = 5;    // bof bof
      int nbbvp = arg(2, stack, 5L);
      int nbgrad = arg(3, stack, 5L);
      double epsfd = arg(4, stack, 1e-5);
      double rho000 = arg(5, stack, 1e-5);
      double epsloc = arg(6, stack, 1e-4);
      double epsij = arg(7, stack, 1e-6);
      int n100 = arg(8, stack, 100L);
      int diagrand = arg(9, stack, 0L);
      R cmin = arg(9, stack, -1000.);
      R cmax = arg(10, stack, 1000.);
      string *datahist = arg(13, stack, (string *)0);
      string *datachist = arg(14, stack, (string *)0);
      int typealgo = arg(15, stack, 1L);

      try {
        Kn &x = *GetAny< Kn * >((*X)(stack));
        const int n = x.N( );
        Kn xmin(n), xmax(n);
        xmin = cmin;
        xmax = cmax;
        Set_arg(11, stack, xmin);
        Set_arg(12, stack, xmax);

        // Kn * para =
        GetAny< KN< double > * >(inittheparam.eval(stack));    // do allocation

        KN_< R > param(x);
        // cout << nbrestart << " ---- \n";
        lgBMO nrj1(stack, n, theparam, JJ, dJJ, nbrestart, nbext1, nbbvp, nbgrad, epsfd, rho000,
                   epsloc, epsij, n100);
        nrj1.diagrand = diagrand;
        nrj1.debug = verbosity;
        nrj1.typealgo = typealgo;
        nrj1.histpath = datahist;
        nrj1.histcpath = datachist;
        double fopt = nrj1.main(x, xmin, xmax);

        if (verbosity) {
          cout << endl << "*** RESULTS SUMMARY ***" << endl;

          if (verbosity > 1) {
            cout << "  The number of call to  J : " << nrj1.nbeval << endl;
            cout << "  The number of call to dJ : " << nrj1.nbevalp << endl;
          }

          if (verbosity) {
            cout << "  Initial J value : " << nrj1.finit << endl;
            cout << "  Final   J  value : " << fopt << endl;
          }
        }
      } catch (...) {
        closetheparam.eval(stack);                // clean memory
        WhereStackOfPtr2Free(stack)->clean( );    // FH mars 2005
        throw;
      }
      closetheparam.eval(stack);                // clean memory
      WhereStackOfPtr2Free(stack)->clean( );    // FH mars 2005

      return 0L;    // SetAny<long>(0);  Modif FH  july 2005
    }

    operator aType( ) const { return atype< long >( ); }
  };

  E_F0 *code(const basicAC_F0 &args) const { return new E_BMO(args, cas); }

  OptimBMO(int c)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< KN< R > * >( )), cas(c) {}

  OptimBMO(int c, int cc)
    : OneOperator(atype< long >( ), atype< Polymorphic * >( ), atype< Polymorphic * >( ),
                  atype< KN< R > * >( )),
      cas(c) {}
};

basicAC_F0::name_and_type OptimBMO::E_BMO::name_param[] = {
  {"eps", &typeid(double)},         {"nbrestart", &typeid(long)}, {"nbbvp", &typeid(long)},
  {"nbgrad", &typeid(long)},        {"epsfd", &typeid(double)},   {"rho000", &typeid(double)},
  {"epsloc", &typeid(double)},      {"epsij", &typeid(double)},   {"n100", &typeid(long)},    // 8
  {"max", &typeid(double)},                                                                   // 9
  {"min", &typeid(double)},                                                                   // 10
  {"vmax", &typeid(double)},                                                                  // 11
  {"vmin", &typeid(double)},                                                                  // 12
  {"histfile", &typeid(string *)},                                                            // 13
  {"histcfile", &typeid(string *)},                                                           // 14
  {"algo", &typeid(long)}                                                                     // 15
};

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  Global.Add("bmo", "(", new OptimBMO(1));       // j + dJ
  Global.Add("bmo", "(", new OptimBMO(1, 1));    // j + dJ
}

LOADFUNC(Load_Init)
