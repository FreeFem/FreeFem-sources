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
// AUTHORS : Sylvain Auliac
// E-MAIL  : auliac@ann.jussieu.fr

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep: cmaes.cpp
/* clang-format on */

/*
 *      This is a freefem interface of the Hansen's CMA-ES C optimizer.
 *      The class CMAES is a quick C++ interface for the contents of the C files
 *      then follow the real FreeFem++ stuff.
 */

#include <iostream>
using namespace std;
#include "ff++.hpp"

#include "cmaes_interface.h"

class CMAES    // Abstract class, because the fitness function prototype may differ from users to
               // users
{
 public:
  // typedef double (*FFT)(double const *); //Fitness Function Type
  CMAES( ) : pop(0), fitvals(0), evo( ) {}

  CMAES(int d, double *xstart, double *stddev, long seed, int lambda,
        const char *ipf = "initials.par")
    : pop(0), fitvals(0), evo( ) {
    fitvals = init(d, xstart, stddev, seed, lambda, ipf);
    cout << SayHello( ) << endl;
  }

  virtual ~CMAES( ) { exit( ); }

  void resume_distribution(char *filename) { return cmaes_resume_distribution(&evo, filename); }

  double *const *&SamplePopulation( ) { return pop = cmaes_SamplePopulation(&evo); }

  double *UpdateDistribution( ) { return cmaes_UpdateDistribution(&evo, fitvals); }

  const char *TestForTermination( ) const { return cmaes_TestForTermination(&evo); }

  double *const *&ReSampleSingle(int index) { return pop = cmaes_ReSampleSingle(&evo, index); }

  double const *ReSampleSingle_old(double *rgx) { return cmaes_ReSampleSingle_old(&evo, rgx); }

  void UpdateEigensystem(int flgforce) { return cmaes_UpdateEigensystem(&evo, flgforce); }

  virtual void PopEval( ) = 0;

  double axisratio( ) const {
    return cmaes_Get(&evo, "axisratio");
  }    // between lengths of longest and shortest principal axis of the distribution ellipsoid

  int eval( ) const { return floor(cmaes_Get(&evo, "eval")); }    // number of function evaluations

  double fitness( ) const {
    return cmaes_Get(&evo, "fitness");
  }    // recent best function evaluation

  double fbestever( ) const { return cmaes_Get(&evo, "fbestever"); }    // ever best function value

  int generation( ) const { return floor(cmaes_Get(&evo, "generation")); }

  int maxeval( ) const {
    return floor(cmaes_Get(&evo, "maxeval"));
  }    // maximal number of function evaluations

  int maxgen( ) const {
    return floor(cmaes_Get(&evo, "maxgen"));
  }    // maximal number of generations

  double maxaxislength( ) const { return cmaes_Get(&evo, "maxaxislength"); }

  double minaxislength( ) const { return cmaes_Get(&evo, "minaxislength"); }

  double maxstddev( ) const { return cmaes_Get(&evo, "maxstddev"); }

  double minstddev( ) const { return cmaes_Get(&evo, "minstddev"); }

  int dimension( ) const { return floor(cmaes_Get(&evo, "dimension")); }

  int popsize( ) const { return floor(cmaes_Get(&evo, "lambda")); }

  double sigma( ) const { return cmaes_Get(&evo, "sigma"); }

  double *diagC( ) const { return const_cast< double * >(cmaes_GetPtr(&evo, "diag(C)")); }

  double *diagD( ) const { return const_cast< double * >(cmaes_GetPtr(&evo, "diag(D)")); }

  double *stddev( ) const { return const_cast< double * >(cmaes_GetPtr(&evo, "stddev")); }

  double *xbestever( ) const { return const_cast< double * >(cmaes_GetPtr(&evo, "xbestever")); }

  double *xbest( ) const { return const_cast< double * >(cmaes_GetPtr(&evo, "xbest")); }

  double *xmean( ) const { return const_cast< double * >(cmaes_GetPtr(&evo, "xmean")); }

  void ReadSignals(char const *filename) const { cmaes_ReadSignals(&evo, filename); }

  char *SayHello( ) const { return cmaes_SayHello(&evo); }

  void WriteToFile(char const *keyword, char const *output_filename) const {
    cmaes_WriteToFile(&evo, keyword, output_filename);
  }

  virtual double *operator( )( ) {
    // ReadSignals("signals.par");
    while (!TestForTermination( )) {
      SamplePopulation( );
      PopEval( );
      UpdateDistribution( );
      // ReadSignals("signals.par");
    }

    cout << "Stop : " << TestForTermination( ) << endl;
    return xmean( );
  }

  cmaes_t &optimizer( ) { return evo; }

 protected:
  double *const *pop;
  double *fitvals;

 private:
  void exit( ) { cmaes_exit(&evo); }

  double *&init(int dimension, double *xstart, double *stddev, long seed, int lambda,
                const char *input_parameter_filename) {
    return fitvals =
             cmaes_init(&evo, dimension, xstart, stddev, seed, lambda, input_parameter_filename);
  }

  mutable cmaes_t evo;
};

/*
 *      Now comes the FreeFem ++ part :
 */

extern Block *currentblock;

typedef double R;

class OptimCMA_ES : public OneOperator {
 public:
  typedef KN< R > Kn;
  typedef KN_< R > Kn_;
  const int cas;

  class ffcalfunc    // to call the freefem function .. J
  {
   public:
    Stack stack;
    Expression JJ, theparame;

    ffcalfunc(Stack s, Expression JJJ, Expression epar) : stack(s), JJ(JJJ), theparame(epar) {}

    double J(Kn_ x) const {
      KN< double > *p = GetAny< KN< double > * >((*theparame)(stack));
      *p = x;
      double ret = GetAny< R >((*JJ)(stack));
      WhereStackOfPtr2Free(stack)->clean( );
      return ret;
    }
  };

  class CMA_ES : public CMAES {
   public:
    typedef KN< double > Rn;
    typedef KN_< double > Rn_;

    CMA_ES( ) : CMAES( ), x(0), fit(0) {}

    /*CMA_ES(ffcalfunc &_ff,int d,Rn &xstart,double *stddev,long seed,int lambda)
     *      : CMAES(d,xstart.n ? xstart:0,stddev,seed,lambda,"non"),x(&xstart),fit(&_ff) {}
     * CMA_ES(ffcalfunc &_ff,int d,Rn &xstart,const Rn &stddev,long seed,int lambda)
     *      : CMAES(d,xstart.n ? xstart:0,stddev,seed,lambda,"non"),x(&xstart),fit(&_ff) {}*/
    CMA_ES(ffcalfunc &_ff, Rn &xstart, const Rn &stddev, long seed, int lambda)
      : CMAES(xstart.n, xstart, stddev, seed, lambda, "non"), x(&xstart), fit(&_ff) {}

    CMA_ES(ffcalfunc &_ff, Rn &xstart, const Rn &stddev, long seed, int lambda, const string &ipf)
      : CMAES(xstart.n, xstart, stddev, seed, lambda, ipf.c_str( )), x(&xstart), fit(&_ff) {}

    void PopEval( ) {
      for (int i = 0; i < popsize( ); ++i) {
        Rn_ popi(pop[i], dimension( ));
        fitvals[i] = fit->J(popi);
      }
    }

    double *operator( )( ) { return *x = Rn(x->n, CMAES::operator( )( )); }

   private:
    ffcalfunc *fit;
    Rn *x;
  };

  class E_CMA_ES : public E_F0mps {
   public:
    const int cas;
    static basicAC_F0::name_and_type name_param[];
    static const int n_name_param = 11;
    Expression nargs[n_name_param];
    Expression X;
    C_F0 inittheparam, theparam, closetheparam;
    Expression JJ;
    long arg(int i, Stack stack, long a) const {
      return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
    }

    R arg(int i, Stack stack, R a) const { return nargs[i] ? GetAny< R >((*nargs[i])(stack)) : a; }

    E_CMA_ES(const basicAC_F0 &args, int cc) : cas(cc) {
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
      if (nbj > 0) {
        opJ = dynamic_cast< const Polymorphic * >(args[0].LeftValue( ));
        assert(opJ);
      }

      JJ = to< R >(C_F0(opJ, "(", theparam));
      closetheparam = C_F0((Expression)Block::snewclose(currentblock), atype< void >( ));
      // closetheparam=currentblock->close(currentblock);   // the cleanning block expression
    }

    virtual AnyType operator( )(Stack stack) const {
      double cost = 1e100;

      WhereStackOfPtr2Free(stack) = new StackOfPtr2Free(stack);    // FH mars 2005
      Kn &x = *GetAny< Kn * >((*X)(stack));
      long n = x.N( );
      long seed = arg(0, stack, 0L);                // The seed for random numbers generation
      double initialStdDev = arg(1, stack, 0.3);    // Initial standard deviation
      KN< double > iSD(n, 1.);
      iSD *= initialStdDev;
      KN< double > initialStdDevs(nargs[2] ? GetAny< KN_< double > >((*nargs[2])(stack))
                                           : (KN_< double >)iSD);
      double stopTolFun = arg(3, stack, 1.E-12);
      double stopTolFunHist = arg(4, stack, 0.);
      double stopTolX = arg(5, stack, 0.);
      double stopTolUpXFactor = arg(6, stack, 1.E3);
      long popsize = arg(7, stack, 4 + (long)floor(3 * log(n)));
      string pipf = nargs[10] ? *GetAny< string * >((*nargs[10])(stack)) : string("");
      ffcalfunc ffJ(stack, JJ, theparam);
      CMA_ES *optim = 0;
      if (pipf.size( ) > 0) {
        cout << "input file : " << pipf << endl;
        optim = new CMA_ES(ffJ, x, initialStdDevs, seed, popsize, pipf);
      } else {
        cout << "no input file " << endl;
        optim = new CMA_ES(ffJ, x, initialStdDevs, seed, popsize);
        long meval = arg(8, stack, static_cast< long >(optim->maxeval( )));
        long mgen = arg(9, stack, static_cast< long >(optim->maxgen( )));
        optim->optimizer( ).sp.stopTolFun = stopTolFun;
        optim->optimizer( ).sp.stopTolFunHist = stopTolFunHist;
        optim->optimizer( ).sp.stopTolX = stopTolX;
        optim->optimizer( ).sp.stopTolUpXFactor = stopTolUpXFactor;
        optim->optimizer( ).sp.stopMaxFunEvals = meval;
        optim->optimizer( ).sp.stopMaxIter = mgen;
      }

      (*optim)( );
      cost = optim->fitness( );
      x = KN_< double >(optim->xbestever( ), optim->dimension( ));
      cout << "Number of fitness evalution(s) : " << optim->eval( ) << endl;

      closetheparam.eval(stack);    // clean memory

      if (optim) {
        delete optim;
        optim = 0;
      }

      WhereStackOfPtr2Free(stack)->clean( );    // FH mars 2005
      return cost;                              // SetAny<long>(0);  Modif FH  july 2005
    }

    operator aType( ) const { return atype< double >( ); }
  };

  E_F0 *code(const basicAC_F0 &args) const { return new E_CMA_ES(args, cas); }

  OptimCMA_ES(int c)
    : OneOperator(atype< double >( ), atype< Polymorphic * >( ), atype< KN< R > * >( )), cas(c) {}
};

basicAC_F0::name_and_type OptimCMA_ES::E_CMA_ES::name_param[] = {
  {"seed", &typeid(long)},
  {"initialStdDev", &typeid(double)},
  {"initialStdDevs", &typeid(KN_< double >)},
  {"stopTolFun", &typeid(double)},
  {"stopTolFunHist", &typeid(double)},
  {"stopTolX", &typeid(double)},
  {"stopTolUpXFactor", &typeid(double)},
  {"popsize", &typeid(long)},
  {"stopMaxFunEval", &typeid(long)},
  {"stopMaxIter", &typeid(long)},
  {"paramFile", &typeid(string *)}};

static void Load_Init( ) { Global.Add("cmaes", "(", new OptimCMA_ES(1)); }

LOADFUNC(Load_Init)
