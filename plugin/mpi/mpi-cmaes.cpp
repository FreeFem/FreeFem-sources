// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Sylvain Auliac
// E-MAIL   : auliac@ann.jussieu.fr
//

/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

// ff-c++-LIBRARY-dep: mpi
// ff-c++-cpp-dep: ../seq/cmaes.cpp -I../seq

/*
        This is a freefem interface of the Hansen's CMA-ES C optimizer.
        The class CMAES is a quick C++ interface for the contents of the C files
        then follow the real FreeFem++ stuff.
*/

#include <iostream>
using namespace std;
#include "ff++.hpp"
#include "mpi.h"
#include "cmaes_interface.h"

template< class T >
struct MPI_TYPE {};
template<>
struct MPI_TYPE< long > {
  static MPI_Datatype TYPE( ) { return MPI_LONG; }
};
template<>
struct MPI_TYPE< int > {
  static MPI_Datatype TYPE( ) { return MPI_INT; }
};
template<>
struct MPI_TYPE< double > {
  static MPI_Datatype TYPE( ) { return MPI_DOUBLE; }
};
template<>
struct MPI_TYPE< char > {
  static MPI_Datatype TYPE( ) { return MPI_BYTE; }
};

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
  //{for(int i=0;i<popsize();++i) fitvals[i] = ff(pop[i]);} //the thing to parralelize

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
  double *&init(int dimension, double *xstart, double *stddev, long seed, int lambda,
                const char *input_parameter_filename) {
    return fitvals =
             cmaes_init(&evo, dimension, xstart, stddev, seed, lambda, input_parameter_filename);
  }

 private:
  void exit( ) { cmaes_exit(&evo); }
  mutable cmaes_t evo;
};

/*
        Now comes the FreeFem ++ part :
*/

extern Block *currentblock;

typedef double R;

class OptimCMA_ES : public OneOperator {
 public:
  typedef KN< R > Kn;
  typedef KN_< R > Kn_;
  const int cas;

  class ffcalfunc    //   to call the freefem function .. J
  {
   public:
    Stack stack;
    Expression JJ, theparame;
    mutable int counter;

    ffcalfunc(Stack s, Expression JJJ, Expression epar)
      : stack(s), JJ(JJJ), theparame(epar), counter(0) {}
    double J(Kn_ x) const {
      ++counter;
      KN< double > *p = GetAny< KN< double > * >((*theparame)(stack));
      *p = x;
      double ret = GetAny< R >((*JJ)(stack));
      WhereStackOfPtr2Free(stack)->clean( );
      return ret;
    }
  };

  class CMA_ES_MPI : public CMAES {
   public:
    typedef KN< double > Rn;
    typedef KN_< double > Rn_;

    // CMA_ES_MPI() : CMAES(),x(0),fit(0) {}
    /*CMA_ES(ffcalfunc &_ff,int d,Rn &xstart,double *stddev,long seed,int lambda)
            : CMAES(d,xstart.n ? xstart:0,stddev,seed,lambda,"non"),x(&xstart),fit(&_ff) {}
    CMA_ES(ffcalfunc &_ff,int d,Rn &xstart,const Rn &stddev,long seed,int lambda)
            : CMAES(d,xstart.n ? xstart:0,stddev,seed,lambda,"non"),x(&xstart),fit(&_ff) {}*/

    CMA_ES_MPI(ffcalfunc &_ff, Rn &xstart, const Rn &stddev, long seed, int lambda, MPI_Comm *_com)
      : CMAES( ), x(0), fit(&_ff), com(_com), myid(0), nproc(1), my_number_of_fitness_eval(0),
        index(0) {
      MPI_Comm_size(*com, &nproc);
      MPI_Comm_rank(*com, &myid);
      x = &xstart;
      init(x->n, xstart.n ? xstart : 0, stddev, seed, lambda, "non");
      my_number_of_fitness_eval = (lambda / nproc) + (myid < (lambda % nproc) ? 1 : 0);
      index = new int[nproc];
      for (int i = 0; i < nproc; ++i) {
        int inofe = lambda / nproc + ((i - 1) < (lambda % nproc) ? 1 : 0);
        index[i] = i > 0 ? index[i - 1] + inofe : 0;
      }
      if (myid == 0) cout << SayHello( ) << endl;
    }
    CMA_ES_MPI(ffcalfunc &_ff, Rn &xstart, const Rn &stddev, long seed, int lambda, MPI_Comm *_com,
               const char *ifname)
      : CMAES( ), x(0), fit(&_ff), com(_com), myid(0), nproc(1), my_number_of_fitness_eval(0),
        index(0) {
      MPI_Comm_size(*com, &nproc);
      MPI_Comm_rank(*com, &myid);
      x = &xstart;
      init(x->n, xstart.n ? xstart : 0, stddev, seed, lambda, ifname);
      my_number_of_fitness_eval = (lambda / nproc) + (myid < (lambda % nproc) ? 1 : 0);
      index = new int[nproc];
      for (int i = 0; i < nproc; ++i) {
        int inofe = lambda / nproc + ((i - 1) < (lambda % nproc) ? 1 : 0);
        index[i] = i > 0 ? index[i - 1] + inofe : 0;
      }
      if (myid == 0) cout << SayHello( ) << endl;
    }
    ~CMA_ES_MPI( ) {
      if (index) delete[] index;
      index = 0;
    }

    void PopEval( ) {
      // cout << endl << "***** in popeval myid = " << myid << "*****" << endl;
      for (int i = 0; i < my_number_of_fitness_eval; ++i) {
        Rn_ popi(pop[i + index[myid]], dimension( ));
        fitvals[i + index[myid]] = fit->J(popi);
      }
    }
    double *operator( )( ) {
      while (!TestForTermination( )) {
        MPI_Barrier(*com);
        SamplePopulation( );    // Every proc samples its own population... but only the one from
                                // proc 0 is used. Fortunately, this is a quick computation regarding
                                // the time spent on fitness function evaluation
        for (int i = 0; i < popsize( ); ++i)
          MPI_Bcast(pop[i], dimension( ), MPI_DOUBLE, 0,
                    *com);    // send the population to every proc
        PopEval( );           // each proc do its work gently
        for (int i = 0; i < nproc; ++i) {
          int nn = i < nproc - 1 ? index[i + 1] - index[i] : popsize( ) / nproc;
          MPI_Bcast(&fitvals[index[i]], nn, MPI_DOUBLE, i, *com);
        }
        UpdateDistribution( );    // There again, only proc 0 would really need to update its
                                  // distribution
      }
      if (myid == 0) cout << "Stop : " << TestForTermination( ) << endl;
      return xmean( );
    }

   private:
    CMA_ES_MPI(const CMA_ES_MPI &);
    CMA_ES_MPI &operator=(const CMA_ES_MPI &);
    ffcalfunc *fit;
    Rn *x;
    MPI_Comm *com;
    int nproc, myid, my_number_of_fitness_eval;
    int *index;
  };

  class E_CMA_ES : public E_F0mps {
   public:
    const int cas;
    static basicAC_F0::name_and_type name_param[];
    static const int n_name_param = 12;
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
      //  the expression to init the theparam of all
      inittheparam =
        currentblock->NewVar< LocalVariable >("the parameter", atype< KN< R > * >( ), X_n);
      theparam = currentblock->Find("the parameter");    //  the expression for the parameter
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
      // cout << "dans le dylib :" << initialStdDevs << endl;
      double stopTolFun = arg(3, stack, 1.E-12);
      double stopTolFunHist = arg(4, stack, 0.);
      double stopTolX = arg(5, stack, 0.);
      double stopTolUpXFactor = arg(6, stack, 1.E3);
      long popsize = arg(7, stack, 4 + (long)floor(3 * log(n)));
      pcommworld vcommworld = 0;
      if (nargs[8]) vcommworld = GetAny< pcommworld >((*nargs[8])(stack));

      MPI_Comm mpiCommWorld = MPI_COMM_WORLD;
      MPI_Comm *commworld = vcommworld ? (MPI_Comm *)vcommworld : &mpiCommWorld;
      // long mu = arg(8,stack,popsize/2);

      long iprint = verbosity;
      ffcalfunc ffJ(stack, JJ, theparam);

      // cout << endl << "ATTENTION :  " << *commworld << endl;

      int myid = 0, nproc = 1;

      MPI_Comm_size(*commworld, &nproc);
      MPI_Comm_rank(*commworld, &myid);

      // cout << endl << "nbr de proc : " << nproc << " -- myid=" << myid << " -- world id=" << wr
      // <<  endl;

      CMA_ES_MPI *optim = 0;
      if (nargs[9])
        optim = new CMA_ES_MPI(ffJ, x, initialStdDevs, seed, popsize, commworld,
                               (GetAny< string * >((*nargs[9])(stack)))->c_str( ));
      else
        optim = new CMA_ES_MPI(ffJ, x, initialStdDevs, seed, popsize, commworld);
      if (!nargs[9]) {
        optim->optimizer( ).sp.stopTolFun = stopTolFun;
        optim->optimizer( ).sp.stopTolFunHist = stopTolFunHist;
        optim->optimizer( ).sp.stopTolX = stopTolX;
        optim->optimizer( ).sp.stopTolUpXFactor = stopTolUpXFactor;
        long meval = arg(10, stack, static_cast< long >(optim->maxeval( )));
        long mgen = arg(11, stack, static_cast< long >(optim->maxgen( )));
        optim->optimizer( ).sp.stopMaxFunEvals = meval;
        optim->optimizer( ).sp.stopMaxIter = mgen;
      }

      (*optim)( );
      cost = optim->fitness( );
      x = KN_< double >(optim->xbestever( ), optim->dimension( ));

      delete optim;
      closetheparam.eval(stack);                // clean memory
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
  {"comm", &typeid(pcommworld)},
  {"paramFile", &typeid(string *)},
  {"stopMaxFunEval", &typeid(long)},
  {"stopMaxIter", &typeid(long)}

  //{"mu",							&typeid(long) }
};

/* --FH:   class Init { public:
  Init();
  };*/

static void Load_Init( ) { Global.Add("cmaesMPI", "(", new OptimCMA_ES(1)); }

LOADFUNC(Load_Init)
