//#include "HashMatrix.hpp"
//#include "lgmat.hpp"

#ifndef FFLANG
#if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)

#define BOOST_NO_CXX17_IF_CONSTEXPR
#include <ff++.hpp>
#include <AFunction_ext.hpp>
#include <lgfem.hpp>
#include <R3.hpp>

#include <htool/htool.hpp>

// include the bemtool library .... path define in where library
//#include <bemtool/operator/block_op.hpp>
#include <bemtool/tools.hpp>
#include <bemtool/fem/dof.hpp>
#include <bemtool/operator/operator.hpp>
#include <bemtool/miscellaneous/htool_wrap.hpp>
//#include "PlotStream.hpp"
#include "HashMatrix.hpp"
#include "common.hpp"

//extern FILE *ThePlotStream;

using namespace std;
//using namespace htool;
//using namespace bemtool;

#include <type_traits>

//typedef   LinearComb<MGauche,C_F0> Finconnue;
//typedef   LinearComb<MDroit,C_F0> Ftest;
//typedef  const Finconnue  finconnue;
//typedef  const Ftest ftest;
//class CDomainOfIntegration;
//class FormBilinear;
#include "common_bem.hpp"
#include "bem.hpp"

#endif
#endif

#include  <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"

//#include "lex.hpp"
#include "HashMatrix.hpp"

#include "SparseLinearSolver.hpp"
#include "Mesh3dn.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include <set>
#include "compositeFESpace.hpp"


string typeFEtoString(int typeFE)
{
  string toto="";
  if(typeFE == 2){
    toto= "FESpace2D Vol";
  }else if(typeFE == 3){
    toto= "FESpace3D Vol";
  }else if(typeFE == 4){
    toto= "FESpace3D Surf";
  }else if(typeFE == 5){
    toto= "FESpace3D Curve";
  }
  else{
    cerr << "func typeFEtoString ::  error in the type of FESpace " << endl;
    ffassert(0);
  }
  return toto;
}


/**
       *  @brief  Builds a new largs  whit each element are included in one block
       *  @param  largs list of argument of the initial Forms 
       *  @param  NpUh  number of FESpace in Uh
       *  @param  NpVh  number of FESpace in Vh
       *  @param  indexBlockUh give the index of the block for a given component of FESpace Uh  
       *  @param  indexBlockVh give the index of the block for a given component of FESpace Vh  
       * 
       */

list<C_F0>  creationLargsForCompositeFESpace( const list<C_F0> & largs, const int &NpUh, const int &NpVh, 
                                        const KN<int> &indexBlockUh, const KN<int> &indexBlockVh ){
  
  // At the end of this function, every element of newlargs is included in one block
  list<C_F0> newlargs; // creation de la nouvelle list largs

  // const_iterator
  list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end(); 
  // loop over largs information 
  if(verbosity>3) cout << "loop over the integral" << endl;

  int count_integral = 0;
  for (ii=ib;ii != ie;ii++) {
    count_integral++;
    if(verbosity>3){
      cout <<"========================================================" << endl;
      cout <<"=                                                      =" << endl;
      cout << "reading the " << count_integral << "-th integral" << endl;
    }
    Expression e=ii->LeftValue();
    aType r = ii->left();
    if(verbosity>3){
      cout << "e=" << e << ", " << "r=" << r << endl;
      cout <<"=                                                      =" << endl;
    }

    // ***************************************
    // Case FormBillinear
    // ***************************************
    if (r==atype<const  FormBilinear *>() ){
      const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      ffassert(bb);
      const CDomainOfIntegration & di= *bb->di;

      BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
      if (Op == NULL) {
        if(mpirank == 0) cout << "dynamic_cast error" << endl; 
        ffassert(0);
      }
      
      size_t Opsize= Op->v.size();
      if( verbosity > 3 ){
        cout << " loop over the term inside the integral" << endl;
        cout << " Number of term in the integral:: Op->v.size()=" << Op->v.size() << endl;
      }
      std::vector< std::pair<int,int> > indexBlock(Opsize);

      for(size_t jj=0; jj<Opsize; jj++){
        indexBlock[jj] = std::pair<int,int>( 
          indexBlockUh[ Op->v[jj].first.first.first ],
          indexBlockVh[ Op->v[jj].first.second.first ] );
      }

      // index to check if a integral is defined on multi block

      int countOP=0;
      BilinearOperator * OpBloc= new BilinearOperator();
      // loop over the block
      for(int ibloc=0; ibloc<NpUh; ibloc++){
        for(int jbloc=0; jbloc<NpVh; jbloc++){

          //BilinearOperator * OpBloc= new BilinearOperator();
          countOP=0;
          for(size_t jj=0; jj<Opsize; jj++){
            if( indexBlock[jj].first == ibloc && indexBlock[jj].second == jbloc){
              if (countOP == 0){
                ffassert( OpBloc->v.empty() );
              } 
              OpBloc->add(Op->v[jj].first, Op->v[jj].second); // Add the billinearOperator to bloc (ibloc,jbloc).
              countOP += 1;
            }
          }
          
          if( countOP > 0 ){   
            // <<  countOP << " voila titi " << "OpBloc->v.size()= " << OpBloc->v.size() << endl; 
            ffassert( OpBloc->v.size() > 0); 
            
            // FormBilinear *titi = new FormBilinear( &di, OpBloc );
            newlargs.push_back( C_F0( new FormBilinear( &di, OpBloc ), r ) ); 
            
            // // Add to the given block
            // for(size_t jj=0; jj<OpBloc->v.size(); jj++){
            //   OpBloc->v[jj].first.first.first = localIndexInTheBlockUh( OpBloc->v[jj].first.first.first ); // finconnue
            //   OpBloc->v[jj].first.second.first = localIndexInTheBlockVh( OpBloc->v[jj].first.second.first ); // ftest
            // }
            // block_largs(ibloc,jbloc).push_back( C_F0( new FormBilinear( &di, OpBloc ), r) );

            OpBloc->v.clear();
            ffassert( (OpBloc->v.empty() == true ) );

          }
          
        }
      }
      delete OpBloc;
    } // end billinear type
    // ***************************************
    // Case linear type 
    // ***************************************
    else if(r == atype< const FormLinear *>()){
      //  if(verbosity>3) cout << "FormLinear in variational form" << endl;
      const FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
      LOperaD * Op = const_cast<LOperaD *>(ll->l);
      if (Op == NULL) {
        if(mpirank == 0) cout << "dynamic_cast error" << endl; 
        ffassert(0);
      }

      size_t Opsize= Op->v.size();
      
      // creation the vector of the indexBlock for each element OpChange
      std::vector< int > indexBlock(Opsize);

      for(size_t jj=0; jj<Opsize; jj++){
        indexBlock[jj] = indexBlockVh[ Op->v[jj].first.first ];
      }

      LOperaD * OpBloc= new LOperaD();

      // put the term inside OpChange in the good block
      for(int jbloc=0; jbloc<NpVh; jbloc++){
        int countOP=0;
        //LOperaD * OpBloc= new LOperaD();
        for(size_t jj=0; jj<Opsize; jj++){
          if( indexBlock[jj] == jbloc ){
              if (countOP == 0){
                ffassert( OpBloc->v.empty() );
              } 
              OpBloc->add(Op->v[jj].first,Op->v[jj].second); // Add the LinearOperator to bloc (ibloc,jbloc).
              countOP += 1;
          }
        }
        if( countOP > 0 ){   
          ffassert( OpBloc->v.size() > 0);   
           
          // // for coonstruction of block_largs
          // for(size_t jj=0; jj<OpBloc->v.size(); jj++){
          //   OpBloc->v[jj].first.first = localIndexInTheBlockVh( OpBloc->v[jj].first.first );
          // }
          newlargs.push_back( C_F0( new FormLinear( (ll->di), OpBloc ), r ) );

          // // for coonstruction of block_largs
          // block_largs(jbloc,jbloc).push_back( C_F0( new FormLinear( (ll->di), OpBloc ), r ) );

          OpBloc->v.clear(); 
          ffassert( ( OpBloc->v.empty() == true)  ); // check if OpBloc is empty after clear();
        }
      }
      delete OpBloc;
    } // end linear type
    // ***************************************
    // Case Boundary condition
    // ***************************************
    else if(r == atype<const  BC_set  *>()){
      if(verbosity>3) cout << " BC in variational form " << endl;
      BC_set * bc=dynamic_cast< BC_set *>(e);

      int sizebc=bc->bc.size();
      std::vector< int > indexBlock(sizebc);
      // calculate the index of the componenent where the bloc
      for (int k=0; k<sizebc; k++)
      {
        pair<int,Expression> &xx=bc->bc[k];
        indexBlock[k] = indexBlockUh[xx.first];
        // xx.first = localIndexInTheBlockVh(xx.first); // change the index corresponding to the block
      }
      bool addBC = false;
      bool *okBC =new bool[NpUh];
      for(int ibloc=0; ibloc<NpUh; ibloc++){
         addBC = false;
        // construction of okBC for ibloc 
        for (int k=0; k<sizebc; k++){
          if(indexBlock[k] == ibloc){
            okBC[k] = true;
            addBC = true;
          }
          else{
            okBC[k] = false;
          }
        }
        // add the BC_set correspond to the ibloc of the composite FESpace 
        if(addBC) newlargs.push_back( C_F0( new BC_set(*bc,okBC), r) ); 
        // if(addBC) block_largs(ibloc,ibloc).push_back( C_F0( new BC_set(*bc,okBC), r) );
      }
      delete [] okBC;
    } // end BC
#ifndef FFLANG
#if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
    // ******************************************
    // Case BemKFormBilinear (KERNEL FORM ONLY)
    // ******************************************
    else if (r==atype<const BemFormBilinear *>() ){
      BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
      ffassert(bbtmp);
      int VVFBEM = bbtmp->type;

      if(VVFBEM ==1){
        BemKFormBilinear * bb= dynamic_cast< BemKFormBilinear *>(e);
        ffassert(bb);
      
        // creation of the new operator
        BemKFormBilinear * bbnew = new BemKFormBilinear( bb->di, FoperatorKBEM(bb->b->kbem, *(bb->b->fi), *(bb->b->ft) ) ); 
        newlargs.push_back( C_F0( bbnew, r ) );
      }
      else if(VVFBEM == 2){
        cerr << " BEM Potential in composite FESpace in construction " << endl;
        ffassert(0);
      }
      else{
        cerr << "VFBEM must be egal to 1(kernel) or 2(potential)" << endl;
        ffassert(0);
      }
    }
#endif
#endif

  } // end iterator oflargs
  return newlargs;
}

KNM< list<C_F0> > computeBlockLargs( const list<C_F0> & largs, const int &NpUh, const int &NpVh, const KN<int> &indexBlockUh, const KN<int> &indexBlockVh ){
  // creation of the return value 
  KNM< list<C_F0> > block_largs( (long)NpUh, (long)NpVh );

  long UhtotalNbItem = indexBlockUh.N();
  long VhtotalNbItem = indexBlockVh.N();    

  // impression des information de la composition largs
  list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end(); 

  // necessaire :: UhtotalNbItem, indexBlockUh

  // Loop to put each term of the varf in each correct block

  // loop over largs information 
  if(verbosity>3) cout << "loop over the integral" << endl;

  int count_integral = 0;
  for (ii=ib;ii != ie;ii++) {
    count_integral++;
    
    Expression e=ii->LeftValue();
    aType r = ii->left();
    if(verbosity>3){
      cout <<"========================================================" << endl;
      cout <<"=                                                      =" << endl;
      cout << "reading the " << count_integral << "-th term of the variational form used to define the matrix" << endl;
      cout << "e=" << e << ", " << "r=" << r << endl;
      cout <<"=                                                      =" << endl;
    }
    // ***************************************
    // Case FormBillinear
    // ***************************************
    if (r==atype<const  FormBilinear *>() ){
      const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      const CDomainOfIntegration & di= *bb->di;

      if(verbosity >3){
        cout << "di.kind=" << di.kind << ", di.dHat=" << di.dHat << endl;
        cout << "di.d=   " << di.d    << ", di.Th=  " << di.Th << endl;
      }
      
      int    d = di.d;
      int dHat = di.dHat;

      // Sert a verifier que "*bb->di->Th" est du bon type ==> A enlever    
      BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
      if (Op == NULL) {
        if(mpirank == 0) cout << "dynamic_cast error" << endl; 
        ffassert(0);
      }
      
      size_t Opsize= Op->v.size();
      if(verbosity >3){
        cout << " loop over the term inside the integral" << endl;
        cout << " Number of term in the integral:: Op->v.size()=" << Op->v.size() << endl;
      }
    
      // index to check if a integral is defined on multi block
      int indexOfBlockUh = -1; // A changer de nom
      int indexOfBlockVh = -1; // A changer de nom
      for(size_t jj=0; jj<Opsize; jj++){
        // attention la fonction test donne la ligne
        //  et la fonction test est en second
        BilinearOperator::K ll = Op->v[jj];
        pair<int,int> finc(ll.first.first), ftest(ll.first.second);
        if(verbosity>3){
          cout << " operateur jj= " << jj << endl;
          cout << " FormBilinear: number of unknown finc=" <<  finc.first << " ,ftest= " << ftest.first << endl;
          cout << " FormBilinear: operator order finc   =" << finc.second << " ,ftest= " << ftest.second << endl; // ordre   only op_id=0
        }
        
        // Fred fait peut être un message après ????
        // verification que la taille des tableaux des fonctions tests et de la fonction inconnue``
        // sont correctes.  
        ffassert( -1  < finc.first  && finc.first < UhtotalNbItem);
        ffassert( -1  < ftest.first && ftest.first < VhtotalNbItem);

        // finc.first : index de component de la fonction inconnue
        // ftest.first: index de component de la fonction test
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // finc.second : renvoie l'index du type de l'operateur: Id, dx(), dy(), dz(), dxx(), dxy(), ...
        //
        // la liste des index des operateurs est definis dans [[????]].
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // exemple vf([u1,u2,..,u30],[v1,v2,...,v30]) = int2d(Th)(dx(u20)*v15)
        //      finc.first  = 20 , ftest.first = 15
        //      finc.second = 1 , ftest.second = 0

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if( jj== 0 ){
          indexOfBlockUh = indexBlockUh[finc.first];
          indexOfBlockVh = indexBlockVh[ftest.first];
        }
        else if( indexOfBlockUh != indexBlockUh[finc.first] ){
          cerr << "The " << count_integral <<"-th integral(s) contains the constribution of two different blocks:" << endl;
          cerr << "the first term correspond to block (" << indexOfBlockUh << " , " <<  indexOfBlockVh << ")" << endl;
          cerr << "the "<<jj <<"-th term correspond to block (" << indexBlockUh(finc.first) << " " <<  indexBlockVh(ftest.first) << ")" << endl;
          cerr << "You need to separate the integral in individual part." << endl;
          cerr << "Remark: scalar product in N dimension correspond to N terms inside the integral" << endl;
          cerr << "A ameliorer Jacques." << endl;
          ffassert(0);
        }
        else if( indexOfBlockVh != indexBlockVh[ftest.first] ){
          cerr << "The " << count_integral <<"-th integral(s) contains the constribution of two different blocks:" << endl;
          cerr << "the first term correspond to block (" << indexOfBlockUh << " , " <<  indexOfBlockVh << ")" <<endl;
          cerr << "the "<<jj <<"-th term correspond to block (" << indexBlockUh(finc.first) << ", " <<  indexBlockVh(ftest.first) << ")" << endl;
          cerr << "You need to separate the integral in individual part." << endl;
          cerr << "Remark: scalar product in N dimension correspond to N terms inside the integral" << endl;
          cerr << "A ameliorer Jacques." << endl;
          ffassert(0);
        }

        ffassert( indexOfBlockUh == indexBlockUh(finc.first) );
        ffassert( indexOfBlockVh == indexBlockVh(ftest.first) );
      }

    ffassert( indexOfBlockUh >= 0 && indexOfBlockVh >= 0);
    
    // A faire :: recuperation des éléments pour chacun des blocs
    // Actuellement, on associe une intégrale par block ==> 
    
    block_largs(indexOfBlockUh,indexOfBlockVh).push_back(*ii);

    if(verbosity>3) cout << "The " << count_integral <<"-th integral(s) is added to the block (" << indexOfBlockUh << " , " <<  indexOfBlockVh << ")" <<endl;

    }
#ifndef FFLANG
#if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
    // ******************************************
    // Case BemKFormBilinear (KERNEL FORM ONLY)
    // ******************************************
    else if (r==atype<const BemFormBilinear *>() ){
      BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
      int VVFBEM = bbtmp->type;

      if(VVFBEM ==1){
        //BemKFormBilinear * bb=new BemKFormBilinear(*dynamic_cast<const BemKFormBilinear *>(e));
        const BemKFormBilinear * bb= dynamic_cast<const BemKFormBilinear *>(e);
        FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
        if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; ffassert(0);}

        int indexOfBlockUh = -1; // A changer de nom
        int indexOfBlockVh = -1; // A changer de nom


        // loop over the index of finconnue
        LOperaG * OpG = const_cast<LOperaG *>(b->fi);
        ffassert( OpG->v.size() == 1);
        size_t jj =0;
        for (LOperaG::const_iterator lop=OpG->v.begin();lop!=OpG->v.end();lop++){

          LOperaG::K lf(*lop);
          pair<int,int> finc(lf.first);
          if(verbosity >3){
            cout << " operateur jj= " << jj << endl;
            cout << " BemFormLinear: number of unknown finc= " << finc.first << endl;
          }
          ffassert( -1  < finc.first && finc.first < UhtotalNbItem);     // check the index 

          if( jj== 0 ){
            indexOfBlockUh = indexBlockUh[finc.first];
          }
          else if( indexOfBlockUh != indexBlockUh[finc.first] ){
            cerr << "The " << count_integral <<"-th term of the varitional form contains the constribution of two different FESpace:" << endl;
            cerr << "This terms correspond to a BEM integral terms" << endl;
            cerr << "the first term correspond to element " << indexOfBlockUh << " of the Composite FESpace (Finconnu)." << endl;
            cerr << "the "<< jj <<"-th term correspond to element " << indexBlockUh(finc.first) << endl;
            cerr << "In a composite FESpace, you need to define a BEM integral for each FESpace individually." << endl;
            cerr << "A ameliorer Jacques." << endl;
            ffassert(0);
          }
          jj+=1;
        }


        // Loop over the index of ftest
        LOperaD * OpD = const_cast<LOperaD *>(b->ft);
        ffassert( OpD->v.size() == 1);
        jj =0; // reinitialisation ton zero
        for (LOperaD::const_iterator lop=OpD->v.begin();lop!=OpD->v.end();lop++){
          
          LOperaD::K lf(*lop);
          pair<int,int> ftest(lf.first);
          if(verbosity>3){
            cout << " operateur jj= " << jj << endl;
            cout << " BemFormLinear: number of unknown ftest= " << ftest.first << endl;
          }
          ffassert( -1  < ftest.first && ftest.first < VhtotalNbItem);    // check the index 

          if( jj== 0 ){
            indexOfBlockVh = indexBlockVh[ftest.first];
          }
          else if( indexOfBlockVh != indexBlockVh[ftest.first] ){
            cerr << "The " << count_integral <<"-th term of the varitional form contains the constribution of two different FESpace:" << endl;
            cerr << "This terms correspond to a BEM integral terms" << endl;
            cerr << "the first term correspond to element " << indexOfBlockVh << " of the Composite FESpace (Ftest)." << endl;
            cerr << "the "<< jj <<"-th term correspond to element " << indexBlockVh(ftest.first) << endl;
            cerr << "In a composite FESpace, you need to define a BEM integral term for each FESpace individually." << endl;
            cerr << "A ameliorer Jacques." << endl;
            ffassert(0);
          }
          jj+=1;
        }
        block_largs(indexOfBlockUh,indexOfBlockVh).push_back(*ii); 

      }else if(VVFBEM == 2){
        BemPFormBilinear * bb=new BemPFormBilinear(*dynamic_cast<const BemPFormBilinear *>(e));
        FoperatorPBEM * b=const_cast<  FoperatorPBEM *>(bb->b);
        if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; }

        cerr << " BEM Potential in composite FESpace in construction " << endl;
        ffassert(0);
      }

    }
#endif
#endif
    else if (r==atype<const  FormLinear *>() ){
      // === FormLinear === //
      const FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
      LOperaD * Op = const_cast<LOperaD *>(ll->l);
      if (Op == NULL) {
        if(mpirank == 0) cout << "dynamic_cast error" << endl; 
        ffassert(0);
      }

      int indexOfBlockVh = -1;
      size_t Opsize= Op->v.size();
      size_t jj =0;
      for (LOperaD::const_iterator lop=Op->v.begin();lop!=Op->v.end();lop++)
      {          
        LOperaD::K lf(*lop);
        pair<int,int> ftest(lf.first);
        if(verbosity>3){
          cout << " operateur jj= " << jj << endl;
          cout << " FormLinear: number of unknown ftest= " << ftest.first << endl;
        }
        ffassert( -1  < ftest.first && ftest.first < VhtotalNbItem);

        if( jj== 0 ){
          indexOfBlockVh = indexBlockVh[ftest.first];
        }
        else if( indexOfBlockVh != indexBlockVh[ftest.first] ){
          cerr << "The " << count_integral <<"-th term of the varitional form contains the constribution of two different FESpace:" << endl;
          cerr << "This terms correspond to an integral terms" << endl;
          cerr << "the first term correspond to element " << indexOfBlockVh << " of the Composite FESpace " << endl;
          cerr << "the "<< jj <<"-th term correspond to element " << indexBlockVh(ftest.first) << endl;
          cerr << "In a composite FESpace, you need to define a BC for each FESpace individually." << endl;
          cerr << "A ameliorer Jacques." << endl;
          ffassert(0);
        }
        jj+=1;
      }
      // Added the boundary condition in the largs block
      block_largs(indexOfBlockVh,indexOfBlockVh).push_back(*ii); 
    }
    else if(r == atype<const  BC_set  *>()){
      if(verbosity >3) cout << " BC in variational form " << endl;
      
      const BC_set * bc=dynamic_cast<const  BC_set *>(e);
      
      // index to check if a integral is defined on multi block
      int indexOfBlockUh = -1;
    
      int kk=bc->bc.size();
      for (int k=0;k<kk;k++)
      {
          pair<int,Expression> xx=bc->bc[k];
          ffassert( -1  < xx.first  && xx.first < UhtotalNbItem); // check the value of index of the component of the varf
        
          if( k == 0) indexOfBlockUh = indexBlockUh[xx.first]; // index of the block Uh
          else if( indexOfBlockUh != indexBlockUh[xx.first] ){
            cerr << "The " << count_integral <<"-th term of the varitional form contains the constribution of two different FESpace:" << endl;
            cerr << "This terms correspond to Boundary condition" << endl;
            cerr << "the first term correspond to element " << indexOfBlockUh << " of the Composite FESpace " << endl;
            cerr << "the "<<kk <<"-th term correspond to element " << indexBlockUh(xx.first) << endl;
            cerr << "In a composite FESpace, you need to define a BC for each FESpace individually." << endl;
            cerr << "A ameliorer Jacques." << endl;
            ffassert(0);
          }
      }
      // Added the boundary condition in the largs block
      block_largs(indexOfBlockUh,indexOfBlockUh).push_back(*ii); 
        
      //ffassert(0);
    }
    else{
      
      cerr << "Composite FESpace only :: bilinear form" << endl;
      cerr << "                       :: linear form" << endl;
      cerr << "                       :: BC" << endl;
      #ifndef FFLANG
      #if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
      cerr << "                       :: BemFormBilinear" << endl;
      #endif
      #endif
      ffassert(0);
    }
  }
  return block_largs;
}


// Info necessaire :: " block_largs, localIndexInTheBlockUh, localIndexInTheBlockVh, NpUh, NpVh  
void changeComponentFormCompositeFESpace( const KN<int> &localIndexInTheBlockUh, const KN<int> &localIndexInTheBlockVh, 
        KNM< list<C_F0> > & block_largs ){
  // put the right number of each component of each block

  long NpUh = block_largs.N();
  long NpVh = block_largs.M();

  for( long i=0; i<NpUh; i++){
      for( long j=0; j<NpVh; j++){
        
        const list<C_F0> *b_largs=&block_largs(i,j); 
        list<C_F0>::const_iterator b_ii,b_ib=b_largs->begin(),b_ie=b_largs->end(); 
        for (b_ii=b_ib;b_ii != b_ie;b_ii++){
          Expression e=b_ii->LeftValue();
          aType r = b_ii->left();
          // Case FormBilinear
          if (r==atype<const  FormBilinear *>() ){
            const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      
            BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
            if (Op == NULL) {
              if(mpirank == 0) cout << "dynamic_cast error" << endl; 
            ffassert(0);
            }
      
            size_t Opsize= Op->v.size();
            if(verbosity>3){
              cout << " loop over the term inside the integral" << endl;
              cout << " Number of term in the integral:: Op->v.size()=" << Op->v.size() << endl;
            }
            KN<size_t> index_operator_finc(Opsize);
            KN<int>    new_index_funct_finc(Opsize);

            KN<size_t> index_operator_ftest(Opsize);
            KN<int>    new_index_funct_ftest(Opsize);

            // index to check if a integral is defined on multi block
            for(size_t jj=0; jj<Opsize; jj++){
              // attention la fonction test donne la ligne
              //  et la fonction test est en second
              BilinearOperator::K ll = Op->v[jj];
              pair<int,int> finc(ll.first.first), ftest(ll.first.second);

              long jj2= jj;

              index_operator_finc[ jj2] = jj;
              new_index_funct_finc[ jj2] = localIndexInTheBlockUh(finc.first);
            
              index_operator_ftest[ jj2]  = jj;
              new_index_funct_ftest[ jj2] = localIndexInTheBlockVh(ftest.first);

            }
            changeIndexFunctionInconnue(*Op, index_operator_finc, new_index_funct_finc );
      
            changeIndexFunctionTest(*Op, index_operator_ftest, new_index_funct_ftest  );          
          }  

          #ifndef FFLANG
          #if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
          // ******************************************
          // Case BemKFormBilinear (KERNEL FORM ONLY)
          // ******************************************
          else if (r==atype<const BemFormBilinear *>() ){
            BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
            int VVFBEM = bbtmp->type;

            if(VVFBEM ==1){
              //BemKFormBilinear * bb=new BemKFormBilinear(*dynamic_cast<const BemKFormBilinear *>(e));
              const BemKFormBilinear * bb= dynamic_cast<const BemKFormBilinear *>(e);
              FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
              if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; exit(0);}

              // loop over the index of finconnue
              LOperaG * OpG = const_cast<LOperaG *>(b->fi);
              ffassert( OpG->v.size() == 1);
              size_t Opsize= OpG->v.size();

              for(size_t jjj=0; jjj<Opsize; jjj++){
                LOperaG::K lf=OpG->v[jjj];
                OpG->v[jjj].first.first = localIndexInTheBlockUh( OpG->v[jjj].first.first );
              
                pair<int,int> finc(lf.first);
                if(verbosity>3){
                  cout << " new value :: block i,j=" << i << ","<< j << ", operateur jj= " << jjj << endl;
                  cout << " BemormBilinear: number of unknown finc = " << finc.first << endl;
                  cout << " BemFormBilinear: operator order   finc = " << finc.second << endl; 
                }
              }


              // Loop over the index of ftest
              LOperaD * OpD = const_cast<LOperaD *>(b->ft);
              ffassert( OpD->v.size() == 1);
              Opsize= OpD->v.size();
              for(size_t jjj=0; jjj<Opsize; jjj++){
                LOperaD::K lf=OpD->v[jjj];
                OpD->v[jjj].first.first = localIndexInTheBlockVh( OpD->v[jjj].first.first );
              
                pair<int,int> ftest(lf.first);
                if(verbosity>3){
                  cout << " new value :: block i,j=" << i << ","<< j << ", operateur jj= " << jjj << endl;
                  cout << " BemormBilinear: number of unknown ftest = " << ftest.first << endl;
                  cout << " BemFormBilinear: operator order   ftest = " << ftest.second << endl; 
                }
              }
            }
            else if(VVFBEM == 2){
              BemPFormBilinear * bb=new BemPFormBilinear(*dynamic_cast<const BemPFormBilinear *>(e));
              FoperatorPBEM * b=const_cast<  FoperatorPBEM *>(bb->b);
              if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; }

              cerr << " BEM Potential in composite FESpace in construction " << endl;
              ffassert(0);
            }

          }
          #endif
          #endif
          // LinearForm
          else if (r==atype<const  FormLinear *>() ){
            const FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
            LOperaD * Op = const_cast<LOperaD *>(ll->l);
            if (Op == NULL) {
              if(mpirank == 0) cout << "dynamic_cast error" << endl; 
              ffassert(0);
            }

            size_t Opsize= Op->v.size();
            for(size_t jjj=0; jjj<Opsize; jjj++){
              Op->v[jjj].first.first = localIndexInTheBlockVh( Op->v[jjj].first.first );
            }
          }
          // case BC_set 
          else if(r == atype<const  BC_set  *>()){
            ffassert( i == j ); // diagonal block 
            BC_set * bc=dynamic_cast< BC_set *>(e); // on ne peut pas utiliser " const BC_set * " ou autrement erreur ce ompilation:  Morice

            //KN<int>  new_index_funct_finc( bc.size() );
            int kk=bc->bc.size();
            //pair<int,Expression>  &bc_ib(bc->bc.begin());
            
            for (int k=0;k<kk;k++)
            {
              pair<int,Expression> &xx2= bc->bc[k];
              //new_index_funct_finc[k] = localIndexInTheBlockUh(bc[k].first);
              // change the index of the component to correspond to the index in the block
              xx2.first = localIndexInTheBlockUh(xx2.first);
              //bc->changeNumberOfComponent(k,localIndexInTheBlockUh(xx.first));
              //bc->bc[k].first = localIndexInTheBlockUh( bc->bc[k].first );
            }
          }
        }
        // listOfComponentBilinearForm(*b_largs);
      }
  }
}


void reverseChangeComponentFormCompositeFESpace(const KN<int>  &beginBlockUh, const KN<int> &beginBlockVh, 
          KNM< list<C_F0> > & block_largs){
  // Info necessaire :: " block_largs, beginBlockUh, beginBlockVh, NpUh, NpVh  

  long NpUh = block_largs.N();
  long NpVh = block_largs.M();

  // put the right number of component of each block
  for( int i=0; i<NpUh; i++){
      for( int j=0; j<NpVh; j++){
        
        const list<C_F0> *b_largs=&block_largs(i,j); 
        list<C_F0>::const_iterator b_ii,b_ib=b_largs->begin(),b_ie=b_largs->end(); 
        for (b_ii=b_ib;b_ii != b_ie;b_ii++){
          Expression e=b_ii->LeftValue();
          aType r = b_ii->left();

          // bilinear case
          if (r==atype<const  FormBilinear *>() ){
            const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      
            BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
            if (Op == NULL) {
              if(mpirank == 0) cout << "dynamic_cast error" << endl; 
            ffassert(0);
            }
      
            size_t Opsize= Op->v.size();
      
            KN<size_t> index_operator_finc(Opsize);
            KN<int>    new_index_funct_finc(Opsize);

            KN<size_t> index_operator_ftest(Opsize);
            KN<int>    new_index_funct_ftest(Opsize);

            // index to check if a integral is defined on multi block
            for(size_t jj=0; jj<Opsize; jj++){
              // attention la fonction test donne la ligne
              //  et la fonction test est en second
              BilinearOperator::K ll = Op->v[jj];
              pair<int,int> finc(ll.first.first), ftest(ll.first.second);

              long jj2= jj;

              index_operator_finc[ jj2] = jj;
              new_index_funct_finc[ jj2] = beginBlockUh[i]+finc.first;
            
              index_operator_ftest[ jj2]  = jj;
              new_index_funct_ftest[ jj2] = beginBlockVh[j]+ftest.first; 
            }
            changeIndexFunctionInconnue(*Op, index_operator_finc, new_index_funct_finc );
      
            changeIndexFunctionTest(*Op, index_operator_ftest, new_index_funct_ftest  );          
          }  
      
#ifndef FFLANG
#if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
          // ******************************************
          // Case BemKFormBilinear (KERNEL FORM ONLY)
          // ******************************************
          else if (r==atype<const BemFormBilinear *>() ){
            BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
            int VVFBEM = bbtmp->type;

            if(VVFBEM ==1){
              // BemKFormBilinear * bb=new BemKFormBilinear(*dynamic_cast<const BemKFormBilinear *>(e));
              const BemKFormBilinear * bb= dynamic_cast<const BemKFormBilinear *>(e);
              FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
              if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; ffassert(0);}

              // loop over the index of finconnue
              LOperaG * OpG = const_cast<LOperaG *>(b->fi);
              ffassert( OpG->v.size() == 1);
              size_t Opsize= OpG->v.size();

              for(size_t jjj=0; jjj<Opsize; jjj++){
                LOperaG::K *lf=&(OpG->v[jjj]);
                OpG->v[jjj].first.first += beginBlockUh[i];
                
              }

              // Loop over the index of ftest
              LOperaD * OpD = const_cast<LOperaD *>(b->ft);
              ffassert( OpD->v.size() == 1);
              Opsize= OpD->v.size();
              for(size_t jjj=0; jjj<Opsize; jjj++){
                LOperaD::K *lf=&(OpD->v[jjj]);
                OpD->v[jjj].first.first += beginBlockVh[j]; 
                
              }
            }
            else if(VVFBEM == 2){
              BemPFormBilinear * bb=new BemPFormBilinear(*dynamic_cast<const BemPFormBilinear *>(e));
              FoperatorPBEM * b=const_cast<  FoperatorPBEM *>(bb->b);
              if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; ffassert(0);}

              cerr << " BEM Potential in composite FESpace in construction " << endl;
              ffassert(0);
            }

          }
#endif
#endif  
          // BC_set
          // case BC_set 
          else if(r == atype<const  BC_set  *>()){
            ffassert( i == j ); // diagonal block 
            BC_set * bc=dynamic_cast<BC_set *>(e);
        
            //KN<int>  new_index_funct_finc( bc.size() );
            int kk=bc->bc.size();
            for (int k=0;k<kk;k++)
            {
              //bc->bc[k].first += beginBlockUh[i];
              pair<int,Expression> &xx=bc->bc[k];
              xx.first += beginBlockUh[i];
            }
          }
        }
      }
  }
}

/*
template<class FESpaceT1,class FESpaceT2>
MatriceMorse<R> * buildInterpolationMatrixT(const FESpaceT1 & Uh,const FESpaceT2 & Vh,void *data);

template< >
MatriceMorse<R> * buildInterpolationMatrixT<FESpaceL,FESpace>(const FESpaceL & Uh,const FESpace & Vh,void *data);
*/

template< class R, class FESpaceT1, class FESpaceT2 >
Matrice_Creuse<R> *  buildMatrixInterpolationForCompositeFESpace(const FESpaceT1 * Uh ,const FESpaceT2 * Vh){
  ffassert(Uh);
  ffassert(Vh);
  int NUh = Uh->N;
  int NVh = Vh->N;
  
  Matrice_Creuse<R> * sparse_mat= new Matrice_Creuse<R>();

  // Remarque pas de U2Vc pour l'instant
  int* data = new int[4 + NUh];
  // default value for the interpolation matrix
  data[0]=false;         // transpose not
  data[1]=(long) op_id;  // get just value
  data[2]=false;         // get just value
  data[3]=0L;            // get just value

  for(int i=0;i<NUh;++i) data[4+i]=i;//

  if(verbosity>3){
    for(int i=0;i<NUh;++i)
    {
      cout << "The Uh componante " << i << " -> " << data[4+i] << "  Componante of Vh  " <<endl;
    }
  }
  for(int i=0;i<NUh;++i){
    if(data[4+i]>=NVh)
    {
      cout << "The Uh componante " << i << " -> " << data[4+i] << " >= " << NVh << " number of Vh Componante " <<endl;
      ExecError("Interpolation incompability between componante ");
    }
  }
  const FESpaceT1 &rUh = *Uh;
  const FESpaceT2 &rVh = *Vh;

  MatriceMorse<R>* titi=buildInterpolationMatrixT<FESpaceT1,FESpaceT2>(rUh,rVh,data);

  sparse_mat->init();
  sparse_mat->typemat=0;//(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparse_mat->A.master( titi );	  //  sparse_mat->A.master(new MatriceMorse<R>(*Uh,*Vh,buildInterpolationMatrix,data));
  if(verbosity>3){
    cout << "sparse_mat->typemat=" << sparse_mat->typemat << endl;
    cout << "N=" << sparse_mat->A->n << endl;
    cout << "M=" << sparse_mat->A->m << endl;
  }
  delete [] data;

  return sparse_mat;
}

// print the information of the list of arguments of the varf 
void listOfComponentBilinearForm(const list<C_F0> & largs){

  list<C_F0>::const_iterator ii,ib=largs.begin(),ie=largs.end(); 

  // loop over largs information 
  cout << "loop over the integral" << endl;

  int count_integral = 0;
  for (ii=ib;ii != ie;ii++) {
    count_integral++;
    Expression e=ii->LeftValue();
    aType r = ii->left();
    
    cout <<"========================================================" << endl;
    cout <<"=                                                      =" << endl;
    cout << "reading the " << count_integral << "-th integral" << endl;
    cout << "e=" << e << ", " << "r=" << r << endl;
    cout <<"=                                                      =" << endl;
    
    // ***************************************
    // Case FormBillinear
    // ***************************************
    if (r==atype<const  FormBilinear *>() ){
      const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      const CDomainOfIntegration & di= *bb->di;
      
      cout << "di.kind=" << di.kind << endl;
      cout << "di.dHat=" << di.dHat << endl;
      cout << "di.d=" << di.d << endl;
      cout << "di.Th=" << di.Th << endl;
      
      int    d = di.d;
      int dHat = di.dHat;
      
      BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
      if (Op == NULL) {
        if(mpirank == 0) cout << "dynamic_cast error" << endl; 
        ffassert(0);
      }
      
      size_t Opsize= Op->v.size();
      
        cout << " loop over the term inside the integral" << endl;
        cout << " Number of term in the integral:: Op->v.size()=" << Op->v.size() << endl;
      
      // index to check if a integral is defined on multi block

      for(size_t jj=0; jj<Opsize; jj++){
        // attention la fonction test donne la ligne
        //  et la fonction test est en second
        BilinearOperator::K ll = Op->v[jj];
        pair<int,int> finc(ll.first.first), ftest(ll.first.second);
        
        cout << " operateur jj= " << jj << endl;
        cout << " FormBilinear: number of unknown finc=" <<  finc.first << " ,ftest= " << ftest.first << endl;
        cout << " FormBilinear: operator order finc   =" << finc.second << " ,ftest= " << ftest.second << endl; // ordre   only op_id=0
        
        // finc.first : index de component de la fonction inconnue
        // ftest.first: index de component de la fonction test
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // finc.second : renvoie l'index du type de l'operateur: Id, dx(), dy(), dz(), dxx(), dxy(), ...
        //
        // la liste des index des operateurs est definis dans [[????]].
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // exemple vf([u1,u2,..,u30],[v1,v2,...,v30]) = int2d(Th)(dx(u20)*v15)
        //      finc.first  = 20 , ftest.first = 15
        //      finc.second = 1 , ftest.second = 0

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
      }
    }
    else if(r == atype<const  BC_set  *>()){
      cout << " BC in variational form " << endl;
      const BC_set * bc=dynamic_cast<const  BC_set *>(e);
  
      int kk=bc->bc.size();
      cout << "bc.size=" << bc->bc.size() << endl;
      for (int k=0;k<kk;k++)
      {
        cout << "bc->bc["<< k << "].first= " << bc->bc[k].first << endl; 
      }
    }
    #ifndef FFLANG
    #if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
    // ******************************************
    // Case BemKFormBilinear (KERNEL FORM ONLY)
    // ******************************************
    else if (r==atype<const BemFormBilinear *>() ){
      BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
      int VVFBEM = bbtmp->type;
      cout << " read index term = "<< count_integral << " VVFBEM=" << VVFBEM << endl;
      if(VVFBEM ==1){
        //BemKFormBilinear * bb=new BemKFormBilinear(*dynamic_cast<const BemKFormBilinear *>(e));
        BemKFormBilinear * bb = dynamic_cast< BemKFormBilinear *>(e);
        FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
        if (b == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; ffassert(0);}

        // loop over the index of finconnue
        LOperaG * OpG = const_cast<LOperaG *>(b->fi);
        ffassert( OpG->v.size() == 1);
        size_t Opsize= OpG->v.size();
        cout << " pointeur  " << OpG << endl;
        for(size_t jjj=0; jjj<Opsize; jjj++){
          LOperaG::K *lf=&(OpG->v[jjj]);
          pair<int,int> finc(lf->first);
          cout << " operateur jj= " << jjj << endl;
          cout << " BemFormBilinear: number of unknown finc = " << finc.first << endl;
          cout << " BemFormBilinear: operator order   finc = " << finc.second << endl; 
        }


        // Loop over the index of ftest
        LOperaD * OpD = const_cast<LOperaD *>(b->ft);
        ffassert( OpD->v.size() == 1);
        Opsize= OpD->v.size();
        for(size_t jjj=0; jjj<Opsize; jjj++){
          LOperaD::K *lf=&(OpD->v[jjj]);
          pair<int,int> ftest(lf->first);
          cout << " operateur jj= " << jjj << endl;
          cout << " BemormBilinear: number of unknown ftest = " << ftest.first << endl;
          cout << " BemFormBilinear: operator order   ftest = " << ftest.second << endl; 
        }
      }
      else if(VVFBEM == 2){
        cerr << " BEM Potential in composite FESpace in construction " << endl;
        ffassert(0);
      }
      else{
        cerr << " VFBEM=1 (kernel) or VFBEM=2 (potential) " << endl;
        ffassert(0);
      }
    }
    #endif
    #endif
    else{
      cerr << "listOfComponentBilinearForm :: vectorial FESpace :: bilinear form only " << endl;
      cerr << "      uniquement terme bilineaire + BC " << endl;
      ffassert(0);
    }
  }

}

/**
       *  @brief  determine if we have BEM bilinear operator in a subblock and the type
       *  @param  largs list of argument of the Bilinear Form
       */

/*
  This function give the good result if only if the FESpace inconnu and  FESpace test are scalar FESpace
  due to this check : 
          if( finc.first==0 && ftest.first==0)      // only first component for finc and ftest
*/

int haveBemSubMatrixBlock(const list<C_F0> & largs, int Uh_NbItem, int Vh_NbItem){

  ffassert( Uh_NbItem == 1 && Vh_NbItem == 1);

  // this function is used to know if a block of the matrix contains BEM Operator
  // return the type of Bem Block
  bool haveBemFormBilinear   = false;
  bool haveMassMatrixforBEMTOOL = false; 
  int nbFB=0; 
  int nbBC=0;
  int nbBEM=0;
  // 
  list<C_F0>::const_iterator ii, b_ib=largs.begin(), b_ie=largs.end(); 
  
#ifndef FFLANG
#if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
  for( ii=b_ib;ii != b_ie;ii++){
    Expression e=ii->LeftValue();
    aType r = ii->left();

    // ******************************************
    // Case BemKFormBilinear (KERNEL FORM ONLY)
    // ******************************************
    if (r==atype<const BemFormBilinear *>() ){
      haveBemFormBilinear = true;
      nbBEM++;
    }
  }
#endif
#endif
  if( nbBEM > 1 ){
    cerr << "Two BEM operator in a sub-matrix defined by a time product of two FESpace is not allowed" << endl;
    ffassert(0);
  }

  for( ii=b_ib;ii != b_ie;ii++){
    Expression e=ii->LeftValue();
    aType r = ii->left();
  
    // ***************************************
    // Case FormBillinear
    // ***************************************
    if (r==atype<const  FormBilinear *>() ){
      nbFB++;
      //
      const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      const CDomainOfIntegration & di= *bb->di;
      // check the integration (keyword)

      BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
      if (Op == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; }

      size_t Opsize= Op->v.size();

      for(size_t jj=0; jj<Opsize; jj++){
        // attention la fonction test donne la ligne
        //  et la fonction test est en second
        BilinearOperator::K ll = Op->v[jj]; //  LinearComb<pair<MGauche,MDroit>,C_F0> BilinearOperator;
        pair<int,int> finc(ll.first.first), ftest(ll.first.second);

        // check if we have a mass matrix that we need to add to BEMTOOL 
        // cf bem.hpp: see pair<BemKernel*,double> getBemKernel(Stack stack, const list<C_F0> & largs).
  
        if( finc.first==0 && ftest.first==0)      // only first component for finc and ftest
          if( finc.second==0 && ftest.second==0 ) // only op_id=0
            if( (di.kind == CDomainOfIntegration::int1d && di.dHat==1) || (di.kind == CDomainOfIntegration::int2d && di.dHat==2) ){
              if( ! haveMassMatrixforBEMTOOL ){
                cerr << " Two mass matrix for BEMTOOl in largs is not allowed " << endl;
                ffassert(0);
              }
              haveMassMatrixforBEMTOOL = true;
            }
      }
    } // atype FormBilinear
    else if(r == atype<const  BC_set  *>()){
      nbBC++;
    } // atype BC
  }

  if( nbBC >0 ){
    if( haveBemFormBilinear ){
      cerr << "Error: Dirichlet Boundary condition and BEM operator in the same variational form." << endl;
      ffassert(0);
    }
  }

  // mixed FEM-BEM terms
  if( haveBemFormBilinear && haveMassMatrixforBEMTOOL  && nbFB > 0 ){
    return 12;    // BEM + mass matrix ==> H-matrix, FEM
  }
  else if( haveBemFormBilinear && nbFB > 0){
    return 11;    // BEM+FEM
  }
  else if( haveBemFormBilinear && haveMassMatrixforBEMTOOL ){
      return 2;   // BEM + mass matrix ==> H-matrix
  }
  else if( haveBemFormBilinear ){
    return 1;    // BEM only
  }
  else{
    return 0;    // FEM only
  }
}

/**
       *  @brief  largs separate in two part :  BEM (H-matrix) and FEM 
       *  @param  largs list of argument of the Bilinear Form
       *  @param  largs_FEM list of argument for the FEM part
       *  @param  largs_BEM list of argument for the BEM part (included sometimes mass matrix ) that be compressed in H-matrix
       */

/*
  This function must be call if we have BEM Bilinear operator in a block.
  
*/
void separateFEMpartBemPart(const list<C_F0> & largs, list<C_F0> &largs_FEM, list<C_F0> &largs_BEM ){
  
  bool newC_F0 = false;
  bool haveBemFormBilinear = false;
  int nbBEM=0;
  // 
  list<C_F0>::const_iterator ii, b_ib=largs.begin(), b_ie=largs.end(); 
  
#ifndef FFLANG
#if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
  for( ii=b_ib;ii != b_ie;ii++){
    Expression e=ii->LeftValue();
    aType r = ii->left();

    // ******************************************
    // Case BemKFormBilinear (KERNEL FORM ONLY)
    // ******************************************
    if (r==atype<const BemFormBilinear *>() ){
      BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
      ffassert(bbtmp);
      int VVFBEM = bbtmp->type;
      haveBemFormBilinear = true;
      nbBEM++;
      if(VVFBEM ==1){
        BemKFormBilinear * bb= dynamic_cast< BemKFormBilinear *>(e);
        ffassert(bb);

        //BemKFormBilinear * bbnew = new BemKFormBilinear( bb->di, FoperatorKBEM(bb->b->kbem, *(bb->b->fi), *(bb->b->ft) ) ); // marche ???
        if(newC_F0){
            BemKFormBilinear * bbnew = new BemKFormBilinear( bb->di, FoperatorKBEM(bb->b->kbem, *(bb->b->fi), *(bb->b->ft) ) ); // marche ???
            largs_BEM.push_back( C_F0( bbnew, r ) ); 
        }
        else largs_BEM.push_back(*ii);

      }else{
        cerr << "case VVFBEM noot coded yet. todo." << endl;
        ffassert(0);
      }  
      // On pourrait retourner l'element de largs. Mais on prefere en creer un nouveau pour etre coherent avec cree deleteNewLargs.
      // largs_BEM.push_back(*ii); 
    }
  }
  if(nbBEM >1){
    cerr << "Two BEM operator in a sub-matrix defined by a time product of the same two FESpace is not allowed" << endl;
    ffassert(0);
  }
#endif
#endif

  for( ii=b_ib;ii != b_ie;ii++){
    Expression e=ii->LeftValue();
    aType r = ii->left();
  
    // ***************************************
    // Case FormBillinear
    // ***************************************
    if (r==atype<const  FormBilinear *>() ){
      
      //
      const FormBilinear * bb=dynamic_cast<const  FormBilinear *>(e);
      const CDomainOfIntegration & di= *bb->di;
      // check the integration (keyword)

      BilinearOperator * Op=const_cast<  BilinearOperator *>(bb->b);
      if (Op == NULL) { if(mpirank == 0) cout << "dynamic_cast error" << endl; }

      size_t Opsize= Op->v.size();
      
      bool  haveBEMmass   = false;
      size_t indexBEMmass = -1;

      for(size_t jj=0; jj<Opsize; jj++){
        // attention la fonction test donne la ligne
        //  et la fonction test est en second
        BilinearOperator::K ll = Op->v[jj]; //  LinearComb<pair<MGauche,MDroit>,C_F0> BilinearOperator;
        pair<int,int> finc(ll.first.first), ftest(ll.first.second);

        // check if we have a mass matrix that we need to add to BEMTOOL 
        // cf bem.hpp: see pair<BemKernel*,double> getBemKernel(Stack stack, const list<C_F0> & largs).

        if( haveBemFormBilinear ){
          // IF BEM VERIFIED THAT WE HAVE A MASS MATRIX
          if( finc.first==0 && ftest.first==0)      // only first component for finc and ftest
            if( finc.second==0 && ftest.second==0 ) // only op_id=0
              if( (di.kind == CDomainOfIntegration::int1d && di.dHat==1) || (di.kind == CDomainOfIntegration::int2d && di.dHat==2) ){   
                
                //haveBEMmass  = true;
                haveBEMmass   = false;
                indexBEMmass  = jj;
              }
        }
      }

      if( !haveBEMmass ){
        if(newC_F0) largs_FEM.push_back( C_F0( new FormBilinear( &di, Op ), r ) );
        else largs_FEM.push_back(*ii) ;
      }
      else{
        if(indexBEMmass <0){ ffassert(0);}
        if(Opsize==1 && indexBEMmass == 0){
          // case one element corresponding the mass matrix adding for BEMTOOOL
          BilinearOperator * OpBEM = new BilinearOperator( Op->v[indexBEMmass].first, Op->v[indexBEMmass].second );
          if(newC_F0){
            BilinearOperator * OpBEM = new BilinearOperator( Op->v[indexBEMmass].first, Op->v[indexBEMmass].second );
            largs_BEM.push_back( C_F0( new FormBilinear( &di, OpBEM ), r ) );
          }
          else largs_BEM.push_back( *ii );
        }
        else{
          ffassert( newC_F0 ); // True only if we have newC_F0
          bool * partFEM = new bool[Opsize];
          for(size_t jj=0; jj<Opsize; jj++){
            partFEM[jj]= true;
          }
          partFEM[indexBEMmass] = false;

          // check 
          if(verbosity>3) cout << "partFEM=" << partFEM << endl;
          
          if(newC_F0){
            BilinearOperator * OpFEM = new BilinearOperator( *Op, partFEM );
            largs_FEM.push_back( C_F0( new FormBilinear( &di, OpFEM ), r ) );
          }
          else largs_FEM.push_back(*ii);

          if(newC_F0){
            BilinearOperator * OpBEM = new BilinearOperator( Op->v[indexBEMmass].first, Op->v[indexBEMmass].second );
            largs_BEM.push_back( C_F0( new FormBilinear( &di, OpBEM ), r ) );
          }
          else largs_BEM.push_back(*ii);

          delete [] partFEM;
        }
      }
      // creation des deux listes

    } // atype FormBilinear
    else if(r == atype<const  BC_set  *>()){
      const BC_set * bc=dynamic_cast<const BC_set *>(e);
      if(newC_F0) largs_FEM.push_back( C_F0( new BC_set(*bc), r) ); 
      else largs_FEM.push_back( *ii ); 
    } // atype BC
    else if(r == atype<const  FormLinear  *>()){
      const FormLinear * ll=dynamic_cast<const  FormLinear *>(e);
      LOperaD * Op = const_cast<LOperaD *>(ll->l);
      if(newC_F0) largs_FEM.push_back( C_F0( new FormLinear( (ll->di), Op ), r ) ); 
      else largs_FEM.push_back( *ii ); 

    } // atype FormLinear
  }
}

/**
       *  @brief  Function to delete element of newlargs = list<C_F0> to avoid memory leak.
       *  @param  largs list of argument of the composite FESpace gene
       */

/* remark: only FormBilinear  and BemFormBilinear is deleted */
void deleteNewLargs(list<C_F0> &newlargs){
  list<C_F0>::iterator ii,ib=newlargs.begin(),ie=newlargs.end(); 
  for (ii=ib;ii != ie;ii++) {
    Expression e=ii->LeftValue();
    aType r = ii->left();
  
    // ***************************************
    // Case FormBillinear
    // ***************************************
    if (r==atype<const  FormBilinear *>() ){
      FormBilinear * bb=dynamic_cast< FormBilinear *>(e);
      //BilinearOperator * Op=dynamic_cast<  BilinearOperator *>(bb->b);
      for(size_t jj=0; jj < bb->b->v.size(); jj++ ){
        bb->b->v[jj].second.Destroy();
      }
      delete bb->b;
      delete bb;
    }
    if (r==atype<const  FormLinear *>() ){
      FormLinear * ll=dynamic_cast< FormLinear *>(e);
      LOperaD * Op = const_cast<LOperaD *>(ll->l);
      if (Op == NULL) {
        if(mpirank == 0) cout << "dynamic_cast error" << endl; 
        ffassert(0);
      }
      
      for(size_t jj=0; jj < Op->v.size(); jj++ ){
        Op->v[jj].second.Destroy();
      }
      delete ll->l;
      delete ll;
    }
    // ****************************************************************
    // BC_set in newlargs correspond to the original BC_set in largs
    //   ==> implies not deleted this element.
    // ****************************************************************
    else if (r==atype<const  BC_set *>() ){
      BC_set * bb=dynamic_cast< BC_set *>(e);
      /*
      //BilinearOperator * Op=dynamic_cast<  BilinearOperator *>(bb->b);
      for(size_t jj=0; jj < bb->b->v.size(); jj++ ){
        bb->b->v[jj].second.Destroy();
      }
      delete bb->b;
      */
      delete bb;
    }
#ifndef FFLANG
#if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
    // ******************************************
    // Case BemKFormBilinear (KERNEL FORM ONLY)
    // ******************************************
    else if (r==atype<const BemFormBilinear *>() ){
      BemFormBilinear * bbtmp= dynamic_cast< BemFormBilinear *>(e);
      ffassert(bbtmp);
      int VVFBEM = bbtmp->type;

      if(VVFBEM ==1){
        BemKFormBilinear * bb= dynamic_cast< BemKFormBilinear *>(e);
        ffassert(bb);
        FoperatorKBEM * b=const_cast<  FoperatorKBEM *>(bb->b);
        LOperaG * OpG = const_cast<LOperaG *>(b->fi);
        for(size_t jj=0; jj < OpG->v.size(); jj++ ){
          OpG->v[jj].second.Destroy();
        }
        LOperaD * OpD = const_cast<LOperaD *>(b->ft);
        for(size_t jj=0; jj < OpD->v.size(); jj++ ){
          OpD->v[jj].second.Destroy();
        }
        delete b->fi;
        delete b->ft;
        delete b;
        delete bb;
        //BemKFormBilinear * bbnew = new BemKFormBilinear( bb->di, FoperatorKBEM(bb->b->kbem, *(bb->b->fi), *(bb->b->ft) ) ); // marche ???
        //newlargs.push_back( C_F0( bbnew, r ) );
      }
      else{
        cerr << "coding only case VFBEM == 1. Other case todo." <<endl;
        ffassert(0);
      }
    }
#endif
#endif

    ii->Destroy();
  }
  newlargs.clear();
}


/**
 * @brief : A décrire 
 * 
 * @param :
 * 
 * */

/* ATTENTION CETTE FONCTION SUPPOSE QUE L'ALLOCATION DE LA MATRICE EST FAIT AVANT */
template<class K,class MMesh, class v_fes1, class v_fes2>
void varfToCompositeBlockLinearSystem_fes(bool initmat, bool initx, 
        /*const typename v_fes1::pfes &pUh, const typename v_fes2::pfes &pVh,*/
        const typename v_fes1::FESpace*& PUh, const typename v_fes2::FESpace*& PVh,
        const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
        KN_<K> *B, KN_<K> *X, MatriceCreuse<K> &A,int *mpirankandsize, bool B_from_varf=false){

  typedef typename  v_fes1::pfes pfes1;
  typedef typename  v_fes2::pfes pfes2;
  typedef typename  v_fes1::FESpace FESpace1;
  typedef typename  v_fes2::FESpace FESpace2;
  typedef typename  FESpace1::Mesh Mesh1;
  typedef typename  FESpace2::Mesh Mesh2;

  varfToCompositeBlockLinearSystem<K,MMesh, FESpace1, FESpace2>( initmat, initx, PUh, PVh, sym, tgv, largs, stack, B, X, A,mpirankandsize,B_from_varf);

}

template<class R,class MMesh, class FESpace1, class FESpace2>
void varfToCompositeBlockLinearSystem(bool initmat, bool initx, const FESpace1 * PUh, const FESpace2 * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<R> *B, KN_<R> *X, MatriceCreuse<R> &A,int *mpirankandsize, bool B_from_varf)
                              {
  typedef typename  FESpace1::Mesh Mesh1;
  typedef typename  FESpace2::Mesh Mesh2;

  // this lines must be defined outside this function

  const FESpace1 & Uh =  *PUh ;
  const FESpace2 & Vh =  *PVh ;
  const MMesh* pTh = (is_same< Mesh1, Mesh2 >::value) ? (MMesh*)&PUh->Th : 0;
  const MMesh &Th= *pTh ;    // integration Th
  bool same=isSameMesh( largs, &Uh.Th, &Vh.Th, stack);


  // PAC(e) :  on retourne le type MatriceCreuse Ici et non Matrice_Creuse<R> dans creationBlockOfMatrixToBilinearForm.
  // Attention pour la generalisation
  if(same){
    if  (AssembleVarForm<R,MatriceCreuse<R>,MMesh,FESpace1,FESpace2 >( stack,Th,Uh,Vh,sym, initmat ? &A:0 , B, largs, mpirankandsize))
    {
      if( B && !B_from_varf ){
        // case problem or solve : b = int2d(Th)( f*v ) but we evaluate before  -int2d(Th)( f*v ) in AssembleVarForm
        *B = - *B; 
        // hach FH
        for (int i=0, n= B->N(); i< n; i++)
        if( abs((*B)[i]) < 1.e-60 ) (*B)[i]=0;
      }
      AssembleBC<R,MMesh,FESpace1,FESpace2> ( stack,Th,Uh,Vh,sym, initmat ? &A:0 , B, initx ? X:0,  largs, tgv, mpirankandsize);   // TODO with problem
    }
    else{
      if( B && !B_from_varf ) *B = - *B; // case problem or solve : b = int2d(Th)( f*v ) but we evaluate before  -int2d(Th)( f*v ) in AssembleVarForm
    }
  }else{
#ifdef V3__CODE
    MatriceMap<R>   AAA;
    cout << "V3__CODE=" << AAA.size() << endl;
    ffassert(0); // code a faire
    MatriceMorse<R> *pMA =   new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,AAA.size(),sym>0);
    bool bc=AssembleVarForm<R,MatriceMap<R>,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,sym>0,initmat ? &AAA:0,B,largs);
    pMA->addMap(1.,AAA);
#else
    MatriceMorse<R> *pMA =  dynamic_cast< HashMatrix<int,R>*>(&A);// new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,0,sym);
    MatriceMap<R>  &  AAA = *pMA;
    bool bc=AssembleVarForm<R,MatriceMap<R>,MMesh,FESpace1,FESpace2>( stack,Th,Uh,Vh,sym>0,initmat ? &AAA:0,B,largs, mpirankandsize);
#endif
    if (bc){
      if( B && !B_from_varf ){
        // case problem or solve : b = int2d(Th)( f*v ) but we evaluate before  -int2d(Th)( f*v ) in AssembleVarForm
        *B = - *B;
        // hach FH
        for (int i=0, n= B->N(); i< n; i++)
        if( abs((*B)[i]) < 1.e-60 ) (*B)[i]=0;
      }
      AssembleBC<R> ( stack,Th,Uh,Vh,sym, initmat ? &A:0 , B, initx ? X:0,  largs, tgv, mpirankandsize);   // TODO with problem
    }
    else{
      if( B && !B_from_varf ) *B = - *B; // case problem or solve : b = int2d(Th)( f*v ) but we evaluate before  -int2d(Th)( f*v ) in AssembleVarForm
    }
  }
}

#ifndef FFLANG
#if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)

/**
 * @brief
 * 
 * @param HMatrix_block  H-matrix corresponding to a local block 
 */

/* pas de initmat */
template< class R>
void varfBemToCompositeBlockLinearSystem_hmat(const int& iUh, const int &jVh, 
                const int &typeUh, const int &typeVh, const generic_v_fes * LLUh, const generic_v_fes * LLVh, 
                const list<C_F0> & b_largs_zz, Stack stack, Data_Bem_Solver &dsbem,
                HMatrixVirt<R>** Hmat){

    //cout << "iUh=" << iUh << " , " << "jVh=" << jVh << endl;
    //cout << LLUh << " " <<  LLVh << endl;

    // reinitialise Hmat
    if( *Hmat)
    delete *Hmat;
    *Hmat =0;

    int VFBEM = typeVFBEM(b_largs_zz,stack); // get type of VFBEM



    // block diagonal matrix
    if( typeUh == 4 && typeVh == 4 ){
      ffassert( iUh==jVh ); // If not a block diagonal not coded yet :: Same mesh for the different FESpace
      // MeshS --- MeshS
      // ==== FESpace 3d Surf: inconnue et test ===
      const FESpaceS * PUh = (FESpaceS *) LLUh->getpVh();
      const FESpaceS * PVh = (FESpaceS *) LLVh->getpVh();

      creationHMatrixtoBEMForm<R, MeshS, FESpaceS, FESpaceS>(PUh, PVh, VFBEM, 
                      b_largs_zz, stack, dsbem, Hmat);

    }
    else if( typeUh == 5 && typeVh == 5  ){
      ffassert( iUh==jVh ); // If not a block diagonal not coded yet :: Same mesh for the different FESpace
      // MeshL --- MeshL
      // ==== FESpace 3d Curve: inconnue et test ===
      const FESpaceL * PUh = (FESpaceL *) LLUh->getpVh();
      const FESpaceL * PVh = (FESpaceL *) LLVh->getpVh();

      creationHMatrixtoBEMForm<R, MeshL, FESpaceL, FESpaceL> ( PUh, PVh, VFBEM, 
                      b_largs_zz, stack, dsbem, Hmat );
    }
    else if( typeUh == 5 && typeVh == 2 ){
      //ffassert( !(iUh==jVh) );
      // case Uh[i] == MeshL et Vh[j] = Mesh2
                    
      if(verbosity>3) cout << " === creation de la matrice BEM pour un bloc non diagonaux === " << endl;
      if(verbosity>3) cout << " ===      hypothesis:: FESpace become FESpaceL for Vh      === " << endl;
      
      const FESpaceL * PUh = (FESpaceL *) LLUh->getpVh();
      creationHMatrixtoBEMForm<R, MeshL, FESpaceL, FESpaceL> ( PUh, PUh, VFBEM, 
                        b_largs_zz, stack, dsbem, Hmat );
    }
    else if( typeUh == 4 && typeVh == 3 ){
      //ffassert( !(iUh==jVh) );
      // case Uh[i] == MeshS et Vh[j] = Mesh3
   
      if(verbosity>3) cout << " === creation de la matrice BEM pour un bloc non diagonaux === " << endl;
      if(verbosity>3) cout << " ===      hypothesis:: FESpace3 become FESpaceS for Vh      === " << endl;
      const FESpaceS * PUh = (FESpaceS *) LLUh->getpVh();
      creationHMatrixtoBEMForm<R, MeshS, FESpaceS, FESpaceS> ( PUh, PUh, VFBEM, 
                        b_largs_zz, stack, dsbem, Hmat );
    }
    else if( typeUh == 2 && typeVh == 5 ){
      // case Uh[i] == Mesh2 et Vh[j] = MeshL
      if(verbosity>3) cout << " === creation de la matrice BEM pour un bloc non diagonaux === " << endl;
      if(verbosity>3) cout << " ===      hypothesis:: FESpace become FESpaceL for Uh      === " << endl;

      const FESpaceL * PVh = (FESpaceL *) LLVh->getpVh();
      creationHMatrixtoBEMForm<R, MeshL, FESpaceL, FESpaceL> ( PVh, PVh, VFBEM,
                        b_largs_zz, stack, dsbem, Hmat );
    }
    else if( typeUh == 3 && typeVh == 4 ){
      // case Uh[i] == Mesh3 et Vh[j] = MeshS
      if(verbosity>3) cout << " === creation de la matrice BEM pour un bloc non diagonaux === " << endl;
      if(verbosity>3) cout << " ===      hypothesis:: FESpace3 become FESpaceS for Uh      === " << endl;
      const FESpaceS * PVh = (FESpaceS *) LLVh->getpVh();
      creationHMatrixtoBEMForm<R, MeshS, FESpaceS, FESpaceS> ( PVh, PVh, VFBEM,
                        b_largs_zz, stack, dsbem, Hmat );
    }
    else{
      cerr << " BEM bilinear form " << endl;
      cerr << " Block ("<< iUh <<" ,"<< jVh << ")" << endl;
      cerr << " =: Pas prise en compte des FESpace inconnue de type := "<< typeFEtoString( typeUh ) << endl;
      cerr << " =:                 avec des FESpace test de type    := "<< typeFEtoString( typeVh ) << endl;
      ffassert(0);
    }
}

/**
 * @brief
 * 
 * @param hm_A_block  matrix creuse corresponding to a local block 
 */
template< class R>
HashMatrix<int,R> * varfBemToBlockDenseMatrix(const int& iUh, const int &jVh, const int &typeUh, const int &typeVh, const generic_v_fes * LLUh, const generic_v_fes * LLVh,
      const list<C_F0> & b_largs_zz, Stack stack, Data_Bem_Solver &dsbem){

  // creation of the H-matrix
  HMatrixVirt<R> ** Hmat = new HMatrixVirt<R> *();

  // construction of the H-matrix

  // IF LLVh is FESpace(Mesh2) ==> LLVh == LLUh (by hypothesis) 
  // Remove this comment in the future.

  //cout << "iUh=" << iUh << " , " << "jVh=" << jVh << endl;
  //cout << LLUh << " " <<  LLVh << endl;

  varfBemToCompositeBlockLinearSystem_hmat( iUh, jVh, typeUh, typeVh,
                LLUh, LLVh, b_largs_zz, stack, dsbem,Hmat);

  // creation of dense matrix
  KNM<R>* M= HMatrixVirtToDense< KNM<R>, R >(Hmat);
  HashMatrix<int,R> * hm_A_block = new HashMatrix<int,R>(*M);

  // delete dense matrix
  M->destroy();
  delete M; 

  // delete Hmat
  if( *Hmat)
    delete *Hmat;
  delete Hmat;

  return hm_A_block;
}

/**
 * @brief computation of matrice of interpolation for BEM 
 *       PAC(e) remove in the future when is not necessary
*/
Matrice_Creuse<double> * buildBlockMatrixInterpolation( const int &i, const int&j, const int&typeUh, const int& typeVh, 
                                             const generic_v_fes * LLUh, const generic_v_fes * LLVh ){
  
  ffassert( !(LLUh == LLVh) );

  if( typeUh == 5 && typeVh == 2 ){
    // case Uh[i] == MeshL et Vh[j] = Mesh2 
    const FESpaceL * PUh = (FESpaceL *) LLUh->getpVh();
    const FESpace * PVh = (FESpace *) LLVh->getpVh();
   
    Matrice_Creuse<double> *  MI_BBB = buildMatrixInterpolationForCompositeFESpace<double,FESpaceL,FESpace>( PUh, PVh );
    return MI_BBB;
    // return  MI_BBB->pHM(); // (HashMatrix return) 
  }
  else if( typeUh == 2 && typeVh == 5 ){
    // case Uh[i] == MeshL et Vh[j] = Mesh2
    const FESpace * PUh = (FESpace *) LLUh->getpVh();
    const FESpaceL * PVh = (FESpaceL *) LLVh->getpVh();

    Matrice_Creuse<double> *  MI_BBB = buildMatrixInterpolationForCompositeFESpace<double,FESpaceL,FESpace>( PVh, PUh );
    return MI_BBB;
  }
  else{
    cerr << "other case in construction" << endl;
    ffassert(0); 
  }
  return nullptr;
}

/**
 * @brief
 * 
 * @param hm_A  global matrix of a Linear System
 */

template< class R>
void varfBemToCompositeBlockLinearSystem(const int& i, const int &j, 
                const int &typeUh, const  int &typeVh,
                const long &sizeUh, const long &sizeVh,
                const long &offsetUh, const long &offsetVh,
                const generic_v_fes *LLUh, const generic_v_fes * LLVh,
                const list<C_F0> & b_largs_zz, Stack stack, Expression const * nargs,
                HashMatrix<int,R> *hm_A, const int &n_name_param){ //=OpCall_FormBilinear_np::n_name_param){

  Data_Bem_Solver dsbem;
  dsbem.factorize=0;
  dsbem.initmat=true;
  SetEnd_Data_Bem_Solver<R>(stack, dsbem, nargs,n_name_param);  // LIST_NAME_PARM_HMAT

  // hm_A_block : is a dense block
  HashMatrix<int,R> * hm_A_block = varfBemToBlockDenseMatrix<R>(i, j, typeUh, typeVh, LLUh, LLVh, b_largs_zz, stack, dsbem );

  // IF LLVh is FESpace(Mesh2) ==> LLVh become LLUh (by hypothesis) 
  // Remove this comment in the future.

  bool isNeedInterpolationMatrix = false;
  bool isNeedInterpolationMatrixLeft = false;

  if( typeUh == 5 && typeVh == 2){
    isNeedInterpolationMatrix = true;
    isNeedInterpolationMatrixLeft = true;
  }
  else if( typeUh == 2 && typeVh == 5){
    isNeedInterpolationMatrix = true;
    isNeedInterpolationMatrixLeft = false;
  }
  if( isNeedInterpolationMatrix ){
    if(verbosity>3) cout << " creation of interpolation matrix " << endl;
    // Creation of the Interpolation Matrix
    Matrice_Creuse<double> * MI = buildBlockMatrixInterpolation( i, j, typeUh, typeVh, LLUh, LLVh );
    MatriceMorse<double> *mMI = MI->pHM(); 

    if(isNeedInterpolationMatrixLeft ){

      // BBB=MI'*BBB;
      MatriceMorse<R> *mAB=new MatriceMorse<R>(MI->M(), hm_A_block->m,0,0);
      AddMul<int,double,R,R>(*mAB,*mMI,*hm_A_block,true,false);

      hm_A->Add(mAB, R(1.0), false, (int) offsetVh,(int) offsetUh); 
      // test function (Vh) are the line and inconnu function (Uh) are the column

      delete mAB;
    }
    else{ 

      // BBB=BBB*MI;
      MatriceMorse<R> *mAB=new MatriceMorse<R>(hm_A_block->m,MI->M(),0,0);
      AddMul<int,R,double,R>(*mAB,*hm_A_block,*mMI,false,false);

      hm_A->Add(mAB, R(1.0), false, (int) offsetVh,(int) offsetUh);
      // test function (Vh) are the line and inconnu function (Uh) are the column

      delete mAB;
    }

    MI->destroy();
    delete MI;
  }
  else{
    hm_A->Add(hm_A_block, R(1.0), false, (int) offsetVh,(int) offsetUh); 
    // test function (Vh) are the line and inconnu function (Uh) are the column
  }
  hm_A_block->destroy();
}
#endif
#endif


#if !defined(PARALLELE) || !defined(WITH_bemtool) || !defined(WITH_htool)

template< class R>
void varfBemToCompositeBlockLinearSystem(const int& i, const int &j, 
                const int &typeUh, const  int &typeVh,
                const long &sizeUh, const long &sizeVh,
                const long &offsetUh, const long &offsetVh,
                const generic_v_fes * LLUh, const generic_v_fes * LLVh,
                const list<C_F0> & b_largs_zz, Stack stack, Expression const * nargs,
                HashMatrix<int,R> *hm_A,const int &n_name_param){
                  #if !defined(WITH_bemtool) || !defined(WITH_htool)
                  cerr << "no BEM library" << endl;
                  #else
                  cerr << "BEM only valid in parallel" << endl;
                  #endif
                  ffassert(0);
                }
#endif


template<class K,class MMesh,class v_fes1, class v_fes2> 
void varfToCompositeBlockLinearSystemALLCASE_pfesT( const int& i, const int &j, 
                const long &offsetUh, const long &offsetVh,
                const typename  v_fes1::FESpace * &PUh, const typename v_fes2::FESpace * &PVh,
                bool initmat, bool initx, const int &sym, const double &tgv, 
                const list<C_F0> & b_largs_zz, Stack stack, 
                KN_<K> *B, KN_<K> *X, HashMatrix<int,K> *hm_A,int *mpirankandsize, bool B_from_varf=false){
     

    typedef typename  v_fes1::pfes pfes1;
    typedef typename  v_fes2::pfes pfes2;
    typedef typename  v_fes1::FESpace FESpace1;
    typedef typename  v_fes2::FESpace FESpace2;

    //const FESpace1 * PUh = (FESpace1*) pfesUh->getpVh(); // update the FESpace1
    //const FESpace2 * PVh = (FESpace2*) pfesVh->getpVh(); // update the FESpace2

    MatriceCreuse<K> * A_block = new MatriceMorse<K>( PVh->NbOfDF, PUh->NbOfDF, 0, sym );
    MatriceCreuse<K>  & BBB(*A_block);

    long N = PVh->NbOfDF;
    long M = PUh->NbOfDF;
    
    ffassert( B && X || (!B && !X) );
    
    if(B && X){

      KN<K> x_block(N);
      x_block=K(0.0);   // initiallise the block to zero ??? 

      KN<K> b_block(M);
      b_block=K(0.0);   // initiallise the block to zero ??? 

      if( i == j ){
          ffassert( PUh->NbOfDF == PVh->NbOfDF); // not coding yet : voir comment prendre les conditions aux limites.
          //ffassert( PUh == PVh); 

          // diagonal block
          // give initial value to X: if necessary
          if(initx){
              for(int ii=0; ii<PUh->NbOfDF; ii++)
                  x_block[ii] = (*X)[ii+offsetUh];
          }
          
          // give initial value to B: if necessary
          for(int ii=0; ii<PVh->NbOfDF; ii++)
            b_block[ii] = (*B)[ii+offsetVh];

          varfToCompositeBlockLinearSystem_fes<K, MMesh, v_fes1, v_fes2>( initmat, initx, PUh, PVh, sym, tgv, b_largs_zz, stack, 
                          &b_block, &x_block, BBB, mpirankandsize, B_from_varf);

          // update X: if necessary
          if(initx){
              for(int ii=0; ii<PUh->NbOfDF; ii++)
                  (*X)[ii+offsetUh] += x_block[ii]; // verification += A faire si on doit le faire ?????
          }

          // update B: if necessary
          for(int ii=0; ii<PVh->NbOfDF; ii++)
              (*B)[ii+offsetVh] += b_block[ii];   
      }
      else{
          // non diagonal block
          // no B : so B_from_varf can be default value
          varfToCompositeBlockLinearSystem_fes<K, MMesh, v_fes1, v_fes2>( initmat, false, PUh, PVh, sym, tgv, b_largs_zz, stack, 
                          0, 0, BBB,mpirankandsize); 
      }
    }
    else{
      // no B : so B_from_varf can be default value
      varfToCompositeBlockLinearSystem_fes<K, MMesh, v_fes1, v_fes2>( initmat, false, PUh, PVh, sym, tgv, b_largs_zz, stack, 
                        0, 0, BBB,mpirankandsize);
    }

    if( hm_A ){
      const HashMatrix<int,K> * hm_block = dynamic_cast<HashMatrix<int,K> *>( &BBB );
      
      if(mpirank ==0 && verbosity > 3) cout << "varfToCompositeBlockLinearSystem_fes::Add local block (" << i << "," << j <<")=> size nnz=" << hm_block->nnz << endl;

      hm_A->Add(hm_block, K(1.0), false, (int) offsetVh,(int) offsetUh); // test function (Vh) are the line and inconnu function (Uh) are the column
    }
    else{
      // only LinearForm and BC
      ffassert( initmat == false );
    }
    delete A_block;
}

template<class R>
void varfToCompositeBlockLinearSystemALLCASE_pfes( const int& i, const int &j, 
                const int &typeUh, const int &typeVh, 
                const long &offsetUh, const long &offsetVh,
                const generic_v_fes * pfesUh, const generic_v_fes * pfesVh,
                bool initmat, bool initx, const int &sym, const double &tgv, 
                const pcommworld &commworld, const list<C_F0> & b_largs_zz, Stack stack,
                KN_<R> *B, KN_<R> *X, HashMatrix<int,R> *hm_A, bool B_from_varf){

  int mpirankandsize[2];
  mpirankandsize[0] = 0;
  mpirankandsize[1] = 1;
  #ifndef FFLANG
  #if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
  #include <mpi.h>
  MPI_Comm comm;
  if(commworld)
      MPI_Comm_dup(*((MPI_Comm*)commworld), &comm);
  else
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_rank(comm, &mpirankandsize[0]);
  MPI_Comm_size(comm, &mpirankandsize[1]);
  if( verbosity >3 && mpirank ==0) cout << "PARALLEL :: mpisize=" << mpirankandsize[1] << endl;
  #endif
  #endif
  if(typeUh == 2 && typeVh == 2){        // Mesh2 -- Mesh2
    const FESpace * PUh = (FESpace*) pfesUh->getpVh(); // update the FESpace1
    const FESpace * PVh = (FESpace*) pfesVh->getpVh(); // update the FESpace2

    varfToCompositeBlockLinearSystemALLCASE_pfesT<R,Mesh,v_fes,v_fes>( i, j, offsetUh, offsetVh, PUh, PVh,
                initmat, initx, sym, tgv, b_largs_zz, stack, B, X, hm_A,mpirankandsize,B_from_varf);
      
  }else if(typeUh == 3 && typeVh == 3){   // Mesh3 -- Mesh3
    const FESpace3 * PUh = (FESpace3*) pfesUh->getpVh(); // update the FESpace1
    const FESpace3 * PVh = (FESpace3*) pfesVh->getpVh(); // update the FESpace2

    varfToCompositeBlockLinearSystemALLCASE_pfesT<R,Mesh3,v_fes3,v_fes3>( i, j, offsetUh, offsetVh, PUh, PVh, initmat, initx, sym, tgv, b_largs_zz, stack, B, X, hm_A,mpirankandsize,B_from_varf);
      
  }else if(typeUh == 4 && typeVh == 4){   // MeshS -- MeshS
    const FESpaceS * PUh = (FESpaceS*) pfesUh->getpVh(); // update the FESpace1
    const FESpaceS * PVh = (FESpaceS*) pfesVh->getpVh(); // update the FESpace2

    varfToCompositeBlockLinearSystemALLCASE_pfesT<R,MeshS,v_fesS,v_fesS>( i, j, offsetUh, offsetVh,
                PUh, PVh, initmat, initx, sym, tgv, b_largs_zz, stack, B, X, hm_A,mpirankandsize,B_from_varf);
      
  }
  else if(typeUh == 5 && typeVh == 5){    // MeshL -- MeshL
    const FESpaceL * PUh = (FESpaceL*) pfesUh->getpVh(); // update the FESpace1
    const FESpaceL * PVh = (FESpaceL*) pfesVh->getpVh(); // update the FESpace2

    varfToCompositeBlockLinearSystemALLCASE_pfesT<R,MeshL,v_fesL,v_fesL>( i, j, offsetUh, offsetVh,   PUh, PVh, initmat, initx, sym, tgv, b_largs_zz, stack, B, X, hm_A,mpirankandsize,B_from_varf);
      
  }
  else if(typeUh == 5 && typeVh == 2){    // MeshL -- Mesh
    const FESpaceL * PUh = (FESpaceL*) pfesUh->getpVh(); // update the FESpace1
    const FESpace * PVh = (FESpace*) pfesVh->getpVh(); // update the FESpace2

    varfToCompositeBlockLinearSystemALLCASE_pfesT<R,MeshL,v_fesL,v_fes>( i, j, offsetUh, offsetVh, PUh, PVh, initmat, initx, sym, tgv, b_largs_zz, stack, B, X, hm_A,mpirankandsize,B_from_varf);
      
  }
  else if(typeUh == 2 && typeVh == 5){    // Mesh -- MeshL
    const FESpace * PUh = (FESpace*) pfesUh->getpVh(); // update the FESpace1
    const FESpaceL * PVh = (FESpaceL*) pfesVh->getpVh(); // update the FESpace2

    varfToCompositeBlockLinearSystemALLCASE_pfesT<R,MeshL,v_fes,v_fesL>( i, j, offsetUh, offsetVh, PUh, PVh, initmat, initx, sym, tgv, b_largs_zz, stack, B, X, hm_A,mpirankandsize,B_from_varf);
  }
  else if(typeUh == 4 && typeVh == 3){    // MeshS -- Mesh3
    const FESpaceS * PUh = (FESpaceS*) pfesUh->getpVh(); // update the FESpace1
    const FESpace3 * PVh = (FESpace3*) pfesVh->getpVh(); // update the FESpace2

    varfToCompositeBlockLinearSystemALLCASE_pfesT<R,MeshS,v_fesS,v_fes3>( i, j, offsetUh, offsetVh, PUh, PVh, initmat, initx, sym, tgv, b_largs_zz, stack, B, X, hm_A,mpirankandsize,B_from_varf);
      
  }
  else if(typeUh == 3 && typeVh == 4){    // Mesh3 -- MeshS
    const FESpace3 * PUh = (FESpace3*) pfesUh->getpVh(); // update the FESpace1
    const FESpaceS * PVh = (FESpaceS*) pfesVh->getpVh(); // update the FESpace2

    varfToCompositeBlockLinearSystemALLCASE_pfesT<R,MeshS,v_fes3,v_fesS>( i, j, offsetUh, offsetVh, PUh, PVh, initmat, initx, sym, tgv, b_largs_zz, stack, B, X, hm_A,mpirankandsize,B_from_varf);
      
  }
  else{
    cerr << "Other type for FESpace in construction, in construction ..." << endl;
    ffassert(0);
  }
}

template<class R>
AnyType OpMatrixtoBilinearFormVG<R>::Op::operator()(Stack stack) const
{
  assert(b && b->nargs);

  pvectgenericfes  * pUh= GetAny<pvectgenericfes *>((*b->euh)(stack));
  pvectgenericfes  * pVh= GetAny<pvectgenericfes *>((*b->evh)(stack));

  ffassert( *pUh && *pVh ); 
  // Update is necessary when we get "pvectgenericfes" to take account a new mesh for example.
  (*pUh)->update();
  (*pVh)->update();

  if( verbosity > 5){
    (*pUh)->printPointer();
    (*pVh)->printPointer();
  }
  int NpUh = (*pUh)->N; // number of fespace in pUh
  int NpVh = (*pVh)->N; // number of fespace in pVh

  KN<int> UhNbOfDf = (*pUh)->vectOfNbOfDF(); // A changer en long
  KN<int> VhNbOfDf = (*pVh)->vectOfNbOfDF();

  KN<int> UhNbItem = (*pUh)->vectOfNbitem();
  KN<int> VhNbItem = (*pVh)->vectOfNbitem();

  const KNM<list<C_F0>> & block_largs=b->block_largs; 
  
  // check if we have a square matrix
  bool A_is_square= (void*)pUh == (void*)pVh || ((*pUh)->totalNbOfDF()) == ( (*pVh)->totalNbOfDF()) ;


  // === simple check if A is symetrical === // 
  // voir avec les autres.
  bool A_is_maybe_sym = (void*)pUh == (void*)pVh; 

  // VF == true => VF type of Matrix
  bool VF=false;
  for(int i=0; i<NpUh;i++){
    for(int j=0; j<NpVh;j++){
      bool block_VF=isVF(b->block_largs(i,j));
      if( !VF ) VF=block_VF;
    }
  }
  //bool VF=isVF(b->largs);    //=== used to set the solver ??? block matrix ??? ===/

  // set parameteer of the matrix :: 
  Data_Sparse_Solver ds;
  ds.factorize=0;
  ds.initmat=true;
  int np_bem = OpCall_FormBilinear_np::n_name_param; // number of parameter with BEM
  int np = OpCall_FormBilinear_np::n_name_param - NB_NAME_PARM_HMAT;
  SetEnd_Data_Sparse_Solver<R>(stack,ds, b->nargs,np);


  // J'ai repris ce qu'il y avait. 
  // PAC(e)     :: Attention peut être pas compatible avec les matrices bloc.
  // A repenser :: surtout pour le parametre symetrique? on le met ce parametre à zéro pour l'instant.
  // set ds.sym = 0 

  ds.sym = 0;
  if(verbosity>3)
    cout << " === we consider the block matrix as a non symetric matrix === (to be change in the future)" << endl; 

  if (! A_is_square )
   {
     if(verbosity>3) cout << " -- the solver  is un set on rectangular matrix  " << endl;
    }

  // A quoi cela correspond?? Gestion du stack + autre
  WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH aout 2007

  Matrice_Creuse<R> & A( * GetAny<Matrice_Creuse<R>*>((*a)(stack)));
  if(init) A.init(); //
  if(verbosity>3){                                         
    cout << " A.N=" <<  A.N() << endl;
    cout << " A.M=" <<  A.M() << endl;
  }
  if( ! pUh || ! pVh) return SetAny<Matrice_Creuse<R>  *>(&A);  //

  // need to define the size of the entire matrix here ==> execution error
  A.resize( (*pVh)->totalNbOfDF(), (*pUh)->totalNbOfDF() ); 
  // test function (Vh) are the line
  // inconnu function (Uh) are the column

  // Assemble the variational form
  int maxJVh=NpVh;
  
  int offsetMatrixUh = 0;
  // loop over the block
  for( int i=0; i<NpUh; i++){
    int offsetMatrixVh = 0;
    if( ds.sym > 0 ){ maxJVh=(i+1); ffassert(maxJVh<NpVh);}
    for( int j=0; j<maxJVh; j++){
      if(verbosity>3) cout << "offsetMatrixUh= " << offsetMatrixUh << ", offsetMatrixVh= " << offsetMatrixVh << endl;
      
      // size of the block
      int N_block = UhNbOfDf[i];
      int M_block = VhNbOfDf[j];
    
      const list<C_F0> &b_largs = block_largs(i,j);
      if(verbosity>2) cout << "size_block =" << b_largs.size() << endl; 
      if( b_largs.size()> 0){
        list<C_F0> largs_FEM;
        list<C_F0> largs_BEM;
        largs_FEM.clear();
        largs_BEM.clear();
        separateFEMpartBemPart( b_largs, largs_FEM, largs_BEM );

        if(verbosity>2){
          cout << " FEM.size()=" << largs_FEM.size() << endl;
          cout << " BEM.size()=" << largs_BEM.size() << endl;
        }
        //largs_BEM.clear();
      // // init matrice creuse of the largs argument
    // Matrice_Creuse<R> *CCC = new Matrice_Creuse<R>() ;
    // CCC->init(); //
    //   CCC->resize(M_block,N_block); // test function (Vh) are the line and inconnu function (Uh) are the column
      // // cout << "block:  i=" << i << "j=" << j <<  " (N,M)=" << M_block << " " << N_block << endl;
      int method =1;
if(method == 1){       
        if( largs_BEM.size() >0 ){
          varfBemToCompositeBlockLinearSystem( i, j, (*pUh)->typeFE[i], (*pVh)->typeFE[j], (long) UhNbOfDf[i], (long) VhNbOfDf[j],
                                        (long) offsetMatrixUh, (long) offsetMatrixVh, (*pUh)->vect[i], (*pVh)->vect[j],
                                        largs_BEM, stack, b->nargs, A.pHM(),np_bem);
        }
}else if(method==2){

        if( largs_BEM.size() >0 ){
#ifndef FFLANG
#if defined(PARALLELE) && defined(WITH_bemtool) && defined(WITH_htool)
          const list<C_F0> & b_largs_zz = largs_BEM;
          
          int VFBEM = typeVFBEM(b_largs_zz,stack);
          if(VFBEM == 2){ cerr << " not implemented with BEM POTENTIAL" << endl; ffassert(0);}
          Data_Bem_Solver dsbem;
          dsbem.factorize=0;
          dsbem.initmat=true;
          SetEnd_Data_Bem_Solver<R>(stack, dsbem, b->nargs,OpCall_FormBilinear_np::n_name_param);  // LIST_NAME_PARM_HMAT

          HMatrixVirt<R> ** Hmat = new HMatrixVirt<R> *();
         
          //
          // a voir dans le futur si la difference entre bloc diagonal et bloc non diagonal a un sens. ::: 08/2022 :::  Morice
          //
          if( i==j ){

            bool samemesh = (void*) (*pUh)->vect[i]->getppTh() == (void*) (*pVh)->vect[j]->getppTh();  // same Fem2D::Mesh     +++ pot or kernel
            if (VFBEM==1)
              ffassert (samemesh);
            if(init)
              *Hmat =0;
            *Hmat =0;
            if( *Hmat)
              delete *Hmat;
            *Hmat =0;


            // block diagonal matrix
            if( (*pUh)->typeFE[i] == 4 && (*pVh)->typeFE[j] == 4 ){
              ffassert( i==j ); // If not a block diagonal not coded yet
              // MeshS --- MeshS
              // ==== FESpace 3d Surf: inconnue et test ===
              const FESpaceS * PUh = (FESpaceS *) (*pUh)->vect[i]->getpVh();
              const FESpaceS * PVh = (FESpaceS *) (*pVh)->vect[j]->getpVh();

              creationHMatrixtoBEMForm<R, MeshS, FESpaceS, FESpaceS>(PUh, PVh, VFBEM, 
                              b_largs_zz, stack, dsbem, Hmat);

            }
            else if( (*pUh)->typeFE[i] == 5 && (*pVh)->typeFE[j] == 5 ){
              ffassert( i==j ); // If not a block diagonal not coded yet
              // MeshL --- MeshL
              // ==== FESpace 3d Curve: inconnue et test ===
              const FESpaceL * PUh = (FESpaceL *) (*pUh)->vect[i]->getpVh();
              const FESpaceL * PVh = (FESpaceL *) (*pVh)->vect[j]->getpVh();

              creationHMatrixtoBEMForm<R, MeshL, FESpaceL, FESpaceL> ( PUh, PVh, VFBEM, 
                              b_largs_zz, stack, dsbem, Hmat );
            }
            else{
              cerr << " BEM bilinear form " << endl;
              cerr << " Block ("<< i <<" ,"<< j << ")" << endl;
              cerr << " =: Pas prise en compte des FESpace inconnue de type := "<< typeFEtoString( (*pUh)->typeFE[i] ) << endl;
              cerr << " =:                 avec des FESpace test de type    := "<< typeFEtoString( (*pVh)->typeFE[j] ) << endl;
              ffassert(0);
            }

            // creation de la matrice dense 
            KNM<R>* M= HMatrixVirtToDense< KNM<R>, R >(Hmat);

            HashMatrix<int,R> *phm= new HashMatrix<int,R>(*M);
            MatriceCreuse<R> *pmc(phm);

            Matrice_Creuse<R> BBB;
            BBB.A=0;
            BBB.A.master(pmc);

            A.pHM()->Add( BBB.pHM(), R(1), false, offsetMatrixVh, offsetMatrixUh ); // test function (Vh) are the line and inconnu function (Uh) are the column

            M->destroy();
            delete M;
            BBB.destroy();
          }
          else{

            bool samemesh = (void*) (*pUh)->vect[i]->getppTh() == (void*) (*pVh)->vect[j]->getppTh();  // same Fem2D::Mesh     +++ pot or kernel
          
            if(init)
              *Hmat =0;
            if( *Hmat)
              delete *Hmat;
            *Hmat =0;
            
            
            // block non diagonal matrix        
            if( (*pUh)->typeFE[i] == 5 && (*pVh)->typeFE[j] == 2 ){
              // case Uh[i] == MeshL et Vh[j] = Mesh2  // Est ce que cela a un sens?
              
              cout << " === creation de la matrice BEM pour un bloc non diagonaux === " << endl;
              //ffassert(0);
              const FESpaceL * PUh = (FESpaceL *) (*pUh)->vect[i]->getpVh();
              creationHMatrixtoBEMForm<R, MeshL, FESpaceL, FESpaceL> ( PUh, PUh, VFBEM, 
                          b_largs_zz, stack, dsbem, Hmat );

            }
             
            else{
              cerr << " BEM bilinear form " << endl;
              cerr << " Block ("<< i <<" ,"<< j << ")" << endl;
              cerr << " =: Pas prise en compte des FESpace inconnue de type := "<< typeFEtoString( (*pUh)->typeFE[i] ) << endl;
              cerr << " =:                 avec des FESpace test de type    := "<< typeFEtoString( (*pVh)->typeFE[j] ) << endl;
              ffassert(0);
            }
            
            // creation de la matrice dense 
            
            KNM<R>* M = HMatrixVirtToDense< KNM<R>, R >(Hmat);
            
            HashMatrix<int,R> *phm= new HashMatrix<int,R>(*M);
            MatriceCreuse<R> *pmc(phm);
            
            Matrice_Creuse<R> *BBB=new Matrice_Creuse<R>();
            BBB->A=0;
            BBB->A.master(pmc);
            //BBB->resize(356,356);

            // BEM matrix is constructed with different FESpace
            ffassert( (*pUh)->vect[i]->getpVh() != (*pVh)->vect[j]->getpVh() ) ;
            
            
            if( (*pUh)->typeFE[i] == 5 && (*pVh)->typeFE[j] == 2 ){
              // case Uh[i] == MeshL et Vh[j] = Mesh2 
              const FESpaceL * PUh = (FESpaceL *) (*pUh)->vect[i]->getpVh();
              const FESpace * PVh = (FESpace *) (*pVh)->vect[j]->getpVh();
              // construction of the matrix of interpolation
              
              //Matrice_Creuse<double> *  MI_BBB=new Matrice_Creuse<double>(); // = buildMatrixInterpolationForCompositeFESpace<double,FESpaceL,FESpace>( PUh, PVh  );
              Matrice_Creuse<double> *  MI_BBB = buildMatrixInterpolationForCompositeFESpace<double,FESpaceL,FESpace>( PUh, PVh  );

              //MI_BBB->resize(356,2922);
              // multiplication matrix*matrix
            
              MatriceMorse<double> *mA= MI_BBB->pHM();
              MatriceMorse<R> *mB= BBB->pHM();            

              //cout << "A=MI " << MI_BBB->N() << " " << MI_BBB->M() << endl;
              //cout << "B=BBB " << BBB->N() << " " << BBB->M() << endl;

              ffassert( MI_BBB->M() >= BBB->M() ); 
              
              ffassert( MI_BBB->N() == BBB->M() );
              
              MatriceMorse<R> *mAB=new MatriceMorse<R>(MI_BBB->M(), BBB->M(),0,0);
              AddMul<int,double,R,R>(*mAB,*mA,*mB,true,false); // BBB=MI_BBB'*BBB;
              
              A.pHM()->Add( mAB, R(1), false, offsetMatrixVh, offsetMatrixUh ); // test function (Vh) are the line and inconnu function (Uh) are the column
              
              delete mAB;
              MI_BBB->destroy();
              delete MI_BBB;
              
            }
            else{
              cerr << "==== to do ==== " << endl;
              ffassert(0);
            }
            M->destroy();
            delete M;
            BBB->destroy();
            delete BBB;
            
          }
          if( *Hmat)
            delete *Hmat;
          delete Hmat;
#endif
#endif
        }
}       
        if( largs_FEM.size() >0 ){
          // computation of the matrix
          const list<C_F0> & b_largs_zz = largs_FEM;
          
          varfToCompositeBlockLinearSystemALLCASE_pfes<R>( i, j, (*pUh)->typeFE[i], (*pVh)->typeFE[j], 
                                                        offsetMatrixUh, offsetMatrixVh, (*pUh)->vect[i], (*pVh)->vect[j],
                                                        true, false, ds.sym, ds.tgv, ds.commworld,
                                                        b_largs_zz, stack, 
                                                        0, 0, A.pHM());
        }

        if(mpirank ==0 && verbosity >3) cout << "Add block (" << i << "," << j <<")=> size nnz=" << A.pHM()->nnz << endl;
      //delete CCC;
      } // end loop bb_largs_ii
    offsetMatrixVh += VhNbOfDf[j];
    } // end loop j
    offsetMatrixUh += UhNbOfDf[i];
  } // end loop i
  
  //pUh = nullptr;
  //pVh = nullptr;  
  //if( *pUh ){ (*pUh)->destroy();}

  A.pHM()->half = ds.sym;
  if (A_is_square)
    SetSolver(stack,VF,*A.A,ds);

  return SetAny<Matrice_Creuse<R>  *>(&A);
}


// Mesh - Mesh
template void varfToCompositeBlockLinearSystem< double, Mesh, FESpace, FESpace>
                              (bool initmat, bool initx, const FESpace * PUh, const FESpace * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<double> *B, KN_<double> *X, MatriceCreuse<double> &A,int *mpirankandsize, bool B_from_varf);

template void varfToCompositeBlockLinearSystem< Complex, Mesh, FESpace, FESpace>
                              (bool initmat, bool initx, const FESpace * PUh, const FESpace * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<Complex> *B, KN_<Complex> *X, MatriceCreuse<Complex> &A,int *mpirankandsize, bool B_from_varf);
// MeshL - MeshL
template void varfToCompositeBlockLinearSystem< double, MeshL, FESpaceL, FESpaceL>
                              (bool initmat, bool initx, const FESpaceL * PUh, const FESpaceL * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<double> *B, KN_<double> *X, MatriceCreuse<double> &A,int *mpirankandsize, bool B_from_varf);

template void varfToCompositeBlockLinearSystem< Complex, MeshL, FESpaceL, FESpaceL>
                              (bool initmat, bool initx, const FESpaceL * PUh, const FESpaceL * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<Complex> *B, KN_<Complex> *X, MatriceCreuse<Complex> &A,int *mpirankandsize, bool B_from_varf);

// Mesh - MeshL
template void varfToCompositeBlockLinearSystem< double, MeshL, FESpace, FESpaceL>
                              (bool initmat, bool initx, const FESpace * PUh, const FESpaceL * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<double> *B, KN_<double> *X, MatriceCreuse<double> &A,int *mpirankandsize, bool B_from_varf);

template void varfToCompositeBlockLinearSystem< Complex, MeshL, FESpace, FESpaceL>
                              (bool initmat, bool initx, const FESpace * PUh, const FESpaceL * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<Complex> *B, KN_<Complex> *X, MatriceCreuse<Complex> &A,int *mpirankandsize, bool B_from_varf);

// MeshL - Mesh
template void varfToCompositeBlockLinearSystem< double, MeshL, FESpaceL, FESpace>
                              (bool initmat, bool initx, const FESpaceL * PUh, const FESpace * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<double> *B, KN_<double> *X, MatriceCreuse<double> &A,int *mpirankandsize, bool B_from_varf);

template void varfToCompositeBlockLinearSystem< Complex, MeshL, FESpaceL, FESpace>
                              (bool initmat, bool initx, const FESpaceL * PUh, const FESpace * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<Complex> *B, KN_<Complex> *X, MatriceCreuse<Complex> &A,int *mpirankandsize, bool B_from_varf);

// FEM block

template void varfToCompositeBlockLinearSystemALLCASE_pfes( const int& i, const int &j, 
                const int &typeUh, const int &typeVh, 
                const long &offsetUh, const long &offsetVh,
                const generic_v_fes * pfesUh, const generic_v_fes * pfesVh,
                bool initmat, bool initx, const int &sym, const double &tgv, 
                const pcommworld &comm, const list<C_F0> & b_largs_zz, Stack stack,
                KN_<double> *B, KN_<double> *X, HashMatrix<int,double> *hm_A, bool B_from_varf=false);

template void varfToCompositeBlockLinearSystemALLCASE_pfes( const int& i, const int &j, 
                const int &typeUh, const int &typeVh, 
                const long &offsetUh, const long &offsetVh,
                const generic_v_fes * pfesUh, const generic_v_fes * pfesVh,
                bool initmat, bool initx, const int &sym, const double &tgv, 
                const pcommworld &comm, const list<C_F0> & b_largs_zz, Stack stack,
                KN_<Complex> *B, KN_<Complex> *X, HashMatrix<int,Complex> *hm_A, bool B_from_varf=false);


// BEM block

template<> void varfBemToCompositeBlockLinearSystem<double>(const int& i, const int &j, 
                const int &typeUh, const  int &typeVh,
                const long &sizeUh, const long &sizeVh,
                const long &offsetUh, const long &offsetVh,
                const generic_v_fes * LLUh, const generic_v_fes *LLVh,
                const list<C_F0> & b_largs_zz, Stack stack, Expression const * nargs,
                HashMatrix<int,double> *hm_A,const int &n_name_param){
                  ffassert(0);
                }

template void varfBemToCompositeBlockLinearSystem(const int& i, const int &j, 
                const int &typeUh, const int &typeVh,
                const long &sizeUh, const long &sizeVh,
                const long &offsetUh, const long &offsetVh,
                const generic_v_fes * LLUh, const generic_v_fes *LLVh,
                const list<C_F0> & b_largs_zz, Stack stack, Expression const *nargs,
                HashMatrix<int,Complex> *hm_A,const int &n_name_param);

template AnyType OpMatrixtoBilinearFormVG<double>::Op::operator()(Stack stack) const;

template AnyType OpMatrixtoBilinearFormVG<Complex>::Op::operator()(Stack stack) const;
