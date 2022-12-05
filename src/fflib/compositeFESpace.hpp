//  composite FESpace
#ifndef COMPOSITE_FESPACE_HPP_
#define COMPOSITE_FESPACE_HPP_

template<class R>  //  to make   A=linearform(x)
struct OpMatrixtoBilinearFormVG
  : public OneOperator
{
  typedef typename Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes>::const_iterator const_iterator;
  int init;
  
  class Op : public E_F0mps {
    public:
      Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes> *b;
      Expression a;
      int init;
      //AnyType operator()(Stack s)  const;
      
      Op(Expression aa,Expression  bb,int initt)
        : b(new Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes>(* dynamic_cast<const Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes> *>(bb))),a(aa),init(initt)
    { 
      assert(b && b->nargs);
      int NN = (int) b->euh->componentNbitem().size();
      int MM = (int) b->evh->componentNbitem().size();

      bool total_iscmplx=false;
      // loop over block
      for(int i=0; i<NN; i++){
        for(int j=0; j<MM; j++){
          // FieldOfForm : optimize the terms (flags -O3) of the variational form and verifies the type of the variational form
          bool iscmplx=FieldOfForm(b->block_largs(i,j),IsComplexType<R>::value)  ;
          // cout<< "FieldOfForm:iscmplx " << iscmplx << " " << IsComplexType<R>::value << " " << ((iscmplx) == IsComplexType<R>::value) << endl;
          ffassert( (iscmplx) == IsComplexType<R>::value);
          if( !total_iscmplx ) total_iscmplx=iscmplx;
        }
      }
    }
    operator aType () const { return atype<Matrice_Creuse<R>  *>();}

    AnyType operator()(Stack s)  const;
  };

  E_F0 * code(const basicAC_F0 & args) const
  { return  new Op(to<Matrice_Creuse<R>*>(args[0]),args[1],init); }
  OpMatrixtoBilinearFormVG(int initt=0) :
    OneOperator(atype<Matrice_Creuse<R>*>(),atype<Matrice_Creuse<R>*>(),atype<const Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes>*>()),
    init(initt){};

};


string typeFEtoString(int typeFE);



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
                                        const KN<int> &indexBlockUh, const KN<int> &indexBlockVh );



KNM< list<C_F0> > computeBlockLargs( const list<C_F0> & largs, const int &NpUh, const int &NpVh, const KN<int> &indexBlockUh, const KN<int> &indexBlockVh );
  


// Info necessaire :: " block_largs, localIndexInTheBlockUh, localIndexInTheBlockVh, NpUh, NpVh  
void changeComponentFormCompositeFESpace( const KN<int> &localIndexInTheBlockUh, const KN<int> &localIndexInTheBlockVh, 
        KNM< list<C_F0> > & block_largs );
  
void reverseChangeComponentFormCompositeFESpace(const KN<int>  &beginBlockUh, const KN<int> &beginBlockVh, 
          KNM< list<C_F0> > & block_largs);

template<class FESpaceT1,class FESpaceT2>
MatriceMorse<R> * buildInterpolationMatrixT(const FESpaceT1 & Uh,const FESpaceT2 & Vh,void *data);

template< >
MatriceMorse<R> * buildInterpolationMatrixT<FESpaceL,FESpace>(const FESpaceL & Uh,const FESpace & Vh,void *data);


template< class R, class FESpaceT1, class FESpaceT2 >
Matrice_Creuse<R> *  buildMatrixInterpolationForCompositeFESpace(const FESpaceT1 * Uh ,const FESpaceT2 * Vh);

// print information of the varf
void listOfComponentBilinearForm(const list<C_F0> & largs);

/**
       *  @brief  determine if we have BEM bilinear operator in a subblock and the type
       *  @param  largs list of argument of the Bilinear Form
       */

/*
  This function give the good result if only if the FESpace inconnu and  FESpace test are scalar FESpace
  due to this check : 
          if( finc.first==0 && ftest.first==0)      // only first component for finc and ftest
*/

int haveBemSubMatrixBlock(const list<C_F0> & largs, int Uh_NbItem, int Vh_NbItem);

/**
       *  @brief  largs separate in two part :  BEM (H-matrix) and FEM 
       *  @param  largs list of argument of the Bilinear Form
       *  @param  largs_FEM list of argument for the FEM part
       *  @param  largs_BEM list of argument for the BEM part (included sometimes mass matrix ) that be compressed in H-matrix
       */

/*
  This function must be call if we have BEM Bilinear operator in a block.
  
*/
void separateFEMpartBemPart(const list<C_F0> & largs, list<C_F0> &largs_FEM, list<C_F0> &largs_BEM );


/**
       *  @brief  Function to delete element of newlargs = list<C_F0> to avoid memory leak.
       *  @param  largs list of argument of the composite FESpace gene
       */

/* remark: only FormBilinear  and BemFormBilinear is deleted */
void deleteNewLargs(list<C_F0> &newlargs);

/**
       *  @brief  Function to create a matrix of a composite FESpace (FreeFem Matrix)
       *  @param  largs list of argument of the composite FESpace gene
       */

template<class R>
void varfToCompositeBlockMatrix( const int &i_Uh, const int &j_Vh, pvectgenericfes  *pUh, pvectgenericfes  *pVh, 
                         //vect_generic_v_fes **pUh, vect_generic_v_fes **pVh, 
                         const int &sym, const double &tgv, const list<C_F0> & b_largs_zz, Stack stack, 
                         Matrice_Creuse<R> &BBB); 
/*
template<class R> 
void computationGlobalMatrix( const KNM<list<C_F0>> &block_largs, Expression *nargs, 
                              pvectgenericfes  * pUh, pvectgenericfes  * pVh ){
}
*/

template<class R,class MMesh, class FESpace1, class FESpace2>
void varfToCompositeBlockLinearSystem(bool initmat, bool initx, const FESpace1 * PUh, const FESpace2 * PVh, 
                              const int &sym, const double &tgv, const list<C_F0> & largs, Stack stack, 
                              KN_<R> *B, KN_<R> *X, MatriceCreuse<R> &A,int *mpirankandsize, bool B_from_varf);


template< class R>
void varfBemToCompositeBlockLinearSystem(const int& i, const int &j, 
                const int &typeUh, const  int &typeVh,
                const long &sizeUh, const long &sizeVh,
                const long &offsetUh, const long &offsetVh,
                const generic_v_fes *LLUh, const generic_v_fes * LLVh,
                const list<C_F0> & b_largs_zz, Stack stack, Expression const * nargs,
                HashMatrix<int,R> *hm_A, const int &n_name_param);

template<class K> 
void varfToCompositeBlockLinearSystemALLCASE_pfes( const int& i, const int &j, 
                const int &typeUh, const int &typeVh,
                const long &offsetUh, const long &offsetVh,
                const generic_v_fes *pfesUh, const generic_v_fes *pfesVh,
                bool initmat, bool initx, const int &sym, const double &tgv, 
                const list<C_F0> & b_largs_zz, Stack stack, 
                KN_<K> *B, KN_<K> *X, HashMatrix<int,K> *hm_A, bool B_from_varf=false);


template<class R> 
void varfToCompositeBlockLinearSystemALLCASE( const int& i, const int &j, 
                const vector< int> &typeUh, const vector< int> &typeVh,
                const vector< long> &sizeUh, const vector< long> &sizeVh,
                const vector< long> &offsetUh, const vector< long> &offsetVh,
                const vector< void *> &LLUh, const vector< void *> &LLVh,
                bool initmat, bool initx, const int &sym, const double &tgv, 
                const list<C_F0> & b_largs_zz, Stack stack, 
                KN_<R> *B, KN_<R> *X, HashMatrix<int,R> *hm_A,int *mpirankandsize, bool B_from_varf=false);

#endif