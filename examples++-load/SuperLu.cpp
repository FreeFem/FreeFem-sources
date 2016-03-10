//   for automatic  compilation with ff-c++
//ff-c++-LIBRARY-dep:   superlu4 blas 
//ff-c++-cpp-dep: 
//  for Super4.0 library   
#include "ff++.hpp"
#include "slu_ddefs.h"
#include "superlu_enum_consts.h"
#define GlobalLU_t GlobalLU_txxxx
#define countnz countnzxxxx
#define fixupL fixupLxxxx
#define print_lu_col print_lu_colxxxx
#define check_tempv check_tempvxxxx
#define PrintPerf PrintPerfxxxx
#define ilu_countnz  ilu_countnzxxxx
#include "slu_zdefs.h"

#undef GlobalLU_t
#undef countnz
#undef fixupL
#undef print_lu_col
#undef check_tempv
#undef PrintPerf
#undef ilu_countnz


template <class R> struct SuperLUDriver
{
    
};

template <> struct SuperLUDriver<double>
{
    /* Driver routines */
    static  Dtype_t R_SLU_T() { return SLU_D;} 
    static void
    gssv(superlu_options_t * p1, SuperMatrix * p2, int * p3, int * p4, SuperMatrix * p5,
	  SuperMatrix * p6, SuperMatrix * p7 , SuperLUStat_t * p8, int * p9)
    { dgssv( p1,p2,p3,p4,p5,p6,p7,p8,p9); }
    
    
    static void
	gssvx(superlu_options_t * p1, SuperMatrix * p2, int * p3, int * p4, int * p5,
	       char * p6, double * p7, double * p8, SuperMatrix * p9, SuperMatrix * p10,
	       void * p11, int p12, SuperMatrix * p13, SuperMatrix * p14,
	       double * p15, double * p16, double * p17, double * p18,
	       mem_usage_t * p19, SuperLUStat_t * p20, int * p21)
    { dgssvx( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,  p11,p12,p13,p14,p15,p16,p17,p18,p19,p20, p21); }
    
    
    
    /* Supernodal LU factor related */
    static void
	Create_CompCol_Matrix(SuperMatrix * p1, int p2 , int p3, int p4, double * p5,
			       int * p6, int * p7, Stype_t p8, Dtype_t p9 , Mtype_t p10)
    {
	    dCreate_CompCol_Matrix( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10);
    }
    
    
    static void
	Create_CompRow_Matrix(SuperMatrix * p1, int p2, int p3, int p4, double * p5,
			       int * p6, int * p7, Stype_t p8, Dtype_t p9, Mtype_t p10)
    {
	dCreate_CompRow_Matrix( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10);
    }
    
    
    static void
	Create_Dense_Matrix(SuperMatrix * p1, int p2, int p3, double * p4, int p5,
			     Stype_t p6, Dtype_t p7, Mtype_t p8)
    {
	dCreate_Dense_Matrix( p1,p2,p3,p4,p5,p6,p7,p8);
    }
    
    
    static void
	Create_SuperNode_Matrix(SuperMatrix * p1, int p2, int p3, int p4, double * p5, 
				 int * p6, int * p7, int * p8, int * p9, int * p10,
				 Stype_t p11, Dtype_t p12, Mtype_t p13)
    {
	    dCreate_SuperNode_Matrix( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,  p11,p12,p13);
    }
    
    static void 
    CompRow_to_CompCol(int p1, int p2, int p3, 
		       double *p4, int *p5, int *p6,
		       double **p7, int **p8, int **p9)
  {
    dCompRow_to_CompCol( p1, p2, p3, p4, p5, p6, p7, p8, p9);
  }

    
};



template <> struct SuperLUDriver<Complex>
{
    /* Driver routines */
  static  Dtype_t R_SLU_T() { return SLU_Z;} 
  static doublecomplex *dc(Complex *p)  { return (doublecomplex *) (void *) p;}
  static doublecomplex **dc(Complex **p)  { return (doublecomplex **) (void *) p;}
  
    static void
    gssv(superlu_options_t * p1, SuperMatrix * p2, int * p3, int * p4, SuperMatrix * p5,
	 SuperMatrix * p6, SuperMatrix * p7 , SuperLUStat_t * p8, int * p9)
    { zgssv( p1,p2,p3,p4,p5,p6,p7,p8,p9); }
    
    
    static void
    gssvx(superlu_options_t * p1, SuperMatrix * p2, int * p3, int * p4, int * p5,
	  char * p6, double * p7, double * p8, SuperMatrix * p9, SuperMatrix * p10,
	  void * p11, int p12, SuperMatrix * p13, SuperMatrix * p14,
	  double * p15, double * p16, double * p17, double * p18,
	  mem_usage_t * p19, SuperLUStat_t * p20, int * p21)
    { zgssvx( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,  p11,p12,p13,p14,p15,p16,p17,p18,p19,p20, p21); }
    
    
    
    /* Supernodal LU factor related */
    static void
    Create_CompCol_Matrix(SuperMatrix * p1, int p2 , int p3, int p4, Complex * p5,
			  int * p6, int * p7, Stype_t p8, Dtype_t p9 , Mtype_t p10)
    {
	zCreate_CompCol_Matrix( p1,p2,p3,p4,dc(p5),p6,p7,p8,p9,p10);
    }
    
    
    static void
    Create_CompRow_Matrix(SuperMatrix * p1, int p2, int p3, int p4, Complex * p5,
			  int * p6, int * p7, Stype_t p8, Dtype_t p9, Mtype_t p10)
    {
	zCreate_CompRow_Matrix( p1,p2,p3,p4,dc(p5),p6,p7,p8,p9,p10);
    }
    
    
    static void
    Create_Dense_Matrix(SuperMatrix * p1, int p2, int p3, Complex * p4, int p5,
			Stype_t p6, Dtype_t p7, Mtype_t p8)
    {
	zCreate_Dense_Matrix( p1,p2,p3,dc(p4),p5,p6,p7,p8);
    }
    
    
    static void
    Create_SuperNode_Matrix(SuperMatrix * p1, int p2, int p3, int p4, Complex * p5, 
			    int * p6, int * p7, int * p8, int * p9, int * p10,
			    Stype_t p11, Dtype_t p12, Mtype_t p13)
    {
	zCreate_SuperNode_Matrix( p1,p2,p3,p4,dc(p5),p6,p7,p8,p9,p10,  p11,p12,p13);
    }
    
  static void 
  CompRow_to_CompCol(int p1, int p2, int p3, Complex *p4, int *p5, 
		     int *p6, Complex **p7, int **p8, int **p9)
  {
    zCompRow_to_CompCol( p1, p2, p3, dc(p4), p5, p6, dc(p7), p8, p9);
  }
    
};

// read options for superlu in freefem++
/*
#ifdef __cpluscplus
int s_(char *ff, ...)
{
  int i = 0;
  while( *(++i+&str) != 0 )
    if( strcmp(str, (char*)*(&str+i)) == 0)
      return i;
  return 0;
}

#else
*/
int s_(char* str, const char* cmp[])
{
  int i = 0;
  while( cmp[i] != 0){
    if( strcmp(str, cmp[i]) == 0){
      //cout << *str << " return" << i << endl;
      return i+1 ;
    }
    i++;
  }
  //cout << *str << " return 0" << endl;
  return 0;
}
//#endif
/*
  static const yes_no_t  enumyes_no_t[2] = {NO, YES};
  static const fact_t  enumfact_t[4] = {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED};
  static const colperm_t  enumcolperm_t[5] = {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD, MY_PERMC};
  static const trans_t  enumtrans_t[3] = {NOTRANS, TRANS, CONJ};
  static const  IterRefine_t enumIterRefine_t[4] = {NOREFINE, SINGLE, DOUBLE, EXTRA};  
  
  static const char*  compyes_no_t[] = {"NO", "YES",0};
  static const char* compfact_t[] = {"DOFACT", "SamePattern", "SamePattern_SameRowPerm", "FACTORED",0};
  static const char* compcolperm_t[] = {"NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD", "MY_PERMC",0};
  static const char* comptrans_t[] = {"NOTRANS", "TRANS", "CONJ",0};
  static const char* compIterRefine_t[] = {"NOREFINE", "SINGLE", "DOUBLE", "EXTRA",0};
  
  static const char* comp[] = {"Fact", "Equil","ColPerm",
  "DiagPivotThresh","Trans","IterRefine",
  "SymmetricMode","PivotGrowth","ConditionNumber",
  "PrintStat",0};
*/

void read_options_freefem(string string_option, superlu_options_t *options){
  static const yes_no_t  enumyes_no_t[2] = {NO, YES};
  static const fact_t  enumfact_t[4] = {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED};
  static const colperm_t  enumcolperm_t[5] = {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD, MY_PERMC};
  static const trans_t  enumtrans_t[3] = {NOTRANS, TRANS, CONJ};
  static const  IterRefine_t enumIterRefine_t[4] = {NOREFINE, SLU_SINGLE, SLU_DOUBLE, SLU_EXTRA};  

  static const char*  compyes_no_t[] = {"NO", "YES",0};
  static const char* compfact_t[] = {"DOFACT", "SamePattern", "SamePattern_SameRowPerm", "FACTORED",0};
  static const char* compcolperm_t[] = {"NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD", "MY_PERMC",0};
  static const char* comptrans_t[] = {"NOTRANS", "TRANS", "CONJ",0};
  static const char* compIterRefine_t[] = {"NOREFINE", "SINGLE", "DOUBLE", "EXTRA",0};
  
  static const char* comp[] = {"Fact", "Equil","ColPerm",
			       "DiagPivotThresh","Trans","IterRefine",
			       "SymmetricMode","PivotGrowth","ConditionNumber",
			       "PrintStat",0};


   /* Set the default values for options argument:
	options.Fact = DOFACT;
	options.Equil = YES;
	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
	options.Trans = NOTRANS;
	options.IterRefine = NOREFINE;
	options.SymmetricMode = NO;
	options.PivotGrowth = NO;
	options.ConditionNumber = NO;
	options.PrintStat = YES;
    */
  //cout << "string_option" <<  *string_option << endl;
    KN<char> kdata(string_option.size()+1);
    
    char * data=kdata;
  strcpy( data, string_option.c_str()); 
  cout << "data=" << data << endl;
  char * tictac;
  tictac = strtok(data," =,\t\n");
  cout << "tictac=" << data << endl;
// #ifdef __cplusplus
//   while(tictac != NULL){
//     int id_option = s_(tictac, "Fact", "Equil","ColPerm",
// 				"DiagPivotThresh","Trans","IterRefine",
// 				"SymmetricMode","PivotGrowth","ConditionNumber",
// 				"PrintStat",0);
//     tictac = strtok(NULL," ,\t\n");
//     int val_options;
//     switch (id_option)
//       { 
//       case 1 : // Fact
// 	val_options= s_(tictac, "DOFACT", "SamePattern", "SamePattern_SameRowPerm", "FACTORED",0);
// 	if( val_options == 0){
// 	  printf("value given for SuperLU for options %s is not correct\n","Fact");
// 	  exit(1);
// 	}
// 	options->Fact= enumfact_t[val_options-1];
// 	break;
//       case 2:  // Equil
// 	val_options= s_(tictac, "NO", "YES", 0);
// 	if( val_options == 0){
// 	  printf("value given for SuperLU for options %s is not correct\n","Equil");
// 	  exit(1);
// 	}
// 	options->Equil= enumyes_no_t[val_options-1];
// 	break;
//       case 3:  // ColPerm
// 	val_options= s_(tictac,"NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD", "MY_PERMC", 0);
// 	if( val_options == 0){
// 	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
// 	  exit(1);
// 	}
// 	options->ColPerm= enumcolperm_t[val_options-1];
//       case 4:  // DiagPivotThresh
// 	options->DiagPivotThresh= strtod(tictac,&tictac);
// 	break;
//       case 5:  // Trans
// 	val_options= s_(tictac, "NOTRANS", "TRANS", "CONJ",0);
// 	if( val_options == 0){
// 	  printf("value given for SuperLU for options %s is not correct\n","Trans");
// 	  exit(1);
// 	}
// 	options->Trans= enumtrans_t[val_options-1];
// 	break;
//       case 6:  // IterRefine
// 	val_options= s_(tictac, "NOREFINE", "SINGLE", "DOUBLE", "EXTRA",0);
// 	if( val_options == 0){
// 	  printf("value given for SuperLU for options %s is not correct\n","IterRefine");
// 	  exit(1);
// 	}
// 	options->IterRefine= enumIterRefine_t[val_options-1];
// 	break;
//       case 7:  // SymmetricMode
// 	val_options= s_(tictac, "NO","YES",0);
// 	if( val_options == 0){
// 	  printf("value given for SuperLU for options %s is not correct\n","SymmetricMode");
// 	  exit(1);
// 	}
// 	options->SymmetricMode= enumyes_no_t[val_options-1];
// 	break;
//       case 8:  // PivotGrowth
// 	val_options= s_(tictac, "NO","YES",0);
// 	if( val_options == 0){
// 	  printf("value given for SuperLU for options %s is not correct\n","PivotGrowth");
// 	  exit(1);
// 	}
// 	options->PivotGrowth= enumyes_no_t[val_options-1];
// 	break;
//       case 9:  // ConditionNumber
// 	val_options= s_(tictac, "NO","YES",0);
// 	if( val_options == 0){
// 	  printf("value given for SuperLU for options %s is not correct\n","ConditionNumber");
// 	  exit(1);
// 	}
// 	options->ConditionNumber = enumyes_no_t[val_options-1];
// 	break;
//       case 10: // PrintStat
// 	val_options= s_(tictac, "NO","YES",0);
// 	if( val_options == 0){
// 	  printf("value given for SuperLU for options %s is not correct\n","PrintStat");
// 	  exit(1);
// 	}
// 	options->PrintStat = enumyes_no_t[val_options-1];
// 	break;
//       case 0: // Equivalent of case default
// 	printf("A false parameter for  SuperLU is given %s \n",tictac);
// 	exit(1);
//       }  
//     tictac = strtok(NULL," ,\t\n");
//   }
// #else
  while(tictac != NULL){
    //char* comp[] = {"Fact", "Equil","ColPerm",
    //"DiagPivotThresh","Trans","IterRefine",
    //"SymmetricMode","PivotGrowth","ConditionNumber",
    //"PrintStat",0 };
    int id_option = s_(tictac, comp);
    tictac = strtok(NULL," =,\t\n");
    int val_options;

    switch (id_option)
      { 
      case 1 : // Fact
	//char* comp1[] = {"DOFACT", "SamePattern", "SamePattern_SameRowPerm", "FACTORED",0};
	val_options= s_(tictac,compfact_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","Fact");
	  exit(1);
	}
	options->Fact = enumfact_t[val_options-1];
	break;
      case 2:  // Equil
	//char* comp2[] = {"NO", "YES", 0};
	val_options= s_(tictac,compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","Equil");
	  exit(1);
	}
	options->Equil = enumyes_no_t[val_options-1];
	break;
      case 3:  // ColPerm
	//char* comp3[] = {"NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD", "MY_PERMC", 0};
	val_options= s_(tictac,compcolperm_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
	  exit(1);
	}
	options->ColPerm = enumcolperm_t[val_options-1];
	break;
      case 4:  // DiagPivotThresh
	options->DiagPivotThresh= strtod(tictac,&tictac);
	break;
      case 5:  // Trans
	//char* comp5[] = {"NOTRANS", "TRANS", "CONJ", 0};
	val_options= s_(tictac, comptrans_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","Trans");
	  exit(1);
	}
	options->Trans = enumtrans_t[val_options-1];
	break;
      case 6:  // IterRefine
	//char* comp6[] = {"NOREFINE", "SINGLE", "DOUBLE", "EXTRA", 0};
	val_options= s_(tictac, compIterRefine_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","IterRefine");
	  exit(1);
	}
	options->IterRefine = enumIterRefine_t[val_options-1];
	break;
      case 7:  // SymmetricMode
	//char* comp7[] = {"NO","YES", 0};
	val_options= s_(tictac, compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","SymmetricMode");
	  exit(1);
	}
	options->SymmetricMode= enumyes_no_t[val_options-1]; 
	break;
      case 8:  // PivotGrowth
	//char* comp8[] = {"NO","YES", 0};
	val_options= s_(tictac,compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","PivotGrowth");
	  exit(1);
	}
	options->PivotGrowth = enumyes_no_t[val_options-1];
	break;
      case 9:  // ConditionNumber
	//char* comp9[] = {"NO","YES", 0};
	val_options= s_(tictac, compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ConditionNumber");
	  exit(1);
	}
	options->ConditionNumber = enumyes_no_t[val_options-1];
	break;
      case 10: // PrintStat
	val_options= s_(tictac, compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","PrintStat");
	  exit(1);
	}
	options->PrintStat = enumyes_no_t[val_options-1];
	break;
      case 0: // Equivalent of case default
	break;
      }  
    tictac = strtok(NULL," =,\t\n");
  }
  //#endif
}

  




template<class R>
class SolveSuperLU :   public MatriceMorse<R>::VirtualSolver, public SuperLUDriver<R>   {
  double eps;
  mutable double  epsr;
  double tgv;
  double tol_pivot_sym,tol_pivot; //Add 31 oct 2005
  
  
  mutable char           equed[1];
  yes_no_t       equil;
  mutable SuperMatrix    A, L, U;
  NCformat       *Astore;
  NCformat       *Ustore;
  SCformat       *Lstore;
  R              *a;
  int            *asub, *xa;
  KN<int>             perm_c; /* column permutation vector */
  KN<int>             perm_r; /* row permutations from partial pivoting */
  string string_option;
  //string *file_option;
  //string *file_perm_r;
  //string *file_perm_c;
  
  KN<int>            etree;
  R         *rhsb, *rhsx, *xact;
  double         *RR, *CC;
  int m, n, nnz;
  
  R         *arow;
  int       *asubrow, *xarow;
  
   
   mutable superlu_options_t options;
   mutable mem_usage_t    mem_usage;
   
public:
  SolveSuperLU(const MatriceMorse<R> &AA,int strategy,double ttgv, double epsilon,
	       double pivot,double pivot_sym, string & param_char, KN<long> pperm_r, 
	       KN<long> pperm_c ) : 
    eps(epsilon),epsr(0),
    tgv(ttgv),
    etree(0),string_option(param_char),perm_r(pperm_r), perm_c(pperm_c),
    RR(0), CC(0), 
    tol_pivot_sym(pivot_sym),tol_pivot(pivot)
  { 
     SuperMatrix    B, X;
     SuperLUStat_t stat;
     void           *work=0;
     int            info, lwork=0, nrhs=1;
     int            i;
     double         ferr[1];
     double         berr[1];
     double          rpg, rcond;
    
     R *bb;
     R *xx;

     A.Store=0;
     B.Store=0;
     X.Store=0;
     L.Store=0;
     U.Store=0;
	
    int status;

    n=AA.n;
    m=AA.m;
    nnz=AA.nbcoef;

    arow=AA.a;
    asubrow=AA.cl;
    xarow=AA.lg;

    /* FreeFem++ use Morse Format */ 
    // FFCS - "this->" required by g++ 4.7
    this->CompRow_to_CompCol(m, n, nnz, arow, asubrow, xarow, 
		       &a, &asub, &xa);

    /* Defaults */
    lwork = 0;
    nrhs  = 0;
    
    /* Set the default values for options argument:
	options.Fact = DOFACT;
	options.Equil = YES;
	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
	options.Trans = NOTRANS;
	options.IterRefine = NOREFINE;
	options.SymmetricMode = NO;
	options.PivotGrowth = NO;
	options.ConditionNumber = NO;
	options.PrintStat = YES;
    */
    set_default_options(&options);
    
    printf(".. default options:\n");
    printf("\tFact\t %8d\n", options.Fact);
    printf("\tEquil\t %8d\n", options.Equil);
    printf("\tColPerm\t %8d\n", options.ColPerm);
    printf("\tDiagPivotThresh %8.4f\n", options.DiagPivotThresh);
    printf("\tTrans\t %8d\n", options.Trans);
    printf("\tIterRefine\t%4d\n", options.IterRefine);
    printf("\tSymmetricMode\t%4d\n", options.SymmetricMode);
    printf("\tPivotGrowth\t%4d\n", options.PivotGrowth);
    printf("\tConditionNumber\t%4d\n", options.ConditionNumber);
    printf("..\n");
    
    if(!string_option.empty()) read_options_freefem(string_option,&options);
    
    printf(".. options:\n");
    printf("\tFact\t %8d\n", options.Fact);
    printf("\tEquil\t %8d\n", options.Equil);
    printf("\tColPerm\t %8d\n", options.ColPerm);
    printf("\tDiagPivotThresh %8.4f\n", options.DiagPivotThresh);
    printf("\tTrans\t %8d\n", options.Trans);
    printf("\tIterRefine\t%4d\n", options.IterRefine);
    printf("\tSymmetricMode\t%4d\n", options.SymmetricMode);
    printf("\tPivotGrowth\t%4d\n", options.PivotGrowth);
    printf("\tConditionNumber\t%4d\n", options.ConditionNumber);
    printf("..\n");

    Dtype_t R_SLU = SuperLUDriver<R>::R_SLU_T(); 

    // FFCS - "this->" required by g++ 4.7
    this->Create_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, R_SLU, SLU_GE);
  
    this->Create_Dense_Matrix(&B, m, 0, (R*) 0, m, SLU_DN, R_SLU, SLU_GE);
    this->Create_Dense_Matrix(&X, m, 0, (R*) 0, m, SLU_DN, R_SLU, SLU_GE);
      

      if ( etree.size() ==0 )   etree.resize(n);   
      if ( perm_r.size() ==0 )   perm_r.resize(n);   
      if ( perm_c.size() ==0 )   perm_c.resize(n);   
     
    if ( !(RR = new double[n]) )
        ABORT("SUPERLU_MALLOC fails for R[].");
    for(int ii=0; ii<n; ii++){
      RR[ii]=1.;
    }
    if ( !(CC = new double[m]) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    for(int ii=0; ii<n; ii++){
      CC[ii]=1.;
    }    
    ferr[0]=0;
    berr[0]=0;
    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    /* ONLY PERFORM THE LU DECOMPOSITION */
    B.ncol = 0;  /* Indicate not to solve the system */
    SuperLUDriver<R>::gssvx(&options, &A, perm_c, perm_r, etree, equed, RR, CC,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &mem_usage, &stat, &info);

   

    if(verbosity>2)
    printf("LU factorization: dgssvx() returns info %d\n", info);
    if(verbosity>3)
    {
    if ( info == 0 || info == n+1 ) {
	
	if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", rpg);
	if ( options.ConditionNumber )
	    printf("Recip. condition number = %e\n", rcond);
        Lstore = (SCformat *) L.Store;
        Ustore = (NCformat *) U.Store;
	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
	       stat.expansions
	       );
	fflush(stdout);
	
    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
    }
    }
    if ( verbosity>5 ) StatPrint(&stat);
    StatFree(&stat);
    if( B.Store)  Destroy_SuperMatrix_Store(&B);
    if( X.Store)  Destroy_SuperMatrix_Store(&X);
    options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */


  }
  void Solver(const MatriceMorse<R> &AA,KN_<R> &x,const KN_<R> &b) const  {
    SuperMatrix    B, X;
    SuperLUStat_t stat;
    void           *work=0;
    int            info=0, lwork=0, nrhs=1;
    int            i;
    double       ferr[1], berr[1];
    double         rpg, rcond;
    double       *xx;

    B.Store=0;
    X.Store=0;
    ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    Dtype_t R_SLU = SuperLUDriver<R>::R_SLU_T(); 

      { 
	  KN_2Ptr<R> xx(x),bb(b);
	  // cout << " xx #### " << xx.c.N() << " "<< xx.ca.N() <<  " " << xx.ca.step << endl;
	  //cout << " bb #### " << bb.c.N() << " "<< bb.ca.N() << " " << bb.ca.step <<endl;
	  // FFCS - "this->" required by g++ 4.7
	  this->Create_Dense_Matrix(&B, m, 1, bb, m, SLU_DN, R_SLU, SLU_GE);
	  this->Create_Dense_Matrix(&X, m, 1, xx, m, SLU_DN, R_SLU, SLU_GE);
	  
	  B.ncol = nrhs;  /* Set the number of right-hand side */
	  
	  /* Initialize the statistics variables. */
	  StatInit(&stat);
	  
	  
	  SuperLUDriver<R>::gssvx(&options, &A, perm_c, perm_r, etree, equed, RR, CC,
				  &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
				  &mem_usage, &stat, &info);
	  
	  
	  
	  if(verbosity>2)
	      printf("Triangular solve: dgssvx() returns info %d\n", info);
	  
      }
  
   

    if(verbosity>3)
    {
    if ( info == 0 || info == n+1 ) {
	
        /* This is how you could access the solution matrix. */
        R *sol = (R*) ((DNformat*) X.Store)->nzval; 
	
	
	if ( options.IterRefine ) {
            printf("Iterative Refinement:\n");
	    printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
	    printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr[0], berr[0]);
	}
	fflush(stdout);
    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
    }
    }
    

    //cout << "   x min max " << x.min() << " " <<x.max() << endl;
    //cout << "=========================================" << endl;
    if( B.Store)  Destroy_SuperMatrix_Store(&B);
    if( X.Store)  Destroy_SuperMatrix_Store(&X);
  }

  ~SolveSuperLU() { 
   if(verbosity>3)
       cout << "~SolveSuperLU S:" << endl;
   //   if (etree)    delete[] etree; 
   //   if (perm_r)   delete[] perm_r; 
   //   if (perm_c)   delete[] perm_c; 
      if (RR)   delete[] RR; 
      if (CC)   delete[] CC; 
      if( A.Store)  Destroy_SuperMatrix_Store(&A);
      if( L.Store)  Destroy_SuperNode_Matrix(&L);
      if( U.Store)  Destroy_CompCol_Matrix(&U);
      
  }
  void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<R> &) (*this) * x; 
  }
     
}; 

MatriceMorse<double>::VirtualSolver *
BuildSolverSuperLU(DCL_ARG_SPARSE_SOLVER(double,A))
{
    if(verbosity>9)
    cout << " BuildSolverSuperLU<double>" << endl;
    return new SolveSuperLU<double>(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym,ds.sparams,ds.perm_r,ds.perm_c);
}

MatriceMorse<Complex>::VirtualSolver *
BuildSolverSuperLU(DCL_ARG_SPARSE_SOLVER(Complex,A))
{
  if(verbosity>9)
    cout << " BuildSolverSuperLU<Complex>" << endl;
  return new SolveSuperLU<Complex>(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym,ds.sparams,ds.perm_r,ds.perm_c);
}


/*  class Init { public:
    Init();
    };*/

//  the 2 default sparse solver double and complex
DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; ;
DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetSuperLU()
{
    if(verbosity>1)
	cout << " SetDefault sparse solver to SuperLU" << endl;
    DefSparseSolver<double>::solver  =BuildSolverSuperLU;
    DefSparseSolver<Complex>::solver =BuildSolverSuperLU;    
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
    return  true;
}

static void Load_Init()
{ 
  
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  
  if(verbosity>1)
    cout << "\n Add: SuperLU,  defaultsolverSuperLU" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  DefSparseSolver<double>::solver =BuildSolverSuperLU;
  DefSparseSolver<Complex>::solver =BuildSolverSuperLU;
  Global.Add("defaulttoSuperLU","(",new OneOperator0<bool>(SetSuperLU));
}

LOADFUNC(Load_Init)
