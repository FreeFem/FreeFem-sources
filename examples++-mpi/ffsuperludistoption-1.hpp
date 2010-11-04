// read options for superlu in freefem++
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

void read_nprow_npcol_freefem(string *string_option, int *nprow, int *npcol){
  
  static const char* comp[] = {"Fact","Equil","ParSymbFact","ColPerm","RowPerm",
			       "DiagPivotThresh","IterRefine","Trans",
			       "ReplaceTinyPivot","SolveInitialized",
			       "RefineInitialized","PrintStat","nprow","npcol",0};

  char data[string_option->size()+1];  
  strcpy( data, string_option->c_str()); 
  char *tictac;
  char *tictac2;
  tictac = strtok(data," =,\t\n");
  
  while(tictac != NULL){
    int id_option = s_(tictac, comp);
    tictac2 = tictac;
    tictac = strtok(NULL," =,\t\n");
    int val_options;

    switch (id_option)
      { 
      case 13: // nprow
	*nprow = atoi(tictac);
	break;
      case 14: // npcol 
	*npcol = atoi(tictac);
	break;
      default: // Equivalent of case default
	if(id_option == 0)
	  {
	    printf("parameter is not valid for superlu_dist %s \n", tictac2 );
	    exit(1);
	  }	  
	break;
      }  
    tictac = strtok(NULL," =,\t\n");
  }
}

void read_nprow_npcol_freefem(string *string_option, int *nprow, int *npcol, int *matrixdist){
  
  static const char* comp[] = {"Fact","Equil","ParSymbFact","ColPerm","RowPerm",
			       "DiagPivotThresh","IterRefine","Trans",
			       "ReplaceTinyPivot","SolveInitialized",
			       "RefineInitialized","PrintStat","nprow","npcol","matrix",0};

  char data[string_option->size()+1];  
  strcpy( data, string_option->c_str()); 
  char *tictac;
  char *tictac2;
  tictac = strtok(data," =,\t\n");
  
  while(tictac != NULL){
    int id_option = s_(tictac, comp);
    tictac2 = tictac;
    tictac = strtok(NULL," =,\t\n");
    int val_options;
    printf("param %s = value %s , id_option %d\n",tictac2,tictac,id_option);
    switch (id_option)
      { 
      case 13: // nprow
	*nprow = atoi(tictac);
	break;
      case 14: // npcol 
	*npcol = atoi(tictac);
	break;
      case 15: // matrix
	printf("parameter matrix \n");
	if(strcmp(tictac,"assembled") == 0)
	  *matrixdist = 0;
	else if(strcmp(tictac,"distributedglobal") == 0) 
	  *matrixdist = 1;
	else if(strcmp(tictac,"distributed") == 0) 
	  *matrixdist = 2;
	else{
	  printf("value of parameter matrix is not correct %s \n", tictac );
	}
	break;
      default: // Equivalent of case default
	if(id_option == 0)
	  {
	    printf("parameter is not valid for superlu_dist %s \n", tictac2 );
	    exit(1);
	  }	  
	break;
      }  
    tictac = strtok(NULL," =,\t\n");
  }
}


void read_options_freefem(string *string_option, superlu_options_t *options, DiagScale_t *diag){
  static const yes_no_t  enumyes_no_t[2] = {NO, YES};
  static const fact_t  enumfact_t[4] = {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED};
  static const colperm_t  enumcolperm_t[6] = {NATURAL, MMD_AT_PLUS_A, MMD_ATA, METIS_AT_PLUS_A,PARMETIS, MY_PERMC};
  static const rowperm_t  enumrowperm_t[3] = {NOROWPERM, LargeDiag, MY_PERMR};
  static const DiagScale_t enumDiagScale_t[4] = {NOEQUIL, ROW, COL, BOTH};
  static const trans_t  enumtrans_t[3] = {NOTRANS, TRANS, CONJ};
  static const IterRefine_t enumIterRefine_t[4] = {NOREFINE, SINGLE, DOUBLE, EXTRA};  
  //static const MemType enumMemType_t[4] = {LUSUP, UCOL, LSUB, USUB};  
  //static const stack_end_t enumstack_end_t[2] = {HEAD, TAIL};
  //static const LU_space_t enumLU_space_t[2] = {SYSTEM, USER};

 
  static const char* compyes_no_t[] = {"NO", "YES",0};
  static const char* compfact_t[] = {"DOFACT", "SamePattern", "SamePattern_SameRowPerm", "FACTORED",0};
  static const char* comprowperm_t[] = {"NOROWPERM", "LargeDiag", "MY_PERMR",0};
  static const char* compcolperm_t[] = {"NATURAL", "MMD_AT_PLUS_A", "MMD_ATA", "METIS_AT_PLUS_A", "PARMETIS", "MY_PERMC",0};
  static const char* compDiagScale_t[] = {"NOEQUIL", "ROW", "COL", "BOTH",0};
  static const char* comptrans_t[] = {"NOTRANS", "TRANS", "CONJ",0};
  static const char* compIterRefine_t[] = {"NOREFINE", "SINGLE", "DOUBLE", "EXTRA",0};
  //static const char* compMemType_t[] = {"LUSUP", "UCOL", "LSUB", "USUB",0};  
  //static const char* compstack_end_t[] = {"HEAD", "TAIL",0};
  //static const char* compLU_space_t[] = {"SYSTEM", "USER",0};

  static const char* comp[] = {"Fact","Equil","ParSymbFact","ColPerm","RowPerm",
			       "DiagPivotThresh","IterRefine","Trans",
			       "ReplaceTinyPivot","SolveInitialized",
			       "RefineInitialized","PrintStat","nprow","npcol","DiagScale","matrix",0};

  char data[string_option->size()+1];  
  strcpy( data, string_option->c_str()); 
  cout << "data=" << data << endl;
  char *tictac;
  char *tictac2;
  tictac = strtok(data," =,\t\n");

  while(tictac != NULL){
    int id_option = s_(tictac, comp);
    tictac2=tictac;
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
      case 3:  // ParSymbFact
	//char* comp3[] = {"NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD", "MY_PERMC", 0};
	val_options= s_(tictac,compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
	  exit(1);
	}
	options->ParSymbFact = enumyes_no_t[val_options-1];
	break;
      case 4:  // ColPerm
	//char* comp3[] = {"NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD", "MY_PERMC", 0};
	val_options= s_(tictac,compcolperm_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
	  exit(1);
	}
	options->ColPerm = enumcolperm_t[val_options-1];
	break;
      case 5:  // RowPerm
	//char* comp3[] = {"NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD", "MY_PERMC", 0};
	val_options= s_(tictac,comprowperm_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
	  exit(1);
	}
	options->RowPerm = enumrowperm_t[val_options-1];
	break;
      case 6:  // DiagPivotThresh
	options->DiagPivotThresh= strtod(tictac,&tictac);
	break;
      case 7:  // IterRefine
	val_options= s_(tictac,compIterRefine_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
	  exit(1);
	}
	options->IterRefine = enumIterRefine_t[val_options-1];
	break;
      case 8:  // Trans
	//char* comp5[] = {"NOTRANS", "TRANS", "CONJ", 0};
	val_options= s_(tictac, comptrans_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","Trans");
	  exit(1);
	}
	options->Trans = enumtrans_t[val_options-1];
	break;
      case 9:  // ReplaceTinyPivot
	//char* comp7[] = {"NO","YES", 0};
	val_options= s_(tictac, compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","SymmetricMode");
	  exit(1);
	}
	options->ReplaceTinyPivot= enumyes_no_t[val_options-1]; 
	break;
      case 10:  // SolveInitialized
	//char* comp8[] = {"NO","YES", 0};
	val_options= s_(tictac,compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","PivotGrowth");
	  exit(1);
	}
	options->SolveInitialized = enumyes_no_t[val_options-1];
	break;
      case 11:  // RefineInitialized
	//char* comp9[] = {"NO","YES", 0};
	val_options= s_(tictac, compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ConditionNumber");
	  exit(1);
	}
	options->RefineInitialized = enumyes_no_t[val_options-1];
	break;
      case 12: // PrintStat
	val_options= s_(tictac, compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","PrintStat");
	  exit(1);
	}
	options->PrintStat = enumyes_no_t[val_options-1];
	break;
	// case 13 nprow
	// case 14 npcol
      case 15: // DiagScale_t
	val_options= s_(tictac, compDiagScale_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","PrintStat");
	  exit(1);
	}
        *diag = enumDiagScale_t[val_options-1];
	break;
      default: // Equivalent of case default
	if(id_option == 0)
	  {
	    printf("parameter is not valid for superlu_dist %s \n", tictac2);
	    exit(1);
	  }	 
	break;
      }  
    tictac = strtok(NULL," =,\t\n");
  }
  
}

// void read_nprow_npcol_freefem(string *string_option, int *nprow, int *npcol, int *matrixdist){
  
//   static const char* comp[] = {"Fact","Equil","ParSymbFact","ColPerm","RowPerm",
// 			       "DiagPivotThresh","IterRefine","Trans",
// 			       "ReplaceTinyPivot","SolveInitialized",
// 			       "RefineInitialized","PrintStat","nprow","npcol","matrix",0};

//   char data[string_option->size()+1];  
//   strcpy( data, string_option->c_str()); 
//   char *tictac;
//   char *tictac2;
//   tictac = strtok(data," =,\t\n");
  
//   while(tictac != NULL){
//     int id_option = s_(tictac, comp);
//     tictac2 = tictac;
//     tictac = strtok(NULL," =,\t\n");
//     int val_options;
//     printf("param %s = value %s , id_option %d\n",tictac2,tictac,id_option);
//     switch (id_option)
//       { 
//       case 13: // nprow
// 	*nprow = atoi(tictac);
// 	break;
//       case 14: // npcol 
// 	*npcol = atoi(tictac);
// 	break;
//       case 15: // matrix
// 	printf("parameter matrix \n");
// 	if(strcmp(tictac,"assembled") == 0)
// 	  *matrixdist = 0;
// 	else if(strcmp(tictac,"distributedglobal") == 0) 
// 	  *matrixdist = 1;
// 	else if(strcmp(tictac,"distributed") == 0) 
// 	  *matrixdist = 2;
// 	else{
// 	  printf("value of parameter matrix is not correct %s \n", tictac );
// 	}
// 	break;
//       default: // Equivalent of case default
// 	if(id_option == 0)
// 	  {
// 	    printf("parameter is not valid for superlu_dist %s \n", tictac2 );
// 	    exit(1);
// 	  }	  
// 	break;
//       }  
//     tictac = strtok(NULL," =,\t\n");
//   }
// }

void read_options_superlu_datafile(string *data_option, superlu_options_t *options, int_t *nprow, int_t *npcol, int *matrixdist, DiagScale_t *diag){
  static const yes_no_t  enumyes_no_t[2] = {NO, YES};
  static const fact_t  enumfact_t[4] = {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED};
  static const colperm_t  enumcolperm_t[6] = {NATURAL, MMD_AT_PLUS_A, MMD_ATA, METIS_AT_PLUS_A,PARMETIS, MY_PERMC};
  static const rowperm_t  enumrowperm_t[3] = {NOROWPERM, LargeDiag, MY_PERMR};
  static const DiagScale_t enumDiagScale_t[4] = {NOEQUIL, ROW, COL, BOTH};
  static const trans_t  enumtrans_t[3] = {NOTRANS, TRANS, CONJ};
  static const IterRefine_t enumIterRefine_t[4] = {NOREFINE, SINGLE, DOUBLE, EXTRA};  
  //static const MemType enumMemType_t[4] = {LUSUP, UCOL, LSUB, USUB};  
  //static const stack_end_t enumstack_end_t[2] = {HEAD, TAIL};
  //static const LU_space_t enumLU_space_t[2] = {SYSTEM, USER};

 
  static const char* compyes_no_t[] = {"NO", "YES",0};
  static const char* compfact_t[] = {"DOFACT", "SamePattern", "SamePattern_SameRowPerm", "FACTORED",0};
  static const char* comprowperm_t[] = {"NOROWPERM", "LargeDiag", "MY_PERMR",0};
  static const char* compcolperm_t[] = {"NATURAL", "MMD_AT_PLUS_A", "MMD_ATA", "METIS_AT_PLUS_A", "PARMETIS", "MY_PERMC",0};
  static const char* compDiagScale_t[] = {"NOEQUIL", "ROW", "COL", "BOTH",0};
  static const char* comptrans_t[] = {"NOTRANS", "TRANS", "CONJ",0};
  static const char* compIterRefine_t[] = {"NOREFINE", "SINGLE", "DOUBLE", "EXTRA",0};
  
  //int_t ffnprow,ffnpcol;
  //int matrixdist;

  char datafile[data_option->size()+1];  
  strcpy( datafile, data_option->c_str()); 
  
  FILE* pfile= fopen( datafile,"rt");
  char data[256];
  char *tictac;
  
  fgets(data,256,pfile);
  cout << "data=" << data << endl;
  tictac = strtok(data," /!#\t\n");
  *nprow = (int) atol(tictac);
  if(verbosity) printf("nprow=%d\n",*nprow);

  fgets(data,256,pfile);
  tictac = strtok(data," /!#\t\n");
  *npcol = (int) atol(tictac);
  if(verbosity) printf("npcol=%d\n",*npcol);

  fgets(data,256,pfile);
  tictac = strtok(data," /!#\t\n");
  if(strcmp(tictac,"assembled") == 0)
    *matrixdist = 0;
  else if(strcmp(tictac,"distributedglobal") == 0) 
    *matrixdist = 1;
  else if(strcmp(tictac,"distributed") == 0) 
    *matrixdist = 2;
  else{
    printf("matrix input %s for superlu_dist is not correct\n", tictac );
    exit(1);
  }

  int id_option=0;
  
  while( !feof(pfile) && id_option<12){
    fgets(data,256,pfile);
    tictac = strtok(data," /!#\t\n");
    id_option++;
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
      case 3:  // ParSymbFact
	//char* comp3[] = {"NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD", "MY_PERMC", 0};
	val_options= s_(tictac,compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
	  exit(1);
	}
	options->ParSymbFact = enumyes_no_t[val_options-1];
	break;
      case 4:  // ColPerm
	val_options= s_(tictac,compcolperm_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
	  exit(1);
	}
	options->ColPerm = enumcolperm_t[val_options-1];
	break;
      case 5:  // RowPerm
	val_options= s_(tictac,comprowperm_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
	  exit(1);
	}
	options->RowPerm = enumrowperm_t[val_options-1];
	break;
      case 6:  // DiagPivotThresh
	options->DiagPivotThresh= strtod(tictac,&tictac);
	break;
      case 7:  // IterRefine
	val_options= s_(tictac,compIterRefine_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ColPerm");
	  exit(1);
	}
	options->IterRefine = enumIterRefine_t[val_options-1];
	break;
      case 8:  // Trans
	//char* comp5[] = {"NOTRANS", "TRANS", "CONJ", 0};
	val_options= s_(tictac, comptrans_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","Trans");
	  exit(1);
	}
	options->Trans = enumtrans_t[val_options-1];
	break;
      case 9:  // ReplaceTinyPivot
	//char* comp7[] = {"NO","YES", 0};
	val_options= s_(tictac, compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","SymmetricMode");
	  exit(1);
	}
	options->ReplaceTinyPivot= enumyes_no_t[val_options-1]; 
	break;
      case 10:  // SolveInitialized
	//char* comp8[] = {"NO","YES", 0};
	val_options= s_(tictac,compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","PivotGrowth");
	  exit(1);
	}
	options->SolveInitialized = enumyes_no_t[val_options-1];
	break;
      case 11:  // RefineInitialized
	//char* comp9[] = {"NO","YES", 0};
	val_options= s_(tictac, compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","ConditionNumber");
	  exit(1);
	}
	options->RefineInitialized = enumyes_no_t[val_options-1];
	break;
      case 12: // PrintStat
	val_options= s_(tictac, compyes_no_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","PrintStat");
	  exit(1);
	}
	options->PrintStat = enumyes_no_t[val_options-1];
	break;
      case 13: // DiagScale_t
	val_options= s_(tictac, compDiagScale_t);
	if( val_options == 0){
	  printf("value given for SuperLU for options %s is not correct\n","PrintStat");
	  exit(1);
	}
        *diag = enumDiagScale_t[val_options-1];
	break;
      default: // Equivalent of case default
	if(id_option == 0 && id_option > 13)
	  {
	    printf("Error in reading data file for superlu_dist %s\n",datafile);
	    exit(1);
	  }	 
	break;
      }  
  }
  fclose(pfile);
}
