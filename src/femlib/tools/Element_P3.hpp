static int TypeOfFE_P3Lagrange::nn[10][3] =  { 
		 { 0 ,0 ,0 } ,
		 { 1 ,1 ,1 } ,
		 { 2 ,2 ,2 } ,
		 { 1 ,1 ,2 } ,
		 { 1 ,2 ,2 } ,
		 { 0 ,2 ,2 } ,
		 { 0 ,0 ,2 } ,
		 { 0 ,0 ,1 } ,
		 { 0 ,1 ,1 } ,
		 { 0 ,1 ,2 } }
;
static int TypeOfFE_P3Lagrange::aa[10][3] =  { 
		 { 0 ,1 ,2 } ,
		 { 0 ,1 ,2 } ,
		 { 0 ,1 ,2 } ,
		 { 0 ,1 ,0 } ,
		 { 0 ,0 ,1 } ,
		 { 0 ,0 ,1 } ,
		 { 0 ,1 ,0 } ,
		 { 0 ,1 ,0 } ,
		 { 0 ,0 ,1 } ,
		 { 0 ,0 ,0 } }
;
static int TypeOfFE_P3Lagrange::ff[10] =  { 6 ,6 ,6 ,2 ,2 ,2 ,2 ,2 ,2 ,1 };
static int TypeOfFE_P3Lagrange::il[10] =  { 3 ,0 ,0 ,0 ,0 ,1 ,2 ,2 ,1 ,1 };
static int TypeOfFE_P3Lagrange::jl[10] =  { 0 ,3 ,0 ,2 ,1 ,0 ,0 ,1 ,2 ,1 };
static int TypeOfFE_P3Lagrange::kl[10] =  { 0 ,0 ,3 ,1 ,2 ,2 ,1 ,0 ,0 ,1 };
