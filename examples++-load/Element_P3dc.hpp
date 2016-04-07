 const int TypeOfFE_P3dcLagrange::nn[10][3] =  { 
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
 const int TypeOfFE_P3dcLagrange::aa[10][3] =  { 
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
 const int TypeOfFE_P3dcLagrange::ff[10] =  { 6 ,6 ,6 ,2 ,2 ,2 ,2 ,2 ,2 ,1 };
 const int TypeOfFE_P3dcLagrange::il[10] =  { 3 ,0 ,0 ,0 ,0 ,1 ,2 ,2 ,1 ,1 };
 const int TypeOfFE_P3dcLagrange::jl[10] =  { 0 ,3 ,0 ,2 ,1 ,0 ,0 ,1 ,2 ,1 };
 const int TypeOfFE_P3dcLagrange::kl[10] =  { 0 ,0 ,3 ,1 ,2 ,2 ,1 ,0 ,0 ,1 };
