//   tools to read ppm file 
/*  use in freefem++ edp
  see :
  real[int,int] ff1("tt.pmm"); // read  image and set to an array. 
  real[int]  ff(ff1.nx*ff1.ny);
  ff=ff1; 
 */
#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;  

#include "RNM.hpp"
#include <cmath>
typedef KNM<double> * pRnm;
typedef KN<double> * pRn;
typedef string * pstring;
#include "ppmimg.h"

#ifdef __cplusplus
extern "C" {
#endif
    
    
    
    PPMimage *loadPPM(const char *imgname,ubyte *type,ubyte quiet) {
	pPPMimage  result;
	FILE      *fp;
	int        i,k,typimg,ret,r,g,b,s,maxval,bitsize;
	const char      *ptr;
	char c,buff[1024],data[256];
	
	/* search for image */
	fprintf(stdout," Loading image: %s\n",imgname);
	ptr = strstr(imgname,".ppm");
	strcpy(data,imgname);
	if ( !ptr ) {
	    ptr = strstr(imgname,".pgm");
	    if ( !ptr )  strcat(data,".ppm");
	    fp = fopen(data,"rb");
	} 
	else
	    fp = fopen(data,"rb");
	if ( !fp ) {
	    fprintf(stderr,"  ## UNABLE TO OPEN FILE %s.\n",data);
	    return(0);
	}
	if ( !quiet )
	    fprintf(stdout,"  opening %s\n",data);
	
	if ( !fgets(buff,sizeof(buff),fp) ) {
	    fprintf(stderr,"  ## INVALID HEADER.\n");
	    return(0);
	}
	
	/* check header file */
	if ( buff[0] != 'P' ) {
	    fprintf(stderr,"  ## INVALID IMAGE FORMAT (MUST BE 'PX').\n");
	    return(0);
	}
	
	switch(buff[1]) {
	    case '2': typimg = P2;  break;
	    case '3': typimg = P3;  break;
	    case '4': typimg = P4;  break;
	    case '5': typimg = P5;  break;
	    case '6': typimg = P6;  break;
	    default:
		fprintf(stderr,"  ## INVALID IMAGE FORMAT (MUST BE 'PX').\n");
		return(0);
	}
	
	/* allocate memory to store imagee */
	result = (PPMimage*) malloc(sizeof(PPMimage));
	assert(result);
	
	do {
	    ret = fscanf(fp,"%s",buff);
	    if ( ret == EOF ) break;
	    /* check and strip comments */
	    if ( buff[0] == '#' )
		do
		    c = getc(fp);
	    while ( c != '\n' );
	    else break;
	}
	while (1);
	
	/* read columns + lines */
	ret  = sscanf(buff,"%d",&s);
	result->sizeX = (short)s;
	ret += fscanf(fp,"%d",&s);
	result->sizeY = (short)s;
	if ( ret != 2 ) {
	    fprintf(stderr,"  ## ERROR LOADING IMAGE.\n");
	    free(result);
	    return(0);
	}
	if ( fscanf(fp,"%d",&maxval) != 1 ) {
	    fprintf(stderr,"  ## INVALID IMAGE SIZE.\n");
	    free(result);
	    return(0);
	}
	
	/* strip line */
	while ( fgetc(fp) != '\n' ) ;
	
	/* size based on type */
	if ( typimg == P2 || typimg == P5 || typimg == P4 )
	    bitsize = result->sizeX*result->sizeY;
	else
	    bitsize = 3*result->sizeX*result->sizeY;
	if ( !quiet )
	    fprintf(stdout,"   image size: %dx%d  %d bytes\n",
		    result->sizeX,result->sizeY,bitsize);
	
	result->data = (ubyte*)malloc(1+bitsize*sizeof(ubyte));
	assert(result->data);
	
	/* read data file */
	switch( typimg ) {
	    case P2:  /* ascii file (grey)  */
	    case P3:  /* ascii file (color) */
		for (i=0; i<bitsize; i++) {
		    int rr=fscanf(fp,"%d",&r);
		    result->data[i] = (ubyte)r;
		}
		break;
		
	    case P5:  /* binary file (grey) */
	    case P6:  /* binary file (color) */
		ret = fread(result->data,sizeof(ubyte),bitsize,fp);
		if ( ret != bitsize ) {
		    fprintf(stderr,"  ## ERROR LOADING IMAGE.\n");
		    free(result->data);
		    free(result);
		    return(0);
		}
		break;
	}
	fclose(fp);
	
	if ( *type == DEFAULT )
	    if ( typimg == P2 || typimg == P5 )
		*type = GREY;
	    else
		*type = COLOR;
	
	/* convert to grey levels */
	    else if ( *type == GREY && (typimg == P3 || typimg == P6) ) {
		fprintf(stdout,"  converting to grey levels\n");
		for (i=0,k=0; i<bitsize; i+=3,k++) {
		    r = (int)result->data[i];
		    g = (int)result->data[i+1];
		    b = (int)result->data[i+2];
		    result->data[k] = (ubyte)(0.3*r+0.59*g+0.11*b);
		}
		result->data = (ubyte*)realloc(result->data,sizeof(ubyte)*bitsize/3+1);
	    }
	
	return(result);
    }
    
    
    int savePPM(const char *imgname,pPPMimage img,int typimg) {
	FILE      *out;
	int        i,c,bitsize;
	
	/* open file */
	out = fopen(imgname,"w");
	if ( !out ) {
	    fprintf(stderr,"  ## UNABLE TO OPEN FILE %s.\n",imgname);
	    return 0;
	}
	
	/* write out image file */
	bitsize = img->sizeX*img->sizeY;
	switch(typimg) {
	    case P2:
		fprintf(out,"P2\n");
		fprintf(out,"# CREATOR: QIZIP Version 1, Rev. 2/2003, (c) INRIA\n");
		fprintf(out,"%d %d\n",img->sizeX,img->sizeY);
		fprintf(out,"255\n");
		c = 0;
		for (i=0; i<img->sizeX*img->sizeY; i++) {
		    fprintf(out,"%3d ",(int)img->data[i]);
		    if ( ++c == 17 ) { 
			c = 0; 
			fprintf(out,"\n");
		    }
		}
		fprintf(out,"\n");
		break;
	    case P5:
		fprintf(out,"P5\n");
		fprintf(out,"# CREATOR: QIZIP Version 1, Rev. 2/2003, (c) INRIA\n");
		fprintf(out,"%d %d\n",img->sizeX,img->sizeY);
		fprintf(out,"255\n");
		fwrite(img->data,sizeof(ubyte),bitsize,out);
		break;
	    case P6:
		fprintf(out,"P6\n");
		fprintf(out,"# CREATOR: QIZIP Version 1, Rev. 2/2003, (c) INRIA\n");
		fprintf(out,"%d %d\n",img->sizeX,img->sizeY);
		fprintf(out,"255\n");
		fwrite(img->data,sizeof(ubyte),3*bitsize,out);
		break;
	}
	fclose(out);
	
	return(1);
    }
    
    /* compute difference image */
    pPPMimage diffImg(pPPMimage bits,pPPMimage img,ubyte itype) {
	pPPMimage  dif;
	double     psnr,dd;
	int        i,bitsize,dmax;
	
	fprintf(stdout,"  Difference image\n");
	bitsize = (int)bits->sizeX*bits->sizeY;
	if ( itype == COLOR )  bitsize *= 3;
	
	dif = (PPMimage *)malloc(sizeof(PPMimage));
	if ( !dif ) {
	    fprintf(stderr,"  Sorry, not enough memory. Bye.\n");
	    return 0;
	}
	dif->sizeX = bits->sizeX; 
	dif->sizeY = bits->sizeY;
	dif->data = (ubyte*)malloc(bitsize*sizeof(ubyte));
	if ( !dif->data ) {
	    fprintf(stderr,"  Sorry, not enough memory. Bye.\n");
	    free(dif);
	    return 0;
	}
	
	dmax = 0;
	psnr = 0.0f;
	for (i=0; i<bitsize; i++) {
	    dd    = abs((int)(bits->data[i]-img->data[i]));
	    dmax  = max(dmax,dd);
	    psnr += (double)dd*dd;
	    dif->data[i] = (ubyte)(255-dd);
	}
	if ( psnr == 0.0f )  fprintf(stderr,"    PSNR problem!");
	else {
	    psnr = 65025.0f / psnr;
	    psnr = 10.0 * log10(bitsize*psnr);
	}
	fprintf(stdout,"    PSNR = %.2f    dmax = %d\n",psnr,dmax);
	
	return(dif);
    }
    
    
#ifdef __cplusplus
}
#endif


pRnm  read_image( pRnm  const & a,const pstring & b)
{
  ubyte type,quiet=1;
  PPMimage * image =loadPPM(b->c_str(),&type, quiet);
  if(!image) {
    std::cerr << " error loadPPM image "<< *b  << endl;
    CompileError("error loadPPM image ");
    return a;
  }
  if(verbosity)   
  cout << " size of image : " << image->sizeX << " x " << image->sizeY << " type =" <<  (int) type << endl; 
  int n = image->sizeX;
  int m = image->sizeY ;
  a->init(n,m);
  ubyte * dd= image->data;

  //  cout << (double) dd[0] / 256. << " "
  //     << (double) dd[250] / 256. << " "
  //     << (double) dd[500] / 256. << "\n "
  //  ;
  int k=0;
  double *mm=*a;
  for(int i=0;i<n;++i)
    for(int j=0;j<m;++j)
      *mm++= (double) dd[k++] / 256. ;
  KN_<double> aa=*a;
  // cout << aa[0] << " "<< aa[250] << "" << aa[500] << endl;
  assert(k==n*m);
  free(image->data);
  free(image);
  return a;
}
pRn  seta( pRn  const & a,const pRnm & b)
{
  *a=*b;
  KN_<double> aa=*a;
  //  cout << aa[0] << " "<< aa[250] << "" << aa[500] << endl;
  return a;
}
class Init { public:
  Init();
};

LOADINIT(Init);
Init::Init(){
  cout << " lood: init ppm2rmn  " << endl;


  TheOperators->Add("<-", 
		    new OneOperator2_<KNM<double> *,KNM<double> *,string*>(&read_image)
		    );
  TheOperators->Add("=", 
		    new OneOperator2_<KN<double> *,KN<double> *,KNM<double>* >(seta)
		    );
  /*
  map_type[typeid(KN<double> ).name()]->AddCast(
						new E_F1_funcT<KN<double>,KNM<double>*>(UnRef<KN<double>,KNM<double> >));
					      //  map_type[typeid(KN<double> ).name()]->AddCast(
						//new E_F1_funcT<KN<double>*,KNM<double>*>(Cast<KN<double>*,KNM<double>*>));
					       
						*/
}
