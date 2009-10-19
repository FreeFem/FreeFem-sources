//   tools to read ppm file 
/*  use in freefem++ edp
  see :
  real[int,int] ff1("tt.pmm"); // read  image and set to an array. 
  real[int]  ff(ff1.nx*ff1.ny);
  ff=ff1; 
 */
#include "pcm.hpp"
#include  <iostream>
#include  <cfloat>
// ff-c++ pcm*cpp -g
//ff-c++-cpp-dep:pcm.cpp
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;  

#include "RNM.hpp"
#include <cmath>


long read1(const long& ,const long&){
    return 1;
}

    
KNM<Complex> * read_pcm(string * filename,KNM<Complex> * p)
    {
      
      PCM pcm(filename->c_str());
      p->resize(pcm.width,pcm.height);
      pcm_complex *pc=pcm.image;
      for(int j=0;j<pcm.height;++j)
      for(int i=0;i<pcm.width;++i,pc++)
	      (*p)(i,j)= Complex(pc->r,pc->i);
	      

      return p;	  
    }
long read_pcm(string * const &filename,KNM<double> * const &u,KNM<double> * const &v)
{
    
    PCM pcm(filename->c_str());
    cout << " pcm  " << filename->c_str() << " : " <<pcm.width << " x " << pcm.height << endl;
    u->resize(pcm.width,pcm.height);
    v->resize(pcm.width,pcm.height);
    pcm_complex *pc;
    float x1=-1e+30,x2=-1e+30;
    for(int j=0;j<pcm.height;++j)
	for(int i=0;i<pcm.width;++i)
	  {
	    pc = pcm.Get(i,j);
	    if(pc)
	      {
		
	      
	    (*u)(i,j)= pc->r;
	    (*v)(i,j)= pc->i;
	    
	    x1 = max(x1,pc->r);
	    x2 = max(x2,pc->i);
	    if(i<0 && j < 0)
	    cout << i << " " << j << " " << pc->r << " " << pc->i << endl;
	       }
	  }
    cout << " max uv : " << x1 << " " << x2 << endl;
    return pcm.width*pcm.height;	  
}

 class Init { public:
  Init();
};

Init init;
Init::Init(){
  cout << " load: init pcm2rmn  " << endl;
 

    Global.Add("readpcm", "(",
    new OneOperator2<KNM<Complex> *,string*,KNM<Complex> * >(&read_pcm),
    new OneOperator3_<long,string*,KNM<double> *,KNM<double> * >(&read_pcm)
		    );
  
}
