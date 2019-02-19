#ifndef WITH_NO_INIT
#include "ff++.hpp"
#endif
#include "AFunction_ext.hpp"

using namespace std;

#include <fstream>

class Marmousi2d{
public:
  string * file;
  KNM<float>* tab;
  
  int mo_file_nx, mo_file_ny;
  double mo_file_xend, mo_file_xstart;
  double mo_file_yend, mo_file_ystart;
  double mo_file_dx, mo_file_dy;
  double cfill;
  
  void init() { file=0; tab=0;}
  void destroy()
  {
    delete file; 
    delete tab;
  } 
};

Marmousi2d * init_Marmousi2d(Marmousi2d * const &a, string * const & s)
{
	if (verbosity)
	  cout << "Reading Marmousi2d Model file " << *s << endl;
    a->file = new string(* s);

        ifstream f((*a->file).c_str(), ios::in | ios::binary);
        if (!f.is_open()) {
          cerr << "Error opening " << (*a->file).c_str() << ": file does not exist." << endl;
          ffassert(f.is_open());
        } 
	
	int sz = 2301*751;
	a->mo_file_nx = 2301;
	a->mo_file_xstart = 0.;
    a->mo_file_xend = 9.2;
    a->mo_file_ny = 751;
    a->mo_file_ystart = -3.;
    a->mo_file_yend = 0.;
	
	float* buff = new float[sz];
	
	f.read((char*) buff, sz*sizeof(float));
	
	f.close();
	
	a->tab = new KNM<float>(2301,751);
		
	int ix,iy,iz;
	
	for (iy = 0; iy < a->mo_file_ny; iy++)
    for (ix = 0; ix < a->mo_file_nx; ix++)
      (*a->tab)(ix,a->mo_file_ny-1-iy) = buff[a->mo_file_ny*ix + iy];
      
    delete [] buff;
         	
	f.close();
} 

double Marmousi2d_eval(Marmousi2d * const & a,const double & xi,const double & yi)
{  
  int ix = a->mo_file_nx*(xi-a->mo_file_xstart)/(a->mo_file_xend-a->mo_file_xstart);
  int iy = a->mo_file_ny*(yi-a->mo_file_ystart)/(a->mo_file_yend-a->mo_file_ystart);
  ix = max(0,min(ix,a->mo_file_nx-1));
  iy = max(0,min(iy,a->mo_file_ny-1));
  
  return (*a->tab)(ix,iy);
}

static void Load_Init()
{ 
  if (verbosity && mpirank == 0)
    cout << " load: Marmousi  " << endl;
    
  Dcl_Type<Marmousi2d *>(InitP<Marmousi2d>,Destroy<Marmousi2d>);
  zzzfff->Add("Marmousi",atype<Marmousi2d*>());
  TheOperators->Add("<-", new OneOperator2_<Marmousi2d *,Marmousi2d* ,string*>(&init_Marmousi2d));
  atype< Marmousi2d * >()->Add("(","",new OneOperator3_<double,Marmousi2d *,double,double>(Marmousi2d_eval));
}
LOADFUNC(Load_Init)
