#ifndef WITH_NO_INIT
#include "ff++.hpp"
#endif
#include "AFunction_ext.hpp"

using namespace std;

#include <fstream>

class Overthrust3d{
public:
  string * file;
  KNMK<float>* tab;
  
  int mo_file_nx, mo_file_ny, mo_file_nz;
  double mo_file_xend, mo_file_xstart;
  double mo_file_yend, mo_file_ystart;
  double mo_file_zend, mo_file_zstart;
  double mo_file_dx, mo_file_dy, mo_file_dz;
  double cfill;
  
  void init() { file=0; tab=0;}
  void destroy()
  {
    delete file; 
    delete tab;
  } 
};

Overthrust3d * init_Overthrust3d(Overthrust3d * const &a, string * const & s)
{
	if (verbosity)
	  cout << "Reading Overthrust3d Model file " << *s << endl;
    a->file = new string(* s);

	
	ifstream f((*a->file).c_str(), ios::in | ios::binary);
        if (!f.is_open()) {
          cerr << "Error opening " << (*a->file).c_str() << ": file does not exist." << endl;
          ffassert(f.is_open());
        }   
	
	int sz = 801*801*187;
	a->mo_file_nx = 801;
	a->mo_file_xstart = 0.;
    a->mo_file_xend = 20.;
    a->mo_file_ny = 801;
    a->mo_file_ystart = 0.;
    a->mo_file_yend = 20.;
    a->mo_file_nz = 187;
    a->mo_file_zstart = -4.65;
    a->mo_file_zend = 0.;
	
	float* buff = new float[sz];
	
	f.read((char*) buff, sz*sizeof(float));
	
	int ix,iy,iz;
	
	/*
	for (iz = 0; iz < a->mo_file_nz; iz++)
	for (iy = 0; iy < a->mo_file_ny; iy++)
	for (ix = 0; ix < a->mo_file_nx; ix++) {
	  f >> buff[a->mo_file_nx*a->mo_file_ny*iz + a->mo_file_nx*iy + ix];
	  buff[a->mo_file_nx*a->mo_file_ny*iz + a->mo_file_nx*iy + ix] *= 1e-3;
	}
	
	ofstream fout("3DOverthrustdata.bin", ios::out | ios::binary);
	
	fout.write((char*) buff, sz*sizeof(float));
	
	fout.close();
	*/
	
	f.close();
	
	a->tab = new KNMK<float>(801,801,187);

	for (iz = 0; iz < a->mo_file_nz; iz++)
	for (iy = 0; iy < a->mo_file_ny; iy++)
	for (ix = 0; ix < a->mo_file_nx; ix++)
      (*a->tab)(ix,iy,a->mo_file_nz-1-iz) = buff[a->mo_file_nx*a->mo_file_ny*iz + a->mo_file_nx*iy + ix];
      
    delete [] buff;
         	
	f.close();
} 

double Overthrust3d_eval(Overthrust3d * const & a,const double & xi,const double & yi,const double & zi)
{  
  int ix = a->mo_file_nx*(xi-a->mo_file_xstart)/(a->mo_file_xend-a->mo_file_xstart);
  int iy = a->mo_file_ny*(yi-a->mo_file_ystart)/(a->mo_file_yend-a->mo_file_ystart);
  int iz = a->mo_file_nz*(zi-a->mo_file_zstart)/(a->mo_file_zend-a->mo_file_zstart);
  ix = max(0,min(ix,a->mo_file_nx-1));
  iy = max(0,min(iy,a->mo_file_ny-1));
  iz = max(0,min(iz,a->mo_file_nz-1));
  
  return (*a->tab)(ix,iy,iz);
}

static void Load_Init()
{ 
  if (verbosity && mpirank == 0)
    cout << " load: Overthrust  " << endl;
    
  Dcl_Type<Overthrust3d *>(InitP<Overthrust3d>,Destroy<Overthrust3d>);
  zzzfff->Add("Overthrust",atype<Overthrust3d*>());
  TheOperators->Add("<-", new OneOperator2_<Overthrust3d *,Overthrust3d* ,string*>(&init_Overthrust3d));
  atype< Overthrust3d * >()->Add("(","",new OneOperator4_<double,Overthrust3d *,double,double,double>(Overthrust3d_eval));
}
LOADFUNC(Load_Init)
