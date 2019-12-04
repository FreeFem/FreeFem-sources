#ifndef WITH_NO_INIT
#include "ff++.hpp"
#endif
#include "AFunction_ext.hpp"

using namespace std;

#include <fstream>

class Marmousi22d {
 public:
  string *file;
  KNM< float > *tab;

  int mo_file_nx, mo_file_ny;
  double mo_file_xend, mo_file_xstart;
  double mo_file_yend, mo_file_ystart;
  double mo_file_dx, mo_file_dy;
  double cfill;

  void init( ) {
    file = 0;
    tab = 0;
  }
  void destroy( ) {
    delete file;
    delete tab;
  }
};

Marmousi22d *init_Marmousi22d(Marmousi22d *const &a, string *const &s) {
  if (verbosity) cout << "Reading Marmousi22d Model file " << *s << endl;
  a->file = new string(*s);

  ifstream f((*a->file).c_str( ), ios::in | ios::binary);
  if (!f.is_open( )) {
    cerr << "Error opening " << (*a->file).c_str( ) << ": file does not exist." << endl;
    ffassert(f.is_open( ));
  }

  int sz = 13601 * 2801;
  a->mo_file_nx = 13601;
  a->mo_file_xstart = 0.;
  a->mo_file_xend = 17.;
  a->mo_file_ny = 2801;
  a->mo_file_ystart = -3.5;
  a->mo_file_yend = 0.;

  float *buff = new float[sz];

  f.read((char *)buff, sz * sizeof(float));

  f.close( );

  a->tab = new KNM< float >(13601, 2801);

  int ix, iy, iz;

  for (iy = 0; iy < a->mo_file_ny; iy++)
    for (ix = 0; ix < a->mo_file_nx; ix++)
      (*a->tab)(ix, a->mo_file_ny - 1 - iy) = buff[a->mo_file_ny * ix + iy];

  delete[] buff;

  f.close( );
}

double Marmousi22d_eval(Marmousi22d *const &a, const double &xi, const double &yi) {
  int ix = a->mo_file_nx * (xi - a->mo_file_xstart) / (a->mo_file_xend - a->mo_file_xstart);
  int iy = a->mo_file_ny * (yi - a->mo_file_ystart) / (a->mo_file_yend - a->mo_file_ystart);
  ix = max(0, min(ix, a->mo_file_nx - 1));
  iy = max(0, min(iy, a->mo_file_ny - 1));

  return (*a->tab)(ix, iy);
}

static void Load_Init( ) {
  if (verbosity && mpirank == 0) cout << " load: Marmousi2  " << endl;

  Dcl_Type< Marmousi22d * >(InitP< Marmousi22d >, Destroy< Marmousi22d >);
  zzzfff->Add("Marmousi2", atype< Marmousi22d * >( ));
  TheOperators->Add("<-",
                    new OneOperator2_< Marmousi22d *, Marmousi22d *, string * >(&init_Marmousi22d));
  atype< Marmousi22d * >( )->Add(
    "(", "", new OneOperator3_< double, Marmousi22d *, double, double >(Marmousi22d_eval));
}
LOADFUNC(Load_Init)
