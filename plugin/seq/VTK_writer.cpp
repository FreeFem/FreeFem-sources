/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Cedric Ody
// E-MAIL  : cedric.listes@gmail.com

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

// from the work of  Sala Lorenzo (Dxwriter)

#include "mode_open.hpp"
#include <iostream>
#include <cfloat>
#include <cmath>
#include <iterator>
using namespace std;
#include "ff++.hpp"
using namespace Fem2D;

class VtkWriter {
  struct tsinfo {
    int imesh;    //!< index of the mesh
    std::string name;
    std::vector< double > vecistant;
  };

 private:
  std::vector< const Fem2D::Mesh * > _vecmesh;
  std::string _nameoffile;

  /*! This string contains the name of data file with \\ where there's a \ in the path*/
  std::string _nameofdatafile;

  //! files containing the data and the timeseries
  std::ofstream _ofdata;

 public:
  VtkWriter( ) { std::cout << "Constructor of VtkWriter" << endl; }

  void openfiles(const std::string &s) {
    _nameoffile = s;
    std::string tmp = s + ".vtu";
    std::cout << tmp << " ";
    _ofdata.open(tmp.c_str( ), std::ios_base::out);
    _nameofdatafile = "";

    for (int i = 0; i < tmp.length( ); ++i) {
      if (tmp.at(i) == '\\') {
        _nameofdatafile.append(1, '\\');
      }

      _nameofdatafile.append(1, tmp.at(i));
    }
  }

  void addmesh(const Fem2D::Mesh *mesh) {
    const Fem2D::Mesh &Th(*mesh);

    _vecmesh.push_back(mesh);
    _ofdata.flags(std::ios_base::scientific);
    _ofdata.precision(15);

    _ofdata << "<?xml version=\"1.0\"?>" << std::endl;
    _ofdata << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">";
    _ofdata << std::endl;
    _ofdata << "<UnstructuredGrid>";
    _ofdata << std::endl;
    _ofdata << "<Piece NumberOfPoints=\"" << Th.nv << "\" NumberOfCells=\"" << Th.nt << "\">";
    _ofdata << std::endl;
    _ofdata << "<Points>" << std::endl;
    _ofdata
      << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">";
    _ofdata << std::endl;

    for (int k = 0; k < Th.nv; ++k) {
      _ofdata << Th(k).x << " " << Th(k).y << " " << 0.0 << std::endl;
    }

    _ofdata << "</DataArray>" << std::endl;
    _ofdata << "</Points>" << std::endl;
    _ofdata << "<Cells>" << std::endl;
    _ofdata << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" "
               "format=\"ascii\">";
    _ofdata << std::endl;

    for (int i = 0; i < Th.nt; ++i) {
      for (int j = 0; j < 3; j++) {
        _ofdata << Th(i, j) << " ";
      }
    }

    _ofdata << std::endl;
    _ofdata << "</DataArray>" << std::endl;
    _ofdata
      << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">";
    _ofdata << std::endl;

    for (int i = 0; i < Th.nt; ++i) {
      _ofdata << 3 + 3 * (i) << " ";
    }

    _ofdata << std::endl;
    _ofdata << "</DataArray>" << std::endl;
    _ofdata
      << "<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">";
    _ofdata << std::endl;

    for (int i = 0; i < Th.nt; ++i) {
      _ofdata << 5 << " ";
    }

    _ofdata << std::endl;
    _ofdata << "</DataArray>" << std::endl;
    _ofdata << "</Cells>" << std::endl;
    _ofdata << "<PointData >" << endl;
  }

  double checkprecision(double val) {
    double tmp;

    if (val >= 0.)
      tmp = max(0., val);
    else
      tmp = min(0., val);

    return tmp;
  }

  /*!Add a field*/
  void addscalar(const string &nameoffield, const Fem2D::Mesh *mesh, const KN< double > &val) {
    _ofdata.flags(std::ios_base::scientific);
    _ofdata.precision(15);

    _ofdata << "<DataArray type=\"Float32\" Name=\"";
    _ofdata << nameoffield << "\" NumberOfComponents=\"1\" format=\"ascii\">";
    _ofdata << std::endl;

    for (int i = 0; i < val.size( ); ++i) {
      _ofdata << checkprecision(val[i]) << std::endl;
    }

    _ofdata << "</DataArray>" << std::endl;

    _ofdata.flush( );
  }

  /*!Add a field*/
  void addvector(const string &nameoffield, const Fem2D::Mesh *mesh, const KN< double > &val,
                 const KN< double > &val2) {
    _ofdata.flags(std::ios_base::scientific);
    _ofdata.precision(15);

    _ofdata << "<DataArray type=\"Float32\" Name=\"";
    _ofdata << nameoffield << "\" NumberOfComponents=\"3\" format=\"ascii\">";
    _ofdata << std::endl;

    for (int i = 0; i < val.size( ); ++i) {
      _ofdata << checkprecision(val[i]) << " " << checkprecision(val2[i]) << " " << 0.0
              << std::endl;
    }

    _ofdata << "</DataArray>" << std::endl;
    _ofdata.flush( );
  }

  /*!Get the mesh associated with the series nameofts*/
  const Fem2D::Mesh *getmeshts(const string &nameofts) { return _vecmesh[0]; }

  void init( ) { new (this) VtkWriter( ); }

  void destroy( ) {
    if (_ofdata.is_open( )) {
      _ofdata << "</PointData>" << endl;
      _ofdata << "<CellData>" << endl;
      _ofdata << "</CellData>" << endl;
      _ofdata << "</Piece>" << endl;
      _ofdata << "</UnstructuredGrid>" << endl;
      _ofdata << "</VTKFile>" << endl;
      _ofdata.close( );
    }
  }
};    // End of class

class Vtkwritesol_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression edx;
  Expression ename;    //!< name of time series or field
  Expression et;       //!< time
  long what;           // 1 scalar, 2 vector, 3 symtensor
  long nbfloat;        // 1 scalar, n vector (3D), n symtensor(3D)
  Expression evct, evct2;

 public:
  Vtkwritesol_Op(const basicAC_F0 &args) : what(0), nbfloat(0) {
    evct = 0;
    evct2 = 0;
    int nbofsol;
    int ddim = 2;
    // There's no named parameter
    args.SetNameParam( );
    if (args.size( ) != 3) {
      CompileError("Vtkwritesol accepts only 4 parameters");
    }

    if (BCastTo< VtkWriter * >(args[0])) {
      edx = CastTo< VtkWriter * >(args[0]);
    }

    if (BCastTo< string * >(args[1])) {
      ename = CastTo< string * >(args[1]);
    }

    if (args[2].left( ) == atype< double >( )) {
      what = 1;
      nbfloat = 1;
      evct = to< double >(args[2]);
    } else if (args[2].left( ) == atype< double * >( )) {
      what = 1;
      nbfloat = 1;
      evct = to< double >(args[2]);
    } else if (BCastTo< pfer >(args[2])) {
      what = 1;
      nbfloat = 1;
      evct = to< double >(args[2]);
    } else if (args[2].left( ) == atype< E_Array >( )) {
      std::cout << "Until now only scalar solution" << std::endl;
      int i = 2;
      const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));

      if (a0->size( ) == ddim) {
        // vector solution
        what = 2;
        nbfloat = a0->size( );
        evct = to< double >((*a0)[0]);
        evct2 = to< double >((*a0)[1]);
      }

      cout << "Passed Until now only scalar solution" << std::endl;
    } else {
      CompileError("savesol in 2D: Sorry no way to save this kind of data");
    }
  }

  // all type
  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< VtkWriter * >( ), atype< string * >( ), true);
  }

  static E_F0 *f(const basicAC_F0 &args) { return new Vtkwritesol_Op(args); }

  AnyType operator( )(Stack stack) const;
};    // end of class

AnyType Vtkwritesol_Op::operator( )(Stack stack) const {
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  VtkWriter &dx = *(GetAny< VtkWriter * >((*edx)(stack)));
  string &name = *(GetAny< string * >((*ename)(stack)));
  const Mesh &Th = *(dx.getmeshts(name));
  int nt = Th.nt;
  int nv = Th.nv;
  int nbsol = nv;
  long longdefault = 0;

  KN< double > valsol(nbsol);
  valsol = 0.;
  KN< int > takemesh(nbsol);
  takemesh = 0;
  MeshPoint *mp3(MeshPointStack(stack));

  for (int it = 0; it < nt; it++) {
    for (int iv = 0; iv < 3; iv++) {
      int i = Th(it, iv);
      mp3->setP(&Th, it, iv);
      valsol[i] = valsol[i] + GetAny< double >((*evct)(stack));
      ++takemesh[i];
    }
  }

  for (int i = 0; i < nbsol; i++) {
    valsol[i] /= takemesh[i];
  }

  // Writes valsol on the file file
  if (what == 1) {
    dx.addscalar(name, &Th, valsol);
  }

  if (what == 2) {
    KN< double > valsol2(nbsol);
    valsol2 = 0.;
    KN< int > takemesh(nbsol);
    takemesh = 0;
    MeshPoint *mp3(MeshPointStack(stack));

    for (int it = 0; it < nt; it++) {
      for (int iv = 0; iv < 3; iv++) {
        int i = Th(it, iv);
        mp3->setP(&Th, it, iv);
        valsol2[i] = valsol2[i] + GetAny< double >((*evct2)(stack));
        ++takemesh[i];
      }
    }

    for (int i = 0; i < nbsol; i++) {
      valsol2[i] /= takemesh[i];
    }

    // Writes valsol on the file file
    dx.addvector(name, &Th, valsol, valsol2);
  }

  return longdefault;
}

// le vrai constructeur est la
VtkWriter *init_VtkWriter(VtkWriter *const &a, string *const &s) {
  std::cout << "start init_VtkWriter" << std::endl;

  a->init( );
  a->openfiles(*s);
  std::cout << "end init_VtkWriter" << std::endl;
  return a;
}

void *call_addmesh(VtkWriter *const &mt, const Fem2D::Mesh *const &pTh) {
  mt->addmesh(pTh);
  return NULL;
}

// Add the function name to the freefem++ table
static void Load_Init( ) {
  Dcl_Type< VtkWriter * >(InitP< VtkWriter >, Destroy< VtkWriter >);
  // declare deux nouveau type pour freefem++  un pointeur et

  zzzfff->Add("VtkWriter", atype< VtkWriter * >( ));    // ajoute le type myType a freefem++
  // constructeur  d'un type myType  dans freefem
  TheOperators->Add("<-", new OneOperator2_< VtkWriter *, VtkWriter *, string * >(&init_VtkWriter));
  Global.Add("Vtkaddmesh", "(",
             new OneOperator2_< void *, VtkWriter *, const Fem2D::Mesh * >(call_addmesh));
  Global.Add("Vtkaddscalar", "(", new OneOperatorCode< Vtkwritesol_Op >);
}

LOADFUNC(Load_Init)
