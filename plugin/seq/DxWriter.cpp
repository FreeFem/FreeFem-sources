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
// AUTHORS : Sala Lorenzo
// E-MAIL  : salallo80@gmail.com

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

#include "mode_open.hpp"
#include <iostream>
#include <cfloat>
#include <cmath>
#include <iterator>
using namespace std;
#include "ff++.hpp"
using namespace Fem2D;

/*!The class DxWriter permits to save in opendx format a "field"
 * (in the dx-language a "field" means the values of a function f(x,y,z) on a grid),
 * a time series (an ordered collection of "fields", so we have field0=f(x,y,z,t0),
 * field1=f(x,y,z,t1) and so on). DxWriter creates two files: one with extension .data where it puts
 * the position of the grid, the connessions, the values; and one with extension.dx where it puts
 * the time series. Now you can save only scalar fields. An example <code> load "DxWriter" mesh
 * Th=square(5,5); DxWriter ff("pippo"); Dxaddmesh(ff, Th); Dxaddtimeseries(ff, "Vx",Th); fespace
 * Vh(Th, P1); Vh vx=x*y; Dxaddsol2ts(ff,"Vx",1.0, vx); vx=0.5*x*y^2+0.2; Dxaddsol2ts(ff,"Vx",2.0,
 * vx); cout<<"Ciao";
 * </code>
 */

class DxWriter {
  struct tsinfo {
    int imesh;    //!< index of the mesh
    std::string name;
    std::vector< double > vecistant;
  };

 private:
  std::vector< const Fem2D::Mesh * > _vecmesh;
  std::vector< tsinfo > _vecofts;
  std::string _nameoffile;
  /*! This string contains the name of data file with \\ where there's a \ in the path*/
  std::string _nameofdatafile;
  //! files containing the data and the timeseries
  std::ofstream _ofdata, _ofts;

  /*!This function is called frequently, so if the main program crashes the files are good
   * and you need only write "end" at the end of data file: echo end>>nameoffile.data and then the
   * files are good
   */
  void save_header( ) {
    std::string s = _nameoffile;

    s.append(".dx");
    _ofts.open(s.c_str( ), std::ios_base::out);

    for (int i = 0; i < _vecofts.size( ); ++i) {
      _ofts << "object \"" << _vecofts[i].name << "\" class series" << std::endl;

      for (int j = 0; j < _vecofts[i].vecistant.size( ); ++j) {
        _ofts << "member " << j << " value file \"" << _nameofdatafile << "\",\""
              << _vecofts[i].name << "_" << j << "\" position " << _vecofts[i].vecistant[j]
              << std::endl;
      }

      _ofts << std::endl;
    }

    _ofts << "end" << std::endl;
    _ofts.close( );
  }

 public:
  DxWriter( ) { std::cout << "Constructor of DxWriter" << endl; }

  void openfiles(const std::string &s) {
    _nameoffile = s;
    std::string tmp = s + ".data";
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
    _ofdata << "object \"pos_" << _vecmesh.size( ) - 1
            << "\" class array type float rank 1 shape 2 items " << Th.nv << " data follows"
            << std::endl;

    for (int k = 0; k < Th.nv; ++k) {    // Scorre tutti i vertici
      _ofdata << Th(k).x << " " << Th(k).y << endl;
    }

    _ofdata << std::endl;
    _ofdata.flags(std::ios_base::fixed);
    _ofdata << "object \"conn_" << _vecmesh.size( ) - 1
            << "\" class array type int rank 1 shape 3 items " << Th.nt << " data follows " << endl;

    for (int i = 0; i < Th.nt; i++) {
      for (int j = 0; j < 3; j++) {
        _ofdata << Th(i, j) << " ";
      }

      _ofdata << endl;
    }

    _ofdata << "attribute \"element type\" string \"triangles\" " << std::endl;
    _ofdata << "attribute \"ref\" string \"positions\" " << std::endl << std::endl;
  }

  /*!Add a new time series, defined on the mesh*/
  void addtimeseries(const string &nameofts, const Fem2D::Mesh *mesh) {
    tsinfo ts;

    ts.name = nameofts;
    std::vector< const Fem2D::Mesh * >::const_iterator first = _vecmesh.begin( ),
                                                       last = _vecmesh.end( );

    if (std::find(first, last, mesh) == last) {
      addmesh(mesh);
      ts.imesh = _vecmesh.size( ) - 1;
    } else {
      ts.imesh = std::distance(first, std::find(first, last, mesh));
    }

    _vecofts.push_back(ts);
  }

  /*!Add an instant to a time series name*/
  void addistant2ts(const string &nameofts, const double t, const KN< double > &val) {
    int jj = -1;

    for (int i = 0; i < _vecofts.size( ); ++i) {
      if (_vecofts[i].name == nameofts) {
        jj = i;
      }
    }

    _vecofts[jj].vecistant.push_back(t);
    _ofdata.flags(std::ios_base::scientific);
    _ofdata.precision(15);
    _ofdata << "object \"" << nameofts << "_data_" << _vecofts[jj].vecistant.size( ) - 1
            << "\" class array type float rank 0 items " << val.size( ) << " data follows"
            << std::endl;

    for (int i = 0; i < val.size( ); ++i) {
      _ofdata << val[i] << std::endl;
    }

    _ofdata << "attribute \"dep\" string \"positions\"" << std::endl << std::endl;
    _ofdata << "object \"" << nameofts << "_" << _vecofts[jj].vecistant.size( ) - 1
            << "\" class field" << std::endl;
    _ofdata << "component \"positions\" value \"pos_" << _vecofts[jj].imesh << "\"" << std::endl;
    _ofdata << "component \"connections\" value \"conn_" << _vecofts[jj].imesh << "\"" << std::endl;
    _ofdata << "component \"data\" value \"" << nameofts << "_data_"
            << _vecofts[jj].vecistant.size( ) - 1 << "\"" << std::endl
            << std::endl;
    _ofdata.flush( );
    save_header( );
  }

  /*!Add a field*/
  void addfield(const string &nameoffield, const Fem2D::Mesh *mesh, const KN< double > &val) {
    std::vector< const Fem2D::Mesh * >::const_iterator first = _vecmesh.begin( ),
                                                       last = _vecmesh.end( );
    int im;

    if (std::find(first, last, mesh) == last) {
      addmesh(mesh);
      im = _vecmesh.size( ) - 1;
    } else {
      im = std::distance(first, std::find(first, last, mesh));
    }

    _ofdata.flags(std::ios_base::scientific);
    _ofdata.precision(15);
    _ofdata << "object \"" << nameoffield << "_data\" class array type float rank 0 items "
            << val.size( ) << " data follows" << std::endl;

    for (int i = 0; i < val.size( ); ++i) {
      _ofdata << val[i] << std::endl;
    }

    _ofdata << "attribute \"dep\" string \"positions\"" << std::endl << std::endl;
    _ofdata << "object \"" << nameoffield << "\" class field" << std::endl;
    _ofdata << "component \"positions\" value \"pos_" << im << "\"" << std::endl;
    _ofdata << "component \"connections\" value \"conn_" << im << "\"" << std::endl;
    _ofdata << "component \"data\" value \"" << nameoffield << "_data\"" << std::endl << std::endl;
    _ofdata.flush( );
  }

  /*!Get the mesh associated with the series nameofts*/
  const Fem2D::Mesh *getmeshts(const string &nameofts) {
    for (int i = 0; i < _vecofts.size( ); ++i) {
      if (_vecofts[i].name == nameofts) {
        return _vecmesh[_vecofts[i].imesh];
      }
    }

    return NULL;
  }

  void init( ) { new (this) DxWriter( ); }

  void destroy( ) {
    if (_ofdata.is_open( )) {
      _ofdata << std::endl << "end" << std::endl;
      _ofdata.close( );
    }
  }
};

class Dxwritesol_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression edx;
  Expression ename;    //!< name of time series or field
  Expression et;       //!< time
  long what;           // 1 scalar, 2 vector, 3 symtensor
  long nbfloat;        // 1 scalar, n vector (3D), n symtensor(3D)
  Expression evct;

 public:
  Dxwritesol_Op(const basicAC_F0 &args) : what(0), nbfloat(0) {
    evct = 0;
    // There's no named parameter
    args.SetNameParam( );
    if (args.size( ) != 4) {
      CompileError("Dxwritesol accepts only 4 parameters");
    }

    if (BCastTo< DxWriter * >(args[0])) {
      edx = CastTo< DxWriter * >(args[0]);
    }

    if (BCastTo< string * >(args[1])) {
      ename = CastTo< string * >(args[1]);
    }

    if (BCastTo< double >(args[2])) {
      et = CastTo< double >(args[2]);
    }

    if (args[3].left( ) == atype< double >( )) {
      what = 1;
      nbfloat = 1;
      evct = to< double >(args[3]);
    } else if (args[3].left( ) == atype< double * >( )) {
      what = 1;
      nbfloat = 1;
      evct = to< double >(args[3]);
    } else if (BCastTo< pfer >(args[3])) {
      what = 1;
      nbfloat = 1;
      evct = to< double >(args[3]);
    } else if (args[3].left( ) == atype< E_Array >( )) {
      CompileError("Until now only scalar solution");
    } else {
      CompileError("savesol in 2D: Sorry no way to save this kind of data");
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< DxWriter * >( ), atype< string * >( ), atype< double >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new Dxwritesol_Op(args); }

  AnyType operator( )(Stack stack) const;
};

AnyType Dxwritesol_Op::operator( )(Stack stack) const {
  MeshPoint *mp(MeshPointStack(stack));
  DxWriter &dx = *(GetAny< DxWriter * >((*edx)(stack)));
  string &name = *(GetAny< string * >((*ename)(stack)));
  double t = GetAny< double >((*et)(stack));
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
  dx.addistant2ts(name, t, valsol);

  return longdefault;
}

// le vrai constructeur est la
DxWriter *init_DxWriter(DxWriter *const &a, string *const &s) {
  std::cout << "start init_DxWriter" << std::endl;

  a->init( );
  a->openfiles(*s);
  std::cout << "end init_DxWriter" << std::endl;
  return a;
}

void *call_addmesh(DxWriter *const &mt, const Fem2D::Mesh *const &pTh) {
  mt->addmesh(pTh);
  return NULL;
}

void *call_addtimeseries(DxWriter *const &mt, string *const &name, const Fem2D::Mesh *const &pTh) {
  mt->addtimeseries(*name, pTh);
  return NULL;
}

// Add the function name to the freefem++ table
static void Load_Init( ) {
  Dcl_Type< DxWriter * >(
    InitP< DxWriter >,
    Destroy< DxWriter >);    // declare deux nouveau type pour freefem++  un pointeur et

  zzzfff->Add("DxWriter", atype< DxWriter * >( ));    // ajoute le type myType a freefem++
  // constructeur  d'un type myType  dans freefem
  TheOperators->Add("<-", new OneOperator2_< DxWriter *, DxWriter *, string * >(&init_DxWriter));

  Global.Add("Dxaddmesh", "(",
             new OneOperator2_< void *, DxWriter *, const Fem2D::Mesh * >(call_addmesh));
  Global.Add("Dxaddtimeseries", "(",
             new OneOperator3_< void *, DxWriter *, std::string *, const Fem2D::Mesh * >(
               call_addtimeseries));

  Global.Add("Dxaddsol2ts", "(", new OneOperatorCode< Dxwritesol_Op >);

  // atype< myType * >()->Add("(","",new OneOperator3_<myType_uv,myType *,double,double
  // >(set_myType_uv));
}

LOADFUNC(Load_Init)
