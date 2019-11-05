// TODO: remove this block as soon as autoconf is removed from FreeFem++
#ifndef CMAKE
#include "config.h"
#endif
extern long verbosity;
#ifdef HAVE_HDF5
#include "write_hdf5.hpp"
using std::max;
using std::min;


WriteHdf5::WriteHdf5(const char *ffname, Int4 nbelem, Int4 nbvertex) : hdf5_filename(ffname),  nbofelem(nbelem), nbofvertex(nbvertex)
{
  //constructeur
  //creation du fichier hdf5
  file_id = H5Fcreate(hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status=0;
}

WriteHdf5::~WriteHdf5()
{
}

void WriteHdf5::WriteHdf5MeshFile2D(float coordinates[][2], int connec[][3])
{
  //ecriture du maillage 2D au format hdf5 pour paraview

  hid_t group_id_connec;
  hid_t group_id_coord;

  hid_t dataset_id_coord, dataspace_id_coord;
  hid_t aid_max_x, aid_min_x, aid_max_y, aid_min_y;
  hid_t attr_max_x, attr_min_x, attr_max_y, attr_min_y;
  hid_t dataset_id_elem2node, dataspace_id_elem2node;
  hid_t dataset_id_nelem, dataspace_id_nelem;
  hid_t dataset_id_nnode, dataspace_id_nnode;
  hid_t dataset_id_elemtype, dataspace_id_elemtype;
  hid_t dataset_id_meshtype, dataspace_id_meshtype;
  hid_t type_id,type_mesh_id;

  hsize_t dims_coord[2];
  hsize_t dims_x_max[1];
  hsize_t dims_x_min[1];
  hsize_t dims_y_max[1];
  hsize_t dims_y_min[1];
  hsize_t dims_elem2node[2];
  hsize_t dims_nelem[1];
  hsize_t dims_nnode[1];

  dims_coord[0] = nbofvertex;
  dims_coord[1] = 2;
  dims_x_max[0] = 1;
  dims_x_min[0] = 1;
  dims_y_max[0] = 1;
  dims_y_min[0] = 1;
  dims_elem2node[0] = nbofelem;
  dims_elem2node[1] = 3;
  dims_nelem[0] = 1;
  dims_nnode[0] = 1;

  float x_max = coordinates[0][0];
  float x_min = coordinates[0][0];
  float y_max = coordinates[0][1];
  float y_min = coordinates[0][1];

  for (int i=0;i<nbofvertex;i++)
    {
      x_max= max(x_max,coordinates[i][0]);
      y_max= max(y_max,coordinates[i][1]);
      x_min= min(x_min,coordinates[i][0]);
      y_min= min(y_min,coordinates[i][1]);
    }

  //ecriture des coordonnees 2D X et Y des noeuds
  group_id_coord = H5Gcreate2(file_id, "/Coordinates", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dataspace_id_coord = H5Screate_simple(2, dims_coord, NULL);
  dataset_id_coord = H5Dcreate2(file_id, "/Coordinates/XY", H5T_IEEE_F32LE, dataspace_id_coord, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id_coord, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, coordinates);

  //ecriture des valeurs max des coordonnees
  aid_max_x =  H5Screate_simple(1, dims_x_max, NULL);
  attr_max_x = H5Acreate2(dataset_id_coord, "X_MAX",  H5T_IEEE_F32LE, aid_max_x, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attr_max_x, H5T_NATIVE_FLOAT, &x_max);
  status = H5Aclose(attr_max_x);
  status = H5Sclose(aid_max_x);

  aid_min_x =  H5Screate_simple(1, dims_x_min, NULL);
  attr_min_x = H5Acreate2(dataset_id_coord, "X_MIN",  H5T_IEEE_F32LE, aid_min_x, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attr_min_x, H5T_NATIVE_FLOAT, &x_min);
  status = H5Aclose(attr_min_x);
  status = H5Sclose(aid_min_x);

  aid_max_y =  H5Screate_simple(1, dims_y_max, NULL);
  attr_max_y = H5Acreate2(dataset_id_coord, "Y_MAX",  H5T_IEEE_F32LE, aid_max_y, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attr_max_y, H5T_NATIVE_FLOAT, &y_max);
  status = H5Aclose(attr_max_y);
  status = H5Sclose(aid_max_y);

  aid_min_y =  H5Screate_simple(1, dims_y_min, NULL);
  attr_min_y = H5Acreate2(dataset_id_coord, "Y_MIN",  H5T_IEEE_F32LE, aid_min_y, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attr_min_y, H5T_NATIVE_FLOAT, &y_min);
  status = H5Aclose(attr_min_y);
  status = H5Sclose(aid_min_y);

  status = H5Dclose(dataset_id_coord);
  status = H5Sclose(dataspace_id_coord);
  status = H5Gclose(group_id_coord);

  //ecriture du tableau de connectivite (3 numeros de noeud definissent 1 triangle)
  group_id_connec = H5Gcreate2(file_id, "/Connectivity", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  dataspace_id_elem2node = H5Screate_simple(2, dims_elem2node, NULL);
  dataset_id_elem2node = H5Dcreate2(file_id, "/Connectivity/ELEM2NODE", H5T_STD_I32LE,
				   dataspace_id_elem2node, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id_elem2node, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, connec);
  status = H5Dclose(dataset_id_elem2node);
  status = H5Sclose(dataspace_id_elem2node);

  //ecriture du type d'elements = triangle, du nombre de triangles et du nombre de noeuds du maillage
  char elemtype[9];
  strcpy(elemtype,"Triangle");
  type_id = H5Tcopy(H5T_C_S1);
  status  = H5Tset_size(type_id, 8);
  dataspace_id_elemtype = H5Screate(H5S_SCALAR);
  dataset_id_elemtype=H5Dcreate2(file_id, "/Connectivity/ELEMTYPE", type_id,
				dataspace_id_elemtype, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id_elemtype, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, elemtype);
  status = H5Dclose(dataset_id_elemtype);
  status = H5Sclose(dataspace_id_elemtype);

  dataspace_id_nelem = H5Screate_simple(1, dims_nelem, NULL);
  dataset_id_nelem = H5Dcreate2(file_id, "/Connectivity/NELEM", H5T_STD_I32LE,
			       dataspace_id_nelem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id_nelem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nbofelem);
  status = H5Dclose(dataset_id_nelem);
  status = H5Sclose(dataspace_id_nelem);

  dataspace_id_nnode = H5Screate_simple(1, dims_nnode, NULL);
  dataset_id_nnode = H5Dcreate2(file_id, "/Connectivity/NNODE", H5T_STD_I32LE,
			       dataspace_id_nnode, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id_nnode, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nbofvertex);
  status = H5Dclose(dataset_id_nnode);
  status = H5Sclose(dataspace_id_nnode);

  //ecriture du type de connectivite
  char meshtype[8];
  strcpy(meshtype,"UNIFORM");
  type_mesh_id = H5Tcopy(H5T_C_S1);
  status  = H5Tset_size(type_mesh_id, 7);
  dataspace_id_meshtype = H5Screate(H5S_SCALAR);
  dataset_id_meshtype=H5Dcreate2(file_id, "/Connectivity/TYPE", type_mesh_id,
				dataspace_id_meshtype, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id_meshtype, type_mesh_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, meshtype);
  status = H5Dclose(dataset_id_meshtype);
  status = H5Sclose(dataspace_id_meshtype);
  status = H5Tclose(type_mesh_id);

  status = H5Gclose(group_id_connec);

  status = H5Fclose(file_id);
  if (status != 0) cout << "HDF5 error" << endl;
  cout<<"save hdf5 file mesh : "<< hdf5_filename << endl;
}

// void WriteHdf5::WriteHdf5MeshFile3D(float coordinates[][3], int connec[][4])
// {
//   //ecriture du maillage 3D au format hdf5 pour paraview
//
//   hid_t group_id_coord;
//   hid_t group_id_connec;
//   hid_t dataset_id_coord, dataspace_id_coord;
//   hid_t aid_max_x, aid_min_x, aid_max_y, aid_min_y, aid_max_z, aid_min_z;
//   hid_t attr_max_x, attr_min_x, attr_max_y, attr_min_y,  attr_max_z, attr_min_z;
//   hid_t dataset_id_elem2node, dataspace_id_elem2node;
//   hid_t dataset_id_elemtype, dataspace_id_elemtype;
//   hid_t dataset_id_nelem, dataspace_id_nelem;
//   hid_t dataset_id_nnode, dataspace_id_nnode;
//   hid_t dataset_id_meshtype, dataspace_id_meshtype;
//   hid_t type_id,type_mesh_id;
//
//   hsize_t dims_coord[2];
//   hsize_t dims_x_max[1];
//   hsize_t dims_x_min[1];
//   hsize_t dims_y_max[1];
//   hsize_t dims_y_min[1];
//   hsize_t dims_z_max[1];
//   hsize_t dims_z_min[1];
//   hsize_t dims_elem2node[2];
//   hsize_t dims_nelem[1];
//   hsize_t dims_nnode[1];
//
//   dims_coord[0] = nbofvertex;
//   dims_coord[1] = 3;
//   dims_x_max[0] = 1;
//   dims_x_min[0] = 1;
//   dims_y_max[0] = 1;
//   dims_y_min[0] = 1;
//   dims_z_max[0] = 1;
//   dims_z_min[0] = 1;
//   dims_elem2node[0] = nbofelem;
//   dims_elem2node[1] = 4;
//   dims_nelem[0] = 1;
//   dims_nnode[0] = 1;
//
//   float x_max = 0 ;
//   float x_min = 0 ;
//   float y_max = 0 ;
//   float y_min = 0 ;
//   float z_max = 0 ;
//   float z_min = 0 ;
//
//   for (int i=0;i<nbofvertex;i++)
//     {
//       x_max= max(x_max,coordinates[i][0]);
//       y_max= max(y_max,coordinates[i][1]);
//       z_max= max(z_max,coordinates[i][2]);
//       x_min= min(x_min,coordinates[i][0]);
//       y_min= min(y_min,coordinates[i][1]);
//       z_min= min(z_min,coordinates[i][2]);
//     }
//
//   //ecriture des coordonnees 3D X,Y,Z des noeuds
//   group_id_coord = H5Gcreate2(file_id, "/Coordinates", H5P_DEFAULT, H5P_DEFAULT,
// 			     H5P_DEFAULT);
//   dataspace_id_coord = H5Screate_simple(2, dims_coord, NULL);
//   dataset_id_coord = H5Dcreate2(file_id, "/Coordinates/XYZ", H5T_IEEE_F32LE,
// 			       dataspace_id_coord, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Dwrite(dataset_id_coord, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
// 		    H5P_DEFAULT, coordinates);
//
//   //ecriture des valeurs max des coordonnees
//   aid_max_x =  H5Screate_simple(1, dims_x_max, NULL);
//   attr_max_x = H5Acreate2(dataset_id_coord, "X_MAX",  H5T_IEEE_F32LE, aid_max_x,
// 			 H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Awrite(attr_max_x, H5T_NATIVE_FLOAT, &x_max);
//   status = H5Aclose(attr_max_x);
//   status = H5Sclose(aid_max_x);
//
//   aid_min_x =  H5Screate_simple(1, dims_x_min, NULL);
//   attr_min_x = H5Acreate2(dataset_id_coord, "X_MIN",  H5T_IEEE_F32LE, aid_min_x,
// 			 H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Awrite(attr_min_x, H5T_NATIVE_FLOAT, &x_min);
//   status = H5Aclose(attr_min_x);
//   status = H5Sclose(aid_min_x);
//
//   aid_max_y =  H5Screate_simple(1, dims_y_max, NULL);
//   attr_max_y = H5Acreate2(dataset_id_coord, "Y_MAX",  H5T_IEEE_F32LE, aid_max_y,
// 			 H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Awrite(attr_max_y, H5T_NATIVE_FLOAT, &y_max);
//   status = H5Aclose(attr_max_y);
//   status = H5Sclose(aid_max_y);
//
//   aid_min_y =  H5Screate_simple(1, dims_y_min, NULL);
//   attr_min_y = H5Acreate2(dataset_id_coord, "Y_MIN",  H5T_IEEE_F32LE, aid_min_y,
// 			 H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Awrite(attr_min_y, H5T_NATIVE_FLOAT, &y_min);
//   status = H5Aclose(attr_min_y);
//   status = H5Sclose(aid_min_y);
//
//   aid_max_z =  H5Screate_simple(1, dims_z_max, NULL);
//   attr_max_z = H5Acreate2(dataset_id_coord, "Z_MAX",  H5T_IEEE_F32LE, aid_max_z,
// 			 H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Awrite(attr_max_z, H5T_NATIVE_FLOAT, &z_max);
//   status = H5Aclose(attr_max_z);
//   status = H5Sclose(aid_max_z);
//
//   aid_min_z =  H5Screate_simple(1, dims_z_min, NULL);
//   attr_min_z = H5Acreate2(dataset_id_coord, "Z_MIN",  H5T_IEEE_F32LE, aid_min_z,
// 			 H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Awrite(attr_min_z, H5T_NATIVE_FLOAT, &z_min);
//   status = H5Aclose(attr_min_z);
//   status = H5Sclose(aid_min_z);
//
//   status = H5Dclose(dataset_id_coord);
//   status = H5Sclose(dataspace_id_coord);
//   status = H5Gclose(group_id_coord);
//
//   //ecriture du tableau de connectivite (4 numeros de noeud definissent 1 tetraedre)
//   group_id_connec = H5Gcreate2(file_id, "/Connectivity", H5P_DEFAULT, H5P_DEFAULT,
// 			      H5P_DEFAULT);
//
//   dataspace_id_elem2node = H5Screate_simple(2, dims_elem2node, NULL);
//   dataset_id_elem2node = H5Dcreate2(file_id, "/Connectivity/ELEM2NODE", H5T_STD_I32LE,
// 				   dataspace_id_elem2node, H5P_DEFAULT, H5P_DEFAULT,
// 				   H5P_DEFAULT);
//   status = H5Dwrite(dataset_id_elem2node, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
// 		    H5P_DEFAULT, connec);
//   status = H5Dclose(dataset_id_elem2node);
//   status = H5Sclose(dataspace_id_elem2node);
//
//   //ecriture du type d'elements = tetraedre, du nombre de tetraedres et du nombre de noeuds du maillage
//   char elemtype[12];
//   strcpy(elemtype,"Tetrahedron");
//   type_id = H5Tcopy(H5T_C_S1);
//   status  = H5Tset_size(type_id, 11);
//   dataspace_id_elemtype = H5Screate(H5S_SCALAR);
//   dataset_id_elemtype=H5Dcreate2(file_id, "/Connectivity/ELEMTYPE", type_id,
// 				dataspace_id_elemtype, H5P_DEFAULT, H5P_DEFAULT,
// 				H5P_DEFAULT);
//   status = H5Dwrite(dataset_id_elemtype, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
// 		    elemtype);
//   status = H5Dclose(dataset_id_elemtype);
//   status = H5Sclose(dataspace_id_elemtype);
//
//   dataspace_id_nelem = H5Screate_simple(1, dims_nelem, NULL);
//   dataset_id_nelem = H5Dcreate2(file_id, "/Connectivity/NELEM", H5T_STD_I32LE,
// 			       dataspace_id_nelem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Dwrite(dataset_id_nelem, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
// 		    &nbofelem);
//   status = H5Dclose(dataset_id_nelem);
//   status = H5Sclose(dataspace_id_nelem);
//
//   dataspace_id_nnode = H5Screate_simple(1, dims_nnode, NULL);
//   dataset_id_nnode = H5Dcreate2(file_id, "/Connectivity/NNODE", H5T_STD_I32LE,
// 			       dataspace_id_nnode, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//   status = H5Dwrite(dataset_id_nnode, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
// 		    &nbofvertex);
//   status = H5Dclose(dataset_id_nnode);
//   status = H5Sclose(dataspace_id_nnode);
//
//   //ecriture du type de connectivite
//   char meshtype[8];
//   strcpy(meshtype,"UNIFORM");
//   type_mesh_id = H5Tcopy(H5T_C_S1);
//   status  = H5Tset_size(type_mesh_id, 7);
//   dataspace_id_meshtype = H5Screate(H5S_SCALAR);
//   dataset_id_meshtype=H5Dcreate2(file_id, "/Connectivity/TYPE", type_mesh_id,
// 				dataspace_id_meshtype, H5P_DEFAULT, H5P_DEFAULT,
// 				H5P_DEFAULT);
//   status = H5Dwrite(dataset_id_meshtype, type_mesh_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
// 		    meshtype);
//   status = H5Dclose(dataset_id_meshtype);
//   status = H5Sclose(dataspace_id_meshtype);
//   status = H5Tclose(type_mesh_id);
//
//   status = H5Gclose(group_id_connec);
//
//   status = H5Fclose(file_id);
//
//   if (status != 0) cout << "HDF5 error" << endl;
//   cout <<"save hdf5 file mesh : "<< hdf5_filename << endl;
// }

void WriteHdf5::WriteHdf5SolFile2DInit()
{
  //creation du groupe /Data contenant toutes les donnees solution
  group_id_data = H5Gcreate2(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void WriteHdf5::WriteHdf5SolFile2DAddField(string * fieldname, int result_order, int trans_dim, int what_type, float *field)
{
  //ajout des champs solution dans le groupe /Data

  string str_data="/Data/";
  size_t size_data = str_data.size() + 1;
  string str_float="FLOAT_";
  size_t size_str_float=str_float.size() + 1;
  string str_underscore="_";
  size_t size_str_underscore=str_underscore.size() + 1;
  const int ldata_type=100;
  char data_type[ldata_type];// to be sure
  hid_t dataset_id_data, dataspace_id_data;
  hid_t type_id;
  hid_t attr_type;
  hid_t aid_type;
  hsize_t dims_data[2];// coorect FH.
  vector<string> type_char(3);
  vector<string> res_char(2);

  type_char[0]= "Scalar";
  type_char[1]= "Vector";
  type_char[2]= "Vector";
  res_char[0]= "Cell";
  res_char[1]= "Node";
  
  size_t size_datafieldname = fieldname->size() + 1;
  char * char_datafieldname = new char[size_datafieldname];
  strncpy(char_datafieldname, fieldname->c_str(), size_datafieldname);
  char * char_datafieldname_tot=new char[size_data+size_datafieldname];
  strncpy(char_datafieldname_tot, str_data.c_str(), size_data);
  strncat(char_datafieldname_tot,char_datafieldname,size_data+size_datafieldname);
  strncpy(data_type, str_float.c_str(), ldata_type);

  if(result_order==0)
    {
      dims_data[0] = nbofelem;
    }
  else
    {
      dims_data[0] = nbofvertex;
    }
  dims_data[1] = trans_dim;

  //ajout des champs fonction de la dimension de l'element (triangle ou noeud) et du type
  //de donnees (vecteur, ...) affectees a l'element
  strncat(data_type,  res_char[result_order].c_str(), ldata_type);
  strncat(data_type, str_underscore.c_str(), ldata_type);
  strncat(data_type, type_char[what_type].c_str(),ldata_type);
  assert(strnlen(data_type,ldata_type)<ldata_type-2); // verif data_type is long enough
  dataspace_id_data = H5Screate_simple(2, dims_data, NULL);
  dataset_id_data = H5Dcreate2(file_id, char_datafieldname_tot, H5T_IEEE_F32LE,
  			      dataspace_id_data, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id_data, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, field);
  type_id = H5Tcopy(H5T_C_S1);
  status  = H5Tset_size(type_id, 17);
  aid_type =  H5Screate(H5S_SCALAR);
  attr_type = H5Acreate2(dataset_id_data, "TYPE",  type_id, aid_type, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attr_type, type_id, data_type);
  status = H5Aclose(attr_type);
  status = H5Sclose(aid_type);
  status = H5Dclose(dataset_id_data);
  status = H5Sclose(dataspace_id_data);
  delete [] char_datafieldname;
  delete [] char_datafieldname_tot;
}

void WriteHdf5::WriteHdf5SolFile2DFinalize()
{
  //fermeture du groupe /Data et du fichier hdf5
  status = H5Gclose(group_id_data);
  status = H5Fclose(file_id);
  if (status != 0) cout << "HDF5 error" << endl;
  if(verbosity>2) cout << "save hdf5 file solution : " << hdf5_filename <<endl;
}

void WriteHdf5::WriteHdf5SolFile3DInit()
{
  //creation du groupe /Data contenant toutes les donnees solution
  group_id_data = H5Gcreate2(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void WriteHdf5::WriteHdf5SolFile3DAddField(string * fieldname, int result_order, int trans_dim, int what_type, float *field)
{
  //ajout des champs solution dans le groupe /Data

  string str_data="/Data/";
  size_t size_data = str_data.size() + 1;
  string str_float="FLOAT_";
  size_t size_str_float=str_float.size() + 1;
  string str_underscore="_";
  size_t size_str_underscore=str_underscore.size() + 1;
  const int ldata_type=100;
  char data_type[ldata_type];
  hid_t dataset_id_data, dataspace_id_data;
  hid_t type_id;
  hid_t attr_type;
  hid_t aid_type;
  hsize_t dims_data[2];// correct FH. 1->2
  vector<string> type_char(3);
  vector<string> res_char(2);

  type_char[0]= "Scalar";
  type_char[1]= "Vector";
  type_char[2]= "Tensor6";
  res_char[0]= "Cell";
  res_char[1]= "Node";

  size_t size_datafieldname = fieldname->size() + 1;
  char * char_datafieldname = new char[size_datafieldname];
  strncpy(char_datafieldname, fieldname->c_str(), size_datafieldname);
  char * char_datafieldname_tot=new char[size_data+size_datafieldname];
  strncpy(char_datafieldname_tot, str_data.c_str(), size_data);
  strncat(char_datafieldname_tot,char_datafieldname,size_data+size_datafieldname);
  strncpy(data_type, str_float.c_str(), ldata_type);

  if(result_order==0)
    {
      dims_data[0] = nbofelem;
    }
  else
    {
      dims_data[0] = nbofvertex;
    }
  dims_data[1] = trans_dim;

  //ajout des champs fonction de la dimension de l'element (tetraedre ou noeud) et du type
  //de donnees (tenseur de 6 composantes, ...) affectees a l'element
  strncat(data_type,  res_char[result_order].c_str(), ldata_type);
  strncat(data_type, str_underscore.c_str(), ldata_type);
  strncat(data_type, type_char[what_type].c_str(),ldata_type);
  assert(strnlen(data_type,ldata_type)< ldata_type-2);
  dataspace_id_data = H5Screate_simple(2, dims_data, NULL);
  dataset_id_data = H5Dcreate2(file_id, char_datafieldname_tot, H5T_IEEE_F32LE,
  			      dataspace_id_data, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id_data, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, field);
  type_id = H5Tcopy(H5T_C_S1);
  status  = H5Tset_size(type_id, 18);
  aid_type =  H5Screate(H5S_SCALAR);
  attr_type = H5Acreate2(dataset_id_data, "TYPE",  type_id, aid_type, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attr_type, type_id, data_type);
  status = H5Aclose(attr_type);
  status = H5Sclose(aid_type);
  status = H5Dclose(dataset_id_data);
  status = H5Sclose(dataspace_id_data);
  delete [] char_datafieldname;
  delete [] char_datafieldname_tot;
}

void WriteHdf5::WriteHdf5SolFile3DFinalize()
{
  //fermeture du groupe /Data et du fichier hdf5
  status = H5Gclose(group_id_data);
  status = H5Fclose(file_id);
  if (status != 0) cout << "HDF5 error" << endl;
  if(verbosity>2) cout << "save hdf5 file solution : " << hdf5_filename <<endl;
}
#endif
