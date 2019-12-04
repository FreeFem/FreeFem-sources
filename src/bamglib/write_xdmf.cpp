#include "write_xdmf.hpp"
#ifdef DEBUG
int debug_w = 1;
#else
int debug_w = 0;
#endif

// La denomination, au niveau des fichiers d'input,
// des fichiers de maillage ou solution sont en partie fixes :
// Le nom des fichiers de maillage doit se finir par .mesh.h5
// Le nom des fichiers solution doit se finir par .sol.h5
// L'interface ne prend en compte que les elements de type triangle et tetraedre.

WriteXdmf::WriteXdmf(const char *ffname, Int4 nbelem, Int4 nbvertex)
  : WXffname(ffname), nbofelem(nbelem), nbofvertex(nbvertex) {
  // constructeur
  // recuperation du nom de fichier passe dans l'input
  // pour creation du nom de fichier xdmf
  int lll = strlen(WXffname);
  string str_xdmf_filename = WXffname;
  string str_xdmf = "xmf";
  str_xdmf_filename.replace(lll - 2, 3, str_xdmf);
  size_t size_xdmf_filename = str_xdmf_filename.size( ) + 1;
  xdmf_filename = new char[size_xdmf_filename];
  strncpy(xdmf_filename, str_xdmf_filename.c_str( ), size_xdmf_filename);
}

WriteXdmf::~WriteXdmf( ) {
  delete[] xdmf_filename;
  delete[] Elemtype;
}

void WriteXdmf::WriteXdmfMeshFile2D( ) {
  // initialisation du type d'element : triangle
  Elemtype = new char[strlen("Triangle") + 1];
  strcpy(Elemtype, "Triangle");
  dimension = 2;
  nbvperelem = 3;

  // creation et ecriture du fichier de description
  // du maillage 2D au format xdfm pour paraview
  xdmf_file.open(xdmf_filename);
  xdmf_file << "<?xml version=\"1.0\" ?>\n";
  xdmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  xdmf_file << "<Xdmf Version=\"2.0\">\n";
  xdmf_file << "  <Domain>\n";
  // pour l'instant nom du maillage par defaut : full mesh
  xdmf_file << "    <Grid Name=\""
            << "full mesh"
            << "\" GridType=\"Uniform\">\n";
  xdmf_file << "      <Topology TopologyType=\"" << Elemtype << "\"  NumberOfElements=\""
            << nbofelem << "\">\n";
  // nb nodes per element, ici classe triangle donc =3
  xdmf_file << "	<DataItem Format=\"HDF\" NumberType=\"Int\" Dimensions=\"" << nbofelem
            << " " << nbvperelem << "\">\n";
  xdmf_file << "           " << WXffname << ":/Connectivity/ELEM2NODE\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Topology>\n";
  xdmf_file << "      <Geometry GeometryType=\"XY\">\n";
  xdmf_file << "	<DataItem NumberType=\"Float\" Precision=\"8\" Dimensions=\"" << nbofvertex
            << " " << dimension << "\" Format=\"HDF\">\n";
  xdmf_file << "           " << WXffname << ":/Coordinates/XY\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Geometry>\n";
  xdmf_file << "    </Grid>\n";
  xdmf_file << "  </Domain>\n";
  xdmf_file << "</Xdmf>\n";
  xdmf_file.close( );

  cout << "save xdmf file mesh : " << xdmf_filename << endl;
}

// void WriteXdmf::WriteXdmfMeshFile3D()
// {
//   //initialisation du type d'element : tetraedre
//   Elemtype= new char[strlen("Tetrahedron")+1];
//   strcpy(Elemtype,"Tetrahedron");
//   dimension=3;
//   nbvperelem=4;
//
//   //creation et ecriture du fichier de description
//   //du maillage 3D au format xdfm pour paraview
//   xdmf_file.open (xdmf_filename);
//   xdmf_file << "<?xml version=\"1.0\" ?>\n";
//   xdmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
//   xdmf_file << "<Xdmf Version=\"2.0\">\n";
//   xdmf_file << "  <Domain>\n";
//   //pour l'instant nom du maillage par defaut : full mesh
//   xdmf_file << "    <Grid Name=\"" << "full mesh" << "\" GridType=\"Uniform\">\n";
//   xdmf_file << "      <Topology TopologyType=\"" << Elemtype
// 	    << "\"  NumberOfElements=\"" << nbofelem << "\">\n";
//   // nb nodes per element, ici classe Tetrahedra donc = 4
//   xdmf_file << "	<DataItem Format=\"HDF\" NumberType=\"Int\" Dimensions=\""
// 	    << nbofelem << " " << nbvperelem << "\">\n";
//   xdmf_file << "           " << WXffname  << ":/Connectivity/ELEM2NODE\n";
//   xdmf_file << "        </DataItem>\n";
//   xdmf_file << "      </Topology>\n";
//   xdmf_file << "      <Geometry GeometryType=\"XYZ\">\n";
//   xdmf_file << "	<DataItem NumberType=\"Float\" Precision=\"8\" Dimensions=\""
// 	    << nbofvertex << " "<< dimension << "\" Format=\"HDF\">\n";
//   xdmf_file << "           " << WXffname << ":/Coordinates/XYZ\n";
//   xdmf_file << "        </DataItem>\n";
//   xdmf_file << "      </Geometry>\n";
//   xdmf_file << "    </Grid>\n";
//   xdmf_file << "  </Domain>\n";
//   xdmf_file << "</Xdmf>\n";
//   xdmf_file.close();
//
//   cout <<"save xdmf file mesh : "<< xdmf_filename << endl;
// }

void WriteXdmf::WriteXdmfSolFile2DInit( ) {
  // initialisation du type d'element : triangle
  Elemtype = new char[strlen("Triangle") + 1];
  strcpy(Elemtype, "Triangle");
  dimension = 2;
  nbvperelem = 3;
  string str_h5_mesh = "mesh.h5";
  string str_h5_mesh_filename = WXffname;
  str_h5_mesh_filename.replace(strlen(WXffname) - 6, 7, str_h5_mesh);

  // creation et debut d'ecriture du fichier de description
  // du fichier solution 2D au format xdfm pour paraview
  // avec le lien vers le fichier de maillage hdf5 correspondant
  xdmf_file.open(xdmf_filename);
  xdmf_file << "<?xml version=\"1.0\" ?>\n";
  xdmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  xdmf_file << "<Xdmf Version=\"2.0\">\n";
  xdmf_file << "  <Domain>\n";
  // pour l'instant nom du maillage par defaut : full mesh
  xdmf_file << "    <Grid Name=\""
            << "full mesh"
            << "\" GridType=\"Uniform\">\n";
  xdmf_file << "      <Topology TopologyType=\"" << Elemtype << "\"  NumberOfElements=\""
            << nbofelem << "\">\n";
  // nb nodes per element, ici classe triangle donc =3
  xdmf_file << "	<DataItem Format=\"HDF\" NumberType=\"Int\" Dimensions=\"" << nbofelem
            << " " << nbvperelem << "\">\n";
  xdmf_file << "           " << str_h5_mesh_filename.c_str( ) << ":/Connectivity/ELEM2NODE\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Topology>\n";
  xdmf_file << "      <Geometry GeometryType=\"XY\">\n";
  xdmf_file << "	<DataItem NumberType=\"Float\" Precision=\"8\" Dimensions=\"" << nbofvertex
            << " " << dimension << "\" Format=\"HDF\">\n";
  xdmf_file << "           " << str_h5_mesh_filename.c_str( ) << ":/Coordinates/XY\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Geometry>\n";
}

void WriteXdmf::WriteXdmfSolFile2DAddField(string *fieldname, int data_type, int result_order,
                                           int trans_dim) {
  vector< string > type_char(3);
  vector< string > res_char(2);

  type_char[0] = "Scalar";
  type_char[1] = "Vector";
  type_char[2] = "Vector";
  res_char[0] = "Cell";
  res_char[1] = "Node";

  // ecriture des descriptions du ou des champs solution 2D
  xdmf_file << "      <Attribute Name=\"" << fieldname->c_str( ) << "\" AttributeType=\""
            << type_char[data_type] << "\" Center=\"" << res_char[result_order] << "\">\n";
  if (result_order == 0) {
    xdmf_file << "        <DataItem Dimensions=\" " << nbofelem << " " << trans_dim
              << " \" NumberType=\""
              << "Float"
              << "\" Format=\"HDF\">\n";
  } else {
    xdmf_file << "        <DataItem Dimensions=\" " << nbofvertex << " " << trans_dim
              << " \" NumberType=\""
              << "Float"
              << "\" Format=\"HDF\">\n";
  }
  xdmf_file << "          " << WXffname << ":/Data/" << fieldname->c_str( ) << "\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Attribute>\n";
}

void WriteXdmf::WriteXdmfSolFile2DFinalize( ) {
  // fin d'ecriture et fermeture du fichier des descriptions solution 2D
  xdmf_file << "    </Grid>\n";
  xdmf_file << "  </Domain>\n";
  xdmf_file << "</Xdmf>\n";
  xdmf_file.close( );
  cout << "save xdmf file solution : " << xdmf_filename << endl;
}

void WriteXdmf::WriteXdmfSolFile3DInit( ) {
  // initialisation du type d'element : tetraedre
  Elemtype = new char[strlen("Tetrahedron") + 1];
  strcpy(Elemtype, "Tetrahedron");
  dimension = 3;
  nbvperelem = 4;
  string str_h5_mesh = "mesh.h5";
  string str_h5_mesh_filename = WXffname;
  str_h5_mesh_filename.replace(strlen(WXffname) - 6, 7, str_h5_mesh);

  // creation et debut d'ecriture du fichier de description
  // du fichier solution 3D au format xdfm pour paraview
  // avec le lien vers le fichier de maillage hdf5 correspondant
  xdmf_file.open(xdmf_filename);
  xdmf_file << "<?xml version=\"1.0\" ?>\n";
  xdmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  xdmf_file << "<Xdmf Version=\"2.0\">\n";
  xdmf_file << "  <Domain>\n";
  // pour l'instant nom du maillage par defaut : full mesh
  xdmf_file << "    <Grid Name=\""
            << "full mesh"
            << "\" GridType=\"Uniform\">\n";
  xdmf_file << "      <Topology TopologyType=\"" << Elemtype << "\"  NumberOfElements=\""
            << nbofelem << "\">\n";
  // nb nodes per element, ici classe tetrahedre donc =4
  xdmf_file << "	<DataItem Format=\"HDF\" NumberType=\"Int\" Dimensions=\"" << nbofelem
            << " " << nbvperelem << "\">\n";
  xdmf_file << "           " << str_h5_mesh_filename.c_str( ) << ":/Connectivity/ELEM2NODE\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Topology>\n";
  xdmf_file << "      <Geometry GeometryType=\"XYZ\">\n";
  xdmf_file << "	<DataItem NumberType=\"Float\" Precision=\"8\" Dimensions=\"" << nbofvertex
            << " " << dimension << "\" Format=\"HDF\">\n";
  xdmf_file << "           " << str_h5_mesh_filename.c_str( ) << ":/Coordinates/XYZ\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Geometry>\n";
}

void WriteXdmf::WriteXdmfSolFile3DAddField(string *fieldname, int data_type, int result_order,
                                           int trans_dim) {
  vector< string > type_char(3);
  vector< string > res_char(2);
  type_char[0] = "Scalar";
  type_char[1] = "Vector";
  type_char[2] = "Tensor6";
  res_char[0] = "Cell";
  res_char[1] = "Node";

  // ecriture des descriptions du ou des champs solution 3D
  xdmf_file << "      <Attribute Name=\"" << fieldname->c_str( ) << "\" AttributeType=\""
            << type_char[data_type] << "\" Center=\"" << res_char[result_order] << "\">\n";
  if (result_order == 0) {
    xdmf_file << "        <DataItem Dimensions=\" " << nbofelem << " " << trans_dim
              << " \" NumberType=\""
              << "Float"
              << "\" Format=\"HDF\">\n";
  } else {
    xdmf_file << "        <DataItem Dimensions=\" " << nbofvertex << " " << trans_dim
              << " \" NumberType=\""
              << "Float"
              << "\" Format=\"HDF\">\n";
  }
  xdmf_file << "          " << WXffname << ":/Data/" << fieldname->c_str( ) << "\n";
  xdmf_file << "        </DataItem>\n";
  xdmf_file << "      </Attribute>\n";
}

void WriteXdmf::WriteXdmfSolFile3DFinalize( ) {
  // fin d'ecriture et fermeture du fichier des descriptions solution 3D
  xdmf_file << "    </Grid>\n";
  xdmf_file << "  </Domain>\n";
  xdmf_file << "</Xdmf>\n";
  xdmf_file.close( );
  cout << "save xdmf file solution : " << xdmf_filename << endl;
}
