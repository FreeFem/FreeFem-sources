/*
 * This file is part of FreeFem++.
 *
 * FreeFem++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FreeFem++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef write_hdf5_h_
#define write_hdf5_h_

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
#endif	// H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;
using std::string;
using std::vector;
typedef long Int4;

class WriteHdf5
{
	private:
		ofstream hdf5_file;
		const char *hdf5_filename;
		Int4 nbofelem;
		Int4 nbofvertex;
		Int4 nbvperelem;
		hid_t file_id;
		hid_t group_id_data;
		herr_t status;
		int dimension;

	public:

		WriteHdf5 (const char *ffname, Int4 nbelem, Int4 nbvertex);
		virtual ~WriteHdf5 ();
		void WriteHdf5MeshFile2D (float coordinates [][2], int connec [][3]);
		void WriteHdf5MeshFile3D (float coordinates [][3], int connec [][4]);
		void WriteHdf5SolFile2DInit ();
		void WriteHdf5SolFile2DAddField (string *fieldname, int result_order, int trans_dim, int what_type, float *field);
		void WriteHdf5SolFile2DFinalize ();
		void WriteHdf5SolFile3DInit ();
		void WriteHdf5SolFile3DAddField (string *fieldname, int result_order, int trans_dim, int what_type, float *field);
		void WriteHdf5SolFile3DFinalize ();
};

#endif

