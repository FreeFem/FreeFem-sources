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

#ifndef write_xdmf_h_
#define write_xdmf_h_

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using std::string;
using std::vector;
typedef long Int4;

class WriteXdmf
{
	private:
		const char *WXffname;
		ofstream xdmf_file;
		char *Elemtype;
		Int4 nbofelem;
		Int4 nbofvertex;
		Int4 nbvperelem;
		char *xdmf_filename;
		int dimension;

	public:

		WriteXdmf (const char *ffname, Int4 nbelem, Int4 nbvertex);
		virtual ~WriteXdmf ();
		void WriteXdmfMeshFile2D ();
		void WriteXdmfMeshFile3D ();
		void WriteXdmfSolFile2DInit ();
		void WriteXdmfSolFile2DAddField (string *fieldname, int data_type, int result_order, int trans_dim);
		void WriteXdmfSolFile2DFinalize ();
		void WriteXdmfSolFile3DInit ();
		void WriteXdmfSolFile3DAddField (string *fieldname, int data_type, int result_order, int trans_dim);
		void WriteXdmfSolFile3DFinalize ();
};

#endif

