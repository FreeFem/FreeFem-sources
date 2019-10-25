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
// AUTHORS : ...
// E-MAIL  : ...

// Example C++ function "myfunction", dynamically loaded into "load.edp"

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "ff++.hpp"

template<class Mesh>
class MatrixEdgeP1:  public E_F0 {
	public:
		typedef Matrice_Creuse<R> *Result;
		Expression emat, expTh;
		MatrixEdgeP1 (const basicAC_F0 &args) {
			args.SetNameParam();
			emat = args[0];
			expTh = to<const Mesh *>(args[1]);
			
		}

		MatrixEdgeP1 ()
		{}

		static ArrayOfaType typeargs () {return ArrayOfaType(atype<Matrice_Creuse<R> *>(), atype<const Mesh *>());}

		static E_F0*f (const basicAC_F0 &args) {return new MatrixEdgeP1(args);}

		AnyType operator () (Stack s) const;
};


template<class Mesh>
AnyType MatrixEdgeP1<Mesh>::operator () (Stack stack) const {
	Matrice_Creuse<R> *sparce_mat = GetAny<Matrice_Creuse<R> *>((*emat)(stack));
	MatriceMorse<R> *amorse = 0;
	MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
         typedef const Mesh *pmesh;
	const Mesh *pTh = GetAny<pmesh>((*expTh)(stack));
	ffassert(pTh);
	const Mesh &Th(*pTh);
	{
		

		// const edge ...
            HashTable<SortArray<int,2>,int> e(Th.nv+Th.nt,Th.nv);
            int ne=0;
            for(int k=0; k<Th.nt;++k)
                for(int i=0;i<3; ++i)
                {
                    SortArray<int,2> eki(Th(k,(i+1)%3),Th(k,(i+2)%3) );
                    if(!e.find(eki)) e.add(eki,ne++);
                }
            if(verbosity>4 && mpirank==0) cout << " ne = " << ne  << endl;
            MatriceMorse<R> * pAij= new MatriceMorse<R>(ne,Th.nv), &Aij = *pAij ;
            //ffassert(Th.ne==ne);
		for (int k = 0; k < ne; k++)
                {
                    
                    int i1 = e.t[k].k.v[0];
                    int i2 = e.t[k].k.v[1];
                    Aij(k,i1) =-1;
                    Aij(k,i2) =+1;
		}

		amorse = pAij;//new MatriceMorse<R>(Th.nv, Th.nv, Aij, false);
	}
	sparce_mat->Uh = UniqueffId();
	sparce_mat->Vh = UniqueffId();
	sparce_mat->A.master(amorse);
         sparce_mat->typemat = 0;//(amorse->n == amorse->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE);// none square matrice (morse)
	*mp = mps;

	if (verbosity > 3) {cout << "  End Build MatEdgeP1 : " << endl;}

	return sparce_mat;
}




static void Load_Init () {
	cout << " lood: init Mat Edge 1 " << endl;
	Global.Add("MatrixEdgeP1", "(", new OneOperatorCode<MatrixEdgeP1<Mesh> >());
        Global.Add("MatrixEdgeP1", "(", new OneOperatorCode<MatrixEdgeP1<Mesh3> >());
	//Global.Add("MatUpWind0", "(", new OneOperatorCode<MatrixUpWind3>());
}

LOADFUNC(Load_Init)
