/****************************************************************************/
/* This file is part of FreeFem++.                                          */
/*                                                                          */
/* FreeFem++ is free software: you can redistribute it and/or modify        */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFem++ is distributed in the hope that it will be useful,             */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFem++. If not, see <http://www.gnu.org/licenses/>.        */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "ff++.hpp"
#include "AddNewFE.h"

// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend
// de l'orientation des aretes
//
/// ---------------------------------------------------------------
namespace  Fem2D {
	// -------------------
	// ttdcnc1_ finite element fully discontinue.
	// -------------------
	class TypeOfFE_P1ttdcnc1_: public TypeOfFE {
		public:
			static int Data [];
			static double Pi_h_coef [];

			TypeOfFE_P1ttdcnc1_ (): TypeOfFE(0, 0, 3, 1, Data, 2, 1, 3, 3, Pi_h_coef) {
				const R2 Pt [] = {R2(0.5, 0.5), R2(0.0, 0.5), R2(0.5, 0.0)};

				for (int i = 0; i < NbDoF; i++) {
					pij_alpha[i] = IPJ(i, i, 0);
					P_Pi_h[i] = Pt[i];
				}
			}

			void FB (const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat, RNMK_ &val) const;

			virtual R operator () (const FElement &K, const R2 &PHat, const KN_<R> &u, int componante, int op) const;
	};

	// on what   nu df on node  node of df
	int TypeOfFE_P1ttdcnc1_::Data [] = {6, 6, 6, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 3};
	;
	double TypeOfFE_P1ttdcnc1_::Pi_h_coef [] = {1., 1., 1.};
	R TypeOfFE_P1ttdcnc1_::operator () (const FElement &K, const R2 &PHat, const KN_<R> &u, int componante, int op) const {
		R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
		R r = 0;

		if (op == 0) {
			R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
			R ll0 = 1 - l0 * 2, ll1 = 1 - l1 * 2, ll2 = 1 - l1 * 2;
			r = u0 * ll0 + u1 * ll1 + ll2 * u2;
		} else {
			const Triangle &T = K.T;
			R2 D0 = T.H(0), D1 = T.H(1), D2 = T.H(2);
			if (op == 1) {
				r = -(D0.x * u0 + D1.x * u1 + D2.x * u2) * 2;
			} else {
				r = -(D0.y * u0 + D1.y * u1 + D2.y * u2) * 2;
			}
		}

		// cout << r << "\t";
		return r;
	}

	void TypeOfFE_P1ttdcnc1_::FB (const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat, RNMK_ &val) const {
		// const Triangle & K(FE.T);
		R2 A(K[0]), B(K[1]), C(K[2]);
		// l1(  cshrink1*(cshrink*((1,0)-G)+G)-G)+G  = 1
		R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;

		if (val.N() < 3) {
			throwassert(val.N() >= 3);
		}

		throwassert(val.M() == 1);
		// throwassert(val.K()==3 );

		val = 0;
		if (whatd[op_id]) {
			RN_ f0(val('.', 0, 0));
			f0[0] = 1 - l0 * 2;
			f0[1] = 1 - l1 * 2;
			f0[2] = 1 - l2 * 2;
		}

		if (whatd[op_dx] || whatd[op_dy]) {
			R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
			if (whatd[op_dx]) {
				RN_ f0x(val('.', 0, op_dx));
				f0x[0] = -Dl0.x * 2;
				f0x[1] = -Dl1.x * 2;
				f0x[2] = -Dl2.x * 2;
			}

			if (whatd[op_dy]) {
				RN_ f0y(val('.', 0, op_dy));
				f0y[0] = -Dl0.y * 2;
				f0y[1] = -Dl1.y * 2;
				f0y[2] = -Dl2.y * 2;
			}
		}
	}

	//
	// end ttdcnc1_
	// ------------------

	/*
	 * static void SetPtPkDCnc(R3 *Pt )
	 * {
	 *    const int d=3	;
	 *
	 *    int n=0;
	 *    R c=1./3.;
	 *    R b[]={c,c,c,c};
	 *    for(int i=0;i<= 4; ++i)
	 *    {   b[i]=0;
	 *        Pt[n++] = R3(b+1) ;
	 *        b[i]=c;
	 *    }
	 *    if(verbosity>9)
	 *        cout << " Pkdc = " << KN_<R3>(Pt,4)<<"\n";
	 *
	 * }
	 *
	 *
	 * class TypeOfFE_LagrangeDC3d: public  GTypeOfFE<Mesh3>
	 * {
	 * //typedef typename  MMesh Mesh;
	 * public:
	 * typedef   Mesh3 Mesh;
	 * typedef   Mesh::Element Element;
	 * typedef   Element::Rd Rd;
	 * typedef   Element::RdHat RdHat;
	 * static const int d=Rd::d;
	 *
	 * const int k;
	 *
	 * struct A4 {
	 *    int dfon[4];
	 *
	 *    A4(int k) {
	 *        //  (k+3)(k+2)(k+1) / 6   // d== 3
	 *
	 *        int ndf = (d== 3) ? ((k+3)*(k+2)*(k+1) / 6) :
	 *        ((d== 2) ? ((k+2)*(k+1) / 2) : k+1 );
	 *
	 *        dfon[0]=dfon[1]=dfon[2]=dfon[3]=0;
	 *        dfon[d]=ndf;
	 *
	 *        if(verbosity>9)
	 *            cout << "A4 "<<   k<< " "   <<dfon[0]<< dfon[1]<<dfon[2]<<dfon[3]<<endl;
	 *    }
	 *    operator const  int  * () const {return dfon;}
	 * };
	 *
	 * RdHat *Pt;
	 * TypeOfFE_LagrangeDC3d(int kk):
	 * //              dfon ,N,nsub(graphique) ,  const mat interpolation , discontinuous
	 * GTypeOfFE<Mesh>(A4(kk),1,Max(kk,1),true,true),k(kk)
	 * {
	 *  int n=this->NbDoF;
	 *  if(verbosity>9)
	 *      cout << "\n +++ Pdc"<<k<<" : ndof : "<< n <<endl;
	 *  SetPtPkDCnc (this->PtInterpolation);
	 *  if(verbosity>9)    cout << this->PtInterpolation<< endl;
	 *  {
	 *      for (int i=0;i<n;i++)
	 *        {
	 *          this->pInterpolation[i]=i;
	 *          this->cInterpolation[i]=0;
	 *          this->dofInterpolation[i]=i;
	 *          this->coefInterpolation[i]=1.;
	 *        }
	 *    }
	 *
	 * }
	 * ~TypeOfFE_LagrangeDC3d(){ } //cout << "TypeOfFE_LagrangeDC3d"<< this->NbDoF<<endl;}
	 *
	 * void FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, RNMK_ & val) const;
	 * virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
	 *
	 * private:
	 * TypeOfFE_LagrangeDC3d( const TypeOfFE_LagrangeDC3d &) ;
	 * void operator=( const TypeOfFE_LagrangeDC3d &) ;
	 * };
	 *
	 * void TypeOfFE_LagrangeDC3d::FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, RNMK_ & val) const
	 * {
	 *    //  1 - 3 l[i]  =>  1 sur la face i et 0 sur b des autre face
	 *  //  const Triangle & K(FE.T);
	 *  R ll[]={1.-P.sum(),P.x,P.y,P.z};
	 *  R l[]={1-ll[0]*3,1-ll[1]*3,1-ll[2]*3,1-ll[3]*3};
	 *  assert(val.N() >=Element::nv);
	 *  assert(val.M()==1 );
	 *
	 *  val=0;
	 *  RN_ f0(val('.',0,op_id));
	 *
	 *  if (whatd & Fop_D0)
	 *    {
	 *      f0[0] = l[0];
	 *      f0[1] = l[1];
	 *      f0[2] = l[2];
	 *      f0[3] = l[3];
	 *    }
	 *  if (whatd & Fop_D1)
	 *    {
	 *      R3 Dl[4];
	 *      K.Gradlambda(Dl);
	 *
	 *      //for(int i=0;i<4;++i)
	 *      //      cout << Dl[i] << endl;
	 *      if (whatd & Fop_dx)
	 *        {
	 *          RN_ f0x(val('.',0,op_dx));
	 *          f0x[0] = Dl[0].x;
	 *          f0x[1] = Dl[1].x;
	 *          f0x[2] = Dl[2].x;
	 *          f0x[3] = Dl[3].x;
	 *          f0x *= -3;
	 *        }
	 *
	 *      if (whatd & Fop_dy) {
	 *          RN_ f0y(val('.',0,op_dy));
	 *          f0y[0] = Dl[0].y;
	 *          f0y[1] = Dl[1].y;
	 *          f0y[2] = Dl[2].y;
	 *          f0y[3] = Dl[3].y;
	 *          f0y *= -3;
	 *
	 *      }
	 *
	 *      if (whatd & Fop_dz) {
	 *          RN_ f0z(val('.',0,op_dz));
	 *          f0z[0] = Dl[0].z;
	 *          f0z[1] = Dl[1].z;
	 *          f0z[2] = Dl[2].z;
	 *          f0z[3] = Dl[3].z;
	 *          f0z *= -3;
	 *
	 *      }
	 *    }
	 *  //  cout << val << endl;
	 * }
	 *
	 * R TypeOfFE_LagrangeDC3d::operator()(const FElement & K,const  R3 & PHat,const KN_<R> & u,int componante,int op) const
	 * {
	 *    R r=0;
	 *    if(k==1)
	 *      {
	 *
	 *       R u0(u(K(0))), u1(u(K(1))), u2(u(K(2))),u3(u(K(3)));
	 *
	 *    if (op==0)
	 *      {
	 *        R ll[4];
	 *        PHat.toBary(ll);
	 *          R l[]={1-ll[0]*3,1-ll[1]*3,1-ll[2]*3,1-ll[3]*3};
	 *        r = u0*l[0]+u1*l[1]+l[2]*u2+l[3]*u3;
	 *      }
	 *    else if(op==op_dx || op==op_dy || op==op_dz)
	 *      {
	 *          const Element & T=K.T;
	 *          R3 D[4];
	 *          T.Gradlambda(D);
	 *          for(int i=0;i<4;++i)
	 *              D[i] *= -3.;
	 *          if (op==op_dx)
	 *              r =  D[0].x*u0 + D[1].x*u1 + D[2].x*u2+ D[3].x*u3 ;
	 *          else if (op==op_dy)
	 *              r =  D[0].y*u0 + D[1].y*u1 + D[2].y*u2+ D[3].y*u3 ;
	 *          else
	 *              r =  D[0].z*u0 + D[1].z*u1 + D[2].z*u2+ D[3].z*u3 ;
	 *      }
	 *      }
	 *    else
	 *        ffassert(0); // to do ..
	 *    return r;
	 * }
	 */

// link with FreeFem++
	static TypeOfFE_P1ttdcnc1_ P1dc1LagrangeP1dc1;
	// static TypeOfFE_LagrangeDC3d TypeOfFE_LagrangeDC3dtt(1);

// a static variable to add the finite element to freefem++
	static AddNewFE P1dcLagrange("P1dcnc", &P1dc1LagrangeP1dc1);
// static AddNewFE3  P1dttLagrange3d("P1dcnc3d",&TypeOfFE_LagrangeDC3dtt,"P1dcnc");
}	// FEM2d namespace

// --- fin --
