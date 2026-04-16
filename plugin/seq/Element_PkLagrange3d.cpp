#include "Element_PkLagrange3d.hpp"
namespace Fem2D {
  class TypeOfFE_PKLagrange_3d : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    typedef GFElement< Mesh3 > FElement;
    int kp;      // PK
    int ndof;    //= (kp + 3) * (kp + 2) * (kp + 1) / 6;
    int dfon[4];
   vector< vector< int > > nl;    //[20][3] for P3; // returns the lambdas used for polynomial
   vector< vector< int > > cl;    //[20][3];  //  polynomial form = prod 1/cp (nl-cl)
   vector< int > cp;              //[20]; //  polynomial factors
   vector< vector< int > > pp;    //[20][4];  //bary. coordinates
    static const int d = Mesh::Rd::d;
    TypeOfFE_PKLagrange_3d(int PK);    // constructor
    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat, RNMK_ &val) const;
    void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf, int *nump) const;
    void OrientDoFOfFace(int odf, const Element &K, int f, int *pp, int inv) const;

   private:
    static int *Preparedfon(int PK) {
      static int locdfon[4];
      // Fill data
      locdfon[0] = 1;                                     // 1 dof per vertice
      locdfon[1] = PK - 1;                                // dof per edge
      locdfon[2] = (PK - 1) * (PK - 2) / 2;               // dof on the surface
      locdfon[3] = (PK - 1) * (PK - 2) * (PK - 3) / 6;    // dof in the volume
      return locdfon;
    }
  };

  TypeOfFE_PKLagrange_3d::TypeOfFE_PKLagrange_3d(int PK) : GTypeOfFE< Mesh >(Preparedfon(PK), 1, PK, false, false) {
    kp = PK;
    ndof = (PK + 1) * (PK + 2) * (PK + 3) / 6;
    nl.resize(ndof, vector< int >(PK));    
    cl.resize(ndof, vector< int >(PK));   
    pp.resize(ndof, vector< int >(4));
    cp.resize(ndof);
    // Filling dfon explicitly 
    dfon[0] = 1;                                     // 1 dof per vertice
    dfon[1] = PK - 1;                                // dof per edge
    dfon[2] = (PK - 1) * (PK - 2) / 2;               // dof on the surface
    dfon[3] = (PK - 1) * (PK - 2) * (PK - 3) / 6;    // dof in the volume*/
    // Fill nl cl cp pp 
    BasisFctPK3D(PK, nl, cl, cp,pp);
    // Fill pp
    //pp = generate_barycoordinates(PK);
    typedef Element E;
    int n = this->NbDoF;
    bool dd = verbosity > 5;
    if (dd) {
      cout << "\n +++ PK  : ndof : " << n << " " << this->PtInterpolation.N( ) << endl;
    }

    R3 *Pt = this->PtInterpolation;
    // construction of interpolation ppoint

    {
      double cc = 1. / PK;

      for (int i = 0; i < ndof; ++i) {
        Pt[i] = R3::KHat[0] * cc * pp[i][0] + R3::KHat[1] * cc * pp[i][1] + R3::KHat[2] * cc * pp[i][2] + R3::KHat[3] * cc * pp[i][3];
      }

      if (dd) {
        cout << this->PtInterpolation << endl;
      }
    }

    for (int i = 0; i < n; i++) {
      this->pInterpolation[i] = i;
      this->cInterpolation[i] = 0;
      this->dofInterpolation[i] = i;
      this->coefInterpolation[i] = 1.;
    }
  }

  // Adapted only to P4 de REF
  void TypeOfFE_PKLagrange_3d::OrientDoFOfFace(int odf, const Element &K, int f, int *pp, int inv) const {
    const int doff = odf + 4 * dfon[0] + 6 * dfon[1] + f * dfon[2];    // first dof of face f 4 vertices + 6 edges+ x* number face
    int np = K.facePermutation(f);
    int i = np / 2, j = np % (2) ? 2 : 1;
    int i0 = i, i1 = (i + j) % 3, i2 = (i + j + j) % 3;
    /*if (inv) {
      pp[doff + i0] = doff + 0;
      pp[doff + i1] = doff + 1;
      pp[doff + i2] = doff + 2;


    } else {
      pp[doff + 0] = doff + i0;
      pp[doff + 1] = doff + i1;
      pp[doff + 2] = doff + i2;
    }*/
    /*cout<<"f="<< f << " Nvface " <<(K.nvface[f][0])<<(K.nvface[f][1])<<(K.nvface[f][2])<<endl;
    cout<<"Np: "<<np<<" i "<<i<<" j "<<j<<endl;
    cout<<"i0: "<<i0<<" i1 "<<i1<<" i2 "<<i2<<endl;*/


    // position de chaque sommet{s,a,s,a,a,s}
    int spos[3]={0,kp-3,dfon[2]-1};
    /*int apos[1][3][3] = {{        // positions des arete
        {-1,  1,  3},   // arete 0-1 → pos1, arete 0-2 : pos3
        { 1, -1,  4},   // arete 1-0 → pos1, arete 1-2 : pos4
        { 3,  4, -1}    // arete 2-0 → pos3, arete 2-1 : pos4
    }};*/
   /*int apos[2][3][3] = {
   {
     {-1,  1,  4},
     { 2, -1,  6},
     { 7,  8, -1}},
 
     {
     {-1,  2,  7},
     { 1, -1,  8},
     { 4,  6, -1}}

  };*/
  /*auto apos=permutationslist ( kp);
  if (kp==7) apos[apos.size()-1]={{-1,10,7},{10,-1,6},{7,6,-1}};
  for (int i0=0; i0<3; i0++){
  for (int i1=0; i1<3; i1++){
  for (int i2=0; i2<3; i2++)
  if((i0!=i1) && (i0!=i2) && (i1!=i2)){*/


//
  ////apos[3]={{-1,8,12},{8,-1,13},{12,13,-1}};
  //  cout<<spos[i0]<<" | "<<spos[0]<<endl;
  //  cout<<spos[i1]<<" | "<<spos[1]<<endl;
  //  cout<<spos[i2]<<" | "<<spos[2]<<endl<<endl;
//
  // for (int k =0;k<2;k++){
  //      cout<<apos[k][i0][i1]<<" | "<<apos[k][0][1]<<endl;
  //      cout<<apos[k][i0][i2]<<" | "<<apos[k][0][2]<<endl;
  //      cout<<apos[k][i1][i2]<<" | "<<apos[k][1][2]<<endl;
  //        cout<<endl;
  //  }
  //exit(0);
  //*/
  
    if (inv) {
        std::vector<int> permut={i0,i1,i2};
        permut[i0]=0;
        permut[i1]=1;
        permut[i2]=2;

        //cout<<"I's "<<i0<<i1<<i2<<endl; 
        for (int i=0;i<=kp-3;i++){
            for (int j=0;j<=kp-3-i;j++){
              int k=(kp-3-j-i);
             std::vector<int> barcoord={i,j,k};
            int oldnumdof=numerotationLocIntern( kp,  barcoord[0] ,  barcoord[1], barcoord[2])   ;
            int newnumdof= numerotationPermutation( kp,barcoord, permut);
            //cout<<" doff perm: "<<oldnumdof <<" | "<<newnumdof<<endl;
            pp[ doff+ newnumdof ]= doff + oldnumdof;
            }}
        
    } 
    //}}}
    
    //exit(0);
    // This was used for the old version of the P4 (We don't need it)
    /*else {
        pp[doff + spos[0]] = doff + spos[i0];
        pp[doff + spos[1]] = doff + spos[i1];
        pp[doff + spos[2]] = doff + spos[i2];
        if(dfon[2]>3){
        pp[doff + apos[0][1]] = doff + apos[i0][i1];
        pp[doff + apos[0][2]] = doff + apos[i0][i2];
        pp[doff + apos[1][2]] = doff + apos[i1][i2];
        }
    }*/

  }



  void TypeOfFE_PKLagrange_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf, int *nump) const {
    //  faux nump don la numerotation des p -> local
    // ne marche que le cas scalaire ??? FH...
    /*
     for (int k=0;k<M.ncoef;++k)
     vdf[M.dofe[k]] += M.coef[k]*vpt(M.p[k],M.comp[k]);
     */
    int Np = M.p.N( );
    int n = this->NbDoF;
    int *p = M.p;

     for (int i = 0; i < n; ++i) {
       M.p[i] = i;
     }

    if (verbosity > 9) {
      cout << " PK  set:" << odf << " : ";
    }

    // int dof = 4+odf;
    int dof = odf;

    int nbrPermAxis = (kp % 2) ? (kp - 1) : (kp - 2);    // Number of permutation per edge
    vector< pair< int, int > > ExId = ExchangeidxVector3D(kp); // indices to exchange

    for (int e = 0; e < 6; ++e) {
      int oe = K.EdgeOrientation(e);
      for (int j = e * nbrPermAxis / 2; j < (e + 1) * nbrPermAxis / 2; j++) {
        int i1 = dof + ExId[j].first;
        int i2 = dof + ExId[j].second;
        // cout << e <<"  "<< i1 << " " << i2 << "  , " << oe <<"  p: " <<  p[i1] << " " <<  p[i2] <<" " << Np << endl;
        ffassert(i1 >= 0 && i2 >= 0);
        ffassert(i1 < Np && i2 < Np);

        if ((oe < 0) && (p[i1] < p[i2])) {
          swap(p[i1], p[i2]);
        } else if ((oe > 0) && (p[i1] > p[i2])) {
          swap(p[i1], p[i2]);
        }
      }
    }

    if (dfon[2] > 1) {
      for (int f = 0; f < 4; ++f) OrientDoFOfFace(odf, K, f, p, 1);
    }
    if (verbosity > 99) {
      cout << " " << M.p << endl;
      ;
    }
  }

  void TypeOfFE_PKLagrange_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat, RNMK_ &val) const {
    assert(val.N( ) >= ndof);    //  degrees of freedom
    assert(val.M( ) == 1);       // 3 components
    // int n = this->NbDoF;
    // -------------
    // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL
    // number (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
    // perm[3] is the local number of the vertex with the biggest global number.)
    // -------------
    R ld[4];
    PHat.toBary(ld);
    ld[0] *= R(kp);
    ld[1] *= R(kp);
    ld[2] *= R(kp);
    ld[3] *= R(kp);

    vector< int > p(ndof, 0); 
    for (int cpt = 0; cpt < p.size( ); cpt++) p[cpt] = cpt;
    {
      // int dof = 4;
      int dof = 0;

      vector< pair< int, int > > ExId = ExchangeidxVector3D(kp);
      int nbrPermAxis = (kp % 2) ? (kp - 1) : (kp - 2);    // Number of permutation per edge
      for (int e = 0; e < 6; ++e) {
        int oe = K.EdgeOrientation(e);
        if (oe < 0) {
          for (int j = e * nbrPermAxis / 2; j < (e + 1) * nbrPermAxis / 2; j++) swap(p[ExId[j].first], p[ExId[j].second]);
        }

        dof += kp - 1;
      }
      
      if (dfon[2] > 1) {
        for (int f = 0; f < 4; ++f) OrientDoFOfFace(0, K, f, p.data( ), 1);
      }
    }
    //    static int ddd = 100;
    //    ddd++;
    val = 0.;
    RN_ f0(val('.', 0, op_id));
    //  if (ddd < 20) {
    //     cout << ld[0] << " " << ld[1] << " " << ld[2] << " " << ld[3] << " ::";
    //  }

    if (whatd & Fop_D0) {
      for (int i = 0; i < ndof; ++i) {
        R fi = 1. / cp[i];

        for (int l = 0; l < kp; ++l) {
          fi *= ld[nl[i][l]] - cl[i][l];
        }

        //      if (ddd < 20) {
        //       cout << " " << fi;
        //    }

        f0[p[i]] = fi;
      }

      ///      if (ddd < 20) {
      ///       cout << endl;
      ///     }
    }

    if (whatd & (Fop_D1 | Fop_D2)) {
      R3 Dld[4], Df[ndof];
      K.Gradlambda(Dld);
      Dld[0] *= R(kp);
      Dld[1] *= R(kp);
      Dld[2] *= R(kp);
      Dld[3] *= R(kp);

      for (int i = 0; i < ndof; ++i) {
        R fi = 1. / cp[i];
        R3 &dfi = Df[p[i]];

        for (int l = 0; l < kp; ++l) {
          double ci = ld[nl[i][l]] - cl[i][l];
          dfi *= ci;
          dfi += fi * Dld[nl[i][l]];
          fi *= ci;
        }

        RN_ f0x(val('.', 0, op_dx));
        RN_ f0y(val('.', 0, op_dy));
        RN_ f0z(val('.', 0, op_dz));
        if (whatd & Fop_dx) {
          for (int i = 0; i < ndof; ++i) {
            f0x[i] = Df[i].x;
          }
        }

        if (whatd & Fop_dy) {
          for (int i = 0; i < ndof; ++i) {
            f0y[i] = Df[i].y;
          }
        }

        if (whatd & Fop_dz) {
          for (int i = 0; i < ndof; ++i) {
            f0z[i] = Df[i].z;
          }
        }

        ffassert(!(whatd & Fop_D2));    // no D2 to do !!!
      }
    }
  }




#include "Element_PkL.hpp"
  // if you need to add a newFE (higher order $ORDER for instance)
  // you can simply do it by following the steps:
  /*
    1) Declare static TypeOfFE_PKLagrange_3d PKLagrangeP$ORDER($ORDER);
    2) GTypeOfFE< Mesh3 > &Elm_PKL$ORDER_3d(PKLagrangeP$ORDER);
    3) then,  you can add your newFE using:
        AddNewFE3("PKLP$ORDER3d", &Elm_PKL$ORDER_3d); in the init fct
    4) finally recompile the file using ff-c++
  */
  static TypeOfFE_PKLagrange_3d PKLagrangeP3(3);
  GTypeOfFE< Mesh3 > &Elm_PKL3_3d(PKLagrangeP3);

  static TypeOfFE_PKLagrange_3d PKLagrangeP4(4);
  GTypeOfFE< Mesh3 > &Elm_PKL4_3d(PKLagrangeP4);

  static TypeOfFE_PKLagrange_3d PKLagrangeP5(5);
  GTypeOfFE< Mesh3 > &Elm_PKL5_3d(PKLagrangeP5);

  static TypeOfFE_PKLagrange_3d PKLagrangeP6(6);
  GTypeOfFE< Mesh3 > &Elm_PKL6_3d(PKLagrangeP6);

  static TypeOfFE_PKLagrange_3d PKLagrangeP7(7);
  GTypeOfFE< Mesh3 > &Elm_PKL7_3d(PKLagrangeP7);

  static TypeOfFE_PKLagrange_3d PKLagrangeP8(8);
  GTypeOfFE< Mesh3 > &Elm_PKL8_3d(PKLagrangeP8);

  static TypeOfFE_PKLagrange_3d PKLagrangeP9(9);
  GTypeOfFE< Mesh3 > &Elm_PKL9_3d(PKLagrangeP9);

  static TypeOfFE_PKLagrange_3d PKLagrangeP10(10);
  GTypeOfFE< Mesh3 > &Elm_PKL10_3d(PKLagrangeP10);

  static TypeOfFE_PKLagrange_3d PKLagrangeP11(11);
  GTypeOfFE< Mesh3 > &Elm_PKL11_3d(PKLagrangeP11);
  
  static void init( ) {
    AddNewFE3("PKLP33d", &Elm_PKL3_3d);
    AddNewFE3("PKLP43d", &Elm_PKL4_3d);
    AddNewFE3("PKLP53d", &Elm_PKL5_3d);
    AddNewFE3("PKLP63d", &Elm_PKL6_3d);
    AddNewFE3("PKLP73d", &Elm_PKL7_3d);
    AddNewFE3("PKLP83d", &Elm_PKL8_3d);
    AddNewFE3("PKLP93d", &Elm_PKL9_3d);
    AddNewFE3("PKLP103d", &Elm_PKL10_3d);
    AddNewFE3("PKLP113d", &Elm_PKL11_3d);

  }
}    // namespace Fem2D
LOADFUNC(Fem2D::init);