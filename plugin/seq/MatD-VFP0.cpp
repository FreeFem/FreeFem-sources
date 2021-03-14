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
//  base on: The VOLNA code for the numerical modeling of tsunami waves: Generation,
// propagation and inundation
// https://hal.archives-ouvertes.fr/hal-00454591v4/document
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
#include "AFunction_ext.hpp"

typedef const Fem2D::Mesh *pMesh;
namespace  Fem2D {
// code Georges Sadaka
typedef double R;

class LeastSq2 {
public:
    R2 l1x,l1y;
    R w2[3];
    int order;
};

void ComputeCoefs (const Mesh &Th,KN<LeastSq2> &geom,long order)
{
    // order : 1 => order 1
    // order : 2 => order 2
    // order : 3= 1|3 => order 1 if one point on boundary, order 2 other
    R d1, d2, d3;
    R w1, w2, w3, sum;
    R2 GH(1./3.,1/3.);
    R2 Pi[3];
    R w[3],di[3];
    int no1=0,no2=0,no3=0;
    //  search vexter on true boundaty
    KN<int> vb(Th.nv,0);
    for (int k=0; k<Th.nt; k++)
    for(int e=0;e<3;++e)
    {
        // search of the edge on the bdy
        int ee=e,ke=Th.ElementAdj(k,ee);
        if(ke==k || ke<0) {//
            vb[Th(k,(e+1)%3)]=1;
            vb[Th(k,(e+2)%3)]=1;
        }
    }
    for (int k=0; k<Th.nt; k++)
    {
        const Triangle &K = Th[k];
        R2  G=K(GH);
        int edgeonbdy=0;
        R sum=0;
        int nvonb =0;
        for(int e=0;e<3;++e)
        {
            nvonb+= vb[Th(k,e)];
            // search of the edge on the bdy
            int ee=e,ke=Th.ElementAdj(k,ee);
            if(ke==k || ke<0)
            {
                edgeonbdy++;
                w[e]=0;
                di[e]=0;
            }
            else {
                const Triangle &Ke=Th[ke];// adjacent triangle ...
                Pi[e]=Ke(GH)-G;
                di[e] = 1./Pi[e].norme();
                sum += di[e];
            }
        }
        bool oo2 = (order==2) || ((order>2) &&!nvonb);
        if(edgeonbdy>=2  || abs(sum) < 1e-100){no3++;oo2=2;} // Bof Bof .... Voir G. Sadaka ...
        if(oo2)
        { // order 2
            no2++;
            R2  P1= K[0]-G, P2=K[1]-G, P3= K[2]-G;
            d1 = P1.norme();
            d2 = P2.norme();
            d3 = P3.norme();
            sum = 1.0/d1 + 1.0/d2 + 1.0/d3;
            w1 = 1.0/(d1*sum);
            w2 = 1.0/(d2*sum);
            w3 = 1.0/(d3*sum);
            R l11 = w1*P1.x*w1*P1.x + w2*P2.x*w2*P2.x + w3*P3.x*w3*P3.x;
            R l12 = w1*P1.x*w1*P1.y + w2*P2.x*w2*P2.y + w3*P3.x*w3*P3.y;
            R l21 = l12;
            R  l22 = w1*P1.y*w1*P1.y + w2*P2.y*w2*P2.y + w3*P3.y*w3*P3.y;
            R det = l11*l22 - l12*l21;
            geom[k].w2[0] = w1*w1;
            geom[k].w2[1] = w2*w2;
            geom[k].w2[2] = w3*w3;
            geom[k].l1x =R2(l22/det,-l12/det);
            geom[k].l1y =R2(-l21/det,l11/det);
            geom[k].order = 2;
        }
        else
        {// order 1
            no1 ++;
            assert(edgeonbdy<2 &&sum > 1e-100);
            w[0] = di[0]/sum;
            w[1] = di[1]/sum;
            w[2] = di[2]/sum;
            if(verbosity>99 &&!mpirank)
                cout << k << " **" << sum << " " << di[0] << " " << di[1] << " " << di[2]
                << " :: " << Pi[0] << " " << Pi[1]<< " " << Pi[2] << " // "
                << w[0] << " " << w[1] << " " << w[2]<< endl;
            
            // K : X_1 is the coordinate of the 3 vertex of the triangle k
            
            R2  WP[3]={ w[0]*Pi[0], w[1]*Pi[1], w[2]*Pi[2]};
            // f = u_ke - u_k
            // L = w .*Pi => L1 = L.x , L2 = L.y,  L 2x3
            //  F = w.*f
            //  L g = F => L'L = LF  ; LF = Pi .*w^2 f
            
            R l11 = WP[0].x*WP[0].x + WP[1].x*WP[1].x + WP[2].x*WP[2].x;
            R l12 = WP[0].x*WP[0].y + WP[1].x*WP[1].y + WP[2].x*WP[2].y;
            R l22 = WP[0].y*WP[0].y + WP[1].y*WP[1].y + WP[2].y*WP[2].y;
            
            R det = l11*l22 - l12*l12;
            geom[k].w2[0] = w[0]*w[0];
            geom[k].w2[1] = w[1]*w[1];
            geom[k].w2[2] = w[2]*w[2];
            geom[k].l1x =R2(l22/det,-l12/det);
            geom[k].l1y =R2(-l12/det,l11/det);
            geom[k].order = 1;
            if(verbosity>99 &&!mpirank)
                cout << k << " : "<<det << " "<< geom[k].l1x <<" ; "<< geom[k].l1y << " : " << geom[k].w2[0] << " " << geom[k].w2[1] << " "<< geom[k].w2[2] <<endl;
        }
    }
    if(verbosity>1&&!mpirank) cout << "ComputeCoefs:: nb Triangle order1 :"<< no1 << " , order2 : "<< no2 << endl;
    if(verbosity &&no3 )
        cout <<" on mpirank: " << mpirank << " ComputeCoefs:: WARNING put order 2:  due to 2 edge border !!! number of case: " << no3 << endl;
}

MatriceMorse<R> *MatVFM01 (const Mesh *pTh )
{
    // U0   entre P0
    // u1 sortie P1
    //  u1 = A *u0
    ffassert(pTh);
    const Mesh &Th=*pTh;
    MatriceMorse<R> &A = *new MatriceMorse<R>(Th.nv,Th.nt,0,0);
    vector<vector<int> > lst(Th.nv);
    vector<vector<double> > w(Th.nv);
    for(int k=0; k<Th.nt; ++k)
        for(int j=0;j<3;++j)
            lst[Th(k,j)].push_back(k);
    R2 barycenter;
    static const double EPS = 1.e-10;
    int degr;
    R ww;
    R Ixx, Iyy, Ixy, ri2;
    R lambda, mu;
    double *sumW=new double[Th.nv];
    R2 GH(1./3.,1/3.);
    for (int i=0; i<Th.nv; ++i) {
        degr=lst[i].size();
        R dg = degr;
        w[i].resize(degr);
        R2 V=Th(i);
        
        // lagrangian multipliers computation
        Ixx = Ixy = Iyy = 0.0;
        R2 RR=R(-dg)*V;
        // cout << " RR =" <<RR << " " << V << " " << R(-degr) << " " << dg << endl;
        for(int l=0;l<degr;++l){
            int k= lst[i][l];
            R2 G = Th[k](GH),GV(V,G);
            
            ri2 = (GV,GV);
            // cout << "ri2=" << ri2 << endl;
            RR += G;
            
            Ixx += GV.x*GV.x/ri2;
            Ixy += GV.x*GV.y/ri2;
            Iyy += GV.y*GV.y/ri2;
        }
        lambda = (RR.y*Ixy - RR.x*Iyy + EPS)/(Ixx*Iyy - Ixy*Ixy + EPS);
        mu = (RR.x*Ixy - RR.y*Ixx + EPS)/(Ixx*Iyy - Ixy*Ixy + EPS);
        // now we compute the weights
        sumW[i] = 0.0;
        for (int l=0; l<degr; l++) {
            int k= lst[i][l];
            R2 G = Th[k](GH),GV(V,G);
            
            ri2 = (GV,GV);
            ww = fabs (1.0 + (lambda*(GV.x)+mu*(GV.y))/ri2);
            sumW[i] += ww;
            w[i][l]=ww;
        }
        for (int l=0; l<degr; l++) {
            int k= lst[i][l];
            w[i][l] /= sumW[i];
            A(i,k) = w[i][l] ;
        }
    }
    delete [] sumW;
    return  &A;
}
MatriceMorse<R> *Matgrads (const Mesh *pTh,long order )
{
    ffassert(pTh);
    const Mesh &Th=*pTh;
    R2 GH(1./3.,1/3.);
    int nt = Th.nt,nv=Th.nv;
    MatriceMorse<R> &A = *new MatriceMorse<R>(Th.nt*2,Th.nt,0,0);
    MatriceMorse<R> *A0 = MatVFM01(pTh);
    KN<LeastSq2> geom(Th.nt);
    vector<vector<int> > lst(Th.nv);
    vector<vector<int> > ltt(Th.nt);
    R2 Pi[3];
    int no1=0,no2=0;
    ComputeCoefs ( Th,geom,order);
    for(int k=0; k<Th.nt; ++k)
    for(int j=0;j<3;++j)
    lst[Th(k,j)].push_back(k);
    
    int icol=0;
    KN<int> col(std::max(nt,nv),icol);
    for(int k=0; k<Th.nt; ++k)
    {
        icol++;
        for(int j=0;j<3;++j)
        {
            int s=Th(k,j);
            for(int l=0; l<lst[s].size(); ++l)
            {
                int kk = lst[s][l];
                if( col[kk] != icol)
                {
                    ltt[k].push_back(kk);
                    col[kk] = icol;
                }
            }
        }
    }
    icol=0;
    col=icol;
    // ltt = pattern de la matrix
    //  calcul des coef ??
    for(int k=0; k<Th.nt; ++k)
    {
        LeastSq2 &gk = geom[k];
        const Triangle &K = Th[k];
        R2  G=K(GH);
        R2  P1[]= { (R2) K[0]-G, (R2) K[1]-G, (R2) K[2]-G};
        int lk[3];
        double cc[3] ={1,1,1};
        int aretesurlebord =0;
        for(int e=0; e<3;++e)
        {
            // recherche des arete de bord
            int ee=e,ke=Th.ElementAdj(k,ee);
            if(ke==k || ke <0) {
                cc[e]=0; // arete sur bord (pas adjacence )
                aretesurlebord++;
                lk[e]=-1;
            }
            else {
                lk[e]= ke;
                const Triangle &Ke=Th[ke];// adjacent triangle ...
                Pi[e]=Ke(GH)-G;
            }
        }
        
        if(gk.order==2)
        {// order 2
            no2++;
            for(int l=0; l<ltt[k].size(); ++l)
            {
                int kk = ltt[k][l];
                // compute coef a(k,kk) => u[i] = delta_i,kk
                // 2 case # kk \cup k == 1 or 2
                //  let s int intersection for k,kk, le numero ks dans k
                // u1[s] = A0(s,kk)
                // dU =
                //  contib k,k
                R dxa_k_kk=0;
                R dya_k_kk=0;
                if( k == kk)
                {
                    R2 B1 = -gk.w2[0]*P1[0]*cc[0]  - gk.w2[1]*P1[1]*cc[1] - gk.w2[2]*P1[2]*cc[2]   ;
                    dxa_k_kk = (B1,gk.l1x);
                    dya_k_kk = (B1,gk.l1y);
                }
                icol++;
                int  s;
                for( int i=0; i< 3; ++i)
                col[Th(kk,i)] = icol;
                for( int ks=0; ks< 3; ++ks)
                if( col[s=Th(k,ks)] == icol)
                {   // s is in k and kk
                    R us =(*A0)(s,kk); // coef for
                    R2 B1 =(us*gk.w2[ks])*P1[ks];
                    dxa_k_kk += (B1,gk.l1x);
                    dya_k_kk += (B1,gk.l1y);
                }
                A(2*k,kk) =dxa_k_kk;
                A(2*k+1,kk) =dya_k_kk;
            }
        }
        else
        {//  order 1 ...
            no1++;
            // f = -1,-1,-1
            R2 B1 = - gk.w2[0]*Pi[0] - gk.w2[1]*Pi[1] - gk.w2[2]*Pi[2]   ;
            A(2*k,k)   = (B1,gk.l1x);
            A(2*k+1,k) = (B1,gk.l1y);
            
            for(int l=0; l<3; ++l)
            {
                
                // f_i = delta_il
                int kk = lk[l];
                if( kk >=0)
                {
                    R2 B1 =  gk.w2[l]*Pi[l]    ;
                    A(2*k,kk)    = (B1,gk.l1x);
                    A(2*k+1,kk)  = (B1,gk.l1y);
                }
            }
        }
    }
    delete A0;
    if(verbosity>1 &&!mpirank) cout << "Matgrads:: nb Triangle order1 :"<< no1 << " , order2 : "<< no2 << endl;
    return &A;
}

}
newpMatrice_Creuse<R> Mat_VFD(Stack stack,const Mesh *const &pTh ,const long &order)
{
    return newpMatrice_Creuse<R>(stack,Matgrads(pTh,order));
}
newpMatrice_Creuse<R> Mat_VFD(Stack stack,const Mesh *const &pTh )
{
    return newpMatrice_Creuse<R>(stack,Matgrads(pTh,3));
}
newpMatrice_Creuse<R> Mat_VFM01(Stack stack,const Mesh *const &pTh )
{
    return newpMatrice_Creuse<R>(stack,MatVFM01(pTh));
}

KN<double>*OrderVF(Stack stack,const Mesh *const &pTh, KN<double>*const &pkb , const long &order)
{
    KN<double> &kb = *pkb;
    const Mesh &Th = *pTh;
    //   double *pkb= Add2StackOfPtr2Free(stack, new double(Th.nt));
    //   KN_<double> kb(pkb,Th.nt);
    ffassert(kb.N()==Th.nt);// verif size
    kb = double(order);
    if( order == 3 || order ==4 )
    {
        //  search vexter on true boundaty
        KN<int> vb(Th.nv,0);
        for (int k=0; k<Th.nt; k++)
            for(int e=0;e<3;++e)
            {
                // search of the edge on the bdy
                int ee=e,ke=Th.ElementAdj(k,ee);
                if(ke==k || ke<0) {//
                    vb[Th(k,(e+1)%3)]=1;
                    vb[Th(k,(e+2)%3)]=1;
                }
            }
        if(order==3)
        {
        double oo[]={2.,1.,1.,2.};
        for (int k=0; k<Th.nt; k++)
            kb[k] =oo[vb[Th(k,0)]+ vb[Th(k,1)]+vb[Th(k,2)]];
        }
        else if(order==4)
         for (int k=0; k<Th.nt; k++)
          kb[k] = (vb[Th(k,0)]+ vb[Th(k,1)]+vb[Th(k,2)])>0;
    }
    return pkb;
}
KN<double> *OrderVF(Stack stack,const Mesh *const &pTh, KN<double>*const &pkb  )
{ return OrderVF(stack,pTh,pkb, 3L);}

KN<double> *SlopeLimiterVF(Stack stack,const Mesh *const &pTh, KN<double>*const &pw,KN<double>*const &pGw, KN<double>*const &palpha  )
{
    KN<double> &w = *pw;
    KN<double> &Gw = *pGw;
    const Mesh &Th = *pTh;
    KN<double> &alpha = *palpha;
    assert( w.N() == Th.nt);
    assert( Gw.N() == Th.nt*2);
    assert( alpha.N() == Th.nt);
    R2 GH(1./3.,1/3.),PE[3]={R2(0.5,0.5),R2(0.,0.5),R2(0.5,0.)};
    for(int k=0; k< Th.nt;++k)
    {
        const Triangle &K = Th[k];
        R2  G=K(GH);

        double wk=w[k], wmink = wk, wmaxk = wk;
        for(int e=0;e<3;++e)
        {
            // search of the edge on the bdy
            int ee=e,ke=Th.ElementAdj(k,ee);
            if( ke>=0 &&ke != k) //  real adj ..
            {
                wmink = min(wmink,w[ke]);
                wmaxk = max(wmaxk,w[ke]);
            }
        }
        R2 Gk(Gw[2*k],Gw[2*k+1]);
        double nGk2=Gk.norme2();
        double a = 1;
        if(nGk2 > K.area*1e-10)
            for(int e=0;e<3;++e)
            {
                double ae=1;
                
                R2 E(K(PE[e])), GE(G,E);// G au bary of edge e
                // calcu de we
                double we = wk + (Gk,GE);
                if (we > wmaxk) ae = (wmaxk-wk)/(we-wk);
                else if(we < wmink ) ae = (wmink-wk)/(we-wk);
                a = min(a,ae);
                if(verbosity>99 &&!mpirank) cout << "      -- "<< e << " ::: " << we << " "<< wk << " " << we-wk << " " << ae << endl;
            }
        if(verbosity>99&&!mpirank) cout << k << " a " << a << " " << wmink << " " << wmaxk << " |G| "
             << G.norme()<< " : " << " : " << Gk.norme() <<  endl;
        alpha[k]=a;
    }
    return palpha;
}

static void Load_Init () {
    if(!mpirank) cout << " load: init MAT_D of VF (G. Sadaka) " << endl;
    Global.Add("MatVFD", "(", new OneOperator1s_<newpMatrice_Creuse<R> , const Mesh *>  (Mat_VFD));
    Global.Add("MatVFD", "(", new OneOperator2s_<newpMatrice_Creuse<R> , const Mesh *,long>  (Mat_VFD));
    Global.Add("MatVFM01", "(", new OneOperator1s_<newpMatrice_Creuse<R> , const Mesh *>  (Mat_VFM01));
    Global.Add("OrderVF", "(", new OneOperator3s_<KN<double>* , const Mesh *,KN<double> *,long >  (OrderVF));
    Global.Add("OrderVF", "(", new OneOperator2s_<KN<double>* , const Mesh *,KN<double> *>  (OrderVF));
    Global.Add("SlopeLimiterVF", "(", new OneOperator4s_<KN<double>* , const Mesh *,KN<double> *,KN<double> *,KN<double> *>  (SlopeLimiterVF));
}

LOADFUNC(Load_Init)
