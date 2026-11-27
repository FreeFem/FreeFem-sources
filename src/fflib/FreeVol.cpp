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
// Under Linux (with FreeFEM++), compile FreeVol.cpp into a plugin with
// ff-c++ -shared -fPIC FreeVol.cpp -I$FREEFEMPP_DIR/include -L$FREEFEMPP_DIR/lib -DFREEVOL_EXPORTS
// or simply ff-c++ -shared -fPIC FreeVol.cpp -DFREEVOL_EXPORTS
// Check that FREEFEMPP_DIR is correctly set (path to the FreeFEM++ installation).
// Make sure the plugins are in the LD_LIBRARY_PATH:
// export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/plugins
// export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.

// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...
//////////////////////////////////////////////////////////////////////
// main file FreeVol to handle finite volume schemes within FreeFem++
/////////////////////////////////////////////////////////////////////

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "FreeVol.hpp"

typedef const Fem2D::Mesh *pMesh;
typedef double R;

// ============================================================================
// ======================== CLASS LeastSq2 ==========================
// ============================================================================
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
    for (long k=0; k<Th.nt; k++)
    for(int e=0;e<3;++e)
    {
        // search of the edge on the bdy
        int ee=e;long ke=Th.ElementAdj(k,ee);
        if(ke==k || ke<0) {//
            vb[Th(k,(e+1)%3)]=1;
            vb[Th(k,(e+2)%3)]=1;
        }
    }
    for (long k=0; k<Th.nt; k++)
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
            int ee=e;long ke=Th.ElementAdj(k,ee);
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
        if(edgeonbdy>=2  || abs(sum) < 1e-30){no3++;oo2=2;} // So-so .... See G. Sadaka ...
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
            assert(edgeonbdy<2 &&sum > 1e-30);
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
    if(verbosity>5&&!mpirank) cout << "ComputeCoefs:: nb Triangle order1 :"<< no1 << " , order2 : "<< no2 << endl;
    if(verbosity &&no3 )
        cout <<" on mpirank: " << mpirank << " ComputeCoefs:: WARNING put order 2:  due to 2 edges on boundary !!! number of cases: " << no3 << endl;
}
// ============================================================================
// ================== UTILITY FUNCTIONS ===========================
// ============================================================================
MatriceMorse<R> *MatVFM01 (const Mesh *pTh, const KN<long> *pvper = 0 )
{
    // U0   input P0
    // u1 sortie P1
    //  u1 = A *u0
    // periodic (pvper: the representative of the periodic class of each vertex): the triangles around all the
    // vertices of the class of the vertex, translated next to it
    ffassert(pTh);
    const Mesh &Th=*pTh;
    const bool per = pvper && pvper->N() == Th.nv;
    MatriceMorse<R> &A = *new MatriceMorse<R>(Th.nv,Th.nt,0,0);
    vector<vector<int> > lst(Th.nv), lsv(Th.nv);// the triangles around the class, and their vertex in the class
    vector<vector<double> > w(Th.nv);
    for(long k=0; k<Th.nt; ++k)
        for(int j=0;j<3;++j) {
            long v = Th(k,j), r = per ? (*pvper)[v] : v;
            lst[r].push_back(k);
            lsv[r].push_back(v);
        }
    R2 barycenter;
    static const double EPS = 1.e-10;
    int degr;
    R ww;
    R Ixx, Iyy, Ixy, ri2;
    R lambda, mu;
    double *sumW=new double[Th.nv];
    R2 GH(1./3.,1/3.);
    for (long i=0; i<Th.nv; ++i) {
        const long ri = per ? (*pvper)[i] : i;
        const vector<int> &lk = lst[ri], &lv = lsv[ri];
        degr=lk.size();
        R dg = degr;
        w[i].resize(degr);
        R2 V=Th(i);
        
        // lagrangian multipliers computation
        Ixx = Ixy = Iyy = 0.0;
        R2 RR=R(-dg)*V;
        // cout << " RR =" <<RR << " " << V << " " << R(-degr) << " " << dg << endl;
        for(int l=0;l<degr;++l){
            long k= lk[l];
            R2 G = Th[k](GH) + (V - (R2)Th(lv[l])),GV(V,G);// translated next to the vertex i (periodic)
            
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
            long k= lk[l];
            R2 G = Th[k](GH) + (V - (R2)Th(lv[l])),GV(V,G);
            
            ri2 = (GV,GV);
            ww = fabs (1.0 + (lambda*(GV.x)+mu*(GV.y))/ri2);
            sumW[i] += ww;
            w[i][l]=ww;
        }
        for (int l=0; l<degr; l++) {
            long k= lk[l];
            w[i][l] /= sumW[i];
            A(i,k) += w[i][l] ;// += : a triangle with 2 vertices of the class (periodic)
        }
    }
    delete [] sumW;
    return  &A;
}
MatriceMorse<R> *Matgrads (const Mesh *pTh,long order, const KN<long> *pvper = 0 )
{
    // periodic (pvper): the vertices of a periodic class are the same vertex (the neighbours on the other side)
    ffassert(pTh);
    const Mesh &Th=*pTh;
    const bool per = pvper && pvper->N() == Th.nv;
#define FREEVOL_REPV(v) (per ? (*pvper)[v] : (v))
    R2 GH(1./3.,1/3.);
    long nt = Th.nt,nv=Th.nv;
    MatriceMorse<R> &A = *new MatriceMorse<R>(Th.nt*2,Th.nt,0,0);
    MatriceMorse<R> *A0 = MatVFM01(pTh, pvper);
    KN<LeastSq2> geom(Th.nt);
    vector<vector<int> > lst(Th.nv);
    vector<vector<int> > ltt(Th.nt);
    R2 Pi[3];
    int no1=0,no2=0;
    ComputeCoefs ( Th,geom,order);
    for(long k=0; k<Th.nt; ++k)
    for(int j=0;j<3;++j)
    lst[FREEVOL_REPV(Th(k,j))].push_back(k);
    
    int icol=0;
    KN<int> col(std::max(nt,nv),icol);
    for(long k=0; k<Th.nt; ++k)
    {
        icol++;
        for(int j=0;j<3;++j)
        {
            long s=FREEVOL_REPV(Th(k,j));
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
    // ltt = pattern of the matrix
    //  computation of the coefs ??
    for(long k=0; k<Th.nt; ++k)
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
            // search for the boundary edges
            int ee=e;long ke=Th.ElementAdj(k,ee);
            if(ke==k || ke <0) {
                cc[e]=0; // edge on the boundary (no adjacency)
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
                //  let s int intersection for k,kk, its number ks in k
                // u1[s] = A0(s,kk)
                // dU =
                //  contib k,k
                R dxa_k_kk=0;
                R dya_k_kk=0;
                if( k == kk)
                {
                    // the 3 vertices of k (index of vertex, not of edge: no cc[e] here)
                    R2 B1 = -gk.w2[0]*P1[0] - gk.w2[1]*P1[1] - gk.w2[2]*P1[2];
                    dxa_k_kk = (B1,gk.l1x);
                    dya_k_kk = (B1,gk.l1y);
                }
                icol++;
                int  s;
                for( int i=0; i< 3; ++i)
                col[FREEVOL_REPV(Th(kk,i))] = icol;
                for( int ks=0; ks< 3; ++ks)
                if( col[FREEVOL_REPV(s=Th(k,ks))] == icol)
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
    if(verbosity>5 &&!mpirank) cout << "Matgrads:: nb Triangle order1 :"<< no1 << " , order2 : "<< no2 << endl;
    return &A;
#undef FREEVOL_REPV
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

KN<double>*OrderVF(const Mesh *const &pTh, KN<double>*const &pkb , const long &order)
{
    KN<double> &kb = *pkb;
    const Mesh &Th = *pTh;
    //   double *pkb= Add2StackOfPtr2Free(stack, new double(Th.nt));
    //   KN_<double> kb(pkb,Th.nt);
    ffassert(kb.N()==Th.nt);// check size
    kb = double(order);
    if( order == 3)
    {
        //  search vexter on true boundaty
        KN<int> vb(Th.nv,0);
        for (long k=0; k<Th.nt; k++)
            for(int e=0;e<3;++e)
            {
                // search of the edge on the bdy
                int ee=e;long ke=Th.ElementAdj(k,ee);
                if(ke==k || ke<0) {//
                    vb[Th(k,(e+1)%3)]=1;
                    vb[Th(k,(e+2)%3)]=1;
                }
            }
        double oo[]={2.,1.,1.,2.};
        for (long k=0; k<Th.nt; k++)
            kb[k] =oo[vb[Th(k,0)]+ vb[Th(k,1)]+vb[Th(k,2)]];
    }
    return pkb;
}
KN<double> *OrderVF(const Mesh *const &pTh, KN<double>*const &pkb  )
{ return OrderVF(pTh,pkb, 3L);}

double sumwtabprod(const KN<double> &aa,const KN<long> &bb,const KN<double> &ww,const long &normalize,KN<double> &cc) {
 //aa: input array
 //bb: local number of relevant values in array
 //ww: weights
 // normalize=1 indicates that we have to normalize so that for each local array the sum of weights is 1
 //cc: result
 long nvp = aa.N();
 long nv = bb.N();
 long np = nvp/nv;
 ffassert((aa.N()==nv*np)&&(cc.N()==nv)&&(ww.N()==nv*np));
 KN<double> sumw(nv);
 if (normalize==1) {
  for (long i=0;i<nv;++i) {
   double sumwloc=0.;
   for (long j=0;j<bb(i);++j) {
    sumwloc+=ww(i*np+j);
   }
   sumw(i)=sumwloc;
  }
 }
 else sumw=1.;
 for (long i=0;i<nv;++i) {
  double sumloc=0.;
  for (long j=0;j<bb(i);++j) {
   sumloc+=aa(i*np+j)*ww(i*np+j);
  }
  cc(i)=sumloc/sumw(i);
 }
 return 0.;
}// end of sumwtabprod

double limiterfunc(double amin,double amax,const long &ilimitertype){
    //ilimitertype = 2: minmod, 3:Sweby, 4: Van Albada 1
    static const double betasweby=1.5;
    double result;
    double asum=amin+amax,absdiff=abs(amax-amin);
    amin=.5*(asum-absdiff);amax=.5*(asum+absdiff);
    if (ilimitertype==2) {
      result=max(amin,0.)+min(amax,0.);//minmod limiter
    } else if (ilimitertype==3) {
      double a1=min(abs(amin),betasweby*amax*amin/(abs(amin)+1.e-30));
      double a2=min(betasweby*abs(amin),amax*amin/(abs(amin)+1.e-30));
      result=amin/(abs(amin)+1.e-30)*max(0.,max(a1,a2));//Sweby limiter
    } else if (ilimitertype==4) {
      result=max(0.,amin*amax)*(amin+amax)/(amin*amin+amax*amax+1.e-30);// Van Albada 1 limiter
    } else {cout << "unknown limiter type" << endl;ffassert(0);}
    return result;
}

// gas dynamics numerical flux (HLLC/Suliciu)
double numfluxgd(const double &gamma,const double* statele,const double* stateri,const double* normalp,double &absvp,double* Fluxesvect)
{
  // statele, stateri: left state, right state
  // normalp: normal
  // absvp: max propagation speed
  // Fluxesvect: numerical flux
  double rhol=statele[0],rhor=stateri[0],uxl=statele[1],uxr=stateri[1],uyl=statele[2],uyr=stateri[2],el=statele[3],er=stateri[3],norx=normalp[0],nory=normalp[1];
  double uxlprime=(norx)*uxl+(nory)*uyl;//rotation in normal/tangent directions
  double uylprime=-(nory)*uxl+(norx)*uyl;
  double uxrprime=(norx)*uxr+(nory)*uyr;
  double uyrprime=-(nory)*uxr+(norx)*uyr;
  double numrhouxprime,numrhouyprime;// num flux in normal/tangent directions

  if ((rhol<=0.)&&(rhor<=0.)) {//test vacuum
   absvp=max(abs(uxlprime),abs(uxrprime));
   Fluxesvect[0]=0.;
   numrhouxprime=0.;
   numrhouyprime=0.;
   Fluxesvect[3]=0.;
   //Fluxesvect[4]=0.;//entropy flux
  }
  else {
   double povrhol=(gamma-1.)*el;
   double povrhor=(gamma-1.)*er;
   double pl=povrhol*rhol;
   double pr=povrhor*rhor;
   double sl=max(1.e-30,sqrt(gamma*povrhol));
   double sr=max(1.e-30,sqrt(gamma*povrhor));
   double al,ar;
   if (pr-pl>=0.) {
    al=sl+.5*(gamma+1.)*max(0.,(pr-pl)/max(1.e-30,rhor)/sr+uxlprime-uxrprime);
    if (rhol>0.) {ar=sr+.5*(gamma+1.)*max(0.,(pl-pr)/max(1.e-30,rhol)/al+uxlprime-uxrprime);}
    else {ar=sr;}
    }
   else {
    ar=sr+.5*(gamma+1.)*max(0.,(pl-pr)/max(1.e-30,rhol)/sl+uxlprime-uxrprime);
    if (rhor>0.) {al=sl+.5*(gamma+1.)*max(0.,(pr-pl)/max(1.e-30,rhor)/ar+uxlprime-uxrprime);}
    else {al=sl;}
   }
   double cl=al*rhol;
   double cr=ar*rhor;
   double povcl,povcr;
   if (al>0.) {povcl=povrhol/al;}
   else {povcl=0.;}
   if (ar>0.) {povcr=povrhor/ar;}
   else {povcr=0.;}
   // eigenvalues
   double vp1=uxlprime-al;
   double us=(cl*uxlprime+cr*uxrprime-(pr-pl))/(cl+cr);
   double vp3=uxrprime+ar;
   absvp=max(abs(vp1),abs(vp3));
   // according to the intermediate state status
   if (vp1>=0.) {
    Fluxesvect[0]=rhol*uxlprime;
    numrhouxprime=rhol*uxlprime*uxlprime+pl;
    numrhouyprime=rhol*uylprime*uxlprime;
    Fluxesvect[3]=(rhol*(uxlprime*uxlprime/2.+uylprime*uylprime/2.+el)+pl)*uxlprime;
    //Fluxesvect[4]=-rhol*log(pow(max(1.e-30,povrhol),(1./gamma))/pow(max(1.e-30,rhol),(1.-1./gamma)))*uxlprime;//entropy flux
   }
   else if (us>=0.) {
    double rhos=cl/(al+(cr*(uxrprime-uxlprime)-(pr-pl))/(cl+cr));
    double pisovcl=(cr*povcl+pr-cr*(uxrprime-uxlprime))/(cl+cr);
    double es=el-.5*povcl*povcl+.5*pisovcl*pisovcl;
    Fluxesvect[0]=rhos*us;
    numrhouxprime=rhos*us*us+pisovcl*cl;
    numrhouyprime=rhos*uylprime*us;
    Fluxesvect[3]=(rhos*(us*us/2.+uylprime*uylprime/2.+es)+pisovcl*cl)*us;
    //Fluxesvect[4]=-rhos*log(pow(max(1.e-30,povrhol),(1./gamma))/pow(max(1.e-30,rhol),(1.-1./gamma)))*us;//entropy flux
    }
   else if (vp3>=0.){
    double rhos=cr/(ar+(cl*(uxrprime-uxlprime)+(pr-pl))/(cl+cr));
    double pisovcr=(pl+cl*povcr-cl*(uxrprime-uxlprime))/(cl+cr);
    double es=er-.5*povcr*povcr+.5*pisovcr*pisovcr;
    Fluxesvect[0]=rhos*us;
    numrhouxprime=rhos*us*us+pisovcr*cr;
    numrhouyprime=rhos*uyrprime*us;
    Fluxesvect[3]=(rhos*(us*us/2.+uyrprime*uyrprime/2.+es)+pisovcr*cr)*us;
    //Fluxesvect[4]=-rhos*log(pow(max(1.e-30,povrhor),(1./gamma))/pow(max(1.e-30,rhor),(1.-1./gamma)))*us;//entropy flux
   }
   else {
    Fluxesvect[0]=rhor*uxrprime;
    numrhouxprime=rhor*uxrprime*uxrprime+pr;
    numrhouyprime=rhor*uyrprime*uxrprime;
    Fluxesvect[3]=(rhor*(uxrprime*uxrprime/2.+uyrprime*uyrprime/2.+er)+pr)*uxrprime;
    //Fluxesvect[4]=-rhor*log(pow(max(1.e-30,povrhor),(1./gamma))/pow(max(1.e-30,rhor),(1.-1./gamma)))*uxrprime;//entropy flux
   }
  }//end of test vacuum
  // recover num flux in x,y directions
  Fluxesvect[1]=(norx)*numrhouxprime-(nory)*numrhouyprime;
  Fluxesvect[2]=(nory)*numrhouxprime+(norx)*numrhouyprime;
  return 0.0;    // dummy return value.
}//end of numfluxgd
// ============================================================================
// ===== FE functions of FreeFEM in 1D (meshL), 2D (mesh), 3D (mesh3)  ========
// ============================================================================
// the fespace of a FE function (FEbase<double,v_fes>, FEbase<double,v_fes3>, FEbase<double,v_fesL>)
template< class FEB > static typename FEB::FESpace &FreeVol_Vh(FEB *b) { return *b; }
// the same code for the FEbase of the dimension of the function
#define FREEVOL_FE_CASE(EXPR) \
    switch (dim) { \
        case 2: { pferbase b = f2.first; EXPR; } \
        case 3: { pf3rbase b = f3.first; EXPR; } \
        case 1: { pfLrbase b = fL.first; EXPR; } \
        case 4: { pfSrbase b = fS.first; EXPR; } \
    }

KN<double> *FreeVol_FE::x() const { FREEVOL_FE_CASE(return b ? b->x() : 0) return 0; }
int FreeVol_FE::comp() const { return dim == 2 ? f2.second : (dim == 3 ? f3.second : (dim == 4 ? fS.second : fL.second)); }
int FreeVol_FE::N() const { FREEVOL_FE_CASE(return FreeVol_Vh(b).N) return 0; }
long FreeVol_FE::NbOfDF() const { FREEVOL_FE_CASE(return FreeVol_Vh(b).NbOfDF) return 0; }
long FreeVol_FE::NbOfNodes() const { FREEVOL_FE_CASE(return FreeVol_Vh(b).NbOfNodes) return 0; }
int FreeVol_FE::MaxNbDFPerElement() const { FREEVOL_FE_CASE(return FreeVol_Vh(b).MaxNbDFPerElement) return 0; }
long FreeVol_FE::dof(long k, int i) const { FREEVOL_FE_CASE(return FreeVol_Vh(b)(k, i)) return 0; }
const void *FreeVol_FE::mesh() const { FREEVOL_FE_CASE(return &FreeVol_Vh(b).Th) return 0; }
const RefCounter *FreeVol_FE::meshrc() const { FREEVOL_FE_CASE(return &FreeVol_Vh(b).Th) return 0; }
KN_<double> FreeVol_FE::values() const {
    KN<double> &pv = *x();
    long nn = NbOfDF() / N();    // number of dof by component
    return pv(SubArray(nn, comp(), N()));
}
#undef FREEVOL_FE_CASE

// ============================================================================
// ========== geometry of the mesh (the same in 1D, 2D, 3D)   =================
// ============================================================================
static void FreeVol_GeomAlloc(FreeVol_Geom &g)
{
    long nkf = g.nt * g.nf;
    g.meas.resize(g.nt); g.G.resize(g.nt);
    g.adj.resize(nkf); g.adjf.resize(nkf); g.fmeas.resize(nkf); g.N.resize(nkf); g.F.resize(nkf);
    g.flab.resize(nkf); g.fvert.resize(nkf * g.nvf);
    g.flab = 0L;
    g.nper = 0; g.pshift.resize(0); g.vper.resize(0); g.plab.resize(0);    // no periodic (SetPeriodic)
}

static void FreeVol_GeomCount(FreeVol_Geom &g)
{
    g.nbdry = 0;
    for (long kf = 0; kf < g.nt * g.nf; ++kf) if (g.adj[kf] < 0) g.nbdry++;
    g.nfaces = (g.nt * g.nf - g.nbdry) / 2 + g.nbdry;
}

// 2D: the mesh of FreeFEM (the computations as before the generic version: same results)
void FreeVol_Geom::build(const Mesh &Th)
{
    static const R2 GH(1./3.,1/3.),PE[3]={R2(0.5,0.5),R2(0.,0.5),R2(0.5,0.)};
    dim = 2; tdim = 2; fedim = 2; surf = false; nf = 3; nvf = 2; nt = Th.nt; nv = Th.nv; pTh = &Th;
    FreeVol_GeomAlloc(*this);
    for (long k = 0; k < nt; ++k) {
        const Triangle &K = Th[k];
        R2 Gk = K(GH);
        G[k] = R3(Gk.x, Gk.y, 0.);
        meas[k] = K.area;
        for (int e = 0; e < 3; ++e) {
            long kf = 3 * k + e;
            int ee = e;
            long ko = Th.ElementAdj(k, ee);
            adj[kf] = (ko == k || ko < 0) ? -1 : ko;
            adjf[kf] = adj[kf] < 0 ? -1 : ee;
            R2 E = K.Edge(e);
            double le = K.lenEdge(e);
            R2 Nn = E.perp() / -le;
            N[kf] = R3(Nn.x, Nn.y, 0.);
            fmeas[kf] = le;
            R2 Fe = K(PE[e]);
            F[kf] = R3(Fe.x, Fe.y, 0.);
            fvert[kf * 2] = Th(k, (e + 1) % 3);
            fvert[kf * 2 + 1] = Th(k, (e + 2) % 3);
        }
    }
    for (long ib = 0; ib < Th.neb; ++ib) {    // the labels of the true boundary (not the internal borders)
        int ie;
        long k = Th.BoundaryElement(ib, ie);
        if (adj[3 * k + ie] < 0) flab[3 * k + ie] = Th.be(ib).lab;
    }
    FreeVol_GeomCount(*this);
}

// 1D (meshL) and 3D (mesh3): the generic meshes of FreeFEM
template< class MeshT >
static void FreeVol_GeomBuildG(FreeVol_Geom &g, const MeshT &Th, int dim)
{
    typedef typename MeshT::Element E;
    g.dim = dim; g.tdim = dim; g.fedim = dim; g.surf = false;
    g.nf = E::nea; g.nvf = E::nva; g.nt = Th.nt; g.nv = Th.nv; g.pTh = &Th;
    FreeVol_GeomAlloc(g);
    for (long k = 0; k < g.nt; ++k) {
        const E &K = Th[k];
        R3 Gk(0., 0., 0.);
        for (int i = 0; i < E::nv; ++i) Gk += (R3)K[i];
        g.G[k] = Gk / E::nv;
        g.meas[k] = K.mesure();
        for (int e = 0; e < g.nf; ++e) {
            long kf = k * g.nf + e;
            int j = e;
            long ko = Th.ElementAdj(k, j);    // < 0 on the boundary
            g.adj[kf] = ko >= 0 ? ko : -1;
            g.adjf[kf] = ko >= 0 ? j : -1;
            R3 Nn = K.N(e);                    // outward, |N| = 2 area (3D), length (1D)
            double ln = Nn.norme();
            g.N[kf] = Nn / ln;
            g.fmeas[kf] = dim == 3 ? 0.5 * ln : 1.;
            R3 Fe(0., 0., 0.);
            for (int l = 0; l < g.nvf; ++l) {
                int iv = E::nvadj[e][l];
                g.fvert[kf * g.nvf + l] = Th(k, iv);
                Fe += (R3)K[iv];
            }
            g.F[kf] = Fe / g.nvf;
        }
    }
    for (long ib = 0; ib < Th.nbe; ++ib) {
        int ie;
        long k = Th.BoundaryElement(ib, ie);
        if (g.adj[k * g.nf + ie] < 0) g.flab[k * g.nf + ie] = Th.be(ib).lab;
    }
    FreeVol_GeomCount(g);
}
void FreeVol_Geom::build(const Mesh3 &Th) { FreeVol_GeomBuildG(*this, Th, 3); }
void FreeVol_Geom::build(const MeshL &Th) { FreeVol_GeomBuildG(*this, Th, 1); }
// surface (meshS): the triangles in 3D, the faces are the edges, N = the unit outward conormal (in the plane of
// the triangle, orthogonal to the edge), Nc = the unit normal of the triangle
void FreeVol_Geom::build(const MeshS &Th)
{
    typedef MeshS::Element E;
    dim = 3; tdim = 2; fedim = 4; surf = true; nf = 3; nvf = 2; nt = Th.nt; nv = Th.nv; pTh = &Th;
    FreeVol_GeomAlloc(*this);
    Nc.resize(nt);
    for (long k = 0; k < nt; ++k) {
        const E &K = Th[k];
        R3 P[3] = {K[0], K[1], K[2]};
        G[k] = (P[0] + P[1] + P[2]) / 3.;
        R3 nk = (P[1] - P[0]) ^ (P[2] - P[0]);
        double ln = nk.norme();
        meas[k] = 0.5 * ln;
        Nc[k] = nk / ln;
        for (int e = 0; e < 3; ++e) {
            long kf = 3 * k + e;
            int j = e;
            long ko = Th.ElementAdj(k, j);    // < 0 on the boundary
            adj[kf] = ko >= 0 ? ko : -1;
            adjf[kf] = ko >= 0 ? j : -1;
            int ia = E::nvadj[e][0], ib = E::nvadj[e][1];
            R3 A = P[ia], B = P[ib], T = B - A;
            R3 Nn = T ^ Nc[k];               // in the plane of the triangle, orthogonal to the edge
            if ((Nn, P[e] - A) > 0) Nn = -Nn;  // outward: away from the opposite vertex
            N[kf] = Nn / Nn.norme();
            fmeas[kf] = T.norme();
            F[kf] = (A + B) * 0.5;
            fvert[kf * 2] = Th(k, ia);
            fvert[kf * 2 + 1] = Th(k, ib);
        }
    }
    for (long ib = 0; ib < Th.nbe; ++ib) {
        int ie;
        long k = Th.BoundaryElement(ib, ie);
        if (adj[3 * k + ie] < 0) flab[3 * k + ie] = Th.be(ib).lab;
    }
    FreeVol_GeomCount(*this);
}
void FreeVol_Geom::build(const FreeVol_FE &fe)
{
    if (fe.dim == 2) build(*(const Mesh *)fe.mesh());
    else if (fe.dim == 3) build(*(const Mesh3 *)fe.mesh());
    else if (fe.dim == 1) build(*(const MeshL *)fe.mesh());
    else if (fe.dim == 4) build(*(const MeshS *)fe.mesh());
    else ExecError("FreeVol: FE function on a mesh, a mesh3, a meshL or a meshS expected");
}

// 1D and 3D: the gradient in each cell by least squares on the centroids of the neighbours by face
// (exact for an affine function when the cell has dim independent neighbours, 0 otherwise)
// surface: the least squares in the plane of the triangle (basis e1, e2), the gradient tangent (3 components)
// the translation of the neighbour by the face kf (periodic: glued face, 0 otherwise): the centroid of the
// neighbour seen from the cell is G[adj[kf]] + FreeVol_PShift(g, kf)
static inline R3 FreeVol_PShift(const FreeVol_Geom &g, long kf) { return g.pshift.N() ? g.pshift[kf] : R3(); }

static MatriceMorse<R> *FreeVol_MatgradsLSQ(const FreeVol_Geom &g)
{
    int d = g.dim;
    MatriceMorse<R> &A = *new MatriceMorse<R>(d * g.nt, g.nt, 0, 0);
    if (g.surf) {
        for (long k = 0; k < g.nt; ++k) {
            R3 e1 = g.F[3 * k] - g.G[k];
            e1 -= g.Nc[k] * (e1, g.Nc[k]);
            e1 /= e1.norme();
            R3 e2 = g.Nc[k] ^ e1;
            double M[2][2] = {{0, 0}, {0, 0}};
            int nb = 0;
            for (int e = 0; e < 3; ++e) {
                long j = g.adj[3 * k + e];
                if (j < 0) continue;
                R3 D = g.G[j] - g.G[k] + FreeVol_PShift(g, 3 * k + e);
                double w2 = 1. / (D, D), d1 = (D, e1), d2 = (D, e2);
                M[0][0] += w2 * d1 * d1; M[0][1] += w2 * d1 * d2; M[1][1] += w2 * d2 * d2;
                nb++;
            }
            M[1][0] = M[0][1];
            double det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
            if (nb < 2 || !(fabs(det) > 1e-10)) continue;    // not enough neighbours: order 1 in this cell
            double I[2][2] = {{M[1][1] / det, -M[0][1] / det}, {-M[1][0] / det, M[0][0] / det}};
            for (int e = 0; e < 3; ++e) {
                long j = g.adj[3 * k + e];
                if (j < 0) continue;
                R3 D = g.G[j] - g.G[k] + FreeVol_PShift(g, 3 * k + e);
                double w2 = 1. / (D, D), d1 = (D, e1), d2 = (D, e2);
                double c1 = w2 * (I[0][0] * d1 + I[0][1] * d2), c2 = w2 * (I[1][0] * d1 + I[1][1] * d2);
                R3 cg = e1 * c1 + e2 * c2;
                for (int a = 0; a < 3; ++a) {
                    A(3 * k + a, j) += cg[a];
                    A(3 * k + a, k) -= cg[a];
                }
            }
        }
        return &A;
    }
    for (long k = 0; k < g.nt; ++k) {
        double M[3][3] = {{0,0,0},{0,0,0},{0,0,0}}, scale = 0.;
        int nb = 0;
        for (int e = 0; e < g.nf; ++e) {
            long j = g.adj[k * g.nf + e];
            if (j < 0) continue;
            R3 D = g.G[j] - g.G[k] + FreeVol_PShift(g, k * g.nf + e);
            double w2 = 1. / (D, D);
            for (int a = 0; a < d; ++a) for (int b = 0; b < d; ++b) M[a][b] += w2 * D[a] * D[b];
            scale = max(scale, sqrt((D, D)));
            nb++;
        }
        if (nb < d) continue;    // not enough neighbours: no slope (order 1 in this cell)
        double I[3][3], det;
        if (d == 1) { det = M[0][0]; I[0][0] = 1. / det; }
        else if (d == 2) {
            det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
            I[0][0] = M[1][1] / det; I[0][1] = -M[0][1] / det; I[1][0] = -M[1][0] / det; I[1][1] = M[0][0] / det;
        } else {
            det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
                int a1 = (b + 1) % 3, a2 = (b + 2) % 3, b1 = (a + 1) % 3, b2 = (a + 2) % 3;
                I[a][b] = (M[a1][b1] * M[a2][b2] - M[a1][b2] * M[a2][b1]) / det;   // inverse = cofactors^T / det
            }
        }
        if (!(fabs(det) > 1e-10)) continue;    // the M matrix is dimensionless (weights 1/|D|^2)
        for (int e = 0; e < g.nf; ++e) {
            long j = g.adj[k * g.nf + e];
            if (j < 0) continue;
            R3 D = g.G[j] - g.G[k] + FreeVol_PShift(g, k * g.nf + e);
            double w2 = 1. / (D, D);
            for (int a = 0; a < d; ++a) {
                double c = 0.;
                for (int b = 0; b < d; ++b) c += I[a][b] * w2 * D[b];
                A(d * k + a, j) += c;
                A(d * k + a, k) -= c;
            }
        }
    }
    return &A;
}

// ============================================================================
// ======================== CLASS FreeVol_model ================================
// ============================================================================
double FreeVol_model::SlopeLimiterVF(const KN<double> &w, const KN<double> &Gw, KN<double> &alpha, const long &ilimiter) const
{
    // ilimiter = 0: no limiter, 1: maximum principle type, 2: minmod, 3: Sweby, 4: Van Albada type, 12=1+2, 13=1+3, 14=1+4
    const FreeVol_Geom &g = geom;
    int d = g.dim;
    assert( w.N() == g.nt);//w: values in cells
    assert( Gw.N() == g.nt*d);//Gw: gradient in cells
    assert( alpha.N() == g.nt);//alpha: computed limitation
    for(long k=0; k< g.nt;++k)//loop on cells
    {
        double wk=w[k], wmink = wk, wmaxk = wk;
        for(int e=0;e<g.nf;++e)//loop on the faces of cell k
        {
            long ke=g.adj[k*g.nf+e];
            if(ke >= 0) // face is not on the boundary
            {
                wmink = min(wmink,w[ke]);
                wmaxk = max(wmaxk,w[ke]);
            }
        }
        double a = 1;
        if (ilimiter!=0){
          if ((ilimiter==1)||(ilimiter>10)){
              for(int e=0;e<g.nf;++e)
              {
                  double ae=1;
                  R3 GE = g.F[k*g.nf+e]-g.G[k];//vector joining the centroid of the cell to the centroid of face e
                  double incwe = 0.;
                  for(int c=0;c<d;++c) incwe += Gw[d*k+c]*GE[c];
                  double we = wk + incwe;//value of w at the centroid of the face
                  if (we > wk) ae = .9*(wmaxk-wk)/(we-wk);
                  else if(we < wk ) ae = .9*(wk-wmink)/(wk-we);
                  a = min(a,ae);
              }
          }
          if (ilimiter>1){
              long iilimiter=ilimiter%10;
              for(int e=0;e<g.nf;++e)
              {
                  double ae=1;
                  R3 GE = g.F[k*g.nf+e]-g.G[k];
                  double incwe = 0.;//increment w(E)-w(G)
                  for(int c=0;c<d;++c) incwe += Gw[d*k+c]*GE[c];
                  long ke=g.adj[k*g.nf+e];
                  if(ke >= 0){ // face is not on the boundary
                   double incwee = 0.;//increment computed with the neighbouring gradient
                   for(int c=0;c<d;++c) incwee += Gw[d*ke+c]*GE[c];
                   ae=limiterfunc(incwe+0.*(incwe-incwee),incwee-0.*(incwee-incwe),iilimiter)*incwe/(incwe*incwe+1.e-20);
                  }
                  a = min(a,ae);
              }
          }
        }
        alpha[k]=a;
    }//end of loop on cells
    return 0.;
}//end of FreeVol_model::SlopeLimiterVF

// the slope on the boundary: 0 (ibdryslope=0) or tangent to the boundary (ibdryslope=1):
// one boundary face: projection on the face; 3D, 2 boundary faces: on their common edge;
// tdim boundary faces or more (corner): 0 (the only slope tangent to all of them)
void FreeVol_model::BdrySlope(KN<double> &Gw, const long &ibdryslope) const
{
    const FreeVol_Geom &g = geom;
    int d = g.dim;
    assert( Gw.N() == g.nt*d);
    for(long k=0; k< g.nt;++k)
    {
        R3 Nb[4];
        int nb=0;// the normals of the physical boundary faces of k
        for(int e=0;e<g.nf;++e)
        {
            long kf=k*g.nf+e;
            if(g.adj[kf] < 0 && g.flab[kf] != -111111) Nb[nb++]=g.N[kf];
        }
        if(nb==0) continue;
        if(ibdryslope==0 || nb>=g.tdim) {// zero slope
            for(int c=0;c<d;++c) Gw[d*k+c] = 0.;
        } else if(nb==1) {// imposing Gw.N=0 on the boundary face
            double GwN=0.;
            for(int c=0;c<d;++c) GwN += Gw[d*k+c]*Nb[0][c];
            for(int c=0;c<d;++c) Gw[d*k+c] = Gw[d*k+c] - GwN*Nb[0][c];
        } else {// 3D, 2 boundary faces: the slope along their common edge t = N1 x N2
            R3 t = Nb[0]^Nb[1];
            double tt=(t,t), Gwt=0.;
            for(int c=0;c<d;++c) Gwt += Gw[d*k+c]*t[c];
            for(int c=0;c<d;++c) Gw[d*k+c] = tt > 1e-20 ? Gwt*t[c]/tt : 0.;
        }
    }
}//end of FreeVol_model::BdrySlope

double FreeVol_model::absvp(long kf,const R3 &N,const double *U) const
{
    double l[Neqmax];
    lambda(kf,N,U,l);
    double maxabsl=0.;
    for (int j=0;j<Neq;++j) maxabsl=max(maxabsl,abs(l[j]));
    return maxabsl;
}

double FreeVol_model::lambdamin(long kf,const R3 &N,const double *U) const
{
    double l[Neqmax];
    lambda(kf,N,U,l);
    double minl=l[0];
    for (int j=0;j<Neq;++j) minl=min(minl,l[j]);
    return minl;
}

double FreeVol_model::lambdamax(long kf,const R3 &N,const double *U) const
{
    double l[Neqmax];
    lambda(kf,N,U,l);
    double maxl=l[0];
    for (int j=0;j<Neq;++j) maxl=max(maxl,l[j]);
    return maxl;
}

// the reconstructed value of the unknowns at the centroid of the face kf (gU = 0: order 1)
void FreeVol_model::GetU(long kf,const KN_d *U,const KN<double> *gU,double *u) const
{
    const FreeVol_Geom &g = geom;
    int d = g.dim;
    long k=kf/g.nf;
    R3 GE = g.F[kf]-g.G[k];
    for(long j=0; j<Neq;++j){
       double ukj = U[j][k];
       if(ord>1 && gU){
          double inc = 0.;
          for(int c=0;c<d;++c) inc += gU[j][d*k+c]*GE[c];
          ukj += inc;
       }
       u[j] = ukj;
    }
}

void FreeVol_model::GetUBC(long kf,const KN_d *U,double *u) const
{
    long k=kf/geom.nf;
    for(int j=0; j<Neq;++j){
       double ukj = U[j][k];//given Dirichlet value
       u[j] = ukj;
    }
}

// ============================================================================
// generic finite volume scheme (all the models, all the dimensions)
// ============================================================================
// the limited slopes of the Neq unknowns in each cell: gU[j][dim*k+c]
void FreeVol_model::LimitedSlopes(KN_d *U, KN<double> *gU) const
{
    const FreeVol_Geom &g = geom;
    int d = g.dim;
    KN<double> u(g.nt), grad_u(d*g.nt), alpha_u(g.nt);
    MatriceMorse<R> &MATGRAD = *(matgrad);
    for(long j=0; j<Neq; ++j){
        u=U[j];
        grad_u = MATGRAD * U[j];
        BdrySlope(grad_u,1L);
        SlopeLimiterVF(u,grad_u,alpha_u,limiter);
        gU[j].resize(d*g.nt);
        for(long k=0; k< g.nt; ++k)
            for(int c=0;c<d;++c) gU[j][d*k+c] = alpha_u[k] * grad_u[d*k+c];
    }
}

// initialization of parameters
void FreeVol_model::set(long n,long ordd,long numfluxa,const FreeVol_FE &fe, KN<long> labelBCC, const char *equationn, long limiterr, KN<double> parameterss, KN<long> intparameterss)
{
    if (n < 1 || n > Neqmax) ExecError("FreeVol: the number of equations of the model must be between 1 and Neqmax = 100");
    Neq=n;
    geom.build(fe);
    ord=ordd;
    numflux=numfluxa;
    if (!(equation)) equation=new string;
    *equation=string(equationn);
    limiter=limiterr;
    parameters.resize(parameterss.N());
    parameters=parameterss;
    intparameters.resize(intparameterss.N());
    intparameters=intparameterss;
    labelBC.resize(labelBCC.n);
    labelBC=labelBCC;
    matgrad=0;
    if(ord>1) {
        if(geom.fedim==2) matgrad=Matgrads((const Mesh *)geom.pTh,2);// 2D: the gradients of FreeVol (vertices)
        else matgrad=FreeVol_MatgradsLSQ(geom);// 1D, 3D, surface: least squares on the neighbours
    }
    if(verbosity>4)
        cout << "FreeVol_model::set " << *equation << " dim=" << geom.dim << " Neq=" << Neq << " ord=" << ord << " numflux=" << numflux <<" limiter="<< limiter << " nt=" << geom.nt <<endl;
}

// periodic boundary conditions, fvscheme FV(..., periodic=[[l1,f1],[l2,f2], ...]) (the syntax of a fespace):
// the boundary faces of the label l1 and of the label l2 with the same values of the functions at their vertices
// are glued: each one becomes an interior face between its cell and the cell on the other side (adj, adjf,
// flab = 0), so all the models compute the flux once (conservative, also the non conservative fluxes left/right of
// the interior faces), the slopes and the limiter see the neighbours on the other side
// the glued faces must be the translation of each other (pshift: the translation of the neighbour, used by the
// least squares gradients of 1D, 3D); their vertices are in the same periodic class (vper: a corner of a domain
// periodic in x and y is a class of 4 vertices), the gradients at the vertices (2D, order 2) are computed on the
// triangles around all the vertices of the class, translated (Matgrads)
// 1D (meshL): periodic=[[l1],[l2], ...], the faces are the ends (a point of label l1 glued to the point of label l2)
// 2D (mesh): periodic=[[l1,f1],[l2,f2], ...], a vertex identified by the value of f (1 function)
// 3D (mesh3): periodic=[[l1,f1,g1],[l2,f2,g2], ...], a vertex identified by the values of (f, g) (2 functions)
// surface (meshS): periodic=[[l1,f1],[l2,f2], ...] as in 2D (the faces are the edges), the vertices in 3D; the glued
// edges must be the translation of each other (the conormals of the 2 triangles are not opposite on a curved
// surface: the flux uses their mean, as on the interior edges)
// 2D, 3D, surface: 2 faces are glued when a permutation of the vertices of one gives the same values at all the vertices
void FreeVol_model::SetPeriodic(int nbc, Expression *per, Stack s)
{
    FreeVol_Geom &g = geom;
    if (nbc <= 0) return;
    if (g.fedim < 1 || g.fedim > 4)
        ExecError("fvscheme: periodic= only on a meshL (1D), a mesh (2D), a mesh3 (3D) or a meshS (surface)");
    if (!s) ExecError("fvscheme periodic: no stack to compute the functions of periodic=");
    const int d = g.fedim, nvf = g.nvf;
    const int nk = d == 4 ? 1 : d - 1;    // nk: the functions which identify a vertex (surface: 1, as in 2D)
    const int nl = nk + 1;                // the expressions by label: [label, f (,g)]
    MeshPoint *mp = MeshPointStack(s), smp = *mp;
    const long nkf = g.nt * g.nf;
    g.pshift.resize(nkf);
    g.pshift = R3();
    g.plab.resize(2 * nbc);
    KN<long> uf(g.nv);    // union-find of the periodic classes of the vertices
    for (long i = 0; i < g.nv; ++i) uf[i] = i;
    struct UF { static long find(KN<long> &u, long i) { while (u[i] != i) { u[i] = u[u[i]]; i = u[i]; } return i; } };
    // the coordinates of the vertex i
    auto vcoord = [&](long i) -> R3 {
        if (d == 2) { const R2 &P = (*(const Mesh *)g.pTh)(i); return R3(P.x, P.y, 0.); }
        if (d == 4) return (*(const MeshS *)g.pTh)(i);
        if (d == 1) return (*(const MeshL *)g.pTh)(i);
        return (*(const Mesh3 *)g.pTh)(i);
    };
    // glue the face kf1 to the face kf2: the vertex l of kf1 is the vertex p[l] of kf2
    auto glue = [&](long kf1, long kf2, const int *p) {
        const long k1 = kf1 / g.nf, k2 = kf2 / g.nf;
        if (k1 == k2) ExecError("fvscheme periodic: a cell is glued to itself (mesh too coarse)");
        bool trans = fabs(g.fmeas[kf1] - g.fmeas[kf2]) <= 1e-6 * g.fmeas[kf1];
        if (g.surf) {    // surface: the vertices of the 2 edges differ by the same translation
            const R3 T = g.F[kf1] - g.F[kf2];
            for (int l = 0; l < nvf && trans; ++l)
                trans = (vcoord(g.fvert[nvf * kf1 + l]) - vcoord(g.fvert[nvf * kf2 + p[l]]) - T).norme() <= 1e-6 * g.fmeas[kf1];
        }
        else trans = trans && (g.N[kf1] + g.N[kf2]).norme() <= 1e-6;
        if (!trans) ExecError("fvscheme periodic: 2 glued faces must be the translation of each other");
        g.adj[kf1] = k2; g.adjf[kf1] = kf2 % g.nf;
        g.adj[kf2] = k1; g.adjf[kf2] = kf1 % g.nf;
        g.flab[kf1] = g.flab[kf2] = 0;
        g.pshift[kf1] = g.F[kf1] - g.F[kf2];
        g.pshift[kf2] = g.F[kf2] - g.F[kf1];
        g.nper += 2;
        for (int l = 0; l < nvf; ++l) {    // the vertices glued are the same vertex
            long r1 = UF::find(uf, g.fvert[nvf * kf1 + l]), r2 = UF::find(uf, g.fvert[nvf * kf2 + p[l]]);
            if (r1 != r2) uf[r1] = r2;
        }
    };
    // the permutations of the vertices of a face (2 in 2D, 6 in 3D)
    static const int perm2[2][3] = {{0, 1, 2}, {1, 0, 2}};
    static const int perm3[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    const int np = nvf == 3 ? 6 : (nvf == 2 ? 2 : 1);
    const int (*perms)[3] = nvf == 3 ? perm3 : perm2;
    const int ne = 2 * nl;    // the expressions by pair: [label1, f1 (,g1)], [label2, f2 (,g2)]
    // domain decomposition (the local mesh has faces of label -111111, the border of the local mesh): macro_ddm with
    // ThPeriodicity (the labels) puts all the elements near the periodic borders on one process, the other processes
    // have no periodic face; a face whose partner is not in the local mesh becomes a face of the border of the local
    // mesh (-111111: its cell is an overlap cell, set again by the exchange)
    bool ddm = false;
    for (long kf = 0; kf < nkf && !ddm; ++kf) ddm = g.adj[kf] < 0 && g.flab[kf] == -111111;
    long ncut = 0;
    for (int ip = 0; ip < nbc; ++ip) {
        Expression *pp = per + ne * ip;
        long lab[2] = {GetAny< long >((*pp[0])(s)), GetAny< long >((*pp[nl])(s))};
        g.plab[2 * ip] = lab[0]; g.plab[2 * ip + 1] = lab[1];
        vector< long > L[2];       // the boundary faces of the 2 labels
        vector< double > key[2];   // the values of the nk functions at the nvf vertices of each face
        double hmn = 1e300;
        for (int side = 0; side < 2; ++side)
            for (long kf = 0; kf < nkf; ++kf)
                if (g.adj[kf] < 0 && g.flab[kf] == lab[side]) {
                    L[side].push_back(kf);
                    if (!nk) continue;
                    double kv[3][2];
                    for (int l = 0; l < nvf; ++l) {
                        R3 P = vcoord(g.fvert[nvf * kf + l]);
                        mp->set(P.x, P.y, P.z);
                        for (int c = 0; c < nk; ++c) kv[l][c] = GetAny< double >((*pp[side * nl + 1 + c])(s));
                        for (int c = 0; c < nk; ++c) key[side].push_back(kv[l][c]);
                    }
                    for (int l = 0; l < nvf; ++l)    // the vertices of a face must have different values
                        for (int m = l + 1; m < nvf; ++m) {
                            double dlm = 0;
                            for (int c = 0; c < nk; ++c) dlm = max(dlm, fabs(kv[l][c] - kv[m][c]));
                            hmn = min(hmn, dlm);
                        }
                }
        if (!ddm && (L[0].empty() || L[0].size() != L[1].size())) {
            cout << "fvscheme periodic: label " << lab[0] << " : " << L[0].size() << " boundary faces, label " << lab[1]
                 << " : " << L[1].size() << " boundary faces" << endl;
            ExecError("fvscheme periodic: the 2 labels of a pair must have the same number of boundary faces (> 0)");
        }
        const long n = L[1].size();
        if (ddm && (L[0].empty() || L[1].empty())) {    // ddm: no pair in the local mesh, the faces are cut faces
            for (int side = 0; side < 2; ++side) for (size_t q = 0; q < L[side].size(); ++q) { g.flab[L[side][q]] = -111111; ncut++; }
            continue;
        }
        if (!nk) {    // 1D: a point of each label
            if (n != 1 || L[0].size() != 1) ExecError("fvscheme periodic (meshL): one end (one boundary point) by label");
            glue(L[0][0], L[1][0], perm2[0]);
        } else {
            if (hmn < 1e-30) ExecError("fvscheme periodic: the functions of periodic= are the same at 2 vertices of a face (they must identify the points of the 2 borders)");
            const double tol = hmn / 8;    // as the periodic fespaces
            const int nkv = nk * nvf;      // the values of a face
            // the faces of the label 2 sorted by the min of the first function on their vertices
            vector< pair< double, long > > s2(n);
            for (long j = 0; j < n; ++j) {
                double m = 1e300;
                for (int l = 0; l < nvf; ++l) m = min(m, key[1][nkv * j + nk * l]);
                s2[j] = make_pair(m, j);
            }
            sort(s2.begin(), s2.end());
            vector< char > used(n, 0);
            const long n0 = L[0].size();
            for (long i = 0; i < n0; ++i) {
                const double *a = &key[0][nkv * i];
                double amin = 1e300;
                for (int l = 0; l < nvf; ++l) amin = min(amin, a[nk * l]);
                long j = -1;
                int jp = -1;
                for (vector< pair< double, long > >::const_iterator it = lower_bound(s2.begin(), s2.end(), make_pair(amin - tol, -1L));
                     j < 0 && it != s2.end() && it->first <= amin + tol; ++it) {
                    long jj = it->second;
                    if (used[jj]) continue;
                    const double *b = &key[1][nkv * jj];
                    for (int q = 0; q < np && j < 0; ++q) {    // a permutation of the vertices of jj with the same values
                        bool ok = true;
                        for (int l = 0; l < nvf && ok; ++l)
                            for (int c = 0; c < nk && ok; ++c) ok = fabs(a[nk * l + c] - b[nk * perms[q][l] + c]) <= tol;
                        if (ok) { j = jj; jp = q; }
                    }
                }
                const long kf1 = L[0][i];
                if (j < 0 && ddm) { g.flab[kf1] = -111111; ncut++; continue; }    // ddm: its partner is not local
                if (j < 0) {
                    cout << "fvscheme periodic: the face of label " << lab[0] << " at " << g.F[kf1].x << " " << g.F[kf1].y << " "
                         << g.F[kf1].z << " has no face of label " << lab[1] << " with the same values" << endl;
                    ExecError("fvscheme periodic: a face has no corresponding face");
                }
                used[j] = 1;
                glue(kf1, L[1][j], perms[jp]);
            }
            if (ddm) for (long j = 0; j < n; ++j) if (!used[j]) { g.flab[L[1][j]] = -111111; ncut++; }
        }
        if (verbosity > 1 && !mpirank)
            cout << " fvscheme periodic: labels " << lab[0] << " and " << lab[1] << " : " << n << " faces glued" << endl;
    }
    if (ddm && ncut && verbosity > 1) cout << " fvscheme periodic (mpirank " << mpirank << "): " << ncut << " periodic faces without partner in the local mesh: border of the local mesh (-111111)" << endl;
    *mp = smp;
    g.vper.resize(g.nv);
    for (long i = 0; i < g.nv; ++i) g.vper[i] = UF::find(uf, i);
    FreeVol_GeomCount(g);
    if (ord > 1) {    // the gradients with the neighbours on the other side
        if (matgrad) delete matgrad;
        if (d == 2) matgrad = Matgrads((const Mesh *)g.pTh, 2, &g.vper);
        else matgrad = FreeVol_MatgradsLSQ(g);    // 1D, 3D, surface: least squares, the neighbours translated (pshift)
    }
}
 FreeVol_model::~FreeVol_model()
{
    Neq=-1;
    ord=0;
    numflux=-1;
    limiter=0;
    if(matgrad) {delete matgrad; matgrad=0;}
    if(equation) {delete equation; equation=0;}
    if ((!mpirank)&&(verbosity>99)) cout << " destroy FreeVol_model " << this << endl;
}

// ============================================================================
// the arrays of the scripts: FE functions in 1D, 2D or 3D
// ============================================================================
// the dimension of an expression of the script: 2 (FE on a mesh), 3 (mesh3), 1 (meshL), 4 (meshS), 0 otherwise
static int FreeVol_DimOf(const C_F0 &e)
{
    if (BCastTo< pfer >(e)) return 2;
    if (BCastTo< pf3r >(e)) return 3;
    if (BCastTo< pfLr >(e)) return 1;
    if (BCastTo< pfSr >(e)) return 4;
    return 0;
}
static Expression FreeVol_CastFE(const C_F0 &e, int dim)
{
    if (dim == 2) return atype< pfer >()->CastTo(e);
    if (dim == 3) return atype< pf3r >()->CastTo(e);
    if (dim == 4) return atype< pfSr >()->CastTo(e);
    return atype< pfLr >()->CastTo(e);
}
static FreeVol_FE FreeVol_EvalFE(Expression e, int dim, Stack s)
{
    FreeVol_FE fe;
    fe.dim = dim;
    if (dim == 2) fe.f2 = GetAny< pfer >((*e)(s));
    else if (dim == 3) fe.f3 = GetAny< pf3r >((*e)(s));
    else if (dim == 4) fe.fS = GetAny< pfSr >((*e)(s));
    else fe.fL = GetAny< pfLr >((*e)(s));
    return fe;
}

// ============================================================================
// ========== OneOperatorFVMOp for computations as "fluxupdate"     ===========
// ============================================================================
template <int NU>
class  OneOperatorFVMOp : public OneOperator {
    public:
    int cas;
    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =0;
    E_F0 * code(const basicAC_F0 & args) const {
        if((verbosity>9)&&(!mpirank)) cout <<" OneOperatorFVMOp "<<n_name_param << endl;
        Expression *nargs = new Expression[n_name_param];
        args.SetNameParam(n_name_param,name_param,nargs);
        array<const E_Array *,NU> peaU;
        for(int i=0; i<NU;++i){
            int j = args.size()-NU+i;
            peaU[i] =  dynamic_cast<const E_Array *>(args[j].LeftValue());
            ffassert(peaU[i]);
            ffassert(peaU[0]->size()==peaU[i]->size());
        }
        return new Op(cas,nargs,args[0].LeftValue(),peaU);
    }

    OneOperatorFVMOp():// CFLUX (FVM,[u],[up],[uBC]);//  NU=3
    OneOperator(atype<R>(),atype<FreeVol *>(), atype<E_Array>(),atype<E_Array>(),atype<E_Array>()),cas(0){ ffassert(NU==3);}

    class Op : public E_F0{
        public:
        Expression *nargs;
        int cas;
        int dim;    // dimension of the arrays: 1 (meshL), 2 (mesh), 3 (mesh3), 4 (meshS)
        Expression eFVM;
        array<const E_Array *,NU> peaU;
        array<vector<Expression>,NU>  eU;

        Op(int kas,Expression *nnargs,Expression FVM, array<const E_Array *,NU> & eaU)
        : cas(kas),nargs(nnargs),eFVM(FVM),peaU(eaU){
            int Neq = peaU[0]->size();
            dim = FreeVol_DimOf((*peaU[0])[0]);
            if (!dim) CompileError("fluxupdate: the arrays must be FE functions on a meshL, a mesh, a mesh3 or a meshS");
            for(int i=0;i<Neq;++i){
                for(int j=0; j<NU;++j){
                    if (FreeVol_DimOf((*peaU[j])[i]) != dim)
                        CompileError("fluxupdate: all the arrays must be FE functions on meshes of the same dimension");
                    eU[j].push_back( FreeVol_CastFE((*peaU[j])[i], dim));
                }
            }
        }
        AnyType operator()(Stack s)  const {
            FreeVol * pVF = GetAny< FreeVol*>((*eFVM)(s));
            // a fvscheme declared without definition: fvscheme FV; then fluxupdate(FV, ...)
            if (!pVF || !pVF->FVM)
                ExecError("fluxupdate: fvscheme not initialized, use fvscheme FV([V1,V2], eq=...) or resetfvscheme");
            if (pVF->FVM->geom.fedim != dim)
                ExecError("fluxupdate: the arrays are not of the dimension of the fvscheme (meshL 1D, mesh 2D, mesh3 3D, meshS surface)");
            {   // the mesh of the unknowns: a new mesh since the definition (adaptmesh, ...)
                FreeVol_FE u0 = FreeVol_EvalFE(eU[0][0], dim, s);
                if (u0.x() && u0.mesh() != pVF->FVM->geom.pTh)
                    FreeVol_MeshChanged(pVF, u0.mesh(), s);
            }
            FreeVol_model *pFVM = pVF->FVM;
            int Neq = pFVM->Neq;
            ffassert( Neq == peaU[0]->size() );
            array<KN<KN_d>,NU> U;
            for(int j=0; j<NU;++j)
              U[j].resize(Neq);
            long dimcomp=0;
            for(int i=0; i<Neq;++i)
            {
                for(int j=0; j<NU;++j){
                    FreeVol_FE uij = FreeVol_EvalFE(eU[j][i], dim, s);
                    if(uij.x() ==0) {
                        cout << " Unset FE "<< j << " maybe you need to initialize ur" << endl;
                        ffassert(0);
                    }
                    if (uij.mesh() != pFVM->geom.pTh)
                        ExecError("fluxupdate: all the arrays ([u], [ur], [uBC]) must be on the mesh of the unknowns (after adaptmesh: ur = 0; uBC = ...;)");
                    if (dimcomp==0) dimcomp=uij.NbOfNodes();
                    else ffassert( uij.NbOfNodes() == dimcomp);
                    //dimcomp=nt for P0
                    U[j][i].set( (*uij.x())(SubArray(dimcomp,uij.comp(),uij.N())));
                }
            }
            // U[0] vector of unknowns (input)
            // U[1] vector of increments of unknows (in fact num fluxes) to be computed
            // U[2] vector of Dirichlet boundary conditions
            return pFVM->computedtFlux(U[0],U[1],U[2]);
        }

        operator aType () const { return atype<R>();}

    };
};

template<int NU>
basicAC_F0::name_and_type  OneOperatorFVMOp<NU>::name_param[]= {};

// ============================================================================
// ========== Registry of the models (eq="advection", ...)          ===========
// ============================================================================
// the table is in a function: built at its first use, whatever the order of the initializations
static map< string, FreeVol_Creator > &FreeVol_Registry() {
    static map< string, FreeVol_Creator > reg;
    return reg;
}

void FreeVol_RegisterModel(const string &eq, FreeVol_Creator creator) {
    FreeVol_Registry()[eq] = creator;
    if((verbosity>1)&&(!mpirank)) cout << " FreeVol: model \"" << eq << "\" registered" << endl;
}

// the plugin of the model eq is "FreeVol_" + eq (ex. FreeVol_advection for eq="advection"):
// loaded at the first use of the equation, so load "FreeVol" is enough
bool load(string ss);    // the load of FreeFEM (load.cpp)

FreeVol *FreeVol_CreateModel(FreeVol *a, const string &eq, const vector< FreeVol_FE > &V,
                             long order, long numflux, KN<long> labelBC, long limiter,
                             KN<double> parameters, KN<long> intparameters, bool rebuild,
                             const void *dual, int nbcperiodic, Expression *periodic, Stack s) {
    map< string, FreeVol_Creator >::const_iterator m = FreeVol_Registry().find(eq);
    if (m == FreeVol_Registry().end()) {
        // not a load during the initialization of a plugin (list of the init functions in use)
        load("FreeVol_" + eq);
        m = FreeVol_Registry().find(eq);
    }
    if (m == FreeVol_Registry().end()) {
        cout << "FreeVol: unknown equation eq=\"" << eq << "\" (plugin FreeVol_" << eq << " not found); known equations:";
        for (m = FreeVol_Registry().begin(); m != FreeVol_Registry().end(); ++m) cout << " \"" << m->first << "\"";
        cout << endl;
        ExecError("FreeVol: unknown equation");
    }
    FreeVol_model *fvm = m->second(V, order, numflux, labelBC, limiter, parameters, intparameters, dual);
    if (nbcperiodic) {
        if (!fvm->PeriodicOK()) { delete fvm; ExecError("fvscheme: periodic= not yet with this model (dual mesh: advectiondm, igddm, swtopodm, eulerdm)"); }
        fvm->SetPeriodic(nbcperiodic, periodic, s);
    }
    FreeVol_Params *par = new FreeVol_Params;
    par->eq = eq;
    par->V = V;
    par->order = order; par->numflux = numflux; par->limiter = limiter;
    par->labelBC.resize(labelBC.N()); par->labelBC = labelBC;
    par->parameters.resize(parameters.N()); par->parameters = parameters;
    par->intparameters.resize(intparameters.N()); par->intparameters = intparameters;
    par->rebuild = rebuild;
    par->dual = dual;
    par->nbcperiodic = nbcperiodic;
    par->periodic = periodic;
    a->destroy();    // a fvscheme can be defined again (resetfvscheme)
    a->FVM = fvm;
    a->par = par;
    a->pThref = V[0].meshrc();
    a->pThref->increment();
    return a;
}

void FreeVol_MeshChanged(FreeVol *a, const void *pThnew, Stack s) {
    if (!a->par->rebuild)
        ExecError("fluxupdate: the mesh has changed since the definition of the fvscheme (adaptmesh?), call resetfvscheme(FV, [V1,V2], eq=..., ...) after the change of mesh (or rebuild=true)");
    // the velocity must be on the new mesh (recomputed after adaptmesh: V1 = ...)
    for (size_t i = 0; i < a->par->V.size(); ++i)
        if (!a->par->V[i].x() || a->par->V[i].mesh() != pThnew)
            ExecError("fluxupdate: the mesh has changed, but the velocity of the fvscheme is not on the new mesh: compute it again after adaptmesh (ex. V1 = ...; V2 = ...;)");
    FreeVol_Params p = *a->par;    // a copy: CreateModel deletes a->par
    FreeVol_CreateModel(a, p.eq, p.V, p.order, p.numflux, p.labelBC, p.limiter,
                        p.parameters, p.intparameters, p.rebuild, p.dual, p.nbcperiodic, p.periodic, s);
    if ((verbosity) && (!mpirank)) cout << " fvscheme rebuilt on the new mesh (nt = " << a->FVM->geom.nt << ")" << endl;
}

// ============================================================================
// ========== OneOperatorFreeVol (for FV <- [V1, V2, ...]) ===========
// ============================================================================
class  OneOperatorFreeVol : public OneOperator {
    // FreeVol *, FreeVol *,E_Array*
public:
    typedef FreeVol * R;
    typedef FreeVol * A;
    typedef E_Array  B;

    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =9;
    aType t0, t1;
    class CODE: public E_F0 {
    public:
        static const int NMax = 100;
        int N ;
        int dim;    // the dimension, given by the FE functions: 1 (meshL), 2 (mesh), 3 (mesh3)
        Expression a0,pfers[NMax];
        Expression edual;    // the dual mesh, fvscheme FV(Thstar, [V1,V2], ...), 0 if not given
        Expression *nargs;
        int nbcperiodic;     // periodic=[[l1,f1],[l2,f2], ...]: the number of pairs, and their expressions
        Expression *periodic;
        long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
        KN<double>* arg(int i,Stack stack,KN<double>* a) const{ return nargs[i] ? GetAny<KN<double>*>( (*nargs[i])(stack) ): a;}
        string * arg(int i,Stack stack,string * a) const{ return nargs[i] ? GetAny<string *>( (*nargs[i])(stack) ): a;}
        KN<long>* arg(int i,Stack stack,KN<long>* a) const{ return nargs[i] ? GetAny<KN<long>*>( (*nargs[i])(stack) ): a;}

        CODE(Expression* nnargs,Expression aa0, const B* ppfers, Expression eedual=0) {
            nargs=nnargs;
            edual=eedual;
            N =ppfers->nbitem();
            a0=aa0;
            ffassert(N<100);// less the 100 equation!!!
            dim = FreeVol_DimOf((*ppfers)[0]);
            for(int i=0; i< N;++i)
                pfers[i]=FreeVol_CastFE((*ppfers)[i], dim);
            // periodic= : as for a fespace, [label, f] in 2D (and on a surface), [label, f1, f2] in 3D, [label] in 1D
            nbcperiodic = 0; periodic = 0;
            GetPeriodic(dim == 1 ? 1 : (dim == 3 ? 3 : 2), nargs[8], nbcperiodic, periodic);
        }
        AnyType operator()(Stack s)  const
        {
            A a =GetAny<A>((*a0)(s));
            long order=arg(0,s,2L);//2 by default
            string *eq=arg(1,s,(string*) 0);//empty string by default
            KN<long> labelBCdef(1);
            labelBCdef[0] = 1;
            KN<long> labelBC = *(arg(2,s,&labelBCdef));//1 by default
            long numflux = arg(3,s,-99L);//-99 means undefined, default value has to be defined in the FreeVol_model separate plugin file
            long limiter = arg(4,s,1L);//1 by default
            KN<double> parametersdef(1);
            parametersdef[0]=0.;
            KN<double> parameters=*(arg(5,s,&parametersdef));//0. by default
            KN<long> intparametersdef(1);
            intparametersdef[0]=0;
            KN<long> intparameters=*(arg(6,s,&intparametersdef));//0 by default
            bool rebuild = nargs[7] ? GetAny<bool>((*nargs[7])(s)) : true;    // after a change of mesh
            if ((!eq)) {cout << "error: equation type needs to be defined by eq=" << endl;ffassert(0);}
            if((verbosity)&&(!mpirank)){
                cout << "FreeVol ";
                cout << "eq : " << *eq;
                cout << ", order " << order << ", dimension " << dim <<endl;
            }
            // the FE functions given to fvscheme (the velocity for the advection)
            vector< FreeVol_FE > V(N);
            for(int i=0; i<N; ++i) V[i] = FreeVol_EvalFE(pfers[i], dim, s);
            // the dual mesh (its type is known by the dimension, checked at the compilation)
            const void *dual = 0;
            if (edual) {
                if (dim == 1) dual = GetAny<dualmeshL *>((*edual)(s));
                else if (dim == 2) dual = GetAny<dualmesh *>((*edual)(s));
                else if (dim == 4) dual = GetAny<dualmeshS *>((*edual)(s));
                else dual = GetAny<dualmesh3 *>((*edual)(s));
                if (!dual) ExecError("fvscheme: the dual mesh is not defined");
            }
            a = FreeVol_CreateModel(a,*eq,V,order,numflux,labelBC,limiter,parameters,intparameters,rebuild,dual,
                                    nbcperiodic,periodic,s);
            return a;
        }
        aType t0,t1; //  return type  type of f,  f(t1, t2)

    }; // End of class CODE
    E_F0 * code(const basicAC_F0 & args) const
    {
        if((verbosity>5)&&(!mpirank)) cout <<" OneOperatorFreeVol "<<n_name_param << endl;
        Expression *nargs = new Expression[n_name_param];   // an array by instruction
        args.SetNameParam(n_name_param,name_param,nargs);
        const B* b= dynamic_cast<const B*>(args[1].LeftValue());// Array of FE functions
        int N = b->nbitem();
        int dim = N ? FreeVol_DimOf((*b)[0]) : 0;
        if (!dim) CompileError("fvscheme: [V1, ...] must be FE functions on a meshL (1D), a mesh (2D), a mesh3 (3D) or a meshS (surface)");
        for(int i=0; i< N;++i)
            if (FreeVol_DimOf((*b)[i]) != dim)
                CompileError("fvscheme: all the FE functions of [V1, ...] must be on meshes of the same dimension");

        return  new CODE(nargs,t0->CastTo(args[0]),b);}

    OneOperatorFreeVol():
      OneOperator(map_type[typeid(R).name()],map_type[typeid(A).name()],map_type[typeid(B).name()]),
      t0( map_type[typeid(A).name()] ),t1(map_type[typeid(B).name()] ) {}

};
basicAC_F0::name_and_type  OneOperatorFreeVol::name_param[]= {
    {  "order", &typeid(long) },
    {  "eq", &typeid(string *) },
    {  "labelBC", &typeid(KN<long>*) },
    {  "numflux", &typeid(long) },
    {  "limiter", &typeid(long) },
    {  "parameters", &typeid(KN<double>*) },
    {  "intparameters", &typeid(KN<long>*) },
    {  "rebuild", &typeid(bool) },    // after a change of mesh: build the model again (default) or stop
    {  "periodic", &typeid(E_Array) }    // periodic boundary conditions (the syntax of a fespace)
};

// fvscheme FV(Thstar, [V1, V2, ...], eq=..., ...): with the dual mesh Thstar built in the script (plugin fvdualmesh:
// meshdualL Thstar = builddualmesh(ThL) in 1D, meshdual in 2D, meshdual3 in 3D), given to the model
// the types of the plugin are not known here (only their names): an operator with an ellipsis, of lower
// priority than OneOperatorFreeVol (FV <- [V1, ...]), the arguments checked at the compilation
class  OneOperatorFreeVolDual : public OneOperator {
public:
    typedef FreeVol * R;
    typedef FreeVol * A;
    typedef E_Array  B;
    aType t0;
    // the dimension of the dual mesh of the type t (1 meshdualL, 2 meshdual, 3 meshdual3, 4 meshdualS: as
    // FreeVol_DimOf), 0 otherwise
    static int DualDim(aType t) {
        const char *nm[4] = {typeid(dualmeshL *).name(), typeid(dualmesh *).name(), typeid(dualmesh3 *).name(), typeid(dualmeshS *).name()};
        for (int d = 0; d < 4; ++d) {
            map<const string, basicForEachType *>::const_iterator i = map_type.find(nm[d]);
            if (i != map_type.end() && i->second == t) return d+1;
        }
        return 0;
    }
    E_F0 * code(const basicAC_F0 & args) const
    {
        const B* b = args.size()==3 ? dynamic_cast<const B*>(args[2].LeftValue()) : 0;// Array of FE functions
        if (!b) CompileError("fvscheme: fvscheme FV([V1, ...], ...) or fvscheme FV(Thstar, [V1, ...], ...)");
        int ddim = DualDim(args[1].left());
        if (!ddim) CompileError("fvscheme FV(Thstar, [V1, ...], ...): Thstar must be a dual mesh (meshdualL, meshdual, meshdual3 or meshdualS, builddualmesh)");
        int N = b->nbitem();
        int dim = N ? FreeVol_DimOf((*b)[0]) : 0;
        if (!dim) CompileError("fvscheme: [V1, ...] must be FE functions on a meshL (1D), a mesh (2D), a mesh3 (3D) or a meshS (surface)");
        for(int i=0; i< N;++i)
            if (FreeVol_DimOf((*b)[i]) != dim)
                CompileError("fvscheme: all the FE functions of [V1, ...] must be on meshes of the same dimension");
        if (ddim != dim) CompileError("fvscheme FV(Thstar, [V1, ...], ...): the dual mesh is not of the dimension of [V1, ...] (meshdualL 1D, meshdual 2D, meshdual3 3D, meshdualS surface)");
        Expression *nargs = new Expression[OneOperatorFreeVol::n_name_param];   // an array by instruction
        args.SetNameParam(OneOperatorFreeVol::n_name_param,OneOperatorFreeVol::name_param,nargs);
        return new OneOperatorFreeVol::CODE(nargs,t0->CastTo(args[0]),b,args[1].left()->CastTo(args[1]));
    }
    OneOperatorFreeVolDual():
      OneOperator(map_type[typeid(R).name()],ArrayOfaType(map_type[typeid(A).name()],true)),
      t0( map_type[typeid(A).name()] ) { pref = -1;}    // FV <- [V1, ...] first
};

// ============================================================================
// ========== registration of the operators in FreeFEM ========================
// ============================================================================
// in the kernel of FreeFEM (no load "FreeVol"): called at the start by msh3_Load_Init (msh3.cpp), after
// Load_Init_fvdualmesh (fvscheme FV(Thstar, ...) uses the types of the dual meshes); the models (FreeVol_<eq>)
// are plugins of plugin/seq, loaded by fvscheme at their first use
void Load_Init_FreeVol() {
    if (verbosity && !mpirank) cout << "FreeVol ";// a kernel module: its name in the line " Load: lg_fem ..." (lg.ypp)

    // 1. utility functions (2D): MatVFD, MatVFM01, OrderVF
    Global.Add("MatVFD", "(", new OneOperator1s_<newpMatrice_Creuse<R> , const Mesh *>  (Mat_VFD));
    Global.Add("MatVFD", "(", new OneOperator2s_<newpMatrice_Creuse<R> , const Mesh *,long>  (Mat_VFD));
    Global.Add("MatVFM01", "(", new OneOperator1s_<newpMatrice_Creuse<R> , const Mesh *>  (Mat_VFM01));
    Global.Add("OrderVF", "(", new OneOperator3_<KN<double>* , const Mesh *,KN<double> *,long >  (OrderVF));

    // 2. the type FreeVol*, the keyword fvscheme
    Dcl_Type< FreeVol * >(InitP< FreeVol >, Destroy< FreeVol >);
    zzzfff->Add("fvscheme",atype<FreeVol*>());
    // 3. definition (fvscheme FV([V1,...], ...)) and new definition (resetfvscheme)
    TheOperators->Add("<-", new OneOperatorFreeVol, new OneOperatorFreeVolDual);
    Global.Add("resetfvscheme","(",new OneOperatorFreeVol, new OneOperatorFreeVolDual);

    // 4. the computations
    Global.Add("fluxupdate","(",new OneOperatorFVMOp<3>());

}
