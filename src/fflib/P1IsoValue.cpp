//
//  P1IsoValue.cpp
//  ff
//
//  Created by Frédéric Hecht on 07/03/2014.
//
//
#include <cassert>
#include <cmath>
using namespace std;
#include "P1IsoValue.hpp"
using namespace Fem2D;

typedef double R;

inline R3 bary(const R3 K[4],R f[4],int i0,int i1,R v)
{
    R d=f[i0]-f[i1];
    assert(fabs(d)>1e-20);
    R l1= (f[i0] - v)/ d;  //  == 1 si v = f[i1]
    R l0 = 1. -l1;
    assert(l0 >=-1e-10 && l1 >= -1e-10);
    return K[i0]*l0 + K[i1]*l1; // == K[i1] si l1 ==1 => v = f[i1]
}

static inline int  signep4(int i0,int i1,int i2,int i3)
{ // calcul du signe dans la permutation
    int s =1;
    if(i0>i1) s=-s,Exchange(i0,i1);
    if(i1>i2) s=-s,Exchange(i1,i2);
    if(i2>i3) s=-s,Exchange(i2,i3); // i3 max
    if(i0>i1) s=-s,Exchange(i0,i1);
    if(i1>i2) s=-s,Exchange(i1,i2); // i2 max < i
    if(i0>i1) s=-s,Exchange(i0,i1);
    return s;
}

int IsoLineK(double *f,Fem2D::R3 *Q,const double eps)
{
    
    static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;
    
    const R3 *K=R3::KHat;
    int kv=0,vk[4],tv[4],kf;
    for(int i=0;i<4;++i)
    {
      if( abs(f[i]) <= eps)
        {
            tv[kv++]=i;
            vk[i]=1;
        }
        else
         vk[kf=i]=0;
    }
    if(kv==3)
    {
        // cout << " full face " << kf << endl;
        // a face complete kf..
      //  warning
        int i1=1,i2=2;
        if(f[kf] <0) i1=2,i2=1;
        Q[0]=K[nvfaceTet[kf][0]];
        Q[1]=K[nvfaceTet[kf][i1]];
        Q[2]=K[nvfaceTet[kf][i2]];
        
        return (f[kf] >0) ? 3:0;// to take one fulL face not to times ...
    }

    R v=0;
    int nP=0;
    int np[4],nm[4],nps[4],nms[4];;
    int km=0,kp=0,kms=0,kps=0;
    
    for (int i=0;i<4;++i)
    {
        if(f[i]<=v+eps) nm[km++]=i;
        if(f[i]>=v-eps) np[kp++]=i;
        // strict ..
        if(f[i]<v-eps) nms[kms++]=i;
        if(f[i]>v+eps) nps[kps++]=i;
    }
    
    //  cout << "IsoLineK: km kp "<< km << " " << kp << endl;
    int h=-1,b[3];
    if(kps==1 && km==3)
    {
        h = nps[0];
        b[0]=nvfaceTet[h][0];
        b[1]=nvfaceTet[h][1];
        b[2]=nvfaceTet[h][2];
    }
    if(kms==1 && kp == 3)
    {
        h = nms[0];
        b[0]=nvfaceTet[h][0];
        b[2]=nvfaceTet[h][1];
        b[1]=nvfaceTet[h][2];
    }
    if(kp==2 && km==2)
    {//  cas quad
        if(signep4(nm[0],nm[1],np[0],np[1]) < 0)
            Exchange(nm[0],nm[1]);
        Q[0]=bary(K,f,nm[0],np[0],v);
        Q[1]=bary(K,f,nm[0],np[1],v);
        Q[2]=bary(K,f,nm[1],np[1],v);
        Q[3]=bary(K,f,nm[1],np[0],v);
        nP=4;
    }
    else if (h>=0)
    { // cas triangle
        Q[0]=bary(K,f,h,b[0],v);
        Q[1]=bary(K,f,h,b[1],v);
        Q[2]=bary(K,f,h,b[2],v);
        nP=3;
    }
    
    return nP;
}

int IsoLineK(double *f,Fem2D::R2 *Q,const double eps)
{
    int debug=0;
    R2 P[3]={ R2(0.,0.),R2(1.,0.),R2(0.,1.)};
    int kv=0,ke=0,e=3;
    int tv[3],te[3],vk[3],i0[3],i1[3];
    for(int i=0;i<3;++i)
    {
        if( abs(f[i]) <= eps) {
            e -= tv[kv++]=i;
            vk[i]=1;
        }
        else
            vk[i]=0;
    }
    if(debug) cout << " ** " <<     kv << endl;
    if(kv>1) //  on 2  vertex on the isoline ....
    {
        if(kv==2)
        {
            if(f[e] > 0.)
            {
                int j0=(e+1)%3;
                int j1=(e+2)%3;
                te[ke]=e+3,i0[ke]=j0,i1[ke]=j0,++ke;
                te[ke]=e,i0[ke]=j1,i1[ke]=j1,++ke;
                // pb d'unicity, need to see the adj triangle ...
                //return 10+e ; // edge number + 10
            }
            else return 0; // skip edge ...
            
        }
        else return 0; //  const funct...
    }
    else // see internal edge ..
        for(int e=0;e<3;++e)
        {
            int j0=(e+1)%3;
            int j1=(e+2)%3;
            if( vk[j0]) //  the intial  point on iso line
            {
                if(0. < f[j1])
                    te[ke]=e,i0[ke]=j0,i1[ke]=j0,++ke;
                else
                    te[ke]=e+3,i0[ke]=j0,i1[ke]=j0,++ke;
            }
            else if (vk[j1]); // skip the final point on iso line
            else if( f[j0] < 0. && 0. < f[j1])  // good  sens
                te[ke]=e,i0[ke]=j0,i1[ke]=j1,++ke;
            else if ( f[j0] > 0. && 0. > f[j1]) // inverse  sens
                te[ke]=e+3,i0[ke]=j1,i1[ke]=j0,++ke;
        }
    if( ke==2)
    {
        // the  K[i1[0]] , Q[0], Q[1] must be direct ...
        // the  K[i0[1]] , Q[0], Q[1] must be direct ...
        // Warning   no trivail case ..  make a plot to see
        //  with is good
        // the first edge must be
        
        if(te[0]<3)  // oriente the line
        {
            assert(te[1] >=3);
            std::swap(te[0],te[1]);
            std::swap(i0[0],i0[1]);
            std::swap(i1[0],i1[1]);
            if(debug) cout << " swap " << endl;
        }
        for(int i=0;i<2;++i)
        {
            int j0=i0[i],j1=i1[i];
            if( j0== j1)
                Q[i] = P[j0];
            else
                Q[i] = (P[j0]*(f[j1]) -  P[j1]*(f[j0]) ) /(f[j1]-f[j0]);
            if(debug) cout << i << " " << j0 << " " << j1 << " : "
                << Q[i] << "***" << endl;
        }
        if(!vk[i1[0]])
            assert( det(P[i1[0]],Q[0],Q[1]) > 0);
        if(!vk[i0[1]])
            assert( det(P[i0[1]],Q[1],Q[0]) > 0);
        return 2;
    }
    // remark, the left of the line is upper .
    return 0;
}
