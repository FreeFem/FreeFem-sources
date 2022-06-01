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

class TypeOfFE_Pk_L : public GTypeOfFE< MeshL > {
public:
    typedef MeshL Mesh;
    typedef MeshL::Element Element;
    typedef GFElement< MeshL > FElement;
    typedef R1 RdHat;
    typedef R3 Rd;
    static const int dHat = 1;
    static const int d = 3;
    const int kp ;
    KN<int> mi ;
    KN<double> ml,mc;
    
    struct A4 {
        int dfon[4];
        A4(int k) {
            dfon[0] = dfon[1] = dfon[2] = dfon[3] = 0;
            ffassert(k>0);
            dfon[0] = 1;
            dfon[1] = k-1;
        }
        operator const int * () const {return dfon;}
    };
    TypeOfFE_Pk_L( int kk)
    :  GTypeOfFE< Mesh >(A4(kk), 1, kk, true, false), kp(kk), mi((kk+1)*kk), ml((kk+1)*kk), mc((kk+1)*kk)
    {
        typedef Element E;
        int n = this->NbDoF;
        bool dd = verbosity > 5;
        if (dd) {
            cout << "\n +++ P"<< kp <<" L  : ndof : " << n << " " << this->PtInterpolation.N( ) << endl;
        }
        ffassert( n == kp+1);
        ffassert(this->PtInterpolation.N( )==kp+1);
        RdHat *Pt = this->PtInterpolation;
        KN<int> pdof(n);
        // construction of interpolation ppoint
        {
            double cc = 1. / kp;
            Pt[0] = R1(0);
            Pt[1] = R1(1.);
            pdof[0] = 0;
            pdof[1]= kp;
            
            for(int i=1; i<kp; ++i)
            {
                Pt[i+1] = R1(cc*i);
                pdof[i+1]= i;
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
        // constructon de
        //                     //  (k*l[li]-i)/(ii-i);;
        int km =0;
        for(int dof=0; dof<n;++dof)
        {
            int ii = pdof[dof]; // point associaed to dof
            // cout << " dof " << dof << "km=" <<km << " ii = " << ii << endl;
            
            for(int i=0; i<= kp; ++i)
            {
                if( i != ii)
                {
                    mi[km]=1;// x
                    ml[km] = double(kp)/(ii-i);
                    mc[km] = -double(i)/(ii-i);
                    km++;
                }
            }
        }
        if(verbosity>9)
        {
            cout << " km == " << km << "kp = " << kp << " " << n << endl;
            
            for(int dof=0,m=0; dof<n;++dof)
            {
                cout.setf(ios::showpos);
                cout << dof << " : " ;
                for(int l=0; l<kp;++l,++m)
                cout << "("<< ml[m]<< " l_" << mi[m] << " +"<<  mc[m] <<")" ;
                cout << endl;
                cout.unsetf(ios::showpos);
            }
        }
        ffassert(km == kp*n);
        
        
        
        
    }
    void FB(const What_d whatd, const Mesh &Th, const Element &K, const RdHat &PHat, RNMK_ &val) const
    {
        const int k =kp;
        double l[dHat+1];
        PHat.toBary(l);
        // l = 0 en i/k , 1 en ii/k  => monome : ( k*l - i )/( ii-i) ok..
        RN_ f0(val('.',0,op_id));
        if (whatd & Fop_D0)
        {  double s=0.;
            for( int dof = 0,m=0; dof < this->NbDoF; ++dof)
            {
                double f=1.;
                for(int i=0; i<k;++i,++m)
                f *= (ml[m]*l[mi[m]]+mc[m]);
                f0[dof]= f;
                s+= f;
            }
            
            ffassert( abs( s-1.)< 1e-7);
        }
        if(whatd & (Fop_D1|Fop_D2))
        {
            bool d2= Fop_D2 && k>1;//  calcul de
            const unsigned int fop[3]={Fop_dx,Fop_dy,Fop_dz};
            const  int op[3]={op_dx,op_dy,op_dz};
            const  int dop[9]={op_dxx,op_dxy,op_dxz, op_dyx,op_dyy,op_dyz, op_dzx,op_dzy,op_dzz};
            
            Rd Dl[dHat+1];
            Rd DDl[dHat+1][d];
            K.Gradlambda(Dl);
            KN<Rd> df(this->NbDoF),ddf(this->NbDoF*d);
            
            for( int dof = 0,m=0; dof < this->NbDoF; ++dof)
            {
                double f=1.;
                Rd Df;
                Rd DDf[d] ;
                for(int i=0; i<k;++i,++m)
                { // f = f*b // df = df*b+db*f; ddf  = ddf*b + 2*df*db + f * ddb
                    int im=mi[m];
                    double b =(ml[m]*l[im]+mc[m]);
                    Rd Db=ml[m]*Dl[im];
                    if(d2)
                        for(int l=0; l<d;++l)
                          DDf[l] = b*DDf[l]+ Db[l]*Df + Db*Df[l] ;
                    Df = b*Df+ (f*Db);
                    f *= b;
                }
                df[dof] = Df;
                if(d2)
                {
                    for(int l=0; l<d;++l)
                      ddf[dof*3+l] = DDf[l];
                }
            }
            // copy data d  D
            for (int dd=0; dd< Rd::d;++dd)
            {
                if (whatd & fop[dd])
                {
                    RN_ dfdd(val('.',0,op[dd]));
                    for(int i=0;i<this->NbDoF;++i)
                      dfdd[i]= df[i][dd];
                }
            }
            // copy data  DD
            if(d2)
            {
                for (int id=0; id< Rd::d;++id)
                for (int jd=0; jd<= id;++jd)
                {
                    int op = dop[id*3+jd];
                    const unsigned int fopij= 1<< op;
                    if (whatd & fopij)
                    {
                        RN_ dfdd(val('.',0,op));
                        for(int i=0;i<this->NbDoF;++i)
                        dfdd[i]= ddf[i*3+id][jd];
                    }
                }
            }
            
        }
    }
};
