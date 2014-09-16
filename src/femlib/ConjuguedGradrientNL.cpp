// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:    29 fev 2000  
// -*- Mode : c++ -*-
//
// SUMMARY  : array modelisation 
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Frederic Hecht
// E-MAIL   : frederic.hecht@ann.jussieu.fr
//

/*
 
 
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 
 
 */

template<class R,class DJ,class M>
R argmin(R  rho,const M & dJ,const KN_<R> &x,KN_<R> &h,KN_<R> &g,KN_<R> &w)
{
  //  find  ro / (dJ(x+roh),h)  =0 
  // remark input: dJ(x)=g 
  int k=0;
  
  R ro0=0, ro=rho,ro1=rho;
  R p0= (g,h),p,p1;
  R ap0=Abs(p0)*0.1; // on arret quand on a divise par 10.
  x += (ro-rold)* h; rold=ro; g =A(x);
  p= ( p1 = (g,h) );
  

  bool loop=true;
  while (k++<100 && loop)
  { //  calcul du  nouveau ro courant
    if (p0*p1 <0) { //  Ok changement de signe 
      ro = (p0*ro1+p1*ro0)/ (p0+p1); 
      x += (ro-rold)* h; rold=ro; g =A(x);
      p = (g,h);
      if (  verbosity >3 )      
      cout << "         ro " << ro << " gh= " << p
           << "; ro0, gh0 = " << ro0 << " " << p0 
           << "; ro1, gh1 = " << ro1 << " " << p1 << endl; 
      
      if(Abs(p) <= ap0 || k>10  ) return ro; 
      if(p0*p<0) { p1=p;ro1=ro;} 
      else {p0=p;ro0=ro;}              
    }
    else { 
      ro *=2; 
      x += (ro-rold)* h; rold=ro; g =A(x);
      p = (g,h);    
      p1=p;
      ro1=ro;    
    }
  }  
}

template<class R,class DJ,class M,class P> 
int ConjuguedGradientNL(const M & dJ,const P & C,KN_<R> &x,const int nbitermax, double &eps,long kprint=1000000000)
{
//  -------------
   throwassert(&x && &b && &A && &C);
   typedef KN<R> Rn;
   int n=b.N();
   if (verbosity>99) kprint=1;
   throwassert(n==x.N());
   R ro=1;
   Rn g(n),h(n),Ah(n), & Cg(Ah);  // on utilise Ah pour stocke Cg  
   g = dJ(x);  
   Cg = C*g; // gradient preconditionne 
   h =-Cg; 
   R g2 = (Cg,g);
   if (g2 < 1e-30) 
    { if(verbosity>1)
       cout << "GCNL  g^2 =" << g2 << " < 1.e-30  Nothing to do " << endl;
     return 2;  }
   if (verbosity>5 ) 
     cout << " 0 GCNL  g^2 =" << g2 << endl;
   R reps2 =eps >0 ?  eps*eps*g2 : -eps; // epsilon relatif 
   eps = reps2;
   for (int iter=0;iter<=nbitermax;iter++)
     { 
       ro = argmin(J,x,h,g,Ah);
       
       Cg = C*g;
       R g2p=g2; 
       g2 = (Cg,g);
       if (  verbosity >1 )
         cout << "CGNL:" <<iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 
       if (g2 < reps2) { 
         if (verbosity )
            cout << "CGNL converge: " << iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 
          return 1;// ok 
          }
       R gamma = g2/g2p;       
       h *= gamma;
       h -= Cg;  //  h = -Cg * gamma* h       
     }
   cout << " CGNL: the method don't converge in " <<nbitermax <<" iterations \n" ;
   return 0; 
}
