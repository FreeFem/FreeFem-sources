
template<class R,class DJ>
R argmin(R  rho,const DJ & dJ, KN_<R> &x,KN_<R> &h,KN_<R> &g,KN_<R> &w)
{
  //  Find  ro such thah (dJ(x+ro h),h)  =0 
  // remark input: dJ(x)=g 
  int k=0;
//  g=dJ*x; //  pour est sure
  R ro0=0, ro=rho,ro1=rho,rold=0;
  R p0= (g,h),p,p1;
  if(p0>0) { h=-g; p0=(g,h);
    cout << "Reset searching directions to gradient! bofbof de F. hecht \n";
  } 
  R ap0=fabs(p0)*0.01; // on arrete quand on a divise par 100.
  
  x += (ro-rold)* h; rold=ro; g=dJ*x;// dJ(x,g);
  p= ( p1 = (g,h) );
  if (  verbosity >=50 )      
    cout << "        ro " <<  ro << " " << p 
         << " rh0= 0 " << p0 << endl;
  
  
  bool loop=true;
  while (k++<100 && loop)
    { //  calcul du  nouveau ro courant
      
      if (p0*p1 <0) { //  Ok changement de signe 
        R lambda= (-p0/(-p0+p1));
        if (lambda>0.8) lambda=0.8;
        if (lambda<0.2) lambda=0.2;
        ro = ro1*lambda + (1-lambda)*ro0 ;
        x += (ro-rold)* h; rold=ro; g=dJ*x;// dJ(x,g);
        assert(ro>1e-30 && ro < 1e+30);
        p = (g,h);
        if ( verbosity >=50 )
	     cout << "         " << ", rho=" << ro << " gh= " << p
             << "; ro0, gh0 = " << ro0 << " " << p0 
             << "; ro1, gh1 = " << ro1 << " " << p1 << " " << lambda ; 
        
        if(fabs(p) <= ap0 || k>100  ) {
          if (  verbosity >=50 )      
            cout << endl << endl;
          return ro; 
        }
        if(p0*p<0) { 
          p1=p;
          ro1=ro;
          if (  verbosity >=50 ) cout << " +\n";} 
        else {
          p0=p;
          ro0=ro;
          if (  verbosity >=50 ) cout <<" -\n";}              
      }
      else 
        { 
          ro *=2; 
          p0=p1;
          x += (ro-rold)* h; rold=ro; g=dJ*x;//dJ(x,g);
          p = (g,h);    
           p1=p;
           ro1=ro;    
           if (  verbosity >=50 ) cout <<p<<" " << ro <<  " 2* " ;
        }
      
    }  
  ExecError("NLCG: ArgMin loop (convexe minimization? )");
  return 0;
}

template<class R,class DJ,class P> 
int NLCG(const DJ & dJ,const P & C,KN_<R> &x,const int nbitermax, double &eps,long kprint=1000000000)
{
  //  -------------
  assert(&x && &dJ && &C);
  typedef KN<R> Rn;
  int n=x.N();
  
  R ro=1;
  Rn g(n),h(n),Ah(n), & Cg(Ah);  // on utilise Ah pour stocke Cg  
  g=dJ*x;// dJ(x,g);  
  Cg = C*g; // gradient preconditionne 
  h =-Cg; 
  R g2 = (Cg,g);
  if (g2 < 1e-30) 
    { if(kprint>1)
      cout << "GCNL  g^2 =" << g2 << " < 1.e-30  Nothing to do " << endl;
    return 2;  }
  if (kprint>5 ) 
    cout << " 0 GCNL  g^2 =" << g2 << endl;
  R reps2 =eps >0 ?  eps*eps*g2 : -eps; // epsilon relatif 
  eps = reps2;
  for (int iter=0;iter<=nbitermax;iter++)
    { 
      ro = argmin(ro,dJ,x,h,g,Ah);
      
      Cg = C*g;
      R g2p=g2; 
      g2 = (Cg,g);
      if (  kprint >1 )
        cout << "CGNL:" <<iter <<  ",  ro = " << ro << " ||g||^2 = " << g2 << endl; 
      if (g2 < reps2) { 
        if (kprint )
          cout << "CGNL converge: " << iter <<",  ro = " << ro << " ||g||^2 = " << g2 << endl; 
        return 1;// ok 
      }
      R gamma = g2/g2p;       
      h *= gamma;
      h -= Cg;  //  h = -Cg * gamma* h       
    }
  if(verbosity)
  cout << " Non convergence de la mŽthode du gradient conjugue NL " <<endl;
  return 0; 
}
