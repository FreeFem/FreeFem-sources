verbosity=2;
mesh3 Th("dodecaedre01");
border cc(t=0,2*pi){x=cos(t);y=sin(t);label=1;}
mesh Th2=buildmesh(cc(50));
fespace Vh2(Th2,P2);
int nbtets=Th.nt;
cout << " Th mes " << Th.measure << " border mes " << Th.bordermeasure << endl;
cout << " nb of Tets = " << nbtets << endl;
if(1) {
  nbtets=2;
  for (int i=0;i<nbtets;i++)
    for (int j=0; j <4; j++)
      cout << i << " " << j << " Th[i][j] = "
	   << Th[i][j] << "  x = "<< Th[i][j].x  << " , y= "<< Th[i][j].y 
	   << ",  label=" << Th[i][j].label << endl;
	    
//   Th(i)   return   the vextex i of Th
//   Th[k]   return   the tet k of Th.

  // get vertices information : 
  int nbvertices=Th.nv;
  //nbvertices=2;
  cout << " nb of vertices = " << nbvertices << endl;
  for (int i=0;i<nbvertices;i++)
	cout << "Th(" <<i  << ") : "   // << endl;	
	     << Th(i).x << " " << Th(i).y  << " " << Th(i).z << " " << Th(i).label // version 2.19 
	  << endl;
 // version >3.4-1
  // --------- new stuff -----------------
  int k=0,l=1,e=1;
  Th.nbe ; // return the number of boundary element \hfilll
  Th.be(k);   // return the boundary element k $\in \{0,...,Th.nbe-1\}$ \hfilll
  Th.be(k)[l];   // return the vertices l $\in \{0,1\}$ of  boundary element k \hfilll
  Th.be(k).Element ;   // return the tet contening the  boundary element k \hfilll
  Th.be(k).whoinElement ;   // return the egde number of triangle contening the  boundary element k \hfilll
  Th[k].adj(e) ; // return adjacent tet to k by face e, and change the value of e to \hfilll
  // the corresponding face in the adjacent tet
  Th[k] == Th[k].adj(e) ;// non adjacent tet return the same 
  Th[k] != Th[k].adj(e) ;// true adjacent tet 
  Th.be(k).N;   // return the Normal  of  boundary element k \hfilll
  
  cout << " print mesh connectivity " << endl;
  int nbelement = Th.nt; 
  for (int k=0;k<nbelement;++k)
    cout << k << " :  " << int(Th[k][0]) << " " << int(Th[k][1]) << " " <<  int(Th[k][2]) 
         << " " <<  int(Th[k][3])
	 << " , label  " << Th[k].label << endl; 
  //  
  
  for (int k=0;k<nbelement;++k)
    for (int e=0,ee;e<4;++e) 
      //  remark FH hack:  set ee to e, and ee is change by method adj, 
      //  in () to make difference with  named parameters. 
      {
	    cout << k <<  " " << e << " <=>  " << int(Th[k].adj((ee=e))) << " " << ee  
	     << "  adj: " << ( Th[k].adj((ee=e)) != Th[k]) << endl;  
      }
      // note :     if k == int(Th[k].adj(ee=e)) not adjacent element 


  int nbboundaryelement = Th.nbe; 
  Th.be;
    for (int k=0;k<nbboundaryelement;++k)
      cout << k << " : " <<  Th.be(k)[0] << " " << Th.be(k)[1] << " , label " << Th.be(k).label 
	   <<  " tet  " << int(Th.be(k).Element) << " " << Th.be(k).whoinElement <<  " N=" << Th.be(k).N << endl; 
    
	  
savemesh(Th,"dd.meshb");
 }
fespace Vh(Th,P23d);
Vh xx=x;
if(xx[].n == Th.nv)
  for(int i=0;i<Th.nv;++i)
    assert(abs(Th(i).x-xx[][i])<1e-6);

func ue =   2*x*x + 3*y*y + 4*z*z + 5*x*y+6*x*z+1;
func f= -18. ;
Vh u=f,b,d,uhe=ue,bc;
cout << " Vh.ndof =  " <<  Vh.ndof << endl;
cout << "  Vh.ndofK " << Vh.ndofK << endl;
cout << Th[5].region << endl;
// cout << Th(0,0,0).region << endl;  a faire ...
cout << Th[5][3] << endl;  // ok.. 


for(int i=0;i<Vh.ndofK;++i )
  cout << Vh(11,i) << " ";
 cout << endl;

cout << ue(0.1,0.2,0.3)<< "  == " << f(0.1,0.2,0.3) << endl; ;
macro Grad3(u) [dx(u),dy(u),dz(u)]  // EOM

varf vbc(u,v) =  on(0,u=1);
varf vlap(u,v) = int3d(Th)(Grad3(v)' *Grad3(u)) + int3d(Th)(f*v) + on(0,u=ue);
varf vBord(u,v,solver=CG) = int2d(Th)(u*v) ;
verbosity=10; 
matrix A=vlap(Vh,Vh);
matrix B=vBord(Vh,Vh);
verbosity=1; 
b[]=vlap(0,Vh);
//bc[]=vbc(0,Vh);
//cout << bc[] <<endl;
{
ofstream fa("A.txt");
ofstream fb("B.txt");
fa << A ;
fb << b[] ;
}


cout << b[]. min << " " << b[].max << endl;
u[]=A^-1*b[];
cout << u[]. min << " " << u[].max << endl;
real err= int3d(Th)( square(u-ue) );
d= ue-u;
cout <<  " err = " << err <<  " " << d[].linfty << endl;
cout << " u (0,0,0) "<< u(0.,0.,0.) << endl;
cout << " dx(u) (0,0,0) "<< dx(u)(0.,0.,0.) << endl;
cout << " dy u (0,0,0) "<< dy(u)(0.,0.,0.) << endl;
cout << " dz u (0,0,0) "<< dz(u)(0.,0.,0.) << endl;
Vh2 u2=u(x,y,0.);
plot(u2,wait=1);
plot(u2,wait=1);
	{ ofstream file("dd.bb"); 
	file << "3 1 1 "<< u[].n << " 2 \n";
	int j;
	for (j=0;j<u[].n ; j++)  
	  file << d[][j] << endl; 
    }  
