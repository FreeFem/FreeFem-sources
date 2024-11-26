---
name: Connectivity
category: Mesh
folder: 3d
---

## Explains how to address the internal data of a mesh
On the mesh of a cube the following data can be accessed
~~~freefem
mesh3 Th=cube(1,1,1);

  // --------- new stuff -----------------
  int k=0,l=1,e=1;  // for example
  Th.nbe ; //   the number of boundary elements
  Th.nt ; // the number of elements
  Th.nv ;  // the number of vertices
  Th.be(k);   //   the boundary element k $\in \{0,...,Th.nbe-1\}$
  Th.be(k)[l];   //   the vertices l $\in \{0,1\}$ of  boundary element k
  Th.be(k).Element ;   //   the triangle contening the  boundary element k
  Th.be(k).whoinElement ; // egde number of triangle containing the  bdy element k
  Th.be(k).N ;   //   the Normal to be(k)   version 4.10.1
  Th.be(k).measure ;   //   the measure of be(k)   version 4.10.1
  Th[k].adj(e) ; //   adjacent triangle to k by edge e
  Th[k].measure ; //    the volume of element k
  Th(l).x ; // first coordinate of node l
  Th[k][0];  // the  number of the first of the 4 vertex  of element k.

// Here are some examples

  cout << " print mesh connectivity " << endl;
  Th[k] == Th[k].adj(e) ;// non adjacent element   the same
  Th[k] != Th[k].adj(e) ;// true adjacent triangle
   int nbelement = Th.nt;
  for (int i=0;i<Th.nv;++i)
  cout << i << " : "  << Th(i).x << " "<< Th(i).y << " " << Th(i).z  << endl; 
 
  for (int k=0;k<nbelement;++k)
    cout << k << " :  " << int(Th[k][0]) << " " << int(Th[k][1]) << " " <<  int(Th[k][2])<< " " <<  int(Th[k][2])
	 << " , label/ region  " << Th[k].label << endl;
  
  for (int k=0;k<nbelement;++k)
    for (int e=0,ee;e<4;++e) 
      //  remark FH hack:  set ee to e, and ee is change by method adj, 
      {
	    cout << k <<  " " << e << " <=>  " << int(Th[k].adj((ee=e))) << " " << ee  
	     << "  adj: " << ( Th[k].adj((ee=e)) != Th[k]) << endl;  
      }
      // note :     if k == int(Th[k].adj(ee=e)) it is not a adjacent element

  int nbboundaryelement = Th.nbe; 
  Th.be;
    for (int k=0;k<nbboundaryelement;++k)
      cout << k << " : " <<  Th.be(k)[0] << " " << Th.be(k)[1] << " " << Th.be(k)[2]  << " , label " << Th.be(k).label 
	   <<  " tet  " << int(Th.be(k).Element) << " " << Th.be(k).whoinElement <<  " N " << Th.be(k).N << endl; 

~~~
The following bounding box bb contains the triangulation. It is useful for limiting the region of a plot

~~~freefem
real[int] bb(4);
	boundingbox(Th,bb); // bb[0] = xmin, bb[1] = xmax, bb[2] = ymin, bb[3] =ymax 
	   cout << " boundingbox  xmin: " << bb[0] << " xmax: " << bb[1] 
	                     << " ymin: " << bb[2] << " ymax: " << bb[3] << endl; 
R3 O(0.5,0.5,0.5);
real ss =0;
 for (int k=0;k<nbboundaryelement;++k)
  ss += solidangle(O,Th.be(k));
 cout << " solid angle = " << ss << " == 4*pi == " << 4*pi << endl;
 assert( abs(ss-4*pi) < 1e-9);
 
 {
 Th = cube(3,3,3);
 func real f(R3 A)
 {
    assert(nuFace>=0); 	 
    cout << "P "<< P << " " << nuTriangle << " " << nuFace << " A = "<< A << endl;
    return solidangle(A,Th[nuTet],nuFace)/area;
 }
 cout << " integral " << int2d(Th,qforder=1)(f(O) )<< " " << 4*pi <<  endl; 
}
~~~
