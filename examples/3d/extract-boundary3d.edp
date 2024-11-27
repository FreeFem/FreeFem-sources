load "medit"
int n= 10;
int nvb = (n+1)^3 - (n-1)^3;// Nb boundary vertices
int ntb = n*n*12; // Nb of Boundary triangle 
mesh3 Th=cube(n,n,n);
int[int] ll=labels(Th);//  get all labels of the mesh .. 
cout << " the labels "<< ll << endl;
cout << Th.nbe << endl; 
meshS ThS=extract(Th,label=ll);// extract boundary of 3d Mesh with given label
cout <<"  mesh nv, nt, nbe" <<  ThS.nv << " - " << ThS. nt << " - " << ThS.nbe << endl; 
cout <<"  computed by hand: " << ntb  << " " << nvb  << endl; 
assert( ThS.nt == ntb) ;
assert( ThS.nv == nvb);
