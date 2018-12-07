  #ifndef WITH_NO_INIT
#include "ff++.hpp"
#endif
#include "AFunction_ext.hpp"

using namespace std;

#include <set>
#include <vector>
#include "splitsimplex.hpp"

using namespace  Fem2D;

Mesh3 * mytruncmesh(const Mesh3 &Th,const long &kksplit,int *split, bool kk, const int newbelabel);

struct Op_mytrunc_mesh3 : public OneOperator {
  typedef const Mesh3 *pmesh3;
  class Op: public E_F0mps   { 
  public:
    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =4;
    Expression nargs[n_name_param];
    
    Expression getmesh,bbb;
    long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
    KN<long> *  arg(int i,Stack stack) const{ return nargs[i] ? GetAny<KN<long> *>( (*nargs[i])(stack) ): 0;}
      
    Op(const basicAC_F0 &  args,Expression t,Expression b) : getmesh(t),bbb(b) 
    { args.SetNameParam(n_name_param,name_param,nargs); }
    AnyType operator()(Stack s)  const ;
  };
  
  E_F0 * code(const basicAC_F0 & args) const 
  { return new Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1])); }
  Op_mytrunc_mesh3() : 
    OneOperator(atype<pmesh3>(),atype<pmesh3>(),atype<bool>()) {};     
};

struct Op_GluMesh3tab : public OneOperator {
  typedef const Mesh3 *pmesh3;
  class Op: public E_F0mps   {
  public:
    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =1;
    Expression nargs[n_name_param];

    Expression getmeshtab;
    long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
 
    Op(const basicAC_F0 &  args,Expression t) : getmeshtab(t)
    { args.SetNameParam(n_name_param,name_param,nargs); }
    AnyType operator()(Stack s)  const ;
  };

  E_F0 * code(const basicAC_F0 & args) const
  { return new Op(args,t[0]->CastTo(args[0])); }
  Op_GluMesh3tab() :
    OneOperator(atype<const pmesh3>(),atype<KN<pmesh3>*>()) {};
};

basicAC_F0::name_and_type Op_GluMesh3tab::Op::name_param[Op_GluMesh3tab::Op::n_name_param] =
 {
   {  "labtodel",             &typeid(long)}
 };

basicAC_F0::name_and_type Op_mytrunc_mesh3::Op::name_param[Op_mytrunc_mesh3::Op::n_name_param] =
 {
   {  "split",             &typeid(long)},
   {  "label",             &typeid(long)},
     { "new2old", &typeid(KN<long>*)},  //  ajout FH pour P. Jovilet jan 2014
     { "old2new", &typeid(KN<long>*)}   //  ajout FH pour P. Jovilet jan 2014

 };


Mesh3 * mytruncmesh(const Mesh3 &Th,const long &kksplit,int *split, bool kk, const int newbelabel)
{
    
    static const int FaceTriangle[4]={3,0,1,2};  //={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}
    
    
    // computation of number of border elements and vertex without split
    int nbe = 0;
    int nt  = 0;
    int nv  = 0;
    int nvtrunc =0;
    int nedge=0;
    int nface=0;
    double hmin=1e100;
    R3 bmin,bmax;
    int nbeee=0,nbfi=0;
    const int kksplit2 = kksplit*kksplit;
    const int kksplit3 = kksplit2*kksplit;
    int ntsplit =0;
    int tagb[4]={1,2,4,8} ;
    KN<int> tagTonB(Th.nt);
    tagTonB=0;
    
    for( int ibe=0; ibe < Th.nbe; ibe++)
    {
        int iff;
        int it=Th.BoundaryElement(ibe,iff);
        tagTonB[it]|= tagb[iff];
        int ifff=iff,itt=Th.ElementAdj(it,ifff);
        if(itt >=0 &&  itt != it)
            tagTonB[itt]|= tagb[ifff];
    }
    
    for (int i=0;i<Th.nt;i++)
        if(split[i])
        {
            ++ntsplit;
            // computation of number of tetrahedrons
            nt=nt+kksplit3;
            // computation of number of border elements
            for (int j=0;j<4;j++)
            {
                int jt=j,it=Th.ElementAdj(i,jt);
                if ( (it==i || it <0)  || ! split[it]) nbeee++;// boundary face ...
                else nbfi++; // internal face count 2 times ...
                if(it==i || it <0) nbe += kksplit2;  //on est sur la frontiere
                else if (!split[it]) nbe += kksplit2; //le voisin ne doit pas etre decoupe
                else if ( (tagTonB[i]&tagb[j] ) != 0 && i<it) nbe += kksplit2; // internal boundary ..
            }
            
            for (int e=0;e<6;e++){
                hmin=min(hmin,Th[i].lenEdge(e));   // calcul de .lenEdge pour un Mesh3
            }
        }
    ffassert( nbfi %2 ==0) ;
    nface = nbeee + nbfi/2;
    double hseuil = (hmin/kksplit)/1000.;
    if(verbosity>5)
        cout << "  number of  not intern boundary faces = " << nbeee << ",  all faces  =  " << nbe << ", hseuil=" << hseuil <<endl;
    
    /* determination de bmin, bmax et hmin */
    
    KN<int> takevertex(Th.nv,-1);
    
    for(int i=0; i<Th.nt; i++){
        if(split[i])
        {
            const Tet &K(Th.elements[i]);
            for(int ii=0; ii<4; ii++)
            {
                int iv= Th.operator()( K[ii] );
                if( takevertex[iv] == -1 )
                {
                    bmin=Minc(Th.vertices[iv],bmin);
                    bmax=Maxc(Th.vertices[iv],bmax);
                    takevertex[iv]=nvtrunc++;
                }
            }
            
        }
    }
    
    if( kksplit > 1 )
    { // compute the number of slip edge ...
        nedge=0;
        HashTable<SortArray<int,2>,int> edges(3*nface,nface);
        for(int i=0; i<Th.nt; i++){
            if(split[i])
            {
                const Tet &K(Th.elements[i]);
                for(int e=0;e<6;++e)
                {
                    
                    int e1 = Th( K[ Th[i].nvedge[e][0] ] );
                    int e2 = Th( K[ Th[i].nvedge[e][1] ] );
                    SortArray<int,2> key(e1,e2);
                    if(!edges.find(key) )
                        edges.add(key,nedge++);
                }
            }
        }
    }
    if(verbosity>10) cout << "    -- nvertex  " << nvtrunc << ", nedges = "<< nedge
        << ", nfaces = " << nface << " ntet =" << ntsplit
        << endl
        << "    -- Euler/Poincare constante = " << nvtrunc-nedge+nface-ntsplit
        << endl;
    
    /* determination des vertex, triangles et tetrahedre obtenue apres splitting dans le Simplex */
    
    int nfacesub = kksplit2;
    int ntetsub = kksplit3;
    int nvsub = (kksplit+1)*(kksplit+2)*(kksplit+3)/6;
    int ntrisub = 4*kksplit2;
    
    R3 *vertexsub; //[nvsub];
    int *tetsub;   //[4*ntetsub];
    int *trisub;   //[4*kksplit*kksplit];
    
    SplitSimplex<R3>( kksplit, nvsub, vertexsub, ntetsub, tetsub);
    SplitSurfaceSimplex( kksplit, ntrisub, trisub);
    
    if(verbosity>3)
        cout << "  -- trunc (3d) : Th.nv= " << Th.nv << "kksplit="<< kksplit << endl;
    
    int ntnosplit  = nt/kksplit3;
    int nbenosplit = nbe/kksplit2;
    int nfacenosplit = (4*ntnosplit+nbenosplit)/2;
    nv = ntnosplit*(nvsub - 4*( (kksplit+1)*(kksplit+2)/2 - 3*(kksplit-1) -3 ) - 6*( kksplit-1 ) - 4);
    if(verbosity>100) cout << "       1) nv= " << nv << endl;
    nv = nv + nfacenosplit*( (kksplit+1)*(kksplit+2)/2 - 3*(kksplit-1) -3 );
    if(verbosity>100) cout << "       2) nv= " << nv << endl;
    nv = nv + nedge*( kksplit-1 );
    if(verbosity>100) cout << "       3) nv= " << nv << endl;
    nv = nv + nvtrunc;
    if(verbosity>100) cout << "       4) nv= " << nv << endl;
    
    
    
    int itt=0;
    int ie=0;
    
    
    Vertex3 *v=new Vertex3[nv];
    Tet *t = new Tet[nt];
    Tet *tt = t;
    
    Triangle3 *b  = new Triangle3[nbe];
    Triangle3 *bb = b;
    R3 hh = (bmax-bmax)/10.;
    EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,bmin-hh,bmax+hh,0);
    const R3 * pP[4];
    int np=0; // nb of new points ..

    {
        KN<R3>  vertextetsub(nvsub);
        KN<int> newindex (nvsub);
        KN<int> neworder(Th.nv,-1);
        int newnv = 0;
        for (int i=0; i<Th.nv; i++)
        	if (takevertex[i] != -1)
          		{
          			neworder[i] = newnv;
            		newnv++; 
          		}

        for(int i=0; i<Th.nt; i++)
            if(split[i])
            {
                const Tet &K(Th.elements[i]); 
                
                for(int ii=0; ii< 4; ii++)
                    pP[ii] = & K[ii];
                
                for( int iv=0; iv<nvsub; iv++)
                    (R3&) vertextetsub[iv]= vertexsub[iv].Bary(pP);     
                
                for( int iv=0; iv<nvsub; iv++)
                {
                    Vertex3 * pvi=gtree->ToClose(vertextetsub[iv],hseuil);
                    
                    if(!pvi)
                    {
                    	int newnum = neworder[Th.operator()( K[iv] )];

                        (R3&) v[newnum]   = vertextetsub[iv];
                        //numv = &K[iv] - Th.vertices;
                        v[newnum].lab = K.lab;
                        newindex[iv] = newnum;
                        gtree->Add( v[newnum] );
                        np++;  
                    }
                    else
                        newindex[iv] = pvi-v;
                    
                    ffassert( np <= nv );
                }    
            	       
                for( int ii=0; ii<ntetsub; ii++)
                {
                    int ivt[4];
                    for( int jj=0; jj< 4; jj++)
                    {
                        ivt[jj] = newindex[tetsub[4*ii+jj]];
                        assert( tetsub[4*ii+jj] < nvsub );
                        /*assert( ivt[jj] < np );*/
                    }
                    (tt++)->set( v, ivt, K.lab);
                    itt++;
                    assert( itt <= nt );
                }
                
                for (int j=0;j<4;j++)
                {
                    int jt=j,it=Th.ElementAdj(i,jt);
                    
                    if ( ( (tagTonB[i]&tagb[j]) ==0 ) &&  !(it==i || it <0)  && !split[it])
                    {
                        // new border not on boundary
                        int ivb[3];
                        
                        for( int ii=0; ii<nfacesub; ii++)
                        {
                            int iface = 3*FaceTriangle[j]*nfacesub+3*ii;
                            
                            for( int jjj=0; jjj<3; jjj++)
                            {
                                ivb[jjj] = newindex[ trisub[iface+jjj] ];
                                assert( trisub[ iface+jjj ] < nvsub );
                                /*assert( ivb[jjj] < np );*/
                            }
                            
                            (bb++)->set( v, ivb, newbelabel);
                            ie++;
                        }
                    }
                    assert( ie <= nbe);
                    
                }
            }
    }
    if(verbosity>10)
        cout  << "    ++ np=" << np << "==  nv=" << nv << endl;
    ffassert( np == nv); 
    if(verbosity>8)
        cout << "   -- Number of new  border face not on Border " << ie << endl;
    delete [] vertexsub; //[nvsub];
    delete [] tetsub;   //[4*ntetsub];
    delete [] trisub;   //[4*kksplit*kksplit];
    
    // split border elements
    int nv2Dsub   = (kksplit+1)*(kksplit+2)/4;
    int ntri2Dsub = kksplit2;
    R2 *vertex2Dsub; //[nvsub];
    int *tri2Dsub;   //[4*kksplit*kksplit];
    
    SplitSimplex<R2>( kksplit, nv2Dsub, vertex2Dsub, ntri2Dsub, tri2Dsub);
    
    
    for( int ibe=0; ibe < Th.nbe; ibe++)
    {
        int iff;
        int it=Th.BoundaryElement(ibe,iff);
        int ifff=iff,itt=Th.ElementAdj(it,ifff);
        if(itt<0) itt=it;
        if( split[it] == 0 && split[itt] == 0) continue; // boundary not on one element
        
        const Triangle3 &K(Th.be(ibe));
        int ivv[3];
        
        ivv[0] = Th.operator()(K[0]);
        ivv[1] = Th.operator()(K[1]);
        ivv[2] = Th.operator()(K[2]);
        
        R3 *vertextrisub = new R3  [nv2Dsub];
        int *newindex = new int[nv2Dsub];
        for( int iv=0; iv<nv2Dsub; iv++)
        {
            double alpha=vertex2Dsub[iv].x;
            double beta=vertex2Dsub[iv].y;
            
            vertextrisub[iv].x = (1-alpha-beta)*Th.vertices[ivv[0]].x + alpha*Th.vertices[ivv[1]].x + beta*Th.vertices[ivv[2]].x;
            vertextrisub[iv].y = (1-alpha-beta)*Th.vertices[ivv[0]].y + alpha*Th.vertices[ivv[1]].y + beta*Th.vertices[ivv[2]].y;
            vertextrisub[iv].z = (1-alpha-beta)*Th.vertices[ivv[0]].z + alpha*Th.vertices[ivv[1]].z + beta*Th.vertices[ivv[2]].z;
            
        }
        
        for( int iv=0; iv<nv2Dsub; iv++)
        {
            const Vertex3 &vi( vertextrisub[iv] );
            Vertex3 * pvi=gtree->ToClose(vi,hseuil);
            assert(pvi);
            newindex[iv] = pvi-v;
        }
        
        for( int ii=0; ii<nfacesub; ii++)
        {
            int ivb[3];
            for( int jjj=0; jjj<3; jjj++)
            {
                ivb[jjj] = newindex[ tri2Dsub[3*ii+jjj] ];
                assert( tri2Dsub[ 3*ii+jjj  ] < nvsub );
                if(verbosity > 199 ) cout << "        " << ivb[jjj] << " np:" << np<< endl;
                assert( ivb[jjj] < np );
            }
            
            (bb++)->set( v, ivb, K.lab);
            ie++;
            assert( ie <= nbe);
        }
        delete [] vertextrisub;
        delete [] newindex;
        
        
    }
    
    delete [] vertex2Dsub;   //[4*ntetsub];
    delete [] tri2Dsub;   //[4*kksplit*kksplit];
    
    
    if(verbosity>99)
    {
        cout << "nbofv initial" << Th.nv << endl;
        cout << "nv=" << nv << " np=" << np << endl;
        cout << "itt=" << itt << " nt=" << nt << endl;
        cout << "ie=" << ie << " nbe=" << nbe << endl;
    }
    ffassert( nv == np );
    ffassert( ie ==nbe);
    ffassert( itt == nt );
    
    //delete gtree;
    
    Mesh3 *Tht = new Mesh3( nv, nt, nbe, v, t, b);
    Tht->BuildGTree(); // Add JM. Oct 2010
    delete gtree;
    
    
    return Tht;
}

AnyType Op_mytrunc_mesh3::Op::operator()(Stack stack)  const {
    
  Mesh3 *pTh = GetAny<Mesh3 *>((*getmesh)(stack));
  Mesh3 &Th = *pTh;
  long kkksplit =arg(0,stack,1L);
  long label =arg(1,stack,2L);
   KN<long> * pn2o =  arg(2,stack);
KN<long> * po2n =  arg(3,stack);

  KN<int> split(Th.nt);
  split=kkksplit;
  MeshPoint *mp= MeshPointStack(stack),mps=*mp;
  long kk=0;
  long ks=kkksplit*kkksplit*kkksplit;
  for (int k=0;k<Th.nt;k++)
    { 
      const Tet & K( Th.elements[k] );
      R3 B(1./4.,1./4.,1./4.);  // 27/09/10 : J.Morice error in msh3.cpp
      mp->set(Th,K(B),B,K,0);
      if (  GetAny<bool>( (*bbb)(stack) )  ) kk++;
      else  split[k]=0  ;    
    }
  // *mp=mps;
  if (verbosity>1) 
    cout << "  -- Trunc mesh: Nb of Tetrahedrons = " << kk << " label=" <<label <<endl;
  Mesh3 * Tht = mytruncmesh(Th,kkksplit,split,false,label);
  
    if(pn2o)
    {
        pn2o->resize(kk*ks);
        KN<long> &n2o(*pn2o);
        int l=0;
        for(int k=0; k< Th.nt; ++k)
            if( split[k] )
                for(int i=0; i< ks; ++i)
                    n2o[l++] = k;
    }
    if(po2n)
    {
        po2n->resize(Th.nt);
        KN<long> &o2n(*po2n);
        int l=0;
        for(int k=0; k< Th.nt; ++k)
            if( split[k] )
            {
                o2n[k] = l;
                l+=ks;
            }
            else o2n[k]=-1;
    }

  Add2StackOfPtr2FreeRC(stack,Tht);//  07/2008 FH
  *mp=mps;
  return Tht;
 };

Mesh3 * GluMesh3tab(KN<pmesh3> * const & tab, long const & lab_delete)
{ 
  int flagsurfaceall = 0;

  int nbt=0;
  int nbe=0;
  int nbex=0;
  int nbv=0;
  int nbvx=0;
  
  double hmin=1e100;
  R3 Pn(1e100,1e100,1e100),Px(-1e100,-1e100,-1e100);
  const Mesh3 * th0=0;

  for(int i = 0;i<tab->n;i++)
    {
      const Mesh3 &Th3(*tab->operator[](i));
      th0=&Th3;
      if(verbosity>1)  cout << " determination of hmin : GluMesh3D + "<< Th3.nv << " " << Th3.nt << " "<< Th3.nbe << endl;
      
      nbt  += Th3.nt;
      nbvx += Th3.nv;
      nbex += Th3.nbe;
      
      for (int k=0;k<Th3.nt;k++){
	for (int e=0;e<6;e++){
	  hmin=min(hmin,Th3[k].lenEdge(e));   // calcul de .lenEdge pour un Mesh3
	}
      }
      
      for (int k=0;k<Th3.nbe;k++){
	for (int e=0;e<3;e++){
	  hmin=min(hmin,Th3.be(k).lenEdge(e));   // calcul de .lenEdge pour un Mesh3
	}
      }
      
      for (int ii=0;ii<Th3.nv;ii++){ 
	R3 P( Th3.vertices[ii].x, Th3.vertices[ii].y, Th3.vertices[ii].z);
	Pn=Minc(P,Pn);
	Px=Maxc(P,Px);     
      }
    } 
  
  if(verbosity > 1) cout << "      - hmin =" <<  hmin << " ,  Bounding Box: " << Pn << " "<< Px << endl;
  
  // probleme memoire
  Vertex3  *v= new Vertex3[nbvx];
  Tet      *t;
  if(nbt!=0) t= new Tet[nbt];
  Tet      *tt=t;
  Triangle3 *b= new Triangle3[nbex];
  Triangle3 *bb= b;
  
  ffassert(hmin>Norme2(Pn-Px)/1e9);
  double hseuil =hmin/10.;

  //int *NumSom= new int[nbvx];

  // VERSION morice
  if(verbosity > 1) cout << " creation of : BuildGTree" << endl;   
  EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,Pn,Px,0);  
  
  nbv=0;
  //int nbv0=0;
    for(int i = 0;i<tab->n;i++)
    {
      const Mesh3 &Th3(*tab->operator[](i));

      if(verbosity>1)  cout << " loop over mesh for create new mesh "<< endl;
      if(verbosity>1)  cout << " GluMesh3D + "<< Th3.nv << " " << Th3.nt <<" " << Th3.nbe << endl;
      //nbv0 =+Th3.nv;
     
      for (int ii=0;ii<Th3.nv;ii++){
	const Vertex3 &vi(Th3.vertices[ii]);
	Vertex3 * pvi=gtree->ToClose(vi,hseuil);

	   
	if(!pvi){
	  v[nbv].x = vi.x;
	  v[nbv].y = vi.y;
	  v[nbv].z = vi.z;
	  v[nbv].lab = vi.lab;
	  //NumSom[ii+nbv0] = nbv;
	  gtree->Add( v[nbv] );
	  nbv++;
	}
	/*
	  else{
	  NumSom[ii+nbv0] = pvi-v;
	  assert(pvi-v <nbv); 
	  }
	*/
      }
	  
      for (int k=0;k<Th3.nt;k++){
	const Tet  &K(Th3.elements[k]);
	int iv[4];
	iv[0]=gtree->ToClose(K[0],hseuil)-v;
	iv[1]=gtree->ToClose(K[1],hseuil)-v;
	iv[2]=gtree->ToClose(K[2],hseuil)-v;  
	iv[3]=gtree->ToClose(K[3],hseuil)-v;  
	(tt++)->set(v,iv,K.lab);
      }
      //nbv0 =+Th3.nv;
    }

  
  if(verbosity > 1) cout << " creation of : BuildGTree for border elements" << endl;
  Vertex3  *becog= new Vertex3[nbex];  
  //Vertex3  becog[nbex]; 
  EF23::GTree<Vertex3> *gtree_be = new EF23::GTree<Vertex3>(becog,Pn,Px,0);
  
  double hseuil_border = hseuil/3;
  //nbv0=0;
  for(int i = 0;i<tab->n;i++)
    {
      const Mesh3 &Th3(*tab->operator[](i));

    for (int k=0;k<Th3.nbe;k++)
      {
	const Triangle3 & K(Th3.be(k));
	if ((K.lab != lab_delete))//&&(K.lab != 3))
	{
	
	int iv[3];
	iv[0]=Th3.operator()(K[0]); 
	iv[1]=Th3.operator()(K[1]); 
	iv[2]=Th3.operator()(K[2]); 

	R cdgx,cdgy,cdgz;
	  
	cdgx = (Th3.vertices[iv[0]].x+ Th3.vertices[iv[1]].x+ Th3.vertices[iv[2]].x)/3.;
	cdgy = (Th3.vertices[iv[0]].y+ Th3.vertices[iv[1]].y+ Th3.vertices[iv[2]].y)/3.;
	cdgz = (Th3.vertices[iv[0]].z+ Th3.vertices[iv[1]].z+ Th3.vertices[iv[2]].z)/3.;
	 
	const R3 r3vi( cdgx, cdgy, cdgz ); 
	const Vertex3 &vi( r3vi);
	    
	Vertex3 * pvi=gtree_be->ToClose(vi,hseuil_border);
	if(!pvi){
	  becog[nbe].x = vi.x;
	  becog[nbe].y = vi.y;
	  becog[nbe].z = vi.z;
	  becog[nbe].lab = vi.lab;
	  gtree_be->Add( becog[nbe++]);
		  
	  int igluv[3];
	  igluv[0]= gtree->ToClose(K[0],hseuil)-v; //NumSom[iv[0]+nbv0];  
	  igluv[1]= gtree->ToClose(K[1],hseuil)-v; //NumSom[iv[1]+nbv0]; 
	  igluv[2]= gtree->ToClose(K[2],hseuil)-v; //NumSom[iv[2]+nbv0]; 
	 
	  (bb++)->set(v,igluv,K.lab);
	}
	}
      }

    //nbv0 =+Th3.nv;
  }
  delete gtree;
  delete gtree_be;
  delete [] becog;
  
  if(verbosity > 2) cout << " nbv="  << nbv  << endl;
  if(verbosity > 2) cout << " nbvx=" << nbvx << endl;
  if(verbosity > 2) cout << " nbt="  << nbt  << endl;
  if(verbosity > 2) cout << " nbe="  << nbe  << endl;
  if(verbosity > 2) cout << " nbex=" << nbex << endl;
  if(verbosity>1)
    {
      cout << "     Nb of glu3D  point " << nbvx-nbv;
      cout << "     Nb of glu3D  Boundary faces " << nbex-nbe << endl;
    }

  if(nbt==0){
    Mesh3 *mpq= new Mesh3(nbv,nbe,v,b);  
    if(flagsurfaceall==1) mpq->BuildBoundaryElementAdj();
    return mpq;
  }
  else{
    Mesh3 *mpq= new Mesh3(nbv,nbt,nbe,v,t,b); 
 /* 
    mpq->BuildBound();
    if(verbosity > 1) cout << "fin de BuildBound" << endl;
    mpq->BuildAdj();
    if(verbosity > 1) cout << "fin de BuildAdj" << endl;
    mpq->Buildbnormalv();  
    if(verbosity > 1) cout << "fin de Buildnormalv()" << endl;
    mpq->BuildjElementConteningVertex();
    if(verbosity > 1) cout << "fin de ConteningVertex()" << endl;
  */
    mpq->BuildGTree();
    if(verbosity > 2) cout << "fin de BuildGTree()" << endl;
    
    //Add2StackOfPtr2FreeRC(stack,mpq);
  
    return mpq;
  }
}

AnyType Op_GluMesh3tab::Op::operator()(Stack stack)  const {

  KN<const Mesh3*> *tab = GetAny<KN<const Mesh3*> *>((*getmeshtab)(stack));
  long labtodel = arg(0,stack,0);

  Mesh3 * Tht = GluMesh3tab(tab,labtodel);

  Add2StackOfPtr2FreeRC(stack,Tht);
  return Tht;
}

#ifndef _ALL_IN_ONE_
static void Load_Init()
{
	
  typedef const Mesh *pmesh;
  typedef const Mesh3 *pmesh3;
  
  if (verbosity && mpirank == 0)
    cout << " load: mymsh3  " << endl;

  Global.Add("truncvord","(", new Op_mytrunc_mesh3);
  Global.Add("gluemesh3","(",new Op_GluMesh3tab); 

}
LOADFUNC(Load_Init)
#endif
