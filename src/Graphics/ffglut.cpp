#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
//#include <pthread.h>
#include <limits>
#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <cmath>
using namespace std;

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <vector>
#include <list>
#include <map>
#include <utility>

#include "rgraph.hpp"
#include "fem.hpp"
#include "RNM.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"

#include "PlotStream.hpp"

extern long verbosity;

// add for the gestion of the endianness of the file.
//PlotStream::fBytes PlotStream::zott; //0123;
//PlotStream::hBytes PlotStream::zottffss; //012345678;
// ---- FH

using namespace Fem2D;
using std::numeric_limits;
const R pi=M_PI;//4*atan(1.); 
using namespace std;

int debug=1;
int casemouse=0,keyact=0;
#include "ffglut.hpp"

#include "ffthreads.hpp"

int version =0;

//Mutex MutexNextPlot;
Thread::Id tidRead=0;
bool NoMorePlot=false;
bool NoMorePlotTilte=false;
ThePlot *currentPlot=0, *nextPlot=0;
bool inThreadRead=false;
FILE *datafile=0;

static  bool TryNewPlot( void );


void LauchNextRead();
void WaitNextRead();
THREADFUNC(ThreadRead,fd);
//void * ThreadRead(void *fd);

int kread=-1;

int Fin(int code)
{
  WaitNextRead();
  if(!NoMorePlot && debug>2)
    cout << " exit before end  " << endl;
  exit(NoMorePlot ? 0  : 1);
}

int   ReadOnePlot(FILE *fp)
{ 
  int err=0;
   if(!fp) return -4; 
  err= feof(fp) ;
  if(err) return -2;
  err= ferror(fp) ;
  if(err) return -3;

  PlotStream f(fp);
  f.set_binary_mode();
  const char *  magic2="#!ffglutdata2..";
  const char *  magic3="#!ffglutdata3..";
  const int lmagic=strlen(magic2);
    char magicxx[32];
  err=0;
  // init ..
  if(kread==-1)
    {
      for(int i=0;i<lmagic;i++)
	{ int c=getc(fp);
	    magicxx[i]=c;
	  //err += c != magic[i];
	  //if(err) break;
	}
	magicxx[lmagic]='\0';
	if( strcmp(magicxx,magic2)==0)  version=2;
	else if( strcmp(magicxx,magic3)==0)  version=3;
	else err =1; 
	
      if(err) {
	if(debug>2)
	 cout << " Err read magic heading " << endl;
	 goto Lreturn;
	 //return err;
	}
      kread++;
      if(debug>2) cout << " Read entete " << version << endl;
	  int c1 =getc(fp);//
      if(c1==13)	  
        int c2 =getc(fp);//	


    }
  long cas; 
  f >> cas;
  err=-1;
  if (feof(fp)) goto Lreturn ;
  if((debug > 2)) cout << " ReadOnePlot " << kread+1<< " cas = " << cas << " " << nextPlot << endl;
  if(cas==PlotStream::dt_newplot)  
    {
      assert(nextPlot==0);
      nextPlot = new ThePlot(f,currentPlot,++kread);
	if(debug>1)	
      cout << " next is build " << nextPlot<< " wait :" << nextPlot->wait << " -> " << kread <<  endl;
      assert(nextPlot);
      err=0;
    }
  else 
    {
      err=1;
      cout << " Error Cas inconnue (skip) " << endl;
    }
 Lreturn:
  f.set_text_mode();
  return err;
}


void TimerNextPlot(int value)
{
  // the routine to  until the end of nextplot.
  // we use gluttimerfunc functionnaly
  //  remark, if we miss we retry.
  // -----
  //  if(debug) cout << " TimeNextPlot  " << endl;
  value=min(1000,(value*3)/2);// try at leat every 1 second (not to heavy computation)
  if(TryNewPlot())
    glutPostRedisplay();
  else 
    glutTimerFunc(value,TimerNextPlot,value);
}

int SendForNextPlot()
{
  //  to send a event to plot the date sheet.
  // and out a timer to wait to the end of read..
  // every 25/ second..  = 1000/25 = 40 ms
  if(NoMorePlot)
    {
    if((debug > 1)) cout << " send signal For Next plot, skip: No More Plot !  " << endl;
    return 0;
    }
  if((debug > 1)) cout << " Try to read read  plot "<< endl;
  //  put a timer for wait to the  end of read
  glutTimerFunc(40,TimerNextPlot,40);
  
  return 1;
}

static  bool TryNewPlot( void )
{
  // the routine to try to see if the next plot is read or not. 
  // -----------------------------------------------------------
  bool ret=false;
  if(debug>2)
    cout << "  TryNewPlot   plot : " << currentPlot << " next = " << nextPlot << endl;;
  if (nextPlot!=0)
    {
      WaitNextRead();
      if(debug>1) cout << " change current plot to: " << nextPlot << " et  Lock Plot . " << endl;;

      AllWindows[glutGetWindow()]->add(nextPlot);
      //if(currentPlot) delete currentPlot; //  a change fait dans add 
      // MutexNextPlot.WAIT();      
      currentPlot=nextPlot;
      nextPlot=0;
      // MutexNextPlot.Free();
      LauchNextRead();  
      ret=true;
    }
  return ret;    
}



int  signep4(int i0,int i1,int i2,int i3)
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
inline R3 bary(const R3 K[4],R f[4],int i0,int i1,R v)
{
  R d=f[i0]-f[i1];
  assert(fabs(d)>1e-20);
  R l1= (f[i0] - v)/ d;  //  == 1 si v = f[i1]  
  R l0 = 1. -l1;
  assert(l0 >=-1e-10 && l1 >= -1e-10);
  return K[i0]*l0 + K[i1]*l1; // == K[i1] si l1 ==1 => v = f[i1] 
}
void drawisoTet(const R3 K[4],R f[4],R v)
{
  static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;//{ {2,1,3},{0,2,3},{1,0,3},{0,1,2} };

  R3 P[4];
  int nP=0;
  int np[4],nm[4];
  int km=0,kp=0;
  for (int i=0;i<4;++i)
    {
      if(f[i]<=v) nm[km++]=i;
      if(f[i]>=v) np[kp++]=i;
    }
  
  //cout << "km kp "<< km << " " << kp << endl;
  int h=-1,b[3];
  if(kp==1 && km==3)
    {
      h = np[0];
      b[0]=nvfaceTet[h][0];
      b[1]=nvfaceTet[h][1];
      b[2]=nvfaceTet[h][2];
    }
  if(km==1 && kp == 3)
    {
      h = nm[0];
      b[0]=nvfaceTet[h][0];
      b[2]=nvfaceTet[h][1];
      b[1]=nvfaceTet[h][2];
    }
  if(kp==2 && km==2)
    {//  cas quad 
      if(signep4(nm[0],nm[1],np[0],np[1]) < 0)
	Exchange(nm[0],nm[1]);
      //  le tet m[0],nm[1],np[0],np[1] est positif
      P[0]=bary(K,f,nm[0],np[0],v);
      P[1]=bary(K,f,nm[0],np[1],v);
      P[2]=bary(K,f,nm[1],np[1],v);
      P[3]=bary(K,f,nm[1],np[0],v);
      nP=4;      
    }
  else if (h>=0)
    { // cas triangle 
      P[0]=bary(K,f,h,b[0],v);
      P[1]=bary(K,f,h,b[1],v);
      P[2]=bary(K,f,h,b[2],v);
      nP=3;
    }
  

  /*
    if(nP)
    {
    cout << "+ " << np[0] << " - " << nm[0] << endl;
    cout << nP << " ;  ";
    for(int i=0;i<nP;++i)
    cout << P[i] << " ;  ";
    cout << endl;
    
    }
  */
    if(nP)
    {
      if(nP>2)
	{
	  R3 N(R3(P[0],P[1])^R3(P[0],P[2]));
	  N /= N.norme();
	  glNormal3d(N.x,N.y,N.z);
	}
      glBegin(GL_POLYGON);
      for(int i=0;i<nP;++i)
	glVertex3f(P[i].x, P[i].y,P[i].z); // 
      glEnd();
    }

  //  verification de l'orientation
  assert(nP < 3 || det(P[0],P[1],P[2],K[np[0]]) >=0)   ;
  assert(nP < 3 || det(P[0],P[1],P[2],K[nm[0]]) <=0)   ;
  
}


int dichotomie(const KN_<double>  &viso,R v) 
{
    int i=0,j=viso.N(),k;
    if  (v <viso[0] || v >viso[j-1]) 
	return -1;  
    while (i<j-1)    
	if ( viso[k=(i+j)/2]> v) j=k;
	else i=k;
    return i;
}


map<int,OneWindow *> AllWindows;
int  ShowGlerror(const char *s)
{
    GLint error = glGetError();
    if ( error != GL_NO_ERROR )
	printf("Attention %s erreur : %x \n",s,error);   
    return error;   
}


//  def des couleurs de la tables 
void DefColor(float & r, float & g, float & b,
              int k,int nb, bool hsv,bool grey,KN<R> colors)
{
    int nbcolors = colors.N()/3;
    if(k<=0) {  r=g=b=1.;} //  white
    else if (k==1)  { r=g=b=0.; } // black
    else if (k >= nb)   {  r=g=b=0.;} // black
    else if (grey) { float gg = 0.1+0.9*float(k-2)/(nb-3); r=g=b=gg;} 
    else if (nbcolors<=1) {  
	float h=float(k-2)/(nb-2),s=1.,v=1.;
	hsvToRgb(h,s,v,r,g,b); 
    return;}     
    else   { //  interpolation dans la table hsv    
	int i= (k-2); 
	int j0= i*(nbcolors-1) / (nb-2);
	int j1=j0+1;
	int i0=  j0*(nb-2)/(nbcolors-1);
	int i1=  j1*(nb-2)/(nbcolors-1);
	int j03=j0*3,j13=j1*3;
	float a=float(i1-i)/(i1-i0),a1=1-a;
	if (hsv)
	  {
	      float h = colors[j03+0]*a + colors[j13+0]*a1;
	      float s = colors[j03+1]*a + colors[j13+1]*a1;
	      float v = colors[j03+2]*a + colors[j13+2]*a1;
	  hsvToRgb(h,s,v,r,g,b); }
	else 
	  {
	      r = colors[j03+0]*a + colors[j13+0]*a1;
	      g = colors[j03+1]*a + colors[j13+1]*a1;
	      b = colors[j03+2]*a + colors[j13+2]*a1;
	  }
    }     
    
}

template<class Mesh>
void Plot(const Mesh & Th,bool fill,bool plotmesh,bool plotborder,ThePlot & plot,GLint gllists,int * lok)
{
    glDisable(GL_DEPTH_TEST);

    ShowGlerror("begin Mesh plot");
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE); 
    R z1= plot.z0;
    R z2= plot.z0;
    
    
    double r=0,g=0,b=0;
    if((debug > 3)) cout<< " OnePlotMesh::Draw " << plotmesh << " " << plotborder << " " <<  Th.nbBrdElmts() << " " << z1 << " "  << z2 << endl;
    // plot.SetColorTable(16) ; 
    bool cc[3]= { plotborder , plotmesh && fill , plotmesh };
    int kk=0;
    //for(int i=0;i<3;i++)
    //  cout << cc[i] << " " << lok[i] << " , ";
    //cout << endl;
    if(cc[kk])
      if(lok[kk])   glCallList(gllists+kk);
      else 
	{ 
	  lok[kk]=1;
	  glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
	  glLineWidth(2); 
	  glBegin(GL_LINES);    
	  for (int i=0;i<Th.nbBrdElmts();i++)
	    {
		const typename  Mesh::BorderElement  & K(Th.be(i)); 
	      plot.color(1+abs(K.lab));
	      glVertex3d(K[0].x,K[0].y,z1);
	      glVertex3d(K[1].x,K[1].y,z1);
	      
	      
	    }
	  glEnd(); 	  
	  glLineWidth(1); 
	  glEndList();  // fin de la list
	}
      else ;
    
    kk++;	
    if(cc[kk])
      if(lok[kk])   glCallList(gllists+kk);
      else 
	{ 
	  lok[kk]=1;
	  glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
	  glPolygonMode(GL_FRONT,GL_FILL);//GL_FILL	
	  glBegin(GL_TRIANGLES);
	  for (int i=0;i<Th.nt;i++)
	    {
	      const typename  Mesh::Element & K(Th[i]);
	      plot.color(K.lab?1+abs(K.lab):0);
	      
	      //glColor3d(r,g,b);
	      int i0= Th(K[0]),  i1= Th(K[1]),   i2= Th(K[2]) ;    		
	      glVertex3d(K[0].x,K[0].y,z2);
	      glVertex3d(K[1].x,K[1].y,z2);
	      glVertex3d(K[2].x,K[2].y,z2);
	      
	    }    
	  glEnd();
	  glEndList();  //
	}
    
    kk++;
    if(cc[kk])
      if(lok[kk])   glCallList(gllists+kk);
      else 
	{ 
	  lok[kk]=1;
	  glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
	  glPolygonMode(GL_FRONT,GL_LINE);
	  glBegin(GL_TRIANGLES);
	  for (int i=0;i<Th.nt;i++)
	    {
		const  typename  Mesh::Element  & K(Th[i]);
		plot.color(fill? 1 : 1+abs(K.lab));
		int i0= Th(K[0]),  i1= Th(K[1]),   i2= Th(K[2]) ;    
		glVertex3d(K[0].x,K[0].y,z1);
		glVertex3d(K[1].x,K[1].y,z1);
		glVertex3d(K[2].x,K[2].y,z1);
		
	    }    
	  
	  glEnd();
	glEndList();  // fin de la list
      }

    ShowGlerror("end Mesh plot");

}


void Plot(const Mesh3 & Th,bool fill,bool plotmesh,bool plotborder,ThePlot & plot,GLint gllists,int * lok)
{
  typedef Mesh3::BorderElement BE;
  typedef Mesh3::Element Tet;
  glEnable(GL_DEPTH_TEST);
  /*
  if(fill)  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE); 
  else glPolygonMode(GL_FRONT_AND_BACK,GL_LINE); 
  */
  ShowGlerror("begin Mesh plot");
  
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE); 
  
  R z1= plot.z0;
  R z2= plot.z0;
  
  double r=0,g=0,b=0;
  
  bool cc[3]= { plotborder , plotborder && fill , plotmesh };
  
  int kk=0;
  if(cc[kk])
    if(lok[kk])   glCallList(gllists+kk);
    else 
      { 
	lok[kk]=1;
	glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
	glLineWidth(1); 
	glAlphaFunc ( GL_GREATER, 0.1 ) ;
	glEnable(GL_ALPHA_TEST) ;
	glLineStipple(1, 0x300C);
	glEnable(GL_LINE_STIPPLE);
	glBegin(GL_TRIANGLES);    
		for (int i=0;i<Th.nbe;i++)
	  {
	    const BE & K(Th.be(i)); 
	    plot.color(1+abs(K.lab),0.25);
	    R3 N(R3(K[0],K[1])^R3(K[0],K[2]));
	    N /= N.norme();
	    glNormal3d(N.x,N.y,N.z);
	    glVertex3d(K[0].x,K[0].y,K[0].z);
	    glVertex3d(K[1].x,K[1].y,K[1].z);
	    glVertex3d(K[2].x,K[2].y,K[2].z);
	  }
	glEnd(); 
	glDisable(GL_LINE_STIPPLE);
	glLineWidth(1); 
	glDisable(GL_ALPHA_TEST) ;

	glEndList();  // fin de la list	  
      }
  
  kk++;
  ShowGlerror("end Mesh plot");
  
}



void OnePlotError::Draw(OneWindow *win)
{
  initlist();
  ThePlot & plot=*win->theplot;
  win->SetScreenView() ;
  glColor3d(0.,0.,0.);
  cout << " Error plot item empty " << item <<  endl;
  int i = 4;
  char s[100];
  sprintf(s,"Warning the item %d fot the plot is empty",item);
  win->Show(s,4+item*2);
  win->SetView() ;
}

template<class Mesh>
void OnePlotMesh<Mesh>::Draw(OneWindow *win)
{
  initlist();
  ThePlot & plot=*win->theplot;
  Plot(*Th,plot.fill,true,true,plot,gllists,oklist);
  ShowGlerror("OnePlotMesh::Draw");
}
void OnePlotMesh3::Draw(OneWindow *win)
{
  initlist();
    ThePlot & plot=*win->theplot;
    Plot(*Th,plot.fill,true,true,plot,gllists,oklist);
    ShowGlerror("OnePlotMesh3::Draw");
}
void OnePlotFE3::Draw(OneWindow *win)
{
  initlist();
  
  ThePlot & plot=*win->theplot;
    ShowGlerror("begin OnePlotFE3 plot");
    ///    plot.SetDefIsoV();
    if(plot.fill && what==6)
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    else
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    
    if(what==6)
      glEnable(GL_DEPTH_TEST);
    else 
      glEnable(GL_DEPTH_TEST);
    win->setLighting();
    if(oklist[0])
      glCallList(gllists+0);
    else
      { 
	oklist[0]=1;
        glNewList(gllists+0,GL_COMPILE_AND_EXECUTE); // save  la list sans affichage
	int nsubV=Psub.N();
	int nsubT=Ksub.N()/4;
	KN<R3> Pn(nsubV);
	int nK=v.N()/ Th->nt;
	for(int k=0,o=0;k<Th->nt;++k, o+= nK)
	  {
	    const Mesh3::Element & K=(*Th)[k];
	    int ii[4];// 
	    R ff[4];
	    R3 Pt[4];
	    for(int i=0;i<nsubV;++i)
	      Pn[i]=K(Psub[i]);
	    int lK=0;
	    if(what==6)
	      for(int sk=0;sk<nsubT;++sk)
		{
		  
		  for(int l=0;l<4;++l)
		    {
		      int iv= Ksub[lK++];
		      ii[l]= iv;
		      Pt[l]=Pn[iv];
		      ff[l]=v[o+iv];
		    }
		  
		  for(int i=0;i< plot.Viso.N();++i)
		    {
		      plot.color(i+4);			    
		      drawisoTet( Pt,ff,plot.Viso[i]);
		    }
		  
		}
	  }
        glEndList();  // fin de la list
      }
    win->unsetLighting();
    ShowGlerror("b mesh  OnePlotFE plot");  
    Plot(*Th,false,plot.drawmeshes,plot.drawborder,plot,gllists+2,&oklist[2]);
    ShowGlerror("OnePlotFE::Draw");
}    



template<class Mesh>
 OnePlotFE<Mesh>::OnePlotFE(const Mesh *T,long w,PlotStream & f)
:OnePlot(w,2,5),Th(T)
{
    R2 P0,P1;
    Th->BoundingBox(P0,P1);
    Pmin=P0;
    Pmax=P1;
    if(version==2)
      {
	  long nsub;
	  
	  f>> nsub;
	  int nsubT=NbOfSubTriangle(nsub);
	  int nsubV=NbOfSubInternalVertices(nsub);
	  
	  Psub.resize(nsubV);
	  Ksub.resize(nsubT*3);
	  for(int i=0,j=0;i<nsubV;++i)
	    Psub[i]=SubInternalVertex(nsub,i);
	   
	  for(int sk=0,p=0;sk<nsubT;++sk)
	      for(int i=0;i<3;++i,++p)
		  Ksub[p]=numSubTriangle(nsub,sk,i);
	  
	  
      }
    else
      {	  f >> Psub ;
	  f >> Ksub ;
	  if(debug>2) {
	  cout << " Psub " << Psub << endl;
	  cout << " Ksub " << Ksub << endl;}

      }
    
    f>> v;
    if(what==1)
      {
	  fmin = min(fmin,v.min());
	  fmax = max(fmax,v.max());
      }
    else if (what==2)
      {  
	  //ffassert(0); // afaire
	  int n= v.N()/2;
	  for (int i=0,j=0;i<n;i++, j+=2)
	    {
		R2 u(v[j],v[j+1]);
		vmax = max(vmax,u.norme());
	    }
	  //cout << " vmax = " << vmax << endl; 
      }
    if(debug>3) cout << "OnePlotFE" << Th <<" " << what<< " " << Psub.N() << " " << Ksub.N()/3 <<" " << v.N() << endl; 
    ffassert(f.good());
    
}

template<class Mesh>
void OnePlotFE<Mesh>::Draw(OneWindow *win)
{
  initlist();
  ThePlot & plot=*win->theplot;
  ShowGlerror("begin OnePlotFE plot");
  //plot.SetDefIsoV();
  win->setLighting();
  //    OneWindow * win=plot.win;// bof bof  la struct est tres mauvaise . 
  assert(win);
  const Mesh & Th(*this->Th);
    int nsubT= Ksub.N()/3;//NbOfSubTriangle(nsub);
  int nsubV=Psub.N();//NbOfSubInternalVertices(nsub);
  int nK=v.N()/ Th.nt;
  if(debug>4)
  cout << "\t\t\tOnePlotMesh::Draw  " <<v.N() << " ,nt " << Th.nt << " " << nK << " " 
       << Psub.N() << " " << what << " ,nv " << Th.nv <<  endl;
  ffassert(v.N()== Th.nt*nK);
  ffassert(nK = nsubV*what);
  int o=0;
  KN<R2> Pn(Psub.N());
  if((debug > 10)) cout << " " <<nsubV  << " " << nsubT << endl;
  
  if(plot.fill && what==1)
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  else
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  
  if(what==2)
    glDisable(GL_DEPTH_TEST);
  else 
    glEnable(GL_DEPTH_TEST);

  R coef = plot.coeff;
  double xmin,xmax,ymin,ymax;
  win->getcadre(xmin,xmax,ymin,ymax);
  double d= Max(ymax-ymin,xmax-xmin);
  R kk = 4*win->hpixel;
  R cc = win->hpixel*40;
  
  int klist=0;
  bool  change=false;
  if( (what==1) ) 
    {
      if ( plot.fill) klist=1;
      change = win->changeiso ;
    }
  else if (what==2)
    change  =  win->changearrow ;
  if(debug>9)
    cout << change << " " << klist << " ... " << oklist[klist] << "  fill = " 
	 <<  plot.fill <<  " " <<  coef << endl;
  if (oklist[klist] && ! change )
    glCallList(gllists+klist);
  else 
    {
      //      R fmn,fmx,vmn,vmx;
            
      // win->theplot->dyn_bfv(win,fmn,fmx,vmn,vmx) ;
      //win->theplot->SetDefIsoV(0,0,fmn,fmx,vmn,vmx) ;
      oklist[klist]=1;
      glNewList(gllists+klist,GL_COMPILE_AND_EXECUTE); // save  la list aevc  affichage
      if(debug>100)
      cout << win->Bmin << ", Bmax:   " << win->Bmax << " Viso: "<< plot.Viso << endl;
      for(int k=0;k<Th.nt;++k, o+= nK)
	{
	  const typename Mesh::Element & K=Th[k];
	  for(int i=0;i<nsubV;++i)
	    Pn[i]=K(Psub[i]);// local to global coord. 
	  if(what==1)
	    for(int sk=0;sk<nsubT;++sk)
	      {
		int i0= Ksub[sk*3+0];//numSubTriangle(nsub,sk,0);
		int i1= Ksub[sk*3+1];//numSubTriangle(nsub,sk,1);
		int i2= Ksub[sk*3+2];//numSubTriangle(nsub,sk,2);
		
		R ff[3]={v[o+i0],v[o+i1],v[o+i2]};
		R2 Pt[3]={Pn[i0],Pn[i1],Pn[i2]};
		if(plot.fill)
		  plot.DrawIsoTfill( Pt, ff, plot.Viso,plot.Viso.N());
		else
		  plot.DrawIsoT( Pt, ff, plot.Viso,plot.Viso.N()); 
	      }
	  else // what ==2
	    for (int i=0,j=0;i<nsubV;++i)
	      {
		R2 P=Pn[i];
		R2 uv(v[o+j],v[o+j+1]);
		j+=2;
		R  l = Max(sqrt((uv,uv)),1e-30) ;
		int col = 2+dichotomie(plot.Varrow,l);
		if(debug>100) 
		  cout << uv << " l= " << l << " " << coef << " " <<col <<  endl;
		
		plot.color(2+col);
		uv = coef*uv;
		l *= coef;
		R2 dd = uv*(-0.01/l);
		R2 dn = dd.perp()*0.5;
		if (l*10000.< kk) continue;
		if (l < kk) 
		  uv = uv*(kk/l);
		else if (l> cc)
		  uv = uv*(cc/l);	   
		glBegin(GL_LINES);          
		
		win->Seg(P,P+uv);
		
		if (10*l>kk) {
		  win->Seg(P+uv,P+uv+dd+dn);
		  win->Seg(P+uv,P+uv+dd-dn);
		}
		glEnd();		
	      }

	}
      glEndList();  // fin de la list
    }
  
  // if(plot.drawmeshes)
  //  if(what==2)
  //  glEnable(GL_DEPTH_TEST);  
  ShowGlerror("b mesh  OnePlotFE plot");  
  win->unsetLighting();
  Plot(Th,false,plot.drawmeshes,plot.drawborder,plot,gllists+2,&oklist[2]);
  ShowGlerror("OnePlotFE::Draw");
}

template<class Mesh>
 void OnePlotFE<Mesh>::dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn,R & vmx) const 
{
  const Mesh & Th(*this->Th);
  int nsubT= Ksub.N()/3;//NbOfSubTriangle(nsub);
  int nsubV=Psub.N();//NbOfSubInternalVertices(nsub);
  int nK=v.N()/ Th.nt;
  ffassert(v.N()== Th.nt*nK);
  ffassert(nK = nsubV*what);
  int o=0;
  KN<R2> Pn(Psub.N());
  KN<R3> P3(Psub.N());
  double xmin,xmax,ymin,ymax;
  KN<int> inCadre(nsubV);
  win->getcadre(xmin,xmax,ymin,ymax);
  bool ccc=false;
  bool ddd=false;
  if(ddd)
  cout << " dyn__ .. " << endl;
  for(int k=0;k<Th.nt;++k, o+= nK)
    {
      const typename Mesh::Element & K=Th[k];
      inCadre=0;
      for(int i=0;i<nsubV;++i)
	Pn[i]=K(Psub[i]);// local to global coord. 
      for(int i=0;i<nsubV;++i)
	{
	  double f=0;
	  if(what==1) f=v[o+i];
	  gluProject(Pn[i].x,Pn[i].y,f,win->modelMatrix,win->projMatrix,win->viewport,
		      &P3[i].x,&P3[i].y,&P3[i].z);
	  if(ddd)
	    cout  <<P3[i]  << ", " ;
	  
	}// local to global coord. 
      if(ddd)
	cout << endl;
      for(int sk=0;sk<nsubT;++sk)
	{
	  int i1=Ksub[sk*3+0], i2=Ksub[sk*3+1], i3=Ksub[sk*3+2];
	  R3 P0= Minc(Minc(P3[i1],P3[i2]),P3[i3]);
	  R3 P1= Maxc(Maxc(P3[i1],P3[i2]),P3[i3]);
	  if( win->InRecScreen(P0,P1))
	      {
		if(debug>100)
		cout << " ???  " << P0 << " " << P1 << " ,  " << win->Bmin 
		     << " , " << win->Bmax << endl;
		inCadre[i1]=2;
		inCadre[i2]=2;
		inCadre[i3]=2;
	      }
	}
      for (int i=0,j=0;i<nsubV;++i)
	if(inCadre[i])
	  {
	    ccc=true;
	    if(what==1)
	      {
		R f=v[o+i];
		fmn=min(f,fmn);
		fmx=max(f,fmx);
		
	      }		    
	    else // what ==2
	      
	      {
		R2 uv(v[o+j],v[o+j+1]);
		j+=2;		
		R  l =(uv,uv) ;
		vmn=min(l,vmn);
		vmx=max(l,vmx);
	      }
	  }
      if(debug>100 && ccc)
	cout << " dny_bfv :  "  << fmn << " " << fmx << " " << vmn << " " << vmx 
	     <<  " : " << Pn[0] << endl;
      
    }
}



void OnePlotCurve::Draw(OneWindow *win)
{
  initlist();

  ThePlot & plot= *win->theplot;
  plot.SetColorTable(16) ;
    double z = plot.z0;
  
  glBegin(GL_LINE_STRIP);    
  plot.color(2);
  // cout << "nePlotCurve::Draw " << xx << " " << yy << endl;
  for (int i=0;i<xx.N();i++)
    {
      glVertex3d(xx[i],yy[i],z);
      
    }
  glEnd(); 
    
}

void OnePlotBorder::Draw(OneWindow *win)
{
  initlist();

  glDisable(GL_DEPTH_TEST);
  ThePlot & plot= *win->theplot;
  R h = 8*win->hpixel;
  
  double z = plot.z0;
  plot.SetColorTable(16) ; 
  
  // vector<vector<pair<long,R2> > > data;
  for(int i=0;i<data.size() ;++i)
    {
      vector<pair<long,R2> > & v=data[i];
      ShowGlerror("end OnePlotBorder::Draw  1");
      
      
      for(int j=1;j<v.size();++j)
	{
	  //	  cout <<v[j].first << endl;
	  plot.color(2+v[j].first);
	  R2 Po(v[j-1].second), Pn(v[j].second);
	  R2 uv(Po,Pn);
          double l = Max(sqrt((uv,uv)),1e-20);
          
          R2 dd = uv*(-h/l);
          R2 dn = dd.perp()*0.5;
	  glLineWidth(2); 
	  glBegin(GL_LINES);    
	  win->Seg(Po,Pn);
	  glEnd();
	  
	  glLineWidth(1);
	  glBegin(GL_LINES);    	  
	  if(j!=1)
	    {
	      win->Seg(Po,Po+dd+dn);
	      win->Seg(Po,Po+dd-dn);
	    }
	  glEnd();
	}
      
      ShowGlerror("end OnePlotBorder::Draw  2");
      
      glPointSize(7); 
      glBegin(GL_POINTS);    
      int l= v.size()-1;
      plot.color(2+v[0].first);
      glVertex3d(v[0].second.x,v[0].second.y,z);
      plot.color(2+v[l].first);
      glVertex3d(v[l].second.x,v[l].second.y,z);
      glEnd();
      glPointSize(1);
      ShowGlerror("end OnePlotBorder::Draw  3");
    }
  ShowGlerror("end OnePlotBorder::Draw");
  
}


OneWindow::OneWindow(int h,int w,ThePlot *p)
  :
  icurrentPlot(lplots.begin()), 
  lplotssize(0),
  height(h),width(w),theplot(0),hpixel(1),
  Bmin(0,0),Bmax(1,1),oBmin(Bmin),oBmax(Bmax),zmin(0),zmax(1),
  windowdump(false),help(false), rapz0(-1.),rapz(1),withlight(false),
  changearrow(true),changeiso(true) 
{
  
  add(p);
}


void OneWindow::set(ThePlot *p)
{
  //ffassert(p);
    bool first = !theplot;
    bool change = theplot != p;
    theplot=p;
    if(p)
      {
	plotdim=p->plotdim;
      }
    rapz0 =-1; // to recompute the defalut rapz
    //    p->win=this;
    //    if(first)
    DefaultView() ;
    
    
}

void OneWindow::add(ThePlot *p)
{
  if(p) {
    lplots.push_back(p);
    lplotssize++;
    ++icurrentPlot;
    if(icurrentPlot==lplots.end())
      --icurrentPlot;// the previous
    if(icurrentPlot != lplots.end())
      set(*icurrentPlot);
    if( lplotssize>10)
      {
	bool isfirst = theplot == *lplots.begin();  
	if(debug >1)
	cout << " delete a plot " << *lplots.begin() << endl;
	delete *lplots.begin();
	lplots.erase(lplots.begin());
        lplotssize--;
	if(isfirst) set(*lplots.begin()); // change to the next plot
      }
  }
  else 
    set(p);
}

void OneWindow::DefaultView() 
{
  if(theplot)
    {
      plotdim=theplot->plotdim;
      R3 A(theplot->Pmin),B(theplot->Pmax);
      R3 D(A,B);
      R dxy= max(D.x,D.y);
      zmax = theplot->fmax;
      zmin = theplot->fmin;
      theta=theplot->theta;
      phi=theplot->phi;
      if(theplot->datadim==3) rapz0=1;
      else   if(rapz0<=0)
	{ //  ( zmax-zmin )*rapz0 =  0.3 dxyy
	  rapz0  =  0.4* dxy/(zmax-zmin) ;
	  if(debug>2)
	    {
	      cout << " rapz0 = " << rapz0 ;
	      cout << " dz = " << zmax-zmin  << " dxy =" << dxy << endl;
	    }
	}
      rapz=rapz0;
      coef_dist=theplot->dcoef;
      focal=theplot->focal;      

      if(theplot->datadim==3)
	{
	  Bmin3=A;
	  Bmax3=B;
	}
      else
	{ // data plot 2d ou 1 d... 
	  if(theplot->boundingbox.size() ==4)
	    {
	      Bmin3.x=theplot->boundingbox[0];
	      Bmin3.y=theplot->boundingbox[1];
	      Bmax3.x=theplot->boundingbox[2];
	      Bmax3.y=theplot->boundingbox[3];	    
	    }
	  else 
	    {
	      Bmin3.x=A.x;
	      Bmin3.y=A.y;
	      Bmax3.x=B.x;
	      Bmax3.y=B.y;
	    }
	  Bmin3.z=theplot->fmin;
	  Bmax3.z=theplot->fmax;
	}
      Pvue3=(Bmin3+Bmax3)/2;
      
      
      
      D *=0.05;      
      if(theplot->boundingbox.size() !=4)
	{
	  A -= D;
	  B += D;
	}
      else
	{
	  R x1=theplot->boundingbox[0],y1=theplot->boundingbox[1];
	  R x2=theplot->boundingbox[2],y2=theplot->boundingbox[3];
	  A = R2(min(x1,x2),min(y1,y2));
	  B = R2(max(x1,x2),max(y1,y2));
	}
      
      if (theplot->aspectratio)
	cadreortho(A.p2(),B.p2());
      else 
	cadre(A.p2(),B.p2());
    }
  hpixel = (Bmax.x-Bmin.x)/width;
  
  // SetView() ;
}

void  OneWindow::SetScreenView() const
{

  glDisable(GL_TEXTURE_2D);
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity(); 
  glOrtho(0,width,0,height,-1,1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void  OneWindow::SetView()
{
  if(plotdim==3 && theplot)
    {
      glViewport(0, 0,width, height);
      
      glMatrixMode(GL_PROJECTION); 
      glLoadIdentity(); 
      R ratio= (double) width / (double)  height; 
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity(); 
      
      
      R aspect=ratio;
      R3 DD(Bmin3,Bmax3);
      DD.z *= rapz;
      R dmax= DD.norme();;
      R dist = 0.5*dmax/sin(focal/2)*coef_dist;
       cam.x=Pvue3.x+cos(phi)*cos(theta)*dist;
       cam.y=Pvue3.y+cos(phi)*sin(theta)*dist;
       cam.z=Pvue3.z*rapz+dist*sin(phi);  
      R znear=max(dist-dmax,1e-30);
      R zfare=dist+dmax;
      gluPerspective(focal*180./M_PI,aspect,znear,zfare);
      /*      
      if (eye)
	{
	  R dmm = -dmax*ceyes;
	  R dx = -dmm*sin(theta);
	  R dy = dmm*cos(theta);
	  camx += dx*eye;
	  camy += dy*eye;
	  }
      */
      if(debug>2)
	{
	  cout <<" setview 3d: rapz " <<  rapz << " cam: ";
	  cout << cam << " Pvue:" ;
	  cout << Pvue3  << " theta " << theta << "  phi = "<<  phi << endl;
	}
      gluLookAt(cam.x,cam.y,cam.z,Pvue3.x,Pvue3.y,Pvue3.z*rapz,0.,0.,1.);
      glScaled(1.,1.,rapz);   

      glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
      ShowGlerror(" Get PM");
      glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
      ShowGlerror(" Get MV");
      glGetIntegerv(GL_VIEWPORT,viewport);
      ShowGlerror(" Get VP");
  

      
    }
  else
    {
      ShowGlerror("Begin SetView");   
      glDisable(GL_DEPTH_TEST);
      glViewport(0, 0,width, height);
      R dz0,dz1,zm=0;
      if(plotdim==3)
	{
	    dz0=Bmin3.z;
	    dz1=Bmax3.z;
	}
      else
	{  
	    R zzmin = Min(zmin,theplot->fminT);
	    R zzmax = Max(zmax,theplot->fmaxT);    
	    R dz = (zzmax-zzmin);
	    zm=(zzmin+zzmax)*0.5;
	    //  to be sur  the the z size is no zero . 
	    dz = max(dz,(Bmax.x-Bmin.x)*0.1);
	    dz = max(dz,(Bmax.y-Bmin.y)*0.1);
	    dz0=-dz;
	    dz1 = dz;
	}
	
      if((debug>3 )) cout << "\t\t\t   SetView " << this << " " << Bmin  << " " 
			  << Bmax  << " dz  " << dz0 << " " << dz1  
			  << " theta " << theta << "  phi = "<<  phi << endl;
      ShowGlerror("0 Set MV");
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      ShowGlerror(" Set MV");
      glMatrixMode(GL_PROJECTION); 
      glLoadIdentity(); 
      ShowGlerror(" Set PM 1");
      glOrtho(Bmin.x,Bmax.x,Bmin.y,Bmax.y,dz0,dz1);
      
      ShowGlerror(" Set PM 2");
      
      R2 M=(Bmin+Bmax)/2.;
      glTranslated(0,0,-zm);
      
      //glLineWidth(1);
      //glColor3d(0.,0.,0.);

      glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
      ShowGlerror(" Get PM");
      glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
      ShowGlerror(" Get MV");
      glGetIntegerv(GL_VIEWPORT,viewport);
      ShowGlerror(" Get VP");
  
      
      ShowGlerror("End SetView ");
  }
    
}
void  OneWindow::resize(int w,int h)
{  double ww=width,hh=height;
    width=w;
    height=h;
    if (theplot->aspectratio)
      {	
	  cadreortho(oBmin,oBmax);
      }
    
}

void OneWindow::zoom(R coef)
{
  coef_dist*=coef;
  R2 M=(oBmin+oBmax)/2.;
  R2 D=(oBmax-oBmin)/2.;
  R2 A=  M - D*coef;
  R2 B=  M + D*coef;
  if (theplot->aspectratio)
    cadreortho(A,B);
  else 
    cadre(A,B);

}
void  OneWindow::zoom(int w,int h,R coef)
{
    GLdouble x=w,y=height-h,z=(zmin+zmax)/2.;
    GLdouble xx,yy,zz;
    

    GLint ok= gluUnProject( x,y,z,modelMatrix,projMatrix,viewport,&xx,&yy,&zz);
    ShowGlerror(" UnPro .. ");
    if(debug>2)
      cout << " ok " << ok << " " << x << " " << y << " " << z 
	   << " -> " << xx << " " << yy << " " << zz << endl;
    R2  oD(oBmin,oBmax);
    R2  D(Bmin,Bmax);
    R2 O(xx,yy);// oBmin.x+D.x*xx/width,oBmin.y+D.y*yy/height); 
    if((debug > 3)) cout<< " zoom : "  << this << " O " << O  
			<< " " << coef << " D = "<<  D<< "as "
			<< theplot->aspectratio <<  endl;
    oD *= 0.5*coef;
    R2 A = O - oD;
    R2 B = O + oD;
    if (theplot->aspectratio)
	cadreortho(A,B);
    else 
	cadre(A,B);
}
void OneWindow::MoveXView(R dx,R dy) 
{
  R3 D(Bmin3,Bmax3);
  R3 dd( dx*D.x*sin(theta),-dx*D.y*cos(theta), - dy*D.z);
  if(debug>2)
  cout << " MoveXView  "<< dx << " " << dy << " " << D << " mn: " << Bmin3 <<"  mx :" << Bmax3 << " d=" << dd << endl;
  Pvue3 += dd/50.;
  // 2d ...  add  FH   july 2009
  R2 D2(-dx*5*hpixel,-dy*5*hpixel);
  oBmin += D2;
  oBmax += D2;
  Bmin += D2;
  Bmax += D2;


  // cout << xm << " " << ym << " " << zm << endl;
}

void OneWindow::cadre(R2 A,R2 B)
{
      
    oBmin=Bmin=A;
    oBmax=Bmax=B;
    hpixel = (Bmax.x-Bmin.x)/width;    
    
    
}

void OneWindow::getcadre(double &xmin,double &xmax,double &ymin,double &ymax)
{
    xmin =  Bmin.x;
    xmax =  Bmax.x;
    ymin = Bmin.y;
    ymax = Bmax.y;
    
}
void OneWindow::Display()
{
  ffassert(this && theplot);
  SetScreenView() ;
  glColor3d(0.,0.,0.);
     
  if(help)
    {
      theplot->DrawHelp(this);
      help=false;
    }
  else
    {
  ShowGlerror("Begin Display");
  
  //  SetView();
  if(theplot)
    theplot->Draw(this);
  ShowGlerror("After Display");
    }
}
void OneWindow::cadreortho(R2 A, R2 B)
{
    R2 D(A,B);
    oBmin=A;
    oBmax=B;
    
    double cxy =  D.y*width/ (D.x*height);
    
    if ( D.y*width < D.x*height)  
	// width -> infty => D.x la ref
	D.y = D.x*(double) height/ width;
    else // height -> infty => D.y la ref
	D.x = D.y*(double) width/height;
    R2 M=(A+B)/2., D2=D/2.;
    
    Bmin= M - D2;
    Bmax= M + D2;
    
    if((debug > 10)) cout << " cadreortho: "<< " :: " << Bmin << " " << Bmax <<" oB " << oBmin << " " << oBmax << endl;
    
    // if((debug > 10)) cout << "cadreortho\n";
}
void OneWindow::setLighting()
{
  if(withlight)
    {
      if(plotdim==3)
	{
	  GLfloat lp0[4] = { cam.x,cam.y,cam.z, 1.0 };
	  glLightfv(GL_LIGHT0,GL_POSITION,lp0);	
	  
	  if(debug>1)  cout << " Light pos  3d:  " << cam << endl;
	}
      else
	{
	  
	  GLfloat position[] = {Pvue3.x,Pvue3.y,Pvue3.z+(Bmax3.z-Bmin3.z)*3,1.f} ;
	  glLightfv(GL_LIGHT0, GL_POSITION, position);
	  
	}
      
      float cca=0.3,ccd=1., ccs=0.8;
      GLfloat ambient[] = {cca,cca,cca,1.0f};//différents paramètres
      GLfloat diffuse[] = {ccd,ccd,ccd,1.0f};
      GLfloat specular_reflexion[] = {ccs,ccs,ccs,1.0f};
      GLubyte shiny_obj = 128;
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);
      glEnable(GL_LIGHTING);//positionnement de la lumière avec
      glLightfv(GL_LIGHT0,GL_AMBIENT,ambient);//les différents paramètres
      glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuse);
      
      glEnable(GL_COLOR_MATERIAL);//spécification de la réflexion sur les matériaux
      glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,ambient);
      glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,diffuse);
      // glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular_reflexion);// on peut le faire avant chaque objet
      //glMateriali(GL_FRONT_AND_BACK,GL_SHININESS,shiny_obj);//si on veut qu'ils aient des caractéristiques #
      glShadeModel(GL_FLAT);  
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0); 
    }
  else
    {
      glDisable(GL_LIGHTING);
      glDisable(GL_LIGHT0);
    }

}

void OneWindow::unsetLighting()
{
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
  //g->lightning=false;
}

OnePlotBorder::OnePlotBorder(PlotStream & f) 
  :OnePlot(4,2,1)
{ 
    long nbd;
    f>> nbd;
    data.resize(nbd);
    for(int i=0;i<nbd;++i)
      {
	  long n; 
	  f>> n;
	  //cout << n << endl;
	  data[i].resize(n+1);
	  for(int j=0;j<=n;++j)
	    {
		long l; 
		double x,y;
		f >> l>> x >> y;
		// cout << x << ' ' << y << ' ' <<  l << endl;
		R2 P(x,y);
		Pmin=Minc(Pmin,P);
		Pmax=Maxc(Pmax,P);
		data[i][j]=make_pair(l,P);
	    }
      }
    ffassert(f.good());
}

void OnePlot::GLDraw(OneWindow *win)
{
  ThePlot & plot= *win->theplot; 
  Draw(win);
  
  win->changeiso=0;
  win->changearrow=0; 
}

void ThePlot::DrawHelp(OneWindow *win) 
{
  int i = 1;
  win->Show("Enter a keyboard character in the FreeFem Graphics window in order to:",i++);
  
  i+=1;
  win->Show("enter) wait next plot",i++);
  win->Show("p)     previous plot (10 plots saved) ",i++);
  win->Show("ESC)   exit from ffglut",i++);
  win->Show("?)  show this help window",i++);
  win->Show("+) -)   zoom in/out  around the cursor 3/2 times ",i++);
  win->Show("=)   reset vue ",i++);
  win->Show("r)   refresh plot ",i++);
  win->Show("up, down, left, right) special keys  to  tanslate   ",i++);
  win->Show("3)   switch 3d/2d plot (in test)  keys : ",i++);
  win->Show("        z) Z) (focal zoom unzoom)  ",i++);
  win->Show("        H) h) switch increase or decrease the Z scale of the plot ",i++);
  win->Show("mouse motion)    ",i++);
  win->Show("   - left button)  rotate    ",i++);
  win->Show("   - right button)       zoom        (ctrl+button on mac) ",i++);
  win->Show("   - right button +alt)  tanslate    (alt+ctrl+button on mac)",i++);

  win->Show("a) A) increase or decrease the arrow size",i++);
  win->Show("B)  switch between show  border meshes or not",i++);
  win->Show("i) I) update or not: the min/max bound of the functions to the window",i++);
  win->Show("n) N) decrease or increase the number of iso value array ",i++);
  win->Show("b)  switch between black and white or color plotting ",i++);
  win->Show("g)  switch between grey or color plotting ",i++);
  win->Show("f)  switch between filling iso or iso line  ",i++);
  win->Show("l)  switch between lighting or not  ",i++);
  win->Show("v)  switch between show or not the numerical value of colors ",i++);
  win->Show("m)  switch between show or not  meshes  ",i++);
  win->Show("w)  window dump in file ffglutXXXX.ppm ",i++);
  win->Show("any other key : nothing ",++i);
}

void ThePlot::dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn,R & vmx) const 
{
  fmn=+1e100;
  fmx=-fmn;
  vmx=0;
  vmn=fmin;
  for (list<OnePlot *>::const_iterator i= plots.begin();i != plots.end(); ++i)
    {
      if(*i)  (*i)->dyn_bfv(win,fmn,fmx,vmn,vmx) ;
    }
  if(debug>4)
    cout << "dyn_bfv  " << fmn << " " << fmx << endl;
  if(fmn>fmx) fmn=fmin,fmx=fmax;
  if(vmn>vmx) vmn=0,vmx=vmax;
}

void ThePlot::Draw(OneWindow *win) 
{
    if((debug>1 )) 
      { 
	cout << "      ThePlot::Plot " << count << " " << this << " win " << win << " " << state ;
	for (list<OnePlot *>::iterator i= plots.begin();i != plots.end(); ++i)
	  cout << (**i).what;
	cout << endl;
	
      }
    if(state==0) {
      state=1;
      win->DefaultView();
    }

    win->SetView();
    for (list<OnePlot *>::iterator i= plots.begin();i != plots.end(); ++i)
	(*i)->Draw(win);
    
    if(cm || value)
      { //  screen plot ...
	win->SetScreenView();
	if(cm)
	  {
	    color(1);
	    win->DrawCommentaire(cm->c_str(),0.1,0.97);
	  }
	if(value)
	  {  
	    int k0=0;
	    if(withiso)
	      {
		win->PlotValue(Viso,k0,"IsoValue");
		k0=1+Viso.N();
	      }
	    if(witharrow)
	      win->PlotValue(Varrow,k0,"Vec Value");
	  }
	// for picking..
	win->SetView();
      }
    changeViso=false;
    changeVarrow=false;
    changeColor=false;
    changeBorder=false;
    changeFill=false;
}

void  ThePlot::SetColorTable(int nb)
{
    tbc.resize(nb);
    for (int i=0;i<nb;++i)
	tbc[i].set(i,nb,this);    
}



ThePlot::ThePlot(PlotStream & fin,ThePlot *old,int kcount)
  :  count(kcount), state(0),
     changeViso(true),changeVarrow(true),changeColor(true),
     changeBorder(true),changeFill(true), withiso(false),witharrow(false),
     plotdim(2),theta(30.*M_PI/180.),phi(20.*M_PI/180.),dcoef(1),focal(20.*M_PI/180.),
     datadim(1)
     
{
  
  hsv=true; // hsv  type   
  coeff=1;
  wait=0;
  value=false;
  fill=false;
  aspectratio=false;
  clean=true;
  uaspectratio=false;
  pViso=false;
  pVarrow=false;
  Niso=0;
  Narrow=20;
  bw=false;
  psfile=0;
  cm=0;
  grey=0;
  if(old) {
    grey=old->grey;
  }
  greyo=grey;
  drawborder=true;
  drawmeshes=false;
  add=false; 
  keepPV=false;
    
  Pmin=R3(+dinfty,+dinfty,+dinfty);
  fmin = +dinfty;    
  fmax = -dinfty;
  Pmax=R3(-dinfty,-dinfty,-dinfty);
  vmax=0;
  
  coefr=1;
  long dimpp=0;
  long cas;
  
  while(1)
    {
	fin >> cas;
	if((debug > 4)) cout << " read cas: " << cas << "  " << PlotStream::dt_endarg << endl;
	if(cas==PlotStream::dt_endarg) break;
	if(version==2)
	    switch (cas) {
		case  0: fin >> coeff; break;
		case 1: fin >> cm; break;
		case 2: fin >> psfile; break;
		case 3: fin >> wait; break;
		case 4: fin >> fill; break;
		case 5: fin >> value; break;
		case 6: fin >> clean; break;
		case 7: fin >> aspectratio;uaspectratio=true; break;
		case 8: fin >> boundingbox; break;
		case 9: fin >> Niso; break;
		case 10: fin >> Narrow; break;
		case 11: fin >> Viso;Niso=Viso.N();pViso=true; break;
		case 12: fin >> Varrow;Narrow=Varrow.N();pVarrow=true; break;
		case 13: fin >> bw; break;
		case 14: fin >> grey; break;
		case 15: fin >> colors; break;
		case 16: fin >> drawborder; break;
		case 17: fin >> dimpp; break;// ajout fevr 2008  v3.0.6
		case 18: fin >> add; break;
		case 19: fin >> keepPV; break;
		default: 
		    cout << "Fatal error: Unknow  case  : " << cas <<endl;
		    ffassert(0);
		    break;
	    }
	else if(version ==3)
	    switch (cas) {
		case  0: fin >= coeff; break;
		case 1: fin >= cm; break;
		case 2: fin >= psfile; break;
		case 3: fin >= wait; break;
		case 4: fin >= fill; break;
		case 5: fin >= value; break;
		case 6: fin >= clean; break;
		case 7: fin >= aspectratio;uaspectratio=true; break;
		case 8: fin >= boundingbox; break;
		case 9: fin >= Niso; break;
		case 10: fin >= Narrow; break;
		case 11: fin >= Viso;Niso=Viso.N();pViso=true; break;
		case 12: fin >= Varrow;Narrow=Varrow.N();pVarrow=true; break;
		case 13: fin >= bw; break;
		case 14: fin >= grey; break;
		case 15: fin >= colors; break;
		case 16: fin >= drawborder; break;
		case 17: fin >= dimpp; break;// ajout fevr 2008  v3.0.6
		case 18: fin >= add; break;
		case 19: fin >= keepPV; break;
		default: 
		    static int nccc=0;
		    if(nccc++<5)
			cout << " Skip Unknow case " << cas <<" (ffglut is too old ?)\n";
		    fin.SkipData();
		    break;
	    }   
	else ffassert(0);
	ffassert(fin.good() && ! fin.eof());
    }
  if(dimpp) plotdim=dimpp; 
  //    if( !uaspectratio) aspectratio= true;
  ffassert(cas==PlotStream::dt_endarg);
  if((debug > 2))
    {
      cout << " coeff " << coeff <<", ";
      if(cm)
	cout << " cm " << *cm <<", ";
      if(wait)	cout << " wait " ;
      if(fill)	cout << " fill " ;
      if(value)	cout << " value " ;
      if(bw)	cout << " bw " ;
      if(grey)	cout << " grey " ;
      if(drawborder)	cout << " drawborder " ;
      if(colors.N()) cout << "\n colors =" << colors;
      if(boundingbox.N()) cout << "\n bb  =" << boundingbox;
      
      cout << endl;
    } 
  fin.GetMeshes();  
  long nbmeshes;
  fin >> nbmeshes;
  if((debug > 2)) cout << " read nb : mesh " << nbmeshes << endl;
 if(version==2)
   {
  Ths.resize(nbmeshes);    
  for(int i=0;i<nbmeshes;++i)
    Ths[i]=0;
   }
   else
     {
	 Ths2.resize(nbmeshes);    
	 for(int i=0;i<nbmeshes;++i)
	     Ths2[i]=0;
     }
    
  for(int i=0;i<nbmeshes;++i)
    { 
      long l;
      fin >> l;
      if(l>=0) 
	{
	if((debug > 3)) cout << " read mesh " << i  << " -> " << l << "  " <<nbmeshes << endl;
	l--;
	ffassert(l>=0 && l < nbmeshes);
	if(version==2)
	  {
	ffassert(Ths[l]==0);
	fin >>Ths[l] ;
	  }
	else
	  {
	      ffassert(Ths2[l]==0);
	      fin >>Ths2[l] ;
	  }
	    
	if((debug > 3))
	  if(version==2)
	    cout << i << " nt/nv " << l << " "  <<Ths[l]->nt << " " << Ths[l]->nv << endl;
	  else
	    cout << i << " nt/nv " << l << " "  <<Ths2[l]->nt << " " << Ths2[l]->nv << endl;
	
	ffassert(fin.good());
	}
      else // Add FH optimisation FH 11/12/2008 (not use to day)
	{// the mesh is already in the previous plot with number ll
	  ffassert(l==-1);
	  long ll;
	  fin >> l>> ll; // read l and ll
	  ffassert(old);
	  if(version==2)
	    {
	      Ths[l]=old->Ths[ll];
	      Ths[l]->increment(); // 
	    }
	  else
	    {
		Ths2[l]=old->Ths2[ll];
		Ths2[l]->increment(); // 
	    }
	  
	}
      
    }
  long nbmeshes3=0;
  if(fin.GetMeshes3())
    { //  il y a des solution 3d; 
      
      fin >> nbmeshes3;
      if((debug > 2)) cout << " read nb : mesh3 " << nbmeshes3 << endl;
      Ths3.resize(nbmeshes3);
      for(int i=0;i<nbmeshes3;++i)
	Ths3[i]=0;
      for(int i=0;i<nbmeshes3;++i)
	{ 
	  long l;
	       fin >> l;
	       if(l>=0) 
		 {
		   if((debug > 3)) cout << " read mesh3 " << i  << " -> " << l 
					<< "  " <<nbmeshes3 << endl;
		   l--;
		   ffassert(l>=0 && l < nbmeshes3);
		   ffassert(Ths3[l]==0);
		   fin >>Ths3[l] ;
		   if((debug > 3))
		     cout << i << " nt/nv " << l << " "  <<Ths3[l]->nt << " " 
			  << Ths3[l]->nv << endl;
		   ffassert(fin.good());
		 }
	       else // Add FH optimisation FH 11/12/2008 (not use to day)
		 {// the mesh is already in the previous plot with number ll
		   ffassert(l==-1);
		   long ll;
		   fin >> l>> ll; // read l and ll
		   ffassert(old);
		     Ths3[l]=old->Ths3[ll];
		     Ths3[l]->increment(); // 
		 }
	       
	}	 
      
      
      
      fin.GetPlots();
    }
  long nbplot;
  int iso3d=0;
  fin >>nbplot;
  if((debug > 2)) cout << " nb item plot " << nbplot << endl;
  for(int i=0;i<nbplot;++i)
    {
      long what;
      OnePlot *p=0;
      fin >> what;
      long imsh;
      if((debug > 2)) cout << "    plot  " << i << " what " << what << endl;
      if(what !=3 && !uaspectratio) aspectratio= true;
      if(what==-1)  // gestion of error (empty plot)
	p = new OnePlotError(fin);
      else if(what==0) 
	{ 
	  
	  fin >> imsh;
	  if(version==2)
	    p=new OnePlotMesh<Mesh>(Ths[imsh-1]);
	  else
	    p=new OnePlotMesh<Mesh2>(Ths2[imsh-1]);
	  
	}
      else if (what==1 || what==2)
	{
	  fin >> imsh;
	  if(what==1) withiso=true;
	  else if (what==2) witharrow=true;
	  if((debug > 10)) cout << " plot : mesh " << imsh << endl;
	  ffassert(imsh>0 && imsh <=nbmeshes);
	  if(version==2)
	    p=new OnePlotFE<Mesh>(Ths[imsh-1],what,fin);
	  else 
	    p=new OnePlotFE<Mesh2>(Ths2[imsh-1],what,fin);
	}
      else if(what==3)
	p=new OnePlotCurve(fin);
      else if(what==4)
	p=new OnePlotBorder(fin);
      else if(what==5) 
	{ 	    
	  fin >> imsh;
	    p=new OnePlotMesh3(Ths3[imsh-1]);
	}
      else if (what==6 )
	{
	  iso3d++;
	  fin >> imsh;
	  if(what==6) withiso=true;
	  
	  if((debug > 10)) cout << " plot : mesh3 " << imsh << endl;
	  ffassert(imsh>0 && imsh <=nbmeshes3);
	  p=new OnePlotFE3(Ths3[imsh-1],what,fin);
	}
      
      else
	{
	  cout << "Bizarre unkown what :  " << what<< endl;
	  ffassert(0);
	}
      ffassert(p);
      plots.push_back(p);
      p->bb(Pmin,Pmax);
      p->bfv(fmin,fmax,vmax);
      plotdim=max(plotdim,p->dim);
      ffassert(fin.good());		      
      datadim=max(datadim,p->dim); 
    }
  if(Niso==0) 
    Niso = iso3d ? 5 : 20;
  
  // cout << "\t\t\t\t  f min, max v max :" << fmin << " " << fmax << " " << vmax << endl;
  
  double ref_f = abs(fmax)+abs(fmin) ; 
  if(fmax < fmin)
    {
      fmax = 1;
      fmin = 0;
    }
  else if( (fmax-fmin) <= 1e-8*ref_f)
    {
      if(ref_f< 1e-20) ref_f=0.5;
      fmax += ref_f/2;
      fmin -= ref_f/2;
    }
  PminT=Pmin;
  PmaxT=Pmax;
  fminT=fmin;
  fmaxT=fmax;
  if(old && 0)
    {
      Pmin= Minc(Pmin,old->PminT);
      Pmax= Maxc(Pmax,old->PmaxT);
      fmax= Max(fmax,old->fmaxT);
      fmin= Min(fmin,old->fminT);
    }
  
  z0= fminT +(fmaxT-fminT)*0.01;
  if((debug > 2)) cout << "               data bound: " << PminT << " " << PmaxT  
		       << " fmin == " << fminT << "  " << fmaxT 
		       << " z0 " << z0 <<  endl;
  fin.GetEndPlot(); 
  Viso.resize(Niso);
  Varrow.resize(Narrow);
  
  SetColorTable(Max(Niso,Narrow)+4) ;           
  SetDefIsoV() ;

}


void ThePlot::SetDefIsoV(int niso,int narr,double fmn,double fmx,double vmn,double vmx)
{
  bool dyn=false; 
  R d,x;
  
  if( fmx>fmn)
    {
      if(debug>3)
	cout << "Set Def dyn_bfv  " << fmn << " " << fmx << endl;

      if(niso>2) 
	Viso.resize(niso); 
      Niso=Viso.N();
      Narrow=narr;
      d =  (fmx-fmn)/(Niso-2) ;
      x =  (fmn+fmx)/2-d*0.5*(Niso-1);
      dyn=true;
    }
  else
    {
      d = 1 ? (fmaxT-fminT)/(Niso-2)  : (fmaxT-fminT)/(Niso-1);       
      x = 1 ? (fminT+fmaxT)/2-d*0.5*(Niso-1) :fminT+d/2;
    }
  if(!pViso || dyn)
    {
      for (int i = 0;i < Niso;i++)
	{Viso[i]=x;x +=d; }
      //if (fill ) {Viso[0]=fminT-d;Viso[Niso-1]=fmaxT+d;}    
    }
  dyn=false;

  if(vmx>vmn)
    {
      if(narr>2) 
	Varrow.resize(niso); 
      Narrow=Varrow.N();
      dyn=true;
      x = sqrt(vmn);
      d = (sqrt(vmx)-x)/(Narrow-1.1);
    }
  else
    {
      x=0; 
      if(debug>10)
	cout << "vmax=  " << vmax << endl;
      d= sqrt(vmax)/(Narrow-1.1);
    }   
  if (!pVarrow || dyn)
    for (int i = 0;i < Narrow;i++)
      {
	Varrow[i]=x;
	x +=d; 
      }
  if(debug>100)
    cout << " Viso ..; " << Viso <<endl;
  SetColorTable(Max(Niso,Narrow)+4) ; 
}

void OneWindow::Show(const char *str,int i)
{
  int hx= 15;
  int ix= width/20;
  int iy= height-hx*i;
  plot(ix,iy,str,3);
}

void  FillRectRasterPos(R x0,R y0,R x1,R y1)
{
  //  if((debug > 10)) cout << "FR Rp:   " << x0 << " " << y0 << " " << x1 << " " << y1 << endl;
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);//GL_FILL
  glBegin(GL_POLYGON);
  glVertex2d(x0, y0);
  glVertex2d(x1, y0);
  glVertex2d(x1, y1);
  glVertex2d(x0, y1);
  glEnd();
  
}

void  OneWindow::FillRect(R x0,R y0,R x1,R y1)
{
  
  double z1=  (zmin+zmax)/2; // line 
  glPolygonMode(GL_FRONT,GL_FILL);//GL_FILL
  glBegin(GL_POLYGON);
  glVertex3d(x0, y0,z1);
  glVertex3d(x1, y0,z1);
  glVertex3d(x1, y1,z1);
  glVertex3d(x0, y1,z1);
  glEnd();
}

void OneWindow::PlotValue(const KN_<double> & Viso,int  k0,const char * cmm)
{
  
  ShowGlerror("PlotValue b");
  //  glRasterPos2f(x,y);
  if((debug > 10)) cout << "PlotValue:" << cmm << " " << k0 << " " << width << " " <<height << endl;
  R xmin=0,xmax=width,ymin=0,ymax=height;
  if((debug > 10)) cout << "PlotValue " << Viso << endl;
  // int ix,iy;
  // GetSizeScreen(ix,iy);   
  
  R dx=(xmax-xmin);
  R dy=(ymax-ymin);
  //  10 points  
  // int kk = Max(30,iy/10);
  R h=10;
  R ho=h*1.1;
  R x0=xmin+dx*0.85;
  R y= ymin+dy*0.97;
  if((debug > 10)) cout << x0 << " " << y << " " << h <<  endl;
  y -= k0* ho;
  this->color(0);
  FillRectRasterPos(x0-h*0.5,y-h*(1.4*Viso.N()+0.3),x0+h*9,y+h*1.5);
  ShowGlerror("PlotValue m");
  this->color(1);

  plot(x0+ho,y,cmm);
  y -=  ho;
  for (int i=0;i<Viso.N();i++)
    {
      if((debug > 10)) cout << " x0 : " << x0<< " " << y << " " << h << " v = " << Viso[i] << endl;
      this->color(i+4);
      FillRectRasterPos(x0,y,x0+h,y+h);
      plot(x0+ho,y+3*h/10,Viso[i]);
      y -=  ho;
      ;                 
    }
  ShowGlerror("PlotValue f");
}


void OneWindow::DrawCommentaire(const char * cm,R x,R y) 
{
  
    R xmin=0,xmax=height,ymin=0,ymax=height;
    float dx=(xmax-xmin);
    float dy=(ymax-ymin);    
    plot(xmin+dx*x,ymin+dy*y,cm);   
}

/*
 void drwstr(R x,R y,char* format, ...) {
 va_list  args;
 char  *s,buffer[1024];
 
 va_start(args,format);
 vsnprintf(buffer,1024,format,args);
 va_end(args);
 }
 */
void  plot(double xx,double yy,const char *cmm,int font)
{
    glRasterPos2f(xx,yy);   
    float x[4];
    glGetFloatv(GL_CURRENT_RASTER_POSITION,x);
    if((debug > 10)) cout<<"avant x : "<<x[0]<<" y : "<<x[1]<<" z : "<<x[2]<< " " << xx <<" " << yy << endl;
    /*
      #define GLUT_BITMAP_9_BY_15((void*)2)
      #define GLUT_BITMAP_8_BY_13((void*)3)
      #define GLUT_BITMAP_TIMES_ROMAN_10((void*)4)
      #define GLUT_BITMAP_TIMES_ROMAN_24((void*)5)
      #define GLUT_BITMAP_HELVETICA_10((void*)6)
      #define GLUT_BITMAP_HELVETICA_12((void*)7)
      #define GLUT_BITMAP_HELVETICA_18((void*)8)
     */
    void * glut_font=GLUT_BITMAP_TIMES_ROMAN_10;    
    switch (font)
      {
      case  0: glut_font=GLUT_STROKE_ROMAN;break;
      case  1: glut_font=GLUT_STROKE_MONO_ROMAN;break;
      case  2: glut_font=GLUT_BITMAP_9_BY_15;break;
      case  3: glut_font=GLUT_BITMAP_8_BY_13;break;
      case  4: glut_font=GLUT_BITMAP_TIMES_ROMAN_10;break;
      case  5: glut_font=GLUT_BITMAP_TIMES_ROMAN_24;break;
      case  6: glut_font=GLUT_BITMAP_HELVETICA_10;break;
      case  7: glut_font=GLUT_BITMAP_HELVETICA_12;break;
      case  8: glut_font=GLUT_BITMAP_HELVETICA_18;break;


	
      }
    for (const char *s=cmm; *s; s++)
      {if((debug > 10)) cout << *s ;
	glutBitmapCharacter(glut_font,*s);
      }
    if((debug > 10)) cout << " ;;; " <<endl;
}
void plot(double x,double y,double i,int font)
{
    char buf[24];
    snprintf(buf,24,"%g",i,font);
    plot(x,y,buf);
}


void hsvToRgb (float h, float s, float v, float & r, float & g, float & b)
{
    int i;
    float aa, bb, cc, f;
    
    if (s == 0) /* Grayscale */
	r = g = b = v;
    else {
	if (h == 1.0) h = 0;
	h *= 6.0;
	i =  int(h);
	f = h - i;
	aa = v * (1 - s);
	bb = v * (1 - (s * f));
	cc = v * (1 - (s * (1 - f)));
	switch (i) {
	    case 0: r = v;  g = cc; b = aa; break;
	    case 1: r = bb; g = v;  b = aa; break;
	    case 2: r = aa; g = v;  b = cc; break;
	    case 3: r = aa; g = bb; b = v;  break;
	    case 4: r = cc; g = aa; b = v;  break;
	    case 5: r = v;  g = aa; b = bb; break;
	}
    }
}

void ThePlot::DrawIsoT(const R2 Pt[3],const R ff[3],const R * Viso,int NbIso, R rapz)
{
    glBegin(GL_LINES);
    R2 PQ[5];
    //int NbIso = Viso.N();
    for(int l=0;l< NbIso;l++)  /*    loop on the level curves */
      {
	  R xf = Viso[l];
	  int im=0;
	  for(int i=0;i<3;i++) // for the  3 edges 
	    {
		int j = (i+1)%3;
		R fi=(ff[i]);
		R fj=(ff[j]);
		
		if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
		  {
		      if (Abs(fi-fj)<=0.1e-10) 	/* one side must be drawn */
			{
			    color(l+4);
			    glVertex3f(Pt[i].x, Pt[i].y,xf*rapz);
			    glVertex3f(Pt[j].x, Pt[j].y,xf*rapz);

			    //MoveTo(Pt[i]);
			    //LineTo(Pt[j]);
			}
		      else
			{
			    R  xlam=(fi-xf)/(fi-fj);
			    PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
			}
		  }
	    }
	  
	  if (im>=2) /*    draw one segment */
	    {
		color(l+4);
		//MoveTo(PQ[0]);
		//LineTo(PQ[1]);
		glVertex3f(PQ[0].x, PQ[0].y,xf*rapz);
		glVertex3f(PQ[1].x, PQ[1].y,xf*rapz);
	    }
      }
    glEnd();
    
} 

void ThePlot::DrawIsoTfill(const R2 Pt[3],const R ff[3],const R * Viso,int NbIso, R rapz)
{
    R2 PQ[10];
    R z[10];
    
    R eps= (Viso[NbIso-1]-Viso[0])*1e-6;
    for(int l=1;l< NbIso;l++)  //   loop on the level curves 
      {
	  R xfb = Viso[l-1];
	  R xfh = Viso[l];
	  assert(xfb < xfh);
	  int im=0;
	  for(int i=0;i<3;i++) // for the  3 edges 
	    {
		int j=(i+1)%3;
		R fi=(ff[i]);
		R fj=(ff[j]);
		R xxfb =  xfb;
		R xxfh =  xfh;
		if (fj<fi ) Exchange(xxfb,xxfh);
		R xf  = xxfb;
		if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
		  {
		      if (Abs(fi-fj)>=0.1e-20)
			{
			    R  xlam=(fi-xf)/(fi-fj);
			    z[im] =  ff[i] * (1.F-xlam)  +  ff[j]* xlam;
			    PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
			    
			}
		  }
		xf = xxfh;  
		if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
		  {
		      if (Abs(fi-fj)>=0.1e-20)
			{
			    R  xlam=(fi-xf)/(fi-fj);
			    z[im] =  ff[i] * (1.F-xlam)  +  ff[j]* xlam;
			    PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
			}
		  }
		if (  xfb-eps <=fj  && fj <= xfh+eps) 
		    z[im]=ff[j],PQ[im++] = Pt[j];
		
	    }
	  if (im>2) 
	    {
	      color(l+4);
	      R3 P[10];
	      for(int i=0;i<im;++i)
		P[i]= R3(PQ[i].x,PQ[i].y,z[i]*rapz);
	      R3 N(R3(P[0],P[1])^R3(P[0],P[2]));
	      N /= N.norme();
	      if(N.z<0) N = -N;
	      glNormal3d(N.x,N.y,N.z);
	      
	      glBegin(GL_POLYGON);
	      //SetColor((xfb+xfh)/2); 
	      for (int i=0;i<im;i++)
		{// if((debug > 10)) cout << i << " \t : " << PQ[i].x << " " <<  PQ[i].y << " " << z[i]*rapz << endl;
		  glVertex3f(P[i].x, P[i].y,P[i].z);
		}
	      glEnd();
	      
	    }
      }
} 


bool WindowDump(int width,int height)
{
    int i,j;
    FILE *fptr;
    static int counter = 0;
    char fname[32];
    unsigned char *image;
    
    /* Allocate our buffer for the image */
    if ((image = new unsigned char[3*width*height]) == NULL) {
	fprintf(stderr,"WindowDump - Failed to allocate memory for image\n");
	return(false);
    }
    
    /* Open the file */
    sprintf(fname,"ffglut_%04d.ppm",counter);
    if ((fptr = fopen(fname, MODE_WRITE_BINARY)) == NULL) {
	fprintf(stderr,"WindowDump - Failed to open file for window dump\n");
	return(false);
    }
    if((debug > 10)) cout << " WindowDump in " << fname << endl;
    /* Copy the image into our buffer */
    glReadBuffer(GL_FRONT);
    glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,image);
    
    /* Write the PPM file */
    fprintf(fptr,"P6\n%d %d\n255\n",width,height);
    for (j=height-1;j>=0;j--) {
	for (i=0;i<width;i++) {
	    fputc(image[3*j*width+3*i+0],fptr);
	    fputc(image[3*j*width+3*i+1],fptr);
	    fputc(image[3*j*width+3*i+2],fptr);
	}
    }
    fclose(fptr);
    
    delete [] image;
    counter++;
    return(true);
}


void Clean() 
{
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

static void Reshape( int width, int height )
{   
    OneWindow * win=CurrentWin();
    if(win)
	win->resize(width,height);
    glutPostRedisplay();
}



void Display(void)
{ 
    OneWindow * win=CurrentWin();
    if(win) 
      {
	  /*    if (win->stereo)
	   { ffassert(0); 
	   
	   glClearColor(1.0, 1.0, 1.0, 0.0);
	   glDrawBuffer(GL_BACK_RIGHT);
	   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	   global->SetView(-1);
	   glCallList(TheDrawList);    
	   glClearColor(1.0, 1.0, 1.0, 0.0);
	   glDrawBuffer(GL_BACK_LEFT);
	   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	   global->SetView(+1);
	   glCallList(TheDrawList);    
	   
	   //win->Display();
	   glFlush();
	   glutSwapBuffers();
	   }
	   else */
	{
	    Clean();
	    win->Display();
	    glFlush();
	    glutSwapBuffers();
	    if ( win->windowdump)
		WindowDump(win->width,win->height);
	    win->windowdump=false;
	}
	  
      }

    if(!win->theplot->wait)
      SendForNextPlot();
    if(!NoMorePlotTilte  &&NoMorePlot)
      {
	NoMorePlotTilte=true;
	glutSetWindowTitle("FreeFem++ / Program ended; enter ESC to exit)");
      }
}

static void Mouse( int button,int state,int x,int y )
{
    // state up or down 
    OneWindow * win=CurrentWin();
     keyact = glutGetModifiers();
    switch(button)
      {
      case GLUT_LEFT_BUTTON:
        casemouse=GLUT_LEFT_BUTTON;
	if(win)
	  {
	    if(win && state == GLUT_DOWN) { win->xold=x,win->yold=y;return;}
	    win->phi += (y-win->yold)/(2.*180.);
	    win->theta -= (x-win->xold)/(2*180.);
	    glutPostRedisplay();
	  }
	break;
      case GLUT_RIGHT_BUTTON:
        casemouse=GLUT_RIGHT_BUTTON;
        if(win)
          {
            if(win && state == GLUT_DOWN) { win->xold=x,win->yold=y;return;}
          }
        break;



      }
}
static void MotionMouse(int x,int y )
{
    OneWindow * win=CurrentWin();
    switch(casemouse)
      {
      case GLUT_LEFT_BUTTON:
	
	if(win)
	  {
	    win->phi += (y-win->yold)/(2.*180.);
	    win->theta -= (x-win->xold)/(2*180.);
	    win->xold=x;
	    win->yold=y;
	    glutPostRedisplay();
	  }
	break;
      case GLUT_RIGHT_BUTTON:
        casemouse=GLUT_RIGHT_BUTTON;
	
        if(win)
          {
	    if(keyact & GLUT_ACTIVE_ALT)
	      {
                int dx = (x-win->xold);
                int dy = (y-win->yold);
		win->MoveXView(dx,-dy);
		glutPostRedisplay();
		{ win->xold=x,win->yold=y;}	      
	      }
	    else {
	      //  zoom en y 
	      R dd= (y-win->yold);
	      
	      { win->xold=x,win->yold=y;}
	      win->zoom(pow(0.99,dd));
	      glutPostRedisplay();

	    }
	  }
        break;
      }

}

static void Key( unsigned char key, int x, int y )
{
    OneWindow * win=CurrentWin();
    int ni=win->theplot->Viso.N();
    int na=win->theplot->Varrow.N();

    switch (key) 
      {
      case 27: // esc char
	Fin(0);
	break;
      case 'w':
	if(win)
		win->windowdump=true;
	break;
      case 'l':
	win->withlight = !win->withlight;
	break;
      case '?' :
	if(win)
	  win->help=true;
	
      case '+':
	win->zoom(x,y,0.7);
	win->coef_dist /= 1.2;
	break;
      case '-':	    
	win->zoom(x,y,1./0.7);
	win->coef_dist *= 1.2;
	break;
      case '3':
	
	win->plotdim=win->plotdim==2?3:2;
	break;
	/*
	  case '2':	    
	  win->plotdim=2;
	  break; */
      case '=':
	win->DefaultView();
	break;
      case 'f':
	win->theplot->fill = ! win->theplot->fill  ;
	break;
      case '@':
	win->dtheta = win->dtheta ? 0 : pi/1800.;
	break;
	
      case 'b':
	win->theplot->grey = ! win->theplot->grey   ;
	break;
      case 'g':
	win->theplot->grey = !  win->theplot->grey   ;
	break;
	
      case 'v':
	win->theplot->value = ! win->theplot->value  ;
	break;
      case 'm':
	win->theplot->drawmeshes = ! win->theplot->drawmeshes  ;
	break;
      case 'B':
	win->theplot->drawborder = ! win->theplot->drawborder  ;
	break;
      case 'H':
	win->rapz *= 1.2;
	break;
      case 'h':
	win->rapz /= 1.2;
	break;
      case 'p':
	if(win->icurrentPlot != win->lplots.begin())
	  win->set(*--win->icurrentPlot);
	break;
      case 'a':
	win->theplot->coeff/= 1.2;
	win->changearrow=true;
	break;
      case 'A':
	win->theplot->coeff*= 1.2;
	win->changearrow=true;
	break;
	
      case 'n':
	na  -=  na < 10  ? 2 : 5;
      ni  -=  ni <10 ? 2 : 5;	
      {
	R fmn,fmx,vmn,vmx;	  
	win->theplot->dyn_bfv(win,fmn,fmx,vmn,vmx) ;
	win->theplot->SetDefIsoV(ni,na,fmn,fmx,vmn,vmx) ;
	
	win->changeiso=true;
	win->changearrow=true;
      }
      break ;
      case 'N':
	na  += 5;
	ni  += 5;	
      case 'i':
	{
	  R fmn,fmx,vmn,vmx;	  
	  win->theplot->dyn_bfv(win,fmn,fmx,vmn,vmx) ;
	  win->theplot->SetDefIsoV(ni,na,fmn,fmx,vmn,vmx) ;
	  
	  win->changeiso=true;
	  win->changearrow=true;
	}
	break;
      case 'I':
	{
	  R fmn,fmx,vmn,vmx;	  
	  win->theplot->SetDefIsoV(ni,na) ;	  
	  win->changeiso=true;
	  win->changearrow=true;
	}
	break;
	
	
      case 'z':
	if(win->focal < M_PI/1.2 ) 
	  {
	    win->coef_dist*=sin(win->focal*1.2/2)/sin(win->focal/2);
	    win->focal *=1.2;
	  }
	break;
      case 'Z':
	if(win->focal > 1e-5)
	  {
	    win->coef_dist*=sin(win->focal/1.2/2)/sin(win->focal/2);
	    win->focal /=1.2;
	  }
	break;
      case '\r':
      case '\n':
	{
	  list<ThePlot*>::iterator ic = win->icurrentPlot;
	  if(++ic == win->lplots.end()) // last plot ->  try new one
	    SendForNextPlot();
	  else
	    win->set(*++win->icurrentPlot);
	  break;
	}
      default:
	if((debug > 10)) cout << " Key Character " << (int) key << " " << key << endl;  
	
      }
    glutPostRedisplay();
}


void SpecialKey(int key, int x, int y)
{
    OneWindow * win=CurrentWin();
    if(win)
      {
    // if((debug > 10)) cout << " SpecialKey " << key << " " << x << " " << y << " : ";
	R dx(0),dy(0);
	switch (key) {
	case  GLUT_KEY_LEFT:   dx = -1; break;
	case  GLUT_KEY_RIGHT:  dx = +1; break;
	case  GLUT_KEY_DOWN:   dy = -1; break;
	case  GLUT_KEY_UP:     dy = +1; break;
	}
	// calcul du deplacement de xm,ym,zm;
	// if((debug > 10)) cout << " " << dx << " " << dy << endl;
	win->MoveXView(dx,dy);
	glutPostRedisplay();
      }
}

void LauchNextRead()
{
  if(!NoMorePlot)
    {
      inThreadRead=true;
      tidRead = Thread::Start(&ThreadRead,(void *) datafile);
    }
}

void WaitNextRead()
{
  if( inThreadRead )
    {
      assert(tidRead!=0);
      Thread::Wait(tidRead);
      tidRead =0;
      inThreadRead=false;;
    }
}

//void * ThreadRead(void *fd)
THREADFUNC(ThreadRead,fd)
{   
  int err=0;
  assert(nextPlot==0);
  //  MutexNextPlot.WAIT(); 
  err=ReadOnePlot((FILE*)fd);
  // MutexNextPlot.Free(); 
  if(debug>1)
    cout << " We Read a plot  : " << kread << " " << nextPlot << " " << err << endl;
  if(err<0)
    NoMorePlot=true; 
  Thread::Exit();
}


int main(int argc,  char** argv)
{
    glutInit(&argc, argv);
    bool stereo=false;
    bool fullscreen = false;

    if(stereo)
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO);
    else  
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );

    if(argc>2) {
      if( strcmp(argv[1],"-nv")==0) debug=0;
      if( strcmp(argv[1],"-v")==0) debug=2,verbosity=2;
      if( strcmp(argv[1],"-vv")==0) debug=5,verbosity=2;
      if( strcmp(argv[1],"-vvv")==0) debug=10, verbosity=1000;
    }
    if(debug>1)		
    cout <<  " mode read = " << MODE_READ_BINARY << endl;
    datafile =0;;
    if(argc>1 && *argv[argc-1] != '-' ) 
     {	
	datafile=fopen(argv[argc-1], "r");
	if(debug >1)
	cout << " fopen :" << argv[argc-1] << " " <<datafile << endl;
     }

    if(datafile==0)
	datafile=stdin;
    if ( !datafile){
	cerr<< " Erreur fdopen stdin in binary " << endl;
	Fin(1);
     }
    int err=ReadOnePlot(datafile);
    if(err) {cout << "Err ReadOnePlot " << err << endl;
      Fin(1);}
   


    if(kread==0) {
      cout << " Error: no graphic data " << endl; 
      Fin(1);
    }
    if(debug>1) 
    cout << "on a lue le premier plot next plot: " << nextPlot << endl;


    int Height = 512;
    int Width = 512*3/2; 
    
    glutInitWindowSize(Width , Height);
    glutInitWindowPosition(100, 100);

    string titre = "FreeFem++: type return key to proceed (or ? for help on other)";
    glutCreateWindow(titre.c_str());
    glutPushWindow();
    if (fullscreen)
	glutFullScreen();
    
    AllWindows[glutGetWindow()]=new OneWindow(Width , Height,currentPlot); 
    TryNewPlot();
    
    glDisable(GL_DEPTH_TEST); 
    glutReshapeFunc( Reshape ); // pour changement de fenetre 
    glutKeyboardFunc( Key );    // pour les evenements clavier
    glutSpecialFunc(SpecialKey);
    glutMouseFunc(Mouse);       // pour les evenements sourie
    glutMotionFunc(MotionMouse); // les mouvements  de la sourie 
    glutDisplayFunc( Display ); // l'affichage
    glutMainLoop(); 

  return 0;
}


