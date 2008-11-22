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
#include "PlotStream.hpp"



using namespace Fem2D;
using std::numeric_limits;
const R pi=M_PI;//4*atan(1.); 
using namespace std;

int debug=2;

#include "ffglut.hpp"

#include "ffthreads.hpp"

//Mutex MutexNextPlot;
Thread::Id tidRead=0;
bool NoMorePlot=false;
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
  if(!NoMorePlot)
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
  const char *  magic="#!ffglutdata...\n";
  err=0;
  // init ..
  if(kread==-1)
    {
      for(int i=0;i<strlen(magic);i++)
	{ int c=getc(fp);
	  err += c != magic[i];
	  if(err) break;
	}
      if(err) return err;
      kread++;
      if(debug>0) cout << " Read entete " << endl;
    }
  err=1;
  
  long cas; 
  f >> cas;
  if (feof(fp)) return -1;
  if((debug > 1)) cout << " ReadOnePlot " << kread+1<< " cas = " << cas << " " << nextPlot << endl;
  if(cas==PlotStream::dt_newplot)  
    {
      assert(nextPlot==0);
      nextPlot = new ThePlot(f,currentPlot,++kread);
      cout << " next is build " << nextPlot<< " wait :" << nextPlot->wait << " -> " << kread <<  endl;
      assert(nextPlot);
      err=0;
    }
  else 
    cout << " Error Cas inconnue (skip) " << endl;
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
    if((debug > 0)) cout << " send signal For Next plot, skip: No More Plot !  " << endl;
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
      AllWindows[glutGetWindow()]->set(nextPlot);
      if(currentPlot) delete currentPlot;
      // MutexNextPlot.WAIT();      
      currentPlot=nextPlot;
      nextPlot=0;
      // MutexNextPlot.Free();
      LauchNextRead();  
      ret=true;
    }
  return ret;    
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

void Plot(const Mesh & Th,bool fill,bool plotmesh,bool plotborder,ThePlot & plot)
{
  glDisable(GL_DEPTH_TEST);
    ShowGlerror("begin Mesh plot");
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE); 
    R z1= plot.z0;
    R z2= plot.z0;
    
    
    double r=0,g=0,b=0;
    if((debug > 3)) cout<< " OnePlotMesh::Draw " << plotmesh << " " << plotborder << " " <<  Th.neb << " " << z1 << " "  << z2 << endl;
    // plot.SetColorTable(16) ; 
    if(plotborder)
      {
	  glLineWidth(2); 
	  glBegin(GL_LINES);    
	  for (int i=0;i<Th.neb;i++)
	    {
		const BoundaryEdge & K(Th.be(i)); 
		plot.color(1+abs(K.lab));
		glVertex3d(K[0].x,K[0].y,z1);
		glVertex3d(K[1].x,K[1].y,z1);
		
		
	    }
	  glEnd(); 
      }
    glLineWidth(1); 
    if(plotmesh)
      {
	  if(fill)
	    {
		glPolygonMode(GL_FRONT,GL_FILL);//GL_FILL	
		glBegin(GL_TRIANGLES);
		for (int i=0;i<Th.nt;i++)
		  {
		      const Triangle & K(Th[i]);
		      plot.color(K.lab?1+abs(K.lab):0);
		      
		      //glColor3d(r,g,b);
		      int i0= Th(K[0]),  i1= Th(K[1]),   i2= Th(K[2]) ;    		
		      glVertex3d(K[0].x,K[0].y,z2);
		      glVertex3d(K[1].x,K[1].y,z2);
		      glVertex3d(K[2].x,K[2].y,z2);
		      
		  }    
		glEnd();
	    }
	  glPolygonMode(GL_FRONT,GL_LINE);
	  glBegin(GL_TRIANGLES);
	  for (int i=0;i<Th.nt;i++)
	    {
		const Triangle & K(Th[i]);
		plot.color(fill? 1 : 1+abs(K.lab));
		
		//glColor3d(r,g,b);
		int i0= Th(K[0]),  i1= Th(K[1]),   i2= Th(K[2]) ;    
		glVertex3d(K[0].x,K[0].y,z1);
		glVertex3d(K[1].x,K[1].y,z1);
		glVertex3d(K[2].x,K[2].y,z1);
		
	    }    
	  
	  glEnd();
      }
    ShowGlerror("end Mesh plot");

}

void OnePlotMesh::Draw(OneWindow *win)
{
  ThePlot & plot=*win->theplot;
  Plot(*Th,plot.fill,true,true,plot);
  ShowGlerror("OnePlotMesh::Draw");
}
void OnePlotFE::Draw(OneWindow *win)
{

  ThePlot & plot=*win->theplot;
  ShowGlerror("begin OnePlotFE plot");
  plot.SetDefIsoV();
  //    OneWindow * win=plot.win;// bof bof  la struct est tres mauvaise . 
  assert(win);
  const Mesh & Th(*this->Th);
  int nsubT=NbOfSubTriangle(nsub);
  int nsubV=NbOfSubInternalVertices(nsub);
  int nK=v.N()/ Th.nt;
  if(debug>4)
  cout << "\t\t\tOnePlotMesh::Draw  " <<v.N() << " ,nt " << Th.nt << " " << nK << " " 
       << nsubV << " " << what << " ,nv " << Th.nv <<  endl;
  ffassert(v.N()== Th.nt*nK);
  ffassert(nK = nsubV*what);
  int o=0;
  KN<R2> Pn(nsubV);
  if((debug > 10)) cout << " " << nsubT << " " << nsubV << endl;
  
  if(plot.fill && what==1)
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  else
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  if(what==2)
    {
      glDisable(GL_DEPTH_TEST);
    }
  R coef = plot.coeff;
  double xmin,xmax,ymin,ymax;
  win->getcadre(xmin,xmax,ymin,ymax);
  double d= Max(ymax-ymin,xmax-xmin);
  R kk = 4*win->hpixel;
  R cc = win->hpixel*40;
  
  {
    for(int k=0;k<Th.nt;++k, o+= nK)
      {
	const Mesh::Element & K=Th[k];
	for(int i=0;i<nsubV;++i)
	  Pn[i]=K(SubInternalVertex(nsub,i));
	if(what==1)
	  for(int sk=0;sk<nsubT;++sk)
	    {
	      int i0=numSubTriangle(nsub,sk,0);
	      int i1=numSubTriangle(nsub,sk,1);
	      int i2=numSubTriangle(nsub,sk,2);
	      
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
	      R2 dd = uv*(-0.005/l);
	      R2 dn = dd.perp()*0.5;
	      if (l*10000.< kk) continue;
	      if (l < kk) 
		uv = uv*(kk/l);
	      else if (l> cc)
		uv = uv*(cc/l);	   
	      glBegin(GL_LINES);          
	      
	      win->Seg(P,P+uv);
	      
	      if (l>kk) {
		win->Seg(P+uv,P+uv+dd+dn);
		win->Seg(P+uv,P+uv+dd-dn);
		}
	      glEnd();		
	    }
	
      }
  }
  // if(plot.drawmeshes)
  if(what==2)
    glEnable(GL_DEPTH_TEST);  
  ShowGlerror("b mesh  OnePlotFE plot");  
  Plot(Th,false,plot.drawmeshes,plot.drawborder,plot);
  ShowGlerror("OnePlotFE::Draw");


  
  
}

void OnePlotCurve::Draw(OneWindow *win)
{
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
  : height(h),width(w),theplot(0),hpixel(1),
    Bmin(0,0),Bmax(1,1),oBmin(Bmin),oBmax(Bmax),zmin(0),zmax(1),
    windowdump(false)
{
    set(p);
}
void OneWindow::set(ThePlot *p)
{
  //ffassert(p);
    bool first = !theplot;
    bool change = theplot != p;
    theplot=p;
    //    p->win=this;
    if(first) DefaultView() ;
    
    
}

void OneWindow::DefaultView() 
{
    if(theplot)
      {
	  R2 A(theplot->Pmin),B(theplot->Pmax);
	  R2 D(A,B);
	  D *=0.05;
	  zmax = theplot->fmax;
	  zmin = theplot->fmin;
	  
	  A -= D;
	  B += D;
	  
	  if (theplot->aspectratio)
	      cadreortho(A,B);
	  else 
	      cadre(A,B);
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
  ShowGlerror("Begin SetView");   
  glDisable(GL_DEPTH_TEST);
  glViewport(0, 0,width, height);
  
  R zzmin = Min(zmin,theplot->fminT);
  R zzmax = Max(zmax,theplot->fmaxT);    
  R dz = (zzmax-zzmin);
  R zm=(zzmin+zzmax)*0.5;
  if((debug > 3)) cout << "\t\t\t   SetView " << this << " " << Bmin  << " " << Bmax << " " << zzmin << " " << zzmax 
		       << "  zm  " << zm << " dz  " << dz << endl;
  ShowGlerror("0 Set MV");
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  ShowGlerror(" Set MV");
  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity(); 
  ShowGlerror(" Set PM 1");
  glOrtho(Bmin.x,Bmax.x,Bmin.y,Bmax.y,-dz,dz);
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
void  OneWindow::resize(int w,int h)
{  double ww=width,hh=height;
    width=w;
    height=h;
    if (theplot->aspectratio)
      {	
	  cadreortho(oBmin,oBmax);
      }
    
}

void  OneWindow::zoom(int w,int h,R coef)
{
    GLdouble x=w,y=h,z=(zmin+zmax)/2.;
    GLdouble xx,yy,zz;
    

    GLint ok= gluUnProject( x,y,z,modelMatrix,projMatrix,viewport,&xx,&yy,&zz);
    ShowGlerror(" UnPro .. ");

    cout << x << " " << y << " " << z << " -> " << xx << " " << yy << " " << zz << endl;
    R2  oD(oBmin,oBmax);
    R2  D(Bmin,Bmax);
    R2 O(xx,yy);// oBmin.x+D.x*xx/width,oBmin.y+D.y*yy/height); 
    if((debug > 3)) cout<< " zoom : "  << this << " O " << O  << " " << coef << " D = "<<  D<< "as "<< theplot->aspectratio <<  endl;
    oD *= 0.5*coef;
    R2 A = O - oD;
    R2 B = O + oD;
    if (theplot->aspectratio)
	cadreortho(A,B);
    else 
	cadre(A,B);
}
void OneWindow::MoveXView(R dx,R dy) 
{}

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
  ShowGlerror("Begin Display");
  
  //  SetView();
  if(theplot)
    theplot->Draw(this);
  ShowGlerror("After Display");
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

OnePlotBorder::OnePlotBorder(PlotStream & f) 
  :OnePlot(4)
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
    //  if(!gllist && plot.Change() )
  {
      if(!gllist) gllist= plot.gllist++;
      glNewList(gllist,GL_COMPILE_AND_EXECUTE); // save  la list sans affichage
      Draw(win);
      glEndList();  // fin de la list
  }
    //  else glCallList(gllist);
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
  :  count(kcount), state(0),gllist(1),
     changeViso(true),changeVarrow(true),changeColor(true),
     changeBorder(true),changeFill(true), withiso(false),witharrow(false)
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
  Niso=20;
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
  Pmin=R2(+dinfty,+dinfty);
  fmin = +dinfty;    
  fmax = -dinfty;
  Pmax=R2(-dinfty,-dinfty);
  vmax=0;
  
  coefr=1;
  
  long cas;
  
  while(1)
    {
      fin >> cas;
      if((debug > 4)) cout << " read cas: " << cas << " " << PlotStream::dt_endarg << endl;
      
      if(cas==PlotStream::dt_endarg) break;
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
      default: 
	cout << " cas : " << cas <<endl;
	ffassert(0);
	break;
      }
      ffassert(fin.good());
    }
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
  Ths.resize(nbmeshes);
  for(int i=0;i<nbmeshes;++i)
    Ths[i]=0;
  for(int i=0;i<nbmeshes;++i)
    { 
      long l;
      fin >> l;
      if((debug > 3)) cout << " read mesh " << i  << " -> " << l << "  " <<nbmeshes << endl;
      l--;
      ffassert(l>=0 && l < nbmeshes);
      ffassert(Ths[l]==0);
      fin >>Ths[l] ;
      if((debug > 3))
	cout << i << " nt/nv " << l << " "  <<Ths[l]->nt << " " << Ths[l]->nv << endl;
      ffassert(fin.good());
    }
  
  fin.GetPlots();
  long nbplot;
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
      if(what==0) 
	{ 
	  
	  fin >> imsh;
	  p=new OnePlotMesh(Ths[imsh-1]);
	}
      else if (what==1 || what==2)
	{
	  fin >> imsh;
	  if(what==1) withiso=true;
	  else if (what==2) witharrow=true;
	  if((debug > 10)) cout << " plot : mesh " << imsh << endl;
	  ffassert(imsh>0 && imsh <=nbmeshes);
	  p=new OnePlotFE(Ths[imsh-1],what,fin);
	}
      else if(what==3)
	p=new OnePlotCurve(fin);
      else if(what==4)
	p=new OnePlotBorder(fin);
      else
	ffassert(0);
      ffassert(p);
      plots.push_back(p);
      p->bb(Pmin,Pmax);
      p->bfv(fmin,fmax,vmax);
      ffassert(fin.good());		      
    }
  //cout << "\t\t\t\t  f min, max v max :" << fmin << " " << fmax << " " << vmax << endl;
    
  if(fmax < fmin)
    {
      fmax = 1;
      fmin = 0;
    }
    
  PminT=Pmin;
  PmaxT=Pmax;
  fminT=fmin;
  fmaxT=fmax;
  if(old)
    {
      Pmin= Minc(Pmin,old->PminT);
      Pmax= Maxc(Pmax,old->PmaxT);
      fmax= Max(fmax,old->fmaxT);
      fmin= Min(fmin,old->fminT);
    }
  
  z0= fminT +(fmaxT-fminT)*0.01;
  if((debug > 2)) cout << "               data bound: " << PminT << " " << PmaxT  << " fmin == " << fminT << "  " << fmaxT << " z0 " << z0 <<  endl;
  fin.GetEndPlot(); 
  Viso.resize(Niso);
  Varrow.resize(Narrow);
  
  SetColorTable(Max(Niso,Narrow)+4) ;           
  
}


void ThePlot::SetDefIsoV()
{
  R d = fill ? (fmaxT-fminT)/(Niso-2)  : (fmaxT-fminT)/(Niso-1);       
  R x = fill ? fminT-d/2 :fminT+d/2;
  if(!pViso)
    {
      for (int i = 0;i < Niso;i++)
	{Viso[i]=x;x +=d; }
      if (fill ) {Viso[0]=fminT-d;Viso[Niso-1]=fmaxT+d;}    
    }
  x=0; 
  if(debug>10)
    cout << "vmax=  " << vmax << endl;
  d= sqrt(vmax)/(Narrow-1.1);   
  if (!pVarrow)
    for (int i = 0;i < Narrow;i++)
      {Varrow[i]=x;x +=d; }
  SetColorTable(Max(Niso,Narrow)+4) ; 
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
void  plot(double xx,double yy,const char *cmm)
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
    //void * glut_font=GLUT_BITMAP_8_BY_13;    
    for (const char *s=cmm; *s; s++)
      {if((debug > 10)) cout << *s ;
	glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10,*s);
      }
    if((debug > 10)) cout << " ;;; " <<endl;
}
void plot(double x,double y,double i)
{
    char buf[24];
    snprintf(buf,24,"%g",i);
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
		glBegin(GL_POLYGON);
		//SetColor((xfb+xfh)/2); 
		for (int i=0;i<im;i++)
		  {// if((debug > 10)) cout << i << " \t : " << PQ[i].x << " " <<  PQ[i].y << " " << z[i]*rapz << endl;
		      glVertex3f(PQ[i].x, PQ[i].y,z[i]*rapz);
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
    sprintf(fname,"L_%04d.ppm",counter);
    if ((fptr = fopen(fname,"w")) == NULL) {
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

}

static void Mouse( int button,int state,int x,int y )
{
    // state up or down 
    OneWindow * win=CurrentWin();
    if(win && state == GLUT_DOWN) { win->xold=x,win->yold=y;return;}
    // if((debug > 10)) cout << "Mouse " << button<< " " << state << " " << x-xold << " " << y-yold << endl;
    //  x gauche -> droitre
    //  y  haut -> bas`
    glutPostRedisplay();
    
}
static void MotionMouse(int x,int y )
{
    OneWindow * win=CurrentWin();
    if(win)
      {
	  win->xold=x;
	  win->yold=y;
	  glutPostRedisplay();
      }
}

static void Key( unsigned char key, int x, int y )
{
    OneWindow * win=CurrentWin();
    switch (key) {
	case 27: // esc char
	    Fin(0);
	    break;
	case 'w':
	    if(win)
		win->windowdump=true;
	    break;
	case '+':
	    win->zoom(x,y,0.7);
	    break;
	case '-':
	    
	    win->zoom(x,y,1./0.7);
	    break;
	case '=':
	    win->DefaultView();
	    break;
	case 'f':
	    win->theplot->fill = ! win->theplot->fill  ;
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
	    
	case 'p':
	    
	    break;
	case 'a':
	    win->theplot->coeff/= 1.2;
	    break;
	case 'A':
	    win->theplot->coeff*= 1.2;
	    break;
	    
	case '\r':
	case '\n':
	  SendForNextPlot();
	    break;
	default:
	    if((debug > 10)) cout << " Key Character " << (int) key << " " << key << endl;  
	    
    }
    glutPostRedisplay();
}


void SpecialKey(int key, int x, int y)
{
    OneWindow * win=CurrentWin();
    
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
  if(debug)
    cout << " We Read a plot  : " << kread << " " << nextPlot << " " << err << endl;
  if(err<0)
    NoMorePlot=true; 
  Thread::Exit();
}


int main(int argc,  char** argv)
{

  try {
    datafile =stdin;
    if(argc>1 && *argv[argc-1] != '-' ) datafile=fopen(argv[argc-1],"r");
    ffassert(datafile);
    int err=ReadOnePlot(datafile);
    if(err) throw string(" no plot in file? ");
    
    bool stereo=false;
    bool fullscreen = false;


    if(kread==0) {
      cout << " Error: no graphic data " << endl; 
      Fin(1);
    }
    
    cout << "on a lue le premier plot next plot: " << nextPlot << endl;

    glutInit(&argc, argv);

    if(stereo)
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO);
    else  
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
    int Height = 512;
    int Width = 512*3/2; 
    
    glutInitWindowSize(Width , Height);
    glutInitWindowPosition(100, 100);

    string titre = "GLUT/ FreeFem++ ";
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
  }
  catch (const StringErr & s)
    {
      cout << " catch error " << s << endl;
      Fin(2);
    }

  return 0;
}


