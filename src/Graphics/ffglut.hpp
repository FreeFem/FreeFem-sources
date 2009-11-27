/*
 *  ffglut.hpp
 *  ff
 *
 *  Created by FrÈdÈric Hecht on 04/11/08.
 *  Copyright 2008 UPMC. All rights reserved.
 *
 */
extern int debug;
static const double dinfty=numeric_limits<double>::max();
void DefColor(float & r, float & g, float & b,
              int k,int nb, bool hsv,bool grey,KN<R> colors);
void hsvToRgb (float h, float s, float v, float & r, float & g, float & b);


class ThePlot;
class OneWindow;

struct OnePlot 
{
  int dim;
  R3 Pmin,Pmax;
  double fmin,fmax;
  double vmax2;
  long what;
  GLsizei ngllists;
  GLint gllists;
  KN<int> oklist;
  int setgllists;
  virtual void Draw(OneWindow *win) =0;
  void bb(R3 & Pmn,R3 &Pmx) const 
  { 
    Pmn=Minc(Pmin,Pmn);
    Pmx=Maxc(Pmax,Pmx);
    }

  void bfv(R & fmn,R &fmx,R & vmx) const 
  { 
    // cout << "\t\t\t\t  f min, max v max :" << fmin << " " << fmax << " " << vmax2 << endl;
    fmn=Min(fmin,fmn);
    fmx=Max(fmax,fmx);
    vmx=Max(vmax2,vmx);
  }

  virtual void dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn,R & vmx) const 
  {//  compute the function bound and arrow bound view ....  
  }
  
  OnePlot(long w,int ddim=2,int nbgllist=0) :
    dim(ddim),
    Pmin(dinfty,dinfty,dinfty),Pmax(-dinfty,-dinfty,-dinfty),
    fmin(dinfty),fmax(-dinfty),vmax2(0),
    what(w),ngllists(nbgllist),gllists(0),
    oklist(nbgllist),setgllists(0){
  }
  
  void initlist()
  {
    if(! setgllists)
      if(ngllists>0) 
	{
	  oklist=0; //  
	  gllists= glGenLists(ngllists);
	  setgllists=1;
	  assert(gllists);
	}
  }
  
  void GLDraw(OneWindow *win);
  virtual ~OnePlot() {
    if(setgllists) glDeleteLists(gllists,ngllists);
  };
private: // pas de copy  car il y a des destructeurs dans les classes derives
  OnePlot(const OnePlot & );
  OnePlot & operator=(const OnePlot & );
};

template<class Mesh>
struct OnePlotMesh : public OnePlot
{
  const Mesh *Th;
  OnePlotMesh(const Mesh *T)
    : OnePlot(0,2,3),Th(T) 
  {
    R2 P0,P1;
    Th->BoundingBox(P0,P1);
    Pmin=P0;
    Pmax=P1;
  }
  void Draw(OneWindow *win);
  
};

struct OnePlotMesh3 : public OnePlot
{
    const Mesh3 *Th;
    OnePlotMesh3(const Mesh3 *T)
      : OnePlot(0,3,3),Th(T) 
    {
	Pmin=Th->Pmin;
	Pmax=Th->Pmax;
	//Th->BoundingBox(Pmin,Pmax);     
    }
    void Draw(OneWindow *win);
};

template<class Mesh>
struct OnePlotFE: public OnePlot 
{
  const Mesh *Th;
 // long nsub;
  KN<double> v;
  KN<R2> Psub;
  KN<int> Ksub;
    
    OnePlotFE(const Mesh *T,long w,PlotStream & f);
    void Draw(OneWindow *win);
    void dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn,R & vmx) const ;
};

struct OnePlotFE3: public OnePlot 
{
    const Mesh3 *Th;
    long nsub;
    KN<double> v;
    KN<R3> Psub;
    KN<int> Ksub;
    OnePlotFE3(const Mesh3 *T,long w,PlotStream & f)
      :OnePlot(w,3,5),Th(T)
    {
	Pmin=Th->Pmin;
	Pmax=Th->Pmax;

	f >> Psub ;
	f >> Ksub ;
	f >>  v;
	if(what==6)
	  {
	      fmin = min(fmin,v.min());
	      fmax = max(fmax,v.max());
	  }
	else if (what==7)
	  {  
	      ffassert(0); // afaire
	      int n= v.N()/2;
	      for (int i=0,j=0;i<n;i++, j+=2)
		{
		    R2 u(v[j],v[j+1]);
		    vmax2 = max(vmax2,u.norme2());
		}
	      //cout << " vmax = " << sqrt(vmax2) << endl; 
	  }
	if(debug>3) cout << "OnePlotFE3" << Th <<" " << what<< " " << nsub <<" " << v.N() << endl
			 << "       Pmin " << Pmin << " Pmax  " << Pmax << endl;
	ffassert(f.good());
	
    }
    void Draw(OneWindow *win);
    
};

struct OnePlotCurve: public OnePlot {
  KN<double> xx,yy;
  OnePlotCurve(PlotStream & f)
    :OnePlot(3,2,1)
  {
    f >> xx>>yy;
    // cout << xx << " " << yy <<endl;
    ffassert(f.good());
    ffassert(xx.N() && yy.N() && xx.N() == yy.N());
    Pmin=Minc(Pmin,R2(xx.min(),yy.min()));
    Pmax=Maxc(Pmax,R2(xx.max(),yy.max()));
    
  }
  void Draw(OneWindow *win);
};

struct OnePlotBorder: public OnePlot {
  vector<vector<pair<long,R2> > > data;
  OnePlotBorder(PlotStream & f);	  
  void Draw(OneWindow *win);
};

// add 11/12/2008 for gestion of error FH  (what -1) 
struct OnePlotError: public OnePlot {
  long item;
  OnePlotError(PlotStream & f) 
    : OnePlot(-1)
  {    f >> item; 
     }	  
  void Draw(OneWindow *win) ;
  
};



class ThePlot { public:
    int count;
    int state; // 0 new 
    GLint gllist;
    KN<R> colors;
    bool hsv; // hsv  type 
    KN<R> boundingbox;
    double coeff;
    bool wait;
    bool value;
    bool fill;
    bool aspectratio;
    bool clean;
    bool uaspectratio;
    bool pViso,pVarrow;
    bool withiso;
    bool witharrow;
    
    long  Niso,Narrow;
    R3 Pmin,Pmax,PminT,PmaxT;//  with R -> true bound
    R  fmin,fmax,fminT,fmaxT; // with  bound with previous plot. 
    R  vmax2;
    KN<R> Viso,Varrow;
    bool bw;
    string * psfile;
    string * cm;
    
    bool grey;
    bool greyo;
    bool drawborder;
    bool drawmeshes;
    bool add,keepPV;
    double echelle;
    vector<Mesh *> Ths;
    vector<Mesh2 *> Ths2;
    vector<Mesh3 *> Ths3;
    list<OnePlot *> plots;
    bool changeViso,changeVarrow,changeColor,changeBorder,changeFill;
    R3 Pvue,Peyes;
    R alpha; 
    R coefr; 
  R z0; //  z pour les objets 2d. 
  //  for 3d plot jan 2009
  int  plotdim;
  R theta, phi, dcoef, focal;
  int datadim;
    // 2D
    
    bool Change() const  { return changeViso||changeVarrow||changeColor||changeBorder||changeFill;}
    ~ThePlot()
    {
	for (list<OnePlot *>::iterator i= plots.begin();i != plots.end(); ++i)
	    if(*i) delete *i;
	for (vector<Mesh *>::iterator i= Ths.begin();i != Ths.end(); ++i)
	    if(*i) delete *i;
	for (vector<Mesh2 *>::iterator i= Ths2.begin();i != Ths2.end(); ++i)
	    if(*i) delete *i;
	for (vector<Mesh3 *>::iterator i= Ths3.begin();i != Ths3.end(); ++i)
	    if(*i) delete *i;

    }
  ThePlot(PlotStream & fin,ThePlot *old , int kcount);
    
    void Draw(OneWindow *win) ;
    void DrawHelp(OneWindow *win) ;
 
    struct RGB  { float r,g,b;
	void set(int k,int nb,ThePlot *theplot ) 
	{
	    DefColor(r,g,b, k,nb,theplot->hsv,theplot->grey,theplot->colors);
	}
    } ;
    vector<RGB>  tbc;
    void color(int i,R alpha=1.) { 
      ffassert(tbc.size());
	RGB c(tbc[min(max(0,i),(const int) tbc.size())]);
	if(alpha<1) {
	    R a1=1.-alpha;
	    c.r=c.r*alpha+a1;
	    c.g=c.g*alpha+a1;
	    c.b=c.b*alpha+a1;
	    
	 // cout << " aaaa" << alpha << endl;
	}
      glColor4f(c.r,c.g,c.b,alpha);
    }
  void  SetColorTable(int nb); 
  void SetDefIsoV(int niso=0,int narr=0,R fmn=1.,R fmx=-1.,R vmn=1.,R vmx=-1.); 
  void DrawIsoT(const R2 Pt[3],const R ff[3],const R * Viso,int NbIso, R rapz=1);
  void DrawIsoTfill(const R2 Pt[3],const R ff[3],const R * Viso,int NbIso, R rapz=1);
  void dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn,R & vmx) const ;

}; 

class OneWindow { 
public:
  ThePlot *theplot;
  list<ThePlot*> lplots;//  to store some plot 
  list<ThePlot*>::iterator icurrentPlot;//=lplots.begin();
  int lplotssize;

  int height,width;
  R2 Bmin,Bmax;
  R2 oBmin,oBmax;// orign box 
  R zmin,zmax;
  R hpixel;// taille pixel en coordonne x,y,z 
  int xold,yold;
  bool windowdump,help;

  GLdouble modelMatrix[16];
  GLdouble projMatrix[16];
  GLint viewport[4];  
  
  int  plotdim;
  R theta, phi, coef_dist, focal, dtheta;
  R  rapz,rapz0;     
  R3 Bmin3,Bmax3,Pvue3;
  R3 cam;
  bool withlight;
  bool changearrow,changeiso;// to rebuild de graphic list if neccessary

  //double  aspx, aspy, echx,echy,ech,rxmin,rxmax,rymin,rymax;
  OneWindow(int h,int w,ThePlot *p);
  void DefaultView() ;
  void  SetView() ;
  void MoveXView(R dx,R dy) ;
  void set(ThePlot *p);
  void add(ThePlot *p);
  void cadre(R2 A,R2 B);
  void cadreortho(R2 A,R2 B);
  void getcadre(double &xmin,double &xmax,double &ymin,double &ymax);
  void Display();
  void resize(int w,int h);
  void zoom(int w,int h,R coef);
  void zoom(R coef);
  float GetHeigthFont(){return 10;}
    void color(int i,R alpha=1) {theplot->color(i,alpha);}
  void FillRect(R x0,R y0,R x1,R y1);
  void PlotValue(const KN_<double> & Viso,int  k0,const char * cmm);
  void DrawCommentaire(const char * cm,R x,R y); 
  void SetScreenView() const ;
  void Show(const char *str,int i);
  void   setLighting();
  void unsetLighting();
  void Seg(R2 A, R2 B) const  { 
    glVertex3d(A.x,A.y,theplot->z0);
    glVertex3d(B.x,B.y,theplot->z0);
  }
 
  int InRecScreen(R2 P1,R2 P2) const 
  {  
    R2 Pmn=Minc(P1,P2),Pmx=Maxc(P1,P2);
    return (Pmx.x >= Bmin.x) && (Pmn.x <= Bmax.x) && (Pmx.y >= Bmin.y) && (Pmn.y <= Bmax.y);
  }
  int InRecScreen(R3 P1,R3 P2) const 
  {  
    R3 Pmn=Minc(P1,P2),Pmx=Maxc(P1,P2);
    return (Pmx.x >= viewport[0]) && (Pmn.x <= viewport[2]) && (Pmx.y >= viewport[1]) && (Pmn.y <= viewport[3]);
  }
    

};

void plot(double x,double y,const char *cmm,int font=-1);
void plot(double x,double y,double i,int fint = -1);

extern map<int,OneWindow*> AllWindows; 

inline OneWindow* CurrentWin() 
 {
  OneWindow* w= AllWindows[glutGetWindow()];
 ffassert(w);
 return w;
 }
