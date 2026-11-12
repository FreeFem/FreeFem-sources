/*
 *  ffglut.hpp
 *  ff
 *
 *  Created by FrÈdÈric Hecht on 04/11/08.
 *  Copyright 2008 UPMC. All rights reserved.
 *
 */
#include <complex>
typedef std::complex<double> Complex;
extern int debug;
static const double dinfty=numeric_limits<double>::max();
void DefColor(float & r, float & g, float & b,
              int k,int nb, bool hsv,bool grey,KN<R> colors);
void hsvToRgb (float h, float s, float v, float & r, float & g, float & b);


class ThePlot;
class OneWindow;

//  a grid of buckets to find quickly the element which contains a point (value under the mouse)
struct FFGrid {
    R3 a,h;
    int n[3];
    vector<int> head,nxt,elem;
    bool built;
    FFGrid() : a(0.,0.,0.), h(1.,1.,1.), built(false) { n[0]=n[1]=n[2]=1;}
    int cell(int i,int j,int k) const { return (k*n[1]+j)*n[0]+i;}
    void ijk(const R3 & P,int c[3]) const
    {
        R p[3]={(P.x-a.x)/h.x,(P.y-a.y)/h.y,(P.z-a.z)/h.z};
        for(int d=0;d<3;++d) c[d]=max(0,min(n[d]-1,(int) floor(p[d])));
    }
    void build(const vector<R3> & bmn,const vector<R3> & bmx)
    { //  the element e is in all the cells of its bounding box [bmn[e],bmx[e]]
        built=true;
        int ne=bmn.size();
        if(!ne) return;
        R3 pmn=bmn[0],pmx=bmx[0];
        for(int e=1;e<ne;++e) { pmn=Minc(pmn,bmn[e]); pmx=Maxc(pmx,bmx[e]);}
        R D[3]={pmx.x-pmn.x,pmx.y-pmn.y,pmx.z-pmn.z};
        R Dm=max(D[0],max(D[1],D[2]));
        int dim=0;
        R vol=1;
        for(int d=0;d<3;++d) if(D[d]>1e-9*Dm) dim++, vol*=D[d];
        R side = dim ? pow(vol/max(1.,ne/2.),1./dim) : 1.;
        R hh[3];
        for(int d=0;d<3;++d)
        {
            n[d] = (D[d]>1e-9*Dm) ? max(1,min(2000,(int) ceil(D[d]/side))) : 1;
            hh[d] = (D[d]>1e-9*Dm) ? D[d]/n[d] : 1.;
        }
        a=pmn;
        h=R3(hh[0],hh[1],hh[2]);
        head.assign(n[0]*n[1]*n[2],-1);
        for(int e=0;e<ne;++e)
        {
            int c0[3],c1[3];
            ijk(bmn[e],c0);
            ijk(bmx[e],c1);
            for(int k=c0[2];k<=c1[2];++k)
                for(int j=c0[1];j<=c1[1];++j)
                    for(int i=c0[0];i<=c1[0];++i)
                    {
                        int c=cell(i,j,k);
                        elem.push_back(e);
                        nxt.push_back(head[c]);
                        head[c]=elem.size()-1;
                    }
        }
    }
};

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
    int num; // number of the curve in the 1d graph (color)
    // true if the plot can be drawn in a 1d graph a la gnuplot, and add the bound of the graph
    //  lx, ly: log scale in x, y (then the bound are the log10 of the positive values)
    virtual bool Graph1d(R2 & Pmn,R2 & Pmx,bool lx=false,bool ly=false) const {return false;}
    virtual void Draw(OneWindow *win) =0;
    virtual bool NextCase() {return false;}
    //  the value of the function at the point P (under the mouse): ncomp = 1 (val) or 2, 3 (vec, val = norm)
    virtual bool Probe(const R3 & P,R & val,R3 & vec,int & ncomp) {return false;}
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
    oklist(nbgllist),setgllists(0),num(0){
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
    //  cut plane: the section of the elements
    bool cutok;
    R3 lastcutN;
    R lastcutc;
    void DrawSection(ThePlot & plot,R3 N,R c);
    OnePlotMesh3(const Mesh3 *T)
    : OnePlot(0,3,3),Th(T),cutok(false),lastcutc(0)
    {
        Pmin=Th->Pmin;
        Pmax=Th->Pmax;
        //Th->BoundingBox(Pmin,Pmax);
    }
    void Draw(OneWindow *win);
};

struct OnePlotMeshS : public OnePlot
{
    const MeshS *Th;
    OnePlotMeshS(const MeshS *T)
    : OnePlot(0,3,3),Th(T)
    {
        Pmin=Th->Pmin;
        Pmax=Th->Pmax;
        //Th->BoundingBox(Pmin,Pmax);
    }
    void Draw(OneWindow *win);
};

struct OnePlotMeshL : public OnePlot
{
    const MeshL *Th;
    OnePlotMeshL(const MeshL *T)
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
    KN<Complex> vc;
    KN<R2> Psub;
    KN<int> Ksub;
    int cas; // in cas of complex  chaage interpertation of complex value
    OnePlotFE(const Mesh *T,long w,PlotStream & f);
    void Draw(OneWindow *win);
    void DrawPretty(OneWindow *win,bool elev);
    FFGrid grid;
    bool Probe(const R3 & P,R & val,R3 & vec,int & ncomp);
    //  pretty rendering: the values of the iso lines (keys n, N), marks on the color bar
    vector<R> isovals;
    bool showisovals;
    void DrawMeshOnSurface(ThePlot & plot,bool border);
    void dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn,R & vmx) const ;
    bool  vc2v();
    bool  NextCase() { cas++; return vc2v();}
};

struct OnePlotFE3: public OnePlot
{
    const Mesh3 *Th;
    long nsub;
    KN<double> v;
    KN<Complex> vc;
    KN<R3> Psub;
    KN<int> Ksub;
    int cas; // in cas of complex  chaage interpertation of complex value
    
    // the boundary faces of the mesh ( 4*k+e : face e of the tet k) to plot the function on the boundary
    vector<int> bfaces;
    vector<int> bedges;// the sharp edges of the boundary (2 vertices by edge), to see the shape of the domain
    //  pretty rendering: the gradient of the function at the vertices (normal of the isosurfaces)
    vector<R3> gradv;
    void ComputeGradV();
    void DrawBoundaryPretty(OneWindow *win,bool plotmesh);
    void DrawBoundaryMesh();
    FFGrid grid;
    bool Probe(const R3 & P,R & val,R3 & vec,int & ncomp);
    bool ProbeRay(const R3 & A,const R3 & B,OneWindow *win,R3 & P,R & val);
    bool bfacesok;
    bool bconstant;// the function is (nearly) constant on the boundary => isosurfaces by default
    R bval;// the value on the boundary if bconstant
    //  pretty rendering: the values of the isosurfaces (marks on the color bar), shown or not
    vector<R> isovals;
    bool showisovals;
    //  p0: the function is constant on each element => one color by element (no iso line, no isosurface),
    //  and the isosurfaces are replaced by the elements with a value greater than a threshold
    bool p0;
    R lastthreshold;
    //  exploded view of a P0 function: the elements with the same value are moved away from the center
    //  (each element if there are too many values), expdir = the direction of the move of each element
    vector<R3> expdir;
    R lastexplode;
    void BuildExplode();
    void DrawThreshold(OneWindow *win,R thr,bool plotmesh,R explode);
    void BoundaryConstant();
    void DrawBoundaryEdges(ThePlot & plot);
    //  cut plane (N,P) = c : the function on the section of the mesh
    int lastcutaxis,lastcutmode;
    R lastcutpos;
    R3 lastcutN;
    void DrawSection(OneWindow *win,R3 N,R c,bool plotmesh,R thr);
    int lastmode;// to rebuild the graphic list if the kind of plot change
    void BuildBoundaryFaces();
    void DrawBoundary(OneWindow *win,bool plotmesh);
    void DrawIsoSurfaces(OneWindow *win);
    void DrawArrows(OneWindow *win);

    OnePlotFE3(const Mesh3 *T,long w,PlotStream & f)
    :OnePlot(w,3,5),Th(T),cas(2),bfacesok(false),bconstant(false),bval(0),showisovals(false),p0(false),lastthreshold(0),lastexplode(0),lastmode(-1),lastcutaxis(0),lastcutmode(-1),lastcutpos(-1)
    {
        Pmin=Th->Pmin;
        Pmax=Th->Pmax;
        
        f >> Psub ;
        f >> Ksub ;
        if(what<10)
            f >>  v;
        else
            f >> vc;
        vc2v();
        if(debug>3) cout << "OnePlotFE3" << Th <<" " << what<< " " << nsub <<" " << v.N() << endl
            << "       Pmin " << Pmin << " Pmax  " << Pmax << endl;
        ffassert(f.good());
        
    }
    void Draw(OneWindow *win);
    bool  vc2v();
    bool  NextCase() { cas++; return vc2v();}
};


struct OnePlotFES: public OnePlot
{
    const MeshS *Th;
    long nsub;
    KN<double> v;
    KN<Complex> vc;
    KN<R2> Psub;
    KN<int> Ksub;
    int cas; // in cas of complex  chaage interpertation of complex value
    
    FFGrid grid;
    bool Probe(const R3 & P,R & val,R3 & vec,int & ncomp);
    OnePlotFES(const MeshS *T,long w,PlotStream & f)
    :OnePlot(w,3,5),Th(T),cas(2)
    {
        Pmin=Th->Pmin;
        Pmax=Th->Pmax;
        
        f >> Psub ;
        f >> Ksub ;
        if(what<10)
            f >>  v;
        else
            f >> vc;
        vc2v();
        if(debug>3) cout << "OnePlotFES :" << Th <<" " << what<< " " << nsub <<" " << v.N() << endl
            << "       Pmin " << Pmin << " Pmax  " << Pmax << endl;
        ffassert(f.good());
        
    }
    void Draw(OneWindow *win);
    bool  vc2v();
    bool  NextCase() { cas++; return vc2v();}
};

struct OnePlotFEL: public OnePlot
{
    const MeshL *Th;
    long nsub;
    KN<double> v;
    KN<Complex> vc;
    KN<R1> Psub;
    KN<int> Ksub;
    int cas; // in cas of complex  chaage interpertation of complex value
    
    OnePlotFEL(const MeshL *T,long w,PlotStream & f)
    :OnePlot(w,3,5),Th(T),cas(2)
    {
        Pmin=Th->Pmin;
        Pmax=Th->Pmax;
        
        f >> Psub ;
        f >> Ksub ;
        if(what<16)
            f >>  v;
        else
            f >> vc;
        vc2v();
        if(debug>3) cout << "OnePlotFEL :" << Th <<" " << what<< " " << nsub <<" " << v.N() << endl
            << "       Pmin " << Pmin << " Pmax  " << Pmax << endl;
        ffassert(f.good());
        
    }
    void Draw(OneWindow *win);
    void DrawGraph(OneWindow *win);
    bool Graph1d(R2 & Pmn,R2 & Pmx,bool lx=false,bool ly=false) const;
    bool  vc2v();
    bool  NextCase() { cas++; return vc2v();}
};

struct OnePlotCurve: public OnePlot {
    KN<double> xx,yy,zz,cc;
    int cas;
    bool Graph1d(R2 & Pmn,R2 & Pmx,bool lx=false,bool ly=false) const;
    OnePlotCurve(PlotStream & f,int nfield=2,ThePlot *theplot=0);
    void Draw(OneWindow *win);
    void  dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn2,R & vmx2) const;
     bool  NextCase() { cas++; return false;}
};

struct OnePlotBorder: public OnePlot {
    vector<vector<pair<long,R3> > > data;
    OnePlotBorder(PlotStream & f);
    void Draw(OneWindow *win);
};


struct OnePlotHMatrix: public OnePlot
{
    int si;
    int sj;
    long nbdense;
    long nblr;
    std::vector<std::pair<int,int>> offsetsdense;
    std::vector<std::pair<int,int>> sizesdense;
    std::vector<std::pair<int,int>> offsetslr;
    std::vector<std::pair<int,int>> sizeslr;
    std::vector<int> rankslr;
    std::vector<double> compression;

    OnePlotHMatrix(long w, PlotStream & f)
    :OnePlot(w)
    {
        //dim = 2;
        Pmin = R3(0,0,0);
        Pmax = R3(1,1,0);

        int offseti, offsetj;
        int sizei, sizej;
        int rank;

        f >> si;
        f >> sj;
        f >> nbdense;
        f >> nblr;

        offsetsdense.resize(nbdense);
        sizesdense.resize(nbdense);
        offsetslr.resize(nblr);
        sizeslr.resize(nblr);
        rankslr.resize(nblr);
        compression.resize(nblr);

        for (int i=0;i<nbdense;i++) {
          f >> offseti;
          f >> offsetj;
          f >> sizei;
          f >> sizej;
          offsetsdense[i] = std::pair<int,int>(offseti,offsetj);
          sizesdense[i] = std::pair<int,int>(sizei,sizej);
        }

        for (int i=0;i<nblr;i++) {
          f >> offseti;
          f >> offsetj;
          f >> sizei;
          f >> sizej;
          f >> rankslr[i];
          f >> compression[i];
          offsetslr[i] = std::pair<int,int>(offseti,offsetj);
          sizeslr[i] = std::pair<int,int>(sizei,sizej);
        }

        ffassert(f.good());
    }
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



class ThePlot
{
public:
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
    int nbN; // intensity for plotting surface normal
    bool bw;
    string * psfile;
    string * cm;
    
    bool grey;
    bool greyo;
    bool drawborder;
    bool drawmeshes;
    bool add,keepPV,pNormalT;
    int periodNormalT;
    double echelle;
    double ArrowSize;
    vector<Mesh *> Ths;
    vector<Mesh2 *> Ths2;
    vector<Mesh3 *> Ths3;
    vector<MeshS *> ThsS;
    vector<MeshL *> ThsL;
    list<OnePlot *> plots;
    bool changeViso,changeVarrow,changeColor,changeBorder,changeFill;
    R3 Pvue,Peyes;
    R alpha;
    R coefr;
    R z0; //  z pour les objets 2d.
    R ZScale;
    //  for 3d plot jan 2009
    int  plotdim;
    bool blockwin, plotNormalT, changePlotdim;
    R theta, phi, dcoef, focal;
    int datadim;
    // 1d graph a la gnuplot (true if the plot contains only 2d curves
    //  or scalar FE functions on a meshL along the x axis)
    bool graph1d;
    int nbcurves;
    bool GraphBound(R2 & Pmn,R2 & Pmx,bool lx=false,bool ly=false) const;
    //  options of the 1d graph given in the comment (cmm=) of the plot:
    //    cmm="title | xlabel=... | ylabel=... | legend=name1,name2,... | logx | logy"
    string gtitle,gxlabel,gylabel;
    vector<string> glegend;
    int glogx,glogy;// -1 not given, 1 log scale
    void GraphOptions();
    // 2D
    long winnum;
    bool NextCase();
    
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
        for (vector<MeshS *>::iterator i= ThsS.begin();i != ThsS.end(); ++i)
            if(*i) delete *i;
        for (vector<MeshL *>::iterator i= ThsL.begin();i != ThsL.end(); ++i)
            if(*i) delete *i;
        
    }
    ThePlot(PlotStream & fin,ThePlot *old , int kcount);
    
    void Draw(OneWindow *win) ;
    void DrawHelp(OneWindow *win) ;
    void DrawHelpAll(OneWindow *win) ;
    void PrettyPlots(OneWindow *win,bool & prettyiso,bool & prettyarrow) const;
    
    struct RGB  { float r,g,b;
        void set(int k,int nb,ThePlot *theplot )
        {
            DefColor(r,g,b, k,nb,theplot->hsv,theplot->grey,theplot->colors);
        }
    } ;
    vector<RGB>  tbc;
    void color(int i,R alpha=1.) {
        ffassert(tbc.size());
        RGB c(tbc[min(max(0,i),(const int) tbc.size()-1)]);
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
    //  pretty rendering: the colors of the iso lines in the colormap (cmap >=0), see colorv
    int cmap;
    R cmapfmn,cmapfmx;
    void colorv(int i,R val);
    void SetDefIsoV(int niso=0,int narr=0,R fmn=1.,R fmx=-1.,R vmn=1.,R vmx=-1.);
    void DrawIsoT(const R2 Pt[3],const R ff[3],const R * Viso,int NbIso, R rapz=1);
    void DrawIsoTfill(const R2 Pt[3],const R ff[3],const R * Viso,int NbIso, R rapz=1);
    void DrawIsoT(const R3 Pt[3],const R ff[3],const R3 Nt[3],const R * Viso,int NbIso,bool changePlotdim, int viewdim,R rapz=1);
    void DrawIsoTfill(const R3 Pt[3],const R ff[3],const R3 Nt[3],const R * Viso,int NbIso,bool changePlotdim,int viewdim, R rapz=1);
    void DrawIsoEfill(const R3 Pt[2],const R ff[2],const R3 Nt[2],const R * Viso,int NbIso,bool changePlotdim,int viewdim, R rapz=1);
    void dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn,R & vmx) const ;
    
};

class OneWindow {
public:
    int countdisplay;// for bug mojave 
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
    bool windowdump,windowdumppng,help;
    bool helpall;// ? once: the help of the current plot, ? twice: all the keys
    int showhx;// spacing of the lines of Show (0: default)
    
    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
    GLint viewport[4];
    
    int  plotdim, viewdim;
    R theta, phi, coef_dist, focal, dtheta;
    R  rapz,rapz0;
    R3 Bmin3,Bmax3,Pvue3;
    R3 cam;
    bool withlight;
    bool changearrow,changeiso;// to rebuild de graphic list if neccessary
    bool keepPV,init,pNormalT,changePlotdim;
    // 1d graph a la gnuplot
    bool graphgrid,graphkey;// show the grid, the key
    //  pretty rendering of the 3d plots (key R): continuous colors (texture of the colormap),
    //  smooth normals, specular light, gradient background, color bar (key v), colormap (key C)
    bool pretty,prettybar;
    int cmapkind;// 0 turbo, 1 viridis, 2 cool to warm, 3 rainbow of FreeFEM, 4 grey
    GLuint cmaptex;
    int cmaptexkind;
    bool Pretty3d() const { return pretty && plotdim==3 && !Graph1d();}
    bool PrettyOn() const { return pretty && !Graph1d();}// 2d and 3d
    bool probe;// show the value of the function under the mouse (key P)
    bool isolines;// pretty rendering 2d with fill: the iso lines on the colors (key o)
    bool axes;// the axes x, y (and z in 3d) with the tics and the values (key j)
    void DrawAxes();
    void DrawAxis(const R3 & A,const R3 & B,R vmin,R vmax,const R2 & C,R gap,const char *name);
    void DrawProbe();
    GLuint ColormapTexture();
    void DrawBackground();
    void DrawColorBar(R fmn,R fmx,R x0,const char *title,const vector<R> *marks=0);
    int graphstyle;// 0 lines, 1 lines and points, 2 points
    bool graphlogx,graphlogy;// log scale
    int userlogx,userlogy;// log scale asked with the keys x, y (-1: no, so the option of the plot)
    //  cut plane of the 3d function: cutaxis = 0 no cut, 1,2,3 plane orthogonal to the axis x,y,z at the beginning
    //  cutN = the normal of the plane (can be rotated)
    //  cutpos = position of the plane in the bounding box of the mesh (between 0 and 1)
    int cutaxis;
    R cutpos;
    R3 cutN;
    void RotateCut(R ah,R av);
    //  the keys of medit: F1 cut plane on/off, F2 edit (the left button rotates the plane),
    //  F3 freeze (the plane is fixed on the object: the mouse turns the object and the plane)
    bool cutedit,cutfreeze;
    R explode;// 3d P0 function: coefficient of the exploded view (0 no)
    bool SetCutPlane(const R3 & Pmn,const R3 & Pmx,R3 & N,R & c);
    GLdouble cuteq[4];// the equation of the clip plane (kept side: >= 0), when cutaxis
    void zoomat(int w,int h,R coef);
    R threshold;// 3d P0 function: plot the elements with a value > min + threshold*(max-min)
    bool isosurf;// 3d scalar function: switch between the isosurfaces and the function on the boundary of the mesh
    int xmouse,ymouse;// mouse position in pixel (to show the coordinates)
    bool Graph1d() const { return theplot && theplot->graph1d && plotdim==2;}
    void GraphBox(int &x0,int &y0,int &x1,int &y1) const;// the box of the axis in pixel
    void DrawGraph();
    //double  aspx, aspy, echx,echy,ech,rxmin,rxmax,rymin,rymax;
    OneWindow(int h,int w,ThePlot *p);
    void DefaultView(int state) ;
    bool  NextCase() ;
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
    void Seg3(R3 A, R3 B) const  {
        glVertex3d(A.x,A.y,A.z);
        glVertex3d(B.x,B.y,B.z);
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

