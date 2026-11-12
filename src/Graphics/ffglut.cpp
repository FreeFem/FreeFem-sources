#ifdef FREEGLUT__
#define GL_SILENCE_DEPRECATION
#include <GL/glut.h>
#else
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#endif
/*
#if !defined(GLUT_WINDOW_SCALE)
#    define NOTGLUT_WINDOW_SCALE
#    define GLUT_WINDOW_SCALE 199
#endif
*/

constexpr double my_pi = 3.141592653589793238462643383279502884;

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
#include <algorithm>
#include <utility>
#include <unistd.h>

#include "rgraph.hpp"
#include "throwassert.hpp"

#include "RNM.hpp"
#include "fem.hpp"

#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"

#include "PlotStream.hpp"
double kscreenscale = 1; // to solve  bug on MacOS Catalina , wrong value in resize
extern long verbosity;
int glutscreenscale=1;//  High Resolution Explained: Features and Benefit
// add for the gestion of the endianness of the file.
//PlotStream::fBytes PlotStream::zott; //0123;
//PlotStream::hBytes PlotStream::zottffss; //012345678;
// ---- FH

using namespace Fem2D;
using std::numeric_limits;
const R pi=my_pi;//4*atan(1.);
using namespace std;

static int   nbSendForNextPlot=0,nbTimerNextPlot=0;
int debug=0;
int casemouse=0,keyact=0;
double gwait=0;//  no wait in second
#include "ffglut.hpp"

#include "ffthreads.hpp"

int version =0;
int subversion =0;

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

int kread=-1;


map<int,OneWindow *> AllWindows;
map<int,int>  Num2Windows;


int Fin(int code)
{
    WaitNextRead();
    if(!NoMorePlot && debug>2)
        cout << " exit before end  " << endl;
    if(NoMorePlot && !code) exit(NoMorePlot ? 0  : 1);
    if(code) exit(NoMorePlot ? 0  : 1);
    return 0;
}

int   ReadOnePlot0(FILE *fp);
int   ReadOnePlot(FILE *fp)
{ //  the stream is locked once for the reading of the plot (see PlotStream::read, getc_unlocked)
    if(!fp) return -4;
#ifndef _WIN32
    flockfile(fp);
#endif
    int r=ReadOnePlot0(fp);
#ifndef _WIN32
    funlockfile(fp);
#endif
    return r;
}

int   ReadOnePlot0(FILE *fp)
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
    const char *  magic3_1="#!ffglutdata3.1";
    const char *  magic3_2="#!ffglutdata3.2";
    const char *  magic4_0="#!ffglutdata4.0";
    const char *  magic4_1="#!ffglutdata4.1";
    /*

     change version 3.2 add colored curve
     */
    const int lmagic=strlen(magic2);
    char magicxx[32];
    err=0;
    // init ..
    if(kread==-1)
    {
        for(int i=0;i<lmagic;i++)
        { int c=getc(fp);
            if(i==0 && c == EOF) return -2; // empty file
            magicxx[i]=c;
        }
        magicxx[lmagic]='\0';
        if( strcmp(magicxx,magic2)==0)  version=2;
        else if( strcmp(magicxx,magic3)==0)  version=3;
        else if( strcmp(magicxx,magic3_1)==0)  version=3;
        else if( strcmp(magicxx,magic3_2)==0)  version=3;
        else if( strcmp(magicxx,magic4_0)==0)  version=4;
        else if( strcmp(magicxx,magic4_1)==0)  version=4;
        else err =1;

        if(err) {
            if(debug>2)
                cout << " Err read magic heading " << endl;
            goto Lreturn;
        }
        kread++;
        if(debug>2) cout << " Read entete " << version << endl;
        int c1 =getc(fp);
        if(c1==13)
            int c2 =getc(fp);

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
            cout << "\n\n\n next is build " << nextPlot<< " wait :" << nextPlot->wait << " -> " << kread <<  " gwait = " << gwait << "\n\n\n"<<endl;
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

    nbTimerNextPlot++;
    // the routine to  until the end of nextplot.
    // we use gluttimerfunc functionnaly
    //  remark, if we miss we retry.
    if(nbTimerNextPlot>1 && debug >2 ) cout << "       ######  Warning more than 1 nbTimerNextPlot" << nbTimerNextPlot << endl;
    value=min(1000,(value*3)/2);// try at leat every 1 second (not to heavy computation)
    if(TryNewPlot())
    {
        if(debug >9) cout << "    ####  TimerNextPlot TryNewPlot:: glutPostRedisplay "<< nbTimerNextPlot << " " << nbSendForNextPlot <<endl;
        nbSendForNextPlot--;// do the plot ..
        glutPostRedisplay();
    }
    else
    {
        if(debug >9) cout << "            ####  TimerNextPlot wait "<< nbTimerNextPlot << " " << nbSendForNextPlot<<endl;
        glutTimerFunc(value,TimerNextPlot,value);
    }
    nbTimerNextPlot--;
}

int SendForNextPlot()
{

    if(nbSendForNextPlot>0 && debug >2 ) cout << "       ######  Warning more than 1 SendForNextPlot" << nbSendForNextPlot+1 << endl;
    if(nbSendForNextPlot) return 0;
    nbSendForNextPlot++;

    //  to send a event to plot the date sheet.
    // and out a timer to wait to the end of read..
    // every 25/ second..  = 1000/25 = 40 ms
    if(NoMorePlot)
    {
        if(gwait )
        {usleep((useconds_t)(1e6*gwait)); Fin(0); }

        if((debug > 1)) cout << " send signal For Next plot, skip: No More Plot !  " <<  endl;
        nbSendForNextPlot--;
        return 0;
    }
    if((debug > 1)) cout << " Try to read read  plot "<< endl;
    //  put a timer for wait to the  end of read
    glutTimerFunc(40,TimerNextPlot,40);
    return 1;
}
static  bool TryNewPlot( void );





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
//  batch of primitives: the OpenGL of macOS (on Metal) is very slow with a lot of small
//  glBegin/glEnd (one by iso piece of each triangle with a fine mesh), so the caller opens
//  one glBegin(GL_TRIANGLES) (ffbatch=1: the polygons are drawn as triangle fans) or
//  one glBegin(GL_LINES) (ffbatch=2: the segments) for all the elements
static int ffbatch=0;
struct FFBatch {
    int old;
    bool on;
    FFBatch(int mode,bool ok=true) : old(ffbatch), on(ok)
    {
        if(on) { ffbatch=mode; glBegin(mode==2 ? GL_LINES : GL_TRIANGLES);}
    }
    ~FFBatch() { if(on) { glEnd(); ffbatch=old;} }
};

static void ffPolygon(const R3 *P,int n)
{
    if(ffbatch==1)
        for(int i=2;i<n;++i)
        {
            glVertex3f(P[0].x,P[0].y,P[0].z);
            glVertex3f(P[i-1].x,P[i-1].y,P[i-1].z);
            glVertex3f(P[i].x,P[i].y,P[i].z);
        }
    else
    {
        glBegin(GL_POLYGON);
        for(int i=0;i<n;++i) glVertex3f(P[i].x,P[i].y,P[i].z);
        glEnd();
    }
}

static void ffPolygonN(const R3 *P,const R3 *N,int n)
{ //  with a normal at each vertex
    if(ffbatch==1)
        for(int i=2;i<n;++i)
        {
            const int ii[3]={0,i-1,i};
            for(int j=0;j<3;++j)
            {
                glNormal3d(N[ii[j]].x,N[ii[j]].y,N[ii[j]].z);
                glVertex3f(P[ii[j]].x,P[ii[j]].y,P[ii[j]].z);
            }
        }
    else
    {
        glBegin(GL_POLYGON);
        for(int i=0;i<n;++i)
        {
            glNormal3d(N[i].x,N[i].y,N[i].z);
            glVertex3f(P[i].x,P[i].y,P[i].z);
        }
        glEnd();
    }
}

void drawisoTet(const R3 K[4],R f[4],R v)
{
    static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;//{ {2,1,3},{0,2,3},{1,0,3},{0,1,2} };

    R3 P[4];
    int nP=0;
    int np[4]={},nm[4]={};
    int km=0,kp=0;
    for (int i=0;i<4;++i)
    {
        if(f[i]<=v) nm[km++]=i;
        if(f[i]>=v) np[kp++]=i;
    }

    int h=-1,b[3]={};
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


    if(nP)
    {
        if(nP>2)
        {
            R3 N(R3(P[0],P[1])^R3(P[0],P[2]));
            N /= N.norme();
            glNormal3d(N.x,N.y,N.z);
        }
        ffPolygon(P,nP);
    }

    //  verification de l'orientation

}


//  pretty rendering: the isosurface in a tet with the normal at the vertices = interpolated gradient G
//  (smooth shading), same triangulation as drawisoTet
void drawisoTetN(const R3 K[4],R f[4],const R3 G[4],R v)
{
    static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;
    int np[4]={},nm[4]={};
    int km=0,kp=0;
    for (int i=0;i<4;++i)
    {
        if(f[i]<=v) nm[km++]=i;
        if(f[i]>=v) np[kp++]=i;
    }
    int ea[4],eb[4],nP=0;
    if(kp==1 && km==3)
    {
        int h=np[0];
        for(int j=0;j<3;++j) ea[j]=h,eb[j]=nvfaceTet[h][j];
        nP=3;
    }
    else if(km==1 && kp == 3)
    {
        int h=nm[0];
        ea[0]=ea[1]=ea[2]=h;
        eb[0]=nvfaceTet[h][0];
        eb[1]=nvfaceTet[h][2];
        eb[2]=nvfaceTet[h][1];
        nP=3;
    }
    else if(kp==2 && km==2)
    {
        if(signep4(nm[0],nm[1],np[0],np[1]) < 0)
            Exchange(nm[0],nm[1]);
        ea[0]=nm[0],eb[0]=np[0];
        ea[1]=nm[0],eb[1]=np[1];
        ea[2]=nm[1],eb[2]=np[1];
        ea[3]=nm[1],eb[3]=np[0];
        nP=4;
    }
    if(!nP) return;
    R3 P[4],N[4];
    for(int i=0;i<nP;++i)
    {
        R d=f[ea[i]]-f[eb[i]];
        R l= fabs(d)>1e-30 ? (f[ea[i]]-v)/d : 0.5;
        P[i]=K[ea[i]]*(1.-l)+K[eb[i]]*l;
        N[i]=G[ea[i]]*(1.-l)+G[eb[i]]*l;
    }
    R3 Ng(R3(P[0],P[1])^R3(P[0],P[2]));
    R ng=Ng.norme();
    if(ng>0) Ng/=ng;
    for(int i=0;i<nP;++i)
    {
        R3 n=N[i];
        if((n,Ng)<0) n=-n;// the normal on the side of the front face
        R nn=n.norme();
        N[i] = nn>0 ? n/nn : Ng;
    }
    ffPolygonN(P,N,nP);
}

//  the colormaps of the pretty rendering: t in [0,1] -> rgb
static void ColormapRGB(int kind,double t,float c[3])
{
    t = max(0.,min(1.,t));
    if(kind==0)
    { // turbo (polynomial approximation, A. Mikhailov 2019)
        c[0]=0.13572138+t*(4.61539260+t*(-42.66032258+t*(132.13108234+t*(-152.94239396+t*59.28637943))));
        c[1]=0.09140261+t*(2.19418839+t*(4.84296658+t*(-14.18503333+t*(4.27729857+t*2.82956604))));
        c[2]=0.10667330+t*(12.64194608+t*(-60.58204836+t*(110.36276771+t*(-89.90310912+t*27.34824973))));
    }
    else if(kind==1)
    { // viridis (polynomial approximation)
        static const double cc[7][3]={{0.2777273272234177,0.005407344544966578,0.3340998053353061},
            {0.1050930431085774,1.404613529898575,1.384590162594685},
            {-0.3308618287255563,0.214847559468213,0.09509516302823659},
            {-4.634230498983486,-5.799100973351585,-19.33244095627987},
            {6.228269936347081,14.17993336680509,56.69055260068105},
            {4.776384997670288,-13.74514537774601,-65.35303263337234},
            {-5.435455855934631,4.645852612178535,26.3124352495832}};
        for(int j=0;j<3;++j)
        {
            double r=cc[6][j];
            for(int k=5;k>=0;--k) r=cc[k][j]+t*r;
            c[j]=r;
        }
    }
    else if(kind==2)
    { // cool to warm (K. Moreland)
        static const double cw[5][3]={{59,76,192},{124,159,249},{221,221,221},{245,156,125},{180,4,38}};
        double x=t*4;
        int i=min(3,(int) x);
        double a=x-i;
        for(int j=0;j<3;++j) c[j]=((1-a)*cw[i][j]+a*cw[i+1][j])/255.;
    }
    else if(kind==3)
    { // the rainbow of ffglut
        float r,g,b;
        hsvToRgb(t*0.98,1.,1.,r,g,b);
        c[0]=r,c[1]=g,c[2]=b;
    }
    else
        c[0]=c[1]=c[2]=0.12+0.85*t;
    for(int j=0;j<3;++j) c[j]=max(0.f,min(1.f,c[j]));
}

inline R det3(const R3 & a,const R3 & b,const R3 & c) { return (a,b^c);}

//  the bound of the data of the plot (see SetDefIsoV)
static void DataRange(const ThePlot & plot,R & fmn,R & fmx)
{
    int n=plot.Viso.N();
    if(plot.pViso || n<3)
    {
        fmn=plot.Viso[0];
        fmx=plot.Viso[n-1];
    }
    else
    {
        R d=(plot.Viso[n-1]-plot.Viso[0])/(n-1);
        fmn=plot.Viso[0]+d/2;
        fmx=plot.Viso[n-1]-d/2;
    }
    if(!(fmx>fmn)) fmx=fmn+1;
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
    bool cc[3]= { plotborder , plotmesh && fill , plotmesh };
    int kk=0;
    if(cc[kk])
        if(lok[kk])   glCallList(gllists+kk);
        else
        {
            lok[kk]=1;
            glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
            glLineWidth(2*kscreenscale);
            glBegin(GL_LINES);
            for (int i=0;i<Th.nbBrdElmts();i++)
            {
                const typename  Mesh::BorderElement  & K(Th.be(i));
                plot.color(1+abs(K.lab));
                glVertex3d(K[0].x,K[0].y,z1);
                glVertex3d(K[1].x,K[1].y,z1);


            }
            glEnd();
            glLineWidth(kscreenscale);
            glEndList();  // fin de la list
        }
        else ;

    kk++;
    glEnable(GL_DEPTH_TEST);

    if(cc[kk])
    {
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
                glVertex3d(K[0].x,K[0].y,z2);
                glVertex3d(K[1].x,K[1].y,z2);
                glVertex3d(K[2].x,K[2].y,z2);

            }
            glEnd();
            glEndList();
        }
    }
    kk++;
    if(cc[kk])
    {
        glEnable(GL_DEPTH_TEST);
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
                glVertex3d(K[0].x,K[0].y,z1);
                glVertex3d(K[1].x,K[1].y,z1);
                glVertex3d(K[2].x,K[2].y,z1);

            }

            glEnd();
            glEndList();  // fin de la list
        }
    }
    ShowGlerror("end Mesh plot");

}


void Plot(const Mesh3 & Th,bool fill,bool plotmesh,bool plotborder,ThePlot & plot,GLint gllists,int * lok)
{
    typedef Mesh3::BorderElement BE;
    typedef Mesh3::Element Tet;
    glEnable(GL_DEPTH_TEST);
    ShowGlerror("begin Mesh3 plot");

    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

    R z1= plot.z0;
    R z2= plot.z0;

    double r=0,g=0,b=0;
    if((debug > 3)) cout<< " OnePlotMesh3::Draw " << plotmesh << " " << plotborder << " " <<  Th.nbBrdElmts() << " " << z1 << " "  << z2 << endl;
    bool cc[3]= { plotborder , plotborder && fill , plotmesh && fill };
    double saturation =plotmesh ? 1 :  0.25;
    int kk=0;
    if(Th.nbe==0)
    {
        if(lok[kk])   glCallList(gllists+kk);
        else
        { // draw all vertex
           glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
           glBegin(GL_POINTS);
            for(int i=0; i<Th.nv; ++i)
            {
             plot.color(1+Th(i).lab,saturation);
              glVertex3d(Th(i).x,Th(i).y,Th(i).z);
            }
            glEnd();
           glEndList();  // fin de la list
            
        }
    }
    else  if(cc[kk])
    {
        if(lok[kk])   glCallList(gllists+kk);
        else
        {
            lok[kk]=1;
            glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
            glLineWidth(kscreenscale);
            glAlphaFunc ( GL_GREATER, 0.1 ) ;
            glEnable(GL_ALPHA_TEST) ;
            if(!plotmesh)
            {
                glLineStipple(1, 0x300C);
                glEnable(GL_LINE_STIPPLE);
            }
            glBegin(GL_TRIANGLES);
            for (int i=0;i<Th.nbe;i++)
            {
                const BE & K(Th.be(i));
                plot.color(1+abs(K.lab),saturation);
                R3 N(R3(K[0],K[1])^R3(K[0],K[2]));
                N /= N.norme();
                glNormal3d(N.x,N.y,N.z);
                glVertex3d(K[0].x,K[0].y,K[0].z);
                glVertex3d(K[1].x,K[1].y,K[1].z);
                glVertex3d(K[2].x,K[2].y,K[2].z);
            }
            glEnd();
            glDisable(GL_LINE_STIPPLE);
            glLineWidth(kscreenscale);
            glDisable(GL_ALPHA_TEST) ;

            glEndList();  // fin de la list
        }
    }
    kk=2;
    if(cc[kk])
    {
        if(lok[kk])   glCallList(gllists+kk);
        else
        {
            lok[kk]=1;
            glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
            glLineWidth(kscreenscale);
            glPolygonMode(GL_FRONT,GL_FILL);//GL_FILL
            glBegin(GL_TRIANGLES);
            for (int i=0;i<Th.nbe;i++)
            {
                const BE & K(Th.be(i));
                plot.color(0,saturation);
                R3 N(R3(K[0],K[1])^R3(K[0],K[2]));
                N /= N.norme();
                glNormal3d(N.x,N.y,N.z);
                glVertex3d(K[0].x,K[0].y,K[0].z);
                glVertex3d(K[1].x,K[1].y,K[1].z);
                glVertex3d(K[2].x,K[2].y,K[2].z);
            }
            glEnd();
            glDisable(GL_LINE_STIPPLE);
            glLineWidth(kscreenscale);
            glDisable(GL_ALPHA_TEST) ;

            glEndList();  // fin de la list
        }
    }
    kk++;
    ShowGlerror("end Mesh3 plot");

}

void Plot(const MeshS & Th,bool fill,bool plotmesh,bool plotborder, ThePlot & plot,GLint gllists,int * lok, OneWindow *win)
{
    //glDisable(GL_DEPTH_TEST);

    ShowGlerror("begin MeshS plot");
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    R z1= plot.z0;
    R z2= plot.z0;
    R coef = plot.coeff;
    bool pNormalT=plot.pNormalT;
    int nbN=plot.nbN;
    
    double r=0,g=0,b=0;
    if((debug > 3)) cout<< " OnePlotMeshS::Draw " << plotmesh << " " << plotborder << " " <<  Th.nbBrdElmts() << " " << z1 << " "  << z2 << endl;
    // plot.SetColorTable(16) ;
    bool cc[3]= { plotborder , plotmesh && fill , plotmesh };
    int kk=0;

    if(cc[kk])
        if(lok[kk])   glCallList(gllists+kk);
        else
        {
            lok[kk]=1;
            glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
            glLineWidth(2*kscreenscale);
            glBegin(GL_LINES);
            for (int i=0;i<Th.nbBrdElmts();i++)
            {
                const MeshS::BorderElement  & K(Th.be(i));
                plot.color(1+abs(K.lab));
                glVertex3d(K[0].x,K[0].y,K[0].z);
                glVertex3d(K[1].x,K[1].y,K[1].z);


            }
            glEnd();
            glLineWidth(kscreenscale);
            glEndList();  // fin de la list
        }
        else ;

    kk++;
    if(cc[kk])
    {
        if(lok[kk])   glCallList(gllists+kk);
        else
        {
            lok[kk]=1;
            glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);//GL_FILL
            glBegin(GL_TRIANGLES);
            for (int i=0;i<Th.nt;i++)
            {
                const MeshS::Element & K(Th[i]);
                plot.color(K.lab?1+abs(K.lab):0);

                //glColor3d(r,g,b);
                glVertex3d(K[0].x,K[0].y,K[0].z);
                glVertex3d(K[1].x,K[1].y,K[1].z);
                glVertex3d(K[2].x,K[2].y,K[2].z);

            }
            glEnd();
            glEndList();
        }
    }
    kk++;
    if(cc[kk])
    {
        if(lok[kk])   glCallList(gllists+kk);
        else
        {
            lok[kk]=1;
            glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
            glPolygonMode(GL_FRONT,GL_LINE);
            glBegin(GL_TRIANGLES);
            for (int i=0;i<Th.nt;i++)
            {
                const MeshS::Element  & K(Th[i]);
                plot.color(fill? 1 : 1+abs(K.lab));
                glVertex3d(K[0].x,K[0].y,K[0].z);
                glVertex3d(K[1].x,K[1].y,K[1].z);
                glVertex3d(K[2].x,K[2].y,K[2].z);

            }

            glEnd();
            glEndList();  // fin de la list
        }
    }
    if (pNormalT) {
       if (debug>5) cout << " plot Normal element at triangles and the tangent to border " << pNormalT << endl;
            R h = 8*win->hpixel;
            coef = coef*sqrt(Th.mes/Th.nt);
            glLineWidth(0.5*kscreenscale);
            glColor3f(1,0,0);
            for (int i=0;i<Th.nbe;  i++ ) {
                typedef MeshS::BorderElement BE;
                typedef BE::Rd Rd;
                typedef BE::RdHat RdHat;
                const MeshS::BorderElement & K(Th.be(i));
                RdHat g(0.5);
                Rd G(K(g));
                Rd dd(K[0],K[1]);
                dd /= dd.norme();
                // fleche ou cone / orientation ???
                R3 nx(1.,0.,0.), ny(0.,1.,0.), nz(0.,0.,1.);
                // R3 dd = uv*(-h/l);
                R3 dnx = (dd^nx)*0.5, dny = (dd^ny)*0.5, dnz = (dd^nz)*0.5;
                dd *= -coef/dd.norme()/5;
                dnx *= -coef/dnx.norme()/5;
                dny *= -coef/dny.norme()/5;
                dnz *= -coef/dnz.norme()/5;

                glLineWidth(kscreenscale);
                glBegin(GL_LINES);
                win->Seg3(G,G+dd+dnx);
                win->Seg3(G,G+dd-dnx);
                win->Seg3(G,G+dd+dny);
                win->Seg3(G,G+dd-dny);
                win->Seg3(G,G+dd+dnz);
                win->Seg3(G,G+dd-dnz);
                glEnd();

                
            }

            for (int i=0;i<Th.nt; i = 0 ? i++ : i+=nbN) {
                const MeshS::Element & K(Th[i]);
                R3 NN=K.NormalTUnitaire();
                NN*=coef;
                R2 PtHat = R2::diag(1. / 3.);
                R3 A(K(PtHat));
                NN+=A;
                R3 dd(A,NN);
                double l = Max(sqrt((dd,dd)),1e-20);
                glBegin(GL_LINES);
                win->Seg3(A,NN);
                glEnd();
                
                // fleche ou cone / orientation ???
                R3 nx(1.,0.,0.), ny(0.,1.,0.), nz(0.,0.,1.);
                // R3 dd = uv*(-h/l);
                R3 dnx = (dd^nx)*0.5, dny = (dd^ny)*0.5, dnz = (dd^nz)*0.5;
                dd *= -coef/dd.norme()/5;
                dnx *= -coef/dnx.norme()/5;
                dny *= -coef/dny.norme()/5;
                dnz *= -coef/dnz.norme()/5;

                glLineWidth(kscreenscale);
                glBegin(GL_LINES);
                win->Seg3(NN,NN+dd+dnx);
                win->Seg3(NN,NN+dd-dnx);
                win->Seg3(NN,NN+dd+dny);
                win->Seg3(NN,NN+dd-dny);
                win->Seg3(NN,NN+dd+dnz);
                win->Seg3(NN,NN+dd-dnz);
                glEnd();
            }
            glLineWidth(kscreenscale);
        }
  

    ShowGlerror("end MeshS plot");

}




void Plot(const MeshL & Th,bool fill,bool plotmesh,bool plotborder,ThePlot & plot,GLint gllists,int * lok, OneWindow *win)
{
    glDisable(GL_DEPTH_TEST);
    
    ShowGlerror("begin MeshL plot");
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    R z1= plot.z0;
    R z2= plot.z0;
    R coef = plot.coeff;
    bool pNormalT=plot.pNormalT;
    int nbN=plot.nbN;

    double r=0,g=0,b=0;
    if((debug > 3)) cout<< " OnePlotMeshL::Draw " << plotmesh << " " << plotborder << " " <<  Th.nbBrdElmts() << " " << z1 << " "  << z2 << endl;
    // plot.SetColorTable(16) ;
    bool cc[3]= { plotborder , plotmesh && fill , plotmesh };
    int kk=0;
    
    if(cc[kk])
        if(lok[kk])
            glCallList(gllists+kk);
        else {
            lok[kk]=1;
            glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
            glPointSize(3);
            glBegin(GL_POINTS);
            for (int i=0;i<Th.nbBrdElmts();i++) {
                const MeshL::BorderElement  & K(Th.be(i));
                plot.color(1+abs(K.lab));
                glVertex3d(K[0].x,K[0].y,K[0].z);
            }
            glEnd();
            glPointSize(1);
            glEndList();  // fin de la list
        }
        else ;
    
    kk++;
    if(cc[kk]) {
        if(lok[kk])
            glCallList(gllists+kk);
        else {
            lok[kk]=1;
            glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
            glLineWidth(2*kscreenscale);
            //glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);//GL_FILL
            glBegin(GL_LINES);
            for (int i=0;i<Th.nt;i++) {
                const MeshL::Element & K(Th[i]);
                plot.color(K.lab?1+abs(K.lab):0);
                glVertex3d(K[0].x,K[0].y,K[0].z);
                glVertex3d(K[1].x,K[1].y,K[1].z);
                
            }
            glEnd();
            glLineWidth(kscreenscale);
            glEndList();
        }
    }
    kk++;
    if(cc[kk]) {
        if(lok[kk])
            glCallList(gllists+kk);
        else {
            lok[kk]=1;
            glNewList(gllists+kk,GL_COMPILE_AND_EXECUTE ); // save  la list sans affichage
            glPolygonMode(GL_FRONT,GL_LINE);
            glBegin(GL_LINES);
            for (int i=0;i<Th.nt;i++) {
                const MeshL::Element  & K(Th[i]);
                plot.color(fill? 1 : 1+abs(K.lab));
                glVertex3d(K[0].x,K[0].y,K[0].z);
                glVertex3d(K[1].x,K[1].y,K[1].z);
                
            }
            glEnd();
            glEndList();  // fin de la list
        }
    }
    if (pNormalT) {
        coef = coef*(Th.mes/Th.nt);
        glLineWidth(0.5*kscreenscale);
        glColor3f(1,0,0);

        for (int i=0;i<Th.nt; i = 0 ? i++ : i+=nbN) {
            const MeshL::Element & K(Th[i]);
            R3 NN=K.NormalTUnitaire();
            NN*=coef;
            R1 PtHat = R1::diag(1. / 2.);
            R3 A(K(PtHat));
            NN+=A;
            R3 dd(A,NN);
            double l = Max(sqrt((dd,dd)),1e-20);
            glBegin(GL_LINES);
            win->Seg3(A,NN);
            glEnd();

            // fleche ou cone / orientation ???
            R3 nx(1.,0.,0.), ny(0.,1.,0.), nz(0.,0.,1.);
            // R3 dd = uv*(-h/l);
            R3 dnx = (dd^nx)*0.5, dny = (dd^ny)*0.5, dnz = (dd^nz)*0.5;
            dd *= -coef/dd.norme()/5;
            dnx *= -coef/dnx.norme()/5;
            dny *= -coef/dny.norme()/5;
            dnz *= -coef/dnz.norme()/5;

            glLineWidth(kscreenscale);
            glBegin(GL_LINES);
            win->Seg3(NN,NN+dd+dnx);
            win->Seg3(NN,NN+dd-dnx);
            win->Seg3(NN,NN+dd+dny);
            win->Seg3(NN,NN+dd-dny);
            win->Seg3(NN,NN+dd+dnz);
            win->Seg3(NN,NN+dd-dnz);
            glEnd();
        }
        glLineWidth(kscreenscale);
    }
    ShowGlerror("end MeshL plot");
    
}



void OnePlotError::Draw(OneWindow *win)
{
    initlist();
    win->SetScreenView() ;
    glColor3d(0.,0.,0.);
    cout << " Error plot item empty " << item <<  endl;
    char s[100];
    snprintf(s,100,"Warning the item %ld fot the plot is empty",item);
    win->Show(s,4+item*2);
    win->SetView() ;
}

static int CutTet(const R3 K[4],const R d[4],const R f[4],R3 P[4],R g[4]);

template<class Mesh>
void OnePlotMesh<Mesh>::Draw(OneWindow *win)
{
    initlist();
    ThePlot & plot=*win->theplot;
    //  if a scalar FE function on this mesh is in the plot, then the mesh is plotted
    //  on the surface z = f(x,y) by OnePlotFE (correct with a 3d view)
    for (list<OnePlot *>::iterator i= plot.plots.begin();i != plot.plots.end(); ++i)
    {
        OnePlotFE<Mesh> * pf=dynamic_cast<OnePlotFE<Mesh> *>(*i);
        if(pf && pf->Th == Th && pf->what%10==1) return;
    }
    Plot(*Th,plot.fill,true,true,plot,gllists,oklist);
    ShowGlerror("OnePlotMesh::Draw");
}

void OnePlotMesh3::Draw(OneWindow *win)
{
    initlist();
    ThePlot & plot=*win->theplot;
    //  if a scalar FE function on this mesh is in the plot, then the mesh is plotted by OnePlotFE3
    for (list<OnePlot *>::iterator i= plot.plots.begin();i != plot.plots.end(); ++i)
    {
        OnePlotFE3 * pf=dynamic_cast<OnePlotFE3 *>(*i);
        if(pf && pf->Th == Th && pf->what%10==6) return;
    }
    //  cut plane: the border is clipped, and the section of the elements is plotted
    R3 N;
    R c=0;
    bool cut=win->SetCutPlane(Th->Pmin,Th->Pmax,N,c);
    if(cut) glEnable(GL_CLIP_PLANE0);
    Plot(*Th,plot.fill,true,true,plot,gllists,oklist);
    glDisable(GL_CLIP_PLANE0);
    if(cut)
    {
        if(cutok && c==lastcutc && N.x==lastcutN.x && N.y==lastcutN.y && N.z==lastcutN.z)
            glCallList(gllists+1);
        else
        {
            cutok=true;
            lastcutc=c;
            lastcutN=N;
            glNewList(gllists+1,GL_COMPILE_AND_EXECUTE);
            DrawSection(plot,N,c);
            glEndList();
        }
    }
    ShowGlerror("OnePlotMesh3::Draw");
}

void OnePlotMesh3::DrawSection(ThePlot & plot,R3 N,R c)
{ //  the section of the elements by the plane (N,P) = c, the color is given by the label of the element
    R3 P[4];
    R g[4];
    glEnable(GL_DEPTH_TEST);
    glLineWidth(kscreenscale);
    for(int pass=0;pass<2;++pass)
    { // pass 0: the filled polygons (pushed back), pass 1: the edges
        if(pass==0)
        {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.,1.);
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        }
        for(int k=0;k<Th->nt;++k)
        {
            const Mesh3::Element & K=(*Th)[k];
            R3 Pt[4];
            R d[4],ff[4]={0.,0.,0.,0.};
            for(int i=0;i<4;++i)
            {
                Pt[i]=K[i];
                d[i]=(N,Pt[i])-c;
            }
            int n=CutTet(Pt,d,ff,P,g);
            if(n==0) continue;
            plot.color(pass ? 1 : (K.lab ? 1+abs(K.lab) : 0));
            glBegin(pass ? GL_LINE_LOOP : GL_POLYGON);
            for(int i=0;i<n;++i)
                glVertex3d(P[i].x,P[i].y,P[i].z);
            glEnd();
        }
        glDisable(GL_POLYGON_OFFSET_FILL);
    }
}

void OnePlotMeshS::Draw(OneWindow *win)
{
    initlist();
    ThePlot & plot=*win->theplot;
    Plot(*Th,plot.fill,true,true,plot,gllists,oklist,win);
    ShowGlerror("OnePlotMeshS::Draw");
}

void OnePlotMeshL::Draw(OneWindow *win)
{
    initlist();
    ThePlot & plot=*win->theplot;
    Plot(*Th,plot.fill,true,true,plot,gllists,oklist,win);
    ShowGlerror("OnePlotMeshL::Draw");
}

void OnePlotFE3::BuildBoundaryFaces()
{ // the boundary faces are the faces which are in only one tet
    if(bfacesok) return;
    bfacesok=true;
    if(Th->TheAdjacencesLink)
    { //  fast: the adjacency of the mesh (-1: no tet on the other side)
        for(int k=0;k<Th->nt;++k)
            for(int e=0;e<4;++e)
                if(Th->TheAdjacencesLink[4*k+e]<0) bfaces.push_back(4*k+e);
    }
    else
    {
    typedef pair<pair<int,int>,pair<int,int> > Face;// the 3 sorted vertices and 4*k+e
    vector<Face> lf;
    lf.reserve(4*Th->nt);
    for(int k=0;k<Th->nt;++k)
        for(int e=0;e<4;++e)
        {
            int i[3],j=0;
            for(int l=0;l<4;++l)
                if(l!=e) i[j++]=(*Th)(k,l);
            sort(i,i+3);
            lf.push_back(Face(make_pair(i[0],i[1]),make_pair(i[2],4*k+e)));
        }
    sort(lf.begin(),lf.end());
    int n=lf.size();
    for(int l=0;l<n;++l)
    {
        bool same0 = l>0 && lf[l-1].first==lf[l].first && lf[l-1].second.first==lf[l].second.first;
        bool same1 = l+1<n && lf[l+1].first==lf[l].first && lf[l+1].second.first==lf[l].second.first;
        if(!same0 && !same1)
            bfaces.push_back(lf[l].second.second);
    }
    }
    //  the sharp edges: angle between the normal of the 2 faces  > 30 degree (or edge in only one face)
    typedef map<pair<int,int>,pair<R3,int> > MapEdge;
    MapEdge me;
    for(size_t l=0;l<bfaces.size();++l)
    {
        int k=bfaces[l]/4, e=bfaces[l]%4;
        const Mesh3::Element & K=(*Th)[k];
        int i[3],j=0;
        for(int ll=0;ll<4;++ll)
            if(ll!=e) i[j++]=ll;
        R3 N(R3(K[i[0]],K[i[1]])^R3(K[i[0]],K[i[2]]));
        N /= N.norme();
        //  outward normal
        if( (N,R3(K[i[0]],K[e])) > 0) N = -N;
        for(j=0;j<3;++j)
        {
            int a=(*Th)(k,i[j]), b=(*Th)(k,i[(j+1)%3]);
            pair<int,int> ab(min(a,b),max(a,b));
            MapEdge::iterator pe=me.find(ab);
            if(pe==me.end())
                me[ab]=make_pair(N,1);
            else
            {
                pe->second.second++;
                if( (pe->second.first,N) < 0.866 )
                {
                    bedges.push_back(a);
                    bedges.push_back(b);
                }
            }
        }
    }
    for(MapEdge::iterator pe=me.begin();pe!=me.end();++pe)
        if(pe->second.second==1)
        {
            bedges.push_back(pe->first.first);
            bedges.push_back(pe->first.second);
        }
}

void OnePlotFE3::BoundaryConstant()
{ // is the function (nearly) constant on the boundary ? (ex. Dirichlet boundary condition)
    BuildBoundaryFaces();
    int nsubV=Psub.N();
    int nK=v.N()/ Th->nt;
    R bmin=dinfty, bmax=-dinfty;
    for(size_t l=0;l<bfaces.size();++l)
    {
        int k=bfaces[l]/4, e=bfaces[l]%4, o=k*nK;
        for(int i=0;i<nsubV;++i)
        {
            R3 Q=Psub[i];
            R lambda[4]={1.-Q.x-Q.y-Q.z,Q.x,Q.y,Q.z};
            if(fabs(lambda[e]) < 1e-6)
            {
                bmin=min(bmin,v[o+i]);
                bmax=max(bmax,v[o+i]);
            }
        }
    }
    bconstant = (bmax-bmin) <= 0.05*(v.max()-v.min());
    bval = (bmin+bmax)/2;
    //  is the function constant on each element ?
    p0=true;
    R eps=1e-10*(v.max()-v.min());
    for(int k=0,o=0;k<Th->nt && p0;++k, o+= nK)
        for(int i=1;i<nK;++i)
            if(fabs(v[o+i]-v[o]) > eps) p0=false;
    if(p0) bconstant=false;// by default all the elements
}

static int ColorOfValue(ThePlot & plot,R val)
{ // the color of the value in the filled iso (see DrawIsoTfill)
    int N=plot.Viso.N();
    if(val<plot.Viso[0]) return 4;
    if(val>plot.Viso[N-1]) return N+4;
    return 5+dichotomie(plot.Viso,val);
}

void OnePlotFE3::BuildExplode()
{ //  exploded view: the direction of the move of each element = center of its group - center of the mesh,
  //  a group = the elements with the same value, or each element if there are too many values
    int nK=v.N()/ Th->nt;
    typedef map<R,pair<R3,R> > MapGroup;
    MapGroup g;
    vector<R3> G(Th->nt);
    R3 C(0.,0.,0.);
    R vol=0;
    for(int k=0;k<Th->nt;++k)
    {
        const Mesh3::Element & K=(*Th)[k];
        G[k]=((R3) K[0]+(R3) K[1]+(R3) K[2]+(R3) K[3])/4.;
        R m=K.mesure();
        pair<R3,R> & gk=g[v[k*nK]];
        if(gk.second==0) gk.first=R3(0.,0.,0.);
        gk.first += G[k]*m;
        gk.second += m;
        C += G[k]*m;
        vol += m;
    }
    C /= vol;
    bool bygroup = g.size() <= 64;
    expdir.resize(Th->nt);
    for(int k=0;k<Th->nt;++k)
    {
        const pair<R3,R> & gk=g[v[k*nK]];
        expdir[k] = (bygroup ? gk.first/gk.second : G[k]) - C;
    }
}

void OnePlotFE3::DrawThreshold(OneWindow *win,R thr,bool plotmesh,R explode)
{ //  P0 function: the elements with a value >= thr, one color by element
  //  explode > 0 : exploded view
    static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;
    ThePlot & plot=*win->theplot;
    int nK=v.N()/ Th->nt;
    if(explode>0) BuildExplode();
    //  the filled faces are pushed back, so the edges are not hidden by the faces
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.,1.);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glBegin(GL_TRIANGLES);
    for(int k=0;k<Th->nt;++k)
    {
        R val=v[k*nK];
        if(val<thr) continue;
        const Mesh3::Element & K=(*Th)[k];
        R3 D = explode>0 ? expdir[k]*explode : R3(0.,0.,0.);
        plot.color(ColorOfValue(plot,val));
        for(int e=0;e<4;++e)
        {
            R3 A=D+K[nvfaceTet[e][0]], B=D+K[nvfaceTet[e][1]], C=D+K[nvfaceTet[e][2]];
            R3 N(R3(A,B)^R3(A,C));
            N /= N.norme();
            glNormal3d(N.x,N.y,N.z);
            glVertex3d(A.x,A.y,A.z);
            glVertex3d(B.x,B.y,B.z);
            glVertex3d(C.x,C.y,C.z);
        }
    }
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);
    if(plotmesh)
    { // the edges of the elements
        glDisable(GL_LIGHTING);
        glLineWidth(kscreenscale);
        plot.color(1);
        glBegin(GL_LINES);
        for(int k=0;k<Th->nt;++k)
        {
            if(v[k*nK]<thr) continue;
            const Mesh3::Element & K=(*Th)[k];
            R3 D = explode>0 ? expdir[k]*explode : R3(0.,0.,0.);
            for(int i=0;i<4;++i)
                for(int j=i+1;j<4;++j)
                {
                    R3 P=D+K[i], Q=D+K[j];
                    glVertex3d(P.x,P.y,P.z);
                    glVertex3d(Q.x,Q.y,Q.z);
                }
        }
        glEnd();
    }
}

void OnePlotFE3::DrawBoundaryEdges(ThePlot & plot)
{ // the sharp edges of the boundary to see the shape of the domain
    BuildBoundaryFaces();
    glDisable(GL_LIGHTING);
    glLineWidth(kscreenscale);
    plot.color(1);
    glBegin(GL_LINES);
    for(size_t l=0;l<bedges.size();++l)
    {
        const R3 & P=(*Th)(bedges[l]);
        glVertex3d(P.x,P.y,P.z);
    }
    glEnd();
}

void OnePlotFE3::ComputeGradV()
{ //  the gradient of the P1 interpolation of the function in each tet, averaged at the vertices
    int nK=v.N()/ Th->nt, nsubV=Psub.N();
    int iv4[4]={-1,-1,-1,-1};// the sub vertices at the vertices of the tet
    for(int i=0;i<nsubV;++i)
    {
        R3 Q=Psub[i];
        R lambda[4]={1.-Q.x-Q.y-Q.z,Q.x,Q.y,Q.z};
        for(int l=0;l<4;++l)
            if(fabs(lambda[l]-1.)<1e-6) iv4[l]=i;
    }
    gradv.assign(Th->nv,R3(0.,0.,0.));
    if(iv4[0]<0 || iv4[1]<0 || iv4[2]<0 || iv4[3]<0) return;
    vector<R> w(Th->nv,0.);
    for(int k=0,o=0;k<Th->nt;++k,o+=nK)
    {
        const Mesh3::Element & K=(*Th)[k];
        R3 GL[4];
        K.Gradlambda(GL);
        R3 g(0.,0.,0.);
        for(int l=0;l<4;++l) g += GL[l]*v[o+iv4[l]];
        R m=K.mesure();
        for(int l=0;l<4;++l)
        {
            int iv=(*Th)(k,l);
            gradv[iv] += g*m;
            w[iv] += m;
        }
    }
    for(int i=0;i<Th->nv;++i)
        if(w[i]>0) gradv[i]/=w[i];
}

void OnePlotFE3::DrawBoundaryPretty(OneWindow *win,bool plotmesh)
{ //  pretty rendering: the function on the boundary with continuous colors (texture of the colormap),
  //  smooth normals (averaged on the faces with an angle < 40 degree, so the sharp edges stay sharp)
    static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;
    ThePlot & plot=*win->theplot;
    BuildBoundaryFaces();
    int nsubV=Psub.N();
    int nsubT=Ksub.N()/4;
    int nK=v.N()/ Th->nt;
    R fmn,fmx;
    DataRange(plot,fmn,fmx);
    //  the sub triangles on the face e of the reference tet
    vector<int> subf[4];
    for(int e=0;e<4;++e)
        for(int sk=0;sk<nsubT;++sk)
            for(int f=0;f<4;++f)
            {
                int i[3],j=0;
                bool on=true;
                for(int l=0;l<4;++l)
                    if(l!=f)
                    {
                        int iv=Ksub[4*sk+l];
                        R3 Q=Psub[iv];
                        R lambda[4]={1.-Q.x-Q.y-Q.z,Q.x,Q.y,Q.z};
                        on = on && fabs(lambda[e]) < 1e-6;
                        i[j++]=iv;
                    }
                if(on)
                    for(int l=0;l<3;++l)
                        subf[e].push_back(i[l]);
            }
    //  the outward normal of the faces, and the faces at each vertex
    int nf=bfaces.size();
    vector<R3> fn(nf);
    map<int,vector<int> > vf;
    for(int l=0;l<nf;++l)
    {
        int k=bfaces[l]/4, e=bfaces[l]%4;
        const Mesh3::Element & K=(*Th)[k];
        R3 A=K[nvfaceTet[e][0]], B=K[nvfaceTet[e][1]], C=K[nvfaceTet[e][2]];
        R3 N(R3(A,B)^R3(A,C));
        if((N,R3(A,(R3) K[e]))>0) N=-N;
        fn[l]=N;// the norm = 2 * area: weight of the average
        for(int j=0;j<3;++j) vf[(*Th)(k,nvfaceTet[e][j])].push_back(l);
    }
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.,1.);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D,win->cmaptex);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
    glColor4f(1.,1.,1.,1.);
    glBegin(GL_TRIANGLES);
    for(int l=0;l<nf;++l)
    {
        int k=bfaces[l]/4, e=bfaces[l]%4, o=k*nK;
        const Mesh3::Element & K=(*Th)[k];
        R3 nl=fn[l]/fn[l].norme();
        R3 nv[3];
        for(int j=0;j<3;++j)
        {
            const vector<int> & fl=vf[(*Th)(k,nvfaceTet[e][j])];
            R3 n(0.,0.,0.);
            for(size_t m=0;m<fl.size();++m)
                if((fn[fl[m]],nl) > 0.766*fn[fl[m]].norme()) n += fn[fl[m]];
            nv[j]=n/n.norme();
        }
        for(size_t j=0;j<subf[e].size();j+=3)
        {
            int ii[3]={subf[e][j],subf[e][j+1],subf[e][j+2]};
            R3 P[3],N[3];
            for(int m=0;m<3;++m)
            {
                R3 Q=Psub[ii[m]];
                R lambda[4]={1.-Q.x-Q.y-Q.z,Q.x,Q.y,Q.z};
                P[m]=K(Q);
                N[m]=nv[0]*lambda[nvfaceTet[e][0]]+nv[1]*lambda[nvfaceTet[e][1]]+nv[2]*lambda[nvfaceTet[e][2]];
                N[m]/=N[m].norme();
            }
            //  the front face is the outward face
            if(((R3(P[0],P[1])^R3(P[0],P[2])),nl)<0)
                Exchange(ii[1],ii[2]),Exchange(P[1],P[2]),Exchange(N[1],N[2]);
            for(int m=0;m<3;++m)
            {
                glNormal3d(N[m].x,N[m].y,N[m].z);
                glTexCoord1d((v[o+ii[m]]-fmn)/(fmx-fmn));
                glVertex3d(P[m].x,P[m].y,P[m].z);
            }
        }
    }
    glEnd();
    glDisable(GL_TEXTURE_1D);
    glDisable(GL_POLYGON_OFFSET_FILL);
    if(!plotmesh)
        DrawBoundaryEdges(plot);
    else
        DrawBoundaryMesh();
}

static bool SubInterp2(const KN<R2> & Psub,const KN<int> & Ksub,const KN<double> & v,int o,int nc,R2 Q,R w3[3]);

bool OnePlotFE3::Probe(const R3 & P,R & val,R3 & vec,int & ncomp)
{ //  3d: the tet which contains P (P is on the visible surface: a small tolerance)
    if(what%10!=6 && what%10!=7) return false;
    if(!grid.built)
    {
        vector<R3> mn(Th->nt),mx(Th->nt);
        for(int k=0;k<Th->nt;++k)
        {
            const Mesh3::Element & K=(*Th)[k];
            mn[k]=mx[k]=K[0];
            for(int l=1;l<4;++l) { mn[k]=Minc(mn[k],(R3) K[l]); mx[k]=Maxc(mx[k],(R3) K[l]);}
        }
        grid.build(mn,mx);
    }
    if(grid.head.empty()) return false;
    int c[3];
    grid.ijk(P,c);
    int nK=v.N()/ Th->nt, nc= what%10==6 ? 1 : 3;
    int best=-1;
    R bmin=-0.05;// tolerance in the barycentric coordinates
    R3 Qbest;
    for(int dk=-1;dk<=1;++dk) for(int dj=-1;dj<=1;++dj) for(int di=-1;di<=1;++di)
    {
        int i=c[0]+di, j=c[1]+dj, kk=c[2]+dk;
        if(i<0 || j<0 || kk<0 || i>=grid.n[0] || j>=grid.n[1] || kk>=grid.n[2]) continue;
        for(int q=grid.head[grid.cell(i,j,kk)];q>=0;q=grid.nxt[q])
        {
            int k=grid.elem[q];
            const Mesh3::Element & K=(*Th)[k];
            R3 GL[4];
            K.Gradlambda(GL);
            R3 AP((R3) K[0],P);
            R l1=(GL[1],AP), l2=(GL[2],AP), l3=(GL[3],AP), l0=1-l1-l2-l3;
            R mn=min(min(l0,l1),min(l2,l3));
            if(mn>bmin) { bmin=mn; best=k; Qbest=R3(l1,l2,l3);}
        }
    }
    if(best<0) return false;
    //  the sub tet
    int o=best*nK, nsubT=Ksub.N()/4, bs=-1;
    R smin=-1e100, mu[4]={0,0,0,0};
    for(int sk=0;sk<nsubT;++sk)
    {
        R3 A=Psub[Ksub[4*sk]], B=Psub[Ksub[4*sk+1]], C=Psub[Ksub[4*sk+2]], D=Psub[Ksub[4*sk+3]];
        R3 AB(A,B),AC(A,C),AD(A,D),AQ(A,Qbest);
        R det=det3(AB,AC,AD);
        if(det==0) continue;
        R m1=det3(AQ,AC,AD)/det, m2=det3(AB,AQ,AD)/det, m3=det3(AB,AC,AQ)/det, m0=1-m1-m2-m3;
        R mn=min(min(m0,m1),min(m2,m3));
        if(mn>smin) { smin=mn; bs=sk; mu[0]=m0; mu[1]=m1; mu[2]=m2; mu[3]=m3;}
    }
    if(bs<0) return false;
    R sw=0;
    for(int j=0;j<4;++j) sw += (mu[j]=max(0.,mu[j]));
    for(int j=0;j<4;++j) mu[j]/=sw;
    R w[3]={0,0,0};
    for(int cc=0;cc<nc;++cc)
        for(int j=0;j<4;++j) w[cc] += mu[j]*v[o+nc*Ksub[4*bs+j]+cc];
    ncomp=nc;
    if(nc==1) val=w[0];
    else { vec=R3(w[0],w[1],w[2]); val=vec.norme();}
    return true;
}

bool OnePlotFE3::ProbeRay(const R3 & A,const R3 & B,OneWindow *win,R3 & P,R & val)
{ //  the first point of the segment [A,B] (from the eye) on an isosurface (the translucent isosurfaces
  //  are not in the depth buffer): the function is sampled along the segment in the mesh,
  //  and the crossing of an iso value is refined by dichotomy
    if(what%10!=6 || isovals.empty()) return false;
    R dummy;
    R3 vec;
    int nc;
    if(!grid.built) Probe(A,dummy,vec,nc);// to build the grid
    //  the segment in the bounding box of the mesh
    R t0=0,t1=1;
    R a[3]={A.x,A.y,A.z}, d[3]={B.x-A.x,B.y-A.y,B.z-A.z};
    R mn[3]={Th->Pmin.x,Th->Pmin.y,Th->Pmin.z}, mx[3]={Th->Pmax.x,Th->Pmax.y,Th->Pmax.z};
    for(int k=0;k<3;++k)
    {
        if(fabs(d[k])<1e-300) { if(a[k]<mn[k] || a[k]>mx[k]) return false; continue;}
        R ta=(mn[k]-a[k])/d[k], tb=(mx[k]-a[k])/d[k];
        if(ta>tb) swap(ta,tb);
        t0=max(t0,ta);
        t1=min(t1,tb);
    }
    if(!(t1>t0)) return false;
    R3 D=B-A;
    R hmin=min(grid.h.x,min(grid.h.y,grid.h.z));
    int ns=min(4000,max(10,(int) ceil((t1-t0)*D.norme()/(0.25*hmin))));
    bool cut = win->cutaxis;
    #define FFINSIDE(Q) (!cut || win->cuteq[0]*(Q).x+win->cuteq[1]*(Q).y+win->cuteq[2]*(Q).z+win->cuteq[3] >= 0)
    R tp=0,fp=0;
    bool okp=false;
    for(int i=0;i<=ns;++i)
    {
        R t=t0+(t1-t0)*i/ns, f;
        R3 Q=A+D*t;
        bool ok = FFINSIDE(Q) && Probe(Q,f,vec,nc);
        if(ok && okp)
            for(size_t l=0;l<isovals.size();++l)
            {
                R v=isovals[l];
                if((fp-v)*(f-v)<=0 && fp!=f)
                { // dichotomy
                    R ta=tp, tb=t, fa=fp;
                    for(int it=0;it<30;++it)
                    {
                        R tm=(ta+tb)/2, fm;
                        if(!Probe(A+D*tm,fm,vec,nc)) break;
                        if((fa-v)*(fm-v)<=0) tb=tm; else ta=tm, fa=fm;
                    }
                    P=A+D*((ta+tb)/2);
                    val=v;
                    #undef FFINSIDE
                    return true;
                }
            }
        okp=ok;
        tp=t;
        fp=f;
    }
    return false;
}

bool OnePlotFES::Probe(const R3 & P,R & val,R3 & vec,int & ncomp)
{ //  3d surface: the triangle nearest to P
    const MeshS & Th(*this->Th);
    if(what%10!=8 && what%10!=9) return false;
    if(!grid.built)
    {
        vector<R3> mn(Th.nt),mx(Th.nt);
        for(int k=0;k<Th.nt;++k)
        {
            const MeshS::Element & K=Th[k];
            mn[k]=mx[k]=K[0];
            for(int l=1;l<3;++l) { mn[k]=Minc(mn[k],(R3) K[l]); mx[k]=Maxc(mx[k],(R3) K[l]);}
            R3 d=(mx[k]-mn[k])*0.1;// the point under the mouse is near the surface
            mn[k]-=d;
            mx[k]+=d;
        }
        grid.build(mn,mx);
    }
    if(grid.head.empty()) return false;
    int c[3];
    grid.ijk(P,c);
    int nK=v.N()/ Th.nt, nc= what%10==8 ? 1 : 3;
    int best=-1;
    R dbest=1e100;
    R2 Qbest;
    for(int q=grid.head[grid.cell(c[0],c[1],c[2])];q>=0;q=grid.nxt[q])
    {
        int k=grid.elem[q];
        const MeshS::Element & K=Th[k];
        R3 A=K[0],B=K[1],C=K[2];
        R3 AB(A,B),AC(A,C),AP(A,P);
        R a11=(AB,AB), a12=(AB,AC), a22=(AC,AC), b1=(AB,AP), b2=(AC,AP);
        R det=a11*a22-a12*a12;
        if(det<=0) continue;
        R l1=(b1*a22-b2*a12)/det, l2=(a11*b2-a12*b1)/det;
        if(l1<-0.05 || l2<-0.05 || l1+l2>1.05) continue;
        R d=(AP-AB*l1-AC*l2).norme()/sqrt(max(a11,a22));
        if(d<dbest) { dbest=d; best=k; Qbest=R2(l1,l2);}
    }
    if(best<0 || dbest>0.1) return false;
    R w[3];
    if(!SubInterp2(Psub,Ksub,v,best*nK,nc,Qbest,w)) return false;
    ncomp=nc;
    if(nc==1) val=w[0];
    else { vec=R3(w[0],w[1],w[2]); val=vec.norme();}
    return true;
}

void OnePlotFE3::DrawBoundaryMesh()
{ //  pretty rendering: the edges of the boundary faces in dark grey
    static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;
    BuildBoundaryFaces();
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_1D);
    glLineWidth(kscreenscale);
    glColor3f(0.15f,0.15f,0.2f);
    glBegin(GL_LINES);
    for(size_t l=0;l<bfaces.size();++l)
    {
        const Mesh3::Element & K=(*Th)[bfaces[l]/4];
        int e=bfaces[l]%4;
        for(int j=0;j<3;++j)
        {
            const R3 & P=K[nvfaceTet[e][j]], & Q=K[nvfaceTet[e][(j+1)%3]];
            glVertex3d(P.x,P.y,P.z);
            glVertex3d(Q.x,Q.y,Q.z);
        }
    }
    glEnd();
}

void OnePlotFE3::DrawBoundary(OneWindow *win,bool plotmesh)
{ //  the function on the boundary of the mesh: iso lines (with hidden faces) or filled iso
    static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;
    ThePlot & plot=*win->theplot;
    if(win->Pretty3d() && !p0)
    {
        DrawBoundaryPretty(win,plotmesh);
        return;
    }
    BuildBoundaryFaces();
    int nsubV=Psub.N();
    int nsubT=Ksub.N()/4;
    int nK=v.N()/ Th->nt;
    KN<R3> Pn(nsubV);
    //  the sub triangles on the face e of the reference tet
    vector<int> subf[4];
    for(int e=0;e<4;++e)
        for(int sk=0;sk<nsubT;++sk)
            for(int f=0;f<4;++f)
            {
                int i[3],j=0;
                bool on=true;
                for(int l=0;l<4;++l)
                    if(l!=f)
                    {
                        int iv=Ksub[4*sk+l];
                        R3 Q=Psub[iv];
                        R lambda[4]={1.-Q.x-Q.y-Q.z,Q.x,Q.y,Q.z};
                        on = on && fabs(lambda[e]) < 1e-6;
                        i[j++]=iv;
                    }
                if(on)
                    for(int l=0;l<3;++l)
                        subf[e].push_back(i[l]);
            }
    //  the filled faces are pushed back, so the lines on the faces are not hidden by the faces
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.,1.);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    if(p0)
    { //  one color by face
        glBegin(GL_TRIANGLES);
        for(size_t l=0;l<bfaces.size();++l)
        {
            int k=bfaces[l]/4, e=bfaces[l]%4;
            const Mesh3::Element & K=(*Th)[k];
            R3 A=K[nvfaceTet[e][0]], B=K[nvfaceTet[e][1]], C=K[nvfaceTet[e][2]];
            R3 N(R3(A,B)^R3(A,C));
            N /= N.norme();
            plot.color(ColorOfValue(plot,v[k*nK]));
            glNormal3d(N.x,N.y,N.z);
            glVertex3d(A.x,A.y,A.z);
            glVertex3d(B.x,B.y,B.z);
            glVertex3d(C.x,C.y,C.z);
        }
        glEnd();
    }
    else if(!plot.fill)
    { //  white faces to hide the iso lines of the back faces
        glDisable(GL_LIGHTING);
        plot.color(0);
        glBegin(GL_TRIANGLES);
        for(size_t l=0;l<bfaces.size();++l)
        {
            const Mesh3::Element & K=(*Th)[bfaces[l]/4];
            int e=bfaces[l]%4;
            for(int j=0;j<3;++j)
            {
                const R3 & P=K[nvfaceTet[e][j]];
                glVertex3d(P.x,P.y,P.z);
            }
        }
        glEnd();
    }
    {
    FFBatch batch(plot.fill ? 1 : 2,!p0);// one glBegin for all the faces
    for(size_t l=0;l<bfaces.size() && !p0;++l)
    {
        int k=bfaces[l]/4, e=bfaces[l]%4, o=k*nK;
        const Mesh3::Element & K=(*Th)[k];
        R3 A=K[nvfaceTet[e][0]], B=K[nvfaceTet[e][1]], C=K[nvfaceTet[e][2]];
        R3 N(R3(A,B)^R3(A,C));
        N /= N.norme();
        R3 Nt[3]={N,N,N};
        for(size_t j=0;j<subf[e].size();j+=3)
        {
            int i0=subf[e][j], i1=subf[e][j+1], i2=subf[e][j+2];
            R ff[3]={v[o+i0],v[o+i1],v[o+i2]};
            R3 Pt[3]={K(Psub[i0]),K(Psub[i1]),K(Psub[i2])};
            if(plot.fill)
                plot.DrawIsoTfill( Pt, ff, Nt, plot.Viso,plot.Viso.N(),false,3);
            else
                plot.DrawIsoT( Pt, ff, Nt, plot.Viso,plot.Viso.N(),false,3);
        }
    }
    }
    glDisable(GL_POLYGON_OFFSET_FILL);
    if(!plotmesh)
        DrawBoundaryEdges(plot);
    else
    { // the edges of the boundary faces
        glDisable(GL_LIGHTING);
        glLineWidth(kscreenscale);
        plot.color(1);
        glBegin(GL_LINES);
        for(size_t l=0;l<bfaces.size();++l)
        {
            const Mesh3::Element & K=(*Th)[bfaces[l]/4];
            int e=bfaces[l]%4;
            for(int j=0;j<3;++j)
            {
                const R3 & P=K[nvfaceTet[e][j]], & Q=K[nvfaceTet[e][(j+1)%3]];
                glVertex3d(P.x,P.y,P.z);
                glVertex3d(Q.x,Q.y,Q.z);
            }
        }
        glEnd();
    }
}

static int CutTet(const R3 K[4],const R d[4],const R f[4],R3 P[4],R g[4])
{ //  the intersection of the tet K with the plane d = 0 : 0, 3 or 4 points P (in the order of the polygon)
  //  and g the value of f at the points P
    int nm[4],np[4],km=0,kp=0;
    for(int i=0;i<4;++i)
        if(d[i]<0) nm[km++]=i;
        else np[kp++]=i;
    if(km==0 || kp==0) return 0;
    //  the edges (a,b) with d[a] < 0 <= d[b]
    int ea[4],eb[4],n=0;
    if(km==1)
        for(int j=0;j<3;++j) ea[n]=nm[0],eb[n++]=np[j];
    else if(kp==1)
        for(int j=0;j<3;++j) ea[n]=nm[j],eb[n++]=np[0];
    else
    {
        ea[0]=nm[0],eb[0]=np[0];
        ea[1]=nm[0],eb[1]=np[1];
        ea[2]=nm[1],eb[2]=np[1];
        ea[3]=nm[1],eb[3]=np[0];
        n=4;
    }
    for(int i=0;i<n;++i)
    {
        int a=ea[i],b=eb[i];
        R l=d[a]/(d[a]-d[b]);
        P[i]=K[a]*(1.-l)+K[b]*l;
        g[i]=f[a]*(1.-l)+f[b]*l;
    }
    return n;
}

void OnePlotFE3::DrawSection(OneWindow *win,R3 N,R c,bool plotmesh,R thr)
{ //  the function on the section of the mesh by the plane (N,P) = c: iso lines or filled iso
  //  P0 function: one color by element, and only the elements with a value >= thr
    ThePlot & plot=*win->theplot;
    int nsubV=Psub.N();
    int nsubT=Ksub.N()/4;
    int nK=v.N()/ Th->nt;
    KN<R3> Pn(nsubV);
    KN<R> dn(nsubV);
    R3 Nt[3]={N,N,N};
    R3 P[4];
    R g[4];
    glDisable(GL_LIGHTING);
    //  the filled polygons are pushed back, so the lines on the section are not hidden
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.,1.);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    for(int k=0,o=0;k<Th->nt;++k, o+= nK)
    {
        const Mesh3::Element & K=(*Th)[k];
        int km=0;
        for(int i=0;i<4;++i)
            km += ((N,(R3) K[i])-c < 0);
        if(km==0 || km==4) continue;// no intersection
        if(p0)
        {
            if(v[o]<thr) continue;
            R3 Pt[4];
            R d[4],ff[4]={0.,0.,0.,0.};
            for(int i=0;i<4;++i)
            {
                Pt[i]=K[i];
                d[i]=(N,Pt[i])-c;
            }
            int n=CutTet(Pt,d,ff,P,g);
            plot.color(ColorOfValue(plot,v[o]));
            glBegin(GL_POLYGON);
            for(int i=0;i<n;++i)
                glVertex3d(P[i].x,P[i].y,P[i].z);
            glEnd();
            continue;
        }
        for(int i=0;i<nsubV;++i)
        {
            Pn[i]=K(Psub[i]);
            dn[i]=(N,Pn[i])-c;
        }
        for(int sk=0;sk<nsubT;++sk)
        {
            R3 Pt[4];
            R ff[4],d[4];
            for(int l=0;l<4;++l)
            {
                int iv= Ksub[4*sk+l];
                Pt[l]=Pn[iv];
                ff[l]=v[o+iv];
                d[l]=dn[iv];
            }
            int n=CutTet(Pt,d,ff,P,g);
            if(n==0) continue;
            if(win->Pretty3d())
            { // pretty rendering: continuous colors (texture of the colormap)
                R fmn,fmx;
                DataRange(plot,fmn,fmx);
                glEnable(GL_TEXTURE_1D);
                glBindTexture(GL_TEXTURE_1D,win->cmaptex);
                glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
                glBegin(GL_POLYGON);
                for(int i=0;i<n;++i)
                {
                    glTexCoord1d((g[i]-fmn)/(fmx-fmn));
                    glVertex3d(P[i].x,P[i].y,P[i].z);
                }
                glEnd();
                glDisable(GL_TEXTURE_1D);
                continue;
            }
            if(!plot.fill)
            { // white polygon to hide what is behind the section
                plot.color(0);
                glBegin(GL_POLYGON);
                for(int i=0;i<n;++i)
                    glVertex3d(P[i].x,P[i].y,P[i].z);
                glEnd();
            }
            for(int i=2;i<n;++i)
            {
                R3 Q[3]={P[0],P[i-1],P[i]};
                R gg[3]={g[0],g[i-1],g[i]};
                if(plot.fill)
                    plot.DrawIsoTfill( Q, gg, Nt, plot.Viso,plot.Viso.N(),false,3);
                else
                    plot.DrawIsoT( Q, gg, Nt, plot.Viso,plot.Viso.N(),false,3);
            }
        }
    }
    glDisable(GL_POLYGON_OFFSET_FILL);
    if(plotmesh)
    { // the section of the elements
        glLineWidth(kscreenscale);
        plot.color(1);
        for(int k=0;k<Th->nt;++k)
        {
            if(p0 && v[k*nK]<thr) continue;
            const Mesh3::Element & K=(*Th)[k];
            R3 Pt[4];
            R d[4],ff[4]={0.,0.,0.,0.};
            for(int i=0;i<4;++i)
            {
                Pt[i]=K[i];
                d[i]=(N,Pt[i])-c;
            }
            int n=CutTet(Pt,d,ff,P,g);
            if(n==0) continue;
            glBegin(GL_LINE_LOOP);
            for(int i=0;i<n;++i)
                glVertex3d(P[i].x,P[i].y,P[i].z);
            glEnd();
        }
    }
}

void OnePlotFE3::DrawIsoSurfaces(OneWindow *win)
{
    ThePlot & plot=*win->theplot;
    //  the values of the isosurfaces: the iso values given by the user (viso=),
    //  otherwise an even number of values strictly inside the bound of the function
    //  (an even number to not get the middle value, often the value on the boundary)
    int N=plot.Viso.N();
    vector<R> vals;
    vector<int> cols;
    if(plot.pViso)
        for(int i=0;i<N;++i)
        {
            vals.push_back(plot.Viso[i]);
            cols.push_back(i+4);
        }
    else if(win->Pretty3d())
    { // pretty rendering: N/5 values (keys n, N: 4 by default) evenly spaced in the range,
      // without the value on the boundary if the function is constant on the boundary
        R fmn,fmx;
        DataRange(plot,fmn,fmx);
        int ns=max(1,N/5);
        for(int k=0;k<ns;++k)
        {
            R val=fmn+(k+0.5)*(fmx-fmn)/ns;
            if(bconstant && fabs(val-bval) < 0.02*(fmx-fmn)) continue;
            vals.push_back(val);
            cols.push_back(4);
        }
    }
    else
    {
        R d=(plot.Viso[N-1]-plot.Viso[0])/(N-1);
        R fmn=plot.Viso[0]+d/2, fmx=plot.Viso[N-1]-d/2;
        int ns=2*max(1,(N+4)/8);
        for(int i=0;i<ns;++i)
        {
            R val=fmn+(i+0.5)*(fmx-fmn)/ns;
            vals.push_back(val);
            cols.push_back(5+dichotomie(plot.Viso,val));// the color of the iso value
        }
    }
    isovals=vals;
    bool pretty=win->Pretty3d();
    glPolygonMode(GL_FRONT_AND_BACK,(plot.fill || pretty) ? GL_FILL : GL_LINE);
    int nsubV=Psub.N();
    int nsubT=Ksub.N()/4;
    KN<R3> Pn(nsubV);
    int nK=v.N()/ Th->nt;
    KN<R3> Gn(nsubV);// pretty: the gradient at the sub vertices (normal of the isosurfaces)
    R fmn=0,fmx=1;
    if(pretty)
    {
        ComputeGradV();
        DataRange(plot,fmn,fmx);
        glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D,win->cmaptex);
        glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
        //  translucent isosurfaces (they do not hide each other)
        R alpha = vals.size()>1 ? 0.6 : 1.;
        glColor4f(1.,1.,1.,alpha);
        if(alpha<1)
        {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
            glDepthMask(GL_FALSE);
        }
    }
    {
    FFBatch batch(1,pretty || plot.fill);// one glBegin for all the isosurfaces (not in wireframe)
    for(int k=0,o=0;k<Th->nt;++k, o+= nK)
    {
        const Mesh3::Element & K=(*Th)[k];
        R ff[4];
        R3 Pt[4],Gt[4];
        for(int i=0;i<nsubV;++i)
            Pn[i]=K(Psub[i]);
        if(pretty)
            for(int i=0;i<nsubV;++i)
            {
                R3 Q=Psub[i];
                R lambda[4]={1.-Q.x-Q.y-Q.z,Q.x,Q.y,Q.z};
                Gn[i]=R3(0.,0.,0.);
                for(int l=0;l<4;++l) Gn[i] += gradv[(*Th)(k,l)]*lambda[l];
            }
        int lK=0;
        for(int sk=0;sk<nsubT;++sk)
        {
            for(int l=0;l<4;++l)
            {
                int iv= Ksub[lK++];
                Pt[l]=Pn[iv];
                ff[l]=v[o+iv];
                Gt[l]=Gn[iv];
            }
            R tmn=min(min(ff[0],ff[1]),min(ff[2],ff[3])), tmx=max(max(ff[0],ff[1]),max(ff[2],ff[3]));
            for(size_t i=0;i<vals.size();++i)
                if(tmn<=vals[i] && vals[i]<=tmx)
                {
                    if(pretty)
                    {
                        glTexCoord1d((vals[i]-fmn)/(fmx-fmn));
                        drawisoTetN(Pt,ff,Gt,vals[i]);
                    }
                    else
                    {
                        plot.color(cols[i]);
                        drawisoTet( Pt,ff,vals[i]);
                    }
                }
        }
    }
    }
    if(pretty)
    {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
    }
}

//  a vector field (nc components) constant on each element (P0): then one arrow by element,
//  at its barycenter, and not one arrow at each vertex of the element
static bool VectorP0(const KN<double> & v,int nt,int nc)
{
    if(nt<=0 || v.N()%nt) return false;
    int nK=v.N()/nt;
    if(nK<=nc) return false;
    R eps=1e-10*max(v.max(),-v.min());
    for(int k=0,o=0;k<nt;++k, o+= nK)
        for(int i=nc;i<nK;i+=nc)
            for(int c=0;c<nc;++c)
                if(fabs(v[o+i+c]-v[o+c]) > eps) return false;
    return true;
}

void OnePlotFE3::DrawArrows(OneWindow *win)
{ //  the vector field: one arrow by vertex, the length is proportional to the norm
  //  and the color is given by the norm
    ThePlot & plot=*win->theplot;
    int nsubV=Psub.N();
    int nK=v.N()/ Th->nt;
    R lmax=sqrt(plot.vmax2);
    if( !(lmax>0) ) return;
    R kk = 4*win->hpixel;// length of the head
    R cc = win->hpixel*40*plot.coeff;// length of the greatest arrow
    if(plot.ArrowSize>0) kk=win->hpixel*max(win->width*plot.ArrowSize/100.,1.);
    KN<bool> done(Th->nv);// to plot only one arrow by vertex of the mesh
    done=false;
    bool pretty=win->Pretty3d();
    glDisable(GL_LIGHTING);
    glLineWidth((pretty ? 1.6 : 1.)*kscreenscale);
    bool vp0=VectorP0(v,Th->nt,3);// P0 field: one arrow at the barycenter of each tetrahedron
    glBegin(GL_LINES);
    for(int k=0,o=0;k<Th->nt;++k, o+= nK)
    {
        const Mesh3::Element & K=(*Th)[k];
        for (int i=0,j=0;i<(vp0 ? 1 : nsubV);++i,j+=3)
        {
            R3 Q= vp0 ? R3(0.25,0.25,0.25) : Psub[i];
            R lambda[4]={1.-Q.x-Q.y-Q.z,Q.x,Q.y,Q.z};
            int iv=-1;
            for(int l=0;l<4;++l)
                if(fabs(lambda[l]-1.)<1e-6) iv=(*Th)(k,l);
            if(iv>=0)
            {
                if(done[iv]) continue;
                done[iv]=true;
            }
            R3 uv(v[o+j],v[o+j+1],v[o+j+2]);
            R l = uv.norme();
            if(l < 1e-4*lmax) continue;
            if(pretty)
            { // the color of the norm in the colormap
                float rgb[3];
                ColormapRGB((plot.grey ? 4 : win->cmapkind),l/lmax,rgb);
                glColor3fv(rgb);
            }
            else
                plot.color(4+dichotomie(plot.Varrow,l));
            R3 d=uv/l;
            R L=max(cc*l/lmax,kk);
            R3 P=K(Q), T=P+d*L;
            //  the head: 4 segments in 2 orthogonal planes
            R h=min(kk,L/2);
            R3 a = fabs(d.x)<0.6 ? R3(1.,0.,0.) : R3(0.,1.,0.);
            R3 n1=d^a;
            n1 /= n1.norme();
            R3 n2=d^n1;
            win->Seg3(P,T);
            win->Seg3(T,T-d*h+n1*(h/2));
            win->Seg3(T,T-d*h-n1*(h/2));
            win->Seg3(T,T-d*h+n2*(h/2));
            win->Seg3(T,T-d*h-n2*(h/2));
        }
    }
    glEnd();
}

void OnePlotFE3::Draw(OneWindow *win)
{
    initlist();

    ThePlot & plot=*win->theplot;
    ShowGlerror("begin OnePlotFE3 plot");
    bool scalar = what%10==6;
    //  scalar: the function on the boundary of the mesh, or the isosurfaces,
    //  by default the isosurfaces if the function is constant on the boundary (key s to switch)
    if(scalar && (win->changeiso || !bfacesok))
        BoundaryConstant();
    bool isosurf = win->isosurf != bconstant;
    showisovals = scalar && isosurf && !p0 && win->Pretty3d();
    //  the mesh is asked with the key m, or with a plot of this mesh: plot(Th,u)
    bool plotmesh = plot.drawmeshes;
    if(scalar)
        for (list<OnePlot *>::iterator i= plot.plots.begin();i != plot.plots.end() && !plotmesh; ++i)
        {
            OnePlotMesh3 * pm=dynamic_cast<OnePlotMesh3 *>(*i);
            plotmesh = pm && pm->Th == Th;
        }
    //  vector: the arrows
    int mode = scalar ? 1+isosurf+2*plot.fill+4*plotmesh+8*plot.drawborder+16*win->Pretty3d() : 100+plotmesh;
    if(scalar && win->Pretty3d()) win->ColormapTexture();// before the display lists
    bool  change= scalar ? win->changeiso : win->changearrow;
    //  P0 function with the key s: only the elements with a value >= thr
    R thr=-dinfty;
    if(scalar && p0 && isosurf)
    {
        int n=plot.Viso.N();
        R d=(plot.Viso[n-1]-plot.Viso[0])/(n-1);
        R fmn=plot.Viso[0]+d/2, fmx=plot.Viso[n-1]-d/2;
        thr = fmn + win->threshold*(fmx-fmn);
    }
    if(thr != lastthreshold) change=true;
    lastthreshold=thr;
    //  P0 function: exploded view (key e)
    R explode = (scalar && p0) ? win->explode : 0;
    if(explode != lastexplode) change=true;
    lastexplode=explode;
    if(mode != lastmode)
    {
        change=true;
        oklist[2]=0;// the border depend of plotmesh
    }
    lastmode=mode;

    glEnable(GL_DEPTH_TEST);
    //  cut plane (N,P) = c: the plot is clipped, we keep the part behind the plane seen from the camera
    R3 N;
    R c=0;
    int cutaxis=win->SetCutPlane(Th->Pmin,Th->Pmax,N,c);
    if(cutaxis) glEnable(GL_CLIP_PLANE0);
    //  the filled isosurfaces are lighted by default (key l to switch)
    bool withlight=win->withlight;
    if(scalar && ((isosurf && (plot.fill || p0)) || explode>0 || (win->Pretty3d() && !p0)) )
        win->withlight = !withlight;
    win->setLighting();
    win->withlight=withlight;
    if(oklist[0] && !change)
        glCallList(gllists+0);
    else
    {
        oklist[0]=1;
        glNewList(gllists+0,GL_COMPILE_AND_EXECUTE); // save  la list sans affichage
        if(!scalar)
            DrawArrows(win);
        else if(isosurf || explode>0)
        {
            if(p0)
                DrawThreshold(win,thr,plotmesh,explode);
            else
                DrawIsoSurfaces(win);
            if((!plotmesh || p0) && plot.drawborder && !(explode>0)) DrawBoundaryEdges(plot);
        }
        else
            DrawBoundary(win,plotmesh);
        glEndList();  // fin de la list
    }
    win->unsetLighting();
    glDisable(GL_CLIP_PLANE0);
    if(scalar && cutaxis && !(explode>0))
    { // the function on the section (not with the exploded view: all the faces of the elements are plotted)
        int cutmode = plot.fill+2*plotmesh+4*isosurf+8*win->Pretty3d();
        if(oklist[1] && !change && cutaxis==lastcutaxis && win->cutpos==lastcutpos && cutmode==lastcutmode
           && N.x==lastcutN.x && N.y==lastcutN.y && N.z==lastcutN.z)
            glCallList(gllists+1);
        else
        {
            oklist[1]=1;
            lastcutaxis=cutaxis;
            lastcutpos=win->cutpos;
            lastcutN=N;
            lastcutmode=cutmode;
            glNewList(gllists+1,GL_COMPILE_AND_EXECUTE);
            DrawSection(win,N,c,plotmesh,thr);
            glEndList();
        }
    }
    ShowGlerror("b mesh  OnePlotFE plot");
    //  the faces of the border: with the arrows, or with the isosurfaces and the mesh
    if(!scalar || (isosurf && plotmesh && !p0) )
    {
        if(cutaxis) glEnable(GL_CLIP_PLANE0);// the border is also clipped
        if(win->Pretty3d())
        { // pretty rendering: the mesh in dark grey, or the outline of the domain
            if(plotmesh) DrawBoundaryMesh();
            else if(plot.drawborder) DrawBoundaryEdges(plot);
        }
        else
            Plot(*Th,false,plotmesh,plot.drawborder,plot,gllists+2,&oklist[2]);
        glDisable(GL_CLIP_PLANE0);
    }
    ShowGlerror("OnePlotFE::Draw");
}

void OnePlotFES::Draw(OneWindow *win)
{

    initlist();
    ThePlot & plot=*win->theplot;
    ShowGlerror("begin OnePlotFE plot");
    //  pretty rendering: continuous colors and smooth shading (lighting by default)
    bool pretty = win->Pretty3d() && what%10==8 && !win->changePlotdim;
    if(pretty) win->ColormapTexture();// before the display lists
    {
        bool wl=win->withlight;
        if(pretty) win->withlight=!wl;
        win->setLighting();
        win->withlight=wl;
    }
    assert(win);
    const MeshS & Th(*this->Th);
    int nsubT= Ksub.N()/3;//NbOfSubTriangle(nsub);
    int nsubV=Psub.N();//NbOfSubInternalVertices(nsub);
    int nK=v.N()/ Th.nt;
    if(debug>4)
        cout << "\t\t\tOnePlotMesh::Draw  " <<v.N() << " ,nt " << Th.nt << " " << nK << " "
        << Psub.N() << " " << what << " ,nv " << Th.nv << " cas=" << cas << endl;
    ffassert(v.N()== Th.nt*nK);
    int o=0;
    KN<R3> Pn(Psub.N());
   
    if((debug > 10)) cout << " " <<nsubV  << " " << nsubT << endl;

    // compute Normal at vertice
    KN<R3> Nv(Th.nv);
    for (int it=0 ; it<Th.nv; it++) //
        Nv[it]=R3(0.,0.,0.); //
    for (int it=0 ; it<Th.nt; it++) {
        const MeshS::Element &K = (Th[it]);
           int iv;
        for (int i = 0; i < 3; i++) {
            iv = Th.operator( )(K[i]);
        Nv[iv]+=K.mesure()*K.NormalT();
        }
    }
    for (int i=0 ; i<Th.nv; i++)
        Nv[i]/=Nv[i].norme();
            
  if(plot.fill && what%10==8)
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    else
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

    if(what%10==9)
        glDisable(GL_DEPTH_TEST);
    else
        glEnable(GL_DEPTH_TEST);

    R coef = plot.coeff;
    double xmin,xmax,ymin,ymax;
    win->getcadre(xmin,xmax,ymin,ymax);
    R cc = win->hpixel*40;
    R kk = 4*win->hpixel;
    if(plot.ArrowSize>0) kk=win->hpixel*max(win->width*plot.ArrowSize/100.,1.);

    int klist=0;
    bool  change=false;
    if( (what%10==8) )
    {
        if ( plot.fill) klist=1;
        change = win->changeiso ;
    }
    else if (what%10==9)
        change  =  win->changearrow ;
    if(debug>9)
        cout << change << " " << klist << " ... " << oklist[klist] << "  fill = "
        <<  plot.fill <<  " " <<  coef << endl;
    if (oklist[klist] && ! change )
        glCallList(gllists+klist);
    else
    {
        oklist[klist]=1;
        glNewList(gllists+klist,GL_COMPILE_AND_EXECUTE); // save  la list aevc  affichage
        if(debug>100)
            cout << win->Bmin << ", Bmax:   " << win->Bmax << " Viso: "<< plot.Viso << endl;
        if(pretty)
        { // the normal at the corners: average on the triangles with an angle < 40 degree
            vector<R3> fn(Th.nt);
            vector<vector<int> > vt(Th.nv);
            for(int k=0;k<Th.nt;++k)
            {
                const MeshS::Element & K=Th[k];
                //  the normal of the order of the vertices (as OpenGL, for the two sided lighting)
                fn[k]=R3((R3) K[0],(R3) K[1])^R3((R3) K[0],(R3) K[2]);
                for(int i=0;i<3;++i) vt[Th(K[i])].push_back(k);
            }
            R fmn,fmx;
            DataRange(plot,fmn,fmx);
            glEnable(GL_TEXTURE_1D);
            glBindTexture(GL_TEXTURE_1D,win->cmaptex);
            glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
            glColor4f(1.,1.,1.,1.);
            glBegin(GL_TRIANGLES);
            for(int k=0;k<Th.nt;++k, o+= nK)
            {
                const MeshS::Element & K=Th[k];
                R3 nk=fn[k]/fn[k].norme();
                R3 nc[3];
                for(int i=0;i<3;++i)
                {
                    const vector<int> & tl=vt[Th(K[i])];
                    R3 n(0.,0.,0.);
                    for(size_t m=0;m<tl.size();++m)
                        if((fn[tl[m]],nk) > 0.766*fn[tl[m]].norme()) n += fn[tl[m]];
                    nc[i]=n/n.norme();
                }
                for(int sk=0;sk<nsubT;++sk)
                    for(int j=0;j<3;++j)
                    {
                        int ii=Ksub[sk*3+j];
                        R2 Q=Psub[ii];
                        R3 n=nc[0]*(1.-Q.x-Q.y)+nc[1]*Q.x+nc[2]*Q.y;
                        n/=n.norme();
                        R3 P=K(Q);
                        glNormal3d(n.x,n.y,n.z);
                        glTexCoord1d((v[o+ii]-fmn)/(fmx-fmn));
                        glVertex3d(P.x,P.y,P.z);
                    }
            }
            glEnd();
            glDisable(GL_TEXTURE_1D);
        }
        else
        {
        FFBatch batch((what%10==8 && plot.fill) ? 1 : 2);// one glBegin for all the elements
        bool vp0 = what%10==9 && VectorP0(v,Th.nt,3);// P0 field: one arrow by triangle
        for(int k=0;k<Th.nt;++k, o+= nK)
        {
            const typename MeshS::Element & K=Th[k];
            int iv[3];
            for (int i = 0; i < 3; i++)
            iv[i] = Th.operator( )(K[i]);
                
            for(int i=0;i<nsubV;++i)
                Pn[i]=K(Psub[i]);// local to global coord.
             
            if(what%10==8)
                for(int sk=0;sk<nsubT;++sk)
                {
                    int i0= Ksub[sk*3+0];//numSubTriangle(nsub,sk,0);
                    int i1= Ksub[sk*3+1];//numSubTriangle(nsub,sk,1);
                    int i2= Ksub[sk*3+2];//numSubTriangle(nsub,sk,2);
               
                    R ff[3]={v[o+i0],v[o+i1],v[o+i2]};
                    R3 Pt[3]={Pn[i0],Pn[i1],Pn[i2]};
                    R3 Nt[3]={Nv[iv[0]],Nv[iv[1]],Nv[iv[2]]};
                    if(plot.fill)
                       plot.DrawIsoTfill( Pt, ff, Nt, plot.Viso,plot.Viso.N(),win->changePlotdim,win->viewdim);
                    else
                       plot.DrawIsoT( Pt, ff, Nt, plot.Viso,plot.Viso.N(),win->changePlotdim,win->viewdim);
                }
            else if (what%10==9)
                for (int i=0,j=0;i<(vp0 ? 1 : nsubV);++i) {
                    R3 P= vp0 ? K(R2(1./3.,1./3.)) : Pn[i];// P0 field: at the barycenter
                    R3 uv(v[o+j],v[o+j+1],v[o+j+2]);
                    j+=3;
                    R l = Max(sqrt((uv,uv)),1e-30) ;
                    int col = 2+dichotomie(plot.Varrow,l);
                    if(debug>100)
                        cout << uv << " l= " << l << " " << coef << " " <<col <<  endl;
                    plot.color(2+col);
                    uv *=coef;
                    
                    l *=coef;
                    R3 dd=uv*(-kk/l);
                    R ndd=dd.norme();
                    dd/=ndd;
                    R3 nx(1.,0.,0.), ny(0.,1.,0.), nz(0.,0.,1.);
                    R3 dnx = (dd^nx)*0.5, dny = (dd^ny)*0.5, dnz = (dd^nz)*0.5;
                   
               
                   dnx *= -coef*ndd/dnx.norme();
                   dny *= -coef*ndd/dny.norme();
                   dnz *= -coef*ndd/dnz.norme();
                    
                    if (l*10000.< kk) continue;
                    if (l < kk)
                        uv *= (kk/l)*coef;
                    else if (l> cc)
                        uv *= (cc/l)*coef;
                   
                    dd*=ndd*coef; // denormalize
                    
                    
                    if(ffbatch!=2) glBegin(GL_LINES);
                    win->Seg3(P,P+uv);

                    if (10*l>kk) {
                        win->Seg3(P+uv,P+uv+dd+dnx);
                        win->Seg3(P+uv,P+uv+dd-dnx);
                        win->Seg3(P+uv,P+uv+dd+dny);
                        win->Seg3(P+uv,P+uv+dd-dny);
                        win->Seg3(P+uv,P+uv+dd+dnz);
                        win->Seg3(P+uv,P+uv+dd-dnz);
                    }
                    if(ffbatch!=2) glEnd();
                }
        }
        }
        glEndList();  // fin de la list

}
    ShowGlerror("b mesh  OnePlotFES plot");
    win->unsetLighting();
    Plot(Th,false,plot.drawmeshes,plot.drawborder,plot,gllists+2,&oklist[2],win);
    ShowGlerror("OnePlotFES::Draw");
}

void OnePlotFEL::Draw(OneWindow *win)
{
    if(win->Graph1d())
    {
        DrawGraph(win);
        return;
    }
    initlist();
    ThePlot & plot=*win->theplot;
    ShowGlerror("begin OnePlotFE plot");
    win->setLighting();
    assert(win);
    const MeshL & Th(*this->Th);
    int nsubT= Ksub.N()/2;
    int nsubV=Psub.N();
    int nK=v.N()/ Th.nt;
    if(debug>4)
        cout << "\t\t\tOnePlotMesh::Draw  " <<v.N() << " ,nt " << Th.nt << " " << nK << " "
        << Psub.N() << " " << what << " ,nv " << Th.nv << " cas=" << cas << endl;
    ffassert(v.N()== Th.nt*nK);
    int o=0;
    KN<R3> Pn(Psub.N());
    if((debug > 10)) cout << " " <<nsubV  << " " << nsubT << endl;
    
    // compute Normal at vertice
    KN<R3> Nv(Th.nv);
    for (int it=0 ; it<Th.nv; it++) //
        Nv[it]=R3(0.,0.,0.); //
    for (int it=0 ; it<Th.nt; it++) {
        const MeshL::Element &K = (Th[it]);
           int iv;
        for (int i = 0; i < 2; i++) {
            iv = Th.operator( )(K[i]);
        Nv[iv]+=K.mesure()*K.NormalT();
        }
    }
    for (int i=0 ; i<Th.nv; i++)
        Nv[i]/=Nv[i].norme();
    
    if(plot.fill && (what == 14 || what == 20 || what == 114 || what == 120))
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    else
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    
    if(what==15 || what==21)
        glDisable(GL_DEPTH_TEST);
    else
        glEnable(GL_DEPTH_TEST);
    
    R coef = plot.coeff;
    double xmin,xmax,ymin,ymax;
    win->getcadre(xmin,xmax,ymin,ymax);
    R kk = 4*win->hpixel;
    if(plot.ArrowSize>0) kk=win->hpixel*max(win->width*plot.ArrowSize/100.,1.);
    
    int klist=0;
    bool  change=false;
    if(what == 14 || what == 20 || what == 114 || what == 120) {
        if ( plot.fill) klist=1;
        change = win->changeiso ;
    }
    else if (what==15 || what==21)
        change  =  win->changearrow ;
    if(debug>9)
        cout << change << " " << klist << " ... " << oklist[klist] << "  fill = "
        <<  plot.fill <<  " " <<  coef << endl;
    if (oklist[klist] && ! change )
        glCallList(gllists+klist);
    else {
        oklist[klist]=1;
        glNewList(gllists+klist,GL_COMPILE_AND_EXECUTE); // save  la list aevc  affichage
        if(debug>100)
            cout << win->Bmin << ", Bmax:   " << win->Bmax << " Viso: "<< plot.Viso << endl;
        {
        FFBatch batch(2);// one glBegin for all the elements
        for(int k=0;k<Th.nt;++k, o+= nK) {
            const typename MeshL::Element & K=Th[k];
            int iv[2];
            for (int i = 0; i < 2; i++)
                iv[i] = Th.operator( )(K[i]);
           for(int i=0;i<nsubV;++i)
                Pn[i]=K(Psub[i]);// local to global coord.
            if(what == 14 || what == 20 || what == 114 || what == 120) {
                for(int sk=0;sk<nsubT;++sk) {
                    int i0= Ksub[sk*2+0], i1=Ksub[sk*2+1];
                    R ff[2]={v[o+i0],v[o+i1]};
                    R3 Pt[2]={Pn[i0],Pn[i1]};
                    R3 Nt[2]={Nv[iv[0]],Nv[iv[1]]};
                    plot.DrawIsoEfill( Pt, ff, Nt, plot.Viso,plot.Viso.N(),win->changePlotdim,win->viewdim);
                }
            }
            else
                ffassert(0);
            
        }
        }
        glEndList();  // fin de la list
        
    }
    ShowGlerror("b mesh  OnePlotFEL plot");
    win->unsetLighting();
    Plot(Th,false,plot.drawmeshes,plot.drawborder,plot,gllists+2,&oklist[2],win);
    ShowGlerror("OnePlotFEL::Draw");
}


template<class Mesh>
OnePlotFE<Mesh>::OnePlotFE(const Mesh *T,long w,PlotStream & f)
:OnePlot(w,2,5),Th(T),cas(2),showisovals(false)
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
    {      f >> Psub ;
        f >> Ksub ;
        if(debug>2) {
            cout << " Psub " << Psub << endl;
            cout << " Ksub " << Ksub << endl;}

    }
    if(what<10)
        f>> v;
    else
    {
        f>> vc;
        vc2v();
    }
    cas =2;
    vc2v();
    if(debug>3) cout << "OnePlotFE" << Th <<" " << what<< " " << Psub.N() << " " << Ksub.N()/3 <<" " << v.N() << " cas " << cas << endl;
    ffassert(f.good());

}



template<class Mesh>
bool  OnePlotFE<Mesh>::vc2v()
{
    bool ret=false;
    if(what>=10)
    {
        ret=true;
        int n = vc.N();
        if(v.size() !=n)
            v.resize(n);
        for(int i=0;i<n;++i)
            if( cas%4== 0) v[i] = vc[i].real();
            else if( cas%4== 1) v[i] = vc[i].imag();
            else if( cas%4== 2) v[i] = abs(vc[i]);
            else if( cas%4== 3) v[i] = arg(vc[i]);
    }

    if(what%10==1)
    {
        fmin = min(fmin,v.min());
        fmax = max(fmax,v.max());
    }
    else if (what%10==2)
    {
        int n= v.N()/2;
        for (int i=0,j=0;i<n;i++, j+=2)
        {
            R2 u(v[j],v[j+1]);
            vmax2 = max(vmax2,u.norme2());
        }
    }
    return ret;
}

bool  OnePlotFE3::vc2v()
{
    bool ret=false;
    if(what>=10)
    {
        ret=true;
        int n = vc.N();
        if(v.size() !=n)
            v.resize(n);
        for(int i=0;i<n;++i)
            if( cas%4== 0) v[i] = vc[i].real();
            else if( cas%4== 1) v[i] = vc[i].imag();
            else if( cas%4== 2) v[i] = abs(vc[i]);
            else if( cas%4== 3) v[i] = arg(vc[i]);
    }


    if(what%10==6)
    {
        fmin = min(fmin,v.min());
        fmax = max(fmax,v.max());
    }
    else if (what%10==7)
    {

        int n= v.N()/3;
        for (int i=0,j=0;i<n;i++, j+=3)
        {
            R3 u(v[j],v[j+1],v[j+2]);
            vmax2 = max(vmax2,u.norme2());
        }
    }
    return ret;
}




bool  OnePlotFES::vc2v()
{
    bool ret=false;
    if(what>=10)
    {
        ret=true;
        int n = vc.N();
        if(v.size() !=n)
            v.resize(n);
        for(int i=0;i<n;++i)
            if( cas%4== 0) v[i] = vc[i].real();
            else if( cas%4== 1) v[i] = vc[i].imag();
            else if( cas%4== 2) v[i] = abs(vc[i]);
            else if( cas%4== 3) v[i] = arg(vc[i]);
    }


    if(what%10==8)
    {
        fmin = min(fmin,v.min());
        fmax = max(fmax,v.max());
    }
    else if (what%10==9)
    {

        int n= v.N()/3;
        for (int i=0,j=0;i<n;i++, j+=3)
        {
            R3 u(v[j],v[j+1],v[j+2]);
            vmax2 = max(vmax2,u.norme2());
        }
    }
    return ret;
}


bool  OnePlotFEL::vc2v()
{
    bool ret=false;
    if(what>=16)
    {
        ret=true;
        int n = vc.N();
        if(v.size() !=n)
            v.resize(n);
        for(int i=0;i<n;++i)
            if( cas%4== 0) v[i] = vc[i].real();
            else if( cas%4== 1) v[i] = vc[i].imag();
            else if( cas%4== 2) v[i] = abs(vc[i]);
            else if( cas%4== 3) v[i] = arg(vc[i]);
    }
    
    
    if(what == 14 || what == 20 || what == 114 || what == 120)
    {
        fmin = min(fmin,v.min());
        fmax = max(fmax,v.max());
    }
    else if (what==15 || what==21)
    {
        
        int n= v.N()/3;
        for (int i=0,j=0;i<n;i++, j+=3)
        {
            R3 u(v[j],v[j+1],v[j+2]);
            vmax2 = max(vmax2,u.norme2());
        }
    }
    return ret;
}
template<class Mesh>
void OnePlotFE<Mesh>::Draw(OneWindow *win)
{
    initlist();
    ThePlot & plot=*win->theplot;
    ShowGlerror("begin OnePlotFE plot");
    //  pretty rendering: continuous colors (fill), smooth shading with the elevation (3d view)
    bool pretty = win->PrettyOn() && what%10==1;
    bool prettyfill = pretty && plot.fill;
    bool elev = win->plotdim==3;
    //  pretty rendering: the values of the iso lines (marks on the color bar), keys n, N:
    //  with fill, N/5 iso lines on the colors (4 by default), otherwise the iso values
    showisovals=pretty && (!prettyfill || win->isolines);
    if(pretty)
    {
        R fmn,fmx;
        DataRange(plot,fmn,fmx);
        isovals.clear();
        int N=plot.Viso.N();
        if(prettyfill)
        {
            int ns=max(1,N/5);
            if(plot.pViso)
                for(int i=0;i<N;++i) isovals.push_back(plot.Viso[i]);
            else
                for(int k=0;k<ns;++k) isovals.push_back(fmn+(k+0.5)*(fmx-fmn)/ns);
        }
        else
            for(int i=0;i<N;++i)
                if(plot.Viso[i]>=fmn && plot.Viso[i]<=fmx) isovals.push_back(plot.Viso[i]);
    }
    if(win->PrettyOn()) win->ColormapTexture();// before the display lists
    {
        bool wl=win->withlight;
        if(prettyfill && elev) win->withlight=!wl;// lighting by default
        win->setLighting();
        win->withlight=wl;
    }
    assert(win);
    const Mesh & Th(*this->Th);
    int nsubT= Ksub.N()/3;//NbOfSubTriangle(nsub);
    int nsubV=Psub.N();//NbOfSubInternalVertices(nsub);
    int nK=v.N()/ Th.nt;
    if(debug>4)
        cout << "\t\t\tOnePlotMesh::Draw  " <<v.N() << " ,nt " << Th.nt << " " << nK << " "
        << Psub.N() << " " << what << " ,nv " << Th.nv << " cas=" << cas << endl;
    ffassert(v.N()== Th.nt*nK);
    ffassert(nK == nsubV*(what%10));
    int o=0;
    KN<R2> Pn(Psub.N());
    if((debug > 10)) cout << " " <<nsubV  << " " << nsubT << endl;

    if(plot.fill && what%10==1)
        glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    else
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

    if(what%10==2)
        glDisable(GL_DEPTH_TEST);
    else
        glEnable(GL_DEPTH_TEST);

    R coef = plot.coeff;
    double xmin,xmax,ymin,ymax;
    win->getcadre(xmin,xmax,ymin,ymax);
    R kk = 4*win->hpixel;
    if(plot.ArrowSize>0) kk=win->hpixel*max(win->width*plot.ArrowSize/100.,1.);
    R cc = win->hpixel*40;

    int klist=0;
    bool  change=false;
    if( (what%10==1) )
    {
        if ( plot.fill) klist=1;
        change = win->changeiso ;
    }
    else if (what%10==2)
        change  =  win->changearrow ;
    if(debug>9)
        cout << change << " " << klist << " ... " << oklist[klist] << "  fill = "
        <<  plot.fill <<  " " <<  coef << endl;
    //  the filled surface is pushed back, so the mesh on the surface is not hidden by the surface
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.,1.);
    if (oklist[klist] && ! change )
        glCallList(gllists+klist);
    else
    {
        oklist[klist]=1;
        glNewList(gllists+klist,GL_COMPILE_AND_EXECUTE); // save  la list aevc  affichage
        if(debug>100)
            cout << win->Bmin << ", Bmax:   " << win->Bmax << " Viso: "<< plot.Viso << endl;
        if(pretty) glLineWidth(1.5*kscreenscale);// the iso lines
        if(prettyfill)
            DrawPretty(win,elev);
        else
        {
        FFBatch batch((what%10==1 && plot.fill) ? 1 : 2);// one glBegin for all the elements
        bool vp0 = what%10==2 && VectorP0(v,Th.nt,2);// P0 field: one arrow by triangle
        for(int k=0;k<Th.nt;++k, o+= nK)
        {
            const typename Mesh::Element & K=Th[k];
            for(int i=0;i<nsubV;++i)
                Pn[i]=K(Psub[i]);// local to global coord.
            if(what%10==1)
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
                for (int i=0,j=0;i<(vp0 ? 1 : nsubV);++i)
                {
                    R2 P= vp0 ? K(R2(1./3.,1./3.)) : Pn[i];// P0 field: at the barycenter
                    R2 uv(v[o+j],v[o+j+1]);
                    j+=2;
                    R  l = Max(sqrt((uv,uv)),1e-30) ;
                    int col = 2+dichotomie(plot.Varrow,l);
                    if(debug>100)
                        cout << uv << " l= " << l << " " << coef << " " <<col <<  endl;

                    if(win->PrettyOn())
                    { // pretty rendering: the color of the norm in the colormap
                        float rgb[3];
                        ColormapRGB(plot.cmap>=0 ? plot.cmap : 0,l/max(sqrt(plot.vmax2),1e-30),rgb);
                        glColor3fv(rgb);
                    }
                    else
                        plot.color(2+col);
                    uv = coef*uv;
                    l *= coef;
                    R2 dd = uv*(-kk/l);
                    R2 dn = dd.perp()*0.5;
                    if (l*10000.< kk) continue;
                    if (l < kk)
                        uv = uv*(kk/l);
                    else if (l> cc)
                        uv = uv*(cc/l);
                    if(ffbatch!=2) glBegin(GL_LINES);

                    win->Seg(P,P+uv);

                    if (10*l>kk) {
                        win->Seg(P+uv,P+uv+dd+dn);
                        win->Seg(P+uv,P+uv+dd-dn);
                    }
                    if(ffbatch!=2) glEnd();
                }

        }
        }
        glEndList();  // fin de la list
    }
    glLineWidth(kscreenscale);
    glDisable(GL_POLYGON_OFFSET_FILL);

    ShowGlerror("b mesh  OnePlotFE plot");
    win->unsetLighting();
    if(what%10==1)
    { // the mesh and the border are plotted on the surface z = f(x,y) (so correct with a 3d view)
        bool plotmesh = plot.drawmeshes;// or a plot of this mesh: plot(Th,u)
        for (list<OnePlot *>::iterator i= plot.plots.begin();i != plot.plots.end() && !plotmesh; ++i)
        {
            OnePlotMesh<Mesh> * pm=dynamic_cast<OnePlotMesh<Mesh> *>(*i);
            plotmesh = pm && pm->Th == this->Th;
        }
        if(change) oklist[2]=oklist[4]=0;
        glEnable(GL_DEPTH_TEST);
        if(plotmesh) DrawMeshOnSurface(plot,false);
        if(plot.drawborder) DrawMeshOnSurface(plot,true);
    }
    else
        Plot(Th,false,plot.drawmeshes,plot.drawborder,plot,gllists+2,&oklist[2]);
    ShowGlerror("OnePlotFE::Draw");
}

template<class Mesh>
void OnePlotFE<Mesh>::DrawMeshOnSurface(ThePlot & plot,bool border)
{ //  plot the edges of the mesh (or of the border) on the surface z = f(x,y),
  //  with the sub triangulation, to follow the surface inside the elements
    int kl = border ? 2 : 4;
    if(oklist[kl])
    {
        glCallList(gllists+kl);
        return;
    }
    const Mesh & Th(*this->Th);
    int nsubT= Ksub.N()/3;
    int nK=v.N()/ Th.nt;
    map<pair<int,int>,int> labe;// the label of the border edges
    if(border)
        for (int i=0;i<Th.nbBrdElmts();i++)
        {
            const typename  Mesh::BorderElement  & E(Th.be(i));
            int a=Th(E[0]), b=Th(E[1]);
            labe[make_pair(min(a,b),max(a,b))]=E.lab;
        }
    oklist[kl]=1;
    glNewList(gllists+kl,GL_COMPILE_AND_EXECUTE);
    glLineWidth((border ? 2 : 1)*kscreenscale);
    plot.color(1);
    glBegin(GL_LINES);
    for(int k=0,o=0;k<Th.nt;++k, o+= nK)
    {
        const typename Mesh::Element & K=Th[k];
        for (int e=0;e<3;++e)
        { // edge e of K (opposite to the vertex e)
            if(border)
            {
                int a=Th(k,(e+1)%3), b=Th(k,(e+2)%3);
                map<pair<int,int>,int>::const_iterator pe=labe.find(make_pair(min(a,b),max(a,b)));
                if(pe==labe.end()) continue;
                plot.color(1+abs(pe->second));
            }
            for(int sk=0;sk<nsubT;++sk)
                for(int j=0;j<3;++j)
                {
                    int i0= Ksub[sk*3+j], i1= Ksub[sk*3+(j+1)%3];
                    R2 Q0=Psub[i0], Q1=Psub[i1];
                    R l0[3]={1.-Q0.x-Q0.y,Q0.x,Q0.y}, l1[3]={1.-Q1.x-Q1.y,Q1.x,Q1.y};
                    if( fabs(l0[e]) > 1e-6 || fabs(l1[e]) > 1e-6 ) continue;// not on the edge e
                    R2 P0=K(Q0), P1=K(Q1);
                    glVertex3d(P0.x,P0.y,v[o+i0]);
                    glVertex3d(P1.x,P1.y,v[o+i1]);
                }
        }
    }
    glEnd();
    glLineWidth(kscreenscale);
    glEndList();
}

//  interpolation in the sub triangles (Psub, Ksub) of the reference triangle at the point Q,
//  nc = number of components (the values of the sub vertex i: v[o+nc*i+c])
static bool SubInterp2(const KN<R2> & Psub,const KN<int> & Ksub,const KN<double> & v,int o,int nc,R2 Q,R w3[3])
{
    int nsubT=Ksub.N()/3, best=-1;
    R bmin=-1e100, mu[3]={0,0,0};
    for(int sk=0;sk<nsubT;++sk)
    {
        R2 A=Psub[Ksub[3*sk]], B=Psub[Ksub[3*sk+1]], C=Psub[Ksub[3*sk+2]];
        R2 AB(A,B),AC(A,C),AQ(A,Q);
        R det=AB.x*AC.y-AB.y*AC.x;
        if(det==0) continue;
        R m1=(AQ.x*AC.y-AQ.y*AC.x)/det, m2=(AB.x*AQ.y-AB.y*AQ.x)/det, m0=1-m1-m2;
        R mn=min(m0,min(m1,m2));
        if(mn>bmin) { bmin=mn; best=sk; mu[0]=m0; mu[1]=m1; mu[2]=m2;}
    }
    if(best<0) return false;
    R sw=0;
    for(int j=0;j<3;++j) sw += (mu[j]=max(0.,mu[j]));
    for(int j=0;j<3;++j) mu[j]/=sw;
    for(int c=0;c<nc;++c)
    {
        w3[c]=0;
        for(int j=0;j<3;++j) w3[c] += mu[j]*v[o+nc*Ksub[3*best+j]+c];
    }
    return true;
}

template<class Mesh>
bool OnePlotFE<Mesh>::Probe(const R3 & P,R & val,R3 & vec,int & ncomp)
{ //  2d: the triangle which contains (P.x,P.y)
    const Mesh & Th(*this->Th);
    if(what%10!=1 && what%10!=2) return false;
    if(!grid.built)
    {
        vector<R3> mn(Th.nt),mx(Th.nt);
        for(int k=0;k<Th.nt;++k)
        {
            const typename Mesh::Element & K=Th[k];
            R2 A=K[0],B=K[1],C=K[2];
            mn[k]=R3(min(A.x,min(B.x,C.x)),min(A.y,min(B.y,C.y)),0.);
            mx[k]=R3(max(A.x,max(B.x,C.x)),max(A.y,max(B.y,C.y)),0.);
        }
        grid.build(mn,mx);
    }
    if(grid.head.empty()) return false;
    int c[3];
    grid.ijk(R3(P.x,P.y,0.),c);
    int nK=v.N()/Th.nt, nc= what%10==1 ? 1 : 2;
    for(int q=grid.head[grid.cell(c[0],c[1],0)];q>=0;q=grid.nxt[q])
    {
        int k=grid.elem[q];
        const typename Mesh::Element & K=Th[k];
        R2 A=K[0],B=K[1],C=K[2];
        R2 AB(A,B),AC(A,C),AP(A,R2(P.x,P.y));
        R det=AB.x*AC.y-AB.y*AC.x;
        if(det==0) continue;
        R l1=(AP.x*AC.y-AP.y*AC.x)/det, l2=(AB.x*AP.y-AB.y*AP.x)/det;
        if(l1<-1e-9 || l2<-1e-9 || l1+l2>1+1e-9) continue;
        R w[3];
        if(!SubInterp2(Psub,Ksub,v,k*nK,nc,R2(l1,l2),w)) return false;
        ncomp=nc;
        if(nc==1) val=w[0];
        else { vec=R3(w[0],w[1],0.); val=vec.norme();}
        return true;
    }
    return false;
}

template<class Mesh>
void OnePlotFE<Mesh>::DrawPretty(OneWindow *win,bool elev)
{ //  pretty rendering of a 2d scalar function: continuous colors (texture of the colormap);
  //  with the elevation (3d view): smooth normals from the gradient averaged at the vertices
    ThePlot & plot=*win->theplot;
    const Mesh & Th(*this->Th);
    int nsubT=Ksub.N()/3, nsubV=Psub.N(), nK=v.N()/Th.nt;
    R fmn,fmx;
    DataRange(plot,fmn,fmx);
    int iv3[3]={-1,-1,-1};// the sub vertices at the vertices of the triangle
    for(int i=0;i<nsubV;++i)
    {
        R2 Q=Psub[i];
        R lambda[3]={1.-Q.x-Q.y,Q.x,Q.y};
        for(int j=0;j<3;++j)
            if(fabs(lambda[j]-1.)<1e-6) iv3[j]=i;
    }
    vector<R2> gv;// the gradient at the vertices
    if(elev && iv3[0]>=0 && iv3[1]>=0 && iv3[2]>=0)
    {
        gv.assign(Th.nv,R2(0.,0.));
        vector<R> w(Th.nv,0.);
        for(int k=0,o=0;k<Th.nt;++k,o+=nK)
        {
            const typename Mesh::Element & K=Th[k];
            R2 A=K[0],B=K[1],C=K[2];
            R2 AB(A,B),AC(A,C);
            R det=AB.x*AC.y-AB.y*AC.x;
            if(fabs(det)<1e-300) continue;
            R da=v[o+iv3[1]]-v[o+iv3[0]], db=v[o+iv3[2]]-v[o+iv3[0]];
            R2 g((da*AC.y-db*AB.y)/det,(db*AB.x-da*AC.x)/det);
            R m=fabs(det)/2;
            for(int j=0;j<3;++j)
            {
                int iv=Th(K[j]);
                gv[iv] += g*m;
                w[iv] += m;
            }
        }
        for(int i=0;i<Th.nv;++i)
            if(w[i]>0) gv[i] = gv[i]/w[i];
    }
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D,win->cmaptex);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,elev ? GL_MODULATE : GL_REPLACE);
    if(elev) glEnable(GL_NORMALIZE);// the normals are scaled by the Z scale of the view
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glColor4f(1.,1.,1.,1.);
    glBegin(GL_TRIANGLES);
    for(int k=0,o=0;k<Th.nt;++k,o+=nK)
    {
        const typename Mesh::Element & K=Th[k];
        R2 g3[3];
        if(gv.size())
            for(int j=0;j<3;++j) g3[j]=gv[Th(K[j])];
        for(int sk=0;sk<nsubT;++sk)
            for(int j=0;j<3;++j)
            {
                int ii=Ksub[sk*3+j];
                R2 Q=Psub[ii], P=K(Q);
                R f=v[o+ii];
                if(gv.size())
                {
                    R2 g=g3[0]*(1.-Q.x-Q.y)+g3[1]*Q.x+g3[2]*Q.y;
                    glNormal3d(-g.x,-g.y,1.);
                }
                glTexCoord1d((f-fmn)/(fmx-fmn));
                glVertex3d(P.x,P.y,f);
            }
    }
    glEnd();
    glDisable(GL_TEXTURE_1D);
    if(elev) glDisable(GL_NORMALIZE);
    //  the iso lines on the colors (keys n, N, on/off with o), in dark grey
    if(!win->isolines) return;
    glDisable(GL_LIGHTING);
    glLineWidth(1.2*kscreenscale);
    glColor3f(0.12f,0.12f,0.15f);
    glBegin(GL_LINES);
    for(int k=0,o=0;k<Th.nt;++k,o+=nK)
    {
        const typename Mesh::Element & K=Th[k];
        for(int sk=0;sk<nsubT;++sk)
        {
            int ii[3]={Ksub[sk*3],Ksub[sk*3+1],Ksub[sk*3+2]};
            R ff[3]={v[o+ii[0]],v[o+ii[1]],v[o+ii[2]]};
            R tmn=Min(ff[0],ff[1],ff[2]), tmx=Max(ff[0],ff[1],ff[2]);
            for(size_t l=0;l<isovals.size();++l)
            {
                R val=isovals[l];
                if(val<tmn || val>tmx || tmn==tmx) continue;
                R2 P[3];
                int n=0;
                for(int e=0;e<3 && n<2;++e)
                {
                    int a=e, b=(e+1)%3;
                    if((ff[a]-val)*(ff[b]-val)<0 || (ff[a]==val && ff[b]!=val))
                    {
                        R t=(val-ff[a])/(ff[b]-ff[a]);
                        P[n++]=K(Psub[ii[a]])*(1.-t)+K(Psub[ii[b]])*t;
                    }
                }
                if(n==2)
                {
                    glVertex3d(P[0].x,P[0].y,val);
                    glVertex3d(P[1].x,P[1].y,val);
                }
            }
        }
    }
    glEnd();
}

template<class Mesh>
void OnePlotFE<Mesh>::dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn2,R & vmx2) const
{
    const Mesh & Th(*this->Th);
    int nsubT= Ksub.N()/3;//NbOfSubTriangle(nsub);
    int nsubV=Psub.N();//NbOfSubInternalVertices(nsub);
    int nK=v.N()/ Th.nt;
    ffassert(v.N()== Th.nt*nK);
    ffassert(nK == nsubV*(what%10));
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
            if(what%10==1) f=v[o+i];
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
                if(what%10==1 )
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
                    vmn2=min(l,vmn2);
                    vmx2=max(l,vmx2);
                }
            }
        if(debug>100 && ccc)
            cout << " dny_bfv :  "  << fmn << " " << fmx << " " << sqrt(vmn2) << " " << sqrt(vmx2)
            <<  " : " << Pn[0] << endl;

    }
}

//  ---  1d graph a la gnuplot (plot with only 2d curves) ---
inline int GraphScale() { return kscreenscale>1.5 ? 2 : 1;}
inline void * GraphFont() { return GraphScale()==2 ? GLUT_BITMAP_HELVETICA_18 : GLUT_BITMAP_HELVETICA_12;}
inline int GraphFontHeight() { return GraphScale()==2 ? 13 : 9;}

static void GraphColor(int k,bool grey)
{ // the default colors of gnuplot
    static const unsigned char c[8][3]={ {0x94,0x00,0xD3},{0x00,0x9E,0x73},{0x56,0xB4,0xE9},{0xE6,0x9F,0x00},
                                         {0xF0,0xE4,0x42},{0x00,0x72,0xB2},{0xE5,0x1E,0x10},{0x00,0x00,0x00}};
    if(grey) glColor3d(0.,0.,0.);
    else glColor3ub(c[k%8][0],c[k%8][1],c[k%8][2]);
}

static int GraphTextWidth(const char *s)
{
    int l=0;
    for (;*s;++s)
        l += glutBitmapWidth(GraphFont(),*s);
    return l;
}

static bool graphbold=true;// the texts of the 1d graph in bold

static void GraphText(R x,R y,const char *s,int align=0)
{ // align: 0 left, 1 center, 2 right
    x -= GraphTextWidth(s)*align/2.;
    //  no bold font in glut => bold = the text is plotted twice with a shift of one pixel
    for (int b=0;b<=graphbold;++b)
    {
        //  +0.5 : center of the pixel, otherwise the 2 positions can be rounded to the same pixel
        glRasterPos2d(floor(x)+b+0.5,floor(y)+0.5);
        for (const char *c=s;*c;++c)
            glutBitmapCharacter(GraphFont(),*c);
    }
}

static void GraphTics(R vmin,R vmax,int nbmax,vector<R> & tics,int & prec)
{ // the tics k*step in [vmin,vmax] with step = 1, 2 or 5 * 10^n, prec: nb of digit to print the tics
    tics.clear();
    prec=6;
    if( !(vmax>vmin) ) return;
    R d=(vmax-vmin)/max(nbmax,1);
    R p=pow(10.,floor(log10(d)));
    R r=d/p;
    R step = p*( r>5 ? 10 : (r>2 ? 5 : (r>1 ? 2 : 1)));
    R k0=ceil(vmin/step-1e-9), k1=floor(vmax/step+1e-9);
    if( !(k1-k0 < 100) ) return;
    for (R k=k0;k<=k1;++k)
        tics.push_back(k*step+0.);// +0. to remove -0
    R vm=max(fabs(vmin),fabs(vmax));
    prec = min(15,max(2,(int) ceil(log10(vm/step))+1));
}

inline R GLog(R v,bool lg) { return lg ? log10(v) : v;}

static void GraphAxisTics(R vmin,R vmax,int nbmax,bool lg,vector<R> & tics,vector<string> & labels,vector<R> & minors,int & prec)
{ //  the tics of an axis with the labels, log scale: vmin, vmax are the log10 of the values
    char buf[64];
    tics.clear();
    labels.clear();
    minors.clear();
    nbmax=max(nbmax,1);
    if(lg && floor(vmax+1e-9)-ceil(vmin-1e-9) >= 1)
    { // a tic by decade (or by step decades), and minor tics 2,3,..,9 if step == 1
        int k0=(int) ceil(vmin-1e-9), k1=(int) floor(vmax+1e-9);
        int step=max(1,(k1-k0+nbmax)/nbmax);
        for(int k=k0;k<=k1;k+=step)
        {
            tics.push_back(k);
            snprintf(buf,sizeof(buf),"%g",pow(10.,k));
            labels.push_back(buf);
        }
        if(step==1)
            for(int k=k0-1;k<=k1;++k)
                for(int m=2;m<10;++m)
                {
                    R t=k+log10((R) m);
                    if(t>=vmin && t<=vmax) minors.push_back(t);
                }
        prec=4;
        return;
    }
    GraphTics(vmin,vmax,nbmax,tics,prec);
    for(size_t i=0;i<tics.size();++i)
    {
        if(lg) snprintf(buf,sizeof(buf),"%.3g",pow(10.,tics[i]));
        else snprintf(buf,sizeof(buf),"%.*g",prec,tics[i]);
        labels.push_back(buf);
    }
    if(lg) prec=4;
}

static void GraphTextRotated(const char *s,vector<unsigned char> & pix,int & w,int & h)
{ //  the bitmap fonts can not be rotated => the text is plotted on the bottom left corner of the window
  //  (must be empty), the pixels are read, rotated of 90 degree in pix (size h x w) and the corner is cleaned
    int fh=GraphFontHeight(), dsc=fh/3+2;
    w=GraphTextWidth(s)+3;
    h=fh+dsc+3;
    vector<unsigned char> a(3*w*h);
    pix.resize(3*w*h);
    glColor3d(0.,0.,0.);
    GraphText(1,dsc,s);
    glPixelStorei(GL_PACK_ALIGNMENT,1);
    glReadBuffer(GL_BACK);
    glReadPixels(0,0,w,h,GL_RGB,GL_UNSIGNED_BYTE,&a[0]);
    glColor3d(1.,1.,1.);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glBegin(GL_QUADS);
    glVertex2d(0,0);
    glVertex2d(w,0);
    glVertex2d(w,h);
    glVertex2d(0,h);
    glEnd();
    for(int j=0;j<h;++j)
        for(int i=0;i<w;++i)
            for(int c=0;c<3;++c)
                pix[3*(i*h+(h-1-j))+c]=a[3*(j*w+i)+c];
}

OnePlotCurve::OnePlotCurve(PlotStream & f,int nfield,ThePlot *theplot)
:OnePlot(3,2,1),cas(0)
{
    f >> xx>>yy;
    if( nfield >=3)
        f >> zz;
    if( nfield ==4)
        f >> cc;
    ffassert(f.good());
    ffassert(xx.N() && yy.N() && xx.N() == yy.N());
    Pmin=Minc(Pmin,R3(xx.min(),yy.min(),Pmin.z));
    Pmax=Maxc(Pmax,R3(xx.max(),yy.max(),Pmax.z));
    if(zz.N())
    {
        Pmin.z = Min(Pmin.z,zz.min());
        Pmax.z = Max(Pmax.z,zz.max());
        if( cc.N() ==0)
        {
            fmax = Max(fmax,zz.max());
            fmin = Min(fmin,zz.min());
            theplot->withiso=true;
        }

    }
    else
        Pmin.z=Pmax.z=0.; // Bof Bof .. FH Juin 2019 ...
    if( cc.N() )
    {
        fmax = Max(fmax,cc.max());
        fmin = Min(fmin,cc.min());
        theplot->withiso=true;
    }
    if(debug>3) cout << " OnePlotCurve nbfield "<< nfield << " max, N= " << xx.max() << " "
        << xx.N() << ", " << yy.max() << " " << yy.N() << " "<< zz.N() << " " << cc.N() << " f " << fmin << " " << fmax << endl;;

}

bool OnePlotCurve::Graph1d(R2 & Pmn,R2 & Pmx,bool lx,bool ly) const
{
    if(zz.N() || cc.N()) return false;
    for (int i=0;i<xx.N();i++)
        if( (!lx || xx[i]>0) && (!ly || yy[i]>0) )
        {
            R2 P(GLog(xx[i],lx),GLog(yy[i],ly));
            Pmn=Minc(Pmn,P);
            Pmx=Maxc(Pmx,P);
        }
    return true;
}

bool OnePlotFEL::Graph1d(R2 & Pmn,R2 & Pmx,bool lx,bool ly) const
{ // scalar function on a mesh along the x axis => graph (x,f(x))
    if( what != 14 && what != 20) return false;
    if( Th->Pmin.y != Th->Pmax.y || Th->Pmin.z != Th->Pmax.z || v.N()==0) return false;
    int nsubV=Psub.N();
    int nK=v.N()/ Th->nt;
    for(int k=0,o=0;k<Th->nt;++k, o+= nK)
    {
        const MeshL::Element & K=(*Th)[k];
        for(int i=0;i<nsubV;++i)
        {
            R x=K(Psub[i]).x, y=v[o+i];
            if( (!lx || x>0) && (!ly || y>0) )
            {
                R2 P(GLog(x,lx),GLog(y,ly));
                Pmn=Minc(Pmn,P);
                Pmx=Maxc(Pmx,P);
            }
        }
    }
    return true;
}

void OnePlotFEL::DrawGraph(OneWindow *win)
{ // 1d graph: the function on each sub segment of each element (so the jumps are not plotted)
    ThePlot & plot=*win->theplot;
    const MeshL & Th(*this->Th);
    double z = plot.z0;
    int nsubT= Ksub.N()/2;
    int nsubV=Psub.N();
    int nK=v.N()/ Th.nt;
    ffassert(v.N()== Th.nt*nK);
    int s=GraphScale();
    glLineWidth(1.5*s);
    glPointSize(5*s);
    GraphColor(num,plot.grey);
    for (int pass=0;pass<2;++pass)
    {
        if(pass==0 && win->graphstyle%3==2) continue;
        if(pass==1 && win->graphstyle%3==0) continue;
        bool lx=win->graphlogx, ly=win->graphlogy;
        glBegin(pass ? GL_POINTS : GL_LINES);
        for(int k=0,o=0;k<Th.nt;++k, o+= nK)
        {
            const MeshL::Element & K=Th[k];
            if(pass)
                for(int i=0;i<nsubV;++i)
                {
                    R x=K(Psub[i]).x, y=v[o+i];
                    if( (lx && x<=0) || (ly && y<=0) ) continue;// log scale
                    glVertex3d(GLog(x,lx),GLog(y,ly),z);
                }
            else
                for(int sk=0;sk<nsubT;++sk)
                {
                    int i0= Ksub[sk*2+0], i1=Ksub[sk*2+1];
                    R x0=K(Psub[i0]).x, y0=v[o+i0], x1=K(Psub[i1]).x, y1=v[o+i1];
                    if( (lx && (x0<=0 || x1<=0)) || (ly && (y0<=0 || y1<=0)) ) continue;// log scale
                    glVertex3d(GLog(x0,lx),GLog(y0,ly),z);
                    glVertex3d(GLog(x1,lx),GLog(y1,ly),z);
                }
        }
        glEnd();
    }
    glLineWidth(kscreenscale);
    glPointSize(1);
}

void OnePlotCurve::dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn2,R & vmx2) const
{
    if( cc.N() == xx.N())
    {
        ; // afaire ????
    }
}
void OnePlotCurve::Draw(OneWindow *win)
{
    initlist();
    ThePlot & plot= *win->theplot;
    double z = plot.z0;
    if(debug>3) cout << " OnePlotCurve::Draw " << (cc.N() != xx.N()) << endl;
    if(win->Graph1d())
    { // 1d graph: style 0 lines, 1 lines and points, 2 points
        int s=GraphScale();
        glLineWidth(1.5*s);
        glPointSize(5*s);
        GraphColor(num,plot.grey);
        for (int pass=0;pass<2;++pass)
        {
            if(pass==0 && win->graphstyle%3==2) continue;
            if(pass==1 && win->graphstyle%3==0) continue;
            //  log scale: only the positive values, so the line can be cut
            bool lx=win->graphlogx, ly=win->graphlogy;
            bool in=false;
            if(pass) glBegin(GL_POINTS);
            for (int i=0;i<xx.N();i++)
            {
                bool ok = (!lx || xx[i]>0) && (!ly || yy[i]>0);
                if(!pass && ok && !in) glBegin(GL_LINE_STRIP), in=true;
                if(!pass && !ok && in) glEnd(), in=false;
                if(ok) glVertex3d(GLog(xx[i],lx),GLog(yy[i],ly),z);
            }
            if(pass || in) glEnd();
        }
        glLineWidth(kscreenscale);
        glPointSize(1);
        return;
    }
    if(cas%2)
        glBegin(GL_POINTS);
    else
        glBegin(GL_LINE_STRIP);
    plot.color(2);
    if( (zz.N()!=xx.N())  && (cc.N()==0)  )
        for (int i=0;i<xx.N();i++)
            glVertex3d(xx[i],yy[i],z);
    else if (cc.N()==0)
        for (int i=0;i<xx.N();i++)
            glVertex3d(xx[i],yy[i],zz[i]);
    else
        for (int i=0;i<xx.N();i++)
        {
            int col = 2+dichotomie(plot.Viso,cc[i]);
            if(col<1) col=1;
            if(col>plot.Viso.N()) col=plot.Viso.N();
            plot.color(col);
            glVertex3d(xx[i],yy[i],zz[i]);
        }

    glEnd();

}

void OnePlotBorder::Draw(OneWindow *win)
{
    initlist();

    glDisable(GL_DEPTH_TEST);
    ThePlot & plot= *win->theplot;
    R h = 10*win->hpixel*kscreenscale;
    double z = plot.z0;
    plot.SetColorTable(16) ;

    if(win->plotdim==3)
    for(int i=0;i<data.size() ;++i) {
        vector<pair<long,R3> > & v=data[i];
        ShowGlerror("end OnePlotBorder::Draw  1");

        for(int j=1;j<v.size();++j) {
            plot.color(2+v[j].first);
            R3 Po(v[j-1].second), Pn(v[j].second);
            R3 uv(Po,Pn);
            double l = Max(sqrt((uv,uv)),1e-20);
            R3 nx(1.,0.,0.), ny(0.,1.,0.), nz(0.,0.,1.);
            R3 dd = uv*(-h/l);
            R3 dnx = (dd^nx)*0.7, dny = (dd^ny)*0.7, dnz = (dd^nz)*0.7;
            
            glLineWidth(2*kscreenscale);
            glBegin(GL_LINES);
            win->Seg3(Po,Pn);
            glEnd();

            glLineWidth(kscreenscale);
            glBegin(GL_LINES);
            if(j!=1) {
              win->Seg3(Po,Po+dd+dnx);
              win->Seg3(Po,Po+dd-dnx);
              win->Seg3(Po,Po+dd+dny);
              win->Seg3(Po,Po+dd-dny);
              win->Seg3(Po,Po+dd+dnz);
              win->Seg3(Po,Po+dd-dnz);
            }
            glEnd();
        }
        ShowGlerror("end OnePlotBorder::Draw  2");

        glPointSize(7);
        glBegin(GL_POINTS);
        int l= v.size()-1;
        plot.color(2+v[0].first);
        glVertex3d(v[0].second.x,v[0].second.y,v[0].second.z);
        plot.color(2+v[l].first);
        glVertex3d(v[l].second.x,v[l].second.y,v[l].second.z);
        glEnd();
        glPointSize(1);
        ShowGlerror("end OnePlotBorder::Draw  3");
    }
    
    else if(win->plotdim==2)
    for(int i=0;i<data.size() ;++i) {
        vector<pair<long,R3> > & v=data[i];
        ShowGlerror("end OnePlotBorder::Draw  1");

        for(int j=1;j<v.size();++j) {
            plot.color(2+v[j].first);
            R2 Po(v[j-1].second.p2()), Pn(v[j].second.p2());
            R2 uv(Po,Pn);
            double l = Max(sqrt((uv,uv)),1e-20);

            R2 dd = uv*(-h/l);
            R2 dn = dd.perp()*0.5;
            glLineWidth(2*kscreenscale);
            glBegin(GL_LINES);
            win->Seg(Po,Pn);
            glEnd();

            glLineWidth(kscreenscale);
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


inline void plotquadx(float y, float iy, float z, float iz, float x1, bool outline = 0){
  glVertex3d(x1, y, z);
  glVertex3d(x1, y+iy, z);
  if (outline) glVertex3d(x1, y+iy, z);
  glVertex3d(x1, y+iy, z+iz);
  glVertex3d(x1, y+iy, z+iz);
  glVertex3d(x1, y, z+iz);
  if (outline) glVertex3d(x1, y, z+iz);
  glVertex3d(x1, y, z);
}

inline void plotquady(float x, float ix, float z, float iz, float y1, bool outline = 0){
  glVertex3d(x, y1, z);
  glVertex3d(x+ix, y1, z);
  if (outline) glVertex3d(x+ix, y1, z);
  glVertex3d(x+ix, y1, z+iz);
  glVertex3d(x+ix, y1, z+iz);
  glVertex3d(x, y1, z+iz);
  if (outline) glVertex3d(x, y1, z+iz);
  glVertex3d(x, y1, z);
}

inline void plotquadz(float x, float ix, float y, float iy, float z1, bool outline = 0){
  glVertex3d(x, y, z1);
  glVertex3d(x+ix, y, z1);
  if (outline) glVertex3d(x+ix, y, z1);
  glVertex3d(x+ix, y+iy, z1);
  glVertex3d(x+ix, y+iy, z1);
  glVertex3d(x, y+iy, z1);
  if (outline) glVertex3d(x, y+iy, z1);
  glVertex3d(x, y, z1);
}

void OnePlotHMatrix::Draw(OneWindow *win)
{
  if (si <= 0 || sj <= 0)
    return;

  R scalez, zeps, dx, dy, dz;

  if (win->plotdim==2) {
    glDisable(GL_DEPTH_TEST);
    scalez = 0;
    zeps = 0;
    dx = 0;
    dy = 0;
    dz = 0;
  }
  else {
    glEnable(GL_DEPTH_TEST);
    scalez = 0.5;
    zeps = 0.001;
    dx = -5e-4*(win->Pvue3.x - win->cam.x);
    dy = -5e-4*(win->Pvue3.y - win->cam.y);
    dz = -5e-4*(win->Pvue3.z - win->cam.z);
  }

  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

  R z1= 0;

  glBegin(GL_TRIANGLES);
  for (int i=0;i < nbdense; i++) {
    std::pair<int,int>& offset = offsetsdense[i];
    std::pair<int,int>& size = sizesdense[i];
    z1 = scalez;
    if (size.first > 1 || size.second > 1) {
      glColor3f(1,0,0); //glColor3f((20.)/255,(80.)/255,(20.)/255);
      plotquadz((float)offset.second/sj, (float)size.second/sj, 1.-(float)offset.first/si, -(float)size.first/si, z1);

      if (win->plotdim==3) {
        plotquadx(1.-(float)offset.first/si, -(float)size.first/si, 0, z1, (float)offset.second/sj);
        plotquadx(1.-(float)offset.first/si, -(float)size.first/si, 0, z1, (float)(offset.second+size.second)/sj);
        plotquady((float)offset.second/sj, (float)size.second/sj, 0, z1, 1.-(float)offset.first/si);
        plotquady((float)offset.second/sj, (float)size.second/sj, 0, z1, 1.-(float)(offset.first+size.first)/si);
        plotquadz((float)offset.second/sj, (float)size.second/sj, 1.-(float)offset.first/si, -(float)size.first/si, 0);
      }
    }
  }

  for (int i=0;i < nblr; i++) {
    std::pair<int,int>& offset = offsetslr[i];
    std::pair<int,int>& size = sizeslr[i];
    z1 = scalez*(1-compression[i]);
    if (size.first > 1 || size.second > 1) {
      glColor3f((20+compression[i]*80)/255,(80+compression[i]*170)/255,(20+compression[i]*80)/255);
      if (win->theplot->fill || win->plotdim==3)
        plotquadz((float)offset.second/sj, (float)size.second/sj, 1.-(float)offset.first/si, -(float)size.first/si, z1);
      else {
        plotquadz((float)offset.second/sj, (float)rankslr[i]/sj, 1.-(float)offset.first/si, -(float)size.first/si, z1);
        plotquadz((float)offset.second/sj, (float)size.second/sj, 1.-(float)offset.first/si, -(float)rankslr[i]/si, z1);
      }

      if (win->plotdim==3) {
        plotquadx(1.-(float)offset.first/si, -(float)size.first/si, 0, z1, (float)offset.second/sj);
        plotquadx(1.-(float)offset.first/si, -(float)size.first/si, 0, z1, (float)(offset.second+size.second)/sj);
        plotquady((float)offset.second/sj, (float)size.second/sj, 0, z1, 1.-(float)offset.first/si);
        plotquady((float)offset.second/sj, (float)size.second/sj, 0, z1, 1.-(float)(offset.first+size.first)/si);
        plotquadz((float)offset.second/sj, (float)size.second/sj, 1.-(float)offset.first/si, -(float)size.first/si, 0);
      }
    }
  }
  glEnd();

  // For sparse matrices
  double vmax = -1e+30;
  for (int i=0;i<nblr;i++)
  if (abs(compression[i]) > vmax && abs(compression[i]) < 1e10)
    vmax = abs(compression[i]);
  glPointSize(2);
  glBegin(GL_POINTS);
  for (int i=0;i < nblr; i++) {
    std::pair<int,int>& offset = offsetslr[i];
    std::pair<int,int>& size = sizeslr[i];
    z1 = 0;
    if (size.first == 1 && size.second == 1) {
      float val = compression[i]/vmax;
      glColor3f((20+val*80)/255,(80+val*170)/255,(20+val*80)/255);
      glVertex3d((float)offset.second/sj, 1.-(float)offset.first/si ,z1);
    }
  }
  glEnd();
  if (win->plotdim==3) {
    glBegin(GL_LINES);
    for (int i=0;i < nblr; i++) {
      std::pair<int,int>& offset = offsetslr[i];
      std::pair<int,int>& size = sizeslr[i];
      z1 = 0;
      if (size.first == 1 && size.second == 1) {
        float val = compression[i]/vmax;
        glColor3f((20+val*80)/255,(80+val*170)/255,(20+val*80)/255);
        glVertex3d((float)offset.second/sj, 1.-(float)offset.first/si ,z1);
        glVertex3d((float)offset.second/sj, 1.-(float)offset.first/si ,z1+(float)val);
      }
    }
    glEnd();
  }

  glLineWidth(kscreenscale);
  glColor3f(0,0,0);
  glBegin(GL_LINES);

  for (int i=0;i < nbdense; i++) {
    std::pair<int,int>& offset = offsetsdense[i];
    std::pair<int,int>& size = sizesdense[i];
    z1 = scalez;
    if (size.first > 1 || size.second > 1) {
      plotquadz((float)offset.second/sj+dx, (float)size.second/sj, 1.-(float)offset.first/si+dy, -(float)size.first/si, dz,1);

      if (win->plotdim==3) {
        plotquadx(1.-(float)offset.first/si+dy, -(float)size.first/si, dz, z1+dz, (float)offset.second/sj+dx,1);
        plotquadx(1.-(float)offset.first/si+dy, -(float)size.first/si, dz, z1+dz, (float)(offset.second+size.second)/sj+dx,1);
        plotquadz((float)offset.second/sj+dx, (float)size.second/sj, 1.-(float)offset.first/si+dy, -(float)size.first/si, z1+dz,1);
      }
    }
  }

  for (int i=0;i < nblr; i++) {
    std::pair<int,int>& offset = offsetslr[i];
    std::pair<int,int>& size = sizeslr[i];
    z1 = scalez*(1-compression[i]);
    if (size.first > 1 || size.second > 1) {
      plotquadz((float)offset.second/sj+dx, (float)size.second/sj, 1.-(float)offset.first/si+dy, -(float)size.first/si, dz,1);
      if (!win->theplot->fill && win->plotdim==2)
        plotquadz((float)(offset.second+rankslr[i])/sj+dx, (float)(size.second-rankslr[i])/sj, 1.-(float)(offset.first+rankslr[i])/si+dy, -(float)(size.first-rankslr[i])/si, dz,1);

      if (win->plotdim==3) {
        plotquadx(1.-(float)offset.first/si+dy, -(float)size.first/si, dz, z1+dz, (float)offset.second/sj+dx,1);
        plotquadx(1.-(float)offset.first/si+dy, -(float)size.first/si, dz, z1+dz, (float)(offset.second+size.second)/sj+dx,1);
        plotquadz((float)offset.second/sj+dx, (float)size.second/sj, 1.-(float)offset.first/si+dy, -(float)size.first/si, z1+dz,1);
      }
    }
  }

  glEnd();
  glColor3f(0,0,0);
  for (int i=0;i < nblr; i++) {
    std::pair<int,int>& offset = offsetslr[i];
    std::pair<int,int>& size = sizeslr[i];
    z1 = scalez*(1-compression[i]);
    if (size.first > 1 || size.second > 1) {
      string s = win->theplot->value ? std::to_string(int(100*compression[i])) : std::to_string(rankslr[i]);
      float scale = 0.005*std::min((float)size.first/si,(float)size.second/sj)/**std::min(mSize.x(),mSize.y())*/;
      if (win->theplot->value)
        scale *= 0.9;

      R textx = -36*scale*s.length()+(float)(offset.second+size.second*0.5)/sj;
      R texty = -48*scale+1.-float(offset.first+size.first*0.5)/si;
      R textz = z1+zeps;

      if (win->theplot->value)
        textx -= 10*scale*s.length();

      bool ptext = 0;
      if (win->plotdim==2 && scale > 0.028/std::min(win->height,win->width)*win->hpixel*win->width)
        ptext = 1;
      if (win->plotdim==3)// && sqrt((win->Pvue3.x-textx)*(win->Pvue3.x-textx)+(win->Pvue3.y-texty)*(win->Pvue3.y-texty)+(win->Pvue3.z-textz)*(win->Pvue3.z-textz)) < 0.1)
        ptext = 1;

      if (ptext) {
        glPushMatrix();
        glTranslatef(textx, texty, textz);
        glScalef(scale,scale,scale);
        for(char& c : s)
          glutStrokeCharacter(GLUT_STROKE_ROMAN,c);
        if (win->theplot->value) {
          glScalef(0.3,0.3,0.3);
          glutStrokeCharacter(GLUT_STROKE_ROMAN,'%');
        }
        glPopMatrix();
      }
    }
  }

}


OneWindow::OneWindow(int h,int w,ThePlot *p)
:
countdisplay(0),
icurrentPlot(lplots.begin()),
lplotssize(0),
height(h),width(w),theplot(0),hpixel(1),
Bmin(0,0),Bmax(1,1),oBmin(Bmin),oBmax(Bmax),zmin(0),zmax(1),
windowdump(false),windowdumppng(false),help(false),helpall(false),showhx(0), rapz0(-1.),rapz(1),withlight(false),
changearrow(true),changeiso(true), keepPV(false),init(false),
graphgrid(false),graphkey(true),graphstyle(0),graphlogx(false),graphlogy(false),pretty(true),probe(true),isolines(true),axes(true),prettybar(true),cmapkind(0),cmaptex(0),cmaptexkind(-1),userlogx(-1),userlogy(-1),explode(0),cutaxis(0),cutpos(0.5),cutN(1.,0.,0.),cutedit(false),cutfreeze(false),threshold(0.5),isosurf(false),xmouse(-1),ymouse(-1)
{

    add(p);
}


void OneWindow::set(ThePlot *p)
{
    theplot=p;
    if(p)
    {
        if( debug) cout << "OneWindow:: set "<< p->keepPV << " -> "<< keepPV << endl;
        plotdim=p->plotdim;
        if(!keepPV)keepPV=p->keepPV;
        pNormalT=p->pNormalT;
        rapz0 = p->ZScale;
        // log scale of the 1d graph asked in the comment of the plot
        //  or with the keys x, y
        if(p->graph1d)
        {
            graphlogx = userlogx>=0 ? userlogx : p->glogx==1;
            graphlogy = userlogy>=0 ? userlogy : p->glogy==1;
        }
    }
    if(!init)
    {
        rapz0 =-1;
        if(p) rapz0 = p->ZScale;
        DefaultView(2) ;
    }
    else if(Graph1d() && !keepPV)
        DefaultView(1);// 1d graph: the axis of an other plot are meaningless

}

void OneWindow::add(ThePlot *p)
{
    if(p) {
        if( debug) cout << "OneWindow:: add "<< p->keepPV << " -> "<< keepPV << endl;
        if(!keepPV) keepPV=p->keepPV;
        lplots.push_back(p);
        lplotssize++;
        ++icurrentPlot;
        if(icurrentPlot==lplots.end())
            --icurrentPlot;// the previous
        if(icurrentPlot != lplots.end())
            set(*icurrentPlot);
        if( lplotssize>20)// pass 10 -> 20 for O. Pironneau 21/12/2020 FH.
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

void OneWindow::DefaultView(int state)
{
    if(debug>1)  cout << "DefaultView " << state << " " <<keepPV << " theplot " << theplot << endl;
    if(keepPV)
    {
        if(state==0 && init) return;
    }
    else {
        if( theplot && theplot->ZScale>0) rapz0=theplot->ZScale;
        else rapz0=-1;
    }

    if(theplot)
    {
        init =1;
        plotdim=theplot->plotdim;
        keepPV=theplot->keepPV;
        R3 A(theplot->Pmin),B(theplot->Pmax);
        R3 D(A,B);
        R dxy= max(D.x,D.y);
        zmax = theplot->fmax;
        zmin = theplot->fmin;
        theta=theplot->theta;
        phi=theplot->phi;
        if(theplot->ZScale>0) rapz0=theplot->ZScale;
        if(theplot->datadim==3) rapz0=1;
        else   if(rapz0<=0)
        {
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



        if(Graph1d())
        { // 1d graph: the bound of the axis
            R2 Gmn(dinfty,dinfty),Gmx(-dinfty,-dinfty);
            theplot->GraphBound(Gmn,Gmx,graphlogx,graphlogy);
            if(Gmn.x>Gmx.x) Gmn.x=0.,Gmx.x=1.;// no value (log scale)
            if(Gmn.y>Gmx.y) Gmn.y=0.,Gmx.y=1.;
            A=R3(Gmn.x,Gmn.y,0.);
            B=R3(Gmx.x,Gmx.y,0.);
            D=R3(A,B);
        }
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

        if(Graph1d())
        {
            if(theplot->boundingbox.size() ==4)
            { // the bounding box is given without the log
                if(graphlogx && A.x>0) A.x=log10(A.x),B.x=log10(B.x);
                if(graphlogy && A.y>0) A.y=log10(A.y),B.y=log10(B.y);
            }
            // no empty axis range
            if( !(B.x>A.x) ) { R d = A.x ? fabs(A.x)*0.05 : 0.5; A.x -= d; B.x += d;}
            if( !(B.y>A.y) ) { R d = A.y ? fabs(A.y)*0.05 : 0.5; A.y -= d; B.y += d;}
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
    glViewport(0, 0,width, height);
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
        gluPerspective(focal*180./my_pi,aspect,znear,zfare);
        if(debug>2)
        {
            cout << " BB " << Bmin3 << " , " << Bmax3 << " rapz " << rapz << " dmax: "<< dmax << endl;
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
        viewdim=3;

    }
    else
    {
        ShowGlerror("Begin SetView");
        glDisable(GL_DEPTH_TEST);
        if(Graph1d())
        { // the curves are plotted (and clipped) in the box of the axis
            int x0,y0,x1,y1;
            GraphBox(x0,y0,x1,y1);
            glViewport(x0,y0,x1-x0,y1-y0);
        }
        else
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
            viewdim=2;
        }

        if((debug>3 )) cout << "\t\t\t   SetView " << this << " " << Bmin  << " "
            << Bmax  << " dz  " << dz0 << " " << dz1
            << " theta " << theta << "  phi = "<<  phi << " // "<< Bmax - Bmin << " w " <<width << " h " <<  height<< endl;
        ShowGlerror("0 Set MV");
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        ShowGlerror(" Set MV");
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        ShowGlerror(" Set PM 1");
        R2 B0=Bmin, B1=Bmax;
        bool pi=false,pa=false;
        if(!Graph1d() && PrettyOn() && prettybar && theplot) theplot->PrettyPlots(this,pi,pa);
        //  margins in pixel: on the right for the color bar (pretty rendering),
        //  on the left and at the bottom for the axes (key j)
        int s=GraphScale();
        R mr = (pi || pa) ? min(120.*s,0.4*width) : 0.;
        R ml = 0., mb = 0.;
        if(axes && !Graph1d())
        {
            ml = min(70.*s + (pi && pa ? 120.*s : 0.),0.3*width);// the bar of the arrows is on the left
            mb = min(50.*s,0.3*height);
        }
        if(mr+ml+mb>0)
        { //  the box Bmin, Bmax in the window without the margins (same scale in x and y if aspectratio)
            R wi=width-ml-mr, hi=height-mb;
            R sx=(B1.x-B0.x)/wi, sy=(B1.y-B0.y)/hi;
            if(theplot && theplot->aspectratio) sx=sy=max(sx,sy);
            R2 M=(B0+B1)/2.;
            B0.x = M.x - sx*(ml+wi/2);
            B1.x = B0.x + sx*width;
            B0.y = M.y - sy*(mb+hi/2);
            B1.y = B0.y + sy*height;
        }
        glOrtho(B0.x,B1.x,B0.y,B1.y,dz0,dz1);

        ShowGlerror(" Set PM 2");

        R2 M=(Bmin+Bmax)/2.;
        glTranslated(0,0,-zm);

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
{
    width=w;
    height=h;
    if (theplot && theplot->aspectratio)
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
    R2 O(xx,yy);
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
    if(Graph1d()) // the scale in x and y are not the same
        D2 = R2(-dx*0.05*(Bmax.x-Bmin.x),-dy*0.05*(Bmax.y-Bmin.y));
    oBmin += D2;
    oBmax += D2;
    Bmin += D2;
    Bmax += D2;

}

void OneWindow::cadre(R2 A,R2 B)
{
    changearrow=true;// the size of the arrows depends on the size of the pixel

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
    if(!theplot) return;
    SetScreenView() ;
    glColor3d(0.,0.,0.);

    if(help)
    {
        theplot->DrawHelp(this);
       // help=false;
    }
    else
    {
        ShowGlerror("Begin Display");

        if(theplot && PrettyOn()) DrawBackground();
        if(theplot)
            theplot->Draw(this);
        DrawProbe();// after the scene: the depth under the mouse
        ShowGlerror("After Display");
    }
}
void OneWindow::cadreortho(R2 A, R2 B)
{
    changearrow=true;// the size of the arrows depends on the size of the pixel
    R2 D(A,B);
    oBmin=A;
    oBmax=B;

    if ( D.y*width < D.x*height)
        // width -> infty => D.x la ref
        D.y = D.x*(double) height/ width;
    else // height -> infty => D.y la ref
        D.x = D.y*(double) width/height;
    R2 M=(A+B)/2., D2=D/2.;

    Bmin= M - D2;
    Bmax= M + D2;
    hpixel = (Bmax.x-Bmin.x)/width;
    if((debug > 10)) cout << " cadreortho: "<< " :: " << Bmin << " " << Bmax <<" oB " << oBmin << " " << oBmax << " w: " <<width
                        << " h: " << height << endl;

}
bool OneWindow::SetCutPlane(const R3 & Pmn,const R3 & Pmx,R3 & N,R & c)
{ //  the cut plane (N,P) = c of a mesh of bounding box Pmn, Pmx, and set the clip plane of OpenGL:
  //  we keep the part behind the plane seen from the camera
    if(!cutaxis) return false;
    N=cutN;
    //  min and max of (N,P) on the bounding box of the mesh
    R cmin=min(N.x*Pmn.x,N.x*Pmx.x)+min(N.y*Pmn.y,N.y*Pmx.y)+min(N.z*Pmn.z,N.z*Pmx.z);
    R cmax=max(N.x*Pmn.x,N.x*Pmx.x)+max(N.y*Pmn.y,N.y*Pmx.y)+max(N.z*Pmn.z,N.z*Pmx.z);
    c = cmin + cutpos*(cmax-cmin);
    R3 CP(cam.x-Pvue3.x,cam.y-Pvue3.y,cam.z-Pvue3.z*rapz);
    R sg = (N,CP) > 0 ? -1 : 1;
    GLdouble eq[4]={sg*N.x,sg*N.y,sg*N.z,-sg*c};
    for(int i=0;i<4;++i) cuteq[i]=eq[i];
    glClipPlane(GL_CLIP_PLANE0,eq);
    return true;
}

void OneWindow::zoomat(int w,int h,R coef)
{ //  zoom, the point under the cursor does not move (wheel of the mouse)
    coef_dist*=coef;
    if(plotdim==3) return;
    GLdouble x=w,y=height-h,z=(zmin+zmax)/2.;
    GLdouble xx,yy,zz;
    gluUnProject( x,y,z,modelMatrix,projMatrix,viewport,&xx,&yy,&zz);
    R2 O(xx,yy);
    R2 A = O - R2(oBmin,O)*coef;
    R2 B = O + R2(O,oBmax)*coef;
    if (theplot && theplot->aspectratio)
        cadreortho(A,B);
    else
        cadre(A,B);
}

void OneWindow::RotateCut(R ah,R av)
{ //  rotation of the normal of the cut plane: angle ah around the vertical axis (z),
  //  and angle av around the horizontal axis of the screen
    if(!cutaxis) return;
    R3 A[2]={R3(0.,0.,1.),R3(-sin(theta),cos(theta),0.)};
    R a[2]={ah,av};
    for(int i=0;i<2;++i)
    { // Rodrigues formula
        R ca=cos(a[i]), sa=sin(a[i]);
        cutN = cutN*ca + (A[i]^cutN)*sa + A[i]*((A[i],cutN)*(1.-ca));
    }
    cutN /= cutN.norme();
}

void OneWindow::DrawProbe()
{ //  the value of the function at the point under the mouse, on the bottom left corner (key P)
    if(!probe || !theplot || Graph1d() || help) return;
    if(xmouse<0 || ymouse<0 || xmouse>=width || ymouse>=height) return;
    GLdouble wx,wy,wz;
    int xp=xmouse, yp=height-1-ymouse;
    R3 P;
    if(plotdim==2)
    {
        gluUnProject(xp,yp,0.5,modelMatrix,projMatrix,viewport,&wx,&wy,&wz);
        P=R3(wx,wy,0.);
    }
    R val=0;
    R3 vec;
    int nc=0;
    bool oniso=false;
    if(plotdim==3)
    { // the point on the visible surface: the depth of the pixel
        GLfloat d=1;
        glReadPixels(xp,yp,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&d);
        //  the translucent isosurfaces are not in the depth buffer: a ray from the eye to this point
        GLdouble ax,ay,az,bx,by,bz;
        gluUnProject(xp,yp,0.,modelMatrix,projMatrix,viewport,&ax,&ay,&az);
        gluUnProject(xp,yp,d<1.f ? d : 1.,modelMatrix,projMatrix,viewport,&bx,&by,&bz);
        for (list<OnePlot *>::iterator i= theplot->plots.begin();i != theplot->plots.end() && !oniso; ++i)
        {
            OnePlotFE3 *f3=dynamic_cast<OnePlotFE3 *>(*i);
            if(f3 && f3->showisovals && f3->ProbeRay(R3(ax,ay,az),R3(bx,by,bz),this,P,val))
                oniso=true, nc=1;
        }
        if(!oniso)
        {
            if(!(d<1.f)) return;
            gluUnProject(xp,yp,d,modelMatrix,projMatrix,viewport,&wx,&wy,&wz);
            P=R3(wx,wy,wz);
        }
    }
    for (list<OnePlot *>::iterator i= theplot->plots.begin();i != theplot->plots.end() && !nc; ++i)
        if(!(*i)->Probe(P,val,vec,nc)) nc=0;
    char buf[256];
    bool d3 = theplot->datadim==3;
    int l = d3 ? snprintf(buf,sizeof(buf),"x = %.4g  y = %.4g  z = %.4g",P.x,P.y,P.z)
               : snprintf(buf,sizeof(buf),"x = %.4g  y = %.4g",P.x,P.y);
    if(nc==1)
        snprintf(buf+l,sizeof(buf)-l,oniso ? "   :  %.6g   (isosurface)" : "   :  %.6g",val);
    else if(nc==2)
        snprintf(buf+l,sizeof(buf)-l,"   :  (%.4g, %.4g)   |.| = %.4g",vec.x,vec.y,val);
    else if(nc==3)
        snprintf(buf+l,sizeof(buf)-l,"   :  (%.4g, %.4g, %.4g)   |.| = %.4g",vec.x,vec.y,vec.z,val);
    int s=GraphScale();
    SetScreenView();
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_1D);
    glColor3d(0.,0.,0.);
    GraphText(6*s,6*s,buf);
    SetView();
}

GLuint OneWindow::ColormapTexture()
{ //  the texture of the colormap (one by window: the OpenGL contexts are not shared)
    int kind = (theplot && theplot->grey) ? 4 : cmapkind;
    if(cmaptex && kind==cmaptexkind) return cmaptex;
    if(!cmaptex) glGenTextures(1,&cmaptex);
    unsigned char data[256*3];
    for(int i=0;i<256;++i)
    {
        float c[3];
        ColormapRGB(kind,i/255.,c);
        for(int j=0;j<3;++j) data[3*i+j]=(unsigned char) (255*c[j]+0.5);
    }
    glBindTexture(GL_TEXTURE_1D,cmaptex);
    glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glTexImage1D(GL_TEXTURE_1D,0,GL_RGB,256,0,GL_RGB,GL_UNSIGNED_BYTE,data);
    cmaptexkind=kind;
    return cmaptex;
}

void OneWindow::DrawBackground()
{ //  pretty rendering: gradient background (light, so the texts in black stay readable)
    SetScreenView();
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_1D);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glBegin(GL_QUADS);
    glColor3f(0.80f,0.84f,0.91f);
    glVertex2d(0,0);
    glVertex2d(width,0);
    glColor3f(0.99f,0.99f,1.f);
    glVertex2d(width,height);
    glVertex2d(0,height);
    glEnd();
}

void OneWindow::DrawColorBar(R fmn,R fmx,R x0,const char *title,const vector<R> *marks)
{ //  pretty rendering: vertical color bar with the values, at the abscissa x0 (pixel)
    if(!theplot || !(fmx>fmn)) return;
    int s=GraphScale(), fh=GraphFontHeight();
    ColormapTexture();
    SetScreenView();
    R x1=x0+18*s, y0=0.22*height, y1=0.78*height;
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D,cmaptex);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glBegin(GL_QUADS);
    glTexCoord1d(0.); glVertex2d(x0,y0);
    glTexCoord1d(0.); glVertex2d(x1,y0);
    glTexCoord1d(1.); glVertex2d(x1,y1);
    glTexCoord1d(1.); glVertex2d(x0,y1);
    glEnd();
    glDisable(GL_TEXTURE_1D);
    glColor3d(0.,0.,0.);
    glLineWidth(s);
    glBegin(GL_LINE_LOOP);
    glVertex2d(x0,y0); glVertex2d(x1,y0); glVertex2d(x1,y1); glVertex2d(x0,y1);
    glEnd();
    if(marks)
    { // the values of the isosurfaces: a line across the bar
        glLineWidth(2*s);
        glBegin(GL_LINES);
        for(size_t i=0;i<marks->size();++i)
        {
            if((*marks)[i]<fmn || (*marks)[i]>fmx) continue;
            R y=y0+((*marks)[i]-fmn)/(fmx-fmn)*(y1-y0);
            glVertex2d(x0,y); glVertex2d(x1,y);
        }
        glEnd();
        glLineWidth(s);
    }
    vector<R> tics;
    int prec;
    GraphTics(fmn,fmx,6,tics,prec);
    char buf[64];
    glBegin(GL_LINES);
    for(size_t i=0;i<tics.size();++i)
    {
        R y=y0+(tics[i]-fmn)/(fmx-fmn)*(y1-y0);
        glVertex2d(x1,y); glVertex2d(x1+5*s,y);
    }
    glEnd();
    for(size_t i=0;i<tics.size();++i)
    {
        R y=y0+(tics[i]-fmn)/(fmx-fmn)*(y1-y0);
        snprintf(buf,sizeof(buf),"%.*g",prec,tics[i]);
        GraphText(x1+8*s,y-fh/2.,buf);
    }
    //  the bound of the data at the ends
    snprintf(buf,sizeof(buf),"max %.4g",fmx);
    GraphText((x0+x1)/2,y1+8*s,buf,1);
    snprintf(buf,sizeof(buf),"min %.4g",fmn);
    GraphText((x0+x1)/2,y0-8*s-fh,buf,1);
    if(title) GraphText((x0+x1)/2,y1+12*s+fh,title,1);
    glLineWidth(kscreenscale);
}

//  the boxes (pixel) of the texts of the axes already drawn: a text over another one is not drawn
static vector<R> axisboxes;

static bool AxisText(R x,R y,const char *t,int align)
{ //  GraphText if the text does not overlap the texts already drawn
    int s=GraphScale();
    R w=GraphTextWidth(t), h=GraphFontHeight();
    R x0=x-w*align/2.-2*s, x1=x0+w+4*s, y0=y-2*s, y1=y+h+2*s;
    for(size_t i=0;i<axisboxes.size();i+=4)
        if(x0<axisboxes[i+2] && axisboxes[i]<x1 && y0<axisboxes[i+3] && axisboxes[i+1]<y1) return false;
    R b[4]={x0,y0,x1,y1};
    axisboxes.insert(axisboxes.end(),b,b+4);
    GraphText(x,y,t,align);
    return true;
}

void OneWindow::DrawAxis(const R3 & A,const R3 & B,R vmin,R vmax,const R2 & C,R gap,const char *name)
{ //  an axis from A to B (world) for the values vmin to vmax, drawn in pixel on the screen:
  //  tics and values on the side opposite to C (the center of the plot on the screen), shifted of gap,
  //  the name of the axis at the middle, after the values
    if(!(vmax>vmin)) return;
    int s=GraphScale(), fh=GraphFontHeight();
    GLdouble x,y,z;
    gluProject(A.x,A.y,A.z,modelMatrix,projMatrix,viewport,&x,&y,&z);
    R2 SA(x,y);
    gluProject(B.x,B.y,B.z,modelMatrix,projMatrix,viewport,&x,&y,&z);
    R2 SB(x,y);
    R2 d(SA,SB);
    R l=d.norme();
    if(l<10*s) return;// the axis is seen end-on
    d/=l;
    R2 o(d.y,-d.x);// outward normal on the screen
    R2 M=(SA+SB)/2.;
    if((o,R2(C,M))<0) o=-o;
    SA += o*gap; SB += o*gap;
    glBegin(GL_LINES);
    glVertex2d(SA.x,SA.y); glVertex2d(SB.x,SB.y);
    glEnd();
    vector<R> tics;
    int prec;
    GraphTics(vmin,vmax,max(2,(int)(l/(70*s))),tics,prec);
    int align = o.x>0.5 ? 0 : (o.x<-0.5 ? 2 : 1);
    char buf[64];
    R ext=0;// the size of the values in the direction o
    for(size_t i=0;i<tics.size();++i)
    {
        R2 P = SA + (SB-SA)*((tics[i]-vmin)/(vmax-vmin));
        if(P.x<0 || P.y<0 || P.x>width || P.y>height) continue;
        R2 Q = P + o*(5*s);
        glBegin(GL_LINES);
        glVertex2d(P.x,P.y); glVertex2d(Q.x,Q.y);
        glEnd();
        snprintf(buf,sizeof(buf),"%.*g",prec,tics[i]);
        Q = P + o*(9*s);
        if(AxisText(Q.x,Q.y-fh/2.+o.y*fh*0.7,buf,align))
            ext=max(ext,fabs(o.x)*GraphTextWidth(buf)+fabs(o.y)*1.2*fh);
    }
    //  the values at the ends of the axis (ex. the max of the function on the axis z),
    //  if no tic is near
    R step = tics.size()>1 ? tics[1]-tics[0] : vmax-vmin;
    R ends[2]={vmin,vmax};
    for(int e=0;e<2;++e)
    {
        bool near=false;
        for(size_t i=0;i<tics.size();++i)
            if(fabs(tics[i]-ends[e])<0.25*step) near=true;
        if(near) continue;
        R2 P = e ? SB : SA;
        if(P.x<0 || P.y<0 || P.x>width || P.y>height) continue;
        R2 Q = P + o*(5*s);
        glBegin(GL_LINES);
        glVertex2d(P.x,P.y); glVertex2d(Q.x,Q.y);
        glEnd();
        snprintf(buf,sizeof(buf),"%.3g",ends[e]);
        Q = P + o*(9*s);
        if(AxisText(Q.x,Q.y-fh/2.+o.y*fh*0.7,buf,align))
            ext=max(ext,fabs(o.x)*GraphTextWidth(buf)+fabs(o.y)*1.2*fh);
    }
    //  the name of the axis at the middle, after the values
    R half=0.5*(fabs(o.x)*GraphTextWidth(name)+fabs(o.y)*fh);
    R2 N = M + o*(gap+9*s+ext+6*s+half);
    AxisText(N.x,N.y-fh/2.,name,1);
}

void OneWindow::DrawAxes()
{ //  the axes x, y (2d) or x, y, z (3d) on the bounding box of the plot, with the tics and the values
    if(!theplot || Graph1d()) return;
    int s=GraphScale();
    R3 A,B;
    if(plotdim==3) A=Bmin3, B=Bmax3;// in 3d view of a 2d function, z = the values of the function
    else A=theplot->Pmin, B=theplot->Pmax;
    if(!(A.x<=B.x && A.y<=B.y)) return;
    SetScreenView();
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_1D);
    glDisable(GL_CLIP_PLANE0);
    glColor3d(0.,0.,0.);
    glLineWidth(s);
    axisboxes.clear();
    GLdouble x,y,z;
    if(plotdim!=3)
    { //  2d: x at the bottom of the plot, y on the left, a few pixels away from the domain
        gluProject((A.x+B.x)/2,(A.y+B.y)/2,0.,modelMatrix,projMatrix,viewport,&x,&y,&z);
        R2 C(x,y);
        DrawAxis(R3(A.x,A.y,0.),R3(B.x,A.y,0.),A.x,B.x,C,6*s,"x");
        DrawAxis(R3(A.x,A.y,0.),R3(A.x,B.y,0.),A.y,B.y,C,6*s,"y");
    }
    else
    { //  3d: for x and y the lowest of the 4 parallel edges of the box on the screen, for z the leftmost
        R3 P[2]={A,B};
        gluProject((A.x+B.x)/2,(A.y+B.y)/2,(A.z+B.z)/2,modelMatrix,projMatrix,viewport,&x,&y,&z);
        R2 C(x,y);
        const char *name[3]={"x","y","z"};
        for(int c=0;c<3;++c)
        {
            int c1=(c+1)%3, c2=(c+2)%3;
            R best=1e100;
            R3 E0,E1;
            for(int i=0;i<2;++i)
                for(int j=0;j<2;++j)
                {
                    R3 Q0=A, Q1=B;
                    Q0[c1]=Q1[c1]=P[i][c1];
                    Q0[c2]=Q1[c2]=P[j][c2];
                    R3 Q=(Q0+Q1)/2.;
                    gluProject(Q.x,Q.y,Q.z,modelMatrix,projMatrix,viewport,&x,&y,&z);
                    R v = c<2 ? y : x;
                    if(v<best) best=v, E0=Q0, E1=Q1;
                }
            DrawAxis(E0,E1,A[c],B[c],C,0.,name[c]);
        }
    }
    glLineWidth(kscreenscale);
}

void OneWindow::GraphBox(int &x0,int &y0,int &x1,int &y1) const
{
    int s=GraphScale();
    x0=85*s;
    y0=45*s;
    //  place for the labels of the axis
    if(theplot && !theplot->gylabel.empty()) x0 += 22*s;
    if(theplot && !theplot->gxlabel.empty()) y0 += 22*s;
    x1=max(width-30*s,x0+10);
    y1=max(height-40*s,y0+10);
}

void OneWindow::DrawGraph()
{ // 1d graph a la gnuplot: grid, curves, box, tics, title, labels, key and mouse coordinates
    int s=GraphScale(), fh=GraphFontHeight();
    int x0,y0,x1,y1;
    GraphBox(x0,y0,x1,y1);
    R ax=(x1-x0)/(Bmax.x-Bmin.x), ay=(y1-y0)/(Bmax.y-Bmin.y);
    vector<R> tx,ty,mx,my;// tics and minor tics
    vector<string> lx,ly;// labels of the tics
    int px,py;
    GraphAxisTics(Bmin.x,Bmax.x,max(2,(x1-x0)/(80*s)),graphlogx,tx,lx,mx,px);
    GraphAxisTics(Bmin.y,Bmax.y,max(2,(y1-y0)/(50*s)),graphlogy,ty,ly,my,py);
    char buf[128];

    //  the label of the y axis is vertical: to do first, when the window is empty
    vector<unsigned char> ypix;
    int yw=0,yh=0;
    if(!theplot->gylabel.empty())
    {
        SetScreenView();
        GraphTextRotated(theplot->gylabel.c_str(),ypix,yw,yh);
    }

    if(graphgrid)
    {
        SetScreenView();
        glLineWidth(s);
        glEnable(GL_LINE_STIPPLE);
        glLineStipple(s,0x3333);
        glColor3d(0.7,0.7,0.7);
        glBegin(GL_LINES);
        for (size_t i=0;i<tx.size();++i)
        {
            R x=x0+(tx[i]-Bmin.x)*ax;
            glVertex2d(x,y0);
            glVertex2d(x,y1);
        }
        for (size_t i=0;i<ty.size();++i)
        {
            R y=y0+(ty[i]-Bmin.y)*ay;
            glVertex2d(x0,y);
            glVertex2d(x1,y);
        }
        glEnd();
        glDisable(GL_LINE_STIPPLE);
    }

    SetView();
    for (list<OnePlot *>::iterator i= theplot->plots.begin();i != theplot->plots.end(); ++i)
        (*i)->Draw(this);

    SetScreenView();
    glLineWidth(s);
    glColor3d(0.,0.,0.);
    glBegin(GL_LINE_LOOP);
    glVertex2d(x0,y0);
    glVertex2d(x1,y0);
    glVertex2d(x1,y1);
    glVertex2d(x0,y1);
    glEnd();
    glBegin(GL_LINES);
    for (int minor=0;minor<2;++minor)
    {
        R lt=(minor ? 3 : 6)*s;// length of the tics
        const vector<R> & vx= minor ? mx : tx, & vy= minor ? my : ty;
        for (size_t i=0;i<vx.size();++i)
        {
            R x=x0+(vx[i]-Bmin.x)*ax;
            glVertex2d(x,y0);
            glVertex2d(x,y0+lt);
            glVertex2d(x,y1);
            glVertex2d(x,y1-lt);
        }
        for (size_t i=0;i<vy.size();++i)
        {
            R y=y0+(vy[i]-Bmin.y)*ay;
            glVertex2d(x0,y);
            glVertex2d(x0+lt,y);
            glVertex2d(x1,y);
            glVertex2d(x1-lt,y);
        }
    }
    glEnd();
    for (size_t i=0;i<tx.size();++i)
        GraphText(x0+(tx[i]-Bmin.x)*ax,y0-8*s-fh,lx[i].c_str(),1);
    for (size_t i=0;i<ty.size();++i)
        GraphText(x0-6*s,y0+(ty[i]-Bmin.y)*ay-fh/2.,ly[i].c_str(),2);
    if(!theplot->gtitle.empty())
        GraphText((x0+x1)/2.,y1+12*s,theplot->gtitle.c_str(),1);
    if(!theplot->gxlabel.empty())
        GraphText((x0+x1)/2.,y0-18*s-2*fh,theplot->gxlabel.c_str(),1);
    if(yw)
    {
        glPixelStorei(GL_UNPACK_ALIGNMENT,1);
        glRasterPos2d(6*s+0.5,floor((y0+y1-yw)/2.)+0.5);
        glDrawPixels(yh,yw,GL_RGB,GL_UNSIGNED_BYTE,&ypix[0]);
    }

    if(graphkey && theplot->nbcurves>1)
    { // the key: name and sample of each curve, on the top right corner
        R xs1=x1-10*s, xs0=xs1-40*s, y=y1-10*s-fh;
        for (list<OnePlot *>::iterator i= theplot->plots.begin();i != theplot->plots.end(); ++i)
        {
            int num=(*i)->num;
            glColor3d(0.,0.,0.);
            if(num < (int) theplot->glegend.size())
                GraphText(xs0-8*s,y,theplot->glegend[num].c_str(),2);
            else
            {
                snprintf(buf,sizeof(buf),"curve %d",num+1);
                GraphText(xs0-8*s,y,buf,2);
            }
            GraphColor(num,theplot->grey);
            R ym=y+fh/2.;
            if(graphstyle%3 != 2)
            {
                glLineWidth(1.5*s);
                glBegin(GL_LINES);
                glVertex2d(xs0,ym);
                glVertex2d(xs1,ym);
                glEnd();
            }
            if(graphstyle%3 != 0)
            {
                glPointSize(5*s);
                glBegin(GL_POINTS);
                glVertex2d((xs0+xs1)/2,ym);
                glEnd();
                glPointSize(1);
            }
            y -= fh+8*s;
        }
    }

    int ym=height-ymouse;
    if(xmouse>=x0 && xmouse<=x1 && ym>=y0 && ym<=y1)
    { // the coordinates of the mouse, on the bottom left corner
        R xv=Bmin.x+(xmouse-x0)/ax, yv=Bmin.y+(ym-y0)/ay;
        if(graphlogx) xv=pow(10.,xv);
        if(graphlogy) yv=pow(10.,yv);
        glColor3d(0.,0.,0.);
        snprintf(buf,sizeof(buf),"%.*g, %.*g",px+2,xv,py+2,yv);
        GraphText(6*s,6*s,buf);
    }
    glLineWidth(kscreenscale);
    SetView();// for picking..
}

void OneWindow::setLighting()
{
    if(withlight)
    {
        if(plotdim==3)
        {
            GLfloat lp0[4] = { (float) cam.x,(float)  cam.y, (float)  cam.z, 1.0F };
            glLightfv(GL_LIGHT0,GL_POSITION,lp0);

            if(debug>1)  cout << " Light pos  3d:  " << cam << endl;
        }
        else
        {
            GLfloat position[] = {(float) Pvue3.x,(float) Pvue3.y,(float) (Pvue3.z+(Bmax3.z-Bmin3.z)*3),1.f} ;
            glLightfv(GL_LIGHT0, GL_POSITION, position);

        }

        float cca=0.3,ccd=1., ccs=0.8;
        GLfloat ambient[] = {cca,cca,cca,1.0f};//differents parametres
        GLfloat diffuse[] = {ccd,ccd,ccd,1.0f};
        GLfloat specular_reflexion[] = {ccs,ccs,ccs,1.0f};
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);
        glEnable(GL_LIGHTING);//positionnement de la lumiere avec
        glLightfv(GL_LIGHT0,GL_AMBIENT,ambient);//les differents parametres
        glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuse);

        glEnable(GL_COLOR_MATERIAL);//specification de la reflexion sur les materiaux
        glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,ambient);
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,diffuse);
        if(Pretty3d())
        { // pretty rendering: smooth shading, specular light
            GLfloat amb[] = {0.32f,0.32f,0.32f,1.f}, dif[] = {0.8f,0.8f,0.8f,1.f}, spec[] = {0.45f,0.45f,0.45f,1.f};
            glLightfv(GL_LIGHT0,GL_AMBIENT,amb);
            glLightfv(GL_LIGHT0,GL_DIFFUSE,dif);
            glLightfv(GL_LIGHT0,GL_SPECULAR,spec);
            glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,spec);
            glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,40.f);
            glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL,GL_SEPARATE_SPECULAR_COLOR);
            glShadeModel(GL_SMOOTH);
        }
        else
        {
            GLfloat nospec[] = {0.f,0.f,0.f,1.f};
            glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,nospec);
            glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL,GL_SINGLE_COLOR);
            glShadeModel(GL_FLAT);
        }
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
}

OnePlotBorder::OnePlotBorder(PlotStream & f)
:OnePlot(4,3,1)
{
    long nbd;
    f>> nbd;
    data.resize(nbd);
    for(int i=0;i<nbd;++i)
    {
        long n;
        f>> n;
        data[i].resize(n+1);
        for(int j=0;j<=n;++j)
        {
            long l;
            double x,y,z;
            f >> l>> x >> y >> z;
            R3 P(x,y,z);
            Pmin=Minc(Pmin,P);
            Pmax=Maxc(Pmax,P);
            data[i][j]=make_pair(l,P);
        }
    }
    ffassert(f.good());
}


void OnePlot::GLDraw(OneWindow *win)
{
    Draw(win);

    win->changeiso=0;
    win->changearrow=0;
}

void ThePlot::DrawHelp(OneWindow *win)
{ //  the help of the current plot: only the keys which act on it (key ?), all the keys with ? again
    if(win->helpall) { DrawHelpAll(win); return; }
    //  what is plotted
    bool graph=win->Graph1d(), view3=win->plotdim==3, data2=datadim!=3;
    bool mesh2=false, fe2s=false, fe2v=false, curve=false, cplx=false, mesh3=false, fe3s=false, fe3p0=false;
    bool fe3v=false, meshS=false, fes=false, fesv=false, fel=false, feld3=false;
    for (list<OnePlot *>::iterator it= plots.begin();it != plots.end(); ++it)
    {
        int w=(*it)->what;
        OnePlotFE3 *f3=dynamic_cast<OnePlotFE3 *>(*it);
        if(w==0) mesh2=true;
        else if(w==3 || w==13) curve=true;
        else if(w==5) mesh3=true;
        else if(w==50) meshS=true;
        else if(w==55) feld3=true;
        else if(w%100==14 || w%100==15 || w%100==20 || w%100==21) fel=true;
        else if(w%10==1) fe2s=true;
        else if(w%10==2) fe2v=true;
        else if(w%10==6 && f3 && f3->p0) fe3p0=true;
        else if(w%10==6) fe3s=true;
        else if(w%10==7) fe3v=true;
        else if(w%10==8) fes=true;
        else if(w%10==9) fesv=true;
        if((w>=10 && w<20 && w!=13 && w!=14 && w!=15) || w>=100) cplx=true;
    }
    bool fe2=fe2s||fe2v, fe3=fe3s||fe3p0||fe3v, cut=mesh3||fe3;
    bool pretty=win->PrettyOn() && (fe2||fe3s||fe3v||fes||fesv);
    string what;
    if(graph) what=" 1d graph";
    else
    {
        if(fe2s) what+=" 2d function";
        if(fe2v) what+=" 2d vector field";
        if(mesh2 && !fe2) what+=" 2d mesh";
        if(curve) what+=" curves";
        if(fe3s) what+=" 3d function";
        if(fe3p0) what+=" 3d P0 function";
        if(fe3v) what+=" 3d vector field";
        if(mesh3 && !fe3) what+=" 3d mesh";
        if(fes||fesv) what+=" surface function";
        if(meshS && !fes && !fesv) what+=" surface mesh";
        if(fel||feld3) what+=" curve mesh function";
        if(cplx) what+=" (complex)";
        if(fe2 && view3) what+=", 3d view";
    }
    int i = 1;
    string t="Keys for this plot:"+what+"   (any key closes this help, ? again: all the keys)";
    win->Show(t.c_str(),i++);
    i+=1;
    win->Show("General:",i++);
    win->Show("   enter) next plot   p) previous plot   ESC) exit ffglut   w) dump PPM   W) dump PNG   *) keep the viewpoint",i++);
    i+=1;
    if(graph)
    {
        win->Show("Graph:",i++);
        win->Show(cplx ? "   k) real part, imaginary part, modulus, argument" : "   k) lines, lines and points, points",i++);
        win->Show("   G) grid   v) key of the curves   F) bold text   x) y) log scale in x, in y",i++);
        win->Show("   +) -) zoom   arrows) translate   =) reset the view   b) g) curves in black",i++);
        win->Show("   mouse: left button) translate   wheel or right button) zoom   cursor) coordinates at the bottom",i++);
        return;
    }
    win->Show("View:",i++);
    win->Show("   +) -) zoom around the cursor   arrows) translate   =) reset the view   r) refresh   j) axes",i++);
    if(view3)
        win->Show("   mouse: left button) rotate   wheel) zoom   right button) zoom (ctrl+button on mac)   alt+right) translate",i++);
    else
        win->Show("   mouse: wheel) zoom   right button) zoom (ctrl+button on mac)   alt+right) translate",i++);
    if(data2 && !curve) win->Show("   3) switch between the 2d view and the 3d view",i++);
    if(view3) win->Show(data2 ? "   z) Z) focal zoom   H) h) increase or decrease the z scale" : "   z) Z) focal zoom",i++);
    i+=1;
    win->Show("This plot:",i++);
    if(cplx) win->Show("   k) real part, imaginary part, modulus, argument",i++);
    if(fe2||fe3s||fe3p0||fes)
        win->Show("   b) g) black and white / grey or colors   i) I) bounds of the colors fitted to the window or not",i++);
    if(fe2s || fes)
        win->Show("   f) fill   n) N) fewer or more iso values   v) values of the colors   m) mesh   B) border",i++);
    if(fe2s && pretty)
        win->Show("   o) iso lines on the colors on/off (with fill)",i++);
    if(fe2v || fe3v || fesv)
        win->Show("   a) A) shorter or longer arrows   v) values of the colors",i++);
    if(mesh2 && !fe2) win->Show("   B) border of the mesh",i++);
    if(view3 && (fe2s||fe3||fes)) win->Show("   l) lighting on/off",i++);
    if(fe3s)
    {
        win->Show("   s) function on the boundary or isosurfaces   f) fill   n) N) fewer or more isosurfaces",i++);
        win->Show("   m) mesh of the boundary   B) edges of the domain   v) values of the colors",i++);
    }
    if(fe3p0)
    {
        win->Show("   s) all the elements or the elements over a threshold   [) ]) lower or raise the threshold",i++);
        win->Show("   e) E) increase or decrease the exploded view   m) edges of the elements",i++);
    }
    if(cut)
    {
        win->Show("   F1) cut plane on/off   F2) edit: the left button rotates the plane   F3) freeze: the object turns with the plane",i++);
        win->Show("   <) >) or ,) .) move the cut plane   x) X) y) Y) or shift + left button) rotate it",i++);
    }
    if(meshS || fes || fesv) win->Show("   T) normals at the surface   q) Q) fewer or more normals",i++);
    if(pretty)
        win->Show("   R) pretty rendering on/off   C) colormap   v) color bar   P) value under the mouse on/off",i++);
    else if(fe2||fe3||fes||fesv)
        win->Show("   P) value under the mouse on/off   R) pretty rendering on/off",i++);
}

void ThePlot::DrawHelpAll(OneWindow *win)
{ //  all the keys of ffglut (about 50 lines: the spacing of the lines is fitted to the window)
    int hx = kscreenscale==2 ? 32 : 20;
    win->showhx = max(kscreenscale==2 ? 25 : 15,min(hx,win->height/50));// at least the old spacing
    struct ResetHx { OneWindow *w; ~ResetHx() { w->showhx=0; } } resethx={win};
    int i = 1;
    win->Show("Enter a keyboard character in the FreeFem Graphics window in order to:",i++);

    i+=1;
    win->Show("enter) wait next plot",i++);
    win->Show("p)     previous plot (10 plots saved) ",i++);
    win->Show("ESC)   exit from ffglut (if termined) ",i++);
    win->Show("^c)    non-ignorable exit from ffglut  ",i++);
    win->Show("?)  show this help window",i++);
    win->Show("+) -)   zoom in/out  around the cursor 3/2 times ",i++);
    win->Show("=)   reset vue ",i++);
    win->Show("r)   refresh plot ",i++);
    win->Show("up, down, left, right) special keys  to  translate   ",i++);
    win->Show("3)   switch 3d/2d plot (in test)  keys : ",i++);
    win->Show(" z) Z) (focal zoom unzoom)  ",i++);
    win->Show(" H) h) switch increase or decrease the Z scale of the plot ",i++);
    win->Show("mouse motion)    ",i++);
    win->Show("   - left button)  rotate    ",i++);
    win->Show("   - wheel)  zoom around the cursor ",i++);
    win->Show("   - right button)       zoom        (ctrl+button on mac) ",i++);
    win->Show("   - right button +alt)  translate    (alt+ctrl+button on mac)",i++);

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
    win->Show("w)  window dump in file ffglutXXXX.ppm    W)  window dump in file ffglut_XXXX.png ",i++);
    win->Show("*)  keep/unkeep viewpoint for next plot",i++);
    win->Show("k)  complex data / change view type ",i++);
    win->Show("T)  show normal at element surface (only meshS)",i++);
    win->Show("q) Q) decrease or increase the number of Normal ploted at element surface (only meshS)",i++);
    win->Show("1d graph (curves or FE function on meshL along x): k) lines/points  G) grid  v) key  F) bold  x) y) log scale  left button) translate",i++);
    win->Show("s)  3d scalar function: switch between the function on the boundary or the isosurfaces",i++);
    win->Show("    3d P0 function: s) switch all the elements or the elements greater than a threshold  [) ]) move the threshold",i++);
    win->Show("    3d P0 function: e) E) increase or decrease the exploded view",i++);
    win->Show("P)  show the value of the function under the mouse (bottom left) on/off",i++);
    win->Show("j)  axes x, y (and z in 3d) with the values on/off",i++);
    win->Show("?)  the help of the current plot only, ? again: this complete help",i++);
    win->Show("o)  2d with fill: the iso lines on the colors on/off   (n, N: number of iso lines)",i++);
    win->Show("R)  pretty rendering on/off (smooth colors and shading)   C) colormap   v) color bar",i++);
    win->Show("F1) 3d function or 3d mesh: cut plane on/off   F2) edit with the mouse   F3) freeze   <) >) move",i++);
    win->Show("    x) X) y) Y) or shift + left button) rotate the cut plane",i++);
    win->Show("any other key : nothing ",++i);
}

void ThePlot::dyn_bfv(OneWindow *win,R & fmn,R &fmx,R & vmn2,R & vmx2) const
{
    fmn=+1e100;
    fmx=-fmn;
    vmx2=0;
    vmn2=+1e100;
    for (list<OnePlot *>::const_iterator i= plots.begin();i != plots.end(); ++i)
    {
        if(*i)  (*i)->dyn_bfv(win,fmn,fmx,vmn2,vmx2) ;
    }
    if(debug>4)
        cout << "dyn_bfv  " << fmn << " " << fmx << " " << vmn2 << " " << vmx2 <<endl;
    if(fmn>fmx) fmn=fmin,fmx=fmax;
    if(vmn2<0) vmn2 = 0;
    if(vmn2>vmx2) vmn2=0,vmx2=vmax2;
}

bool ThePlot::GraphBound(R2 & Pmn,R2 & Pmx,bool lx,bool ly) const
{ // true if all the plots can be drawn in a 1d graph
    bool ret= !plots.empty();
    for (list<OnePlot *>::const_iterator i= plots.begin();i != plots.end(); ++i)
        ret = (**i).Graph1d(Pmn,Pmx,lx,ly) && ret;
    return ret;
}

static string GraphTrim(const string & s)
{
    size_t a=s.find_first_not_of(" \t"), b=s.find_last_not_of(" \t");
    return a==string::npos ? string() : s.substr(a,b-a+1);
}

void ThePlot::GraphOptions()
{ //  the options of the 1d graph in the comment of the plot, the fields are separated by |
  //    cmm="title | xlabel=... | ylabel=... | legend=name1,name2,... | logx | logy"
    glogx=glogy=-1;
    if(!cm) return;
    string c=*cm+"|";
    for(size_t a=0,b;(b=c.find('|',a))!=string::npos;a=b+1)
    {
        string f=GraphTrim(c.substr(a,b-a));
        if(f.compare(0,7,"xlabel=")==0) gxlabel=GraphTrim(f.substr(7));
        else if(f.compare(0,7,"ylabel=")==0) gylabel=GraphTrim(f.substr(7));
        else if(f.compare(0,7,"legend=")==0)
        {
            string l=f.substr(7)+",";
            for(size_t la=0,lb;(lb=l.find(',',la))!=string::npos;la=lb+1)
                glegend.push_back(GraphTrim(l.substr(la,lb-la)));
        }
        else if(f=="logx") glogx=1;
        else if(f=="logy") glogy=1;
        else if(f=="loglog") glogx=glogy=1;
        else if(!f.empty())
            gtitle += (gtitle.empty() ? "" : " | ") + f;// the title
    }
}

bool ThePlot::NextCase()
{
    bool ret=false;
    for (list<OnePlot *>::iterator i= plots.begin();i != plots.end(); ++i)
        ret =  (**i).NextCase() || ret; //  alway  run NextCase ..
    return  ret;
}

void ThePlot::PrettyPlots(OneWindow *win,bool & prettyiso,bool & prettyarrow) const
{ //  the plots drawn with the colormap of the pretty rendering
    prettyiso=prettyarrow=false;
    if(!win->PrettyOn()) return;
    for (list<OnePlot *>::const_iterator i= plots.begin();i != plots.end(); ++i)
    {
        const OnePlotFE3 *f3=dynamic_cast<const OnePlotFE3 *>(*i);
        const OnePlotFES *fs=dynamic_cast<const OnePlotFES *>(*i);
        const OnePlotFE<Mesh2> *f2=dynamic_cast<const OnePlotFE<Mesh2> *>(*i);
        const OnePlotFE<Mesh> *f1=dynamic_cast<const OnePlotFE<Mesh> *>(*i);
        int w2 = f2 ? f2->what%10 : (f1 ? f1->what%10 : 0);
        if(win->Pretty3d())
        {
            if(f3 && f3->what%10==6 && !f3->p0) prettyiso=true;
            if(f3 && f3->what%10==7) prettyarrow=true;
            if(fs && fs->what%10==8 && !win->changePlotdim) prettyiso=true;
        }
        if(w2==1) prettyiso=true;
        if(w2==2) prettyarrow=true;
    }
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
        win->DefaultView(0);
    }

    if(win->Graph1d())
    {
        win->DrawGraph();
        changeViso=false;
        changeVarrow=false;
        changeColor=false;
        changeBorder=false;
        changeFill=false;
        return;
    }
    //  pretty rendering: the iso lines with the colormap
    cmap = win->PrettyOn() ? (grey ? 4 : win->cmapkind) : -1;
    DataRange(*this,cmapfmn,cmapfmx);
    win->SetView();
    for (list<OnePlot *>::iterator i= plots.begin();i != plots.end(); ++i)
        (*i)->Draw(win);
    //  the display lists are up to date: they are only called at the next display (fast with a
    //  fine mesh), until a key or the zoom asks to rebuild them (see OnePlot::GLDraw, never called)
    win->changeiso=false;
    win->changearrow=false;

    //  the axes with the tics and the values (key j)
    if(win->axes)
    {
        win->DrawAxes();
        win->SetView();
    }

    //  pretty rendering: the color bars of the functions drawn with the colormap (key v)
    bool prettyiso=false, prettyarrow=false;
    PrettyPlots(win,prettyiso,prettyarrow);
    bool colorbar = win->prettybar && (prettyiso || prettyarrow);
    if(colorbar)
    {
        R xr = win->width-95*GraphScale();
        if(prettyiso)
        {
            R fmn,fmx;
            DataRange(*this,fmn,fmx);
            const vector<R> *marks=0;// the values of the isosurfaces
            for (list<OnePlot *>::iterator i= plots.begin();i != plots.end() && !marks; ++i)
            {
                OnePlotFE3 *f3=dynamic_cast<OnePlotFE3 *>(*i);
                OnePlotFE<Mesh2> *f2=dynamic_cast<OnePlotFE<Mesh2> *>(*i);
                OnePlotFE<Mesh> *f1=dynamic_cast<OnePlotFE<Mesh> *>(*i);
                if(f3 && f3->showisovals) marks=&f3->isovals;
                if(f2 && f2->showisovals) marks=&f2->isovals;
                if(f1 && f1->showisovals) marks=&f1->isovals;
            }
            win->DrawColorBar(fmn,fmx,xr,0,marks);
        }
        if(prettyarrow)
            win->DrawColorBar(0.,sqrt(vmax2),prettyiso ? 30.*GraphScale() : xr,"|vector|");
        win->SetView();
    }
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
            if(withiso && !(colorbar && prettyiso))
            {
                win->PlotValue(Viso,k0,"IsoValue");
                k0=1+Viso.N();
            }
            if(witharrow && !(colorbar && prettyarrow))
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

void ThePlot::colorv(int i,R val)
{ //  pretty rendering: the color of the value in the colormap, otherwise the color i of the table
    if(cmap>=0)
    {
        float c[3];
        ColormapRGB(cmap,(val-cmapfmn)/(cmapfmx-cmapfmn),c);
        glColor3fv(c);
    }
    else
        color(i);
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
plotdim(2),changePlotdim(false),theta(30.*my_pi/180.),phi(20.*my_pi/180.),dcoef(1),focal(20.*my_pi/180.),
datadim(1), winnum(0), keepPV(0), pNormalT(0)

{
    string * pdffile=0,*svgfile=0;// add for  fujiwara used today !!!
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
    ZScale=-1; // auto
    ArrowSize=-1;
    if(old) {
        grey=old->grey;
    }
    greyo=grey;
    drawborder=true;
    drawmeshes=false;
    add=false;
    keepPV=false;
    pNormalT=false;
    nbN=1; // arbitrary initialisation
    cmap=-1;// the classic colors
    cmapfmn=0;
    cmapfmx=1;
    echelle=1.;

    Pmin=R3(+dinfty,+dinfty,+dinfty);
    fmin = +dinfty;
    fmax = -dinfty;
    Pmax=R3(-dinfty,-dinfty,-dinfty);
    vmax2=0;

    coefr=1;
    long dimpp=0, dimHatpp=0;
    long cas;
#define READ_VTK_PARAM(index,type)                    \
case 20+index: {type dummy; fin >= dummy;} break;

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
                case 20: fin >> echelle;break;
                default:
                    cout << "Fatal error: Unknown  case  : " << cas <<endl;
                    ffassert(0);
                    break;
            }
        else if(version == 3 || version == 4)
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
                case 20: fin >= echelle;break;
                case 21: fin >= ZScale;break;

                    // unsed parameter ...
                    //READ_VTK_PARAM(1,double); // ZScale
                    READ_VTK_PARAM(2,bool); // WhiteBackground
                    READ_VTK_PARAM(3,bool); // OpaqueBorders
                    READ_VTK_PARAM(4,bool); // BorderAsMesh
                    READ_VTK_PARAM(5,bool); // ShowMeshes
                    READ_VTK_PARAM(6,long); // ColorScheme
                    READ_VTK_PARAM(7,long); // ArrowShape
                    //       READ_VTK_PARAM(8,double); // ArrowSize
                case 28: fin >= ArrowSize; break;
                    READ_VTK_PARAM(9,long); // ComplexDisplay
                    READ_VTK_PARAM(10,bool); // LabelColors
                    READ_VTK_PARAM(11,bool); // ShowAxes
                    READ_VTK_PARAM(12,bool); // CutPlane
                    READ_VTK_PARAM(13,KN<double>); // CameraPosition
                    READ_VTK_PARAM(14,KN<double>); // CameraFocalPoint
                    READ_VTK_PARAM(15,KN<double>); // CameraViewUp
                    READ_VTK_PARAM(16,double); // CameraViewAngle
                    READ_VTK_PARAM(17,KN<double>); // CameraClippingRange
                    READ_VTK_PARAM(18,KN<double>); // CutPlaneOrigin
                    READ_VTK_PARAM(19,KN<double>); // CutPlaneNormal
                 case 40: fin >= winnum; break;
      
                //case 43: fin >= winnum; break;
                case 43: fin >= pNormalT;break;
                case 44: fin >= pdffile;break;
                case 45: fin >= svgfile;break;

                default:
                    static int nccc=0;
                    if(nccc++<5)
                        cout << " Skip Unknown case " << cas <<" (ffglut is too old ?)\n";
                    fin.SkipData();
                    break;
            }
        else ffassert(0);
        ffassert(fin.good() && ! fin.eof());
    }
    if(dimpp) plotdim=dimpp;
    ffassert(cas==PlotStream::dt_endarg);
    if((debug > 2))
    {  cout << "    ***** get ::: ";
        cout << " Window num " << winnum << ", ";
        cout << " coeff " << coeff <<", ";
        if(cm)
            cout << " cm " << *cm <<", ";
        if(wait)    cout << " wait " ;
        if(fill)    cout << " fill " ;
        if(value)    cout << " value " ;
        if(bw)    cout << " bw " ;
        if(grey)    cout << " grey " ;
        if(drawborder)    cout << " drawborder " ;
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
            {
                if(version==2)
                    cout << i << " nt/nv " << l << " "  <<Ths[l]->nt << " " << Ths[l]->nv << endl;
                else
                    cout << i << " nt/nv " << l << " "  <<Ths2[l]->nt << " " << Ths2[l]->nv << endl;
            }
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
                Ths[l]->increment();
            }
            else
            {
                Ths2[l]=old->Ths2[ll];
                Ths2[l]->increment();
            }

        }

    }
    long nbmeshes3=0;
    long nbmeshesS=0;
    long nbmeshesL=0;
    
    int getMesh3Type = fin.GetMeshes3(); // 0 : dt_mesh3, 1 : dt_meshS, 2 : dt_meshL, 3 : dt_plots
    if (getMesh3Type!=3) //  read GetPlots if false ...
    {
     //  There are 3D volume solution
     if(getMesh3Type==0)
     {

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
                     cout << "mesh3 i: "<< i << " nt/nv " << l << " "  <<Ths3[l]->nt << " "
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
         getMesh3Type=fin.GetMeshes3();
     }
     //  There are 3D surface solution
     if(getMesh3Type==1) {
         fin >> nbmeshesS;
         if((debug > 2)) cout << " read nb : meshS " << nbmeshesS << endl;
         ThsS.resize(nbmeshesS);
         for(int i=0;i<nbmeshesS;++i)
             ThsS[i]=0;
         for(int i=0;i<nbmeshesS;++i)
         {
             long l;
             fin >> l;
             if(l>=0)
             {
                 if((debug > 3)) cout << "meshS i: "<< " read meshS " << i  << " -> " << l
                     << "  " <<nbmeshesS << endl;
                 l--;
                 ffassert(l>=0 && l < nbmeshesS);
                 ffassert(ThsS[l]==0);
                 fin >>ThsS[l] ;
                 if((debug > 3))
                     cout << "meshS i: "<< i << " nt/nv " << l << " "  <<ThsS[l]->nt << " "
                     << ThsS[l]->nv << endl;
                 ffassert(fin.good());
             }
             else // Add FH optimisation FH 11/12/2008 (not use to day)
             {// the mesh is already in the previous plot with number ll
                 ffassert(l==-1);
                 long ll;
                 fin >> l>> ll; // read l and ll
                 ffassert(old);
                 ThsS[l]=old->ThsS[ll];
                 ThsS[l]->increment(); //
             }
        }
        getMesh3Type=fin.GetMeshes3();
     }
     //  There are 3D line solution
     if(getMesh3Type==2) {
            fin >> nbmeshesL;
            if((debug > 2)) cout << " read nb : meshL " << nbmeshesL << endl;
            ThsL.resize(nbmeshesL);
            for(int i=0;i<nbmeshesL;++i)
                ThsL[i]=0;
            for(int i=0;i<nbmeshesL;++i)
            {
                long l;
                fin >> l;
                if(l>=0)
                {
                    if((debug > 3)) cout << " read meshL " << i  << " -> " << l
                        << "  " <<nbmeshesL << endl;
                    l--;
                    ffassert(l>=0 && l < nbmeshesL);
                    ffassert(ThsL[l]==0);
                    fin >>ThsL[l] ;
                    if((debug > 3))
                        cout << "meshL i: "<< i << " nt/nv " << l << " "  <<ThsL[l]->nt << " "
                        << ThsL[l]->nv << endl;
                    ffassert(fin.good());
                }
                else // Add FH optimisation FH 11/12/2008 (not use to day)
                {// the mesh is already in the previous plot with number ll
                    ffassert(l==-1);
                    long ll;
                    fin >> l>> ll; // read l and ll
                    ffassert(old);
                    ThsL[l]=old->ThsL[ll];
                    ThsL[l]->increment(); //
                }
            }
            getMesh3Type=fin.GetMeshes3();
        }
     else if (getMesh3Type==4)

      fin.GetPlots();
    }
  
    long nbplot;
    int iso3d=0, iso3dB=0;
    fin >>nbplot;
    if((debug > 2)) cout << " nb item plot " << nbplot << endl;
    for(int i=0;i<nbplot;++i)
    {
        long what;
        OnePlot *p=0;
        fin >> what;
        long imsh;
        if((what !=3 && what != 13  )&& !uaspectratio) aspectratio= true;
        if((debug > 2)) cout << "    plot  " << i << " what " << what << " as : " << aspectratio << endl;

        if(what==-1)  // gestion of error (empty plot)
            p = new OnePlotError(fin);
        else if(what==0)
        {
            fin >> imsh;
             if((debug > 9)) cout << " plot : mesh " << imsh << endl;
            if(version==2)
                p=new OnePlotMesh<Mesh>(Ths[imsh-1]);
            else
                p=new OnePlotMesh<Mesh2>(Ths2[imsh-1]);

        }
        else if (what==1 || what==2 || what==11 || what==12 )
        {
            fin >> imsh;
            if(what%10==1) withiso=true;
            else if (what%10==2) witharrow=true;
            if((debug >= 10)) cout << " plot : mesh " << imsh << endl;
            ffassert(imsh>0 && imsh <=nbmeshes);
            if(version==2)
                p=new OnePlotFE<Mesh>(Ths[imsh-1],what,fin);
            else
                p=new OnePlotFE<Mesh2>(Ths2[imsh-1],what,fin);
        }
        else if(what==3)
            p=new OnePlotCurve(fin);
        else if(what==13)
            p=new OnePlotCurve(fin,4,this);// add zz and color ..
        else if(what==4)
            p=new OnePlotBorder(fin);
        
        // volume 3D meshes (mesh3 or mesh3 only in case with pointer on meshS or/and meshL
        else if (what == 5)
        {
            fin >> imsh;
            p=new OnePlotMesh3(Ths3[imsh-1]);if((debug > 9)) cout << " plot : mesh3 (volume 3D) " << imsh << endl;
        }
        // surface 3D meshes
        else if(what==50)
        {
            fin >> imsh;if((debug > 9)) cout << " plot : meshS (surface 3D) " << imsh << endl;
            p=new OnePlotMeshS(ThsS[imsh-1]);
        }
        // line 3D meshes
        else if (what==55)
        {
            fin >> imsh;if((debug > 9)) cout << " plot : meshL (line 3D) " << imsh << endl;
            p=new OnePlotMeshL(ThsL[imsh-1]);
        }
        else if (what==6  || what==7|| what==16  || what==17)
        {
            iso3d++;
            fin >> imsh;
            if(what==6||what==16 ) withiso=true;
            else witharrow=true;

            if((debug > 10)) cout << " plot : mesh3 (volume) " << imsh << endl;
            ffassert(imsh>0 && imsh <=nbmeshes3);
            p=new OnePlotFE3(Ths3[imsh-1],what,fin);
        }
        else if (what==8  || what==9 || what==18  || what==19)
        {
            iso3dB++;
            fin >> imsh;
            if(what==8||what==18 ) withiso=true;
            else if (what%10==9) witharrow=true;
            if((debug > 10)) cout << " plot : meshS (3D surface) " << imsh << endl;
            ffassert(imsh>0 && imsh <=nbmeshesS);

            p=new OnePlotFES(ThsS[imsh-1],what,fin);
        }
        else if (what==14  || what==15 || what==20  || what==21)
        {
            iso3dB++;
            fin >> imsh;
            if(what==14||what==20 ) withiso=true;
            else if (what==15||what==21 ) witharrow=true;
            if((debug > 10)) cout << " plot : meshL (3D curve) " << imsh << endl;
            ffassert(imsh>0 && imsh <=nbmeshesL);
            p=new OnePlotFEL(ThsL[imsh-1],what,fin);
        }
        else if (what == 31)
        {
            p=new OnePlotHMatrix(what,fin);
        }
        else
        {
            cout << "Bizarre unknown what :  " << what<< endl;
            ffassert(0);
        }
        ffassert(p);
        plots.push_back(p);
        p->bb(Pmin,Pmax);
        p->bfv(fmin,fmax,vmax2);
        plotdim=max(plotdim,p->dim);
        ffassert(fin.good());
        datadim=max(datadim,p->dim);
    }
    if(dimpp) plotdim=dimpp;
    // 1d graph a la gnuplot if the plot contains only 2d curves
    //  or scalar FE functions on a meshL along the x axis (use dim=3 to get the 3d plot)
    {
        R2 Gmn(dinfty,dinfty),Gmx(-dinfty,-dinfty);
        graph1d = GraphBound(Gmn,Gmx);
    }
    nbcurves=0;
    for (list<OnePlot *>::iterator i= plots.begin();i != plots.end(); ++i)
        (**i).num=nbcurves++;
    glogx=glogy=-1;
    if(graph1d)
    {
        GraphOptions();
        if(!dimpp) plotdim=2;
        if(plotdim==2 && !uaspectratio) aspectratio=false;
    }
   
    if(Niso==0)
        Niso = 20;// also in 3d (function on the boundary), the isosurfaces are a part of the iso values
    if(iso3dB) Niso=20;
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

    z0= fminT -(fmaxT-fminT)*0.01;
    if((debug > 2)) cout << "               data bound: " << PminT << " " << PmaxT
        << " fmin == " << fminT << "  " << fmaxT
        << " z0 " << z0 <<  endl;
    fin.GetEndPlot();
    Viso.resize(Niso);
    Varrow.resize(Narrow);

    SetColorTable(Max(Niso,Narrow)+4) ;
    SetDefIsoV(Niso,Narrow,fmin,fmax) ;
    
}


void ThePlot::SetDefIsoV(int niso,int narr,double fmn,double fmx,double vmn,double vmx)
{
    niso = max(niso,2);
    narr = max(narr,2);

    bool dyni=niso!=Viso.N() , dyna = narr != Varrow.N();
    R d,x;
    if(debug>3 && !(fmx>fmn) )
        cout << " SetDefIsoV  (not) " << endl;
    if(debug>3)
        cout << " SetDefIsoV  " << dyni << " " <<  fmn << " " << fmx << "pViso= " << pViso << " / "<< dyna << " "  << vmn << " " << vmx << endl;

    // resize !!!
    Viso.resize(niso);
    Niso=Viso.N();
   
    Varrow.resize(narr);
    Narrow=Varrow.N();

    if( fmx>fmn)
    {
        dyni = niso !=Viso.N() ;
        //       pViso = pViso && niso ==Viso.N(); FH:  debile !!!!!
        d =  (fmx-fmn)/(max(Niso-2,1L))  ;
        x =  (fmn+fmx)/2-d*0.5*(Niso-1);
    }
    else
    {
        d = 1 ? (fmaxT-fminT)/(max(Niso-2,1L))  : (fmaxT-fminT)/(Niso-1);
        x = 1 ? (fminT+fmaxT)/2-d*0.5*(Niso-1) :fminT+d/2;
    }
    if(!pViso || dyni)
    {
        for (int i = 0;i < Niso;i++)
        {Viso[i]=x;x +=d; }
    }
  

    if(vmx>vmn)
    {
        dyna=true;
        x = sqrt(vmn);
        d = (sqrt(vmx)-x)/(Narrow-1.1);
    }
    else
    {
        x=0;
        if(debug>10)
            cout << "vmax2=  " << vmax2 << endl;
        d= sqrt(vmax2)/(Narrow-1.1);
    }
    if (!pVarrow || dyna)
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
    int hx= 20;// 15 was too small for the font: the lines overlapped
    if(kscreenscale==2) hx = 32; // OK ...
    if(showhx>0) hx=showhx;// the complete help: the lines must fit in the window
    int ix= width/20;
    int iy= height-hx*i;
    plot(ix,iy,str);
}

void  FillRectRasterPos(R x0,R y0,R x1,R y1)
{
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
    if((debug > 10)) cout << "PlotValue:" << cmm << " " << k0 << " " << width << " " <<height << " " << kscreenscale << endl;
    R xmin=0,xmax=width,ymin=0,ymax=height;
    if((debug > 10)) cout << "PlotValue " << Viso << endl;

    R dx=(xmax-xmin);
    R dy=(ymax-ymin);
    //  10 points
    R h=14*kscreenscale;
    if( kscreenscale==2) h = 24;// BofBof F.H
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
    y -=  ho*2;
    for (int i=0;i<Viso.N();i++)
    {
        if((debug > 10)) cout << " x0 : " << x0<< " " << y << " " << h << " v = " << Viso[i] << endl;
        this->color(i+4);
        FillRectRasterPos(x0,y,x0+h,y+h);
        plot(x0+ho,y+3*h/10,Viso[i]);
        y -=  h;
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

void  plot(double xx,double yy,const char *cmm,int font)
{
    glRasterPos2f(xx,yy);
    float x[4];
    glGetFloatv(GL_CURRENT_RASTER_POSITION,x);
    if((debug > 10)) cout<<"avant x : "<<x[0]<<" y : "<<x[1]<<" z : "<<x[2]<< " " << xx <<" " << yy << endl;
    void * glut_font=GLUT_BITMAP_HELVETICA_18;
    if(kscreenscale==2) glut_font=GLUT_BITMAP_TIMES_ROMAN_24;
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
        h = h - floorf(h);
        h*=6.;
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
    int nbl =0;
     R2 PQ[5];
    R  eps2= Min(R2(Pt[0],Pt[1]).norme2(),R2(Pt[0],Pt[2]).norme2(),R2(Pt[1],Pt[2]).norme2() )*1e-8;
    R tmn=Min(ff[0],ff[1],ff[2]), tmx=Max(ff[0],ff[1],ff[2]);
    for(int l=0;l< NbIso;l++)  /*    loop on the level curves */
    {
        if(Viso[l] < tmn || Viso[l] > tmx) continue;// fast: only the iso lines in the triangle
        R xf = Viso[l];
        int im=0;
        for(int i=0;i<3;i++) // for the  3 edges
        {
            int j = (i+1)%3;
            R fi=(ff[i]);
            R fj=(ff[j]);

            if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))
            {
                if (Abs(fi-fj)<=0.1e-10)     /* one side must be drawn */
                {
                    colorv(l+4,Viso[l]);
                    glVertex3f(Pt[i].x, Pt[i].y, xf*rapz);
                    glVertex3f(Pt[j].x, Pt[j].y, xf*rapz);
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
            colorv(l+4,Viso[l]);
            if( R2(PQ[0],PQ[1]).norme2() > eps2 )
            {
                if(nbl++==0 && ffbatch!=2)    glBegin(GL_LINES);
                glVertex3f(PQ[0].x, PQ[0].y, xf*rapz);
                glVertex3f(PQ[1].x, PQ[1].y, xf*rapz);
            }
        }
    }
    if(nbl && ffbatch!=2) glEnd();

}

// draw iso values for FE surface
void ThePlot::DrawIsoT(const R3 Pt[3],const R ff[3],const R3 Nt[3],const R * Viso,int NbIso,bool changePlotdim,int viewdim,R rapz)
{
    int nbl=0;
    R3 PQ[5];
    R3 NQ[5];
    R  eps2= Min(R3(Pt[0],Pt[1]).norme2(),R3(Pt[0],Pt[2]).norme2(),R3(Pt[1],Pt[2]).norme2() )*1e-8;
    R tmn=Min(ff[0],ff[1],ff[2]), tmx=Max(ff[0],ff[1],ff[2]);
    for(int l=0;l< NbIso;l++) { /*    loop on the level curves */
        if(Viso[l] < tmn || Viso[l] > tmx) continue;// fast: only the iso lines in the triangle
        R xf = Viso[l];
        int im=0;
        for(int i=0;i<3;i++) { // for the  3 edges
            int j = (i+1)%3;
            R fi=(ff[i]);
            R fj=(ff[j]);

            if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf))) {
                if (Abs(fi-fj)<=0.1e-10) {    /* one side must be drawn */
                    colorv(l+4,Viso[l]);
                    if(nbl++==0 && ffbatch!=2)  glBegin(GL_LINES);
                    if(viewdim==3)
                        if(changePlotdim){
                            glVertex3f(Pt[i].x+NQ[i].x*xf*rapz, Pt[i].y+NQ[i].y*xf*rapz, Pt[i].z+NQ[i].z*xf*rapz);
                            glVertex3f(Pt[j].x+NQ[j].x*xf*rapz, Pt[j].y+NQ[j].y*xf*rapz, Pt[j].z+NQ[j].z*xf*rapz);
                        }
                        else {
                            glVertex3f(Pt[i].x, Pt[i].y, Pt[i].z);
                            glVertex3f(Pt[j].x, Pt[j].y, Pt[j].z);
                        }
                     else if(viewdim==2) {
                            glVertex3f(Pt[i].x, Pt[i].y, xf*rapz);
                            glVertex3f(Pt[j].x, Pt[j].y, xf*rapz);
                        }
                }
                else {
                    R  xlam=(fi-xf)/(fi-fj);
                    NQ[im]   = Nt[i] * (1.F-xlam)  +  Nt[j]* xlam;
                    PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
                }
            }
        }

        if (im>=2)
        {
            colorv(l+4,Viso[l]);
            if( R3(PQ[0],PQ[1]).norme2() > eps2 ) {
                 if(nbl++==0 && ffbatch!=2)  glBegin(GL_LINES);
                 if(viewdim==3)
                     if(changePlotdim){
                         glVertex3f(PQ[0].x+NQ[0].x*xf*rapz, PQ[0].y+NQ[0].y*xf*rapz, PQ[0].z+NQ[0].z*xf*rapz);
                         glVertex3f(PQ[1].x+NQ[1].x*xf*rapz, PQ[1].y+NQ[1].y*xf*rapz, PQ[1].z+NQ[1].z*xf*rapz);
                     }
                     else {
                         glVertex3f(PQ[0].x, PQ[0].y, PQ[0].z);
                         glVertex3f(PQ[1].x, PQ[1].y, PQ[1].z);
                     }
                     else if(viewdim==2) {
                         glVertex3f(PQ[0].x, PQ[0].y, xf*rapz);
                         glVertex3f(PQ[1].x, PQ[1].y, xf*rapz);
                     }
            }
        }
    }
    if( nbl && ffbatch!=2) glEnd();
}




// draw filling values for 2d FE
void ThePlot::DrawIsoTfill(const R2 Pt[3],const R ff[3],const R * Viso,int NbIso, R rapz)
{
    R2 PQ[10];
    R z[10];

    R eps= (Viso[NbIso-1]-Viso[0])*1e-6;
    R tmn=Min(ff[0],ff[1],ff[2]), tmx=Max(ff[0],ff[1],ff[2]);
    for(int l=1;l< NbIso;l++)  //   loop on the level curves
    {
        //  fast: only the bands which cut the triangle (a lot of time with a fine mesh)
        if(Viso[l] < tmn-eps || Viso[l-1] > tmx+eps) continue;
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
                    if(verbosity>999) cout << im-1 << " PQ " << PQ[im-1] << endl;
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
            colorv(l+4,0.5*(Viso[l-1]+Viso[l]));
            R3 P[10];
            for(int i=0;i<im;++i)
                P[i]= R3(PQ[i].x,PQ[i].y,z[i]*rapz);
            R3 N(R3(P[0],P[1])^R3(P[0],P[2]));
            N /= N.norme();
            if(N.z<0) N = -N;
            glNormal3d(N.x,N.y,N.z);
            ffPolygon(P,im);

        }
    }
}

// draw filling values for 3d FE surface
void ThePlot::DrawIsoTfill(const R3 Pt[3],const R ff[3],const R3 Nt[3],const R * Viso,int NbIso,bool changePlotdim,int viewdim,R rapz)
{
    R3 PQ[10];
    R3 NQ[10];
    R z[10];
    R eps= (Viso[NbIso-1]-Viso[0])*1e-6;
    R tmn=Min(ff[0],ff[1],ff[2]), tmx=Max(ff[0],ff[1],ff[2]);
    for(int l=1;l< NbIso;l++)  //   loop on the level curves
    {
        //  fast: only the bands which cut the triangle (a lot of time with a fine mesh)
        if(Viso[l] < tmn-eps || Viso[l-1] > tmx+eps) continue;
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
            if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf))) {
                if (Abs(fi-fj)>=0.1e-20) {
                    R  xlam=(fi-xf)/(fi-fj);
                    z[im] =  fi * (1.F-xlam)  +  fj* xlam;
                    NQ[im] =  Nt[i] * (1.F-xlam)  +  Nt[j] * xlam;
                    PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
               }
            }
            xf = xxfh;
            if(((fi<=xf)&&(fj>=xf))||((fi>=xf)&&(fj<=xf)))  {
                if (Abs(fi-fj)>=0.1e-20) {
                    R  xlam=(fi-xf)/(fi-fj);
                    z[im] =  fi * (1.F-xlam)  +  fj* xlam;
                    NQ[im] =  Nt[i] * (1.F-xlam)  +  Nt[j] * xlam;
                    PQ[im++]   = Pt[i] * (1.F-xlam)  +  Pt[j]* xlam;
                }
            }
            if (  xfb-eps <=fj  && fj <= xfh+eps)
            {z[im]=fj,NQ[im]=Nt[j]; PQ[im++] = Pt[j];}
        }
        if (im>2) {
            colorv(l+4,0.5*(Viso[l-1]+Viso[l]));
            R3 P[10];
            for(int i=0;i<im;++i)
             if(viewdim==3)
                 if(changePlotdim)
                    P[i]= R3(PQ[i].x+rapz*z[i]*NQ[i].x,PQ[i].y+rapz*z[i]*NQ[i].y,PQ[i].z+rapz*z[i]*NQ[i].z);
                 else
                     P[i]= R3(PQ[i].x,PQ[i].y,PQ[i].z);
               else if(viewdim==2)
                P[i]= R3(PQ[i].x,PQ[i].y,rapz*z[i]);
            
            R3 N(R3(P[0],P[1])^R3(P[0],P[2]));
            N /= N.norme();
            if(N.z<0) N = -N;
            glNormal3d(N.x,N.y,N.z);

            ffPolygon(P,im);
        }
    }
}
// draw filling values for 3d FE curve
void ThePlot::DrawIsoEfill(const R3 Pt[2],const R ff[2],const R3 Nt[2],const R * Viso,int NbIso,bool changePlotdim,int viewdim, R rapz)
{
    R3 PQ[10], NQ[10];
    R z[10];
    R eps= (Viso[NbIso-1]-Viso[0])*1e-6;
    for(int l=1;l< NbIso;l++)  //   loop on the level curves
    {
        R xfb = Viso[l-1];
        R xfh = Viso[l];
        assert(xfb < xfh);
        int im=0;
        R fi=(ff[0]), fj=(ff[1]), xxfb=xfb, xxfh=xfh;
        int sens = 0;
        if (fj<fi ) {Exchange(xxfb,xxfh); sens=1;}
        
        if (Abs(fi-fj)>=0.1e-20) {
            if ((xxfh>=fi && xxfb<=fj) || (xxfh<=fi && xxfb>=fj)) {
                PQ[im] = Pt[0];
                NQ[im] = Nt[0];
                z[im++] = fi;
                if ((!sens && (xxfb >= fi)) || (sens && (xxfb <= fi))) {
                    R  xlam=(fi-xxfb)/(fi-fj);
                    PQ[im-1] = Pt[0] * (1.F-xlam)  +  Pt[1]* xlam;
                    NQ[im-1] = Nt[0] * (1.F-xlam)  +  Nt[1]* xlam;
                    z[im-1] = fi * (1.F-xlam)  +  fj* xlam;
                }
                PQ[im] = Pt[1];
                NQ[im] = Nt[1];
                z[im++] = fj;
                if ((!sens && (xxfh <= fj)) || (sens && (xxfh >= fj))) {
                    R  xlam=(fi-xxfh)/(fi-fj);
                    PQ[im-1] = Pt[0] * (1.F-xlam)  +  Pt[1]* xlam;
                    NQ[im-1] = Nt[0] * (1.F-xlam)  +  Nt[1]* xlam;
                    z[im-1] = fi * (1.F-xlam)  +  fj* xlam;
                }
            }
        }
        else if (xfb-eps <=fj  && fj <= xfh+eps) {
            PQ[im] = Pt[0]; NQ[im] = Nt[0]; z[im++] =fi;
            PQ[im] = Pt[1]; NQ[im] = Nt[1]; z[im++] =fj;
        }
        
        if (im>1) {
            color(l+4);
            R3 P[10];
            for(int i=0;i<im;++i)
                if(viewdim==3)
                    if(changePlotdim)
                        P[i]= R3(PQ[i].x+rapz*z[i]*NQ[i].x,PQ[i].y+rapz*z[i]*NQ[i].y,PQ[i].z+rapz*z[i]*NQ[i].z);
                    else
                        P[i]= R3(PQ[i].x,PQ[i].y,PQ[i].z);
                else if(viewdim==2)
                    P[i]= R3(PQ[i].x,PQ[i].y,rapz*z[i]);
      
            glLineWidth(3*kscreenscale);
            glPolygonMode(GL_FRONT,GL_LINE);
           
            if(ffbatch!=2) glBegin(GL_LINES);
            for (int i=0;i<im;i++) {
                if((debug > 100)) cout << i << " \t : " << Pt[i].x << " " <<  Pt[i].y << " " << Pt[i].z << endl;
                glVertex3f(P[i].x, P[i].y, P[i].z);
            }
            if(ffbatch!=2) glEnd();
            glLineWidth(kscreenscale);
        }
    }
}


bool WindowDump (int width, int height) {
  int i, j;
  FILE *fptr;
  static int counter = 0;
  char fname[32];
  unsigned char *image;

  /* Allocate our buffer for the image */
  if ((image = new unsigned char[3*width*height]) == NULL) {
    fprintf(stderr, "WindowDump - Failed to allocate memory for image\n");
    return false;
  }

  /* Open the file */
  snprintf(fname,32, "ffglut_%04d.ppm", counter);
  if ((fptr = fopen(fname, MODE_WRITE_BINARY)) == NULL) {
    fprintf(stderr, "WindowDump - Failed to open file for window dump\n");
    delete[] image;
    return false;
  }
  if((debug > 10)) cout << " WindowDump in " << fname << endl;

  /* Copy the image into our buffer */
  glReadBuffer(GL_FRONT);
  glPixelStorei(GL_PACK_ALIGNMENT,1);// no padding at the end of the lines
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);

  /* Write the PPM file */
  fprintf(fptr, "P6\n%d %d\n255\n", width, height);
  for (j = height-1; j >= 0; j--) {
    for (i = 0; i < width; i++) {
      fputc(image[3*j*width+3*i+0], fptr);
      fputc(image[3*j*width+3*i+1], fptr);
      fputc(image[3*j*width+3*i+2], fptr);
    }
  }
  fclose(fptr);

  delete[] image;
  counter++;
  return true;
}


static void PNGChunk(FILE *f,const char *type,const vector<unsigned char> & data)
{ //  a chunk of a png file: length, type, data, crc of type and data
    static unsigned int crctable[256];
    static bool init=false;
    if(!init)
    {
        for(unsigned int n=0;n<256;++n)
        {
            unsigned int c=n;
            for(int k=0;k<8;++k)
                c = (c&1) ? 0xedb88320U ^ (c>>1) : c>>1;
            crctable[n]=c;
        }
        init=true;
    }
    unsigned int len=data.size(), crc=0xffffffffU;
    unsigned char b[4]={(unsigned char)(len>>24),(unsigned char)(len>>16),(unsigned char)(len>>8),(unsigned char) len};
    fwrite(b,1,4,f);
    fwrite(type,1,4,f);
    for(int i=0;i<4;++i)
        crc = crctable[(crc^(unsigned char) type[i])&0xff]^(crc>>8);
    for(size_t i=0;i<data.size();++i)
        crc = crctable[(crc^data[i])&0xff]^(crc>>8);
    crc ^= 0xffffffffU;
    if(len) fwrite(&data[0],1,len,f);
    unsigned char e[4]={(unsigned char)(crc>>24),(unsigned char)(crc>>16),(unsigned char)(crc>>8),(unsigned char) crc};
    fwrite(e,1,4,f);
}

bool WindowDumpPNG (int width, int height)
{ //  window dump in a png file, without library: the data are not compressed (stored blocks of zlib)
    static int counter = 0;
    char fname[32];
    snprintf(fname,32, "ffglut_%04d.png", counter);
    FILE *f=fopen(fname, MODE_WRITE_BINARY);
    if(!f)
    {
        fprintf(stderr, "WindowDumpPNG - Failed to open file for window dump\n");
        return false;
    }
    vector<unsigned char> image(3*width*height);
    glReadBuffer(GL_FRONT);
    glPixelStorei(GL_PACK_ALIGNMENT,1);// no padding at the end of the lines
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &image[0]);
    //  the lines from the top, with the filter byte 0 at the beginning of each line
    vector<unsigned char> raw;
    raw.reserve((3*width+1)*height);
    for (int j = height-1; j >= 0; j--)
    {
        raw.push_back(0);
        raw.insert(raw.end(),image.begin()+3*j*width,image.begin()+3*(j+1)*width);
    }
    //  zlib stream: header, stored blocks ( <= 65535 bytes), adler32
    vector<unsigned char> z;
    z.push_back(0x78);
    z.push_back(0x01);
    unsigned int s1=1,s2=0;
    for(size_t i=0;i<raw.size();i+=65535)
    {
        size_t n=min((size_t) 65535,raw.size()-i);
        z.push_back(i+n==raw.size());// last block
        z.push_back(n&0xff);
        z.push_back(n>>8);
        z.push_back((~n)&0xff);
        z.push_back(((~n)>>8)&0xff);
        z.insert(z.end(),raw.begin()+i,raw.begin()+i+n);
        for(size_t k=i;k<i+n;++k)
        {
            s1=(s1+raw[k])%65521;
            s2=(s2+s1)%65521;
        }
    }
    unsigned int adler=(s2<<16)|s1;
    z.push_back(adler>>24);
    z.push_back((adler>>16)&0xff);
    z.push_back((adler>>8)&0xff);
    z.push_back(adler&0xff);
    static const unsigned char sig[8]={137,80,78,71,13,10,26,10};
    fwrite(sig,1,8,f);
    vector<unsigned char> hd(13,0);
    hd[0]=width>>24, hd[1]=(width>>16)&0xff, hd[2]=(width>>8)&0xff, hd[3]=width&0xff;
    hd[4]=height>>24, hd[5]=(height>>16)&0xff, hd[6]=(height>>8)&0xff, hd[7]=height&0xff;
    hd[8]=8;// 8 bits
    hd[9]=2;// RGB
    PNGChunk(f,"IHDR",hd);
    PNGChunk(f,"IDAT",z);
    PNGChunk(f,"IEND",vector<unsigned char>());
    fclose(f);
    cout << " window dump in " << fname << endl;
    counter++;
    return true;
}

void Clean()
{
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

static void Reshape( int width, int height )
{
     GLint pView[4];
 
    width*= kscreenscale;
    height*= kscreenscale;
  //  cout << " width " << width << " "<<height << endl;
    OneWindow * win=CurrentWin();
    if(win)
    {
        win->resize(width,height);
    //    glGetIntegerv(GL_VIEWPORT, pView);
    //     cout << "pView " << pView[0]<< " "<< pView[1]<< " "<< pView[2]<< " "<< pView[3]<< endl;

    }
    glutPostRedisplay();
}



void Display(void)
{
    OneWindow * win=CurrentWin();
    if(win)
    {
        if (win->countdisplay++<0) { // HAsh for MAC OS Mojave ???? FH oct 2018 ...
            if(debug>4  || 1) cout << "\n **** Hack Mojove :: glutReshapeWindow " << win->width << " " << win->height <<" \n" << endl;
            win->width--;
            win->height--;
#ifdef GLUT_WINDOW_SCALE
            glutscreenscale= glutGet(GLUT_WINDOW_SCALE);
#endif
            glutReshapeWindow(win->width,win->height);
            win->resize(win->width,win->height);
            glutPostRedisplay();
            return;
        }

        if(debug>9)
        {
            cout << "\n\n     Display "<< win->theplot  << "  true wiat: " << (!win->theplot || !win->theplot->wait || gwait) << endl;
        }
        {
            Clean();
            win->Display();
            glFlush();
            glutSwapBuffers();
            if ( win->windowdump)
                WindowDump(win->width,win->height);
            win->windowdump=false;
            if ( win->windowdumppng)
                WindowDumpPNG(win->width,win->height);
            win->windowdumppng=false;
        }

    }
    if(debug>9)
    {
        cout << "  Send  SendForNextPlot : win->theplot " << win->theplot << "   " << (!win->theplot || !win->theplot->wait || gwait);;
        if(win->theplot) cout << "    wait" << win->theplot->wait << endl;
        cout <<"\n\n\n"<< endl;
    }
    if(!win->theplot || !win->theplot->wait || gwait )
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
    x *= kscreenscale;
    y *= kscreenscale;

    OneWindow * win=CurrentWin();
    keyact = glutGetModifiers();
    if(button==3 || button==4)
    { // wheel of the mouse (button 3 and 4): zoom, the point under the cursor does not move
        if(win && win->theplot && state == GLUT_DOWN)
        {
            win->zoomat(x,y,button==3 ? 0.9 : 1./0.9);
            glutPostRedisplay();
        }
        return;
    }
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
    x *= kscreenscale;
    y *= kscreenscale;
    if(win) win->xmouse=x,win->ymouse=y;
    switch(casemouse)
    {
        case GLUT_LEFT_BUTTON:

            if(win && win->Graph1d())
            { // 1d graph: translate
                int x0,y0,x1,y1;
                win->GraphBox(x0,y0,x1,y1);
                R2 D2(-(x-win->xold)*(win->Bmax.x-win->Bmin.x)/(x1-x0),(y-win->yold)*(win->Bmax.y-win->Bmin.y)/(y1-y0));
                win->oBmin += D2;
                win->oBmax += D2;
                win->Bmin += D2;
                win->Bmax += D2;
                win->xold=x;
                win->yold=y;
                glutPostRedisplay();
            }
            else if(win && win->cutaxis && (win->cutedit || ((keyact & GLUT_ACTIVE_SHIFT) && !win->cutfreeze)))
            { // F2 (edit) or shift (not with F3): rotation of the cut plane
                win->RotateCut((x-win->xold)/180.,(y-win->yold)/180.);
                win->xold=x;
                win->yold=y;
                glutPostRedisplay();
            }
            else if(win)
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

static void PassiveMotionMouse(int x,int y )
{
    OneWindow * win=CurrentWin();
    if(win)
    {
        x *= kscreenscale;
        y *= kscreenscale;
        bool move = x != win->xmouse || y != win->ymouse;
        win->xmouse=x;
        win->ymouse=y;
        if(move && (win->Graph1d() || win->probe)) // to show the coordinates / value under the mouse
            glutPostRedisplay();
    }
}

static void Key( unsigned char key, int x, int y )
{
    x *= kscreenscale;
    y *= kscreenscale;
    OneWindow * win=CurrentWin();
    if(!win->theplot)return;
    if(debug>1) cout << "Key winnum:  " <<win->theplot->winnum << endl;
    int ni=win->theplot->Viso.N();
    int na=win->theplot->Varrow.N();
    int nbN=win->theplot->nbN;
    
    bool washelp = win->help;
    win->help= false;
    switch (key)
    {
        case 27: // esc char
            Fin(0);
            break;
        case 3: // esc char
            Fin(9);
            break;
        case 'w':
            if(win)
                win->windowdump=true;
            break;
        case 'l':
            win->withlight = !win->withlight;
            break;
        case '?' : // the help of the current plot, again: all the keys
            if(win)
            {
                win->helpall = washelp && !win->helpall;
                win->help= true;
            }
            break;
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
            win->changeiso=win->changearrow=true;// the lists depend on the dimension of the view
            break;
        case '2':
            win->changePlotdim=win->changePlotdim==0 ? 1 : 0;
            win->changeiso=win->changearrow=true;
            break;
        case '=':
            win->DefaultView(1);
            break;
        case 'f':
            win->theplot->fill = ! win->theplot->fill  ;
            break;
        case '@':
            win->dtheta = win->dtheta ? 0 : pi/1800.;
            break;
        case '*': // add FH  mars 2013 ..
            win->keepPV = ! win->keepPV;
            if(debug) cout << "  ... win->keepPV "<< win->keepPV << endl;
            break;
        case 'k':
            if(win->theplot->NextCase())
            {
                R fmn,fmx,vmn,vmx;
                win->theplot->dyn_bfv(win,fmn,fmx,vmn,vmx) ;
                win->theplot->SetDefIsoV(ni,na,fmn,fmx,vmn,vmx) ;
                if(fmx>=fmn && win->theplot->datadim!=3)
                { //  the new part of the complex function (real, imag, modulus, arg) has other values:
                  //  the z range of the view (the elements are drawn at z = value, outside they were
                  //  clipped by the near and far planes: white areas)
                    win->zmin=fmn;
                    win->zmax=fmx;
                    win->Bmin3.z=fmn;
                    win->Bmax3.z=fmx;
                    win->Pvue3.z=(fmn+fmx)/2;
                }
                win->changeiso=true;
                win->changearrow=true;
                if(win->Graph1d()) win->DefaultView(1);// new bound of the axis
            }
            else if(win->Graph1d())
                win->graphstyle++;
            break;

        case 'b':
            win->theplot->grey = ! win->theplot->grey   ;
            win->changeiso=true;
            win->changearrow=true;

            break;
        case 'g':
            win->theplot->grey = !  win->theplot->grey   ;
            win->changeiso=true;
            win->changearrow=true;

            break;

        case 'v':
            bool pi,pa;
            if(win->Graph1d())
                win->graphkey = ! win->graphkey;
            else if(win->PrettyOn() && (win->theplot->PrettyPlots(win,pi,pa),pi||pa))
                win->prettybar = ! win->prettybar;
            else
                win->theplot->value = ! win->theplot->value  ;
            break;
        case 'G':
            win->graphgrid = ! win->graphgrid;
            break;
        case 'o': // the iso lines on the colors (pretty rendering 2d with fill)
            win->isolines = ! win->isolines;
            win->changeiso=true;
            break;
        case 'P': // the value under the mouse
            win->probe = ! win->probe;
            break;
        case 'j': // the axes x, y (and z) with the tics and the values
            win->axes = ! win->axes;
            win->changearrow=true;// the size of the pixel changes with the margins
            break;
        case 'R': // pretty rendering of the 3d plots
            win->pretty = ! win->pretty;
            win->changeiso=true;
            win->changearrow=true;
            break;
        case 'C': // colormap of the pretty rendering
            win->cmapkind = (win->cmapkind+1)%5;
            win->changeiso=win->changearrow=true;// the iso lines and the arrows have the colors in the lists
            {
                static const char *names[5]={"turbo","viridis","cool to warm","rainbow","grey"};
                cout << " ffglut: colormap " << names[win->cmapkind] << endl;
            }
            break;
        case 's':
            win->isosurf = ! win->isosurf;
            break;
        case 'e': // exploded view of the 3d P0 function
            win->explode = min(win->explode+0.1,3.);
            break;
        case 'E':
            win->explode = max(win->explode-0.1,0.);
            break;
        case 'W':
            if(win)
                win->windowdumppng=true;
            break;
        case '[': // threshold of the 3d P0 function
            win->threshold = max(win->threshold-0.05,0.);
            break;
        case ']':
            win->threshold = min(win->threshold+0.05,1.);
            break;
        case 'x': // 1d graph: log scale in x, otherwise rotation of the cut plane of 5 degree
            if(win->Graph1d())
            {
                win->graphlogx = win->userlogx = ! win->graphlogx;
                win->DefaultView(1);
            }
            else
                win->RotateCut(0.,5*pi/180.);
            break;
        case 'X':
            win->RotateCut(0.,-5*pi/180.);
            break;
        case 'y':
            if(win->Graph1d())
            {
                win->graphlogy = win->userlogy = ! win->graphlogy;
                win->DefaultView(1);
            }
            else
                win->RotateCut(5*pi/180.,0.);
            break;
        case 'Y':
            win->RotateCut(-5*pi/180.,0.);
            break;
        case '<':
        case ',':
            win->cutpos = max(win->cutpos-0.02,0.);
            break;
        case '>':
        case '.':
            win->cutpos = min(win->cutpos+0.02,1.);
            break;
        case 'F':
            graphbold = ! graphbold;
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
        {
            na  -=  na < 10  ? 2 : 5;
            ni  -=  ni <10 ? 2 : 5;
            na = max(na,2);
            ni = max(ni,2);
            R fmn,fmx,vmn,vmx;
            win->theplot->dyn_bfv(win,fmn,fmx,vmn,vmx) ;
            win->theplot->SetDefIsoV(ni,na,fmn,fmx,vmn,vmx) ;

            win->changeiso=true;
            win->changearrow=true;
        }
            break ;
        case 'N':
        {
            na  += na < 10  ? 2 : 5;
            ni  += ni < 10  ? 2 : 5;
        }
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
            if(win->focal < my_pi/1.2 )
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
        case 'T':
            win->theplot->pNormalT = ! win->theplot->pNormalT  ;
            break;
        case 'Q':
        {
            win->theplot->nbN  /= win->theplot->nbN < 0  ? 0 : 2;
            break;
        }
        case 'q':
        {
            win->theplot->nbN  *= win->theplot->nbN < 0  ? 0 : 2;
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
        R dx(0),dy(0);
        //  the cut plane as in medit: F1 on/off, F2 edit with the mouse, F3 freeze on the screen
        if(key==GLUT_KEY_F1 || key==GLUT_KEY_F2 || key==GLUT_KEY_F3)
        {
            win->help=false;
            if(key==GLUT_KEY_F1)
            {
                win->cutaxis = win->cutaxis ? 0 : 1;
                win->cutN = R3(1.,0.,0.);
                win->cutedit = win->cutfreeze = false;
            }
            else if(win->cutaxis && key==GLUT_KEY_F2)
            {
                win->cutedit = !win->cutedit;
                win->cutfreeze = false;
            }
            else if(win->cutaxis && key==GLUT_KEY_F3)
            { //  freeze: the plane is fixed on the object, the mouse turns the object with the plane
                win->cutfreeze = !win->cutfreeze;
                win->cutedit = false;
            }
            glutPostRedisplay();
            return;
        }
        switch (key) {
            case  GLUT_KEY_LEFT:   dx = -1; break;
            case  GLUT_KEY_RIGHT:  dx = +1; break;
            case  GLUT_KEY_DOWN:   dy = -1; break;
            case  GLUT_KEY_UP:     dy = +1; break;
        }
        // calcul du deplacement de xm,ym,zm;
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

THREADFUNC(ThreadRead,fd)
{
    int err=0;
    assert(nextPlot==0);
    if(gwait) usleep((useconds_t) (gwait*1.e6) );
    err=ReadOnePlot((FILE*)fd);
    if(debug>1)
        cout << " We Read a plot  : " << kread << " " << nextPlot << " " << err << endl;
    if(err<0)
    {
        NoMorePlot=true;

    }
    Thread::Exit();
    return  0;
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
        int iwnp= nextPlot-> winnum;
        if(debug>1) cout << " change current plot to: " << nextPlot << " et  Lock Plot . winnum  " << iwnp << endl;;
        if(Num2Windows[iwnp]==0)
        {
            ostringstream titre;
            titre <<  "W "<< iwnp << " /FreeFem++: type return key to proceed (or ? for help on other)";
            int Height = 512;
            int Width = 512*3/2;

            glutInitWindowSize(Width , Height);
            glutInitWindowPosition(100+iwnp*50, 100+iwnp*50);
            int iw0=glutCreateWindow(titre.str().c_str());
            if(debug>1)  cout << " ** glutCreateWindow  " << iw0 << endl;
            glDisable(GL_DEPTH_TEST);
            glutReshapeFunc( Reshape ); // pour changement de fenetre
            glutKeyboardFunc( Key );    // pour les evenements clavier
            glutSpecialFunc(SpecialKey);
            glutMouseFunc(Mouse);       // pour les evenements sourie
            glutMotionFunc(MotionMouse); // les mouvements  de la sourie
            glutPassiveMotionFunc(PassiveMotionMouse);
            glutDisplayFunc( Display ); // l'affichage
            glutPushWindow();
            glutShowWindow();
            AllWindows[iw0]= new OneWindow(Width+1 , Height+1,nextPlot);
            Num2Windows[iwnp]=iw0;
        }
        else
            AllWindows[Num2Windows[iwnp]]->add(nextPlot);
        glutSetWindow(Num2Windows[iwnp]);
        currentPlot=nextPlot;
        nextPlot=0;
        LauchNextRead();
        ret=true;
    }
    if(gwait &&  NoMorePlot )
    {usleep((useconds_t)(1e6*gwait)); Fin(0); }

    return ret;
}
const char * Index(const char * p, const char c)
{
    int k=0;
    while(k++<1000000)
        if(!*p) return 0;
        else if(  *p==c) return p;
        else ++p;
    return 0;
}
const char * rIndex(const char * p, const char c)
{
    int k=0;
    const char *pp=0;
    while(k++<1000000)
        if(!*p) break;
        else if(  *p==c) pp= p;
        else ++p;
    return pp;
}
void     SetDefWin(const char *p,int & iii0,int & jjj0,int & Width,int &Height)
{
    // syntax
    //   1024x1024+100+100
    // or
    //  1024x1024
    const char  *bx = p;
    const char *by = Index(p,'x');
    const char *ox = Index(p,'+');
    const char *oy = rIndex(p,'+');
    if(by ==0) return;
    Width= atoi(bx);
    Height= atoi(by+1);
    if(ox && (ox != oy))
    {
        iii0= atoi(ox+1);
        jjj0=atoi(oy+1);
    }
    if(debug>1)
        cout << "  position = "<< Width << "x" << Height << "+"<< iii0 << "+" << jjj0 << endl;
    assert(Width >0 &&  Width < 3000);
    assert(Height >0 &&  Height < 3000);
    assert(iii0 >0 &&  iii0 < 3000);
    assert(jjj0 >0 &&  jjj0 < 3000);

}
int main(int argc,  char** argv)
{
    lockOrientation=0;// to get bad mesh !! FH mai 2021 
    ffapi::init();
    glutInit(&argc, argv);
    bool stereo=false;
    bool fullscreen = false;

    if(stereo)
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO);
    else
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);// antialiasing
    int i1=1;
    int W = glutGet(GLUT_SCREEN_WIDTH);
    int H = glutGet(GLUT_SCREEN_HEIGHT);
    glutscreenscale = 1;
 #ifdef GLUT_WINDOW_SCALE
    glutscreenscale= glutGet(GLUT_WINDOW_SCALE);
    if(glutscreenscale < 0) // unsupported by original GLUT
        glutscreenscale = 1;
#endif
    int Height = glutscreenscale*2*H/3,Width = glutscreenscale*W*2/3, iii0=100,jjj0=100;
    
    string titre = "W0/FreeFem++: type return key to proceed (or ? for help on other)";
    int eerr=0;
    if(argc>1)
    {
        if( (i1 < argc) &&strcmp(argv[i1],"-?")==0) i1++,eerr=-1;
        if (i1 < argc)
        {
            if( strcmp(argv[i1],"-nv")==0) i1++,debug=0;
            else if( strcmp(argv[i1],"-v")==0) i1++,debug=2,verbosity=2;
            else if( strcmp(argv[i1],"-vv")==0) i1++,debug=5,verbosity=2;
            else if( strcmp(argv[i1],"-vvv")==0) i1++,debug=10, verbosity=1000;
            else if( strcmp(argv[i1],"-vvvv")==0) i1++,debug=99999, verbosity=1000;
            else if( strcmp(argv[i1],"-0")==0) return 0; // do nothing
        }
        if( (i1+1 < argc) && (strcmp(argv[i1],"-wait")==0)) { i1++; gwait=atof(argv[i1++]); }
        if( (i1+1 < argc) && (strcmp(argv[i1],"-g")==0)) {
            i1++;
            SetDefWin(argv[i1++], iii0,jjj0,Width,Height);
        }
        if( (i1+1 < argc) && (strcmp(argv[i1],"-t")==0)) {
            i1++;
            titre = argv[i1++];
        }
        if( (i1 -  argc > 1) )
        {
            eerr = 1;
            cout << " error ming args " << i1 -  argc  << endl;
        }
    }
    if(debug>1)
    {
        cout << " glutscreenscale " << glutscreenscale << " " << Width << " " << Height <<  " " << W << " " << H << endl;
    }
    datafile =0;
    if(debug>1)
        cout << "ffglut : " << (argc>i1) << eerr << endl;
    if(argc>i1 && (eerr==0))
    {
        datafile=fopen(argv[argc-1], "r");
        if(debug >1)
            cout << " fopen :" << argv[argc-1] << " " <<datafile << endl;
        if(datafile==0)
            eerr=100;

    }
    if(eerr)
    {
        cerr << " Erreur ffglut  [-nv|-v|-vv|-vvv] [-wait 0.5] [-g 512x300+10+10] [-t title] [file]" << endl;
        cerr << " err number " << eerr << endl;
        abort();
    }


    if(debug>1)
        cout <<  " mode read = " << MODE_READ_BINARY << endl;
    if(datafile==0)
        datafile=stdin;
    if ( !datafile){
        cerr<< " Erreur fdopen stdin in binary " << endl;
        Fin(1);
    }
    int err=ReadOnePlot(datafile);
    if(err) {
        Fin(1);}



    if(kread==0) {
        Fin(1);
    }
    if(debug>1)
        cout << "on a lue le premier plot next plot: " << nextPlot << endl;

    
  
    glutInitWindowSize((Width+1) , (Height+1));
    glutInitWindowPosition(iii0*glutscreenscale,jjj0*glutscreenscale);

    int iw0=glutCreateWindow(titre.c_str());
    GLint pView[4];
    glGetIntegerv(GL_VIEWPORT, pView);
 // hack to correct bug in Catalina , incompatibility Width and pView (glutscreenscale???)
    kscreenscale = (double) ( pView[2]- pView[0])/ (double) (Width+1) ;
    if (kscreenscale < 0.1) kscreenscale= 1;
    if(kscreenscale > 10.) kscreenscale=1;
    if( debug>1 )
    cout << "pView " << pView[0]<< " "<< pView[1]<< " "<< pView[2]<< " "<< pView[3]<< " kscreenscale=  " << kscreenscale << endl;
    
    Num2Windows[0]=iw0;
    glDisable(GL_DEPTH_TEST);
    glutReshapeFunc( Reshape ); // pour changement de fenetre
    glutKeyboardFunc( Key );    // pour les evenements clavier
    glutSpecialFunc(SpecialKey);
    glutMouseFunc(Mouse);       // pour les evenements sourie
    glutMotionFunc(MotionMouse); // les mouvements  de la sourie
    glutPassiveMotionFunc(PassiveMotionMouse);
    glutDisplayFunc( Display ); // l'affichage
    glutSetWindow(iw0);
    AllWindows[iw0]=new OneWindow(Width+1 , Height+1 ,currentPlot);
    TryNewPlot();
    glutMainLoop();

    return 0;
}
