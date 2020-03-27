#ifndef _DataFindBoundary_hpp__
#define _DataFindBoundary_hpp__

#include <RNM.hpp>
#include <fstream>

template<typename Mesh>
struct GenericDataFindBoundary
{
    typedef typename Mesh::Vertex Vertex;
    
    typedef typename Mesh::Element Element;
    typedef typename Mesh::BorderElement BorderElement;
    typedef  typename Mesh::Rd Rd;
    typedef  typename Mesh::RdHat RdHat;
    //using   EF23::GTree;

    static const int d = Rd::d;
    static const int dHat = RdHat::d;
    static const bool bborder = d == dHat; //  build border .??
    const Mesh *pTh;
    EF23::GTree<Vertex> *tree;
    mutable long nbfind, nbelement;
    KN<Vertex> P;// the barycentre of element
    KN<double> delta;// de dist of recheche
    KN<Vertex*> lp;// buffer of recheche ..
    long debug ;
    GenericDataFindBoundary(Mesh const * pTh,int ddebug=0);
    ~GenericDataFindBoundary() ;//{delete tree;}
    int Find(Rd P,double *l,int & outside) const ;
    void gnuplot(const string & fn);
};


#endif

