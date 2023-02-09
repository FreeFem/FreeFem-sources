// -*- C++ -*-
// Time-stamp: "2022-10-04 18:04:46 fujiwara"
//
// SUMMARY : Plot 2D Numerical Results, Mesh, and Indices into PDF file
// ORG     : Graduate School of Informatics, Kyoto University, Japan
// AUTHOR  : Hiroshi Fujiwara
// E-MAIL  : fujiwara@acs.i.kyoto-u.ac.jp
// 
// The newest version is avalilable at:
// http://www-an.acs.i.kyoto-u.ac.jp/~fujiwara/ff
//
//ff-c++-LIBRARY-dep: [zlib]
//----------------------------------------------------------------------
// Usage:
//
// (1) bool plotPDF( PDFfilename, mesh Th [, options] );
//     
//     Drawing 'mesh Th' only.
//
// (2) bool plotPDF( PDFfilename, mesh Th, Vh u [, options] );
// 
//     Drawing both 'mesh Th' and 'FE function u'.
//
// (3) bool plotPDF( PDFfilename, mesh Th, Vh<complex> u [, options] );
// 
//     Drawing both 'mesh Th' and 'FE<complex> function u'.
//
// (4) bool plotPDF( PDFfilename, mesh Th, [Vh u, Vh v] [, options] );
//
//     Drawing both 'mesh Th' and 'FE vector field [u,v]'
//
// Options and Default Values:
//
//  --------------------
//  page control options
//  --------------------
//   int size      = 512   : figure pane sizes in pixel
//   real ar       = auto  : aspect ratio, x/y
//
//   bool meshpage = true  : true if show mesh page
//   bool index    = false : true if show mesh with index page
//   bool belabel  = false : true if show boundary edge labels defined
//                         : by "border ... { label=...; }" in .edp
//   bool isoline  = true  : true if show isoline of profile
//   bool fill     = true  : true if show profile with fill-style
//
//  --------------------------
//  plot-style control options
//  --------------------------
//   bool gray        = false : true if monochrome (black white)
//   bool bw          = false : equivalent to gray
//   bool value       = true  : true if show legend in isoline and fill pages
//   bool logscale    = false : true if plot functions in logarithmic scale
//   real withmesh    = 0.0   : thickness over-layed mesh, 0=white, 1=black
//   real fontscale   = 1.0   : ratio of font size for mesh index and belabel
//   real[int] viso   = []    : values to plot isolines
//   int nbiso        = 12    : number of isolines
//   real[int] frange = []    : value range in plotting with fill-style
//   int nbfill       = 32    : number of colors (levels) in fill-style
//   string fetype    = "P1"  : finite element
//                            : P0,P1,P2,P1nc(Crouzeix-Raviart) are available
//                            : If not specified, piecewise linear interpolated
//                            : shape is drawn
//   string title             : title of the output PDF file
//   string cmm               : comment shown on the graph
//   real[int] fmargin = [0,0,0,0] : margin (left, bottom, right, and top)
//   int prec         = 3     : number of digits in legends
//   real lw          = 1.0   : line width of contours and arrows
//   real[int,int] palette    : color spec. for functions
//                              (default : Color Universal Design ver.4.)
//      palette = [ [r0,g0,b0,], [r1,g1,b1], ..., [rN-1,gN-1,bN-1] ];
//
//      N (>= 2) triplets of integers (0 <= r,g,b <= 255).
//      [r0,g0,b0,] => fmin, [r_{N-1},g_{N-1},b_{N-1}] => fmax
//      Inteval between fmin and fmax is divided into N-1 subintervals.
//      RGB is interpolated in each sub-interval.
//
// If both viso and nbiso are specified, viso is prior to nbiso.
// If fetype="P0", then isoline is ignored.
//
//  -------------------------
//  index page
//  -------------------------
//   bool idcell   = true  : true if show cell (triangle) ID
//   bool idvert   = true  : true if show vertex ID
//   bool idedge   = false : true if show edge ID
//
//  -------------------------
//  vector field (vector-valued function)
//  -------------------------
//   real ArrowSize   = 1.0 : scale for arrow head size
//   real coef        = 1.0 : scale for arrow length
//   bool unitarrow   = false : true if draw arrows with same length
//
//  * fill is ignored.
//
//  -------------------------
//  complex-valued function
//  -------------------------
//   bool zreal = true  : true if show real part
//   bool zimag = true  : true if show imaginary part
//   bool zabs  = false : true if show modulus (absolue value)
//   bool zarg  = false : true if show argument
//
//  * fetype is ignored (isoline, fill are interpolated as P1 type).
//----------------------------------------------------------------------
// Examples:
//
// mesh Th = ...;
// fespace Vh(Th,P1);   // Only P1-FE is accepted by this module
// Vh u;
//
// load "plotPDF"
// bool ret = plotPDF( "sample_mesh", Th );  // Save to "sample_mesh.pdf"
//                                           // File extension ".pdf" is
//                                           // automatically appended.
// assert( ret == true );  // Error check.
//
// plotPDF("sample_u.pdf", Th, u ); // save solution and mesh
//
//
// plotPDF( "sample", Th, u, meshpage=false,belabel=false,fill=false); // save only isoline
// plotPDF( "sample", Th,    index=true, fontscale=0.7); // include mesh,label, and index
//                                                       // small font
// plotPDF( "sample", Th, u, value=false ); // value is omitted in solution pages
//
//  
// real[int] visoarray=[-0.02,-0.01,0,0.01,0.02];    // specify isoline-values
// plotPDF("sample_with_viso",Th,u, viso=visoarray);
//
// plotPDF("sample_with_viso",Th,u, nbiso=10); // show 10 isolines, values are
//                                            // determined, equi-spaced
//                                            // between min(u) and max(u).
// real[int] range=[-0.02,0.02];  // specify range for plot (min and max)
// plotPDF("sample_with_range",Th,u, frange=range);
//
// plotPDF("sample_fill16",Th,u, nbfill=16);
// plotPDF("sample_fill64",Th,u, nbfill=64);
//
// Vh v;
// plotPDF("vector_field",Th,[u,v], isoline=false, withmesh=0.25); // draw vector field [u,v]
// plotPDF("vector_field",Th,[u,v], isoline=false, withmesh=0.25,coef=2.0,ArrowSize=0.8);
//
// fespace Vh0(Th,P0);
// Vh0 u0;
//
// plotPDF( "sampleP0", Th, u0, fetype="P0"); // plot u0 as P0-FE
// plotPDF( "sampleP0", Th, u0, title="sample P0"); // save with document title
//
// Vh0 uasP0 = u; // interpolation inside FreeFEM
// plotPDF( "cast_fromP1", Th, uasP0,  );
// plotPDF( "asP0-FE",     Th, u, fetype="P0" );
// plotPDF( "solution",    Th, u, cmm="solution" );
//----------------------------------------------------------------------
// TODO:
//  (1) fetype : RT0(Raviart-Thomas), P1b, P2b, P1dc, P2dc
//  (2) append to existing PDF file (e.g., for time evolutional problem)
//  (3) brush up source using pair in <utility> and tuple in <tuple>
//  (4) Vector field options:
//        int nbarrow      = 0   : number of colors (0 : unlimited)
//        real[int] varrow = []  : value of colors
//  (5) Symbol font
//  (6) paper size, axes
//  (7) Non-triangle element
//----------------------------------------------------------------------

#ifndef WITH_NO_INIT
#include "ff++.hpp"
#include "AFunction_ext.hpp"
#endif

#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <queue>

//----------------------------------------------------------------------

namespace
{
    const long   DEFAULT_PAGESIZE    = 512;
    const double DEFAULT_ASPECTRATIO = 0; // "0" means auto
    const char   AppName[] = "FreeFEM plotPDF module";

    const int PLOTPDF_NOPTIONS = 34;

    const double PADDING = 20;
    const int DEFAULT_MARGIN[4] = { 0, 0, 0, 0 }; // left, bottom, right, top

    const double DEFAULT_INDEX_FONTSIZE  = 16;
    const double DEFAULT_LEGEND_FONTSIZE = 12;
    const double DEFAULT_COMMENT_FONTSIZE = 16;
    const double DEFAULT_COMMENT_BASE = DEFAULT_INDEX_FONTSIZE*2; // avoid overlap

    const int  NUM_LABELS = 12; // number of labels (values) in fill-style
    const long DEFAULT_ISOLINES    = 12;
    const long DEFAULT_FILL_COLORS = 32;
    const double DEFAULT_LINEWIDTH = 1;

    const bool DEFAULT_MONOCHROME    = false;
    const bool DEFAULT_SHOW_LEGEND   = true;
    const bool DEFAULT_LOGSCALE      = false;
    const double DEFAULT_WITHMESH    = 0; // gray scale, 0=white, 1=black

    // Complex-valued function
    const bool DEFAULT_ZREAL = true; // draw real part
    const bool DEFAULT_ZIMAG = true; // draw imaginary part
    const bool DEFAULT_ZABS  = false;  // draw modulus of complex-valued function
    const bool DEFAULT_ZARG  = false; // draw argument

    // Vector field
    const double DEFAULT_ARROW_SCALE = 1.0;
    const double DEFAULT_AHEAD_SCALE = 1.0;
    const bool   DEFAULT_UNIT_ARROW = false; // draw arrow with its length
    const long   DEFAULT_ARROW_COLORS  = 0; // 0 : unlimited
    const double MAX_ARROW_LENGTH     = 50;   // [pixel]
    const double DEF_ARROW_HEAD_SIZE  = 8;    // [pixel]
    const double ARROW_HEAD_ANGLE     = 0.23; // 13/180*PI = 0.226893 [rad]

    const bool DEFAULT_SHOW_MESH    = true;
    const bool DEFAULT_SHOW_INDEX   = false;
    const bool DEFAULT_SHOW_BELABEL = false;
    const bool DEFAULT_SHOW_ISOLINE = true;
    const bool DEFAULT_SHOW_FILL    = true;

    const bool DEFAULT_SHOW_IDCELL = true;
    const bool DEFAULT_SHOW_IDVERT = true;
    const bool DEFAULT_SHOW_IDEDGE = false;
    const float DEFAULT_IDCELL_MONO[] = { 0, 0, 0 }; // RGB : black
    const float DEFAULT_IDCELL_RGB[] = { 0, 0, 0 }; // RGB : black
    const float DEFAULT_IDVERT_MONO[] = { 0.5, 0.5, 0.5 }; // RGB
    const float DEFAULT_IDVERT_RGB[] = { 0.0/255.0, 90.0/255.0, 255.0/255.0 }; // RGB : blue
    const float DEFAULT_IDEDGE_MONO[] = { 0.8, 0.8, 0.8 }; // RGB
    const float DEFAULT_IDEDGE_RGB[] = { 255.0/255.0, 75.0/255.0, 0.0/25.0 }; // RGB : red

    const long DEFAULT_PRECISION_LEGEND = 3; // #digits in legend
    const double LEGEND_FONTWIDTH = static_cast<double>(PADDING)/DEFAULT_PRECISION_LEGEND;    // depend on font family and size

    const char DEFAULT_FETYPE[] = "P1";
    const int DEFAULT_P2_INERVALS = 5;

    const int DEFAULT_PALETTE_NCOLORS = 5;
    const double DEFAULT_PALETTE[ DEFAULT_PALETTE_NCOLORS ][3] =
    {
	// color universal design (ver.3, Aug. 2013)
	// {   0,  65, 255 }, // RGB at rf=0.09 (blue)
	// {   0, 161, 255 }, // RGB at rf=0.25
	// {  53, 161, 107 }, // RGB at rf=0.50 (green)
	// { 255, 161,   0 }, // RGB at rf=0.75
	// { 255,  40,   0 }  // RGB at rf=1.00 (red)

	// color universal design (ver.4, Apr. 2018)
	{   0,  90, 255 }, // RGB at rf=0.00 (blue)
	{   0, 175, 255 }, // RGB at rf=0.25
	{   3, 175, 122 }, // RGB at rf=0.50 (green)
	{ 255, 175,   0 }, // RGB at rf=0.75
	{ 255,  75,   0 }  // RGB at rf=1.00 (red)

	// if you like vivid color
	// {   0,   0, 255 }, // RGB at rf=0.00 (blue)
	// {   0, 255, 255 }, // RGB at rf=0.25
	// {   0, 255,   0 }, // RGB at rf=0.50 (green)
	// { 255, 255,   0 }, // RGB at rf=0.75
	// { 255,   0,   0 }  // RGB at rf=1.00 (red)
    };
}

//----------------------------------------------------------------------
// Simple PDF class
//----------------------------------------------------------------------
#include <iomanip>
#include <sstream>
#include <list>
#include <ctime>
#include <cstring> // strlen
#include <cstdlib> // exit
// Change F. Hecht 
//#define HAVE_ZLIB   // need zlib : CC src.cpp -lz
//#if !defined(NO_ZLIB) && !defined(DISABLE_ZLIB) && !defined(WITHOUT_ZLIB) 
//#define HAVE_ZLIB
//#endif
#if defined(WITH_zlib)
#define HAVE_ZLIB
#endif

class SimplePDFModule
{
    int byte_offset;
    std::list<int> xref;

    struct OutlineItem {
	int headPageObjectNumber;
	char *label;
    };
    std::list<OutlineItem> outline;

    const std::string filename;
    const std::string DocumentTitle;
    int page_obj_offset, page;

    std::string get_datetime() const
    {
	std::time_t now = std::time( NULL );
	const std::tm* lt = std::localtime(&now);

	std::stringstream s;

	s << std::setfill('0') << std::right // valid for all operations <<
	  << "20" << std::setw(2) << lt->tm_year-100
	  << std::setw(2) << lt->tm_mon+1    // setw() is reset at each ops <<
	  << std::setw(2) << lt->tm_mday
	  << std::setw(2) << lt->tm_hour
	  << std::setw(2) << lt->tm_min
	  << std::setw(2) << lt->tm_sec;

	return s.str();
    }

#ifdef HAVE_ZLIB
    int deflate_compress( char* &buf, const std::string &Stream ) const;
#endif

public:

    SimplePDFModule( const char *const filename, const char *const title = "" );

    ~SimplePDFModule();

    void addPage( const std::stringstream &ContentStream, const int WIDTH, const int HEIGHT, const int *const MARGIN );

    void addBookmark( const char *const BookmarkLabel );
};

SimplePDFModule::SimplePDFModule( const char *const PDFfilename, const char *const title )
    : filename( PDFfilename ), DocumentTitle( title ), page(0)
{
    std::ofstream fout( filename.c_str(), std::ios::binary );

    if( !fout ){
	std::cerr << "plotPDF() : Cannot open the file: " << filename << std::endl;
	return;
    }

    //--------------------------------------------------
    // Header Section
    //--------------------------------------------------
    std::stringstream header;
#if 1
    // match to pdflatex & dvipdfmx in texlive 2019
    header << "%PDF-1.5\n"
           << '%' << char(0xd0) << char(0xd4) << char(0xc5) << char(0xd8) << '\n';
#else
    header << "%PDF-1.7\n"
           << '%' << char(0xe2) << char(0xe3) << char(0xcf) << char(0xd3) << '\n';
#endif

    fout << header.str();
    byte_offset = header.str().length();

    fout.close();

    //--------------------------------------------------
    // Body Section : Objects
    //--------------------------------------------------
    std::list<const std::string*> obj;

    std::stringstream strDocumentInfo;
    strDocumentInfo << "1 0 obj\n"
		    << "<<\n";
    
    if( strlen( DocumentTitle.c_str() ) > 0 )
	strDocumentInfo << "  /Title (" << DocumentTitle << ")\n";

    strDocumentInfo << "  /Creator (" << AppName << ")\n"
		    << "  /CreationDate (D:" << get_datetime() << ")\n"
		    << ">>\n"
		    << "endobj\n";
    const std::string DocumentInfo = strDocumentInfo.str();

    // /PageLayout \in { /SinglePage (default), /OneColumn,
    //                   /TwoColumnLeft, /TwoColumnRight, /TwoPageLeft, /TwoPageRight }
    // /PageMode \in { /UseNone (default), /UseOutlines, /UseThumbs, /FullScreen, /UseOC, /UseAttachements }
    const std::string DocumentCatalog =
	"2 0 obj\n"
	"<<\n"
	"  /Pages 3 0 R\n"
	"  /Type /Catalog\n"
	"  /PageLayout /SinglePage\n"
	"  /PageMode 4 0 R\n"
	"  /Outlines 5 0 R\n"
	">>\n"
	"endobj\n";

    //--------------------------------------------------
    // Object 3 is reserved for PageTree
    // Object 4 is reserved for PageMode
    // Object 5 is reserved for Outlines
    // Object 6 is reserved for Outlines (document title)
    //--------------------------------------------------

    //--------------------------------------------------
    // Roman Fonts
    //--------------------------------------------------
    const std::string FontEH =
	"7 0 obj\n"
	"<<\n"
	"  /BaseFont /Helvetica\n"
	"  /Subtype /Type1\n"
	"  /Type /Font\n"
	">>\n"
	"endobj\n";

    fout.open( filename.c_str(), std::ios::app|std::ios::binary );

    obj.push_back( &DocumentInfo );
    obj.push_back( &DocumentCatalog );
    obj.push_back( &FontEH );

    for(std::list<const std::string*>::const_iterator itr = obj.begin(); itr != obj.end(); itr++){
	fout << **itr;
	xref.push_back( byte_offset );
	byte_offset += (*itr)->length();
    }

    fout.close();

    page_obj_offset = obj.size()+5;

    return;
}

SimplePDFModule::~SimplePDFModule()
{
    std::ofstream fout( filename.c_str(), std::ios::app|std::ios::binary );

    //----------------------------------------
    // PageTree
    //----------------------------------------
    std::stringstream strPageTree;
    strPageTree << "3 0 obj\n"
		<< "<<\n"
		<< "  /Kids [";

    for(int i = 0; i < page; i++)
	strPageTree << ((i == 0)? "": " ") << page_obj_offset+i*2 << " 0 R";

    strPageTree	<< "]\n"
		<< "  /Type /Pages\n"
		<< "  /Count " << page << "\n"
		<< ">>\n"
		<< "endobj\n";

    fout << strPageTree.str();

    std::list<int>::iterator itr = xref.begin(); // obj 1 (DocumentInfo)
    itr++; // obj 2 (DocumentCatalog)
    itr++; // obj 3 (FontEH at present)
    itr = xref.insert( itr, byte_offset ); // byte_offset is inserted before itr (FontEH)
                                           // and itr indicates byte_offset

    byte_offset += strPageTree.str().length();

    //----------------------------------------
    // /PageMode \in { /UseNone (default), /UseOutlines, /UseThumbs, /FullScreen, /UseOC, /UseAttachements }
    //----------------------------------------
    std::stringstream strPageMode;
    strPageMode << "4 0 obj\n"
		<< "  " << (outline.empty()? "/UseNone\n": "/UseOutlines\n")
		<< "endobj\n";

    fout << strPageMode.str();

    itr++; // obj 3 (FontEH at present)
    itr = xref.insert( itr, byte_offset ); // byte_offset is inserted before itr (FontEH)
                                           // and itr indicates byte_offset
    byte_offset += strPageMode.str().length();

    //----------------------------------------
    // Outlines
    //----------------------------------------
    std::stringstream strOutlines;
    strOutlines << "5 0 obj\n"
		<< "<<\n"
		<< "  /Type /Outlines\n"
		<< "  /Count " << 1+outline.size() << "\n"  // 1 : Document Title
		<< "  /First 6 0 R\n"
		<< "  /Last  6 0 R\n"
		<< ">>\n"
		<< "endobj\n";

    fout << strOutlines.str();

    itr++; // obj 3 (FontEH at present)
    itr = xref.insert( itr, byte_offset ); // byte_offset is inserted before itr (FontEH)
                                           // and itr indicates byte_offset
    byte_offset += strOutlines.str().length();

    //----------------------------------------
    // the top of outline
    //----------------------------------------
    strOutlines.str(""); // clear string buffer
    strOutlines.clear( std::stringstream::goodbit ); // clear status
    
    strOutlines << "6 0 obj\n"
		<< "<<\n"
		<< "  /Title (" << DocumentTitle << ")\n"
		<< "  /Parent 5 0 R\n"
		<< "  /Count " << outline.size() << "\n";

    if( !outline.empty() ){

	const int OutlineItemIdFirst = xref.size()+2;
	const int OutlineItemIdLast  = xref.size()+outline.size()+1;

	strOutlines << "  /First " << OutlineItemIdFirst << " 0 R\n"
		    << "  /Last  " << OutlineItemIdLast  << " 0 R\n";
    }

    strOutlines	<< "  /Dest [" << outline.front().headPageObjectNumber << " 0 R /Fit]\n"
		<< ">>\n"
		<< "endobj\n";

    fout << strOutlines.str();

    itr++; // obj 3 (FontEH at present)
    itr = xref.insert( itr, byte_offset ); // byte_offset is inserted before itr (FontEH)
                                           // and itr indicates byte_offset
    byte_offset += strOutlines.str().length();

    //----------------------------------------
    // outline items
    //----------------------------------------
    for(std::list<OutlineItem>::iterator itr = outline.begin(); itr != outline.end(); itr++){

	std::stringstream strOutlineItem;

	const int object_id = xref.size()+1;

	strOutlineItem << object_id << " 0 obj\n"
		       << "<<\n"
		       << "  /Title (" << itr->label << ")\n"
		       << "  /Parent 6 0 R\n";

	if( itr != outline.begin() )
	    strOutlineItem << "  /Prev " << object_id-1 << " 0 R\n";

	std::list<OutlineItem>::const_iterator succ = itr;
	succ++;
	if( succ != outline.end() )
	    strOutlineItem << "  /Next " << object_id+1 << " 0 R\n";

	strOutlineItem << "  /Count 0\n"
		       << "  /Dest [" << itr->headPageObjectNumber << " 0 R /Fit]\n"
		       << ">>\n"
		       << "endobj\n";

	fout << strOutlineItem.str();

	delete [] itr->label;

	xref.push_back( byte_offset );
	byte_offset += strOutlineItem.str().length();
    }

    //----------------------------------------
    // Corss-reference Table Section
    //----------------------------------------
    const int nobjs = xref.size()+1; // including 0th object (0 0 obj)
    
    fout << "xref\n";
    fout << "0 " << nobjs << std::endl;
    fout << "0000000000 65535 f \n"; // 0 0 obj

    for(std::list<int>::const_iterator itr = xref.begin(); itr != xref.end(); itr++)
	fout << std::setw(10) << std::setfill('0') << *itr << " 00000 n \n";

    //----------------------------------------
    // Trailer Section
    //----------------------------------------
    fout << "trailer\n"
	 << "<<\n"
	 << "  /Info 1 0 R\n"
	 << "  /Root 2 0 R\n"
	 << "  /Size " << nobjs << "\n"
	 << ">>\n"
	 << "startxref\n"
	 << byte_offset << std::endl    // starting position of 'xref'
	 << "%%EOF\n";

    fout.close();

    return;
}
//

#ifdef HAVE_ZLIB

#include <zlib.h>

int SimplePDFModule::deflate_compress( char* &outbuf, const std::string &Stream ) const
{
    char *inbuf = new char [ Stream.length()+1 ];

    for(size_t i = 0; i < Stream.length(); i++) // strcpy or strncpy
	inbuf[i] = Stream[i];
    inbuf[ Stream.length() ] = '\0';

    z_stream z;
    z.zalloc = Z_NULL;
    z.zfree  = Z_NULL;
    z.opaque = Z_NULL;

    if( deflateInit(&z, Z_DEFAULT_COMPRESSION) != Z_OK ) {
	std::cerr << "zlib: deflateInit(): " << ((z.msg)? z.msg: "Error") << std::endl;
	std::exit(1);
    }

    outbuf = new char [ Stream.length()+1 ];

    z.next_in   = reinterpret_cast<unsigned char*>(inbuf);
    z.avail_in  = Stream.length();
    z.next_out  = reinterpret_cast<unsigned char*>(outbuf);
    z.avail_out = Stream.length();

    if( deflate(&z, Z_FINISH) != Z_STREAM_END ){
	std::cerr << "zlib: deflate(): " << ((z.msg)? z.msg: "Error") << std::endl;
	std::exit(1);
    }

    if( deflateEnd(&z) != Z_OK ){
	std::cerr << "zlib: deflateEnd(): " << ((z.msg)? z.msg: "Error") << std::endl;
	std::exit(1);
    }

    delete [] inbuf;

    const int compressed_buf_length = Stream.length() - z.avail_out;

    outbuf[compressed_buf_length] = '\n';

    return compressed_buf_length+1;
}

#endif

void SimplePDFModule::addPage( const std::stringstream &ContentStream, const int WIDTH, const int HEIGHT, const int *const MARGIN )
{
    const int &marginl = MARGIN[0]; // left
    const int &marginb = MARGIN[1]; // bottom
    const int &marginr = MARGIN[2]; // right
    const int &margint = MARGIN[3]; // top

    //----------------------------------------
    // PageObject
    //----------------------------------------
    std::stringstream strPageObject;
    
    strPageObject << page_obj_offset+2*page << " 0 obj\n"
		  << "<<\n"
		  << "  /Type /Page\n"
		  << "  /Parent 3 0 R\n"
		  << "  /Resources << /Font << /F1 7 0 R >> >>\n"
		  << "  /MediaBox [0 0 " << WIDTH+marginl+marginr << ' ' << HEIGHT+marginb+margint << "]\n"
		  << "  /Contents " << page_obj_offset+2*page+1 << " 0 R\n"
		  << ">>\n"
		  << "endobj\n";

    const std::string PageObject = strPageObject.str();

    //----------------------------------------
    // pageContent (stream)
    //----------------------------------------
#ifdef HAVE_ZLIB
    char *buf;
    const int compressed_buf_length = deflate_compress( buf, ContentStream.str() );
#endif

    //----------------------------------------
    // pageContent (object)
    //----------------------------------------
    std::stringstream strContent;
    strContent << page_obj_offset+2*page+1 << " 0 obj\n"
	       << "<< /Length ";
#ifdef HAVE_ZLIB
    strContent << compressed_buf_length << " /Filter /FlateDecode";
#else
    strContent << ContentStream.str().length();
#endif
    strContent << " >>\n"
	       << "stream\n";

#ifdef HAVE_ZLIB
    // Following works well even if buf[] has '\0'
    strContent << std::string( buf+0, buf+compressed_buf_length );
    delete [] buf;
#else
    strContent << ContentStream.str();
#endif

    strContent << "endstream\n"
	       << "endobj\n";

    const std::string Content = strContent.str();

    //----------------------------------------
    // Add Page Objects
    //----------------------------------------
    std::ofstream fout( filename.c_str(), std::ios::app|std::ios::binary );

    const std::string*const PageObj[2] = { &PageObject, &Content };

    for(int i = 0; i < 2; i++){
    	xref.push_back( byte_offset );
    	fout << *PageObj[i];
    	byte_offset += PageObj[i]->length();
    }

    fout.close();
    page++;

    return;
}

void SimplePDFModule::addBookmark( const char *const BookmarkLabel )
{
    OutlineItem item;

    item.headPageObjectNumber = page_obj_offset+2*page;

    item.label = new char [ strlen(BookmarkLabel)+1 ];
    strcpy( item.label, BookmarkLabel );
    
    outline.push_back( item );

    return;
}

//----------------------------------------------------------------------
// End of SimplePDFModule
//----------------------------------------------------------------------

void GaussElimination(double *const x, double *const *const a, const int N)
{
    const double EPS = 1e-10;

    // forward elimination
    for(int k = 0; k < N-1; k++){

	{ // begin pivoting
	    int row = k;
	    
	    for(int r = k+1; r < N; r++)
		if(fabs(a[row][k]) < fabs(a[r][k]))
		    row = r;

	    if(fabs(a[row][k]) < EPS){
		std::cerr << "singular matrix : " << row << std::endl;
		exit(1);
	    }

	    // k == row: no pivoting
	    if(k != row){
		for(int j = 0; j <= N; j++){
		    const double tmp = a[k][j];
		    a[k][j] = a[row][j];
		    a[row][j] = tmp;
		}
	    }
	} // end pivoting

	const double d = 1 / a[k][k];
        for(int i = k+1; i < N; i++){
            for(int j = k+1; j <= N; j++)
                a[i][j] -= a[i][k] * a[k][j] *d;
            a[i][k] = 0;
        }
    }

    // backward substitution
    for(int i = N-1; i >= 0; i--){
	for(int j = i+1; j < N; j++)
	    a[i][N] -= a[i][j]*a[j][N];
	a[i][N] /= a[i][i];
    }

    for(int i = 0; i < N; i++)
	x[i] = (fabs(a[i][N]) < EPS)? 0: a[i][N];

    return;
}

std::string get_string(Stack stack, Expression e, const char *const DEFAULT)
{
    const size_t length = 128;
    char *const carg = new char [ length ];

    if( !e )
	strcpy(carg, DEFAULT);
    else
	strncpy(carg, GetAny<string*>((*e)(stack))->c_str(), length);

    return std::string( carg );
}

void setrgbcolor( std::stringstream &st, const double f_,
		  const KNM<double> &palette,
		  const double fmin_, const double fmax_,
		  const bool monochrome, const bool logscale )
{
    const double EPS = 1e-3;
    const double THRESHOLD = -1e10;

    if( (logscale && (fmin_ <= 0)) || (logscale && (f_ <= 0)) )
	cout << "plotPDF(): logscale for non-positive values.\n";

    const double f    = (logscale)? log(fabs(f_)): f_;
    const double fmin = (logscale)? log(fabs(fmin_)): fmin_;
    const double fmax = (logscale)? log(fabs(fmax_)): fmax_;

    if( fabs(fmax-fmin) < EPS ){
	st << 0.5 << ' ' << 0.5 << ' ' << 0.5 << ' '; // middle
	return;
    }

    double rf = (f - fmin) / (fmax - fmin);

    if( (1 <= rf) && (rf <= 1+EPS) ){
	rf = 1;
    } else if( (-EPS <= rf) && (rf <= 0) ){
	rf = 0;
    } else if( (rf > 1+EPS) || (rf < -EPS) ){
	st << 1 << ' ' << 1 << ' ' << 1 << ' '; // white
	return;
    }
		
    double r(palette(0,0)), g(palette(0,1)), b(palette(0,2));

    if( palette.N() == 1 ){

	if(monochrome)
	    r = g = b = rf;

	if(f <= THRESHOLD)
	    r = b = g = 0;

	st << r << ' ' << g << ' ' << b << ' ';
	return;
    }

    const double step = static_cast<double>(1) / (palette.N()-1);

    for(int i = 0; i < palette.N()-1; i++){

	const double l = i*step;
	const double u = (i == palette.N()-1)? 1: (i+1)*step;

	if( rf > u )
	    continue;

	// linear interpolation
	const double p = (rf - l)/step;

	r = (1-p)*palette(i,0) + p*palette(i+1,0);
	g = (1-p)*palette(i,1) + p*palette(i+1,1);
	b = (1-p)*palette(i,2) + p*palette(i+1,2);

	break;
    }

    if(monochrome)
	r = g = b = rf;

    if(f <= THRESHOLD)
	r = b = g = 0;

    st << r << ' ' << g << ' ' << b << ' ';

    return;
}

//----------------------------------------------------------------------
void addComment( std::stringstream &Content,
		 const double posy, const int marginl, const int marginb,
		 const double fontscale, const string &comment )
{
    std::stringstream &st = Content;

    st << "q\n";
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";
    //st << "1 w\n"; // setlinewidth
    st << "0 0 0 RG\n"; // black

    const int posx = 0;

    st << "BT\n";
    st << "/F1 " << DEFAULT_COMMENT_FONTSIZE*fontscale << " Tf\n";
    st << "1 0 0 1 " << posx << ' ' << posy+DEFAULT_COMMENT_BASE << " Tm "
       << "(" << comment << ") Tj\n";
    st << "ET\n";
    st << "Q\n";

    return;
}

//----------------------------------------------------------------------

void plot_mesh( std::stringstream &Content, const Fem2D::Mesh &Th,
		const double scale, const double ar, const double x0, const double y0,
		const int marginl, const int marginb,
		const double textfontsize, const bool monochrome,
		const double withmesh, const double linewidth,
		const bool idcell, const bool idvert, const int idedge,
		const int mode = 0)
{
    // mode == 0 (default) : without index
    // mode == 1 : FreeFEM index
    // mode == 2 : Boundary Label
    enum { MODE_NOLABEL = 0, MODE_INDEX = 1, MODE_BELABEL = 2 };

    const int nVertices = Th.nv;
    const int nElements = Th.nt;
    const int nEdges    = Th.neb;
    const double &r     = scale;

    std::stringstream &st = Content;
    st.str("");

    st << "q\n";
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";
    st << linewidth << " w\n"; // setlinewidth
    
    //------------------------------
    // vertices
    //------------------------------
    if( (mode == MODE_INDEX) && idvert ){

	if(monochrome)
	    st << DEFAULT_IDVERT_MONO[0] << ' ' << DEFAULT_IDVERT_MONO[1] << ' '
	       << DEFAULT_IDVERT_MONO[2] << " rg\n" << std::endl;
	else
	    st << DEFAULT_IDVERT_RGB[0] << ' ' << DEFAULT_IDVERT_RGB[1] << ' '
	       << DEFAULT_IDVERT_RGB[2] << " rg\n" << std::endl;

	int id = mode;

	st << "BT\n";
	st << "/F1 " << textfontsize << " Tf\n";
	for(int n = 0; n < nVertices; n++){

	    st << "1 0 0 1 " << r*ar*(Th(n).x-x0) << ' ' << r*(Th(n).y-y0) << " Tm "
	       << "(" << id << ") Tj\n";
	    id++;
	}
	st << "ET\n";
    }

    //------------------------------
    // element (triangle)
    //------------------------------
    if( mode == MODE_BELABEL ){

	const double grayscale0 = (withmesh < 0)? 0: withmesh;
	const double grayscale1 = (withmesh < 1)? withmesh: 1;
	const double grayscale = 1-grayscale1;

	st << grayscale << ' ' << grayscale << ' ' << grayscale << " RG\n";

    } else {
	st << "0 0 0 RG\n"; // black
    }
    for(int n = 0; n < nElements; n++){

	const int &v0 = Th(n,0);
	const int &v1 = Th(n,1);
	const int &v2 = Th(n,2);

	st << r*ar*(Th(v0).x-x0) << ' ' << r*(Th(v0).y-y0) << " m ";
	st << r*ar*(Th(v1).x-x0) << ' ' << r*(Th(v1).y-y0) << " l ";
	st << r*ar*(Th(v2).x-x0) << ' ' << r*(Th(v2).y-y0) << " l ";
	st << "s" << std::endl;
    }

    if( (mode == MODE_INDEX) && idcell ){
	
	if(monochrome)
	    st << DEFAULT_IDCELL_MONO[0] << ' ' << DEFAULT_IDCELL_MONO[1] << ' '
	       << DEFAULT_IDCELL_MONO[2] << " rg\n" << std::endl;
	else
	    st << DEFAULT_IDCELL_RGB[0] << ' ' << DEFAULT_IDCELL_RGB[1] << ' '
	       << DEFAULT_IDCELL_RGB[2] << " rg\n" << std::endl;

	int id = mode;
	st << "BT\n";
	st << "/F1 " << textfontsize << " Tf\n";

	for(int n = 0; n < nElements; n++){

	    const Fem2D::Triangle &t = Th[n];
	    const Fem2D::R2 g = ( t[0] + t[1] + t[2] ) / 3;
	    
	    st << "1 0 0 1 " << r*ar*(g.x-x0) << ' ' << r*(g.y-y0) << " Tm "
	       << "(" << id << ") Tj\n"; // show
	    id++;
	}
	st << "ET\n";
    }

    //------------------------------
    // sides (edges) on the boundary
    //------------------------------
    if( (mode == MODE_INDEX) || (mode == MODE_BELABEL) ){

	if(monochrome)
	    st << DEFAULT_IDEDGE_MONO[0] << ' ' << DEFAULT_IDEDGE_MONO[1] << ' '
	       << DEFAULT_IDEDGE_MONO[2] << " RG\n" << std::endl;
	else
	    st << DEFAULT_IDEDGE_RGB[0] << ' ' << DEFAULT_IDEDGE_RGB[1] << ' '
	       << DEFAULT_IDEDGE_RGB[2] << " RG\n" << std::endl;
    }

    for(int n = 0; n < nEdges; n++){

	const int &v0 = Th( Th.bedges[n][0] );
	const int &v1 = Th( Th.bedges[n][1] );

	st << r*ar*(Th(v0).x-x0) << ' ' << r*(Th(v0).y-y0) << " m "
	   << r*ar*(Th(v1).x-x0) << ' ' << r*(Th(v1).y-y0) << " l S\n";
    }

    if( ((mode == MODE_INDEX) && idedge) || (mode == MODE_BELABEL) ){

	if(monochrome)
	    st << DEFAULT_IDEDGE_MONO[0] << ' ' << DEFAULT_IDEDGE_MONO[1] << ' '
	       << DEFAULT_IDEDGE_MONO[2] << " rg\n" << std::endl;
	else
	    st << DEFAULT_IDEDGE_RGB[0] << ' ' << DEFAULT_IDEDGE_RGB[1] << ' '
	       << DEFAULT_IDEDGE_RGB[2] << " rg\n" << std::endl;

	int id = mode;

	st << "BT\n";
	st << "/F1 " << textfontsize << " Tf\n";

	for(int n = 0; n < nEdges; n++){

	    const int &v0 = Th( Th.bedges[n][0] );
	    const int &v1 = Th( Th.bedges[n][1] );

	    const double mx = ( Th(v0).x + Th(v1).x ) / 2;
	    const double my = ( Th(v0).y + Th(v1).y ) / 2;

	    if(mode == MODE_INDEX){

		st << "1 0 0 1 " << r*ar*(mx-x0) << " " << r*(my-y0) << " Tm "
		   << "(" << id << ") Tj\n";

	    } else if(mode == MODE_BELABEL){

		st << "1 0 0 1 " << r*ar*(mx-x0) << " " << r*(my-y0) << " Tm "
		   << "(" << Th.be(n) << ") Tj\n"; // same as Th.bedges[n]
	    }
	    id++;
	}
	st << "ET\n";
    }

    st << "Q\n";
    return;
} //

//---------------------------------------------------------------------------
// P1 Finite Element (also used in P1nc : Non-Conforming Finite Element)
//---------------------------------------------------------------------------

void trackP1isoline( std::vector<double> &px, std::vector<double> &py,
		     const double *const vx, const double *const vy,
		     const double value, const double *const vf )
{
    const double TOLERANCE = 1e-12;

    const int nVertices_in_each_Element = 3;

    for(int i = 0; i < nVertices_in_each_Element; i++){

	const double &x0 = vx[i];
	const double &y0 = vy[i];
	const double &f0 = vf[i];
		    
	const double &x1 = vx[(i+1)%3];
	const double &y1 = vy[(i+1)%3];
	const double &f1 = vf[(i+1)%3];

	if( (value < f0) && (value < f1) ) // no cossing point
	    continue;
		
	if( (f0 < value) && (f1 < value) ) // no crossing point
	    continue;

	if( (fabs(f0-f1) < TOLERANCE) && (fabs(f0-value) < TOLERANCE) ){ // coinside: edge is isoline
	    px.push_back( x0 );
	    py.push_back( y0 );
	    px.push_back( x1 );
	    py.push_back( y1 );
	    continue;
	}

	const double t = (value-f0)/(f1-f0); // (1-t)*f0 + t*f1 = value
	const double x = (1-t)*x0 + t*x1;
	const double y = (1-t)*y0 + t*y1;

	px.push_back( x );
	py.push_back( y );
    }

    if( px.size() == 3 ){   // if f(x,y)==value on either vertex

	if( (px[0] == px[1]) && (py[0] == py[1]) ){
	    px[1] = px[2];
	    py[1] = py[2];
	}

	//if( (px[1] == px[2]) && (py[1] == py[2]) ) // node 0 is differ from node 1
	//if( (px[2] == px[0]) && (py[2] == py[0]) ) // node 0 is differ from node 1
    }

    return;
}

void plot_P1_isoline( std::stringstream &Content, const Fem2D::Mesh &Th, const KN<double> &f_P1,
		      const KNM<double> &palette,
		      const int sizex, const int sizey, const double scale, const double ar,
		      const double x0, const double y0, const double y1,
		      const int marginl, const int marginb,
		      const double textfontsize, const bool monochrome,
		      const bool legend, const int prec, const bool logscale,
		      const double withmesh,
		      const int NISOLINES, const KN<double>*const viso,
		      const double linewidth )
{
    const double EPS = 1e-10;

    //------------------------------
    // values in plot
    //------------------------------
    const double fmax = (viso)? viso->max(): f_P1.max();
    const double fmin = (viso)? viso->min(): f_P1.min();
    std::vector<double> isoline_val;

    if( viso ){

	for(int m = 0; m < viso->size(); m++)
	    isoline_val.push_back( (*viso)[m] );

    } else if( logscale && (fmin > 0) ) {

	// fmin * step^N = fmax <=> step^N = fmax/fmin
	// <=> N = log_{step}(fmax/fmin) = (log(fmax/fmin))/log(step)
	// <=> log(step) = (1/N)(log(fmax/fmin))
	// <=> step = exp( (1/N)(log(fmax/fmin)) )
	const double df = exp( (static_cast<double>(1)/NISOLINES)*(log(fmax/fmin)) );

	isoline_val.push_back( fmin*sqrt(df) );
	for(int m = 1; m < NISOLINES; m++)
	    isoline_val.push_back( isoline_val[m-1] * df );

    } else {
	
	if( logscale )
	    std::cout << "plotPDF(): logscale for non-positive values.\n";

#if 1
	const double df = (fmax - fmin) / NISOLINES;
	for(int m = 0; m < NISOLINES; m++)
	    isoline_val.push_back( fmin + df/2 + m*df );
#else
	const double df = (fmax - fmin) / (NISOLINES+1);
	for(int m = 0; m < NISOLINES; m++)
	    isoline_val.push_back( fmin + (m+1)*df );
#endif
    }

    // color/gray-scale range to plot
#if 1
    // If user specifies an irrelevant color range (in viso array),
    const double cmax = *max_element( isoline_val.begin(), isoline_val.end() );
    const double cmin = *min_element( isoline_val.begin(), isoline_val.end() );
#else
    // map color/grayscale range to function values
    const double cmax = fmax;
    const double cmin = fmin;
#endif

    std::stringstream &st = Content;
    st.str("");

    //------------------------------
    // element (triangle)
    //------------------------------
    if( 0 < withmesh ){

	st << "q\n";
	st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm ";
	st << "1 w\n"; // setlinewidth

	const double grayscale1 = (withmesh < 1)? withmesh: 1;
	const double grayscale = 1-grayscale1;
	st << grayscale << ' ' << grayscale << ' ' << grayscale << " RG\n";

	for(int k = 0; k < Th.nt; k++){

	    const int &v0 = Th(k,0);
	    const int &v1 = Th(k,1);
	    const int &v2 = Th(k,2);

	    st << scale*ar*(Th(v0).x-x0) << ' ' << scale*(Th(v0).y-y0) << " m ";
	    st << scale*ar*(Th(v1).x-x0) << ' ' << scale*(Th(v1).y-y0) << " l ";
	    st << scale*ar*(Th(v2).x-x0) << ' ' << scale*(Th(v2).y-y0) << " l ";
	    st << "s" << std::endl;
	}

	st << "Q\n";

    } // withmesh

    //------------------------------
    // main routine
    //------------------------------
    st << "q\n";
    st << linewidth << " w\n"; // setlinewidth
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";
	
    const int &nTriangles = Th.nt;
    for(int k = 0; k < nTriangles; k++){

	const int &v0 = Th(k,0);
	const int &v1 = Th(k,1);
	const int &v2 = Th(k,2);

	const double vx[] = { Th(v0).x,     Th(v1).x,     Th(v2).x };
	const double vy[] = { Th(v0).y,     Th(v1).y,     Th(v2).y };
	const double vf[] = { f_P1[3*k+0],  f_P1[3*k+1],  f_P1[3*k+2] };

	for(size_t m = 0; m < isoline_val.size(); m++){

	    const double &value = isoline_val[m];

	    std::vector<double> px, py;
	    trackP1isoline( px, py, vx, vy, value, vf );

	    assert( px.size() == py.size() );

	    if( px.size() == 0 ) continue; // goto next isoline_val

	    setrgbcolor(st, value, palette, cmin, cmax, monochrome, logscale);

	    if( px.size() > 3 ){ 

		// f(x,y)==value on all vertices, i.e. f(x,y) \equiv value on the triangle
		st << "rg\n";
		st << scale*ar*(vx[0] - x0) << ' ' << scale*(vy[0] - y0) << " m "
		   << scale*ar*(vx[1] - x0) << ' ' << scale*(vy[1] - y0) << " l "
		   << scale*ar*(vx[2] - x0) << ' ' << scale*(vy[2] - y0) << " l f\n";
		
	    } else {

		// assert( (px.size() == 2) || (px.size() == 3) );
		st << "RG\n";
		st << scale*ar*(px[0] - x0) << ' ' << scale*(py[0] - y0) << " m "
		   << scale*ar*(px[1] - x0) << ' ' << scale*(py[1] - y0) << " l\n";
		st << "S\n";
	    }
	}
    }

    st << "Q\n";

    if( legend ){

	st << "q\n";
	st << "1 w\n"; // setlinewidth
	st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";

	for(size_t m = 0; m < isoline_val.size(); m++){

	    const double &f = isoline_val[m];

	    setrgbcolor(st, f, palette, cmin, cmax, monochrome, logscale);
	    st << "rg\n";

	    st << "BT /F1 " << textfontsize << " Tf "
	       << "1 0 0 1 " << sizex-PADDING << " " << (m+1)*(scale*(y1-y0)-textfontsize)/(isoline_val.size()+1) << " Tm "
	       << "(" << ((f >= 0)? "\\ ": "");

	    if( (fabs(f) > 1e-3) || (fabs(f) < 1e-12) ){
		st << std::setprecision(prec) << std::setfill('0') << f << ") Tj ET\n";
	    } else {
		st << std::resetiosflags(std::ios::fixed) << std::setiosflags(std::ios::scientific)
		   << std::setprecision(prec) << f << resetiosflags(std::ios::scientific)
		   << std::setiosflags(std::ios::fixed) << ") Tj ET\n";
	    }
	    // debug: st << std::setiosflags(std::ios::fixed)
	}

	st << "Q\n";

    } // legend

    //------------------------------
    // edges
    //------------------------------
    const int &nEdges = Th.neb;

    st << "q\n";
    st << "1 w\n"; // setlinewidth
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";
    st << "0 0 0 RG\n";

    for(int k = 0; k < nEdges; k++){

	const int &v0 = Th( Th.bedges[k][0] );
	const int &v1 = Th( Th.bedges[k][1] );

	st << scale*ar*(Th(v0).x-x0) << ' ' << scale*(Th(v0).y-y0) << " m ";
	st << scale*ar*(Th(v1).x-x0) << ' ' << scale*(Th(v1).y-y0) << " l s" << std::endl;
    }

    st << "Q\n";
    return;
}

//----------------------------------------------------------------------

void plot_P1_fill( std::stringstream &Content, const Fem2D::Mesh &Th, const KN<double> &f_P1_,
		   const KNM<double> &palette,
		   const int sizex, const int sizey, const double scale, const double ar,
		   const double x0, const double y0, const double y1,
		   const int marginl, const int marginb,
		   const double textfontsize, const bool monochrome,
		   const bool legend, const int prec, const bool logscale,
		   const double withmesh,
		   const long nbfill, const KN<double> *const frange )
{
    const double EPS = 1e-10;

    const int &nVertices  = Th.nv;
    const int &nTriangles = Th.nt;
    const int &nEdges     = Th.neb;
    const double &r       = scale;

    const double fmax = (frange)? (*frange)[1]: f_P1_.max();
    const double fmin = (frange)? (*frange)[0]: f_P1_.min();
    const double df   = (logscale)?
	pow( fmax/fmin, static_cast<double>(1)/nbfill ):
	(fmax - fmin)/nbfill;

    std::stringstream &st = Content;
    st.str("");

    st << "q\n";
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";

    //------------------------------
    // element(triangle)-wise process
    //------------------------------
    double *const f_P1 = new double [ Th.nv ];
    for(int k = 0; k < nTriangles; k++){

	const int &v0 = Th(k,0); const double &f0 = f_P1[v0];
	const int &v1 = Th(k,1); const double &f1 = f_P1[v1];
	const int &v2 = Th(k,2); const double &f2 = f_P1[v2];

	f_P1[v0] = f_P1_[3*k+0]; f_P1[v1] = f_P1_[3*k+1]; f_P1[v2] = f_P1_[3*k+2];

	// find minimum
	int v_min = v2;
	int end1 = v0; // destination point
	int end2 = v1; // destination point

	if( (f0 <= f1) && (f0 <= f2) ){
	    v_min = v0;
	    end1 = v1;
	    end2 = v2;
	} else if( (f1 <= f0) && (f1 <= f2) ){
	    v_min = v1;
	    end1 = v2;
	    end2 = v0;
	}

	if( frange && (fmax < f_P1[v_min]) )
	    continue;

	if( frange && (f_P1[end1] < fmin) && (f_P1[end2] < fmin) )
	    continue;

	int beg1 = v_min; // current point
	int beg2 = v_min; // current point
	double t1 = 0;    // secting parameters on edge beg1--end1
	double t2 = 0;    // secting parameters on edge beg2--end2

	int level = (logscale)?
	    static_cast<int>( log(f_P1[v_min]/fmin) / log(df) ):
	    static_cast<int>( (f_P1[v_min]-fmin) / df );
	    
	double f = (logscale)? (pow(df,level) * fmin) : (level*df + fmin);
	    
	if( frange && (f_P1[v_min] < fmin) ){

	    if( f_P1[end1] < fmin ){
		beg1 = end1;
		end1 = end2;
	    } if( f_P1[end2] < fmin ){
		beg2 = end2;
		end2 = end1;
	    }

	    t1 = (fmin - f_P1[beg1]) / (f_P1[end1] - f_P1[beg1]);
	    t2 = (fmin - f_P1[beg2]) / (f_P1[end2] - f_P1[beg2]);
	    level = 0;
	    f = fmin;
	}

	do {
	    const double EPS = 1e-10;
	    if( (end1 == end2) && (t1 >= 1-EPS) && (t2 >= 1-EPS) )
		break;

	    // find points which divide edges by t:(1-t1), t2:(1-t2) respectively
	    const double p1x = (1-t1)*Th(beg1).x + t1*Th(end1).x;
	    const double p1y = (1-t1)*Th(beg1).y + t1*Th(end1).y;
	    const double p2x = (1-t2)*Th(beg2).x + t2*Th(end2).x;
	    const double p2y = (1-t2)*Th(beg2).y + t2*Th(end2).y;

	    const double f_next = (logscale)? f * df: f + df;

	    if( frange && (fmax-EPS <= f) )
		break;

	    if( level <= 0 ){
		setrgbcolor(st, fmin, palette, fmin, fmax, monochrome, logscale);
	    } else if( level >= nbfill-1 ){
		
		// exceptional handling : if local minimum (in triangle) is
		// already global maximum, then level == nbfill.
		// Thus the condition is not level == nbfill-1,
		// but level >= nbfill-1.
		setrgbcolor(st, fmax, palette, fmin, fmax, monochrome, logscale);

	    } else {
		const double c = (logscale)?
		    (pow(df,level+0.5)*fmin): ((level+0.5)*df + fmin);
		setrgbcolor(st, c, palette, fmin, fmax, monochrome, logscale);
	    }
	    st << "rg\n";
		
	    if( (f_next >= f_P1[end1]) && (f_next >= f_P1[end2]) ){

		st << r*ar*(p1x-x0)        << ' ' << r*(p1y-y0)        << " m "
		   << r*ar*(Th(end1).x-x0) << ' ' << r*(Th(end1).y-y0) << " l "
		   << r*ar*(Th(end2).x-x0) << ' ' << r*(Th(end2).y-y0) << " l "
		   << r*ar*(p2x-x0)        << ' ' << r*(p2y-y0)        << " l f\n";

		break;
	    }

	    if( f_next >= f_P1[end1] ){
		  
		st << r*ar*(p1x-x0)        << ' ' << r*(p1y-y0)        << " m "
		   << r*ar*(Th(end1).x-x0) << ' ' << r*(Th(end1).y-y0) << " l "
		   << r*ar*(p2x-x0)        << ' ' << r*(p2y-y0)        << " l f\n";

		beg1 = end1;
		end1 = end2;
		t1 = 0;
		continue;
	    } 

	    if( f_next >= f_P1[end2] ){

		st << r*ar*(p1x-x0) << ' ' << r*(p1y-y0) << " m "
		   << r*ar*(Th(end2).x-x0) << ' ' << r*(Th(end2).y-y0) << " l "
		   << r*ar*(p2x-x0) << ' ' << r*(p2y-y0) << " l f\n";
		  
		beg2 = end2;
		end2 = end1;
		t2 = 0;
		continue;
	    }

	    const double t1_next = (f_next - f_P1[beg1]) / (f_P1[end1] - f_P1[beg1]);
	    const double t2_next = (f_next - f_P1[beg2]) / (f_P1[end2] - f_P1[beg2]);

	    const double q1x = (1-t1_next)*Th(beg1).x + t1_next*Th(end1).x;
	    const double q1y = (1-t1_next)*Th(beg1).y + t1_next*Th(end1).y;
	    const double q2x = (1-t2_next)*Th(beg2).x + t2_next*Th(end2).x;
	    const double q2y = (1-t2_next)*Th(beg2).y + t2_next*Th(end2).y;

	    st << r*ar*(p1x-x0) << ' ' << r*(p1y-y0) << " m "
	       << r*ar*(q1x-x0) << ' ' << r*(q1y-y0) << " l "
	       << r*ar*(q2x-x0) << ' ' << r*(q2y-y0) << " l "
	       << r*ar*(p2x-x0) << ' ' << r*(p2y-y0) << " l f\n";

	    // update
	    f = f_next;
	    level++;
	    t1 = t1_next;
	    t2 = t2_next;
	} while(level <= nbfill);

    } // element(triangle)-wise process

    delete [] f_P1;

    //------------------------------
    // legend
    //------------------------------
    if( legend ){

	const double dy = r*(y1-y0)/nbfill;

	for(int m = 0; m < nbfill; m++){

	    if( m == 0 ){
		setrgbcolor(st, fmin, palette, fmin, fmax, monochrome, logscale);
	    } else if( m == nbfill-1 ){
		setrgbcolor(st, fmax, palette, fmin, fmax, monochrome, logscale);
	    } else {
		const double f = (logscale)? fmin * pow(df,m+0.5): fmin + (m+0.5)*df;
		setrgbcolor(st, f, palette, fmin, fmax, monochrome, logscale);
	    }
	    st << "rg\n";

	    st << sizex-PADDING   << " " << m*dy << " m "
               << sizex-PADDING/2 << " " << m*dy << " l "
               << sizex-PADDING/2 << " " << (m+1)*dy << " l "
	       << sizex-PADDING   << " " << (m+1)*dy << " l f\n";
	}

	const double EPS = 1e-10;
        const double dl = (logscale)?
	    pow( fmax/fmin, static_cast<double>(1)/(NUM_LABELS-1) ):
	    (fmax-fmin)/(NUM_LABELS-1);

        for(int m = 0; m < NUM_LABELS; m++){

	    const double f = (logscale)? fmin * pow(dl,m): fmin + m*dl;

	    if( logscale ){

		if( f <= fmin*df ){
		    setrgbcolor(st, fmin, palette, fmin, fmax, monochrome, logscale);
		} else if( f >= fmax/df ){
		    setrgbcolor(st, fmax, palette, fmin, fmax, monochrome, logscale);
		} else {
		    const double dc = pow( fmax/fmin, static_cast<double>(1)/nbfill );
		    const int mc = static_cast<int>( log(f/fmin) / log(dc) );
		    const double c = fmin * pow(dc, mc+0.5);
		    setrgbcolor(st, c, palette, fmin, fmax, monochrome, logscale);
		}

	    } else {

		if( f <= fmin+df ){
		    setrgbcolor(st, fmin, palette, fmin, fmax, monochrome, logscale);
		} else if( f >= fmax-df ){
		    setrgbcolor(st, fmax, palette, fmin, fmax, monochrome, logscale);
		} else {
		    const double dc = (fmax-fmin) / nbfill;
		    const int mc = static_cast<int>( (f-fmin)/dc );
		    const double c = fmin + (mc+0.5)*dc;
		    setrgbcolor(st, c, palette, fmin, fmax, monochrome, logscale);
		}
	    }
	    st << " rg\n";

            st << "BT /F1 " << textfontsize << " Tf "
	       << "1 0 0 1 " << sizex << " " << m*(r*(y1-y0)-textfontsize)/(NUM_LABELS-1) << " Tm "
	       << "(" << ((f >= 0)? "\\ ": "");

	    if( (fabs(f) > 1e-3) || (fabs(f) < 1e-12) ){
		st << std::setprecision(prec) << std::setfill('0') << f << ") Tj ET\n";
	    } else {
		st << std::resetiosflags(std::ios::fixed) << std::setiosflags(std::ios::scientific)
		   << std::setprecision(prec) << f << resetiosflags(std::ios::scientific)
		   << std::setiosflags(std::ios::fixed) << ") Tj ET\n";
	    }
	    // debug: st << std::setiosflags(std::ios::fixed)
        }
    } // legend

    //------------------------------
    // element (triangle)
    //------------------------------
    if( 0 < withmesh ){

	const double grayscale1 = (withmesh < 1)? withmesh: 1;
	const double grayscale = 1-grayscale1;
	st << grayscale << ' ' << grayscale << ' ' << grayscale << " RG\n";

	for(int n = 0; n < Th.nt; n++){

	    const int &v0 = Th(n,0);
	    const int &v1 = Th(n,1);
	    const int &v2 = Th(n,2);

	    st << r*ar*(Th(v0).x-x0) << ' ' << r*(Th(v0).y-y0) << " m ";
	    st << r*ar*(Th(v1).x-x0) << ' ' << r*(Th(v1).y-y0) << " l ";
	    st << r*ar*(Th(v2).x-x0) << ' ' << r*(Th(v2).y-y0) << " l ";
	    st << "s" << std::endl;
	}
    } // withmesh

    //------------------------------
    // edges
    //------------------------------
    st << "0 0 0 RG\n";
    for(int k = 0; k < nEdges; k++){

	const int &v0 = Th( Th.bedges[k][0] );
	const int &v1 = Th( Th.bedges[k][1] );

	st << r*ar*(Th(v0).x-x0) << ' ' << r*(Th(v0).y-y0) << " m "
	   << r*ar*(Th(v1).x-x0) << ' ' << r*(Th(v1).y-y0) << " l S\n";
    }

    st << "Q\n";
    return;
} //

//----------------------------------------------------------------------
// P0 Finite Element
//----------------------------------------------------------------------

void plot_P0_fill( std::stringstream &Content, const Fem2D::Mesh &Th, const KN<double> &f_P0,
		   const KNM<double> &palette,
		   const int sizex, const int sizey, const double scale, const double ar,
		   const double x0, const double y0, const double y1,
		   const int marginl, const int marginb,
		   const double textfontsize, const bool monochrome,
		   const bool legend, const int prec, const bool logscale,
		   const double withmesh,
		   const long nbfill, const KN<double> *const frange )
{
    const double EPS = 1e-10;

    const int &nVertices  = Th.nv;
    const int &nTriangles = Th.nt;
    const int &nEdges     = Th.neb;
    const double &r       = scale;

    const double fmax = (frange)? (*frange)[1]: f_P0.max();
    const double fmin = (frange)? (*frange)[0]: f_P0.min();
    const double df   = (logscale)?
	exp( (static_cast<double>(1)/nbfill)*(log(fmax/fmin)) ):
	(fmax - fmin)/nbfill;

    std::stringstream &st = Content;
    st.str("");

    st << "q\n";
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";

    //------------------------------
    // element(triangle)-wise process
    //------------------------------
    for(int k = 0; k < nTriangles; k++){

	const int &v0 = Th(k,0);
	const int &v1 = Th(k,1);
	const int &v2 = Th(k,2);

	const double &f = f_P0[k];

	if( frange && (f < fmin) )
	    continue;
	if( frange && (f > fmax) )
	    continue;

	int level = (logscale)?
	    static_cast<int>( log(f/fmin) / log(df) ):
	    static_cast<int>( (f-fmin) / df );
	    
	if( level == 0 ){
	    setrgbcolor(st, fmin, palette, fmin, fmax, monochrome, logscale);
	} else if( level >= nbfill-1 ){
	    setrgbcolor(st, fmax, palette, fmin, fmax, monochrome, logscale);
	} else {
	    const double c = (logscale)? (pow(df,level+0.5)*fmin): ((level+0.5)*df + fmin);
	    setrgbcolor(st, c, palette, fmin, fmax, monochrome, logscale);
	}
	st << "rg\n";
		
	st << r*ar*(Th(v0).x-x0) << ' ' << r*(Th(v0).y-y0) << " m "
	   << r*ar*(Th(v1).x-x0) << ' ' << r*(Th(v1).y-y0) << " l "
	   << r*ar*(Th(v2).x-x0) << ' ' << r*(Th(v2).y-y0) << " l f\n";

    } // element(triangle)-wise process

    //------------------------------
    // legend
    //------------------------------
    if( legend ){

	const int &LEVELS  = nbfill;
	const double dy = r*(y1-y0)/LEVELS;

	for(int m = 0; m < LEVELS; m++){

	    if( m == 0 ){
		setrgbcolor(st, fmin, palette, fmin, fmax, monochrome, logscale);
	    } else if( m == LEVELS-1 ){
		setrgbcolor(st, fmax, palette, fmin, fmax, monochrome, logscale);
	    } else {
		const double f = (logscale)? fmin * pow(df,m+0.5): fmin + (m+0.5)*df;
		setrgbcolor(st, f, palette, fmin, fmax, monochrome, logscale);
	    }
	    st << "rg\n";

	    st << sizex-PADDING   << " " << m*dy << " m "
	       << sizex-PADDING/2 << " " << m*dy << " l "
	       << sizex-PADDING/2 << " " << (m+1)*dy << " l "
	       << sizex-PADDING   << " " << (m+1)*dy << " l f\n";
	}

	const double EPS = 1e-10;
        const double dl = (logscale)?
	    pow( fmax/fmin, static_cast<double>(1)/(NUM_LABELS-1) ):
	    (fmax-fmin)/(NUM_LABELS-1);

	for(int m = 0; m < NUM_LABELS; m++){

	    const double f = (logscale)? fmin*pow(dl,m): fmin + m*dl;

	    if( logscale ){

		if( f <= fmin*df ){
		    setrgbcolor(st, fmin, palette, fmin, fmax, monochrome, logscale);
		} else if( f >= fmax/df ){
		    setrgbcolor(st, fmax, palette, fmin, fmax, monochrome, logscale);
		} else {
		    const double dc = pow( fmax/fmin, static_cast<double>(1)/LEVELS );
		    const int mc = static_cast<int>( log(f/fmin) / log(dc) );
		    const double c = fmin * pow(dc, mc+0.5);
		    setrgbcolor(st, c, palette, fmin, fmax, monochrome, logscale);
		}

	    } else {

		if( f <= fmin+df ){
		    setrgbcolor(st, fmin, palette, fmin, fmax, monochrome, logscale);
		} else if( f >= fmax-df ){
		    setrgbcolor(st, fmax, palette, fmin, fmax, monochrome, logscale);
		} else {
		    const double dc = (fmax-fmin) / LEVELS;
		    const int mc = static_cast<int>( (f-fmin)/dc );
		    const double c = fmin + (mc+0.5)*dc;
		    setrgbcolor(st, c, palette, fmin, fmax, monochrome, logscale);
		}
	    }
	    st << " rg\n";
	    
	    st << "BT /F1 " << textfontsize << " Tf "
	       << "1 0 0 1 " << sizex << " " << m*(r*(y1-y0)-textfontsize)/(NUM_LABELS-1) << " Tm "
	       << "(" << ((f >= 0)? "\\ ": "");

	    if( (fabs(f) > 1e-3) || (fabs(f) < 1e-12) ){
		st << std::setprecision(prec) << std::setfill('0') << f << ") Tj ET\n";
	    } else {
		st << std::resetiosflags(std::ios::fixed) << std::setiosflags(std::ios::scientific)
		   << std::setprecision(prec) << f << resetiosflags(std::ios::scientific)
		   << std::setiosflags(std::ios::fixed) << ") Tj ET\n";
	    }
	    // debug: st << std::setiosflags(std::ios::fixed)
	}
    } // legend

    //------------------------------
    // element (triangle)
    //------------------------------
    if( 0 < withmesh ){

	const double grayscale1 = (withmesh < 1)? withmesh: 1;
	const double grayscale = 1-grayscale1;
	st << grayscale << ' ' << grayscale << ' ' << grayscale << " RG\n";

	for(int n = 0; n < Th.nt; n++){
	    
	    const int &v0 = Th(n,0);
	    const int &v1 = Th(n,1);
	    const int &v2 = Th(n,2);

	    st << r*ar*(Th(v0).x-x0) << ' ' << r*(Th(v0).y-y0) << " m ";
	    st << r*ar*(Th(v1).x-x0) << ' ' << r*(Th(v1).y-y0) << " l ";
	    st << r*ar*(Th(v2).x-x0) << ' ' << r*(Th(v2).y-y0) << " l ";
	    st << "s" << std::endl;
	}
    } // withmesh

    //------------------------------
    // edges
    //------------------------------
    st << "0 0 0 RG\n";
    for(int k = 0; k < nEdges; k++){

	const int &v0 = Th( Th.bedges[k][0] );
	const int &v1 = Th( Th.bedges[k][1] );

	st << r*ar*(Th(v0).x-x0) << ' ' << r*(Th(v0).y-y0) << " m "
	   << r*ar*(Th(v1).x-x0) << ' ' << r*(Th(v1).y-y0) << " l S\n";
    }

    st << "Q\n";
    return;
} //

//----------------------------------------------------------------------
// P2 Finite Element
//----------------------------------------------------------------------

void findQuadraticPolynomial(double *const phi, const double *const vx, const double *const vy, const double *const f_P2)
{
    const int NQ = 6; // number of unknowns in quadratic polynomial

    double *a[NQ];
    a[0] = new double [NQ*(NQ+1)];  // NQ*(NQ+1) matrix for Gauss Elimination
    for(int i = 1; i < NQ; i++)
	a[i] = a[i-1] + (NQ+1);

    const double ex[] = { (vx[1]+vx[2])/2, (vx[2]+vx[0])/2, (vx[0]+vx[1])/2 };
    const double ey[] = { (vy[1]+vy[2])/2, (vy[2]+vy[0])/2, (vy[0]+vy[1])/2 };

    const double x[] = { vx[0], vx[1], vx[2], ex[0], ex[1], ex[2] };
    const double y[] = { vy[0], vy[1], vy[2], ey[0], ey[1], ey[2] };


    for(int i = 0; i < NQ; i++){
	a[i][0] = x[i] * x[i];
	a[i][1] = x[i] * y[i];
	a[i][2] = y[i] * y[i];
	a[i][3] = x[i];
	a[i][4] = y[i];
	a[i][5] = 1;
	a[i][6] = f_P2[i];
    }

    GaussElimination(phi, a, NQ);

    delete [] a[0];
    
    return;
}

void findCanonicalForm( double *const PHI, const double *const phi )
{
    const double &a = phi[0]; const double &b = phi[1]; const double &c = phi[2];
    const double &d = phi[3]; const double &e = phi[4]; const double &f = phi[5];
    
    // phi(x,y) = a*x*x + b*x*y + c*y*y + d*x + e*y + f
    //          = [x,y][ a, b/2; b/2, c ][x;y] + [d,e][x;y] + f
    //          = [x,y]P P^T [ a, b/2; b/2, c ]P P^T[x;y] + [d,e]P P^T[x;y] + f
    //          = [X,Y] [lambda1,0;0,lambda2] [X;Y] + [D,E][X;Y] + f
    //          = lambda1*X*X + lambda2*Y*Y + D*X + E*Y + f
    //          = lambda1 * (X + D/(2*lambda1))^2 + lambda2 * (Y + E/(2*lambda2))^2
    //              + ( -D*D/(4*lambda1) - E*E/(4*lambda2) + f )
    //          = PHI(X,Y)
    // where
    //   v1 = [v1x;v1y] (resp. v2 = [v2x;v2y]) is an eigenvector of lambda1 (resp. lambda2)
    //   P = [v1x, v2x; v1y, v2y],  satisfying P^T P = I <=> || v1 || = || v2 || = 1, v1 \perp v2 = 0
    //   [X;Y] = P^T[x;y] <=> [x;y] = P[X;Y]
    //   [D;E] = P^T[d;e] <=> [d;e] = P[D;E]

    const double det   = (a-c)*(a-c)+b*b;
    const double sqdet = sqrt( det );

    // eigenvalues of [ a, b /2; b/2, c]
    double &lambda1 = PHI[0]; double &lambda2 = PHI[1];
    lambda1 = ( (a+c) + sqdet ) / 2;
    lambda2 = ( (a+c) - sqdet ) / 2;

    // corresponding eivenvectors
    double &v1x = PHI[2]; double &v1y= PHI[3];
    double &v2x = PHI[4]; double &v2y= PHI[5];

    if( a < c ){

	const double n = sqrt( 2*det - 2*(a-c)*sqdet );
	v1x = b / n;
	v1y = (-(a-c)+sqdet) / n;
	v2x = (a-c-sqdet) / n;
	v2y = b / n;

    } else if( a > c ){

	const double n = sqrt( 2*det + 2*(a-c)*sqdet );
	v1x = (a-c+sqdet) / n;
	v1y = b / n;
	v2x = b / n;
	v2y = (-(a-c)-sqdet) / n;

    } else { // a == c

	lambda1 = (2*a+b)/2;
	lambda2 = (2*a-b)/2;
	v1x = v1y = v2x = v2y = 1/sqrt(static_cast<double>(2));  v2y *= -1;
    }

    double &D = PHI[6]; double &E = PHI[7]; double &F = PHI[8];

    // [D;E] = P^T[d;e]
    const double p[2][2]  = { { v1x, v2x }, { v1y, v2y } };
    const double pT[2][2] = { { p[0][0], p[1][0] }, { p[0][1], p[1][1] } };
    D = pT[0][0] * d + pT[0][1] * e;
    E = pT[1][0] * d + pT[1][1] * e;

    const double EPS = 1e-10;

    F = f;
    if( fabs(lambda1) > EPS )
	F -= D*D/(4*lambda1);

    if( fabs(lambda2) > EPS )
	F -= E*E/(4*lambda2);

    return;
}

void transformTriangle( double *const Vx, double *const Vy,
			const double *const vx, const double *const vy, const double *const PHI )
{
    const double &lambda1 = PHI[0]; const double &lambda2 = PHI[1];
    const double &D = PHI[6]; const double &E = PHI[7];

    const double &ev1x = PHI[2]; const double &ev1y = PHI[3];
    const double &ev2x = PHI[4]; const double &ev2y = PHI[5];
    const double  P[2][2] = { { ev1x, ev2x }, { ev1y, ev2y } };
    const double PT[2][2] = { { P[0][0], P[1][0] }, { P[0][1], P[1][1] } };

    const int nVertices = 3;
    
    const double EPS = 1e-10;

    for(int i = 0; i < nVertices; i++){
	Vx[i] = PT[0][0] * vx[i] + PT[0][1] * vy[i];
	if( fabs(lambda1) > EPS )
	    Vx[i] += D/(2*lambda1);
	Vy[i] = PT[1][0] * vx[i] + PT[1][1] * vy[i];
	if( fabs(lambda2) > EPS )
	    Vy[i] += E/(2*lambda2);
    }

    return;
}

void findZeros( std::vector<double> &zx, std::vector<double> &zy,
		const double v1x, const double v1y, const double v2x, const double v2y,
		const double *const phi, const double value )
{
    // crosss points of ax^2 + bxy + cy^2 + dx + ey + f = value
    // and segment : [x;y] = (1-t)*[v1x;v1y] + t*[v2x;v2y] with 0 <= t <= 1
    const double &a = phi[0]; const double &b = phi[1]; const double &c = phi[2];
    const double &d = phi[3]; const double &e = phi[4]; const double &f = phi[5];

    const double EPS = 1e-10;

    // If phi(x,y) == value at both v1 and v2, then return
    const double phi1 = a*v1x*v1x + b*v1x*v1y + c*v1y*v1y + d*v1x + e*v1y + f - value;
    const double phi2 = a*v2x*v2x + b*v2x*v2y + c*v2y*v2y + d*v2x + e*v2y + f - value;

    if( fabs(phi1) + fabs(phi2) < EPS ){

	// examine phi at mid point of v1 and v2
	const double v3x = (v1x + v2x) / 2;
	const double v3y = (v1y + v2y) / 2;
	const double phi3 = a*v3x*v3x + b*v3x*v3y + c*v3y*v3y + d*v3x + e*v3y + f - value;
	
	// If phi(mid point) == 0, then phi(x,y) may vanish on the segment (line) v1--v2.
	if( fabs(phi3) < EPS )
	    return;

	zx.push_back( v1x ); zy.push_back( v1y );
	zx.push_back( v2x ); zy.push_back( v2y );
	return;
    }

    // At^2 + Bt + C = 0
    const double A = a*(v1x-v2x)*(v1x-v2x) + b*(v1x-v2x)*(v1y-v2y) + c*(v1y-v2y)*(v1y-v2y);
    const double B = 2*a*v1x*(-v1x+v2x) + b*v1y*(-v1x+v2x) + b*v1x*(-v1y+v2y) + 2*c*v1y*(-v1y+v2y)
	+ d*(-v1x+v2x) + e*(-v1y+v2y);
    const double C = a*v1x*v1x + b*v1x*v1y + c*v1y*v1y + d*v1x + e*v1y + f - value;

    // If A == 0, then solve Bt + C = 0
    if( fabs(A) < EPS ){

	if(fabs(B) < EPS){ // A == B == 0
	    //if(fabs(C) > EPS ){
	    //    cerr << "findZeros(): no solutions" << endl;
	    //}
	    return;
	}

	const double t = -C/B;

	if( (-EPS < t) && (t < 1+EPS) ){ // 0 <= t <= 1
	    const double x = (1-t)*v1x + t*v2x;
	    const double y = (1-t)*v1y + t*v2y;
	    zx.push_back(x);
	    zy.push_back(y);
	}
	return;
    }

    // solve At^2 + Bt + C = 0
    const double D_ = B*B - 4*A*C;
    const double D = (fabs(D_) < EPS)? 0: D_;

    if( D < 0 ) return; // no cross points
    
    const double t01[2] = { (-B + sqrt(D))/(2*A), (-B - sqrt(D))/(2*A) };
    for(int i = 0; i < 2; i++){

        const double &t = t01[i]; 

        if( (-EPS < t) && (t < 1+EPS) ){ // 0 <= t <= 1
	
	    const double x = (1-t)*v1x + t*v2x;
	    const double y = (1-t)*v1y + t*v2y;
            zx.push_back( x );
	    zy.push_back( y );
        }
    }
    return;
}

void invTransformCubicBzeirs( std::vector< std::vector<double> > &Cxs, std::vector< std::vector<double> > &Cys,
			      const double *const PHI )
{
    const double EPS = 1e-10;

    const double &lambda1 = PHI[0]; const double &lambda2 = PHI[1];
    const double &D = PHI[6]; const double &E = PHI[7];

    const double &ev1x = PHI[2]; const double &ev1y = PHI[3];
    const double &ev2x = PHI[4]; const double &ev2y = PHI[5];
    const double  P[2][2] = { { ev1x, ev2x }, { ev1y, ev2y } };

    assert( Cxs.size() == Cys.size() );

    if( Cxs.size() == 0 ) return; // nothing to transform

    for(size_t k = 0; k < Cxs.size(); k++){

	std::vector<double> &Cx = Cxs[k];
	std::vector<double> &Cy = Cys[k];

	assert( Cx.size() == Cy.size() );

	for(size_t j = 0; j < Cx.size(); j++){ // j=0 : starting point

	    double Px = Cx[j];
	    double Py = Cy[j];

	    if( fabs(lambda1) > EPS )
		Px -= D/(2*lambda1);
	    if( fabs(lambda2) > EPS )
		Py -= E/(2*lambda2);

	    Cx[j] = P[0][0] * Px + P[0][1] * Py;
	    Cy[j] = P[1][0] * Px + P[1][1] * Py;
	}
    } // k
    return;
}

void drawCubicBeziers( std::stringstream &Content,
		       const std::vector< std::vector<double> > &Cxs,
		       const std::vector< std::vector<double> > &Cys,
		       const double scale, const double ar, const double x0, const double y0 )
{
    const double EPS = 1e-10;

    std::stringstream &st = Content;

    assert( Cxs.size() == Cys.size() );

    if( Cxs.size() == 0 ) return; // nothing to draw

    for(size_t k = 0; k < Cxs.size(); k++){

	const std::vector<double> &px = Cxs[k];
	const std::vector<double> &py = Cys[k];

	assert( px.size() == py.size() );
	assert( ( px.size()-1 ) % 3 == 0 ); // ignore the starting point

	st << scale*ar*(px[0]-x0) << ' ' << scale*(py[0]-y0) << " m\n";

	// sequential control points
	for(size_t j = 1; j < px.size(); j += 3){ // j=0 : starting point

	    for(int i = 0; i < 3; i++)
		st << scale*ar*(px[j+i]-x0) << ' ' << scale*(py[j+i]-y0) << ' ';
	    st << "c\n";
	}
	st << "S\n";
    } // k
    return;
}

bool isInsideTriangle( const double px, const double py, const double *const vx, const double *const vy )
{
    const double EPS = 1e-10;	
    // Find a, b \in \Real such that [px;py]-v0 =  a (v1-v0) + b (v2-v0).
    // if (0 < a < 1) && (0 < b < 1) && (0 < a+b < 1), then p is inside triangle v0--v1--v2

    // solve (v1x-v0x) * a + (v2x-v0x) * b = px-v0x
    //       (v1y-v0y) * a + (v2y-v0y) * b = py-v0y

    const double m00 = vx[1] - vx[0];
    const double m01 = vx[2] - vx[0];
    const double rhs0 = px - vx[0];

    const double m10 = vy[1] - vy[0];
    const double m11 = vy[2] - vy[0];
    const double rhs1 = py - vy[0];

    const double det = m00*m11 - m01*m10;
    assert( fabs(det) > EPS );

    const double a = ( m11 * rhs0 - m01 * rhs1) / det;
    const double b = (-m10 * rhs0 + m00 * rhs1) / det;

    return (0 < a) && (a < 1) && (0 < b) && (b < 1) && (0 < a+b) && (a+b < 1);
}

void trackParabolaCore( std::vector< std::vector<double> > &Cx, std::vector< std::vector<double> > &Cy,
			const double a, const double b, std::vector<double> &x,
			const double *const Vx, const double *const Vy )
{
    // y = a * x^2 + b
    std::sort( x.begin(), x.end() );

    // control points: Px[4], Py[4]
    // Bezier curve : P(t) = (1-t)^3 P0 + 3(1-t)^2 t P1 + 3(1-t)t^2 P2 + t^3 P3
    //                     = (-P0+3P1-3P2+P3) t^3 + (3P0-6P1+3P2) t^2 + (-3P0+3P1) t + P0
    //
    // parabola between x0 <= x <= x1 : 
    // P(t) = [ (x1-x0)t + x0 ; a*( (x1-x0)t + x0 )^2 + b*( (x1-x0)t + x0 ) + c ]
    //      = [ (x1-x0)t + x0 ; a(x1-x0)^2 t^2 + (2*a*x0+b)*(x1-x0) t + (a x_0^2 + bx_0 + c) ]
    //
    // Therefore
    //   P0 = [ x0 ; ax_0^2 + bx_0 + c ]
    //   -3P0+3P1 = [ x1-x0 ; (2*a*x0+b)*(x1-x0) ]
    //   3P0-6P1+3P2 = [ 0 ; a(x1-x0)^2 ]
    //   -P0+3P1-3P2+P3 = [ 0 ; 0 ]
    // <=>
    //   P0 = [ x0 ; ax_0^2 + bx_0 + c ]
    //   P1 = P0 + [ x1-x0 ; (2*a*x0+b)*(x1-x0) ]/3
    //   P2 = -P0 + 2*P1 + [ 0 ; a(x1-x0)^2 ]/3
    //   P3 = P0-3P1+3P2
    // From this, P0x = x0, P1x = x0+(x1-x0)/3, P2x = P0x + 2/3*(x1-x0), P3x = x1
	
    for(int i = 0; i+1 < x.size(); i++){

	const double &X0 = x[i];
	const double &X1 = x[i+1];
	const double h   = X1 - X0;

	const double EPS = 1e-10;
	if( h < EPS)
	    continue;

	// examine the arc X[i]--X[i+1] is inside triangle or not
	const double mX0 = X0 + h/100;
	const double mY0 = a*mX0*mX0 + b;

	const double mX1 = X1 - h/100;
	const double mY1 = a*mX1*mX1 + b;
	if( !(isInsideTriangle( mX0, mY0, Vx, Vy ) && isInsideTriangle( mX1, mY1, Vx, Vy )) )
	    continue;

	const double P0y = a*X0*X0 + b;
	const double P1y =  P0y + (2*a*X0)*h/3;
	const double P2y = -P0y + 2*P1y + a*h*h/3;
	const double P3y =  P0y - 3*P1y + 3*P2y;

	Cx.push_back( std::vector<double> { X0,  X0+h/3, X1-h/3, X1  } );
	Cy.push_back( std::vector<double> { P0y, P1y,    P2y,    P3y } );
    }
    return;
}

void trackParabola( std::vector< std::vector<double> > &Cx, std::vector< std::vector<double> > &Cy,
		    const double *const PHI, const std::vector<double> &zx, const std::vector<double> &zy,
		    const double *const Vx, const double *const Vy )

{
    const double EPS = 1e-10;

    // PHI(X,Y) = lambda1*X*X + lambda2*Y*Y + D*X + E*Y + f
    //          = lambda1 * (X + D/(2*lambda1))^2 + lambda2 * (Y + E/(2*lambda2))^2
    //            + ( -D*D/(4*lambda1) - E*E/(4*lambda2) + f)
    //          = lambda1*(X + D/(2*lambda1))^2 + lambda2*(Y + E/(2*lambda2))^2 + F
    // X' = X + D/(2*lambda1), Y' = Y + E/(2*lambda2)
    // They are new variables, (not derivatives of X and Y)
    const double &lambda1 = PHI[0]; const double &lambda2 = PHI[1];
    const double &D = PHI[6]; const double &E = PHI[7]; const double &F = PHI[8];

    const double &ev1x = PHI[2]; const double &ev1y = PHI[3];
    const double &ev2x = PHI[4]; const double &ev2y = PHI[5];
    const double  P[2][2] = { { ev1x, ev2x }, { ev1y, ev2y } };
    const double PT[2][2] = { { P[0][0], P[1][0] }, { P[0][1], P[1][1] } };

    assert( zx.size() == zy.size() );
#if 1
    std::vector<double> Zx, Zy;
    for(size_t i = 0; i < zx.size(); i++){
	Zx.push_back( PT[0][0]*zx[i] + PT[0][1]*zy[i] );
	Zy.push_back( PT[1][0]*zx[i] + PT[1][1]*zy[i] );
    }
#endif
    if( fabs(lambda1) > EPS ){

	// lambda1*X^2 + D*X + E*Y + F = 0
        // <=> lambda1*( X + D/(2*lambda1) )^2 + E*Y + F - D*D/(4*lambda1) = 0
        // <=> Y = (-lambda1/E) * (X')^2 - F'/E
	if( fabs(E) < EPS ) return;

	for(std::vector<double>::iterator itr = Zx.begin(); itr != Zx.end(); itr++)
	    *itr += D/(2*lambda1);

	const double a = -lambda1/E; // Y' = a*(X')^2 + b
	const double b = -F/E;

	trackParabolaCore( Cx, Cy, a, b, Zx, Vx, Vy );

    } else {

	assert( fabs(lambda2) > EPS );

	// D*X + lambda2*Y^2 + E*Y + F = 0 <=> X = (-lambda2/D)*Y^2 - (E/D)*Y - F/D
	// <=> D*X + lambda2*(Y - E/(2*lambda2))^2 + F - E*E/(4*lambda2) = 0
        // <=> X = (-lambda2/D)*(Y')^2 - F'/D
	if( fabs(D) < EPS ) return;

	for(std::vector<double>::iterator itr = Zy.begin(); itr != Zy.end(); itr++)
	    *itr += E/(2*lambda2);

	const double a = -lambda2/D; // X' = a*(Y')^2 + b
	const double b = -F/D;

	trackParabolaCore( Cy, Cx, a, b, Zy, Vy, Vx );
    }

    return;
}

void trackEllipse( std::vector< std::vector<double> > &Cxs, std::vector< std::vector<double> > &Cys,
		   const double *const PHI, const double *const Vx, const double *const Vy )
{
    const double &lambda1 = PHI[0]; const double &lambda2 = PHI[1];
    const double &D = PHI[6]; const double &E = PHI[7]; const double &F = PHI[8];

    assert( lambda1*lambda2 > 0 );

    // lambda1*(X')^2 + lambda2*(Y')^2 + F = 0
    // Y' = 0 => X'= sqrt( -F/lambda1 ), it means that -F/lambda1 > 0
    if( -F/lambda1 <= 0 )
	return;

    // lambda1 * X*X + lambda2 * Y*Y + F = 0
    // <=> (-lambda1/F) * X*X + (-lambda2/F) * Y*Y = 1
    // <=> (X/a)^2 + (Y/b)^2 = 1, 1/a = sqrt(-lambda1/F), 1/b = sqrt(-lambda2/F)
    const double a = sqrt( -F / lambda1 );
    const double b = sqrt( -F / lambda2 );

    // Ellipse is tangent to an edge of the triangle element, and localtes outside of the triangle
    // Examine both opposite sides (a,0) and (-a,0) are belong to inside of the triangle element
    if( !isInsideTriangle( -a, 0, Vx, Vy ) && !isInsideTriangle( -a, 0, Vx, Vy ) )
	return;

    const double PI = atan(static_cast<double>(1)) * 4;
    const double c  = 35*(32/(PI*PI*PI) - 96/(PI*PI*PI*PI)) - static_cast<double>(13)/12;
    const double p1 = c*b;
    const double p2 = c*a;

    // X' = X+D/(2*lambda1), Y' = Y+E/(2*lambda2) <=> X = X'-D/(2*lambda1), Y = Y'-E/(2*lambda2)
    // [x;y] = P[X;Y]

    std::vector<double> Cx, Cy;

    // starting point
    Cx.push_back( a );
    Cy.push_back( 0 );

    // quater 1
    Cx.push_back( a );  Cy.push_back( p1 );
    Cx.push_back( p2 ); Cy.push_back( b );
    Cx.push_back( 0 );  Cy.push_back( b );

    // quater 2
    Cx.push_back( -p2 ); Cy.push_back( b );
    Cx.push_back( -a );  Cy.push_back( p1 );
    Cx.push_back( -a );  Cy.push_back( 0 );

    // quater 3
    Cx.push_back( -a );  Cy.push_back( -p1 );
    Cx.push_back( -p2 ); Cy.push_back( -b );
    Cx.push_back( 0 );   Cy.push_back( -b );

    // quater 4
    Cx.push_back( p2 ); Cy.push_back( -b );
    Cx.push_back( a );  Cy.push_back( -p1 );
    Cx.push_back( a );  Cy.push_back( 0 );

    Cxs.push_back( Cx ); Cys.push_back( Cy );

    return;
}

void trackEllipse( std::vector< std::vector<double> > &Cx, std::vector< std::vector<double> > &Cy,
		   const double *const PHI, const std::vector<double> &zx, const std::vector<double> &zy,
		   const double *const Vx, const double *const Vy )
{
    const int INTERVALS = DEFAULT_P2_INERVALS;
    const double EPS = 1e-10;    
    const double &lambda1 = PHI[0]; const double &lambda2 = PHI[1];
    const double &D = PHI[6]; const double &E = PHI[7]; const double &F = PHI[8];

    assert( lambda1 * lambda2 > 0 );
    assert( fabs(lambda1) + fabs(lambda2) > EPS );
    assert( lambda1 * F < 0 );

    const double &ev1x = PHI[2]; const double &ev1y = PHI[3];
    const double &ev2x = PHI[4]; const double &ev2y = PHI[5];
    const double  P[2][2] = { { ev1x, ev2x }, { ev1y, ev2y } };
    const double PT[2][2] = { { P[0][0], P[1][0] }, { P[0][1], P[1][1] } };

    assert( zx.size() == zy.size() );
#if 1
    std::vector<double> Zx, Zy;
    for(size_t i = 0; i < zx.size(); i++){
	Zx.push_back( PT[0][0]*zx[i] + PT[0][1]*zy[i] + D/(2*lambda1) );
	Zy.push_back( PT[1][0]*zx[i] + PT[1][1]*zy[i] + E/(2*lambda2) );
    }
#endif
    // lambda1 * (X')^2 + lambda2 * (Y')^2 = -F, with lambda1*F < 0 <=> lambda1*(-F) > 0
    // <=> (X')^2 / lambda2 + (Y')^2 / lambda1 = -F/(lambda1*lambda2)
    // <=> X' = sqrt(-F/(lambda1*lambda2))*sqrt(lambda2)*cos(t) = sqrt(-F/lambda1) * cos(t)
    //     Y' = sqrt(-F/)lambda1*lambda2))*sqrt(lambda1)*sin(t) = sqrt(-F/lambda2) * sin(t)
    // <=> cos(t) = X'*sqrt(-lambda1/F), sin(t) = Y'*sqrt(-lambda2/F)

    // lambda1 * X*X + lambda2 * Y*Y + F = 0
    // <=> (-lambda1/F) * X*X + (-lambda2/F) * Y*Y = 1
    // <=> (X/a)^2 + (Y/b)^2 = 1, 1/a = sqrt(-lambda1/F), 1/b = sqrt(-lambda2/F)
    const double a = sqrt( -F / lambda1 );
    const double b = sqrt( -F / lambda2 );

    std::vector<double> theta;
    assert( Zx.size() == Zy.size() );
    for(size_t i = 0; i < Zx.size(); i++)
	theta.push_back( atan2( Zy[i]/b, Zx[i]/a ) );

    std::sort( theta.begin(), theta.end() );

 
    const double PI = atan(static_cast<double>(1))*4;
    theta.push_back( theta[0] + 2*PI );

    for(size_t i = 0; i < theta.size(); i++){

	const double &t0 = theta[i];
	const double &t1 = theta[i+1];

	if( t1-t0 < EPS )
	    continue;

	// examine the arc theta[i]--theta[i+1] is inside triangle or not
	const double mt0 = t0 + (t1-t0)/100;
	const double mx0 = a * cos(mt0);
	const double my0 = b * sin(mt0);

	const double mt1 = t1 - (t1-t0)/100;
	const double mx1 = a * cos(mt1);
	const double my1 = b * sin(mt1);
	if( !(isInsideTriangle( mx0, my0, Vx, Vy ) && isInsideTriangle( mx1, my1, Vx, Vy )) )
	    continue;

	std::vector<double> Cx_local, Cy_local;

	const double dt = (t1 - t0) / INTERVALS;

	Cx_local.push_back( a*cos(t0) );
	Cy_local.push_back( b*sin(t0) );

	for(int k = 0; k < INTERVALS; k++){

	    const double eta1 = t0 + k*dt;
	    const double eta2 = t0 + (k+1)*dt;

	    const double P0x = a * cos( eta1 );
	    const double P0y = b * sin( eta1 );

	    const double P3x = a * cos( eta2 );
	    const double P3y = b * sin( eta2 );
#if 0
	    // (X/a)^2 + (Y/b)^2 = 1
	    // tangent at (X0,Y0) : X0*X/a^2 + Y0*Y/b^2 = 1, whose directional vector is (a*a*Y0,-b*b*X0)
	    const double alpha = 0; // alpha = 0 => segment
#else
	    // L.Maisonobe (2003) Drawing an elliptical arc using polylines, quadratic or cubic Bezier curves
	    // http://www.spaceroots.org/documents/ellipse/elliptical-arc.pdf
	    // or
	    // https://mortoray.com/2017/02/16/rendering-an-svg-elliptical-arc-as-bezier-curves
	    const double te = tan( (eta2-eta1)/2 );
	    const double alpha = sin(eta2-eta1) * (sqrt(4+3*te*te) - 1) / 3;
#endif
	    const double P1x = P0x + alpha * ( -a*sin(eta1) );
	    const double P1y = P0y + alpha * (  b*cos(eta1) );

	    const double P2x = P3x - alpha * ( -a*sin(eta2) );
	    const double P2y = P3y - alpha * (  b*cos(eta2) );

	    Cx_local.push_back( P1x ); Cy_local.push_back( P1y );
	    Cx_local.push_back( P2x ); Cy_local.push_back( P2y );
	    Cx_local.push_back( P3x ); Cy_local.push_back( P3y );
	}
	Cx.push_back( Cx_local );
	Cy.push_back( Cy_local );
    }

    return;
}

void trackHyperbolaCore( std::vector< std::vector<double> > &Cx, std::vector< std::vector<double> > &Cy,
			 const double sign, const double a, const double b, std::vector<double> &x,
			 const double *const Vx, const double *const Vy )
{
    const int INTERVALS = DEFAULT_P2_INERVALS;

    // y = sign * sqrt( a * x^2 + b )

    std::sort( x.begin(),  x.end() );

    std::vector<double> Y;
    for(size_t i = 0; i < x.size(); i++)
	Y.push_back( sign * sqrt( a*x[i]*x[i] + b ) );
	    
    // algorithm in "Mathematical Illustrations" by Bill Casselman,
    // Chapter 7, (2005) Cambridge University Press.
    // https://www.cambridge.org/jp/academic/subjects/mathematics/geometry-and-topology/mathematical-illustrations-manual-geometry-and-postscript?format=PB&isbn=9780521547888
    // See also http://www.math.ubc.ca/~cass/graphics/text/www/

    for(size_t i = 0; i+1 < x.size(); i++){

	std::vector<double> Cx_local, Cy_local;

	const double &X0 = x[i];
	const double &X1 = x[i+1];
	const double h   = (X1 - X0) / INTERVALS;

	const double EPS = 1e-10;
	if( h < EPS )
	    continue;

	// examine the arc X[i]--X[i+1] is inside triangle or not
	const double mX0 = X0 + h/100;
	const double mY0 = sign * sqrt( a*mX0*mX0 + b );

	const double mX1 = X1 - h/100;
	const double mY1 = sign * sqrt( a*mX1*mX1 + b );
	if( !(isInsideTriangle( mX0, mY0, Vx, Vy ) && isInsideTriangle( mX1, mY1, Vx, Vy )) )
	    continue;

	Cx_local.push_back( X0 );
	Cy_local.push_back( Y[i] );
		
	for(int k = 0; k < INTERVALS; k++){

	    const double P0x = X0 + k*h;
	    const double P3x = X0 + (k+1)*h;
	    const double P1x = P0x + h/3;
	    const double P2x = P3x - h/3;

	    const double P0y = sign * sqrt( a*P0x*P0x + b );
	    const double P3y = sign * sqrt( a*P3x*P3x + b );
	    // y = f(x) = sqrt(a*x*x + b)
	    // => f' = 2*a*x / 2 / sqrt(a*x*x+b) = a*x/y
	    const double dY0 = a * P0x / P0y;
	    const double dY3 = a * P3x / P3y;

	    const double P1y = P0y + h*dY0/3;
	    const double P2y = P3y - h*dY3/3;

	    Cx_local.push_back( P1x ); Cy_local.push_back( P1y );
	    Cx_local.push_back( P2x ); Cy_local.push_back( P2y );
	    Cx_local.push_back( P3x ); Cy_local.push_back( P3y );
	}

	Cx.push_back( Cx_local );
	Cy.push_back( Cy_local );
    }

    return;
}

void trackHyperbola( std::vector< std::vector<double> > &Cx, std::vector< std::vector<double> > &Cy,
		     const double *const PHI, const std::vector<double> &zx, const std::vector<double> &zy,
		     const double *const Vx, const double *const Vy )
{
    const double EPS = 1e-10;
    const double &lambda1 = PHI[0]; const double &lambda2 = PHI[1];
    const double &D = PHI[6]; const double &E = PHI[7]; const double &F = PHI[8];

    assert( lambda1 * lambda2 < 0 );
    assert( fabs(lambda1) + fabs(lambda2) > EPS );

    const double &ev1x = PHI[2]; const double &ev1y = PHI[3];
    const double &ev2x = PHI[4]; const double &ev2y = PHI[5];
    const double  P[2][2] = { { ev1x, ev2x }, { ev1y, ev2y } };
    const double PT[2][2] = { { P[0][0], P[1][0] }, { P[0][1], P[1][1] } };

    assert( zx.size() == zy.size() );
#if 1
    std::vector<double> Zx, Zy;
    for(size_t i = 0; i < zx.size(); i++){
	Zx.push_back( PT[0][0]*zx[i] + PT[0][1]*zy[i] + D/(2*lambda1) );
	Zy.push_back( PT[1][0]*zx[i] + PT[1][1]*zy[i] + E/(2*lambda2) );
    }
#endif

    if( lambda1*F > 0 ){

	// lambda1*(X')^2 + lambda2*(Y')^2 + F = 0 <=> Y' = \pm sqrt( (-lambda1*(X')^2 - F)/lambda2 )
	const double a = -lambda1/lambda2; // Y' = \pm sqrt( a(X')^2 + b )
	const double b = -F/lambda2;
	
	std::vector<double> Zx_plus, Zx_minus;
	for(size_t i = 0; i < Zy.size(); i++){
	    if( Zy[i] > 0 ){
		Zx_plus.push_back( Zx[i] );
	    } else {
		Zx_minus.push_back( Zx[i] );
	    }
	}

	// Y' = + sqrt( a(X')^2 + b )
	trackHyperbolaCore( Cx, Cy, +1, a, b, Zx_plus, Vx, Vy );

	// Y' = - sqrt( a(X')^2 + b )
	trackHyperbolaCore( Cx, Cy, -1, a, b, Zx_minus, Vx, Vy );

    } else {

	assert( lambda2*F > 0 );

	// lambda1*(X')^2 + lambda2*(Y')^2 + F = 0 <=> X' = \pm sqrt( (-lambda2*(Y')^2 - F) / lambda1 )
	const double a = -lambda2/lambda1; // X' = \pm sqtt( a(Y')^2 + b )
	const double b = -F/lambda1;
	
	std::vector<double> Zy_plus, Zy_minus;
	for(size_t i = 0; i < Zx.size(); i++){
	    if( Zx[i] > 0 ){
		Zy_plus.push_back( Zy[i] );
	    } else {
		Zy_minus.push_back( Zy[i] );
	    }
	}

	// X' = + sqrt( a(Y')^2 + b )
	trackHyperbolaCore( Cy, Cx, +1, a, b, Zy_plus, Vy, Vx );

	// X' = - sqrt( a(Y')^2 + b )
	trackHyperbolaCore( Cy, Cx, -1, a, b, Zy_minus, Vy, Vx );

    }
    return;
}

void plot_P2_isoline( std::stringstream &Content, const Fem2D::Mesh &Th, const KN<double> &f_P2_,
		      const KNM<double> &palette,
		      const int sizex, const int sizey, const double scale, const double ar,
		      const double x0, const double y0, const double y1,
		      const int marginl, const int marginb,
		      const double textfontsize, const bool monochrome,
		      const bool legend, const int prec, const bool logscale,
		      const double withmesh,
		      const int NISOLINES, const KN<double>*const viso,
		      const double linewidth )
{
    const double EPS = 1e-10;

    //------------------------------
    // values to plot
    //------------------------------
    const double fmax = (viso)? viso->max(): f_P2_.max();
    const double fmin = (viso)? viso->min(): f_P2_.min();
    std::vector<double> isoline_val;

    if( viso ){

	for(int m = 0; m < viso->size(); m++)
	    isoline_val.push_back( (*viso)[m] );

    } else if( logscale && (fmin > 0) ) {
	
	const double df = exp( (static_cast<double>(1)/NISOLINES)*(log(fmax/fmin)) );

	isoline_val.push_back( fmin*sqrt(df) );
	for(int m = 1; m < NISOLINES; m++)
	    isoline_val.push_back( isoline_val[m-1] * df );

    } else {

	if( logscale )
	    std::cout << "plotPDF(): logscale for non-positive values.\n";

#if 1
	const double df = (fmax - fmin) / NISOLINES;
	for(int m = 0; m < NISOLINES; m++)
	    isoline_val.push_back( fmin + df/2 + m*df );
#else
	const double df = (fmax - fmin) / (NISOLINES+1);
	for(int m = 0; m < NISOLINES; m++)
	    isoline_val.push_back( fmin + (m+1)*df );
#endif
    }

    std::stringstream &st = Content;
    st.str("");

    st << "q\n";
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";
    st << "1 w\n"; // setlinewidth

    //------------------------------
    // element (triangle)
    //------------------------------
    if( 0 < withmesh ){

	const double grayscale1 = (withmesh < 1)? withmesh: 1;
	const double grayscale = 1-grayscale1;
	st << grayscale << ' ' << grayscale << ' ' << grayscale << " RG\n";

	for(int n = 0; n < Th.nt; n++){

	    const int &v0 = Th(n,0);
	    const int &v1 = Th(n,1);
	    const int &v2 = Th(n,2);

	    st << scale*ar*(Th(v0).x-x0) << ' ' << scale*(Th(v0).y-y0) << " m ";
	    st << scale*ar*(Th(v1).x-x0) << ' ' << scale*(Th(v1).y-y0) << " l ";
	    st << scale*ar*(Th(v2).x-x0) << ' ' << scale*(Th(v2).y-y0) << " l ";
	    st << "s" << std::endl;
	}

    } // withmesh
    
    st << "Q\n";

    // color/gray-scale range to plot
#if 1
    // If user specifies an irrelevant color range (in viso array),
    // we do not take care of it.
    const double cmax = *max_element( isoline_val.begin(), isoline_val.end() );
    const double cmin = *min_element( isoline_val.begin(), isoline_val.end() );
#else
    // map color/grayscale range to function values
    const double cmax = fmax;
    const double cmin = fmin;
#endif

    //------------------------------
    // main routine
    //------------------------------
    const int NQ = 6; // number of coefficients in quadratic polynomial
    const int &nTriangles = Th.nt;

    for(int k = 0; k < nTriangles; k++){

	st << "q\n";
	st << linewidth << " w\n"; // setlinewidth
	st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";
	
	const int &v0 = Th(k,0); const int &v1 = Th(k,1); const int &v2 = Th(k,2);

	const double vx[] = { Th(v0).x,  Th(v1).x,  Th(v2).x };
	const double vy[] = { Th(v0).y,  Th(v1).y,  Th(v2).y };

	// f_P2[i] = phi(v[i]), f_P2[i+3] = phi(e[i]), i=0,1,2.
	const double *const f_P2 = f_P2_ + k*NQ;

	// quadratic polynomials : ax^2 + bxy + cy^2 + dx + ey + f
	double phi[ NQ ];
	const double &a = phi[0]; const double &b = phi[1]; const double &c = phi[2];
	const double &d = phi[3]; const double &e = phi[4]; const double &f = phi[5];

	findQuadraticPolynomial(phi, vx, vy, f_P2);
	
	const bool isLinear = (fabs(a) + fabs(b) + fabs(c) < EPS * (fabs(d)+fabs(e) + fabs(f)));

	// normalize Quadratic Polynomial
	double PHI[9];
	findCanonicalForm( PHI, phi );

	// PHI(X,Y) = lambda1*X*X + lambda2*Y*Y + D*X + E*Y + f
        //          = lambda1 * (X + D/(2*lambda1))^2 + lambda2 * (Y + E/(2*lambda2))^2
        //            + ( -D*D/(4*lambda1) - E*E/(4*lambda2) + f)
        //          = lambda1*(X + D/(2*lambda1))^2 + lambda2*(Y + E/(2*lambda2))^2 + F
	const double &lambda1 = PHI[0]; const double &lambda2 = PHI[1];
	const double &D = PHI[6];   const double &E = PHI[7];   double &F = PHI[8];

	// lambda1*X^2 + lambda2*Y^2 + D*X + E*Y + f = value,
	// if lambda2 == E == 0, then lambda1*X^2 + D*X + f = value, i.e. X=const
	// if lambda1 == D == 0, then lambda2*Y^2 + E*Y + f = value, i.e. Y=const
	const bool isParallelY = (fabs(lambda2) + fabs(E) < EPS * (fabs(lambda1)+fabs(D) + fabs(F)));
	const bool isParallelX = (fabs(lambda1) + fabs(D) < EPS * (fabs(lambda2)+fabs(E) + fabs(F)));

	// Question: If PHI - isovalue = (a1 X + b1 Y + c1)(a2 X + b2 Y + c2) = 0, what happes? 
 	if( isLinear || isParallelX || isParallelY ){    // phi(x,y) is linear

	    for(size_t m = 0; m < isoline_val.size(); m++){

		const double &value = isoline_val[m];
		
		std::vector<double> px, py;
		trackP1isoline( px, py, vx, vy, value, f_P2 );

		assert( px.size() == py.size() );

		if( px.size() == 0 ) continue; // goto next isoline_val

		setrgbcolor(st, value, palette, cmin, cmax, monochrome, logscale);

		if( px.size() > 3 ){ 

		    // f(x,y)==value on all vertices, i.e. f(x,y) \equiv value on the triangle
		    st << "rg\n";
		    st << scale*ar*(vx[0] - x0) << ' ' << scale*(vy[0] - y0) << " m "
		       << scale*ar*(vx[1] - x0) << ' ' << scale*(vy[1] - y0) << " l "
		       << scale*ar*(vx[2] - x0) << ' ' << scale*(vy[2] - y0) << " l f\n";
		
		} else {

		    // assert( (px.size() == 2) || (px.size() == 3) );
		    st << "RG\n";
		    st << scale*ar*(px[0] - x0) << ' ' << scale*(py[0] - y0) << " m "
		       << scale*ar*(px[1] - x0) << ' ' << scale*(py[1] - y0) << " l S\n";
		}
	    }

	    st << "Q\n";
	    continue; // goto next Triangle
	}

	double Vx[3], Vy[3];
	transformTriangle( Vx, Vy, vx, vy, PHI );

	const bool isParabolic  = (fabs(lambda1) < EPS) || (fabs(lambda2) < EPS);
	const bool isElliptic   = (!isParabolic) && (lambda1*lambda2 > 0);
	const bool isHyperbolic = (!isParabolic) && (lambda1*lambda2 < 0);

	const double &ev1x = PHI[2]; const double &ev1y = PHI[3];
	const double &ev2x = PHI[4]; const double &ev2y = PHI[5];
	const double p[2][2] = { { ev1x, ev2x }, { ev1y, ev2y } };

	for(size_t m = 0; m < isoline_val.size(); m++){

	    const double &value = isoline_val[m];

	    // examine values of phi at all vertices
	    double phi_vertices[3];
	    for(int i = 0; i < 3; i++)
		phi_vertices[i] = a*vx[i]*vx[i] + b*vx[i]*vy[i] + c*vy[i]*vy[i] + d*vx[i] + e*vy[i] + f - value;

	    // If phi == value at all three vertices, then skip
	    // the curve phi==value is ellipse outside the triangle,
	    // or two lines (factrization of hyperbola), which coinside with segments.
	    // BUG: there is exception: two lines case
	    if( fabs(phi_vertices[0]) + fabs(phi_vertices[1]) + fabs(phi_vertices[2]) < EPS )
		continue;

	    std::vector<double> zx, zy;
	    findZeros( zx, zy, vx[0], vy[0], vx[1], vy[1], phi, value );
	    findZeros( zx, zy, vx[1], vy[1], vx[2], vy[2], phi, value );
	    findZeros( zx, zy, vx[2], vy[2], vx[0], vy[0], phi, value );

	    assert( zx.size() == zy.size() );

	    F -= value; // modify constant term in the canonical form
	    
	    setrgbcolor(st, value, palette, cmin, cmax, monochrome, logscale);
	    st << "RG\n";

	    std::vector< std::vector<double> > Cx, Cy; // control points of Bezier curves PHI=0

	    if( isParabolic && (zx.size() >= 2) ){

		trackParabola(Cx, Cy, PHI, zx, zy, Vx, Vy);

	    } else if( isElliptic ){

		if( zx.size() >= 2 ){
		    trackEllipse(Cx, Cy, PHI, zx, zy, Vx, Vy);
		} else {
		    // Ellipse is included inside triangle (might tangent to an edge)
		    // draw whole ellipse
		    trackEllipse(Cx, Cy, PHI, Vx, Vy);
		}
		
	    } else if( isHyperbolic && (zx.size() >= 2) ){

		trackHyperbola(Cx, Cy, PHI, zx, zy, Vx, Vy);
	    }

	    if( Cx.size() > 0 ) {
		invTransformCubicBzeirs( Cx, Cy, PHI );
		drawCubicBeziers( st, Cx, Cy, scale, ar, x0, y0 );
	    }

	    F += value; // recover the canonical form

	} // for isoline_val
	
	st << "Q\n";

    } // for triangle

    st << "q\n";
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";

    if( legend ){

	for(size_t m = 0; m < isoline_val.size(); m++){

	    const double &f = isoline_val[m];

	    setrgbcolor(st, f, palette, cmin, cmax, monochrome, logscale);
	    st << "rg\n";

	    st << "BT /F1 " << textfontsize << " Tf "
	       << "1 0 0 1 " << sizex-PADDING << " " << (m+1)*(scale*(y1-y0)-textfontsize)/(isoline_val.size()+1) << " Tm "
	       << "(" << ((f >= 0)? "\\ ": "");

	    if( (fabs(f) > 1e-3) || (fabs(f) < 1e-12) ){
		st << std::setprecision(prec) << std::setfill('0') << f << ") Tj ET\n";
	    } else {
		st << std::resetiosflags(std::ios::fixed) << std::setiosflags(std::ios::scientific)
		   << std::setprecision(prec) << f << resetiosflags(std::ios::scientific)
		   << std::setiosflags(std::ios::fixed) << ") Tj ET\n";
	    }
	    // debug: st << std::setiosflags(std::ios::fixed)
	}
    } // legend

    //------------------------------
    // edges
    //------------------------------
    const int &nEdges = Th.neb;

    st << "0 0 0 RG\n";

    for(int k = 0; k < nEdges; k++){

	const int &v0 = Th( Th.bedges[k][0] );
	const int &v1 = Th( Th.bedges[k][1] );

	st << scale*ar*(Th(v0).x-x0) << ' ' << scale*(Th(v0).y-y0) << " m ";
	st << scale*ar*(Th(v1).x-x0) << ' ' << scale*(Th(v1).y-y0) << " l s" << std::endl;
    }

    st << "Q\n";
    return;
} //

bool isSegment( const std::vector<double> &cx, const std::vector<double> &cy, const int i )
{
    // examine a part of cubic Bezier curve (i--i+1--i+2--i+3) is a segment or not.
    return (cx[i] == cx[i+1]) && (cy[i] == cy[i+1]) && (cx[i+2] == cx[i+3]) && (cy[i+2] == cy[i+3]);
}

int findSegment( const double x, const double y,
		 const std::vector<double> &cx, const std::vector<double> &cy )
{
    const double EPS = 1e-10;

    // cx and cy are control points of cubic Bezier curve.
    // (cx[0],cy[0]) -- (cx[1],cy[1]) -- (cx[2],cy[2]) -- (cx[3],cy[3])  : 0
    //               -- ...
    //               -- (cx[3i+1],cy[3i+1]) -- (cx[3i+2],cy[3i+2]) -- (cx[3i+3],cy[3i+3]) : 3i

    for(size_t i = 0; i+3 < cx.size(); i += 3){

	if( !isSegment( cx, cy, i ) )
	    continue;

	const double &x0 = cx[i];   const double &y0 = cy[i];
	const double &x1 = cx[i+3]; const double &y1 = cy[i+3];
	
	// if (x,y) \in (x0,y0) -- (x1,y1), then
	// x = (1-t)*x0 + t*x1 and y = (1-t)*y0 + t*y1, 0 <= t <= 1
	// <=> x = x0 - t*x0 + t*x1, y = y0 - t*y0 + t*y1
	// <=> t = (x-x0) / (x1-x0) = (y-y0) / (y1-y0)

	if( fabs( (x-x0)*(y1-y0) - (x1-x0)*(y-y0) ) > EPS )
	    continue;

	if( fabs( x1-x0 ) > EPS ){
	    const double t = (x-x0) / (x1-x0);

	    if( (-EPS < t) && (t < 1+EPS) )
		return i;
	}

	if( fabs( y1-y0 ) > EPS ){
	    const double t = (y-y0) / (y1-y0);
	    
	    if( (-EPS < t) && (t < 1+EPS) )
		return i;
	}
    }

    return -1;
}

void splitByBorder( std::vector< std::vector<double> > &partition_x, std::vector< std::vector<double> > &partition_y,
		    const std::vector< std::vector<double> > &cxs, const std::vector< std::vector<double> > &cys )
{
    assert( partition_x.size() == partition_y.size() );
    assert( cxs.size() == cys.size() );

    for(size_t i = 0; i < cxs.size(); i++){

	std::vector<double> cx = cxs[i];
	std::vector<double> cy = cys[i];

	assert( cx.size() == cy.size() );

	const double &beg_x = cx[0];
	const double &beg_y = cy[0];
    
	const double &end_x = cx.back(); // cx[ cx.size()-1 ]
	const double &end_y = cy.back(); // cy[ cy.size()-1 ]

	for(size_t j = 0; j < partition_x.size(); j++){

	    const std::vector<double> &px = partition_x[j];
	    const std::vector<double> &py = partition_y[j];

	    assert( px.size() == py.size() );

	    int beg_id = findSegment( beg_x, beg_y, px, py );

	    if( beg_id < 0 ) continue;

	    int end_id = findSegment( end_x, end_y, px, py );

	    if( end_id < 0 ) continue;

	    // findSegment returns 3n ( multiple of 3 )

	    if( end_id < beg_id ){
		int tmp_id = beg_id;
		beg_id = end_id;
		end_id = tmp_id;
		std::reverse(cx.begin(), cx.end());
		std::reverse(cy.begin(), cy.end());
	    }

	    if( beg_id == end_id ){

		if( (px[beg_id]-beg_x)*(px[beg_id]-beg_x) + (py[beg_id]-beg_y)*(py[beg_id]-beg_y)
		    > (px[beg_id]-end_x)*(px[beg_id]-end_x) + (py[beg_id]-end_y)*(py[beg_id]-end_y) ){

		    std::reverse(cx.begin(), cx.end());
		    std::reverse(cy.begin(), cy.end());
		}
	    }

	    // begin and end points are on segment beg_id-component and end_id-component respectively,
	    // which are parts of boundary of partition.
	    // partition_x, partition_y are represented as a cubic Bezier curve

	    std::vector<double> p0x, p0y;

	    for(size_t k = 0; k <= beg_id+1; k++){
		p0x.push_back( px[k] );
		p0y.push_back( py[k] );
	    }

	    p0x.push_back( cx[0] ); p0y.push_back( cy[0] );

	    for(size_t k = 0; k < cx.size(); k++){
		p0x.push_back( cx[k] );
		p0y.push_back( cy[k] );
	    }

	    p0x.push_back( cx.back() ); p0y.push_back( cy.back() );

	    for(size_t k = end_id+2; k < px.size(); k++){
		p0x.push_back( px[k] );
		p0y.push_back( py[k] );
	    }

	    std::vector<double> p1x, p1y;

	    p1x.push_back( cx[0] ); p1y.push_back( cy[0] );
	    p1x.push_back( cx[0] ); p1y.push_back( cy[0] );

	    for(size_t k = beg_id+2; k <= end_id+1; k++){
		p1x.push_back( px[k] );
		p1y.push_back( py[k] );
	    }

	    p1x.push_back( cx.back() ); p1y.push_back( cy.back() );

	    for(int k = cx.size()-1; k >= 0; k--){
		p1x.push_back( cx[k] );
		p1y.push_back( cy[k] );
	    }

	    partition_x.erase( partition_x.begin() + j );
	    partition_y.erase( partition_y.begin() + j );

	    partition_x.push_back( p0x ); partition_y.push_back( p0y );
	    partition_x.push_back( p1x ); partition_y.push_back( p1y );

	    break;
	}
    }

    return;
}

double findFillValue( const std::vector<double> &cx, const std::vector<double> &cy, const double *const phi )
{
    const double &a = phi[0]; const double &b = phi[1]; const double &c = phi[2];
    const double &d = phi[3]; const double &e = phi[4]; const double &f = phi[5];

    int npoints = 0;

    // find average of phi at nodal points
    double sum = 0;
    for(size_t i = 0; i < cx.size(); i += 3){

	const double &x = cx[i];
	const double &y = cy[i];

	sum += (a*x*x + b*x*y + c*y*y + d*x + e*y + f);
	npoints++;

	if( (i+3 < cx.size()) && isSegment(cx,cy,i) ){

	    const double mx = ( cx[i]+cx[i+3] ) / 2;
	    const double my = ( cy[i]+cy[i+3] ) / 2;
	    sum += (a*mx*mx + b*mx*my + c*my*my + d*mx + e*my + f);
	    npoints++;
	}
    }

    return sum / npoints;
}

void P2_fill_linear( std::stringstream &Content,
		     const Fem2D::Mesh &Th, const int triangle_id,
		     const double *const vf,
		     const KNM<double> &palette,
		     const double fmax, const double fmin, const double df,
		     const int nbfill, const KN<double> *const frange,
		     const double scale, const double ar, const double x0, const double y0,
		     const bool monochrome, const bool logscale )
{
    std::stringstream &st = Content;

    const int &v0 = Th(triangle_id,0);
    const int &v1 = Th(triangle_id,1);
    const int &v2 = Th(triangle_id,2);
    const double x[] = { Th(v0).x, Th(v1).x, Th(v2).x };
    const double y[] = { Th(v0).y, Th(v1).y, Th(v2).y };
    // phi(x,y) = a*x + b*y + c
    const double *const &phi = vf;

    // find minimum
    int v_min = 2;
    int end1 = 0; // destination point
    int end2 = 1; // destination point

    if( (phi[0] <= phi[1]) && (phi[0] <= phi[2]) ){
	v_min = 0;
	end1 = 1;
	end2 = 2;
    } else if( (phi[1] <= phi[0]) && (phi[1] <= phi[2]) ){
	v_min = 1;
	end1 = 2;
	end2 = 0;
    }

    if( frange && (fmax < phi[v_min]) )
	return;

    if( frange && (phi[end1] < fmin) && (phi[end2] < fmin) )
	return;

    int beg1 = v_min; // current point
    int beg2 = v_min; // current point
    double t1 = 0;    // secting parameters on edge beg1--end1
    double t2 = 0;    // secting parameters on edge beg2--end2

    int level = (logscale)?
	static_cast<int>( log(phi[v_min]/fmin) / log(df) ):
	static_cast<int>( (phi[v_min]-fmin) / df );
    
    double f = (logscale)? (pow(df,level) * fmin) : (level*df + fmin);

    if( frange && (phi[v_min] < fmin) ){

	if( phi[end1] < fmin ){
	    beg1 = end1;
	    end1 = end2;
	} if( phi[end2] < fmin ){
	    beg2 = end2;
	    end2 = end1;
	}

	t1 = (fmin - phi[beg1]) / (phi[end1] - phi[beg1]);
	t2 = (fmin - phi[beg2]) / (phi[end2] - phi[beg2]);
	level = 0;
	f = fmin;
    }

    do {
	const double EPS = 1e-10;
	if( (end1 == end2) && (t1 >= 1-EPS) && (t2 >= 1-EPS) )
	    break;

	if( (end1 == end2) && (t1 >= 1-EPS) && (t2 >= 1-EPS) )
	    break;

	// find points which divide edges by t:(1-t1), t2:(1-t2) respectively
	const double p1x = (1-t1)*x[beg1] + t1*x[end1];
	const double p1y = (1-t1)*y[beg1] + t1*y[end1];
	const double p2x = (1-t2)*x[beg2] + t2*x[end2];
	const double p2y = (1-t2)*y[beg2] + t2*y[end2];

	const double f_next = (logscale)? f * df: f + df;
	
	if( frange && (fmax-EPS <= f) )
	    break;

	if( level == 0 ){
	    setrgbcolor(st, fmin, palette, fmin, fmax, monochrome, logscale);
	} else if( level >= nbfill-1 ){
		
	    // exceptional handling : if local minimum (in triangle) is
	    // already global maximum, then level == nbfill.
	    // Thus the condition is not level == nbfill-1,
	    // but level >= nbfill-1.
	    setrgbcolor(st, fmax, palette, fmin, fmax, monochrome, logscale);

	} else {
	    const double c = (logscale)?
		(pow(df,level+0.5)*fmin): ((level+0.5)*df + fmin);
	    setrgbcolor(st, c, palette, fmin, fmax, monochrome, logscale);
	}
	st << "rg\n";
		
	if( (f_next >= phi[end1]) && (f_next >= phi[end2]) ){

	    st << scale*ar*(p1x-x0)        << ' ' << scale*(p1y-y0)        << " m "
	       << scale*ar*(x[end1]-x0) << ' ' << scale*(y[end1]-y0) << " l "
	       << scale*ar*(x[end2]-x0) << ' ' << scale*(y[end2]-y0) << " l "
	       << scale*ar*(p2x-x0)        << ' ' << scale*(p2y-y0)        << " l f\n";

	    break;
	}

	if( f_next >= phi[end1] ){
		  
	    st << scale*ar*(p1x-x0)        << ' ' << scale*(p1y-y0)        << " m "
	       << scale*ar*(x[end1]-x0) << ' ' << scale*(y[end1]-y0) << " l "
	       << scale*ar*(p2x-x0)        << ' ' << scale*(p2y-y0)        << " l f\n";

	    beg1 = end1;
	    end1 = end2;
	    t1 = 0;
	    continue;
	} 

	if( f_next >= phi[end2] ){

	    st << scale*ar*(p1x-x0) << ' ' << scale*(p1y-y0) << " m "
	       << scale*ar*(x[end2]-x0) << ' ' << scale*(y[end2]-y0) << " l "
	       << scale*ar*(p2x-x0) << ' ' << scale*(p2y-y0) << " l f\n";
		  
	    beg2 = end2;
	    end2 = end1;
	    t2 = 0;
	    continue;
	}

	const double t1_next = (f_next - phi[beg1]) / (phi[end1] - phi[beg1]);
	const double t2_next = (f_next - phi[beg2]) / (phi[end2] - phi[beg2]);

	const double q1x = (1-t1_next)*x[beg1] + t1_next*x[end1];
	const double q1y = (1-t1_next)*y[beg1] + t1_next*y[end1];
	const double q2x = (1-t2_next)*x[beg2] + t2_next*x[end2];
	const double q2y = (1-t2_next)*y[beg2] + t2_next*y[end2];

	st << scale*ar*(p1x-x0) << ' ' << scale*(p1y-y0) << " m "
	   << scale*ar*(q1x-x0) << ' ' << scale*(q1y-y0) << " l "
	   << scale*ar*(q2x-x0) << ' ' << scale*(q2y-y0) << " l "
	   << scale*ar*(p2x-x0) << ' ' << scale*(p2y-y0) << " l f\n";

	// update
	f = f_next;
	level++;
	t1 = t1_next;
	t2 = t2_next;

    } while(level <= nbfill);

    return;
}

void plot_P2_fill( std::stringstream &Content, const Fem2D::Mesh &Th, const KN<double> &f_P2_,
		   const KNM<double> &palette,
		   const int sizex, const int sizey, const double scale, const double ar,
		   const double x0, const double y0, const double y1,
		   const int marginl, const int marginb,
		   const double textfontsize, const bool monochrome,
		   const bool legend, const int prec, const bool logscale,
		   const double withmesh,
		   const long nbfill, const KN<double> *const frange )
{
    const double EPS = 1e-10;

    //------------------------------
    // borders to plot
    //------------------------------
    const double tmp_fmax = (frange)? (*frange)[1]: f_P2_.max();
    const double tmp_fmin = (frange)? (*frange)[0]: f_P2_.min();
    const double tmp_df = (logscale)?
	pow( tmp_fmax/tmp_fmin, static_cast<double>(1)/nbfill ):
	(tmp_fmax - tmp_fmin)/nbfill;

    std::vector<double> border_val;
    if( logscale && (tmp_fmin > 0) ){

	border_val.push_back( tmp_fmin );
	for(int m = 1; m <= nbfill; m++)
	    border_val.push_back( border_val[m-1] * tmp_df );

    } else {

	if( logscale )
	    std::cout << "plotPDF(): logscale for non-positive values.\n";

	for(int m = 0; m <= nbfill; m++)
	    border_val.push_back( m*tmp_df + tmp_fmin );
    }

    // considering wider range, since min/max of P2 element are not
    // attained at vertices of triangles.
    const double fmax = (logscale)? tmp_fmax*tmp_df: tmp_fmax + tmp_df;
    const double fmin = (logscale)? tmp_fmin/tmp_df: tmp_fmin - tmp_df;
    
    std::stringstream &st = Content;
    st.str("");

    st << "q\n";
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";

    //------------------------------
    // element(triangle)-wise process
    //------------------------------
    const int NQ = 6; // number of unknowns in quadratic polynomial
    const int &nTriangles = Th.nt;
    
    for(int k = 0; k < nTriangles; k++){

	const int &v0 = Th(k,0); const int &v1 = Th(k,1); const int &v2 = Th(k,2);

	const double vx[] = { Th(v0).x,  Th(v1).x,  Th(v2).x };
	const double vy[] = { Th(v0).y,  Th(v1).y,  Th(v2).y };

	// f_P2[i] = phi(v[i]), f_P2[i+3] = phi(e[i]), i=0,1,2.
	const double *const f_P2 = f_P2_ + k*NQ;

	// quadratic polynomials : ax^2 + bxy + cy^2 + dx + ey + f
	double phi[ NQ ];
	const double &a = phi[0]; const double &b = phi[1]; const double &c = phi[2];
	const double &d = phi[3]; const double &e = phi[4]; const double &f = phi[5];

	findQuadraticPolynomial(phi, vx, vy, f_P2);
	
	const bool isLinear = (fabs(a) + fabs(b) + fabs(c) < EPS * (fabs(d)+fabs(e) + fabs(f)));

	// normalize Quadratic Polynomial
	double PHI[9];
	findCanonicalForm( PHI, phi );

	// PHI(X,Y) = lambda1*X*X + lambda2*Y*Y + D*X + E*Y + f
        //          = lambda1 * (X + D/(2*lambda1))^2 + lambda2 * (Y + E/(2*lambda2))^2
        //            + ( -D*D/(4*lambda1) - E*E/(4*lambda2) + f)
        //          = lambda1*(X + D/(2*lambda1))^2 + lambda2*(Y + E/(2*lambda2))^2 + F
	const double &lambda1 = PHI[0]; const double &lambda2 = PHI[1];
	const double &D = PHI[6];   const double &E = PHI[7];   double &F = PHI[8];

	// lambda1*X^2 + lambda2*Y^2 + D*X + E*Y + F = value,
	// if lambda2 == E == 0, then lambda1*X^2 + D*X + F = value, i.e. X=const
	// if lambda1 == D == 0, then lambda2*Y^2 + E*Y + F = value, i.e. Y=const
	const bool isParallelY = (fabs(lambda2) + fabs(E) < EPS * (fabs(lambda1)+fabs(D) + fabs(F)));
	const bool isParallelX = (fabs(lambda1) + fabs(D) < EPS * (fabs(lambda2)+fabs(E) + fabs(F)));

	// Question: If PHI - isovalue = (a1 X + b1 Y + c1)(a2 X + b2 Y + c2) = 0, what happes? 
	if( isLinear || isParallelX || isParallelY ){    // phi(x,y) is linear

	    P2_fill_linear( st, Th, k, f_P2, palette,
			    tmp_fmax, tmp_fmin, tmp_df, nbfill, frange,
			    scale, ar, x0, y0, monochrome, logscale );

	    continue; // goto next Triangle
	}

	double Vx[3], Vy[3];
	transformTriangle( Vx, Vy, vx, vy, PHI );

	const bool isParabolic  = (fabs(lambda1) < EPS) || (fabs(lambda2) < EPS);
	const bool isElliptic   = (!isParabolic) && (lambda1*lambda2 > 0);
	const bool isHyperbolic = (!isParabolic) && (lambda1*lambda2 < 0);

	const double &ev1x = PHI[2]; const double &ev1y = PHI[3];
	const double &ev2x = PHI[4]; const double &ev2y = PHI[5];
	const double p[2][2] = { { ev1x, ev2x }, { ev1y, ev2y } };

	const std::vector<double> Tx { Vx[0], Vx[0], Vx[1], Vx[1],
					      Vx[1], Vx[2], Vx[2],
					      Vx[2], Vx[0], Vx[0] };
	const std::vector<double> Ty { Vy[0], Vy[0], Vy[1], Vy[1],
					      Vy[1], Vy[2], Vy[2],
					      Vy[2], Vy[0], Vy[0] };

	std::vector< std::vector<double> > partition_x { Tx };
	std::vector< std::vector<double> > partition_y { Ty };

	// divide T by isolines
	for(size_t m = 0; m < border_val.size(); m++){

	    const double &value = border_val[m];

	    // examine values of phi at all vertices
	    double phi_vertices[3];
	    for(int i = 0; i < 3; i++)
		phi_vertices[i] = a*vx[i]*vx[i] + b*vx[i]*vy[i] + c*vy[i]*vy[i] + d*vx[i] + e*vy[i] + f - value;

	    // If phi == value at all three vertices, then skip
	    // the curve phi==value is ellipse outside the triangle,
	    // or two lines (factrization of hyperbola), which coinside with segments.
	    // BUG: there is exception: two lines case
	    if( fabs(phi_vertices[0]) + fabs(phi_vertices[1]) + fabs(phi_vertices[2]) < EPS )
		continue;

	    std::vector<double> zx, zy;
	    findZeros( zx, zy, vx[0], vy[0], vx[1], vy[1], phi, value );
	    findZeros( zx, zy, vx[1], vy[1], vx[2], vy[2], phi, value );
	    findZeros( zx, zy, vx[2], vy[2], vx[0], vy[0], phi, value );

	    assert( zx.size() == zy.size() );

	    F -= value; // modify constant term

	    std::vector< std::vector<double> > Cx, Cy; // control points of Bezier curves PHI=0

	    if( isParabolic && (zx.size() >= 2) ){

		trackParabola(Cx, Cy, PHI, zx, zy, Vx, Vy);

	    } else if( isElliptic ){

		if( zx.size() >= 2 ){
		    trackEllipse(Cx, Cy, PHI, zx, zy, Vx, Vy);
		} else {
		    // Ellipse is included inside triangle (might tangent to an edge)
		    // draw whole ellipse
		    trackEllipse(Cx, Cy, PHI, Vx, Vy);
		}
		
	    } else if( isHyperbolic && (zx.size() >= 2) ){

		trackHyperbola(Cx, Cy, PHI, zx, zy, Vx, Vy);
	    }

	    if( Cx.size() > 0 ) 
		splitByBorder( partition_x, partition_y, Cx, Cy );

	    F += value; // correct constant term

	} // for each border value

	assert( partition_x.size() == partition_y.size() );

	invTransformCubicBzeirs( partition_x, partition_y, PHI );

	// draw divided partitions
	for(size_t i = 0; i < partition_x.size(); i++){

	    const std::vector<double> &cx = partition_x[i];
	    const std::vector<double> &cy = partition_y[i];

	    const double fillvalue = findFillValue( cx, cy, phi );
	    double fm = (logscale)?
		tmp_fmin*pow(tmp_df, static_cast<int>( log(fillvalue/tmp_fmin) / log(tmp_df) ) + 0.5):
		floor((fillvalue - fmin)/tmp_df) * tmp_df + fmin + tmp_df/2;

	    // Since max/min of P2 approx cannot be obtained from only vertex values,
	    // the next is not enough for P2 element.
	    //if( (value < fmin) || (fmax < value) )
	    //    continue; // goto next sub-triangle
	    // Following solves the issue.
	    if( logscale ) {

		if( fillvalue < tmp_fmin ){
		    if(frange) continue; // No color, goto next partition
		    fm = tmp_fmin; // if range is not specified, truncate fill color
		} else if( fillvalue < tmp_fmin*tmp_df ){
		    fm = tmp_fmin;
		}
		
		if( fillvalue > tmp_fmax ){
		    if(frange) continue; // No color, goto next partition
		    fm = tmp_fmax; // if range is not specified, truncate fill color
		} else if( fillvalue > tmp_fmax/tmp_df ){
		    fm = tmp_fmax;
		}

	    } else { // if( logscale )

		if( fillvalue < tmp_fmin ) {
		    if(frange) continue; // No color, goto next partition
		    fm = tmp_fmin; // if range is not specified, truncate fill color
		} else if( fillvalue < tmp_fmin+tmp_df ){
		    fm = tmp_fmin;
		}
		
		if( fillvalue > tmp_fmax ){
		    if(frange) continue; // No color, goto next partition
		    fm = tmp_fmax; // if range is not specified, truncate fill color
		} else if( fillvalue > tmp_fmax-tmp_df ){
		    fm = tmp_fmax;
		}

	    } // if(logscale)

 	    setrgbcolor(st, fm, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
#if 1
	    st << "rg\n";
	    
	    // assume : curve begin point and end point are conisnde.
	    st << scale*ar*(cx[0]-x0) << ' ' << scale*(cy[0]-y0) << " m ";
	    for(size_t i = 1; i < cx.size(); i += 3){
		st << scale*ar*(cx[i+0]-x0) << ' ' << scale*(cy[i+0]-y0) << ' '
		   << scale*ar*(cx[i+1]-x0) << ' ' << scale*(cy[i+1]-y0) << ' '
		   << scale*ar*(cx[i+2]-x0) << ' ' << scale*(cy[i+2]-y0) << " c ";
	    }
	    // f: closing path and fill (without boundary), b: unclosing path and fill (with boundary)
	    st << "f\n";
#else
	    // debug : drawing the borders
	    st << "RG\n";

	    st << scale*ar*(cx[0]-x0) << ' ' << scale*(cy[0]-y0) << " m ";
	    for(size_t i = 1; i < cx.size(); i += 3){
		st << scale*ar*(cx[i+0]-x0) << ' ' << scale*(cy[i+0]-y0) << ' '
		   << scale*ar*(cx[i+1]-x0) << ' ' << scale*(cy[i+1]-y0) << ' '
		   << scale*ar*(cx[i+2]-x0) << ' ' << scale*(cy[i+2]-y0) << " c ";
	    }
	    // s: closing path and stroke, S: unclosing path and stroke
	    st << "S\n";
#endif
	}

    } // element(triangle)-wise process

    //------------------------------
    // legend
    //------------------------------
    if( legend ){

	const double dy = scale*(y1-y0)/nbfill;

	for(int m = 0; m < border_val.size()-1; m++){
#if 0
	    const double f = (logscale)?
		sqrt( border_val[m] * border_val[m+1] ):
		(border_val[m] + border_val[m+1]) / 2;

	    setrgbcolor(st, f, palette, fmin, fmax, monochrome, logscale);
#else
	    if( m == 0 ){
		setrgbcolor(st, tmp_fmin, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
	    } else if( m == nbfill-1 ){
		setrgbcolor(st, tmp_fmax, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
	    } else {
		const double f = (logscale)? tmp_fmin * pow(tmp_df,m+0.5): tmp_fmin + (m+0.5)*tmp_df;
		setrgbcolor(st, f, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
	    }
#endif
	    st << "rg\n";

	    st << sizex-PADDING   << " " << m*dy << " m "
               << sizex-PADDING/2 << " " << m*dy << " l "
               << sizex-PADDING/2 << " " << (m+1)*dy << " l "
	       << sizex-PADDING   << " " << (m+1)*dy << " l f\n";
	}

	const double EPS = 1e-10;

        const double dl = (logscale)?
	    pow( tmp_fmax/tmp_fmin, static_cast<double>(1)/(NUM_LABELS-1) ):
	    (tmp_fmax-tmp_fmin)/(NUM_LABELS-1);

        for(int m = 0; m < NUM_LABELS; m++){
#if 0
            const double f = (logscale)? tmp_fmin * pow(dl,m): tmp_fmin + m*dl;
	    if( nbfill < NUM_LABELS ){
		const double ff = (logscale)?
		    tmp_fmin * pow(tmp_df, static_cast<int>(log(f/tmp_fmin)/log(tmp_df)) ):
		    static_cast<int>((f-tmp_fmin)/tmp_df) * tmp_df + tmp_fmin;
		setrgbcolor(st, ff, palette, fmin, fmax, monochrome, logscale);
	    } else {
		setrgbcolor(st, f, palette, fmin, fmax, monochrome, logscale);
	    }
#else
	    const double f = (logscale)? tmp_fmin * pow(dl,m): tmp_fmin + m*dl;
	    if( logscale ){

		if( f <= tmp_fmin*tmp_df ){
		    setrgbcolor(st, tmp_fmin, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
		} else if( f >= tmp_fmax/tmp_df ){
		    setrgbcolor(st, tmp_fmax, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
		} else {
		    const double dc = pow( tmp_fmax/tmp_fmin, static_cast<double>(1)/nbfill );
		    const int mc = static_cast<int>( log(f/tmp_fmin) / log(dc) );
		    const double c = tmp_fmin * pow(dc, mc+0.5);
		    setrgbcolor(st, c, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
		}

	    } else {

		if( f <= tmp_fmin+tmp_df ){
		    setrgbcolor(st, tmp_fmin, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
		} else if( f >= tmp_fmax-tmp_df ){
		    setrgbcolor(st, tmp_fmax, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
		} else {
		    const double dc = (tmp_fmax-tmp_fmin) / nbfill;
		    const int mc = static_cast<int>( (f-tmp_fmin)/dc );
		    const double c = tmp_fmin + (mc+0.5)*dc;
		    setrgbcolor(st, c, palette, tmp_fmin, tmp_fmax, monochrome, logscale);
		}
	    }
#endif
            st << " rg\n";

            st << "BT /F1 " << textfontsize << " Tf "
	       << "1 0 0 1 " << sizex << " " << m*(scale*(y1-y0)-textfontsize)/(NUM_LABELS-1) << " Tm "
	       << "(" << ((f >= 0)? "\\ ": "");

	    if( (fabs(f) > 1e-3) || (fabs(f) < 1e-12) ){
		st << std::setprecision(prec) << std::setfill('0') << f << ") Tj ET\n";
	    } else {
		st << std::resetiosflags(std::ios::fixed) << std::setiosflags(std::ios::scientific)
		   << std::setprecision(prec) << f << resetiosflags(std::ios::scientific)
		   << std::setiosflags(std::ios::fixed) << ") Tj ET\n";
	    }
	    // debug: st << std::setiosflags(std::ios::fixed)
        }
    } // legend

    //------------------------------
    // element (triangle)
    //------------------------------
    if( 0 < withmesh ){

	const double grayscale1 = (withmesh < 1)? withmesh: 1;
	const double grayscale = 1-grayscale1;
	st << grayscale << ' ' << grayscale << ' ' << grayscale << " RG\n";

	for(int n = 0; n < Th.nt; n++){

	    const int &v0 = Th(n,0);
	    const int &v1 = Th(n,1);
	    const int &v2 = Th(n,2);

	    st << scale*ar*(Th(v0).x-x0) << ' ' << scale*(Th(v0).y-y0) << " m ";
	    st << scale*ar*(Th(v1).x-x0) << ' ' << scale*(Th(v1).y-y0) << " l ";
	    st << scale*ar*(Th(v2).x-x0) << ' ' << scale*(Th(v2).y-y0) << " l ";
	    st << "s" << std::endl;
	}
    } // withmesh

    //------------------------------
    // edges
    //------------------------------
    const int &nEdges = Th.neb;

    st << "0 0 0 RG\n";

    for(int k = 0; k < nEdges; k++){

	const int &v0 = Th( Th.bedges[k][0] );
	const int &v1 = Th( Th.bedges[k][1] );

	st << scale*ar*(Th(v0).x-x0) << ' ' << scale*(Th(v0).y-y0) << " m "
	   << scale*ar*(Th(v1).x-x0) << ' ' << scale*(Th(v1).y-y0) << " l S\n";
    }

    st << "Q\n";
    return;
} //

//----------------------------------------------------------------------
// vector field
//----------------------------------------------------------------------

void plot_vector( std::stringstream &Content,
		  const double arrow_origin_x, const double arrow_origin_y,
		  const double vector_x, const double vector_y, const double vector_norm2,
		  const double alength_scale, const double ahead_scale,
		  const double x0, const double y0, const double scale, const double ar,
		  const double fmin, const double fmax,
		  const bool unit_arrow, const bool logscale,
		  const KNM<double> &palette, const bool monochrome )
{
    std::stringstream &st = Content;
    const double &r       = scale;

    const double &ox = arrow_origin_x;
    const double &oy = arrow_origin_y;

    const double AH_SIZE = (alength_scale > 0)? (DEF_ARROW_HEAD_SIZE * ahead_scale): (DEF_ARROW_HEAD_SIZE * (-ahead_scale));
    const double &AH_ANGLE = ARROW_HEAD_ANGLE;

    // what happen if (fmin == 0) && logscale ?
    const double favg = logscale? (sqrt(fmax*fmin)): (fmax+fmin)/2;

    // scaled arrow length is arrow_scale * cf2, which is truncated if greater than fmax
    const double &as = alength_scale;
    const double alength = (unit_arrow)? (as*favg)/fmax*MAX_ARROW_LENGTH:
	( (logscale)? as*(log(vector_norm2/fmin))/(log(fmax/fmin))*MAX_ARROW_LENGTH: as*vector_norm2/fmax*MAX_ARROW_LENGTH);

    // In logscale, vector_norm2 = fmin*(fmax/fmin)^r <=> r = (log(vector_norm2)-log(fmin))/(log(fmax)-log(fmin))
    // Then, alength = alength_scale * r * MAX_ARROW_LENGTH.

    const double arrow_head_x = r*ar*(ox-x0) + alength * ar*vector_x/vector_norm2;
    const double arrow_head_y = r*(oy-y0) + alength * vector_y/vector_norm2;

    setrgbcolor(st, vector_norm2, palette, fmin, fmax, monochrome, logscale);
    st << "RG\n";
    st << r*ar*(ox-x0) << ' ' << r*(oy-y0) << " m "; // arrow origin
    st << arrow_head_x << ' ' << arrow_head_y << " l S" << std::endl; // arrow head

    // Need fabs? Yes, if coef = alength < 0, arrows go reverse direction.
    if( fabs(alength) > AH_SIZE/2 ){
	    
	const double theta = atan2( -vector_y, -vector_x ); // reverse of arrow direction
	
	st << arrow_head_x + AH_SIZE*cos(theta-AH_ANGLE) << ' ' << arrow_head_y + AH_SIZE*sin(theta-AH_ANGLE) << " m "
	   << arrow_head_x << ' ' << arrow_head_y << " l "
	   << arrow_head_x + AH_SIZE*cos(theta+AH_ANGLE) << ' ' << arrow_head_y + AH_SIZE*sin(theta+AH_ANGLE) << " l S" << std::endl;
    }

    return;
}

void plot_vector2flow( std::stringstream &Content, const Fem2D::Mesh &Th,
		       const KN<double> &fx, const KN<double> &fy, const KN<double> &f2, 
		       const bool fromVertex, const KNM<double> &palette,
		       const double alength_scale, const bool unit_arrow,
		       const double linewidth, const double ahead_scale,
		       const int nbarrow, const KN<double> *const varrow,
		       const int sizex, const int sizey, const double scale, const double ar,
		       const double x0, const double y0, const double y1,
		       const int marginl, const int marginb,
		       const double textfontsize, const bool monochrome,
		       const bool legend, const int prec, const bool logscale,
		       const double withmesh,
		       const long nbfill, const KN<double> *const frange, const string &fetype,
		       const bool isoline, const int NISOLINES, const KN<double>*const viso )
{
    const double EPS = 1e-10;

    const int &nVertices  = Th.nv;
    const int &nTriangles = Th.nt;
    const int &nEdges     = Th.neb;
    const double &r       = scale;

    const double fmax = (frange)? (*frange)[1]: f2.max();
    const double fmin = (frange)? (*frange)[0]: f2.min();

    const double df   = (logscale)?
	exp( (static_cast<double>(1)/nbfill)*(log(fmax/fmin)) ):
	(fmax - fmin)/nbfill;

    if( varrow )
	std::cout << "plotPDF(): Option 'varrow' is not implmeneted yet." << endl;
    if( nbarrow != 0 )
	std::cout << "plotPDF(): Option 'nbarrow' is not implmeneted yet." << endl;

    std::stringstream &st = Content;
    st.str("");

    //------------------------------
    // element (triangle)
    //------------------------------
    if( 0 < withmesh ){

	st << "q\n";
	st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";

	const double grayscale1 = (withmesh < 1)? withmesh: 1;
	const double grayscale = 1-grayscale1;
	st << grayscale << ' ' << grayscale << ' ' << grayscale << " RG\n";

	for(int n = 0; n < Th.nt; n++){

	    const int &v0 = Th(n,0);
	    const int &v1 = Th(n,1);
	    const int &v2 = Th(n,2);

	    st << r*ar*(Th(v0).x-x0) << ' ' << r*(Th(v0).y-y0) << " m ";
	    st << r*ar*(Th(v1).x-x0) << ' ' << r*(Th(v1).y-y0) << " l ";
	    st << r*ar*(Th(v2).x-x0) << ' ' << r*(Th(v2).y-y0) << " l ";
	    st << "s" << std::endl;
	}

	st << "Q\n";

    } // withmesh

    //------------------------------
    // edges
    //------------------------------
    st << "q\n";
    st << "1 w\n"; // setlinewidth
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";
    st << "0 0 0 RG\n";
    for(int k = 0; k < nEdges; k++){

	const int &v0 = Th( Th.bedges[k][0] );
	const int &v1 = Th( Th.bedges[k][1] );

	st << r*ar*(Th(v0).x-x0) << ' ' << r*(Th(v0).y-y0) << " m "
	   << r*ar*(Th(v1).x-x0) << ' ' << r*(Th(v1).y-y0) << " l S\n";
    }

    st << "Q\n";


    //------------------------------
    // legend
    //------------------------------
    if( legend ){

	st << "q\n";
	st << "1 w\n"; // setlinewidth
	st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";

	const double dy = r*(y1-y0)/nbfill;

	for(int m = 0; m < nbfill; m++){

	    const double f = (logscale)? fmin*pow(df,m): fmin + m*df;

	    setrgbcolor(st, f, palette, fmin, fmax, monochrome, logscale);
	    st << "rg\n";

	    st << sizex-PADDING   << " " << m*dy << " m "
               << sizex-PADDING/2 << " " << m*dy << " l "
               << sizex-PADDING/2 << " " << (m+1)*dy << " l "
	       << sizex-PADDING   << " " << (m+1)*dy << " l f\n";
	}

	const double EPS = 1e-10;
        const double dl = (logscale)?
	    exp( static_cast<double>(1)/(NUM_LABELS-1)*log(fmax/fmin) ):
	    (fmax-fmin)/(NUM_LABELS-1);
        for(int m = 0; m < NUM_LABELS; m++){

            const double f = (logscale)? fmin*pow(dl,m): fmin + m*dl;

	    if( nbfill < NUM_LABELS ){
		const double ff = (logscale)?
		    fmin * pow(df, static_cast<int>(log(f*0.999)/log(df)) ):
		    static_cast<int>((f*0.999-fmin)/df) * df + fmin;
		setrgbcolor(st, ff, palette, fmin, fmax, monochrome, logscale);
	    } else {
		setrgbcolor(st, f, palette, fmin, fmax, monochrome, logscale);
	    }
	    st << " rg\n";

            st << "BT /F1 " << textfontsize << " Tf "
	       << "1 0 0 1 " << sizex << " " << m*(r*(y1-y0)-textfontsize)/(NUM_LABELS-1) << " Tm "
	       << "(" << ((f >= 0)? "\\ ": "");

	    if( (fabs(f) > 1e-3) || (fabs(f) < 1e-12) ){
		st << std::setprecision(prec) << std::setfill('0') << f << ") Tj ET\n";
	    } else {
		st << std::resetiosflags(std::ios::fixed) << std::setiosflags(std::ios::scientific)
		   << std::setprecision(prec) << f << resetiosflags(std::ios::scientific)
		   << std::setiosflags(std::ios::fixed) << ") Tj ET\n";
	    }
	    // debug: st << std::setiosflags(std::ios::fixed)
        }

	st << "Q\n";

    } // legend

    //------------------------------
    // isoline (main routine in plot_P1_isoline)
    //------------------------------
    st << "q\n";
    st << linewidth << " w\n"; // setlinewidth
    st << "1 0 0 1 " << PADDING+marginl << " " << PADDING+marginb << " cm\n";

    if( isoline ){

	if( (fetype != "P1") && (fetype != "P1nc") ){
	    std::cout << "plotPDF() : isoline for vector filed is interpolated as P1 type" << std::endl;
	}

	std::vector<double> isoline_val;

	if( viso ){

	    for(int m = 0; m < viso->size(); m++)
		isoline_val.push_back( (*viso)[m] );

	} else if( logscale && (fmin > 0) ) {

	    // fmin * step^N = fmax <=> step^N = fmax/fmin
	    // <=> N = log_{step}(fmax/fmin) = (log(fmax/fmin))/log(step)
	    // <=> log(step) = (1/N)(log(fmax/fmin))
	    // <=> step = exp( (1/N)(log(fmax/fmin)) )
	    const double df = exp( (static_cast<double>(1)/NISOLINES)*(log(fmax/fmin)) );

	    isoline_val.push_back( fmin*sqrt(df) );
	    for(int m = 1; m < NISOLINES; m++)
		isoline_val.push_back( isoline_val[m-1] * df );

	} else {
	
	    if( logscale )
		std::cout << "plotPDF(): logscale for non-positive values.\n";

	    const double df = (fmax - fmin) / NISOLINES;
	    for(int m = 0; m < NISOLINES; m++)
		isoline_val.push_back( fmin + df/2 + m*df );
	}

	const KN<double> &f_P1 = f2;

	for(int k = 0; k < nTriangles; k++){

	    const int &v0 = Th(k,0);
	    const int &v1 = Th(k,1);
	    const int &v2 = Th(k,2);

	    const double vx[] = { Th(v0).x,     Th(v1).x,     Th(v2).x };
	    const double vy[] = { Th(v0).y,     Th(v1).y,     Th(v2).y };
	    const double vf[] = { f_P1[3*k+0],  f_P1[3*k+1],  f_P1[3*k+2] };

	    for(size_t m = 0; m < isoline_val.size(); m++){

		const double &value = isoline_val[m];

		std::vector<double> px, py;
		trackP1isoline( px, py, vx, vy, value, vf );

		assert( px.size() == py.size() );

		if( px.size() == 0 ) continue; // goto next isoline_val

		setrgbcolor(st, value, palette, fmin, fmax, monochrome, logscale);

		if( px.size() > 3 ){ 

		    // f(x,y)==value on all vertices, i.e. f(x,y) \equiv value on the triangle
		    st << "rg\n";
		    st << scale*ar*(vx[0] - x0) << ' ' << scale*(vy[0] - y0) << " m "
		       << scale*ar*(vx[1] - x0) << ' ' << scale*(vy[1] - y0) << " l "
		       << scale*ar*(vx[2] - x0) << ' ' << scale*(vy[2] - y0) << " l f\n";
		
		} else {

		    // assert( (px.size() == 2) || (px.size() == 3) );
		    st << "RG\n";
		    st << scale*ar*(px[0] - x0) << ' ' << scale*(py[0] - y0) << " m "
		       << scale*ar*(px[1] - x0) << ' ' << scale*(py[1] - y0) << " l\n";
		    st << "S\n";
		}
	    }
	}
    }

    //------------------------------
    // element(triangle)-wise process
    // This is final so that the legend does not hide arrows.
    //------------------------------
    for(int k = 0; k < nTriangles; k++){

	const int &v0 = Th(k,0);
	const int &v1 = Th(k,1);
	const int &v2 = Th(k,2);

	const double vx[] = { Th(v0).x, Th(v1).x, Th(v2).x };
	const double vy[] = { Th(v0).y, Th(v1).y, Th(v2).y };
	const double vfx[] = { fx[3*k+0],  fx[3*k+1],  fx[3*k+2] };
	const double vfy[] = { fy[3*k+0],  fy[3*k+1],  fy[3*k+2] };

	if( fromVertex ) {

	    // arrow from triangle vertices
	    const double vf2[] = { f2[3*k+0],  f2[3*k+1],  f2[3*k+2] };

	    for(int i = 0; i < 3; i++){

		if( (vf2[i] < fmin) || (fmax < vf2[i]) )
		    continue;

		plot_vector( st, vx[i],vy[i], vfx[i],vfy[i],vf2[i], alength_scale, ahead_scale,
			     x0,y0,scale,ar, fmin,fmax, unit_arrow,logscale, palette,monochrome );
	    }

	} else {

	    // arrow from center of triangles
	    const double cx = (vx[0] + vx[1] + vx[2])/3;
	    const double cy = (vy[0] + vy[1] + vy[2])/3;

	    const double cfx = (vfx[0]+vfx[1]+vfx[2])/3;
	    const double cfy = (vfy[0]+vfy[1]+vfy[2])/3;
	    const double cf2 = sqrt( cfx*cfx + cfy*cfy );

	    if( (cf2 < fmin) || (fmax < cf2) )
		continue;

	    plot_vector( st, cx,cy, cfx,cfy,cf2, alength_scale, ahead_scale,
			 x0,y0,scale,ar, fmin,fmax, unit_arrow,logscale, palette,monochrome );

	}

    } // element(triangle)-wise process

    st << "Q\n";

    return;
} //

//----------------------------------------------------------------------
// Interface
//----------------------------------------------------------------------

class PLOTPDF_Op : public E_F0mps
{
public:
    Expression eTh, ef, efilename, efx, efy, ez;
    static const int n_name_param = PLOTPDF_NOPTIONS;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    
    double arg(int i, Stack stack, double defvalue) const { return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): defvalue; }
    long arg(int i, Stack stack, long defvalue) const { return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): defvalue; }
    KN<double>* arg(int i, Stack stack, KN<double>* defvalue) const { return nargs[i] ? GetAny<KN<double>*>( (*nargs[i])(stack) ): defvalue; }
    bool arg(int i, Stack stack, bool defvalue) const { return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): defvalue; }
    
public:
    // mesh only
    PLOTPDF_Op(const basicAC_F0 &args, Expression filename, Expression th )
    : eTh(th), ef(0), efilename(filename), efx(0), efy(0), ez(0)
    {
        args.SetNameParam( n_name_param, name_param, nargs );
    }

    // real valued function
    PLOTPDF_Op(const basicAC_F0 &args, Expression filename, Expression th, Expression f )
    : eTh(th), ef(f), efilename(filename), efx(0), efy(0), ez(0)
    {
        args.SetNameParam( n_name_param, name_param, nargs );
    }

    // vector field
    PLOTPDF_Op(const basicAC_F0 &args, Expression filename, Expression th, Expression f, Expression fx, Expression fy)
    : eTh(th), ef(0), efilename(filename), efx(fx), efy(fy), ez(0)
    {
        args.SetNameParam( n_name_param, name_param, nargs );
    }

    // complex-valued function
    PLOTPDF_Op(const basicAC_F0 &args, Expression filename, Expression th, Expression f, const Complex z)
    : eTh(th), ef(0), efilename(filename), efx(0), efy(0), ez(f)
    {
        args.SetNameParam( n_name_param, name_param, nargs );
    }

    AnyType operator()(Stack stack) const;
};

basicAC_F0::name_and_type PLOTPDF_Op::name_param[] =
{
    {  "size",      &typeid(long)},
    {  "ar",        &typeid(double)},
    {  "fontscale", &typeid(double)},
    {  "viso",      &typeid(KN<double>*)},
    {  "nbiso",     &typeid(long)},
    {  "nbfill",    &typeid(long)},
    {  "frange",    &typeid(KN<double>*)},
    {  "gray",      &typeid(bool)},
    {  "bw",        &typeid(bool)},
    {  "value",     &typeid(bool)},
    {  "withmesh",  &typeid(double)},
    {  "meshpage",  &typeid(bool)},
    {  "index",     &typeid(bool)},
    {  "belabel",   &typeid(bool)},
    {  "isoline",   &typeid(bool)},
    {  "fill",      &typeid(bool)},
    {  "fetype",    &typeid(string*)},
    {  "title",     &typeid(string*)},
    {  "cmm",       &typeid(string*)},
    {  "fmargin",   &typeid(KN<double>*)},
    {  "prec",      &typeid(long)},
    {  "logscale",  &typeid(bool)},
    {  "zabs",      &typeid(bool)},
    {  "zreal",     &typeid(bool)},
    {  "zimag",     &typeid(bool)},
    {  "zarg",      &typeid(bool)},
    {  "coef",      &typeid(double)}, // arrow length scale
    {  "ArrowSize", &typeid(double)},
    {  "unitarrow", &typeid(bool)},
    {  "idcell",    &typeid(bool)},
    {  "idvert",    &typeid(bool)},
    {  "idedge",    &typeid(bool)},
    {  "lw",        &typeid(double)},
    {  "palette",   &typeid(KNM<double>*)}
    //----------------------------------------------------------------------
    // If you add new options, modify PLOTPDF_NOPTIONS in global namspace
    //----------------------------------------------------------------------
    //{  "nbarrow",   &typeid(long)}, // not implemented
    //{  "varrow",    &typeid(KN<double>*)}, // not implemented
};

AnyType PLOTPDF_Op::operator()(Stack stack) const
{
    // options
    const long draw_pane_size = arg(0,stack,DEFAULT_PAGESIZE);
    const double aspectratio = arg(1,stack,DEFAULT_ASPECTRATIO);
    const double fontscale = arg(2,stack,1.0);
    const KN<double> *const viso = arg(3,stack,reinterpret_cast<KN<double>*>(0));
    const long nbiso    = arg(4,stack,DEFAULT_ISOLINES);
    const long nbfill   = arg(5,stack,DEFAULT_FILL_COLORS);
    const KN<double> *const frange = arg(6,stack,reinterpret_cast<KN<double>*>(0));
    const bool gray     = arg(7,stack,DEFAULT_MONOCHROME);
    const bool bw       = arg(8,stack,DEFAULT_MONOCHROME);
    const bool legend   = arg(9,stack,DEFAULT_SHOW_LEGEND);
    const double withmesh = arg(10,stack,DEFAULT_WITHMESH);
    const bool mesh     = arg(11,stack,DEFAULT_SHOW_MESH);
    const bool index    = arg(12,stack,DEFAULT_SHOW_INDEX);
    const bool belabel  = arg(13,stack,DEFAULT_SHOW_BELABEL);
    const bool isoline  = arg(14,stack,DEFAULT_SHOW_ISOLINE);
    const bool fill     = arg(15,stack,DEFAULT_SHOW_FILL);
    const std::string fetype  = get_string(stack,nargs[16],DEFAULT_FETYPE);
    const std::string title   = get_string(stack,nargs[17],AppName);
    const std::string comment = get_string(stack,nargs[18],"");
    const KN<double> *const fmargin = arg(19,stack,reinterpret_cast<KN<double>*>(0));
    const long prec     = arg(20,stack,DEFAULT_PRECISION_LEGEND);
    const bool logscale = arg(21,stack,DEFAULT_LOGSCALE);
    const bool zabs = arg(22,stack,DEFAULT_ZABS);
    const bool zreal = arg(23,stack,DEFAULT_ZREAL);
    const bool zimag = arg(24,stack,DEFAULT_ZIMAG);
    const bool zarg = arg(25,stack,DEFAULT_ZARG);
    const double arrow_scale = arg(26,stack,DEFAULT_ARROW_SCALE);
    const double ahead_scale = arg(27,stack,DEFAULT_AHEAD_SCALE);
    const bool unit_arrow = arg(28,stack,DEFAULT_UNIT_ARROW);
    const bool idcell = arg(29,stack,DEFAULT_SHOW_IDCELL);
    const bool idvert = arg(30,stack,DEFAULT_SHOW_IDVERT);
    const bool idedge = arg(31,stack,DEFAULT_SHOW_IDEDGE);
    const double linewidth = arg(32,stack,DEFAULT_LINEWIDTH);
    //const KNM<double> *const RGBpalette = arg(33,stack,reinterpret_cast<KNM<double>*>(nullptr));
    const KNM<double> *const RGBpalette = nargs[33]?
	GetAny<KNM<double>*>((*nargs[33])(stack)): nullptr;
    //const long nbarrow = arg(34,stack,DEFAULT_ARROW_COLORS);
    //const KN<double> *const varrow = arg(35,stack,reinterpret_cast<KN<double>*>(0));
    const long nbarrow = DEFAULT_ARROW_COLORS;
    const KN<double> *const varrow = reinterpret_cast<KN<double>*>(0);

    if( verbosity >= 90 ){

	// invoked by FreeFem++-nw -v 90 users.edp

	std::cout << "plotPDF:options -------------------------" << std::endl;
	std::cout << "plotPDF:option: pagesize=" << draw_pane_size << std::endl;
	std::cout << "plotPDF:option: aspectratio=" << aspectratio << std::endl;
	std::cout << "plotPDF:option: fontscale=" << fontscale << std::endl;
	if( viso )
	    std::cout << "plotPDF:option: viso->size()=" << viso->size() << std::endl;
	else
	    std::cout << "plotPDF:option: viso[] is empty" << std::endl;
	std::cout << "plotPDF:option: nbiso=" << nbiso << std::endl;
	std::cout << "plotPDF:option: nbfill=" << nbfill << std::endl;
	if( frange )
	    std::cout << "plotPDF:option: frange->size()=" << frange->size() << std::endl;
	else
	    std::cout << "plotPDF:option: frange[] is empty" << std::endl;
	std::cout << "plotPDF:option: gray=" << gray << std::endl;
	std::cout << "plotPDF:option: bw=" << bw << std::endl;
	std::cout << "plotPDF:option: value=" << legend << std::endl;
	std::cout << "plotPDF:option: withmesh=" << withmesh << std::endl;
	std::cout << "plotPDF:option: mesh=" << mesh << std::endl;
	std::cout << "plotPDF:option: belabel=" << belabel << std::endl;
	std::cout << "plotPDF:option: index=" << index << std::endl;
	std::cout << "plotPDF:option: isoline=" << isoline << std::endl;
	std::cout << "plotPDF:option: fill=" << fill << std::endl;
	std::cout << "plotPDF:option: fetype=" << fetype << std::endl;
	std::cout << "plotPDF:option: title=" << title << std::endl;
	std::cout << "plotPDF:option: cmm=" << comment << std::endl;
	if( fmargin ){
	    std::cout << "plotPDF:option: fmargin->size()=" << fmargin->size() << ", [";
	    for( int i = 0; i < fmargin->size(); i++)
		std::cout << (*fmargin)[i] << ',';
	    std::cout << "]" << std::endl;
	} else {
	    std::cout << "plotPDF:option: fmargin[] is empty" << std::endl;
	}
	std::cout << "plotPDF:option: prec=" << prec << std::endl;
	std::cout << "plotPDF:option: logscale=" << logscale << std::endl;
	std::cout << "plotPDF:option: zabs=" << zabs << std::endl;
	std::cout << "plotPDF:option: zreal=" << zreal << std::endl;
	std::cout << "plotPDF:option: zimag=" << zimag << std::endl;
	std::cout << "plotPDF:option: zarg=" << zarg << std::endl;
	std::cout << "plotPDF:option: coef=" << arrow_scale << std::endl;
	std::cout << "plotPDF:option: ArrowSize=" << ahead_scale << std::endl;
	std::cout << "plotPDF:option: unitarrow=" << unit_arrow << std::endl;
#if 0
	// not implemented yet
	std::cout << "plotPDF:option: nbarrow=" << nbarrow << std::endl;
	if( varrow ){
	    std::cout << "plotPDF:option: varrow->size()=" << varrow->size() << ", [";
	    for( int i = 0; i < varrow->size(); i++)
		std::cout << (*varrow)[i] << ',';
	    std::cout << "]" << std::endl;
	} else {
	    std::cout << "plotPDF:option: varrow[] is empty" << std::endl;
	}
#endif
	std::cout << "plotPDF:option: idcell=" << idcell << std::endl;
	std::cout << "plotPDF:option: idvert=" << idvert << std::endl;
	std::cout << "plotPDF:option: idedge=" << idedge << std::endl;
	std::cout << "plotPDF:option: lw=" << linewidth << std::endl;

	if( RGBpalette ){
	    std::cout << "plotPDF:option: palette["
		      << RGBpalette->N() << ',' << RGBpalette->M() << "]=";
	    // RGBpalette->size()
	    std::cout << "[ ";
	    for(int i = 0; i < RGBpalette->N(); i++){
		std::cout << '[';
		for(int j = 0; j < RGBpalette->M(); j++){
		    std::cout << RGBpalette->operator()(i,j);
		    if( j != RGBpalette->M()-1 )
			std::cout << ',';
		    else
			std::cout << ']';
		}
		if( i != RGBpalette->N()-1 )
		    std::cout << ", ";
	    }
	    std::cout << " ]" << std::endl;
	} else {
	    std::cout << "plotPDF:option: palette[] is empty" << std::endl;
	}

    }

    //----------------------------------------
    // detecting filename ending with ".pdf"
    //----------------------------------------
    const std::string *const filename = GetAny<std::string*>((*efilename)(stack));
    ffassert(filename);

    const char PDFextension[] = ".pdf";
    const size_t pos = filename->rfind( PDFextension );
    const bool hasPDFextension = ( (pos != std::string::npos) && (pos == (filename->length() - strlen(PDFextension))) );

    const std::string filename_with_extension = (hasPDFextension)? *filename: *filename + PDFextension;

    if(! (mesh || index || belabel || (isoline && (nbiso > 0)) || (fill && (nbfill > 0)) || ez ) ){
	std::cerr << "plotPDF() : No file output : " << filename_with_extension << std::endl;
	return false;
    }

    if(! (mesh || index || belabel || ef || (!ef && efx && efy) || ez ) ){
	std::cerr << "plotPDF() : No file output" << filename_with_extension << std::endl;
	return false;
    }

    const bool monochrome = bw || gray;

    //const std::string PDFtitle(*filename);
    const std::string PDFtitle(title.c_str());
    SimplePDFModule mesh_figure( filename_with_extension.c_str(), PDFtitle.c_str() );

    //----------------------------------------
    // Lower & Upper Bound of Domain
    //----------------------------------------
    const Mesh *const pTh = GetAny<const Mesh *const>((*eTh)(stack));
    ffassert(pTh);
    const Fem2D::Mesh & Th(*pTh);
    const int nVertices  = Th.nv;
    const int nTriangles = Th.nt;

    R2 Pmin, Pmax;
    Th.BoundingBox(Pmin, Pmax);

    const double &x0 = Pmin.x;
    const double &y0 = Pmin.y;
    
    const double &x1 = Pmax.x;
    const double &y1 = Pmax.y;

    //--------------------------------------------------
    // comment
    //--------------------------------------------------
    const int comment_height
	= (comment.length() == 0)? 0: (DEFAULT_COMMENT_FONTSIZE + DEFAULT_COMMENT_BASE);

    //----------------------------------------
    // user specific margin (by option)
    //----------------------------------------
    const int margin[] = {
	( (fmargin) && (fmargin->size() >= 1) )? static_cast<int>((*fmargin)[0]): DEFAULT_MARGIN[0], // left
	( (fmargin) && (fmargin->size() >= 2) )? static_cast<int>((*fmargin)[1]): DEFAULT_MARGIN[1], // bottom
	( (fmargin) && (fmargin->size() >= 3) )? static_cast<int>((*fmargin)[2]): DEFAULT_MARGIN[2], // right
	( (fmargin) && (fmargin->size() >= 4) )? static_cast<int>((*fmargin)[3]): DEFAULT_MARGIN[3]  // top
    };
    const int &marginl = margin[0];
    const int &marginb = margin[1];
    const int &marginr = margin[2];
    const int &margint = margin[3];

    const int legend_width = 3*PADDING + static_cast<int>(LEGEND_FONTWIDTH*prec);

    //--------------------------------------------------
    // scaling factor, [x0,x1]*[y0,y1] -> [0,PAGESIZE]^2
    //--------------------------------------------------
    if( aspectratio < 0 ){
	std::cout << "plotPDF(): ar should be positive." << std::endl;
    }
    const double ar = (aspectratio <= 0)? 1: (y1-y0)/(x1-x0)*aspectratio;
    const double rx = draw_pane_size / ((x1 - x0) * ar);
    const double ry = draw_pane_size / (y1 - y0);
    const double scale = (rx < ry)? rx: ry;   // min(rx,ry)

    const int sizex = static_cast<int>(2*PADDING + scale*(x1-x0)*ar);
    const int sizey = static_cast<int>(2*PADDING + scale*(y1-y0)) + comment_height;

    const double index_fontsize = DEFAULT_INDEX_FONTSIZE * fontscale;

    //------------------------------
    // color palette
    //------------------------------
    // compile error
    // KNM<double> default_palette_array
    //  = KNM_<double>(DEFAULT_PALETTE, DEFAULT_PALETTE_NCOLORS, 3);

    const KNM<double> default_palette_array( DEFAULT_PALETTE_NCOLORS, 3 );

    // why can we modify const KNM<double> ?
    for(int i = 0; i < DEFAULT_PALETTE_NCOLORS; i++)
	for(int j = 0; j < 3; j++)
	    default_palette_array(i,j) = DEFAULT_PALETTE[i][j];

    bool validPalette = RGBpalette;
    if( RGBpalette && ((RGBpalette->N() == 0) || (RGBpalette->M() < 3)) ){
	std::cout << "plotPDF(): palette is given with illeagal form." << std::endl;
	std::cout << "plotPDF(): using default color palette." << std::endl;
	validPalette = false;
    }
    const KNM<double> palette = (RGBpalette && validPalette)? (*RGBpalette): default_palette_array;

    for(int i = 0; i < palette.N(); i++)
	for(int j = 0; j < 3; j++)
	    palette(i,j) /= 255; // why can we modify const KNM<double> ?

    //--------------------------------------------------
    // PDF Contents
    //--------------------------------------------------
    std::stringstream Content;
    Content.setf( std::ios::fixed );  // PDF does not support floating-point format
    Content.precision(3);

    if( mesh ){
	mesh_figure.addBookmark( "Mesh" );
	plot_mesh( Content, Th, scale, ar, x0, y0, marginl, marginb, index_fontsize, monochrome, withmesh, linewidth, idcell, idvert, idedge );
	if( comment.length() > 0 )
	    addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	mesh_figure.addPage( Content, sizex, sizey, margin );
    }

    if( index ){
	mesh_figure.addBookmark( "Mesh with Index" );
	plot_mesh( Content, Th, scale, ar, x0, y0, marginl, marginb, index_fontsize, monochrome, withmesh, linewidth, idcell, idvert, idedge, 1 );
	if( comment.length() > 0 )
	    addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	mesh_figure.addPage( Content, sizex, sizey, margin );
    }

    if( belabel ){
	mesh_figure.addBookmark( "Mesh with Boundary Label" );
	plot_mesh( Content, Th, scale, ar, x0, y0, marginl, marginb, index_fontsize, monochrome, withmesh, linewidth, idcell, idvert, idedge, 2 );
	if( comment.length() > 0 )
	    addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	mesh_figure.addPage( Content, sizex, sizey, margin );
    }

    if( !ef && !efx && !efy && !ez )
	return true;

    //----------------------------------------------------------------------
    // vector-valued function
    //----------------------------------------------------------------------
    const double legend_fontsize = DEFAULT_LEGEND_FONTSIZE;

    if( efx && efy ){

	if( verbosity >= 5 )
	    std::cout << "plotPDF() : Vector filed is interpolated as P1 type." << std::endl;

	const bool fromVertex = (fetype == "P1") || (fetype == "P2");

	const int &nTriangles = Th.nt;
	KN<double> fx( 3*nTriangles ), fy( 3*nTriangles ), f2( 3*nTriangles );

	// find function values on vertices
	for(int it = 0; it < nTriangles; it++){
	
	    for(int iv = 0; iv < 3; iv++){
		MeshPointStack(stack)->setP(pTh,it,iv); // at the iv-th vertex
		fx[3*it+iv] = GetAny<double>( (*efx)(stack) ); // Expression ef is atype<double>()
		fy[3*it+iv] = GetAny<double>( (*efy)(stack) ); // Expression ef is atype<double>()
		f2[3*it+iv] = sqrt( fx[3*it+iv]*fx[3*it+iv] + fy[3*it+iv]*fy[3*it+iv] ); // Euclid norm
	    }
	}

        mesh_figure.addBookmark( "Profile : Vector-valued Function" );
	plot_vector2flow( Content, Th, fx, fy, f2, fromVertex, palette, arrow_scale, unit_arrow, 
			  linewidth, ahead_scale, nbarrow, varrow,
			  sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			  legend_fontsize, monochrome, legend, prec, logscale,
			  withmesh, nbfill, frange, fetype, isoline, nbiso, viso );
	if( comment.length() > 0 )
	    addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	return true;
    }

    //----------------------------------------------------------------------
    // complex-valued function (Interpolate Input Array as P1 function)
    //----------------------------------------------------------------------
    if( ez ){

	if( verbosity >= 5 )
	    std::cout << "plotPDF() : Complex-valued function is interpolated as P1 type." << std::endl;

	if( (fetype != "P1") && (fetype != "P1nc") ){
	    std::cout << "plotPDF() : Unknown fetype : " << fetype << std::endl;
	    std::cout << "plotPDF() : Interpolated as piecewise-linear" << std::endl;
	}

	const int &nTriangles = Th.nt;
	KN< Complex > fz( 3*nTriangles );
	KN<double> fr( 3*nTriangles ), fi( 3*nTriangles ), fm( 3*nTriangles ), fa( 3*nTriangles );

	// find function values on vertices
	for(int it = 0; it < nTriangles; it++){
	
	    for(int iv = 0; iv < 3; iv++){
		MeshPointStack(stack)->setP(pTh,it,iv); // at the iv-th vertex
		fz[3*it+iv] = GetAny<Complex>( (*ez)(stack) ); // Expression ez is atype<Complex>()
		fr[3*it+iv] = fz[3*it+iv].real();
		fi[3*it+iv] = fz[3*it+iv].imag();
		fm[3*it+iv] = abs( fz[3*it+iv] );
		fa[3*it+iv] = atan2( fi[3*it+iv], fr[3*it+iv] );
	    }
	}

	// P1, P1nc and Other Interpolated Functions
	//------------------------------
	// isoline
	//------------------------------
	if( isoline && (nbiso > 0) && zreal ){
	    mesh_figure.addBookmark( "Isoline : Real Part" );
	    plot_P1_isoline( Content, Th, fr, // real part
			     palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			     legend_fontsize, monochrome, legend, prec, logscale,
			     withmesh, nbiso, viso, linewidth );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}

	if( isoline && (nbiso > 0) && zimag ){
	    mesh_figure.addBookmark( "Isoline : Imaginary Part" );
	    plot_P1_isoline( Content, Th, fi, // imaginary part
			     palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			     legend_fontsize, monochrome, legend, prec, logscale,
			     withmesh, nbiso, viso, linewidth );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}
	
	if( isoline && (nbiso > 0) && zabs ){
	    mesh_figure.addBookmark( "Isoline : Modulus of Complex-valued Function" );
	    plot_P1_isoline( Content, Th, fm, // modulus
			     palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			     legend_fontsize, monochrome, legend, prec, logscale,
			     withmesh, nbiso, viso, linewidth );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}

	if( isoline && (nbiso > 0) && zarg ){
	    mesh_figure.addBookmark( "Isoline : Argument" );
	    plot_P1_isoline( Content, Th, fa, // argument
			     palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			     legend_fontsize, monochrome, legend, prec, logscale,
			     withmesh, nbiso, viso, linewidth );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}
	
	//------------------------------
	// fill-style
	//------------------------------
	if( fill && (nbfill > 0) && zreal ){
	    mesh_figure.addBookmark( "Real Part of Complex-valued Function" );
	    plot_P1_fill( Content, Th, fr, // real part
			  palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			  legend_fontsize, monochrome, legend, prec, logscale,
			  withmesh, nbfill, frange );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}

	if( fill && (nbfill > 0) && zimag ){
	    mesh_figure.addBookmark( "Imaginary Part of Complex-valued Function" );
	    plot_P1_fill( Content, Th, fi, // imaginary part
			  palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			  legend_fontsize, monochrome, legend, prec, logscale,
			  withmesh, nbfill, frange );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}

	if( fill && (nbfill > 0) && zabs ){
	    mesh_figure.addBookmark( "Modulus of Complex-valued Function" );
	    plot_P1_fill( Content, Th, fm, // modulus
			  palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			  legend_fontsize, monochrome, legend, prec, logscale,
			  withmesh, nbfill, frange );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}

	if( fill && (nbfill > 0) && zarg ){
	    mesh_figure.addBookmark( "Argument of Complex-valued Function" );
	    plot_P1_fill( Content, Th, fa, // argument
			  palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			  legend_fontsize, monochrome, legend, prec, logscale,
			  withmesh, nbfill, frange );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}

	return true;
    }

    //----------------------------------------------------------------------
    // Output Specific FE-type (scalar function)
    //----------------------------------------------------------------------
    if( fetype == "P0" ){ 

	const int &nTriangles = Th.nt;
	KN<double> f_FE( nTriangles );

	for(int it = 0; it < nTriangles; it++){

	    const int &v0 = Th(it,0);
	    const int &v1 = Th(it,1);
	    const int &v2 = Th(it,2);

	    const double x = (Th(v0).x+Th(v1).x+Th(v2).x)/3;
	    const double y = (Th(v0).y+Th(v1).y+Th(v2).y)/3;
	    MeshPointStack(stack)->set(x,y); // barycenter

	    f_FE[it] = GetAny<double>( (*ef)(stack) ); // Expression ef is atype<double>()
	}

        mesh_figure.addBookmark( "P0 Profile" );
	plot_P0_fill( Content, Th, f_FE,
		      palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
		      legend_fontsize, monochrome, legend, prec, logscale,
		      withmesh, nbfill, frange );
	if( comment.length() > 0 )
	    addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	return true;
    } 

    if( fetype == "P2" ){

	KN<double> f_FE( 6*nTriangles );

	// loop over triangle
	for(int it = 0; it < nTriangles; it++){

	    const int &v0 = Th(it,0);
	    const int &v1 = Th(it,1);
	    const int &v2 = Th(it,2);

	    const double x[] = { Th(v0).x, Th(v1).x, Th(v2).x };
	    const double y[] = { Th(v0).y, Th(v1).y, Th(v2).y };

	    // f_FE[i] = func(v[i]), f_FE[i+3] = func(e[i]), i=0,1,2.
	    for(int iv = 0; iv < 3; iv++){
		
		MeshPointStack(stack)->setP(pTh,it,iv); // at the iv-th vertex
		f_FE[ 6*it + iv ] = GetAny<double>( (*ef)(stack) ); // Expression ef is atype<double>()

		const double ex = (x[(iv+1)%3]+x[(iv+2)%3]) / 2;  // mid-point of iv-th edge
		const double ey = (y[(iv+1)%3]+y[(iv+2)%3]) / 2;

		MeshPointStack(stack)->set( ex, ey );
		f_FE[ 6*it + 3 + iv ] = GetAny<double>( (*ef)(stack) ); // Expression ef is atype<double>()
	    }
	}

	if( isoline && (nbiso > 0) ){
	    mesh_figure.addBookmark( "P2 Isoline" );
	    plot_P2_isoline( Content, Th, f_FE,
			     palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			     legend_fontsize, monochrome, legend, prec, logscale,
			     withmesh, nbiso, viso, linewidth );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}

	if( fill && (nbfill > 0) ){
	    mesh_figure.addBookmark( "P2 Profile" );
	    plot_P2_fill( Content, Th, f_FE,
			  palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			  legend_fontsize, monochrome, legend, prec, logscale,
			  withmesh, nbfill, frange );
	    if( comment.length() > 0 )
		addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	    mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
	}

	return true;
    }

    //----------------------------------------------------------------------
    // Interpolate Input Array as P1 (or P1nc) function
    //----------------------------------------------------------------------
    KN<double> f_FE( 3*nTriangles );
    
    if( (fetype != "P1") && (fetype != "P1nc") ){
	std::cout << "plotPDF() : Unknown fetype : " << fetype << std::endl;
	std::cout << "plotPDF() : Interpolated as piecewise-linear" << std::endl;
    }

    // find function values on vertices
    for(int it = 0; it < nTriangles; it++){
	
	for(int iv = 0; iv < 3; iv++){
	    MeshPointStack(stack)->setP(pTh,it,iv); // at the iv-th vertex
	    f_FE[3*it+iv] = GetAny<double>( (*ef)(stack) ); // Expression ef is atype<double>()
	}
    }

#if 0
    {   // output for GNUPLOT
	// usage in GNUPLOT:
	// gnuplot> splot 'FreeFEM2Gnuplot.dat'
	std::ofstream fout("FreeFEM2Gnuplot.dat");
	for(int it = 0; it < nTriangles; it++){

	    const int &v0 = Th(it,0);
	    const int &v1 = Th(it,1);
	    const int &v2 = Th(it,2);

	    const double vx[] = { Th(v0).x, Th(v1).x, Th(v2).x };
	    const double vy[] = { Th(v0).y, Th(v1).y, Th(v2).y };

	    fout << vx[0] << ' ' << vy[0] << ' ' << f_FE[ 3*it+0 ] << std::endl;
	    fout << vx[1] << ' ' << vy[1] << ' ' << f_FE[ 3*it+1 ] << std::endl << std::endl;
	    fout << vx[2] << ' ' << vy[2] << ' ' << f_FE[ 3*it+2 ] << std::endl;
	    fout << vx[2] << ' ' << vy[2] << ' ' << f_FE[ 3*it+2 ] << std::endl;
	    fout << std::endl << std::endl;
	}
	fout.close();
    }
#endif

    // P1, P1nc and Other Interpolated Functions
    if( isoline && (nbiso > 0) ){
	mesh_figure.addBookmark( "Isoline" );
	plot_P1_isoline( Content, Th, f_FE,
			 palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
			 legend_fontsize, monochrome, legend, prec, logscale,
			 withmesh, nbiso, viso, linewidth );
	if( comment.length() > 0 )
	    addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
    }
	
    if( fill && (nbfill > 0) ){
	mesh_figure.addBookmark( "Profile" );
	plot_P1_fill( Content, Th, f_FE,
		      palette, sizex, sizey, scale, ar, x0, y0, y1, marginl, marginb,
		      legend_fontsize, monochrome, legend, prec, logscale,
		      withmesh, nbfill, frange );
	if( comment.length() > 0 )
	    addComment( Content, scale*(y1-y0), marginl, marginb, fontscale, comment );
	mesh_figure.addPage( Content, (legend? (sizex+legend_width): sizex), sizey, margin );
    }
    return true;
}

class PLOTPDF: public OneOperator
{
    const int argc;
public:
    PLOTPDF()     : OneOperator( atype<long>(), atype<std::string*>(), atype<const Mesh*>() ), argc(2) {}
    PLOTPDF(int)  : OneOperator( atype<long>(), atype<std::string*>(), atype<const Mesh*>(), atype<double>() ), argc(3) {}
    //PLOTPDF(int,int) : OneOperator( atype<long>(), atype<std::string*>(), atype<const Mesh*>(), atype<E_Array>() ), argc(4) {}
    // vector field
    PLOTPDF(int,int) : OneOperator( atype<long>(), atype<std::string*>(), atype<const Mesh*>(), atype<E_Array>() ), argc(3) {}

    // complex-valued function
    PLOTPDF(int,int,int) : OneOperator( atype<long>(), atype<std::string*>(), atype<const Mesh*>(), atype< Complex >() ), argc(5) {}
    
    E_F0 * code(const basicAC_F0 & args) const
    {
        if(argc == 2) {

	    // mesh
            return  new PLOTPDF_Op( args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]) );

        } else if(argc == 3) {
	    
	    if( BCastTo<double>( args[2] ) ){

		// scalar (real-valued) function
		return  new PLOTPDF_Op( args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]) );

	    } else if( args[2].left() == atype<E_Array>() ) {

		// vector-valued function
		const int DimVector = 2;

		const E_Array *const a2 = dynamic_cast<const E_Array*>( args[2].LeftValue() );
		if( a2->size() != DimVector ){
		    std::cerr << "plotPDF() : Error: The size of vector-valued function is not valid." << std::endl;
		    ffassert(0);
		}
		Expression fx = to<double>( (*a2)[0] );
		Expression fy = to<double>( (*a2)[1] );
		return  new PLOTPDF_Op( args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), fx, fy ); 
	    }

	} else if(argc == 5) {

	    // complex-valued function
	    return new PLOTPDF_Op( args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), static_cast<Complex>(0) );
	}
	ffassert(0);
    }
};

//----------------------------------------------------------------------

static void Load_Init()
{
    if ( verbosity && (mpirank == 0) ) { std::cout << " load: plotPDF " << std::endl; }

    Global.Add("plotPDF", "(", new PLOTPDF);        // mesh only
    Global.Add("plotPDF", "(", new PLOTPDF(0));     // real valued
    Global.Add("plotPDF", "(", new PLOTPDF(0,0));   // vector valued
    Global.Add("plotPDF", "(", new PLOTPDF(0,0,0)); // complex-valued
}

LOADFUNC( Load_Init );

//----------------------------------------------------------------------
// End of file
//----------------------------------------------------------------------
