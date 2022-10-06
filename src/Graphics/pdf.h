// -*- C++ -*-
// Time-stamp: "2019-12-29 10:38:37 fujiwara"
//
// SUMMARY : Class in order to generate Portable Document Format
// ORG     : Graduate School of Informatics, Kyoto University, Japan
// AUTHOR  : Hiroshi Fujiwara
// E-MAIL  : fujiwara@acs.i.kyoto-u.ac.jp
// 
// The newest version is avalilable at:
// http://www-an.acs.i.kyoto-u.ac.jp/~fujiwara/ff++-programs
//
//----------------------------------------------------------------------
#ifndef PDF_H
#define PDF_H

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <list>
#include <ctime>
#include <cstring> // strlen
#include <cstdlib> // exit

//#define HAVE_ZLIB   // need zlib : CC ***.cpp -lz
//#undef HAVE_ZLIB

//----------------------------------------------------------------------
// Simple PDF class
//----------------------------------------------------------------------
#include <iomanip>
#include <sstream>
#include <list>
#include <ctime>
#include <cstring> // strlen
#include <cstdlib> // exit

class SimplePDF_FF
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

    SimplePDF_FF( const char *const filename, const char *const title = "",
	       const char *const AppName = "" );

    ~SimplePDF_FF();

    void addPage( const std::stringstream &ContentStream, const int WIDTH, const int HEIGHT );

    void addBookmark( const char *const BookmarkLabel );
};

SimplePDF_FF::SimplePDF_FF( const char *const PDFfilename, const char *const title,
                      const char *const AppName ) : filename( PDFfilename ), DocumentTitle( title ), page(0)
{
    std::ofstream fout( filename.c_str(), std::ios::binary );

    if( !fout ){
	std::cerr << "Cannot open the file: " << filename << std::endl;
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

    fout.open( filename.c_str(), std::ios::app );

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

SimplePDF_FF::~SimplePDF_FF()
{
    std::ofstream fout( filename.c_str(), std::ios::app );

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

int SimplePDF_FF::deflate_compress( char* &outbuf, const std::string &Stream ) const
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

void SimplePDF_FF::addPage( const std::stringstream &ContentStream,
		         const int WIDTH, const int HEIGHT )
{
    //----------------------------------------
    // PageObject
    //----------------------------------------
    std::stringstream strPageObject;

    strPageObject << page_obj_offset+2*page << " 0 obj\n"
		  << "<<\n"
		  << "  /Type /Page\n"
		  << "  /Parent 3 0 R\n"
		  << "  /Resources << /Font << /F1 7 0 R >> >>\n"
		  << "  /MediaBox [0 0 " << WIDTH << ' ' << HEIGHT << "]\n"
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
    std::ofstream fout( filename.c_str(), std::ios::app );

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

void SimplePDF_FF::addBookmark( const char *const BookmarkLabel )
{
    OutlineItem item;

    item.headPageObjectNumber = page_obj_offset+2*page;

    item.label = new char [ strlen(BookmarkLabel)+1 ];
    strcpy( item.label, BookmarkLabel );
    
    outline.push_back( item );

    return;
}

#endif // PDF_H

//----------------------------------------------------------------------
// End of file
//----------------------------------------------------------------------
