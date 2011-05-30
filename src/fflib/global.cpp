// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*
 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "config-wrapper.h"

#include <complex>
#include "AFunction.hpp"
#include "error.hpp"
#include "lex.hpp"
#include "RNM.hpp"
#include <queue>
#include "environment.hpp"

#define  FF_GRAPH_PTR_DCL
#include "rgraph.hpp"
//#include "PlotStream.hpp"


 long verbosity = 1;
 long searchMethod = 0; //pichon
 FILE *ThePlotStream=0; //  Add for new plot. FH oct 2008

  KN<String> *pkarg;//  for the list of argument  mars 2010 
 Map_type_of_map map_type_of_map ; //  to store te type 
Map_type_of_map map_pair_of_type ; //  to store te type 

 basicForEachType *  typevarreal,  * typevarcomplex;  //  type of real and complex variable


 mylex *zzzfff;
  bool lexdebug;
#include "lg.tab.hpp"
 YYSTYPE *plglval;

 int TheCurrentLine=-1; // unset: by default
 long mpisize=0,mpirank=0;

 
   C_F0 *pOne=0,*pZero=0,*pminusOne=0;
// const C_F0 & One(*pOne), &Zero(*pZero);
 
 Polymorphic * TheOperators=0, //=new Polymorphic(), 
             * TheRightOperators=0; //=new Polymorphic();

 TableOfIdentifier Global;

 long E_Border::Count =0;

 typedef list<TableOfIdentifier *> ListOfTOfId;

  ListOfTOfId tables_of_identifier;

  const int AC_F0::MaxSize=1024; // maximal number of parameters



map<const string,basicForEachType *> map_type;
bool showCPU= false;


size_t CodeAlloc::nb=0, CodeAlloc::lg=0,CodeAlloc::nbpx=0,CodeAlloc::chunk=2048; 
size_t CodeAlloc::nbt,CodeAlloc::nbdl=0;
CodeAlloc ** CodeAlloc::mem=0;
bool CodeAlloc::sort=true;
bool  CodeAlloc::cleanning=false;
bool echo_edp=true; // add F.H of remove script dump 

//  add F. Hecht 
EnvironmentData  ffenvironment;


