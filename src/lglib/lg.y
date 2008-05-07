%{ 
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
#define eflval yylval 
#include <iostream>
#include  <complex>
#include <string>
#include "error.hpp"
class Iden;
#include "strversionnumber.hpp"

#ifdef __MWERKS__
#ifdef __INTEL__
#include <malloc.h>
#else
#include <alloca.h>
#endif
#endif
#include "AFunction.hpp"
//  to reserve space to graphical pointer function
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpacen.hpp" 
#include "FESpace.hpp" 
#include "lgfem.hpp" 
#include "lex.hpp"
#include "environment.hpp"

class Routine;
bool load(string s);

template <class R> class FE;
template <class R,int i> class FE_;
extern mylex *zzzfff;
#ifdef PARALLELE
  void initparallele(int &, char **&);
  void init_lgparallele();
  void end_parallele();
#endif
#ifdef HAVE_LIBARPACK
  void init_eigenvalue();
#endif
   
  aType dcltype;
const int nbembtype=10;
aType rettype[nbembtype];
Block * routineinblock[nbembtype]; // Add FH july 2005 pb clean on return 
int kkembtype=-1;
int inloopcount=0;
Block *currentblock;
// Add FH july 2005 
//  problem clean variable after break,continue and return.
const int sizeStackOfLoop=100; 
Block * StackOfLoop[sizeStackOfLoop];
// end ADD
double CPUcompileInit =0;
//class pfes;
C_F0  fespacetype;
bool fespacecomplex;

int ShowAlloc(const char *s,size_t &);
inline int yylex()  {return zzzfff->scan();}
inline int lineno() {return zzzfff->lineno();}

extern bool withrgraphique;
inline void fingraphique()
 { if(withrgraphique) 
   { withrgraphique=false;
    rattente(1);
    closegraphique();
  }}

void lgerror (const char* s) ;

%}

%union{ 
 double dnum;
 long lnum;
 char * str;
 char oper[8];
 CC_F0 cexp;
 Routine   *routine;
 AC_F0 args;
 aType type;
 CListOfInst cinst;
 Block * block; 
 ListOfId *clist_id;
/* ListCatch * clist_Catchs;*/
}

/* BISON Declarations */

%type <cinst>   input
%type <cinst>   instructions
%type <cexp>   instruction
%type  <cexp>  try
%type  <cexp>  catchs
/* %type  <cexp>  throw */
%type <cexp>  declaration 
%type <cexp>  declaration_for
%type <cexp>  list_of_dcls
%type <cexp>  fespace_def
%type <cexp>  fespace_def_list
%type <cexp>   Expr
%type <cexp>   no_comma_expr
%type <cexp>   sub_script_expr
%type <cexp>   no_ternary_expr
%type <cexp>   no_set_expr
%type <cexp>   unary_expr
%type <cexp>   pow_expr
%type <cexp>   primary
%type <oper>   unop
%type <args>  parameters
%type <args>  array 
%type <args>  parameters_list
%type  <cexp> begin
%type  <cexp> end
%type  <clist_id> list_of_id_args
%type  <clist_id> list_of_id1
%type <cexp>  spaceIDs
%type <cexp>  spaceIDa
%type <cexp>  spaceIDb
%type <cexp> ID_space
%type <cexp>  ID_array_space
%type <args> bornes;
%type <args>   border_expr;
%type <type>   type_of_dcl;
%type <str> id;
/* Add precedence rules to solve dangling else s/r conflict */

%nonassoc IF
%nonassoc ELSE

%left <oper> ','
%right <oper> '=' SET
%left  <oper> LTLT GTGT
%left  <oper> OR '|' 
%left  <oper> AND '&'  
%left  <oper> EQ NE
%left  <oper>  '<' '>' LE GE 
%left  <oper>  '+' '-'
%left  <oper> '*' '/'  '%' DOTSTAR DOTSLASH
%right <oper> UNARY PLUSPLUS MOINSMOINS '!'
%right <oper>  '^' '\''
%right <oper>  '_' 
%left  <oper>  '(' '[' '.'

%token <oper> ')'  ']' 
 
%token <lnum> LNUM
%token <dnum> DNUM
%token <dnum> CNUM
%token <str> ID
%token <str> FESPACEID
%token <str> IDPARAM
%token <str> STRING

%token ENDOFFILE
%token INCLUDE
%token LOAD
%token BIDON 

%token FOR
%token WHILE
%token IF
%token ELSE
%token BREAK
%token CONTINUE
%token RETURN
%token TRY
%token CATCH
%token THROW

%token <type> TYPE
%token <type> FUNCTION
%token <str> FESPACE

%token DOTSTAR
%token DOTSLASH
%token AND
%token OR
%token EQ
%token NE
%token LE
%token GE
%token PLUSPLUS
%token MOINSMOINS
%token SET
%token LTLT
%token PLUSEQ
%token MOINSEQ
%token MULEQ
%token DIVEQ
%token GTGT
%token ARROW
%token BORDER
%token CURVE
%token SOLVE

%% 

start:   input ENDOFFILE {
	            
                        size_t sizestack = currentblock->size()+1024 ; //  before close 
                        $1+=currentblock->close(currentblock);
                        if(verbosity) cout << " sizestack + 1024 =" << sizestack << "  ( " << sizestack-1024 <<" )\n" ;   
                        size_t lg0,lg1;                       
                        int NbPtr = ShowAlloc("init execution ",lg0); // number of un delele ptr
                        if(verbosity) cout << endl;  
                        { Stack stack = newStack(sizestack);
                        double CPUcompile= CPUtime();
                        try {                  
                          $1.eval(stack);}
                        catch ( E_exception & e)  {
                          cerr << e.what() << endl;
                          return 1; }
                        catch( Error & err) {
                          cerr << err.what() << endl;
			  cerr << " err code " << err.errcode() << endl;
                          return err.errcode();
                        }
                         catch( ...) { cerr << "Strange catch exception ???\n"; 
                          cerr << " at exec line  " << TheCurrentLine << endl;
                          return 1; 
                         }

                        if(verbosity)  cout << "times: compile "<< CPUcompile-CPUcompileInit <<"s, execution " <<  CPUtime()-CPUcompile << "s\n";
                        deleteStack(stack);
                        //debugstack.clear() 
                        } 
                        fingraphique();
                        NbPtr = ShowAlloc("end execution -- ",lg1) - NbPtr;
                        
                        if (NbPtr) { cout << " ######## We forget of deleting   " << NbPtr << " Nb pointer,   " <<  lg1-lg0 << "Bytes\n" ;}
  return 0;}
;

input:   instructions 
;
         
instructions:  instruction   {$$=$1;;;}
        | instructions  instruction   { $$= ($1+=$2) }         
        ;

list_of_id_args:   { $$=new ListOfId();}
            | id                      { $$ = new ListOfId(); $$->push_back(UnId($1))}
            | id '=' no_comma_expr    { $$ = new ListOfId(); $$->push_back(UnId($1,$3)) }
            | FESPACE id              { $$ = new ListOfId(); $$->push_back(UnId($2,Find($1),atype<FE<double> **>()))}
            | FESPACE '&' id              { $$ = new ListOfId(); $$->push_back(UnId($3,Find($1),atype<FE<double> **>(),true))}
            | type_of_dcl id          { $$ = new ListOfId(); $$->push_back(UnId($2,C_F0(),$1->right())) }
            | type_of_dcl '&' id      { $$ = new ListOfId(); $$->push_back(UnId($3,C_F0(),$1,true)) }
            | '[' list_of_id_args ']' { $$ = new ListOfId(); $$->push_back(UnId($2)) }
            | list_of_id_args ',' id                     { $$ = $1; $$->push_back(UnId($3)) }
            | list_of_id_args ',''[' list_of_id_args ']' { $$ = $1; $$->push_back(UnId($4)) }
            | list_of_id_args ',' id '=' no_comma_expr   { $$ = $1; $$->push_back(UnId($3,$5)) }
            | list_of_id_args ',' FESPACE id             { $$ = $1; $$->push_back(UnId($4,Find($3),atype<FE<double> **>())) }
            | list_of_id_args ',' FESPACE '&' id             { $$ = $1; $$->push_back(UnId($5,Find($3),atype<FE<double> **>(),true)) }
            | list_of_id_args ',' type_of_dcl id         { $$ = $1; $$->push_back(UnId($4,C_F0(),$3->right())) }
            | list_of_id_args ',' type_of_dcl '&' id     { $$ = $1; $$->push_back(UnId($5,C_F0(),$3,true)) }
;

list_of_id1:  id                      { $$ = new ListOfId(); $$->push_back(UnId($1)); }
            | list_of_id1 ',' id      { $$=$1  ; $$->push_back(UnId($3)); }
;
         
id: ID | FESPACE; 

list_of_dcls:    ID                         {$$=currentblock->NewVar<LocalVariable>($1,dcltype)}   
              |  ID '='   no_comma_expr     {$$=currentblock->NewVar<LocalVariable>($1,dcltype,$3)}
              |  ID  '(' parameters_list ')' {$$=currentblock->NewVar<LocalVariable>($1,dcltype,$3);
                                              $3.destroy()}
              |  list_of_dcls ',' list_of_dcls  {$$=C_F0($1,$3)}
                           
              
;


parameters_list:
	   no_set_expr {$$=$1} 
	|  FESPACE  ID  {$$=Find($1)} 
	|  ID '=' no_set_expr { $$=make_pair<const char *,const C_F0>($1,$3)} 	
	| parameters_list ',' no_set_expr { $$ = ($1 += $3) }
	| parameters_list ',' id '=' no_set_expr { $$= ($1+= make_pair<const char *,const C_F0>($3,$5))}
; 

type_of_dcl:   TYPE 
             | TYPE '[' TYPE ']' {$$=TypeArray($1,$3)}
             | TYPE '[' TYPE ',' TYPE ']' {$$=TypeArray($1,$3,$5)}
             | TYPE '<' TYPE '>' {$$=TypeTemplate($1,$3)}
             
             
;


ID_space:
    ID                                  { $$ =  NewFEvariable($1,currentblock,fespacetype,fespacecomplex); }
 |  ID '[' no_set_expr ']'              { $$ =  NewFEarray($1,currentblock,fespacetype,$3,fespacecomplex); }
 |  ID '=' no_set_expr                  { $$ =  NewFEvariable($1,currentblock,fespacetype,$3,fespacecomplex) }
 |  '[' list_of_id1 ']'                 { $$ =  NewFEvariable($2,currentblock,fespacetype,fespacecomplex) }
 |  '[' list_of_id1 ']' '[' no_set_expr ']'  { $$ =  NewFEarray($2,currentblock,fespacetype,$5,fespacecomplex) }
 |  '[' list_of_id1 ']' '=' no_set_expr { $$ =  NewFEvariable($2,currentblock,fespacetype,$5,fespacecomplex) }
; 
ID_array_space:
    ID '(' no_set_expr ')'              { $$ =  NewFEarray($1,currentblock,fespacetype,$3,fespacecomplex); }
 |  '[' list_of_id1 ']' '(' no_set_expr ')'  { $$ =  NewFEarray($2,currentblock,fespacetype,$5,fespacecomplex) }

;

fespace:  FESPACE {fespacecomplex=false;  fespacetype = Find($1);}
        | FESPACE '<' TYPE '>' {
             if ($3 != typevarreal && $3 != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=($3==typevarcomplex);
             fespacetype = Find($1);}
;        
spaceIDa  :      ID_array_space {  $$ = $1  } 
            |    spaceIDa ',' ID_array_space { $$=C_F0($1,$3);} ;
            
spaceIDb  :      ID_space {  $$ = $1  } 
            |    spaceIDb ',' ID_space { $$=C_F0($1,$3);} ;

spaceIDs :    fespace               spaceIDb    { $$=0;  $$ = $2} 
           |  fespace '[' TYPE ']'  spaceIDa    { $$=0;  $$ = $5}
;

fespace_def:
   ID '(' parameters_list ')' 
     {$$=currentblock->NewVar<LocalVariableFES,size_t>($1,atype<pfes*>(),$3,dimFESpaceImage($3));
     $3.destroy(); };
     
fespace_def_list:  fespace_def
                 | fespace_def_list ',' fespace_def {$$=C_F0($1,$3)}
;
    
declaration:   type_of_dcl {dcltype=$1} list_of_dcls ';' {$$=$3} 
             | FESPACEID  fespace_def_list    ';' {$$=$2}  
             | spaceIDs ';'{ $$=$1} 
             | FUNCTION ID '=' Expr ';' {$$=currentblock->NewID($1,$2,$4);} 
             | FUNCTION type_of_dcl ID  '(' list_of_id_args ')' 
                   {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = $2->right();
                      routineinblock[kkembtype] = currentblock;
                      $<routine>5=new Routine($1,$2->right(),$3,$5,currentblock);
                     // cout << " \n after new routine \n " << endl;                      
                      }
                    '{' instructions'}' 
                     { currentblock=$<routine>5->Set($9);
                       currentblock->Add($3,"(",$<routine>5);
                       kkembtype--;
                       $$=0;
                    
                        }
             | FUNCTION ID '(' list_of_id_args ')' 
                      {Block::open(currentblock); $1->SetArgs($4);}
                   '='   no_comma_expr  ';' 
                      {  $<cinst>$=currentblock->close(currentblock);
                         $$=currentblock->NewID($1,$2,$8,*$4);
                         delete $4; //  FH 23032005
                         } 
;              

begin: '{'  {  Block::open(currentblock)};
end:   '}'  {  $$=currentblock->close(currentblock)};

for_loop:  FOR {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;};
while_loop:  WHILE {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;}; // modif FH july 2005
                

declaration_for: 
    type_of_dcl {dcltype=$1; Block::open(currentblock);  }
               list_of_dcls {$$=$3};

try: TRY { Block::open(currentblock) };

instruction:   ';' {$$=0;} 
         | INCLUDE  STRING  {zzzfff->input($2);$$= 0; }
         | LOAD  STRING  {load($2);$$= 0; }
         |  try  '{' instructions '}' catchs {$$=Try($3,$5,currentblock->close(currentblock));}
         |  Expr  ';' {$$=$1}  
         |  declaration  {$$=$1} 
         |  for_loop  '(' Expr ';' Expr ';' Expr ')' instruction {inloopcount--; $$=For($3,$5,$7,$9)} 
         |  for_loop  '(' declaration_for ';' Expr ';' Expr ')' instruction 
                {inloopcount--; 
                $$=C_F0(For($3,$5,$7,$9),currentblock->close(currentblock))}
                 
         |  while_loop '(' Expr ')' instruction {inloopcount--;$$=While($3,$5)}
         |  IF '(' Expr ')'   instruction  {$$=FIf($3,$5)}
         |  IF '(' Expr ')'   instruction  ELSE instruction {$$=FIf($3,$5,$7)}
         |  begin  instructions end { 
                      $$=C_F0(new E_block($2,$3),atype<void>()) }
         |  BORDER  ID   border_expr {
                      $$=0;currentblock->NewID(atype<const E_Border *>(),$2,C_F0(TheOperators,"[border]",$3))} 
         |  BORDER  ID   '['  array ']' ';' {
                      $$=0;currentblock->NewID(atype<const E_Border *>(),$2,C_F0(TheOperators,"[border]",$4))} 
                               
         |  BREAK ';' {
                    if(inloopcount) 
                      $$= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") }
         |  CONTINUE ';' { 
                    if(inloopcount)
                        $$= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop")}
         |  RETURN  Expr ';' { 
                    if (kkembtype>=0)
                      $$= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo($2)) ,atype<void>());
                     else lgerror(" return not in routine ") }

;
catchs:  
  CATCH '(' '.' '.' '.'  ')' instruction {$$ =  $7; }
;

bornes: '(' ID '=' Expr ',' Expr ')' { 
   Block::open(currentblock);
   $$ = currentblock->NewVar<LocalVariable>($2,atype<double*>());
   $$+= $4;
   $$+= $6 }
;
border_expr:   bornes instruction {   
   $$ = ($1 += $2);
   currentblock->close(currentblock)} 
 ;

Expr:	 
         no_comma_expr 
       | Expr ',' Expr {$$=C_F0(TheOperators,$2,$1,$3);}
;

	
unop:
	  '-' 
	| '+' 
	| '!' 
	| PLUSPLUS 
	| MOINSMOINS 
;

no_comma_expr:
      no_set_expr 
	| no_set_expr '=' no_comma_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_set_expr PLUSEQ no_comma_expr {$$=C_F0(TheOperators,"+=",$1,$3)}
	| no_set_expr MOINSEQ no_comma_expr {$$=C_F0(TheOperators,"-=",$1,$3)}
	| no_set_expr MULEQ no_comma_expr {$$=C_F0(TheOperators,"*=",$1,$3)}
	| no_set_expr DIVEQ no_comma_expr {$$=C_F0(TheOperators,"/=",$1,$3)}
;


no_set_expr:
	no_ternary_expr
	| no_ternary_expr '?' no_set_expr ':' no_set_expr {$$=C_F0(TheOperators,"?:",$1,$3,$5)}
;
no_ternary_expr:
	  unary_expr 
	| no_ternary_expr '*' no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr DOTSTAR no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr DOTSLASH no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr '/' no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr '%' no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr '+' no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr '-' no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr LTLT no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr GTGT no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr '&' no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr AND no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr '|' no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr OR no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr '<' no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr LE no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr '>' no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr GE no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr EQ no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	| no_ternary_expr NE no_ternary_expr {$$=C_F0(TheOperators,$2,$1,$3)}
	
;

sub_script_expr:  
	    no_set_expr {$$=$1} 
    |   ':' {$$=C_F0(TheOperators,":")}
	|   no_set_expr ':' no_set_expr {$$=C_F0(TheOperators,":",$1,$3)}
	|   no_set_expr ':' no_set_expr ':' no_set_expr {$$=C_F0(TheOperators,":",$1,$3,$5)} 	
;
  
parameters:  {$$=0} 
	|   FESPACE {$$=Find($1)} 
	|  id '=' no_set_expr { $$=make_pair<const char *,const C_F0>($1,$3)} 	
	|   sub_script_expr {$$=$1} 
	| parameters ',' FESPACE { $$ = ($1 += Find($3)) }
	| parameters ',' sub_script_expr { $$ = ($1 += $3) }
	| parameters ',' id '=' no_set_expr { $$= ($1+= make_pair<const char *,const C_F0>($3,$5)) } 
; 

array:   no_comma_expr {$$=$1} 
       | array ',' no_comma_expr {$$ = ($1 += $3) };
     
    
unary_expr:
    pow_expr   
  | unop  pow_expr %prec UNARY {$$=C_F0(TheOperators,$1,$2)} 
;   

pow_expr: primary
  |      primary  '^' unary_expr   {$$=C_F0(TheOperators,$2,$1,$3)} 
  |      primary  '_' unary_expr   {$$=C_F0(TheOperators,$2,$1,$3)} 
  |      primary '\''              {$$=C_F0(TheOperators,$2,$1)} 
;

primary:  
           ID           {$$=Find($1);}
  |        LNUM         {$$= CConstant($1)}
  |        DNUM         {$$= CConstant($1)}
  |        CNUM         {$$= CConstant(complex<double>(0,$1))}
  |        STRING {$$= CConstant<const char *>($1)}
  |        primary '('  parameters ')'  {$$=C_F0($1,$2,$3);}
  |        primary '[' sub_script_expr ']'    {$$=C_F0($1,$2,$3)}
  |        primary '[' sub_script_expr ',' sub_script_expr ']'  {$$=C_F0($1,$2,$3,$5)}
  |        primary '['  ']'        {$$=C_F0($1,"[]")}
  |        primary '.'  ID       { $$=C_F0($1,$3) ;}
  |        FESPACE '.'  ID       { $$=C_F0(Find($1),$3) ;}
  |        FESPACE '(' parameters ')'     { $$=C_F0(Find($1),$2,$3) ;}
  |        primary PLUSPLUS      {$$=C_F0(TheRightOperators,$2,$1)} 
  |        primary MOINSMOINS    {$$=C_F0(TheRightOperators,$2,$1)} 
  |        TYPE '('  Expr ')' {
             if ($1->right()->CastingFrom($3.left()) ) 
                $$=$1->right()->CastTo($3)  ;
             else { $$=$1->right()->Find("<--",basicAC_F0_wa($3));
             if (!$$.left()) { cerr << " no wait to change " << $3.left()->right()->name() << " in " << 
                                        $1->right()->name() << endl;
                                CompileError(" Error in type(exp) "); }
             }
            }
  |        '(' Expr ')' {$$=$2}
  |        '[' array  ']' { $$=C_F0(TheOperators,"[]",$2)} 

;


%% 


#include <fstream>
using namespace std;
// bool lgdebug;
// bool lexdebug;
void ForDebug();
void ForDebug()
{
  int i=0;
  i++;
}
//extern void ShowAlloc(const char *s, size_t lg);
//extern void ShowNbAlloc(const char *s);
void init_lgfem() ;
void init_lgmesh() ;
void init_lgmesh3() ;
void init_algo();
bool withrgraphique = false;
//string  StrVersionNumber();

int Compile()
{
  extern   YYSTYPE *plglval;  // modif FH 
  plglval = &lglval;
  int retvalue=0;
  //  int ok;
  
  currentblock=0;
  Block::open(currentblock);  
  try {
    retvalue=yyparse(); //  compile
    if(retvalue==0)  
      if(currentblock) 
	retvalue=1,cerr <<  "Error:a block is not close" << endl;   
      else {
	  if( verbosity) {
	      cerr << " CodeAlloc : nb ptr  "<< CodeAlloc::nb << ",  size :"  <<  CodeAlloc::lg << endl;
	      cerr <<  "Bien: On a fini Normalement" << endl; }
	}
  }

  catch (Error & e) 
    {
      retvalue=e.errcode();
      cerr << "error " << e.what() 
	   << "\n code = "<<  retvalue << endl;
    }
  catch(std::ios_base::failure & e)
    {
     cerr << "std  catch io failure \n what : " << e.what() << endl;; 
     cerr << " at exec line  " << TheCurrentLine << endl; 
    }
  catch(std::exception & e)
    {
     cerr << "std  catch exception \n what : " << e.what() << endl;; 
     cerr << " at exec line  " << TheCurrentLine << endl; 
    
    }
  catch(...)
   {
     cerr << "Strange catch exception ???\n"; 
     cerr << " at exec line  " << TheCurrentLine << endl; 
    }
  return retvalue; 
}

int mainff (int  argc, char **argv)
{
GetEnvironment();     
  //  size_t lg000;
 // ShowAlloc("begin main ",lg000);
  int retvalue=0;
#ifdef PARALLELE
   initparallele(argc,argv);
#endif
  CPUcompileInit= CPUtime();
  withrgraphique = false;
   atexit(ForDebug);
//  AllFunctions::maptype  xlocal;
//  local=&xlocal;
  lexdebug = false;
  lgdebug = false;

  char *  cc= new char [1024];
  //  istream * ccin=0;
  if ( ! (getprog(cc,argc,argv)>0) ) 
    return 1; 
  if(verbosity) { 
      cout << "-- FreeFem++ v" << StrVersionNumber() << endl;
      if(verbosity>1) cout << "   file :" << cc << " " << " verbosity= " << verbosity << endl;
  }
  
  zzzfff = Newlex(cout);
    
  
/*  
  ccin= new ifstream(cc);
  if (argc >1 && (ccin!=0) )  
     ccin= new ifstream(argv[1]),throwassert(ccin);
  if (ccin!=0) 
    zzzfff = new  mylex(*ccin,cout) ;
  else 
    zzzfff = new  mylex(cin,cout) ;
*/    
//  les motsclefs    
   zzzfff->Add("include",INCLUDE);
   zzzfff->Add("load",LOAD);
   zzzfff->Add("while",WHILE);
   zzzfff->Add("for",FOR);
   zzzfff->Add("if",IF);
   zzzfff->Add("else",ELSE);
   zzzfff->Add("end",ENDOFFILE);
   zzzfff->Add("break",BREAK);
   zzzfff->Add("continue",CONTINUE);
   zzzfff->Add("return",RETURN);
   zzzfff->Add("border",BORDER);
   zzzfff->Add("fespace",FESPACEID);
   zzzfff->Add("try",TRY);
   zzzfff->Add("catch",CATCH);
   zzzfff->Add("throw",THROW);
   Init_map_type();
   if(verbosity) cout << " Load: ";
   init_lgfem() ;
   init_lgmesh() ;
   init_lgmesh3() ;
   init_algo();
   
#ifdef HAVE_LIBARPACK
   init_eigenvalue();
#endif   
#ifdef PARALLELE
   init_lgparallele(); 
#endif 
#ifdef HAVE_LIBUMFPACK   
  if(verbosity)  cout << " UMFPACK ";  
#endif
 // callInitsFunct(); Pb opimisation 
  if(verbosity)  cout << endl;
  zzzfff->input(cc);
  EnvironmentLoad(); // just before compile     
  retvalue= Compile(); 
      
#ifdef PARALLELE
  end_parallele();
#endif
  //  currentblock->close(currentblock).eval(thestack);
  fingraphique();
  Destroylex( zzzfff);
  
   // ClearMem();
  return retvalue;
}


 
