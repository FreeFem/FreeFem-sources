
/*  A Bison parser, made from lg.y
    by GNU Bison version 1.28  */

#define YYBISON 1  /* Identify Bison output.  */

#define yyparse lgparse
#define yylex lglex
#define yyerror lgerror
#define yylval lglval
#define yychar lgchar
#define yydebug lgdebug
#define yynerrs lgnerrs
#define	IF	257
#define	ELSE	258
#define	SET	259
#define	LTLT	260
#define	GTGT	261
#define	OR	262
#define	AND	263
#define	EQ	264
#define	NE	265
#define	LE	266
#define	GE	267
#define	DOTSTAR	268
#define	DOTSLASH	269
#define	UNARY	270
#define	PLUSPLUS	271
#define	MOINSMOINS	272
#define	LNUM	273
#define	DNUM	274
#define	CNUM	275
#define	ID	276
#define	FESPACEID	277
#define	IDPARAM	278
#define	STRING	279
#define	ENDOFFILE	280
#define	INCLUDE	281
#define	LOAD	282
#define	BIDON	283
#define	FOR	284
#define	WHILE	285
#define	BREAK	286
#define	CONTINUE	287
#define	RETURN	288
#define	TRY	289
#define	CATCH	290
#define	THROW	291
#define	TYPE	292
#define	FUNCTION	293
#define	FESPACE	294
#define	PLUSEQ	295
#define	MOINSEQ	296
#define	MULEQ	297
#define	DIVEQ	298
#define	ARROW	299
#define	BORDER	300
#define	CURVE	301
#define	SOLVE	302

#line 1 "lg.y"
 

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

int ShowAlloc(char *s,size_t &);
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


#line 77 "lg.y"
typedef union{ 
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
} YYSTYPE;
#ifndef YYDEBUG
#define YYDEBUG 1
#endif

#include <stdio.h>

#ifndef __cplusplus
#ifndef __STDC__
#define const
#endif
#endif



#define	YYFINAL		360
#define	YYFLAG		-32768
#define	YYNTBASE	74

#define YYTRANSLATE(x) ((unsigned)(x) <= 302 ? yytranslate[x] : 117)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    30,     2,     2,     2,    24,    13,    32,    34,
    37,    22,    20,     5,    21,    36,    23,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,    73,    69,    16,
     6,    17,    72,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
    35,     2,    38,    31,    33,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    70,    11,    71,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     1,     3,     4,     7,     8,
     9,    10,    12,    14,    15,    18,    19,    25,    26,    27,
    28,    29,    39,    40,    41,    42,    43,    44,    45,    46,
    47,    48,    49,    50,    51,    52,    53,    54,    55,    56,
    57,    58,    59,    60,    61,    62,    63,    64,    65,    66,
    67,    68
};

#if YYDEBUG != 0
static const short yyprhs[] = {     0,
     0,     3,     5,     7,    10,    11,    13,    17,    20,    24,
    27,    31,    35,    39,    45,    51,    56,    62,    67,    73,
    75,    79,    81,    83,    85,    89,    94,    98,   100,   103,
   107,   111,   117,   119,   124,   131,   136,   138,   143,   147,
   151,   158,   164,   169,   176,   178,   183,   185,   189,   191,
   195,   198,   204,   209,   211,   215,   216,   221,   225,   228,
   234,   235,   246,   247,   257,   259,   261,   263,   265,   266,
   270,   272,   274,   277,   280,   286,   289,   291,   301,   311,
   317,   323,   331,   335,   339,   346,   349,   352,   356,   364,
   372,   375,   377,   381,   383,   385,   387,   389,   391,   393,
   397,   401,   405,   409,   413,   415,   421,   423,   427,   431,
   435,   439,   443,   447,   451,   455,   459,   463,   467,   471,
   475,   479,   483,   487,   491,   495,   499,   501,   503,   507,
   513,   514,   516,   520,   522,   526,   530,   536,   538,   542,
   544,   547,   549,   553,   557,   560,   562,   564,   566,   568,
   570,   575,   580,   587,   591,   595,   599,   604,   607,   610,
   615,   619
};

static const short yyrhs[] = {    75,
    46,     0,    76,     0,   102,     0,    76,   102,     0,     0,
    79,     0,    79,     6,   108,     0,    60,    79,     0,    60,
    13,    79,     0,    82,    79,     0,    82,    13,    79,     0,
    35,    77,    38,     0,    77,     5,    79,     0,    77,     5,
    35,    77,    38,     0,    77,     5,    79,     6,   108,     0,
    77,     5,    60,    79,     0,    77,     5,    60,    13,    79,
     0,    77,     5,    82,    79,     0,    77,     5,    82,    13,
    79,     0,    79,     0,    78,     5,    79,     0,    42,     0,
    60,     0,    42,     0,    42,     6,   108,     0,    42,    34,
    81,    37,     0,    80,     5,    80,     0,   109,     0,    60,
    42,     0,    42,     6,   109,     0,    81,     5,   109,     0,
    81,     5,    79,     6,   109,     0,    58,     0,    58,    35,
    58,    38,     0,    58,    35,    58,     5,    58,    38,     0,
    58,    16,    58,    17,     0,    42,     0,    42,    35,   109,
    38,     0,    42,     6,   109,     0,    35,    78,    38,     0,
    35,    78,    38,    35,   109,    38,     0,    35,    78,    38,
     6,   109,     0,    42,    34,   109,    37,     0,    35,    78,
    38,    34,   109,    37,     0,    60,     0,    60,    16,    58,
    17,     0,    84,     0,    86,     5,    84,     0,    83,     0,
    87,     5,    83,     0,    85,    87,     0,    85,    35,    58,
    38,    86,     0,    42,    34,    81,    37,     0,    89,     0,
    90,     5,    89,     0,     0,    82,    92,    80,    69,     0,
    43,    90,    69,     0,    88,    69,     0,    59,    42,     6,
   106,    69,     0,     0,    59,    82,    42,    34,    77,    37,
    93,    70,    76,    71,     0,     0,    59,    42,    34,    77,
    37,    94,     6,   108,    69,     0,    70,     0,    71,     0,
    50,     0,    51,     0,     0,    82,   100,    80,     0,    55,
     0,    69,     0,    47,    45,     0,    48,    45,     0,   101,
    70,    76,    71,   103,     0,   106,    69,     0,    91,     0,
    97,    34,   106,    69,   106,    69,   106,    37,   102,     0,
    97,    34,    99,    69,   106,    69,   106,    37,   102,     0,
    98,    34,   106,    37,   102,     0,     3,    34,   106,    37,
   102,     0,     3,    34,   106,    37,   102,     4,   102,     0,
    95,    76,    96,     0,    66,    42,   105,     0,    66,    42,
    35,   113,    38,    69,     0,    52,    69,     0,    53,    69,
     0,    54,   106,    69,     0,    56,    34,    36,    36,    36,
    37,   102,     0,    34,    42,     6,   106,     5,   106,    37,
     0,   104,   102,     0,   108,     0,   106,     5,   106,     0,
    21,     0,    20,     0,    30,     0,    28,     0,    29,     0,
   109,     0,   109,     6,   108,     0,   109,    61,   108,     0,
   109,    62,   108,     0,   109,    63,   108,     0,   109,    64,
   108,     0,   110,     0,   110,    72,   109,    73,   109,     0,
   114,     0,   110,    22,   110,     0,   110,    25,   110,     0,
   110,    26,   110,     0,   110,    23,   110,     0,   110,    24,
   110,     0,   110,    20,   110,     0,   110,    21,   110,     0,
   110,     8,   110,     0,   110,     9,   110,     0,   110,    13,
   110,     0,   110,    12,   110,     0,   110,    11,   110,     0,
   110,    10,   110,     0,   110,    16,   110,     0,   110,    18,
   110,     0,   110,    17,   110,     0,   110,    19,   110,     0,
   110,    14,   110,     0,   110,    15,   110,     0,   109,     0,
    73,     0,   109,    73,   109,     0,   109,    73,   109,    73,
   109,     0,     0,    60,     0,    79,     6,   109,     0,   111,
     0,   112,     5,    60,     0,   112,     5,   111,     0,   112,
     5,    79,     6,   109,     0,   108,     0,   113,     5,   108,
     0,   115,     0,   107,   115,     0,   116,     0,   116,    31,
   114,     0,   116,    33,   114,     0,   116,    32,     0,    42,
     0,    39,     0,    40,     0,    41,     0,    45,     0,   116,
    34,   112,    37,     0,   116,    35,   111,    38,     0,   116,
    35,   111,     5,   111,    38,     0,   116,    35,    38,     0,
   116,    36,    42,     0,    60,    36,    42,     0,    60,    34,
   112,    37,     0,   116,    28,     0,   116,    29,     0,    58,
    34,   106,    37,     0,    34,   106,    37,     0,    35,   113,
    38,     0
};

#endif

#if YYDEBUG != 0
static const short yyrline[] = { 0,
   203,   239,   242,   243,   246,   247,   248,   249,   250,   251,
   252,   253,   254,   255,   256,   257,   258,   259,   260,   263,
   264,   267,   267,   269,   270,   271,   273,   279,   281,   282,
   283,   284,   287,   288,   289,   290,   296,   298,   299,   300,
   301,   302,   304,   306,   310,   311,   316,   317,   319,   320,
   322,   323,   326,   331,   332,   335,   335,   336,   337,   338,
   339,   347,   354,   356,   363,   364,   366,   368,   372,   374,
   376,   378,   379,   380,   381,   382,   383,   384,   385,   389,
   390,   391,   392,   394,   396,   399,   403,   407,   413,   417,
   423,   428,   430,   434,   436,   437,   438,   439,   442,   444,
   445,   446,   447,   448,   452,   454,   456,   458,   459,   460,
   461,   462,   463,   464,   465,   466,   467,   468,   469,   470,
   471,   472,   473,   474,   475,   476,   480,   482,   483,   484,
   487,   488,   489,   490,   491,   492,   493,   496,   497,   500,
   502,   505,   506,   507,   508,   511,   513,   514,   515,   516,
   517,   518,   519,   520,   521,   522,   523,   524,   525,   526,
   535,   536
};
#endif


#if YYDEBUG != 0 || defined (YYERROR_VERBOSE)

static const char * const yytname[] = {   "$","error","$undefined.","IF","ELSE",
"','","'='","SET","LTLT","GTGT","OR","'|'","AND","'&'","EQ","NE","'<'","'>'",
"LE","GE","'+'","'-'","'*'","'/'","'%'","DOTSTAR","DOTSLASH","UNARY","PLUSPLUS",
"MOINSMOINS","'!'","'^'","'\\''","'_'","'('","'['","'.'","')'","']'","LNUM",
"DNUM","CNUM","ID","FESPACEID","IDPARAM","STRING","ENDOFFILE","INCLUDE","LOAD",
"BIDON","FOR","WHILE","BREAK","CONTINUE","RETURN","TRY","CATCH","THROW","TYPE",
"FUNCTION","FESPACE","PLUSEQ","MOINSEQ","MULEQ","DIVEQ","ARROW","BORDER","CURVE",
"SOLVE","';'","'{'","'}'","'?'","':'","start","input","instructions","list_of_id_args",
"list_of_id1","id","list_of_dcls","parameters_list","type_of_dcl","ID_space",
"ID_array_space","fespace","spaceIDa","spaceIDb","spaceIDs","fespace_def","fespace_def_list",
"declaration","@1","@2","@3","begin","end","for_loop","while_loop","declaration_for",
"@4","try","instruction","catchs","bornes","border_expr","Expr","unop","no_comma_expr",
"no_set_expr","no_ternary_expr","sub_script_expr","parameters","array","unary_expr",
"pow_expr","primary", NULL
};
#endif

static const short yyr1[] = {     0,
    74,    75,    76,    76,    77,    77,    77,    77,    77,    77,
    77,    77,    77,    77,    77,    77,    77,    77,    77,    78,
    78,    79,    79,    80,    80,    80,    80,    81,    81,    81,
    81,    81,    82,    82,    82,    82,    83,    83,    83,    83,
    83,    83,    84,    84,    85,    85,    86,    86,    87,    87,
    88,    88,    89,    90,    90,    92,    91,    91,    91,    91,
    93,    91,    94,    91,    95,    96,    97,    98,   100,    99,
   101,   102,   102,   102,   102,   102,   102,   102,   102,   102,
   102,   102,   102,   102,   102,   102,   102,   102,   103,   104,
   105,   106,   106,   107,   107,   107,   107,   107,   108,   108,
   108,   108,   108,   108,   109,   109,   110,   110,   110,   110,
   110,   110,   110,   110,   110,   110,   110,   110,   110,   110,
   110,   110,   110,   110,   110,   110,   111,   111,   111,   111,
   112,   112,   112,   112,   112,   112,   112,   113,   113,   114,
   114,   115,   115,   115,   115,   116,   116,   116,   116,   116,
   116,   116,   116,   116,   116,   116,   116,   116,   116,   116,
   116,   116
};

static const short yyr2[] = {     0,
     2,     1,     1,     2,     0,     1,     3,     2,     3,     2,
     3,     3,     3,     5,     5,     4,     5,     4,     5,     1,
     3,     1,     1,     1,     3,     4,     3,     1,     2,     3,
     3,     5,     1,     4,     6,     4,     1,     4,     3,     3,
     6,     5,     4,     6,     1,     4,     1,     3,     1,     3,
     2,     5,     4,     1,     3,     0,     4,     3,     2,     5,
     0,    10,     0,     9,     1,     1,     1,     1,     0,     3,
     1,     1,     2,     2,     5,     2,     1,     9,     9,     5,
     5,     7,     3,     3,     6,     2,     2,     3,     7,     7,
     2,     1,     3,     1,     1,     1,     1,     1,     1,     3,
     3,     3,     3,     3,     1,     5,     1,     3,     3,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     3,     3,     3,     3,     3,     1,     1,     3,     5,
     0,     1,     3,     1,     3,     3,     5,     1,     3,     1,
     2,     1,     3,     3,     2,     1,     1,     1,     1,     1,
     4,     4,     6,     3,     3,     3,     4,     2,     2,     4,
     3,     3
};

static const short yydefact[] = {     0,
     0,    95,    94,    97,    98,    96,     0,     0,   147,   148,
   149,   146,     0,   150,     0,     0,    67,    68,     0,     0,
     0,    71,    33,     0,    45,     0,    72,    65,     0,     2,
    56,     0,     0,    77,     0,     0,     0,     0,     3,     0,
     0,    92,    99,   105,   107,   140,   142,     0,     0,     0,
     0,   138,     0,     0,    54,     0,    73,    74,    86,    87,
     0,     0,     0,     0,     0,    33,     0,     0,   131,     0,
     0,     1,     4,     0,     0,    37,    49,    51,    59,     0,
     0,     0,     0,     0,    76,   141,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,   158,   159,     0,   145,     0,   131,     0,     0,     0,
   161,     0,   162,     0,     0,    58,    88,     0,     0,     0,
     0,     5,     0,     0,   146,   132,   128,     0,   127,   134,
     0,   156,     0,     0,     0,    84,    24,     0,    22,     0,
    23,     0,    20,     0,     0,     0,    66,    83,    69,     0,
     0,     0,     0,    93,   100,   101,   102,   103,   104,   115,
   116,   120,   119,   118,   117,   125,   126,   121,   123,   122,
   124,   113,   114,   108,   111,   112,   109,   110,     0,   143,
   144,     0,   154,     0,   155,     0,   139,   146,     0,     0,
    28,    55,    36,   160,     0,    34,     0,     5,    23,     0,
     6,     0,     5,    46,     0,     0,     0,   157,     0,     0,
    91,     0,     0,     0,    57,     0,     0,    40,    39,     0,
     0,    50,     0,     0,     0,     0,     0,     0,   151,     0,
   152,    81,     0,    29,     0,    53,     0,    60,     0,     0,
     8,     0,    63,     0,     0,    10,     0,   133,   129,   135,
     0,   136,     0,     0,    25,     0,    27,     0,     0,    47,
    52,    21,     0,     0,    38,    70,     0,     0,    80,     0,
    75,   106,     0,     0,    30,    23,     0,    31,    35,    12,
     9,     5,    23,    13,     0,     0,     7,    11,    61,     0,
     0,     0,    85,    26,     0,     0,     0,    42,     0,     0,
     0,     0,   153,    82,     0,     0,     0,    16,     0,     0,
    18,     0,     0,   130,   137,     0,     0,     0,    48,    41,
     0,     0,     0,    32,    14,    17,    15,    19,     0,     0,
    93,     0,    43,     0,     0,     0,    64,     0,    90,     0,
    79,    78,     0,    62,    44,     0,    89,     0,     0,     0
};

static const short yydefgoto[] = {   358,
    29,    30,   210,   152,   211,   148,   200,    31,    77,   270,
    32,   271,    78,    33,    55,    56,    34,    74,   323,   296,
    35,   158,    36,    37,   160,   233,    38,    39,   281,   145,
   146,    40,    41,    42,    43,    44,   140,   141,    53,    45,
    46,    47
};

static const short yypact[] = {   462,
   -10,-32768,-32768,-32768,-32768,-32768,   586,   586,-32768,-32768,
-32768,-32768,    70,-32768,    -9,    77,-32768,-32768,    58,    78,
   586,-32768,    10,    68,    87,   123,-32768,-32768,    -7,   462,
-32768,    -1,   107,-32768,   462,   143,   148,   118,-32768,     3,
   182,-32768,    31,   275,-32768,-32768,   682,   586,   164,   165,
    23,-32768,    17,   176,-32768,     4,-32768,-32768,-32768,-32768,
     5,   153,   586,   155,   112,   146,   173,   160,   114,   183,
   136,-32768,-32768,   184,   126,    19,-32768,   226,-32768,   303,
   613,   586,   462,   586,-32768,-32768,   586,   586,   586,   586,
   586,   586,   586,   586,   586,   586,   586,   586,   586,   586,
   586,   586,   586,   586,   586,   586,   586,   586,   586,   586,
   586,-32768,-32768,   586,-32768,   586,   114,   505,   190,    42,
-32768,   586,-32768,   640,    70,-32768,-32768,   216,    91,    18,
   586,   154,   204,   222,   237,   157,-32768,   238,   175,-32768,
    92,-32768,   207,   586,   462,-32768,   151,     8,-32768,   208,
-32768,    24,-32768,   586,   586,    74,-32768,-32768,-32768,   187,
    11,    93,   356,-32768,-32768,-32768,-32768,-32768,-32768,   737,
   737,   752,   752,   765,   765,   242,   242,   772,   772,   772,
   772,   256,   256,-32768,-32768,-32768,-32768,-32768,   179,-32768,
-32768,    95,-32768,    26,-32768,   462,-32768,   264,   158,    99,
-32768,-32768,-32768,-32768,   213,-32768,    12,   154,    -2,   100,
   267,     1,   154,-32768,   586,   586,   532,-32768,   269,    27,
-32768,   586,   640,   184,-32768,   125,   -12,    96,-32768,   265,
   -12,-32768,   184,   586,   586,   462,   220,   586,-32768,   559,
-32768,   298,   586,-32768,   667,-32768,   271,-32768,    28,   -12,
-32768,   195,-32768,   586,   -12,-32768,   101,-32768,   234,   157,
   304,-32768,   586,   243,-32768,   102,-32768,   -12,   277,-32768,
   308,-32768,   586,   586,-32768,   312,    13,    14,-32768,   287,
-32768,-32768,   284,   462,-32768,   165,   319,-32768,-32768,-32768,
-32768,   154,     7,   321,    57,   322,-32768,-32768,-32768,   586,
   586,   324,-32768,-32768,    30,   586,   125,-32768,   292,   586,
   586,   299,-32768,-32768,   586,    33,   -12,-32768,   586,   -12,
-32768,   586,   266,-32768,-32768,   586,   300,   302,-32768,-32768,
   103,   108,   305,-32768,-32768,-32768,-32768,-32768,   280,   462,
   315,   586,-32768,   462,   462,   328,-32768,   409,-32768,   323,
-32768,-32768,   329,-32768,-32768,   462,-32768,   340,   365,-32768
};

static const short yypgoto[] = {-32768,
-32768,   -32,  -201,   110,   -48,  -113,   144,   -18,   212,    63,
-32768,-32768,-32768,-32768,   246,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,   -30,-32768,-32768,
-32768,    -6,-32768,    -3,   -65,   636,  -116,   258,   235,    89,
   339,-32768
};


#define	YYLAST		798


static const short yytable[] = {    73,
    51,   194,    80,   139,    52,    67,   249,    84,   125,    84,
   250,   257,   224,   255,    61,    84,    84,    84,    84,   317,
   138,   122,   205,    48,   154,    62,   153,    84,   227,   149,
   240,   122,   252,    75,   227,    57,    87,   252,    72,   149,
    76,   120,   149,    63,    64,   189,    84,   151,   149,    73,
   163,   139,   139,   155,   123,   206,   129,   151,   201,   121,
   151,   228,   159,   241,   264,   290,   151,   327,   138,   320,
   335,    85,   126,   127,   161,   162,   225,   164,   196,   235,
   248,   310,   311,   165,   166,   167,   168,   169,   229,   230,
   316,    88,    89,    90,    91,    84,   217,    84,   149,   217,
   262,   273,    68,   245,   252,   252,   245,    84,   231,    65,
   267,    54,    84,   212,   221,    76,   151,   131,   197,   276,
    69,    58,    70,   283,   207,    66,    59,   204,   218,   236,
   274,   239,    73,     2,     3,   246,   253,   299,   304,   344,
    52,     4,     5,     6,   345,   132,    60,     7,     8,   258,
   259,   139,     9,    10,    11,   135,   222,   201,    14,   268,
   251,    62,   -23,   256,    71,   242,   269,   149,   261,   143,
   144,    49,   282,   136,   139,    79,    81,   285,   272,   288,
    64,    82,   153,   150,   223,   151,   137,    83,   208,   212,
    69,    69,    70,    70,   212,   149,   287,    63,    69,   244,
    70,   291,   190,   294,   191,   279,   298,   308,   309,   124,
   128,    66,   130,   209,   133,     7,     8,   134,   265,   153,
     9,    10,    11,    12,   142,   147,    14,   277,   278,   292,
   156,   195,   203,   295,   324,   325,   149,   213,   214,    49,
   328,    50,   -22,   215,   318,   226,   321,   216,   219,   334,
   297,   238,    66,   314,   293,   234,   302,   100,   101,   102,
   103,   104,   105,   106,   107,   108,   109,   110,   336,   243,
   247,   338,   254,   212,   263,   280,   350,   106,   107,   108,
   109,   110,    92,    93,    94,    95,    96,    97,    98,    99,
   100,   101,   102,   103,   104,   105,   106,   107,   108,   109,
   110,   284,   275,   331,   332,     1,   300,   348,   289,   301,
   306,   303,   307,   351,   352,   337,   224,    73,   339,   341,
   312,   313,     2,     3,   315,   357,   319,   322,   326,   330,
     4,     5,     6,   342,   333,   340,     7,     8,   343,   359,
   346,     9,    10,    11,    12,    13,   111,    14,   347,    15,
    16,   349,    17,    18,    19,    20,    21,    22,     1,   355,
    23,    24,    25,   353,   360,   356,   266,   232,    26,   329,
   202,    27,    28,   157,   192,     2,     3,   305,   220,    86,
     0,     0,     0,     4,     5,     6,     0,     0,     0,     7,
     8,     0,     0,     0,     9,    10,    11,    12,    13,     0,
    14,     0,    15,    16,     0,    17,    18,    19,    20,    21,
    22,     1,     0,    23,    24,    25,     0,     0,     0,     0,
     0,    26,     0,     0,    27,    28,   237,     0,     2,     3,
     0,     0,     0,     0,     0,     0,     4,     5,     6,     0,
     0,     0,     7,     8,     0,     0,     0,     9,    10,    11,
    12,    13,     0,    14,     0,    15,    16,     0,    17,    18,
    19,    20,    21,    22,     1,     0,    23,    24,    25,     0,
     0,     0,     0,     0,    26,     0,     0,    27,    28,   354,
     0,     2,     3,     0,     0,     0,     0,     0,     0,     4,
     5,     6,     0,     0,     0,     7,     8,     0,     0,     0,
     9,    10,    11,    12,    13,     0,    14,     0,    15,    16,
     0,    17,    18,    19,    20,    21,    22,     0,     0,    23,
    24,    25,     0,     0,     2,     3,     0,    26,     0,     0,
    27,    28,     4,     5,     6,     0,     0,     0,     7,     8,
     0,     0,   193,     9,    10,    11,    12,     0,     0,    14,
     0,     2,     3,     0,     0,     0,     0,     0,     0,     4,
     5,     6,    49,     0,    50,     7,     8,     0,     0,     0,
     9,    10,    11,   135,     0,     0,    14,   137,     2,     3,
     0,     0,     0,     0,     0,     0,     4,     5,     6,    49,
     0,   260,     7,     8,     0,     0,     0,     9,    10,    11,
    12,     0,     0,    14,   137,     2,     3,     0,     0,     0,
     0,     0,     0,     4,     5,     6,    49,     0,    50,     7,
     8,     0,     0,     0,     9,    10,    11,    12,     0,     0,
    14,   137,     2,     3,     0,     0,     0,     0,     0,     0,
     4,     5,     6,    49,     0,    50,     7,     8,     0,     0,
     0,     9,    10,    11,    12,     0,     0,    14,     0,     2,
     3,     0,     0,     0,     0,     0,     0,     4,     5,     6,
    23,     0,    50,     7,     8,     0,     0,     0,     9,    10,
    11,   198,     0,     0,    14,     0,     2,     3,     0,     0,
     0,     0,     0,     0,     4,     5,     6,    49,     0,   199,
     7,     8,     0,     0,     0,     9,    10,    11,   135,   112,
   113,    14,   114,   115,   116,   117,   118,   119,     0,     0,
     0,     0,     0,     0,    49,     0,   286,   170,   171,   172,
   173,   174,   175,   176,   177,   178,   179,   180,   181,   182,
   183,   184,   185,   186,   187,   188,    94,    95,    96,    97,
    98,    99,   100,   101,   102,   103,   104,   105,   106,   107,
   108,   109,   110,    96,    97,    98,    99,   100,   101,   102,
   103,   104,   105,   106,   107,   108,   109,   110,    98,    99,
   100,   101,   102,   103,   104,   105,   106,   107,   108,   109,
   110,   104,   105,   106,   107,   108,   109,   110
};

static const short yycheck[] = {    30,
     7,   118,    35,    69,     8,    24,   208,     5,     5,     5,
    13,   213,     5,    13,    21,     5,     5,     5,     5,    13,
    69,     5,     5,    34,     6,    16,    75,     5,     5,    42,
     5,     5,     5,    35,     5,    45,     6,     5,    46,    42,
    42,    48,    42,    34,    35,   111,     5,    60,    42,    80,
    83,   117,   118,    35,    38,    38,    63,    60,   124,    37,
    60,    38,    81,    38,    38,    38,    60,    38,   117,    13,
    38,    69,    69,    69,    81,    82,    69,    84,    37,    69,
    69,    69,    69,    87,    88,    89,    90,    91,   154,   155,
   292,    61,    62,    63,    64,     5,     5,     5,    42,     5,
   217,     6,    16,     5,     5,     5,     5,     5,    35,    42,
   224,    42,     5,   132,   145,    42,    60,     6,   122,   233,
    34,    45,    36,   240,   131,    58,    69,    37,    37,    37,
    35,    37,   163,    20,    21,    37,    37,    37,    37,    37,
   144,    28,    29,    30,    37,    34,    69,    34,    35,   215,
   216,   217,    39,    40,    41,    42,     6,   223,    45,    35,
   209,    16,     6,   212,    42,   196,    42,    42,   217,    34,
    35,    58,   238,    60,   240,    69,    34,   243,   227,   245,
    35,    34,   231,    58,    34,    60,    73,    70,    35,   208,
    34,    34,    36,    36,   213,    42,   245,    34,    34,    42,
    36,   250,   114,   252,   116,   236,   255,   273,   274,    34,
    58,    58,    58,    60,    42,    34,    35,    58,   222,   268,
    39,    40,    41,    42,    42,    42,    45,   234,   235,    35,
     5,    42,    17,   252,   300,   301,    42,    34,    17,    58,
   306,    60,     6,     6,   293,    38,   295,    73,    42,   315,
   254,    73,    58,   284,    60,    69,   263,    16,    17,    18,
    19,    20,    21,    22,    23,    24,    25,    26,   317,     6,
    58,   320,     6,   292,     6,    56,   342,    22,    23,    24,
    25,    26,     8,     9,    10,    11,    12,    13,    14,    15,
    16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
    26,     4,    38,   310,   311,     3,    73,   340,    38,     6,
    34,    69,     5,   344,   345,   319,     5,   348,   322,   326,
    34,    38,    20,    21,     6,   356,     6,     6,     5,    38,
    28,    29,    30,    34,    36,    70,    34,    35,    37,     0,
    36,    39,    40,    41,    42,    43,    72,    45,    69,    47,
    48,    37,    50,    51,    52,    53,    54,    55,     3,    37,
    58,    59,    60,    36,     0,    37,   223,   156,    66,   307,
   125,    69,    70,    71,   117,    20,    21,   268,   144,    41,
    -1,    -1,    -1,    28,    29,    30,    -1,    -1,    -1,    34,
    35,    -1,    -1,    -1,    39,    40,    41,    42,    43,    -1,
    45,    -1,    47,    48,    -1,    50,    51,    52,    53,    54,
    55,     3,    -1,    58,    59,    60,    -1,    -1,    -1,    -1,
    -1,    66,    -1,    -1,    69,    70,    71,    -1,    20,    21,
    -1,    -1,    -1,    -1,    -1,    -1,    28,    29,    30,    -1,
    -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
    42,    43,    -1,    45,    -1,    47,    48,    -1,    50,    51,
    52,    53,    54,    55,     3,    -1,    58,    59,    60,    -1,
    -1,    -1,    -1,    -1,    66,    -1,    -1,    69,    70,    71,
    -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    -1,    28,
    29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
    39,    40,    41,    42,    43,    -1,    45,    -1,    47,    48,
    -1,    50,    51,    52,    53,    54,    55,    -1,    -1,    58,
    59,    60,    -1,    -1,    20,    21,    -1,    66,    -1,    -1,
    69,    70,    28,    29,    30,    -1,    -1,    -1,    34,    35,
    -1,    -1,    38,    39,    40,    41,    42,    -1,    -1,    45,
    -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    -1,    28,
    29,    30,    58,    -1,    60,    34,    35,    -1,    -1,    -1,
    39,    40,    41,    42,    -1,    -1,    45,    73,    20,    21,
    -1,    -1,    -1,    -1,    -1,    -1,    28,    29,    30,    58,
    -1,    60,    34,    35,    -1,    -1,    -1,    39,    40,    41,
    42,    -1,    -1,    45,    73,    20,    21,    -1,    -1,    -1,
    -1,    -1,    -1,    28,    29,    30,    58,    -1,    60,    34,
    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
    45,    73,    20,    21,    -1,    -1,    -1,    -1,    -1,    -1,
    28,    29,    30,    58,    -1,    60,    34,    35,    -1,    -1,
    -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,    20,
    21,    -1,    -1,    -1,    -1,    -1,    -1,    28,    29,    30,
    58,    -1,    60,    34,    35,    -1,    -1,    -1,    39,    40,
    41,    42,    -1,    -1,    45,    -1,    20,    21,    -1,    -1,
    -1,    -1,    -1,    -1,    28,    29,    30,    58,    -1,    60,
    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    28,
    29,    45,    31,    32,    33,    34,    35,    36,    -1,    -1,
    -1,    -1,    -1,    -1,    58,    -1,    60,    92,    93,    94,
    95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
   105,   106,   107,   108,   109,   110,    10,    11,    12,    13,
    14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
    24,    25,    26,    12,    13,    14,    15,    16,    17,    18,
    19,    20,    21,    22,    23,    24,    25,    26,    14,    15,
    16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
    26,    20,    21,    22,    23,    24,    25,    26
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/share/bison.simple"
/* This file comes from bison-1.28.  */

/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

#ifndef YYSTACK_USE_ALLOCA
#ifdef alloca
#define YYSTACK_USE_ALLOCA
#else /* alloca not defined */
#ifdef __GNUC__
#define YYSTACK_USE_ALLOCA
#define alloca __builtin_alloca
#else /* not GNU C.  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi) || (defined (__sun) && defined (__i386))
#define YYSTACK_USE_ALLOCA
#include <alloca.h>
#else /* not sparc */
/* We think this test detects Watcom and Microsoft C.  */
/* This used to test MSDOS, but that is a bad idea
   since that symbol is in the user namespace.  */
#if (defined (_MSDOS) || defined (_MSDOS_)) && !defined (__TURBOC__)
#if 0 /* No need for malloc.h, which pollutes the namespace;
	 instead, just don't use alloca.  */
#include <malloc.h>
#endif
#else /* not MSDOS, or __TURBOC__ */
#if defined(_AIX)
/* I don't know what this was needed for, but it pollutes the namespace.
   So I turned it off.   rms, 2 May 1997.  */
/* #include <malloc.h>  */
 #pragma alloca
#define YYSTACK_USE_ALLOCA
#else /* not MSDOS, or __TURBOC__, or _AIX */
#if 0
#ifdef __hpux /* haible@ilog.fr says this works for HPUX 9.05 and up,
		 and on HPUX 10.  Eventually we can turn this on.  */
#define YYSTACK_USE_ALLOCA
#define alloca __builtin_alloca
#endif /* __hpux */
#endif
#endif /* not _AIX */
#endif /* not MSDOS, or __TURBOC__ */
#endif /* not sparc */
#endif /* not GNU C */
#endif /* alloca not defined */
#endif /* YYSTACK_USE_ALLOCA not defined */

#ifdef YYSTACK_USE_ALLOCA
#define YYSTACK_ALLOC alloca
#else
#define YYSTACK_ALLOC malloc
#endif

/* Note: there must be only one dollar sign in this file.
   It is replaced by the list of actions, each action
   as one case of the switch.  */

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0
#define YYACCEPT	goto yyacceptlab
#define YYABORT 	goto yyabortlab
#define YYERROR		goto yyerrlab1
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL		goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(token, value) \
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    { yychar = (token), yylval = (value);			\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { yyerror ("syntax error: cannot back up"); YYERROR; }	\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

#ifndef YYPURE
#define YYLEX		yylex()
#endif

#ifdef YYPURE
#ifdef YYLSP_NEEDED
#ifdef YYLEX_PARAM
#define YYLEX		yylex(&yylval, &yylloc, YYLEX_PARAM)
#else
#define YYLEX		yylex(&yylval, &yylloc)
#endif
#else /* not YYLSP_NEEDED */
#ifdef YYLEX_PARAM
#define YYLEX		yylex(&yylval, YYLEX_PARAM)
#else
#define YYLEX		yylex(&yylval)
#endif
#endif /* not YYLSP_NEEDED */
#endif

/* If nonreentrant, generate the variables here */

#ifndef YYPURE

int	yychar;			/*  the lookahead symbol		*/
YYSTYPE	yylval;			/*  the semantic value of the		*/
				/*  lookahead symbol			*/

#ifdef YYLSP_NEEDED
YYLTYPE yylloc;			/*  location data for the lookahead	*/
				/*  symbol				*/
#endif

int yynerrs;			/*  number of parse errors so far       */
#endif  /* not YYPURE */

#if YYDEBUG != 0
int yydebug;			/*  nonzero means print parse trace	*/
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif

/*  YYINITDEPTH indicates the initial size of the parser's stacks	*/

#ifndef	YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif

/* Define __yy_memcpy.  Note that the size argument
   should be passed with type unsigned int, because that is what the non-GCC
   definitions require.  With GCC, __builtin_memcpy takes an arg
   of type size_t, but it can handle unsigned int.  */

#if __GNUC__ > 1		/* GNU C and GNU C++ define this.  */
#define __yy_memcpy(TO,FROM,COUNT)	__builtin_memcpy(TO,FROM,COUNT)
#else				/* not GNU C or C++ */
#ifndef __cplusplus

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_memcpy (to, from, count)
     char *to;
     char *from;
     unsigned int count;
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#else /* __cplusplus */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_memcpy (char *to, char *from, unsigned int count)
{
  register char *t = to;
  register char *f = from;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#endif
#endif

#line 217 "/usr/share/bison.simple"

/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
#ifdef __cplusplus
#define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else /* not __cplusplus */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
#endif /* not __cplusplus */
#else /* not YYPARSE_PARAM */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif /* not YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
#ifdef YYPARSE_PARAM
int yyparse (void *);
#else
int yyparse (void);
#endif
#endif

int
yyparse(YYPARSE_PARAM_ARG)
     YYPARSE_PARAM_DECL
{
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YYSTYPE *yyvsp;
  int yyerrstatus;	/*  number of tokens to shift before error messages enabled */
  int yychar1 = 0;		/*  lookahead token as an internal (translated) token number */

  short	yyssa[YYINITDEPTH];	/*  the state stack			*/
  YYSTYPE yyvsa[YYINITDEPTH];	/*  the semantic value stack		*/

  short *yyss = yyssa;		/*  refer to the stacks thru separate pointers */
  YYSTYPE *yyvs = yyvsa;	/*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YYLSP_NEEDED
  YYLTYPE yylsa[YYINITDEPTH];	/*  the location stack			*/
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;

#define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
#define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  int yystacksize = YYINITDEPTH;
  int yyfree_stacks = 0;

#ifdef YYPURE
  int yychar;
  YYSTYPE yylval;
  int yynerrs;
#ifdef YYLSP_NEEDED
  YYLTYPE yylloc;
#endif
#endif

  YYSTYPE yyval;		/*  the variable used to return		*/
				/*  semantic values from the action	*/
				/*  routines				*/

  int yylen;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Starting parse\n");
#endif

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YYLSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state, which is found in  yystate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
yynewstate:

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      YYSTYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;
#ifdef YYLSP_NEEDED
      YYLTYPE *yyls1 = yyls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = yyssp - yyss + 1;

#ifdef yyoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
#ifdef YYLSP_NEEDED
      /* This used to be a conditional around just the two extra args,
	 but that might be undefined if yyoverflow is a macro.  */
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yyls1, size * sizeof (*yylsp),
		 &yystacksize);
#else
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yystacksize);
#endif

      yyss = yyss1; yyvs = yyvs1;
#ifdef YYLSP_NEEDED
      yyls = yyls1;
#endif
#else /* no yyoverflow */
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	{
	  yyerror("parser stack overflow");
	  if (yyfree_stacks)
	    {
	      free (yyss);
	      free (yyvs);
#ifdef YYLSP_NEEDED
	      free (yyls);
#endif
	    }
	  return 2;
	}
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;
#ifndef YYSTACK_USE_ALLOCA
      yyfree_stacks = 1;
#endif
      yyss = (short *) YYSTACK_ALLOC (yystacksize * sizeof (*yyssp));
      __yy_memcpy ((char *)yyss, (char *)yyss1,
		   size * (unsigned int) sizeof (*yyssp));
      yyvs = (YYSTYPE *) YYSTACK_ALLOC (yystacksize * sizeof (*yyvsp));
      __yy_memcpy ((char *)yyvs, (char *)yyvs1,
		   size * (unsigned int) sizeof (*yyvsp));
#ifdef YYLSP_NEEDED
      yyls = (YYLTYPE *) YYSTACK_ALLOC (yystacksize * sizeof (*yylsp));
      __yy_memcpy ((char *)yyls, (char *)yyls1,
		   size * (unsigned int) sizeof (*yylsp));
#endif
#endif /* no yyoverflow */

      yyssp = yyss + size - 1;
      yyvsp = yyvs + size - 1;
#ifdef YYLSP_NEEDED
      yylsp = yyls + size - 1;
#endif

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

  goto yybackup;
 yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY)
    {
#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Reading a token: ");
#endif
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0)		/* This means end of input. */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more */

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      yychar1 = YYTRANSLATE(yychar);

#if YYDEBUG != 0
      if (yydebug)
	{
	  fprintf (stderr, "Next token is %d (%s", yychar, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise meaning
	     of a token, for further debugging info.  */
#ifdef YYPRINT
	  YYPRINT (stderr, yychar, yylval);
#endif
	  fprintf (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    goto yydefault;

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrlab;

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting token %d (%s), ", yychar, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (yyerrstatus) yyerrstatus--;

  yystate = yyn;
  goto yynewstate;

/* Do the default action for the current state.  */
yydefault:

  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
yyreduce:
  yylen = yyr2[yyn];
  if (yylen > 0)
    yyval = yyvsp[1-yylen]; /* implement default value of the action */

#if YYDEBUG != 0
  if (yydebug)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
	fprintf (stderr, "%s ", yytname[yyrhs[i]]);
      fprintf (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif


  switch (yyn) {

case 1:
#line 203 "lg.y"
{
	            
                        size_t sizestack = currentblock->size()+1024 ; //  before close 
                        yyvsp[-1].cinst+=currentblock->close(currentblock);
                        cout << " sizestack + 1024 =" << sizestack << "  ( " << sizestack-1024 <<" )\n" ;   
                        size_t lg0,lg1;                       
                        int NbPtr = ShowAlloc("init execution ",lg0); // number of un delele ptr
                        cout << endl;  
                        { Stack stack = newStack(sizestack);
                        double CPUcompile= CPUtime();
                        try {                  
                          yyvsp[-1].cinst.eval(stack);}
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

                        cout << "times: compile "<< CPUcompile-CPUcompileInit <<"s, execution " <<  CPUtime()-CPUcompile << "s\n";
                        deleteStack(stack);
                        //debugstack.clear() 
                        } 
                        fingraphique();
                        NbPtr = ShowAlloc("end execution -- ",lg1) - NbPtr;
                        
                        if (NbPtr) { cout << " ######## We forget of deleting   " << NbPtr << " Nb pointer,   " <<  lg1-lg0 << "Bytes\n" ;}
  return 0;;
    break;}
case 3:
#line 242 "lg.y"
{yyval.cinst=yyvsp[0].cexp;;;;
    break;}
case 4:
#line 243 "lg.y"
{ yyval.cinst= (yyvsp[-1].cinst+=yyvsp[0].cexp) ;
    break;}
case 5:
#line 246 "lg.y"
{ yyval.clist_id=new ListOfId();;
    break;}
case 6:
#line 247 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str));
    break;}
case 7:
#line 248 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 8:
#line 249 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>()));
    break;}
case 9:
#line 250 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true));
    break;}
case 10:
#line 251 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 11:
#line 252 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 12:
#line 253 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 13:
#line 254 "lg.y"
{ yyval.clist_id = yyvsp[-2].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str)) ;
    break;}
case 14:
#line 255 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 15:
#line 256 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 16:
#line 257 "lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>())) ;
    break;}
case 17:
#line 258 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true)) ;
    break;}
case 18:
#line 259 "lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 19:
#line 260 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 20:
#line 263 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 21:
#line 264 "lg.y"
{ yyval.clist_id=yyvsp[-2].clist_id  ; yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 24:
#line 269 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[0].str,dcltype);
    break;}
case 25:
#line 270 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-2].str,dcltype,yyvsp[0].cexp);
    break;}
case 26:
#line 271 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-3].str,dcltype,yyvsp[-1].args);
                                              yyvsp[-1].args.destroy();
    break;}
case 27:
#line 273 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 28:
#line 280 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 29:
#line 281 "lg.y"
{yyval.args=Find(yyvsp[-1].str);
    break;}
case 30:
#line 282 "lg.y"
{ yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);
    break;}
case 31:
#line 283 "lg.y"
{ yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 32:
#line 284 "lg.y"
{ yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp));
    break;}
case 34:
#line 288 "lg.y"
{yyval.type=TypeArray(yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 35:
#line 289 "lg.y"
{yyval.type=TypeArray(yyvsp[-5].type,yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 36:
#line 290 "lg.y"
{yyval.type=TypeTemplate(yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 37:
#line 297 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[0].str,currentblock,fespacetype,fespacecomplex); ;
    break;}
case 38:
#line 298 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;
    break;}
case 39:
#line 299 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-2].str,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;
    break;}
case 40:
#line 300 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-1].clist_id,currentblock,fespacetype,fespacecomplex) ;
    break;}
case 41:
#line 301 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;
    break;}
case 42:
#line 302 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-3].clist_id,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;
    break;}
case 43:
#line 305 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;
    break;}
case 44:
#line 306 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;
    break;}
case 45:
#line 310 "lg.y"
{fespacecomplex=false;  fespacetype = Find(yyvsp[0].str);;
    break;}
case 46:
#line 311 "lg.y"
{
             if (yyvsp[-1].type != typevarreal && yyvsp[-1].type != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=(yyvsp[-1].type==typevarcomplex);
             fespacetype = Find(yyvsp[-3].str);;
    break;}
case 47:
#line 316 "lg.y"
{  yyval.cexp = yyvsp[0].cexp  ;
    break;}
case 48:
#line 317 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 49:
#line 319 "lg.y"
{  yyval.cexp = yyvsp[0].cexp  ;
    break;}
case 50:
#line 320 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 51:
#line 322 "lg.y"
{ yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;
    break;}
case 52:
#line 323 "lg.y"
{ yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;
    break;}
case 53:
#line 328 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariableFES,size_t>(yyvsp[-3].str,atype<pfes*>(),yyvsp[-1].args,dimFESpaceImage(yyvsp[-1].args));
     yyvsp[-1].args.destroy(); ;
    break;}
case 55:
#line 332 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 56:
#line 335 "lg.y"
{dcltype=yyvsp[0].type;
    break;}
case 57:
#line 335 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 58:
#line 336 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 59:
#line 337 "lg.y"
{ yyval.cexp=yyvsp[-1].cexp;
    break;}
case 60:
#line 338 "lg.y"
{yyval.cexp=currentblock->NewID(yyvsp[-4].type,yyvsp[-3].str,yyvsp[-1].cexp);;
    break;}
case 61:
#line 340 "lg.y"
{   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = yyvsp[-4].type->right();
                      routineinblock[kkembtype] = currentblock;
                      yyvsp[-1].routine=new Routine(yyvsp[-5].type,yyvsp[-4].type->right(),yyvsp[-3].str,yyvsp[-1].clist_id,currentblock);
                     // cout << " \n after new routine \n " << endl;                      
                      ;
    break;}
case 62:
#line 348 "lg.y"
{ currentblock=yyvsp[-5].routine->Set(yyvsp[-1].cinst);
                       currentblock->Add(yyvsp[-7].str,"(",yyvsp[-5].routine);
                       kkembtype--;
                       yyval.cexp=0;
                    
                        ;
    break;}
case 63:
#line 355 "lg.y"
{Block::open(currentblock); yyvsp[-4].type->SetArgs(yyvsp[-1].clist_id);;
    break;}
case 64:
#line 357 "lg.y"
{  yyval.cinst=currentblock->close(currentblock);
                         yyval.cexp=currentblock->NewID(yyvsp[-8].type,yyvsp[-7].str,yyvsp[-1].cexp,*yyvsp[-5].clist_id);
                         delete yyvsp[-5].clist_id; //  FH 23032005
                         ;
    break;}
case 65:
#line 363 "lg.y"
{  Block::open(currentblock);
    break;}
case 66:
#line 364 "lg.y"
{  yyval.cexp=currentblock->close(currentblock);
    break;}
case 67:
#line 366 "lg.y"
{ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;
    break;}
case 68:
#line 368 "lg.y"
{ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;
    break;}
case 69:
#line 373 "lg.y"
{dcltype=yyvsp[0].type; Block::open(currentblock);  ;
    break;}
case 70:
#line 374 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 71:
#line 376 "lg.y"
{ Block::open(currentblock) ;
    break;}
case 72:
#line 378 "lg.y"
{yyval.cexp=0;;
    break;}
case 73:
#line 379 "lg.y"
{zzzfff->input(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 74:
#line 380 "lg.y"
{load(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 75:
#line 381 "lg.y"
{yyval.cexp=Try(yyvsp[-2].cinst,yyvsp[0].cexp,currentblock->close(currentblock));;
    break;}
case 76:
#line 382 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 77:
#line 383 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 78:
#line 384 "lg.y"
{inloopcount--; yyval.cexp=For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 79:
#line 386 "lg.y"
{inloopcount--; 
                yyval.cexp=C_F0(For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp),currentblock->close(currentblock));
    break;}
case 80:
#line 389 "lg.y"
{inloopcount--;yyval.cexp=While(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 81:
#line 390 "lg.y"
{yyval.cexp=FIf(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 82:
#line 391 "lg.y"
{yyval.cexp=FIf(yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 83:
#line 392 "lg.y"
{ 
                      yyval.cexp=C_F0(new E_block(yyvsp[-1].cinst,yyvsp[0].cexp),atype<void>()) ;
    break;}
case 84:
#line 394 "lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-1].str,C_F0(TheOperators,"[border]",yyvsp[0].args));
    break;}
case 85:
#line 396 "lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-4].str,C_F0(TheOperators,"[border]",yyvsp[-2].args));
    break;}
case 86:
#line 399 "lg.y"
{
                    if(inloopcount) 
                      yyval.cexp= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;
    break;}
case 87:
#line 403 "lg.y"
{ 
                    if(inloopcount)
                        yyval.cexp= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");
    break;}
case 88:
#line 407 "lg.y"
{ 
                    if (kkembtype>=0)
                      yyval.cexp= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo(yyvsp[-1].cexp)) ,atype<void>());
                     else lgerror(" return not in routine ") ;
    break;}
case 89:
#line 414 "lg.y"
{yyval.cexp =  yyvsp[0].cexp; ;
    break;}
case 90:
#line 417 "lg.y"
{ 
   Block::open(currentblock);
   yyval.args = currentblock->NewVar<LocalVariable>(yyvsp[-5].str,atype<double*>());
   yyval.args+= yyvsp[-3].cexp;
   yyval.args+= yyvsp[-1].cexp ;
    break;}
case 91:
#line 423 "lg.y"
{   
   yyval.args = (yyvsp[-1].args += yyvsp[0].cexp);
   currentblock->close(currentblock);
    break;}
case 93:
#line 430 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 100:
#line 444 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 101:
#line 445 "lg.y"
{yyval.cexp=C_F0(TheOperators,"+=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 102:
#line 446 "lg.y"
{yyval.cexp=C_F0(TheOperators,"-=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 103:
#line 447 "lg.y"
{yyval.cexp=C_F0(TheOperators,"*=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 104:
#line 448 "lg.y"
{yyval.cexp=C_F0(TheOperators,"/=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 106:
#line 454 "lg.y"
{yyval.cexp=C_F0(TheOperators,"?:",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 108:
#line 458 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 109:
#line 459 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 110:
#line 460 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 111:
#line 461 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 112:
#line 462 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 113:
#line 463 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 114:
#line 464 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 115:
#line 465 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 116:
#line 466 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 117:
#line 467 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 118:
#line 468 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 119:
#line 469 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 120:
#line 470 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 121:
#line 471 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 122:
#line 472 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 123:
#line 473 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 124:
#line 474 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 125:
#line 475 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 126:
#line 476 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 127:
#line 481 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 128:
#line 482 "lg.y"
{yyval.cexp=C_F0(TheOperators,":");
    break;}
case 129:
#line 483 "lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 130:
#line 484 "lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 131:
#line 487 "lg.y"
{yyval.args=0;
    break;}
case 132:
#line 488 "lg.y"
{yyval.args=Find(yyvsp[0].str);
    break;}
case 133:
#line 489 "lg.y"
{ yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);
    break;}
case 134:
#line 490 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 135:
#line 491 "lg.y"
{ yyval.args = (yyvsp[-2].args += Find(yyvsp[0].str)) ;
    break;}
case 136:
#line 492 "lg.y"
{ yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 137:
#line 493 "lg.y"
{ yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 138:
#line 496 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 139:
#line 497 "lg.y"
{yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 141:
#line 502 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[0].cexp);
    break;}
case 143:
#line 506 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 144:
#line 507 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 145:
#line 508 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 146:
#line 512 "lg.y"
{yyval.cexp=Find(yyvsp[0].str);;
    break;}
case 147:
#line 513 "lg.y"
{yyval.cexp= CConstant(yyvsp[0].lnum);
    break;}
case 148:
#line 514 "lg.y"
{yyval.cexp= CConstant(yyvsp[0].dnum);
    break;}
case 149:
#line 515 "lg.y"
{yyval.cexp= CConstant(complex<double>(0,yyvsp[0].dnum));
    break;}
case 150:
#line 516 "lg.y"
{yyval.cexp= CConstant<const char *>(yyvsp[0].str);
    break;}
case 151:
#line 517 "lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].args);;
    break;}
case 152:
#line 518 "lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].cexp);
    break;}
case 153:
#line 519 "lg.y"
{yyval.cexp=C_F0(yyvsp[-5].cexp,yyvsp[-4].oper,yyvsp[-3].cexp,yyvsp[-1].cexp);
    break;}
case 154:
#line 520 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,"[]");
    break;}
case 155:
#line 521 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].str) ;;
    break;}
case 156:
#line 522 "lg.y"
{ yyval.cexp=C_F0(Find(yyvsp[-2].str),yyvsp[0].str) ;;
    break;}
case 157:
#line 523 "lg.y"
{ yyval.cexp=C_F0(Find(yyvsp[-3].str),yyvsp[-2].oper,yyvsp[-1].args) ;;
    break;}
case 158:
#line 524 "lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 159:
#line 525 "lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 160:
#line 526 "lg.y"
{
             if (yyvsp[-3].type->right()->CastingFrom(yyvsp[-1].cexp.left()) ) 
                yyval.cexp=yyvsp[-3].type->right()->CastTo(yyvsp[-1].cexp)  ;
             else { yyval.cexp=yyvsp[-3].type->right()->Find("<--",basicAC_F0_wa(yyvsp[-1].cexp));
             if (!yyval.cexp.left()) { cerr << " no wait to change " << yyvsp[-1].cexp.left()->right()->name() << " in " << 
                                        yyvsp[-3].type->right()->name() << endl;
                                CompileError(" Error in type(exp) "); }
             }
            ;
    break;}
case 161:
#line 535 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 162:
#line 536 "lg.y"
{ yyval.cexp=C_F0(TheOperators,"[]",yyvsp[-1].args);
    break;}
}
   /* the action file gets copied in in place of this dollarsign */
#line 543 "/usr/share/bison.simple"

  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;

#ifdef YYLSP_NEEDED
  yylsp++;
  if (yylen == 0)
    {
      yylsp->first_line = yylloc.first_line;
      yylsp->first_column = yylloc.first_column;
      yylsp->last_line = (yylsp-1)->last_line;
      yylsp->last_column = (yylsp-1)->last_column;
      yylsp->text = 0;
    }
  else
    {
      yylsp->last_line = (yylsp+yylen-1)->last_line;
      yylsp->last_column = (yylsp+yylen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;

yyerrlab:   /* here on detecting error */

  if (! yyerrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++yynerrs;

#ifdef YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  /* Start X at -yyn if nec to avoid negative indexes in yycheck.  */
	  for (x = (yyn < 0 ? -yyn : 0);
	       x < (sizeof(yytname) / sizeof(char *)); x++)
	    if (yycheck[x + yyn] == x)
	      size += strlen(yytname[x]) + 15, count++;
	  msg = (char *) malloc(size + 15);
	  if (msg != 0)
	    {
	      strcpy(msg, "parse error");

	      if (count < 5)
		{
		  count = 0;
		  for (x = (yyn < 0 ? -yyn : 0);
		       x < (sizeof(yytname) / sizeof(char *)); x++)
		    if (yycheck[x + yyn] == x)
		      {
			strcat(msg, count == 0 ? ", expecting `" : " or `");
			strcat(msg, yytname[x]);
			strcat(msg, "'");
			count++;
		      }
		}
	      yyerror(msg);
	      free(msg);
	    }
	  else
	    yyerror ("parse error; also virtual memory exceeded");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror("parse error");
    }

  goto yyerrlab1;
yyerrlab1:   /* here on error raised explicitly by an action */

  if (yyerrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (yychar == YYEOF)
	YYABORT;

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Discarding token %d (%s).\n", yychar, yytname[yychar1]);
#endif

      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3;		/* Each real token shifted decrements this */

  goto yyerrhandle;

yyerrdefault:  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (yyn) goto yydefault;
#endif

yyerrpop:   /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss) YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

yyerrhandle:

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yyerrdefault;

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    goto yyerrdefault;

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrpop;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrpop;

  if (yyn == YYFINAL)
    YYACCEPT;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;

 yyacceptlab:
  /* YYACCEPT comes here.  */
  if (yyfree_stacks)
    {
      free (yyss);
      free (yyvs);
#ifdef YYLSP_NEEDED
      free (yyls);
#endif
    }
  return 0;

 yyabortlab:
  /* YYABORT comes here.  */
  if (yyfree_stacks)
    {
      free (yyss);
      free (yyvs);
#ifdef YYLSP_NEEDED
      free (yyls);
#endif
    }
  return 1;
}
#line 541 "lg.y"
 


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
void init_algo();
bool withrgraphique = false;
//string  StrVersionNumber();

int Compile()
{
  extern   YYSTYPE *plglval;  // modif FH 
  plglval = &lglval;
  int retvalue=0;
  int ok;
  
  currentblock=0;
  Block::open(currentblock);  
  try {
    retvalue=yyparse(); //  compile
    if(retvalue==0)  
      if(currentblock) 
	retvalue=1,cerr <<  "Error:a block is not close" << endl;   
      else {
        cerr << " CodeAlloc : nb ptr  "<< CodeAlloc::nb << ",  size :"  <<  CodeAlloc::lg << endl;
	cerr <<  "Bien: On a fini Normalement" << endl; 
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

  size_t lg000;
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

  cout << "-- FreeFem++ v" << StrVersionNumber() << endl;
  char *  cc= new char [1024];
  istream * ccin=0;
  if ( ! getprog(cc,argc,argv)>0) 
    return 1; 
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
   cout << " Load: ";
   init_lgfem() ;
   init_lgmesh() ;
   init_algo();
   
#ifdef HAVE_LIBARPACK
   init_eigenvalue();
#endif   
#ifdef PARALLELE
   init_lgparallele(); 
#endif 
#ifdef HAVE_LIBUMFPACK   
  cout << " UMFPACK ";  
#endif
 // callInitsFunct(); Pb opimisation 
   cout << endl;
  zzzfff->input(cc);
  GetEnvironment();     
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


 
