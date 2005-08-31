
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
#define	TYPE	289
#define	FUNCTION	290
#define	FESPACE	291
#define	PLUSEQ	292
#define	MOINSEQ	293
#define	MULEQ	294
#define	DIVEQ	295
#define	ARROW	296
#define	BORDER	297
#define	CURVE	298
#define	SOLVE	299

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


#line 76 "lg.y"
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



#define	YYFINAL		347
#define	YYFLAG		-32768
#define	YYNTBASE	71

#define YYTRANSLATE(x) ((unsigned)(x) <= 299 ? yytranslate[x] : 112)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    30,     2,     2,     2,    24,    13,    32,    34,
    37,    22,    20,     5,    21,    36,    23,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,    70,    66,    16,
     6,    17,    69,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
    35,     2,    38,    31,    33,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    67,    11,    68,     2,     2,     2,     2,     2,
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
    57,    58,    59,    60,    61,    62,    63,    64,    65
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
   270,   272,   275,   278,   281,   283,   293,   303,   309,   315,
   323,   327,   331,   338,   341,   344,   348,   356,   359,   361,
   365,   367,   369,   371,   373,   375,   377,   381,   385,   389,
   393,   397,   399,   405,   407,   411,   415,   419,   423,   427,
   431,   435,   439,   443,   447,   451,   455,   459,   463,   467,
   471,   475,   479,   483,   485,   487,   491,   497,   498,   500,
   504,   506,   510,   514,   520,   522,   526,   528,   531,   533,
   537,   541,   544,   546,   548,   550,   552,   554,   559,   564,
   571,   575,   579,   583,   588,   591,   594,   599,   603
};

static const short yyrhs[] = {    72,
    46,     0,    73,     0,    98,     0,    73,    98,     0,     0,
    76,     0,    76,     6,   103,     0,    57,    76,     0,    57,
    13,    76,     0,    79,    76,     0,    79,    13,    76,     0,
    35,    74,    38,     0,    74,     5,    76,     0,    74,     5,
    35,    74,    38,     0,    74,     5,    76,     6,   103,     0,
    74,     5,    57,    76,     0,    74,     5,    57,    13,    76,
     0,    74,     5,    79,    76,     0,    74,     5,    79,    13,
    76,     0,    76,     0,    75,     5,    76,     0,    42,     0,
    57,     0,    42,     0,    42,     6,   103,     0,    42,    34,
    78,    37,     0,    77,     5,    77,     0,   105,     0,    57,
    42,     0,    42,     6,   105,     0,    78,     5,   105,     0,
    78,     5,    76,     6,   105,     0,    55,     0,    55,    35,
    55,    38,     0,    55,    35,    55,     5,    55,    38,     0,
    55,    16,    55,    17,     0,    42,     0,    42,    35,   105,
    38,     0,    42,     6,   105,     0,    35,    75,    38,     0,
    35,    75,    38,    35,   105,    38,     0,    35,    75,    38,
     6,   105,     0,    42,    34,   105,    37,     0,    35,    75,
    38,    34,   105,    37,     0,    57,     0,    57,    16,    55,
    17,     0,    81,     0,    83,     5,    81,     0,    80,     0,
    84,     5,    80,     0,    82,    84,     0,    82,    35,    55,
    38,    83,     0,    42,    34,    78,    37,     0,    86,     0,
    87,     5,    86,     0,     0,    79,    89,    77,    66,     0,
    43,    87,    66,     0,    85,    66,     0,    56,    42,     6,
   101,    66,     0,     0,    56,    79,    42,    34,    74,    37,
    90,    67,    73,    68,     0,     0,    56,    42,    34,    74,
    37,    91,     6,   103,    66,     0,    67,     0,    68,     0,
    50,     0,    51,     0,     0,    79,    97,    77,     0,    66,
     0,    47,    45,     0,    48,    45,     0,   101,    66,     0,
    88,     0,    94,    34,   101,    66,   101,    66,   101,    37,
    98,     0,    94,    34,    96,    66,   101,    66,   101,    37,
    98,     0,    95,    34,   101,    37,    98,     0,     3,    34,
   101,    37,    98,     0,     3,    34,   101,    37,    98,     4,
    98,     0,    92,    73,    93,     0,    63,    42,   100,     0,
    63,    42,    35,   108,    38,    66,     0,    52,    66,     0,
    53,    66,     0,    54,   101,    66,     0,    34,    42,     6,
   101,     5,   101,    37,     0,    99,    98,     0,   103,     0,
   101,     5,   101,     0,    21,     0,    20,     0,    30,     0,
    28,     0,    29,     0,   104,     0,   104,     6,   103,     0,
   104,    58,   103,     0,   104,    59,   103,     0,   104,    60,
   103,     0,   104,    61,   103,     0,   105,     0,   105,    69,
   104,    70,   104,     0,   109,     0,   105,    22,   105,     0,
   105,    25,   105,     0,   105,    26,   105,     0,   105,    23,
   105,     0,   105,    24,   105,     0,   105,    20,   105,     0,
   105,    21,   105,     0,   105,     8,   105,     0,   105,     9,
   105,     0,   105,    13,   105,     0,   105,    12,   105,     0,
   105,    11,   105,     0,   105,    10,   105,     0,   105,    16,
   105,     0,   105,    18,   105,     0,   105,    17,   105,     0,
   105,    19,   105,     0,   105,    14,   105,     0,   105,    15,
   105,     0,   105,     0,    70,     0,   105,    70,   105,     0,
   105,    70,   105,    70,   105,     0,     0,    57,     0,    76,
     6,   105,     0,   106,     0,   107,     5,    57,     0,   107,
     5,   106,     0,   107,     5,    76,     6,   105,     0,   103,
     0,   108,     5,   103,     0,   110,     0,   102,   110,     0,
   111,     0,   111,    31,   109,     0,   111,    33,   109,     0,
   111,    32,     0,    42,     0,    39,     0,    40,     0,    41,
     0,    45,     0,   111,    34,   107,    37,     0,   111,    35,
   106,    38,     0,   111,    35,   106,     5,   106,    38,     0,
   111,    35,    38,     0,   111,    36,    42,     0,    57,    36,
    42,     0,    57,    34,   107,    37,     0,   111,    28,     0,
   111,    29,     0,    55,    34,   101,    37,     0,    34,   101,
    37,     0,    35,   108,    38,     0
};

#endif

#if YYDEBUG != 0
static const short yyrline[] = { 0,
   195,   231,   234,   235,   238,   239,   240,   241,   242,   243,
   244,   245,   246,   247,   248,   249,   250,   251,   252,   255,
   256,   259,   259,   261,   262,   263,   265,   271,   273,   274,
   275,   276,   279,   280,   281,   282,   288,   290,   291,   292,
   293,   294,   296,   298,   302,   303,   308,   309,   311,   312,
   314,   315,   318,   323,   324,   327,   327,   328,   329,   330,
   331,   339,   346,   348,   355,   356,   358,   360,   364,   366,
   368,   369,   370,   371,   372,   373,   374,   378,   379,   380,
   381,   383,   385,   388,   392,   396,   404,   410,   415,   417,
   421,   423,   424,   425,   426,   429,   431,   432,   433,   434,
   435,   438,   440,   442,   444,   445,   446,   447,   448,   449,
   450,   451,   452,   453,   454,   455,   456,   457,   458,   459,
   460,   461,   462,   466,   468,   469,   470,   473,   474,   475,
   476,   477,   478,   479,   482,   483,   486,   488,   491,   492,
   493,   494,   497,   499,   500,   501,   502,   503,   504,   505,
   506,   507,   508,   509,   510,   511,   512,   521,   522
};
#endif


#if YYDEBUG != 0 || defined (YYERROR_VERBOSE)

static const char * const yytname[] = {   "$","error","$undefined.","IF","ELSE",
"','","'='","SET","LTLT","GTGT","OR","'|'","AND","'&'","EQ","NE","'<'","'>'",
"LE","GE","'+'","'-'","'*'","'/'","'%'","DOTSTAR","DOTSLASH","UNARY","PLUSPLUS",
"MOINSMOINS","'!'","'^'","'\\''","'_'","'('","'['","'.'","')'","']'","LNUM",
"DNUM","CNUM","ID","FESPACEID","IDPARAM","STRING","ENDOFFILE","INCLUDE","LOAD",
"BIDON","FOR","WHILE","BREAK","CONTINUE","RETURN","TYPE","FUNCTION","FESPACE",
"PLUSEQ","MOINSEQ","MULEQ","DIVEQ","ARROW","BORDER","CURVE","SOLVE","';'","'{'",
"'}'","'?'","':'","start","input","instructions","list_of_id_args","list_of_id1",
"id","list_of_dcls","parameters_list","type_of_dcl","ID_space","ID_array_space",
"fespace","spaceIDa","spaceIDb","spaceIDs","fespace_def","fespace_def_list",
"declaration","@1","@2","@3","begin","end","for_loop","while_loop","declaration_for",
"@4","instruction","bornes","border_expr","Expr","unop","no_comma_expr","no_ternary_expr",
"no_set_expr","sub_script_expr","parameters","array","unary_expr","pow_expr",
"primary", NULL
};
#endif

static const short yyr1[] = {     0,
    71,    72,    73,    73,    74,    74,    74,    74,    74,    74,
    74,    74,    74,    74,    74,    74,    74,    74,    74,    75,
    75,    76,    76,    77,    77,    77,    77,    78,    78,    78,
    78,    78,    79,    79,    79,    79,    80,    80,    80,    80,
    80,    80,    81,    81,    82,    82,    83,    83,    84,    84,
    85,    85,    86,    87,    87,    89,    88,    88,    88,    88,
    90,    88,    91,    88,    92,    93,    94,    95,    97,    96,
    98,    98,    98,    98,    98,    98,    98,    98,    98,    98,
    98,    98,    98,    98,    98,    98,    99,   100,   101,   101,
   102,   102,   102,   102,   102,   103,   103,   103,   103,   103,
   103,   104,   104,   105,   105,   105,   105,   105,   105,   105,
   105,   105,   105,   105,   105,   105,   105,   105,   105,   105,
   105,   105,   105,   106,   106,   106,   106,   107,   107,   107,
   107,   107,   107,   107,   108,   108,   109,   109,   110,   110,
   110,   110,   111,   111,   111,   111,   111,   111,   111,   111,
   111,   111,   111,   111,   111,   111,   111,   111,   111
};

static const short yyr2[] = {     0,
     2,     1,     1,     2,     0,     1,     3,     2,     3,     2,
     3,     3,     3,     5,     5,     4,     5,     4,     5,     1,
     3,     1,     1,     1,     3,     4,     3,     1,     2,     3,
     3,     5,     1,     4,     6,     4,     1,     4,     3,     3,
     6,     5,     4,     6,     1,     4,     1,     3,     1,     3,
     2,     5,     4,     1,     3,     0,     4,     3,     2,     5,
     0,    10,     0,     9,     1,     1,     1,     1,     0,     3,
     1,     2,     2,     2,     1,     9,     9,     5,     5,     7,
     3,     3,     6,     2,     2,     3,     7,     2,     1,     3,
     1,     1,     1,     1,     1,     1,     3,     3,     3,     3,
     3,     1,     5,     1,     3,     3,     3,     3,     3,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     3,     3,     1,     1,     3,     5,     0,     1,     3,
     1,     3,     3,     5,     1,     3,     1,     2,     1,     3,
     3,     2,     1,     1,     1,     1,     1,     4,     4,     6,
     3,     3,     3,     4,     2,     2,     4,     3,     3
};

static const short yydefact[] = {     0,
     0,    92,    91,    94,    95,    93,     0,     0,   144,   145,
   146,   143,     0,   147,     0,     0,    67,    68,     0,     0,
     0,    33,     0,    45,     0,    71,    65,     0,     2,    56,
     0,     0,    75,     0,     0,     0,     3,     0,     0,    89,
    96,   102,   104,   137,   139,     0,     0,     0,     0,   135,
     0,     0,    54,     0,    72,    73,    84,    85,     0,     0,
     0,     0,     0,    33,     0,     0,   128,     0,     0,     1,
     4,     0,     0,    37,    49,    51,    59,     0,     0,     0,
     0,    74,   138,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,   155,   156,
     0,   142,     0,   128,     0,     0,     0,   158,     0,   159,
     0,     0,    58,    86,     0,     0,     0,     0,     5,     0,
     0,   143,   129,   125,     0,   124,   131,     0,   153,     0,
     0,     0,    82,    24,     0,    22,     0,    23,     0,    20,
     0,     0,     0,    66,    81,    69,     0,     0,     0,    90,
    97,    98,    99,   100,   101,   112,   113,   117,   116,   115,
   114,   122,   123,   118,   120,   119,   121,   110,   111,   105,
   108,   109,   106,   107,     0,   140,   141,     0,   151,     0,
   152,     0,   136,   143,     0,     0,    28,    55,    36,   157,
     0,    34,     0,     5,    23,     0,     6,     0,     5,    46,
     0,     0,     0,   154,     0,     0,    88,     0,     0,     0,
    57,     0,     0,    40,    39,     0,     0,    50,     0,     0,
     0,     0,     0,   148,     0,   149,    79,     0,    29,     0,
    53,     0,    60,     0,     0,     8,     0,    63,     0,     0,
    10,     0,   130,   126,   132,     0,   133,     0,     0,    25,
     0,    27,     0,     0,    47,    52,    21,     0,     0,    38,
    70,     0,     0,    78,   103,     0,     0,    30,    23,     0,
    31,    35,    12,     9,     5,    23,    13,     0,     0,     7,
    11,    61,     0,     0,     0,    83,    26,     0,     0,     0,
    42,     0,     0,     0,   150,    80,     0,     0,     0,    16,
     0,     0,    18,     0,     0,   127,   134,     0,     0,     0,
    48,    41,     0,     0,    32,    14,    17,    15,    19,     0,
     0,    90,     0,    43,     0,     0,    64,     0,    87,     0,
    77,    76,    62,    44,     0,     0,     0
};

static const short yydefgoto[] = {   345,
    28,    29,   206,   149,   207,   145,   196,    30,    75,   265,
    31,   266,    76,    32,    53,    54,    33,    72,   315,   289,
    34,   155,    35,    36,   157,   229,    37,   142,   143,    38,
    39,    40,    41,    42,   137,   138,    51,    43,    44,    45
};

static const short yypact[] = {   436,
    11,-32768,-32768,-32768,-32768,-32768,   605,   605,-32768,-32768,
-32768,-32768,    47,-32768,    -5,    91,-32768,-32768,   128,   130,
   605,   234,    69,    99,   158,-32768,-32768,   160,   436,-32768,
    81,   141,-32768,   436,   177,   179,-32768,     3,   283,-32768,
    46,   496,-32768,-32768,   298,   605,   180,    19,    55,-32768,
     6,   191,-32768,     4,-32768,-32768,-32768,-32768,     5,   157,
   605,   162,   104,    27,   212,   201,   503,   215,   156,-32768,
-32768,   225,     9,   186,-32768,   270,-32768,   340,   635,   605,
   605,-32768,-32768,   605,   605,   605,   605,   605,   605,   605,
   605,   605,   605,   605,   605,   605,   605,   605,   605,   605,
   605,   605,   605,   605,   605,   605,   605,   605,-32768,-32768,
   605,-32768,   605,   503,   140,   235,    85,-32768,   605,-32768,
   665,    47,-32768,-32768,   261,   122,    41,   605,    77,   253,
   271,   286,   165,-32768,   288,     7,-32768,   125,-32768,   260,
   605,   436,-32768,   192,    29,-32768,   257,-32768,    43,-32768,
   605,   605,   213,-32768,-32768,-32768,   238,    30,   126,-32768,
-32768,-32768,-32768,-32768,-32768,   884,   884,   899,   899,   912,
   912,   756,   756,   331,   331,   331,   331,   259,   259,-32768,
-32768,-32768,-32768,-32768,   231,-32768,-32768,   168,-32768,    56,
-32768,   436,-32768,   299,   237,   171,   867,-32768,-32768,-32768,
   251,-32768,    31,    77,     0,   172,   307,    34,    77,-32768,
   605,   605,   541,-32768,   308,    62,-32768,   605,   665,   225,
-32768,   218,    -1,   187,   867,   745,    -1,-32768,   225,   605,
   605,   436,   605,-32768,   573,-32768,   311,   605,-32768,   695,
-32768,   297,-32768,    63,    -1,-32768,   217,-32768,   605,    -1,
-32768,   178,   867,   220,   165,   315,-32768,   605,   276,-32768,
   181,-32768,    -1,   302,-32768,   332,-32768,   605,   605,-32768,
   341,    32,    33,-32768,-32768,   309,   436,   867,    19,   342,
   867,-32768,-32768,-32768,    77,    36,   343,    45,   353,-32768,
-32768,-32768,   605,   605,   358,-32768,-32768,    70,   605,   218,
   867,   776,   605,   605,-32768,-32768,   605,    75,    -1,-32768,
   605,    -1,-32768,   605,   300,   867,   867,   605,   330,   807,
-32768,-32768,   182,   183,   867,-32768,-32768,-32768,-32768,   305,
   436,   329,   605,-32768,   436,   436,-32768,   395,-32768,   837,
-32768,-32768,-32768,-32768,   372,   373,-32768
};

static const short yypgoto[] = {-32768,
-32768,   -32,  -197,   113,    53,   -92,   159,   -20,   224,    86,
-32768,-32768,-32768,-32768,   267,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,   -28,-32768,-32768,    -7,
-32768,    -2,  -104,    51,  -110,   285,   263,   -48,   361,-32768
};


#define	YYLAST		938


static const short yytable[] = {    49,
    71,    78,    65,   185,   190,    50,   244,    81,   122,    81,
   119,   252,   245,    59,    89,    90,    91,    92,    93,    94,
    95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
   105,   106,   107,   220,    81,    81,    81,    81,   117,    55,
   146,   146,    60,   120,    46,   201,   250,   223,   309,    71,
   146,    84,    67,   126,    68,   148,   148,   312,   156,    81,
   235,    62,   186,   147,   187,   148,   119,   247,    82,   123,
   124,   158,   159,   160,   223,   146,   212,   146,   202,   247,
   224,   161,   162,   163,   164,   165,   146,   308,    52,    81,
   148,   118,   148,   236,   221,   231,   243,   303,   304,   259,
   283,   148,   257,    85,    86,    87,    88,   319,   208,   128,
    63,   204,   326,   217,    66,    73,   193,   136,   146,   135,
   203,   192,    74,    64,   276,   150,    81,   262,   275,   213,
    81,    64,    67,   205,    68,    56,   271,   129,    50,   166,
   167,   168,   169,   170,   171,   172,   173,   174,   175,   176,
   177,   178,   179,   180,   181,   182,   183,   184,   200,     2,
     3,   214,   232,   237,   136,   136,   135,     4,     5,     6,
   -23,   197,   213,     7,     8,   240,   247,   189,     9,    10,
    11,    12,   247,   208,    14,   240,    81,    81,   208,   140,
   141,   151,   268,    57,    47,    58,    48,   218,    67,    69,
    68,   225,   226,   274,   234,    70,    77,   241,   248,   134,
    79,   125,    80,    61,   292,   260,   127,   297,   335,   336,
   152,   269,   272,   273,   121,   219,   288,    89,    90,    91,
    92,    93,    94,    95,    96,    97,    98,    99,   100,   101,
   102,   103,   104,   105,   106,   107,   290,   227,   306,    60,
   295,   285,   263,   130,    74,   131,   139,   246,   146,   264,
   251,   253,   254,   136,   208,   256,   144,    61,    62,   197,
    67,    64,    68,   286,   153,   267,   191,   199,   239,   150,
   103,   104,   105,   106,   107,   136,   209,   210,   278,   293,
   281,   -22,   280,   211,   222,   323,   324,   284,   338,   287,
   233,   215,   291,   230,   238,   242,   341,   342,   328,    71,
   332,   330,   249,   258,   277,   150,     7,     8,   301,   302,
   294,     9,    10,    11,    12,   109,   110,    14,   111,   112,
   113,   114,   115,   116,   282,   299,   300,    47,   310,    48,
   313,   296,     1,   316,   317,   220,   305,   307,   311,   320,
   101,   102,   103,   104,   105,   106,   107,   325,   314,     2,
     3,   327,   318,   333,   329,   339,   331,     4,     5,     6,
   337,   346,   347,     7,     8,   298,   228,   261,     9,    10,
    11,    12,    13,   340,    14,   321,    15,    16,   198,    17,
    18,    19,    20,    21,    22,    23,    24,     1,   188,    83,
     0,     0,    25,   216,     0,    26,    27,   154,     0,     0,
     0,     0,     0,     0,     2,     3,     0,     0,     0,     0,
     0,     0,     4,     5,     6,     0,     0,     0,     7,     8,
     0,     0,     0,     9,    10,    11,    12,    13,     1,    14,
     0,    15,    16,     0,    17,    18,    19,    20,    21,    22,
    23,    24,     0,     0,     0,     2,     3,    25,     0,     0,
    26,    27,   343,     4,     5,     6,     0,     0,     0,     7,
     8,     0,     0,     0,     9,    10,    11,    12,    13,     0,
    14,     0,    15,    16,     0,    17,    18,    19,    20,    21,
    22,    23,    24,     0,     0,     0,     0,     0,    25,     0,
     0,    26,    27,    89,    90,    91,    92,    93,    94,    95,
    96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
   106,   107,     2,     3,     0,     0,     0,     0,     0,     0,
     4,     5,     6,     0,     0,     0,     7,     8,     0,     0,
     0,     9,    10,    11,   132,     0,     0,    14,     0,     0,
     0,     0,     0,     0,     0,     0,     0,    47,     0,   133,
     2,     3,     0,     0,   108,     0,     0,     0,     4,     5,
     6,     0,   134,     0,     7,     8,     0,     0,     0,     9,
    10,    11,   132,     0,     0,    14,     0,     0,     0,     0,
     0,     0,     2,     3,     0,    47,     0,   255,     0,     0,
     4,     5,     6,     0,     0,     0,     7,     8,     0,     0,
   134,     9,    10,    11,    12,     0,     0,    14,     0,     0,
     0,     0,     0,     0,     2,     3,     0,    47,     0,    48,
     0,     0,     4,     5,     6,     0,     0,     0,     7,     8,
     0,     0,   134,     9,    10,    11,    12,     0,     0,    14,
     0,     0,     0,     0,     2,     3,     0,     0,     0,    47,
     0,    48,     4,     5,     6,     0,     0,     0,     7,     8,
     0,     0,     0,     9,    10,    11,    12,     0,     0,    14,
     0,     0,     0,     0,     2,     3,     0,     0,     0,    22,
     0,    48,     4,     5,     6,     0,     0,     0,     7,     8,
     0,     0,     0,     9,    10,    11,   194,     0,     0,    14,
     0,     0,     0,     0,     2,     3,     0,     0,     0,    47,
     0,   195,     4,     5,     6,     0,     0,     0,     7,     8,
     0,     0,     0,     9,    10,    11,   132,     0,     0,    14,
     0,     0,     0,     0,     0,     0,     0,     0,     0,    47,
     0,   279,    89,    90,    91,    92,    93,    94,    95,    96,
    97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
   107,    97,    98,    99,   100,   101,   102,   103,   104,   105,
   106,   107,   270,    89,    90,    91,    92,    93,    94,    95,
    96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
   106,   107,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,   322,    89,    90,    91,    92,    93,    94,
    95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
   105,   106,   107,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,   334,    89,    90,    91,    92,    93,    94,
    95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
   105,   106,   107,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,   344,    89,    90,    91,    92,    93,    94,
    95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
   105,   106,   107,    91,    92,    93,    94,    95,    96,    97,
    98,    99,   100,   101,   102,   103,   104,   105,   106,   107,
    93,    94,    95,    96,    97,    98,    99,   100,   101,   102,
   103,   104,   105,   106,   107,    95,    96,    97,    98,    99,
   100,   101,   102,   103,   104,   105,   106,   107
};

static const short yycheck[] = {     7,
    29,    34,    23,   108,   115,     8,   204,     5,     5,     5,
     5,   209,    13,    21,     8,     9,    10,    11,    12,    13,
    14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
    24,    25,    26,     5,     5,     5,     5,     5,    46,    45,
    42,    42,    16,    38,    34,     5,    13,     5,    13,    78,
    42,     6,    34,    61,    36,    57,    57,    13,    79,     5,
     5,    35,   111,    55,   113,    57,     5,     5,    66,    66,
    66,    79,    80,    81,     5,    42,    70,    42,    38,     5,
    38,    84,    85,    86,    87,    88,    42,   285,    42,     5,
    57,    37,    57,    38,    66,    66,    66,    66,    66,    38,
    38,    57,   213,    58,    59,    60,    61,    38,   129,     6,
    42,    35,    38,   142,    16,    35,   119,    67,    42,    67,
   128,    37,    42,    55,   235,    73,     5,   220,   233,     5,
     5,    55,    34,    57,    36,    45,   229,    34,   141,    89,
    90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
   100,   101,   102,   103,   104,   105,   106,   107,    37,    20,
    21,    37,    37,   192,   114,   115,   114,    28,    29,    30,
     6,   121,     5,    34,    35,     5,     5,    38,    39,    40,
    41,    42,     5,   204,    45,     5,     5,     5,   209,    34,
    35,     6,     6,    66,    55,    66,    57,     6,    34,    42,
    36,   151,   152,   232,    37,    46,    66,    37,    37,    70,
    34,    55,    34,    34,    37,   218,    55,    37,    37,    37,
    35,    35,   230,   231,    34,    34,   247,     8,     9,    10,
    11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
    21,    22,    23,    24,    25,    26,   249,    35,   277,    16,
   258,    35,    35,    42,    42,    55,    42,   205,    42,    42,
   208,   211,   212,   213,   285,   213,    42,    34,    35,   219,
    34,    55,    36,    57,     5,   223,    42,    17,    42,   227,
    22,    23,    24,    25,    26,   235,    34,    17,   238,    70,
   240,     6,   240,     6,    38,   303,   304,   245,   331,   247,
    70,    42,   250,    66,     6,    55,   335,   336,   311,   338,
   318,   314,     6,     6,     4,   263,    34,    35,   268,   269,
     6,    39,    40,    41,    42,    28,    29,    45,    31,    32,
    33,    34,    35,    36,    38,    34,     5,    55,   286,    57,
   288,    66,     3,   293,   294,     5,    38,     6,     6,   299,
    20,    21,    22,    23,    24,    25,    26,   307,     6,    20,
    21,   309,     5,    34,   312,    37,    67,    28,    29,    30,
    66,     0,     0,    34,    35,   263,   153,   219,    39,    40,
    41,    42,    43,   333,    45,   300,    47,    48,   122,    50,
    51,    52,    53,    54,    55,    56,    57,     3,   114,    39,
    -1,    -1,    63,   141,    -1,    66,    67,    68,    -1,    -1,
    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    -1,
    -1,    -1,    28,    29,    30,    -1,    -1,    -1,    34,    35,
    -1,    -1,    -1,    39,    40,    41,    42,    43,     3,    45,
    -1,    47,    48,    -1,    50,    51,    52,    53,    54,    55,
    56,    57,    -1,    -1,    -1,    20,    21,    63,    -1,    -1,
    66,    67,    68,    28,    29,    30,    -1,    -1,    -1,    34,
    35,    -1,    -1,    -1,    39,    40,    41,    42,    43,    -1,
    45,    -1,    47,    48,    -1,    50,    51,    52,    53,    54,
    55,    56,    57,    -1,    -1,    -1,    -1,    -1,    63,    -1,
    -1,    66,    67,     8,     9,    10,    11,    12,    13,    14,
    15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
    25,    26,    20,    21,    -1,    -1,    -1,    -1,    -1,    -1,
    28,    29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,
    -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,    -1,    57,
    20,    21,    -1,    -1,    69,    -1,    -1,    -1,    28,    29,
    30,    -1,    70,    -1,    34,    35,    -1,    -1,    -1,    39,
    40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,    -1,
    -1,    -1,    20,    21,    -1,    55,    -1,    57,    -1,    -1,
    28,    29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,
    70,    39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,
    -1,    -1,    -1,    -1,    20,    21,    -1,    55,    -1,    57,
    -1,    -1,    28,    29,    30,    -1,    -1,    -1,    34,    35,
    -1,    -1,    70,    39,    40,    41,    42,    -1,    -1,    45,
    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    55,
    -1,    57,    28,    29,    30,    -1,    -1,    -1,    34,    35,
    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,
    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    55,
    -1,    57,    28,    29,    30,    -1,    -1,    -1,    34,    35,
    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,
    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    55,
    -1,    57,    28,    29,    30,    -1,    -1,    -1,    34,    35,
    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    55,
    -1,    57,     8,     9,    10,    11,    12,    13,    14,    15,
    16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
    26,    16,    17,    18,    19,    20,    21,    22,    23,    24,
    25,    26,    38,     8,     9,    10,    11,    12,    13,    14,
    15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
    25,    26,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    38,     8,     9,    10,    11,    12,    13,
    14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
    24,    25,    26,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    37,     8,     9,    10,    11,    12,    13,
    14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
    24,    25,    26,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    37,     8,     9,    10,    11,    12,    13,
    14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
    24,    25,    26,    10,    11,    12,    13,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
    12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
    22,    23,    24,    25,    26,    14,    15,    16,    17,    18,
    19,    20,    21,    22,    23,    24,    25,    26
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
#line 195 "lg.y"
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
#line 234 "lg.y"
{yyval.cinst=yyvsp[0].cexp;;;;
    break;}
case 4:
#line 235 "lg.y"
{ yyval.cinst= (yyvsp[-1].cinst+=yyvsp[0].cexp) ;
    break;}
case 5:
#line 238 "lg.y"
{ yyval.clist_id=new ListOfId();;
    break;}
case 6:
#line 239 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str));
    break;}
case 7:
#line 240 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 8:
#line 241 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>()));
    break;}
case 9:
#line 242 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true));
    break;}
case 10:
#line 243 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 11:
#line 244 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 12:
#line 245 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 13:
#line 246 "lg.y"
{ yyval.clist_id = yyvsp[-2].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str)) ;
    break;}
case 14:
#line 247 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 15:
#line 248 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 16:
#line 249 "lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>())) ;
    break;}
case 17:
#line 250 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true)) ;
    break;}
case 18:
#line 251 "lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 19:
#line 252 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 20:
#line 255 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 21:
#line 256 "lg.y"
{ yyval.clist_id=yyvsp[-2].clist_id  ; yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 24:
#line 261 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[0].str,dcltype);
    break;}
case 25:
#line 262 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-2].str,dcltype,yyvsp[0].cexp);
    break;}
case 26:
#line 263 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-3].str,dcltype,yyvsp[-1].args);
                                              yyvsp[-1].args.destroy();
    break;}
case 27:
#line 265 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 28:
#line 272 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 29:
#line 273 "lg.y"
{yyval.args=Find(yyvsp[-1].str);
    break;}
case 30:
#line 274 "lg.y"
{ yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);
    break;}
case 31:
#line 275 "lg.y"
{ yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 32:
#line 276 "lg.y"
{ yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp));
    break;}
case 34:
#line 280 "lg.y"
{yyval.type=TypeArray(yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 35:
#line 281 "lg.y"
{yyval.type=TypeArray(yyvsp[-5].type,yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 36:
#line 282 "lg.y"
{yyval.type=TypeTemplate(yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 37:
#line 289 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[0].str,currentblock,fespacetype,fespacecomplex); ;
    break;}
case 38:
#line 290 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;
    break;}
case 39:
#line 291 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-2].str,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;
    break;}
case 40:
#line 292 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-1].clist_id,currentblock,fespacetype,fespacecomplex) ;
    break;}
case 41:
#line 293 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;
    break;}
case 42:
#line 294 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-3].clist_id,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;
    break;}
case 43:
#line 297 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;
    break;}
case 44:
#line 298 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;
    break;}
case 45:
#line 302 "lg.y"
{fespacecomplex=false;  fespacetype = Find(yyvsp[0].str);;
    break;}
case 46:
#line 303 "lg.y"
{
             if (yyvsp[-1].type != typevarreal && yyvsp[-1].type != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=(yyvsp[-1].type==typevarcomplex);
             fespacetype = Find(yyvsp[-3].str);;
    break;}
case 47:
#line 308 "lg.y"
{  yyval.cexp = yyvsp[0].cexp  ;
    break;}
case 48:
#line 309 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 49:
#line 311 "lg.y"
{  yyval.cexp = yyvsp[0].cexp  ;
    break;}
case 50:
#line 312 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 51:
#line 314 "lg.y"
{ yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;
    break;}
case 52:
#line 315 "lg.y"
{ yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;
    break;}
case 53:
#line 320 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariableFES,size_t>(yyvsp[-3].str,atype<pfes*>(),yyvsp[-1].args,dimFESpaceImage(yyvsp[-1].args));
     yyvsp[-1].args.destroy(); ;
    break;}
case 55:
#line 324 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 56:
#line 327 "lg.y"
{dcltype=yyvsp[0].type;
    break;}
case 57:
#line 327 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 58:
#line 328 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 59:
#line 329 "lg.y"
{ yyval.cexp=yyvsp[-1].cexp;
    break;}
case 60:
#line 330 "lg.y"
{yyval.cexp=currentblock->NewID(yyvsp[-4].type,yyvsp[-3].str,yyvsp[-1].cexp);;
    break;}
case 61:
#line 332 "lg.y"
{   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = yyvsp[-4].type->right();
                      routineinblock[kkembtype] = currentblock;
                      yyvsp[-1].routine=new Routine(yyvsp[-5].type,yyvsp[-4].type->right(),yyvsp[-3].str,yyvsp[-1].clist_id,currentblock);
                     // cout << " \n after new routine \n " << endl;                      
                      ;
    break;}
case 62:
#line 340 "lg.y"
{ currentblock=yyvsp[-5].routine->Set(yyvsp[-1].cinst);
                       currentblock->Add(yyvsp[-7].str,"(",yyvsp[-5].routine);
                       kkembtype--;
                       yyval.cexp=0;
                    
                        ;
    break;}
case 63:
#line 347 "lg.y"
{currentblock = new Block(currentblock); yyvsp[-4].type->SetArgs(yyvsp[-1].clist_id);;
    break;}
case 64:
#line 349 "lg.y"
{  yyval.cinst=currentblock->close(currentblock);
                         yyval.cexp=currentblock->NewID(yyvsp[-8].type,yyvsp[-7].str,yyvsp[-1].cexp,*yyvsp[-5].clist_id);
                         delete yyvsp[-5].clist_id; //  FH 23032005
                         ;
    break;}
case 65:
#line 355 "lg.y"
{  currentblock = new Block(currentblock);
    break;}
case 66:
#line 356 "lg.y"
{  yyval.cexp=currentblock->close(currentblock);
    break;}
case 67:
#line 358 "lg.y"
{ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;
    break;}
case 68:
#line 360 "lg.y"
{ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;
    break;}
case 69:
#line 365 "lg.y"
{dcltype=yyvsp[0].type;currentblock = new Block(currentblock);
    break;}
case 70:
#line 366 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 71:
#line 368 "lg.y"
{yyval.cexp=0;;
    break;}
case 72:
#line 369 "lg.y"
{zzzfff->input(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 73:
#line 370 "lg.y"
{load(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 74:
#line 371 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 75:
#line 372 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 76:
#line 373 "lg.y"
{inloopcount--; yyval.cexp=For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 77:
#line 375 "lg.y"
{inloopcount--; 
                yyval.cexp=C_F0(For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp),currentblock->close(currentblock));
    break;}
case 78:
#line 378 "lg.y"
{inloopcount--;yyval.cexp=While(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 79:
#line 379 "lg.y"
{yyval.cexp=FIf(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 80:
#line 380 "lg.y"
{yyval.cexp=FIf(yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 81:
#line 381 "lg.y"
{ 
                      yyval.cexp=C_F0(new E_block(yyvsp[-1].cinst,yyvsp[0].cexp),atype<void>()) ;
    break;}
case 82:
#line 383 "lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-1].str,C_F0(TheOperators,"[border]",yyvsp[0].args));
    break;}
case 83:
#line 385 "lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-4].str,C_F0(TheOperators,"[border]",yyvsp[-2].args));
    break;}
case 84:
#line 388 "lg.y"
{
                    if(inloopcount) 
                      yyval.cexp= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;
    break;}
case 85:
#line 392 "lg.y"
{ 
                    if(inloopcount)
                        yyval.cexp= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");
    break;}
case 86:
#line 396 "lg.y"
{ 
                    if (kkembtype>=0)
                      yyval.cexp= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo(yyvsp[-1].cexp)) ,atype<void>());
                     else lgerror(" return not in routine ") ;
    break;}
case 87:
#line 404 "lg.y"
{ 
   currentblock = new Block(currentblock);
   yyval.args = currentblock->NewVar<LocalVariable>(yyvsp[-5].str,atype<double*>());
   yyval.args+= yyvsp[-3].cexp;
   yyval.args+= yyvsp[-1].cexp ;
    break;}
case 88:
#line 410 "lg.y"
{   
   yyval.args = (yyvsp[-1].args += yyvsp[0].cexp);
   currentblock->close(currentblock);
    break;}
case 90:
#line 417 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 97:
#line 431 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 98:
#line 432 "lg.y"
{yyval.cexp=C_F0(TheOperators,"+=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 99:
#line 433 "lg.y"
{yyval.cexp=C_F0(TheOperators,"-=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 100:
#line 434 "lg.y"
{yyval.cexp=C_F0(TheOperators,"*=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 101:
#line 435 "lg.y"
{yyval.cexp=C_F0(TheOperators,"/=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 103:
#line 440 "lg.y"
{yyval.cexp=C_F0(TheOperators,"?:",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 105:
#line 444 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 106:
#line 445 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 107:
#line 446 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 108:
#line 447 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 109:
#line 448 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 110:
#line 449 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 111:
#line 450 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 112:
#line 451 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 113:
#line 452 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 114:
#line 453 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 115:
#line 454 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 116:
#line 455 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 117:
#line 456 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 118:
#line 457 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 119:
#line 458 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 120:
#line 459 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 121:
#line 460 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 122:
#line 461 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 123:
#line 462 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 124:
#line 467 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 125:
#line 468 "lg.y"
{yyval.cexp=C_F0(TheOperators,":");
    break;}
case 126:
#line 469 "lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 127:
#line 470 "lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 128:
#line 473 "lg.y"
{yyval.args=0;
    break;}
case 129:
#line 474 "lg.y"
{yyval.args=Find(yyvsp[0].str);
    break;}
case 130:
#line 475 "lg.y"
{ yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);
    break;}
case 131:
#line 476 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 132:
#line 477 "lg.y"
{ yyval.args = (yyvsp[-2].args += Find(yyvsp[0].str)) ;
    break;}
case 133:
#line 478 "lg.y"
{ yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 134:
#line 479 "lg.y"
{ yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 135:
#line 482 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 136:
#line 483 "lg.y"
{yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 138:
#line 488 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[0].cexp);
    break;}
case 140:
#line 492 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 141:
#line 493 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 142:
#line 494 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 143:
#line 498 "lg.y"
{yyval.cexp=Find(yyvsp[0].str);;
    break;}
case 144:
#line 499 "lg.y"
{yyval.cexp= CConstant(yyvsp[0].lnum);
    break;}
case 145:
#line 500 "lg.y"
{yyval.cexp= CConstant(yyvsp[0].dnum);
    break;}
case 146:
#line 501 "lg.y"
{yyval.cexp= CConstant(complex<double>(0,yyvsp[0].dnum));
    break;}
case 147:
#line 502 "lg.y"
{yyval.cexp= CConstant<const char *>(yyvsp[0].str);
    break;}
case 148:
#line 503 "lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].args);;
    break;}
case 149:
#line 504 "lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].cexp);
    break;}
case 150:
#line 505 "lg.y"
{yyval.cexp=C_F0(yyvsp[-5].cexp,yyvsp[-4].oper,yyvsp[-3].cexp,yyvsp[-1].cexp);
    break;}
case 151:
#line 506 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,"[]");
    break;}
case 152:
#line 507 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].str) ;;
    break;}
case 153:
#line 508 "lg.y"
{ yyval.cexp=C_F0(Find(yyvsp[-2].str),yyvsp[0].str) ;;
    break;}
case 154:
#line 509 "lg.y"
{ yyval.cexp=C_F0(Find(yyvsp[-3].str),yyvsp[-2].oper,yyvsp[-1].args) ;;
    break;}
case 155:
#line 510 "lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 156:
#line 511 "lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 157:
#line 512 "lg.y"
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
case 158:
#line 521 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 159:
#line 522 "lg.y"
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
#line 527 "lg.y"
 


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
  currentblock = new Block(currentblock);  
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

int mymain (int  argc, char **argv)
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
  zzzfff = new  mylex(cout);
  zzzfff->input(cc);
    
  
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
   
  retvalue= Compile(); 
      
#ifdef PARALLELE
  end_parallele();
#endif
  //  currentblock->close(currentblock).eval(thestack);
  fingraphique();
  delete zzzfff;
  
   // ClearMem();
  return retvalue;
}


 
