
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

#ifdef __MWERKS__
#ifdef __INTEL__
#include <malloc.h>
#else
#include <alloca.h>
#endif
#endif
#include "AFunction.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp" 
#include "lgfem.hpp" 
#include "lex.hpp"
#include "rgraph.hpp"
class Routine;
bool load(string s);

template <class R> class FE;
template <class R,int i> class FE_;
mylex *zzzfff;
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
int kkembtype=-1;
int inloopcount=0;
Block *currentblock;
static double CPUcompileInit =0;
//class pfes;
C_F0  fespacetype;
bool fespacecomplex;

int ShowAlloc(char *s);
inline int yylex()  {return zzzfff->scan();}
inline int lineno() {return zzzfff->lineno();}
void ShowKeyWord(ostream & f ) 
 {
   zzzfff->dump(f);
 
 }

extern bool withrgraphique;
inline void fingraphique()
 { if(withrgraphique) 
   { withrgraphique=false;
    rattente(1);
    closegraphique();
  }}

void yyerror (const char* s) 
{
  cerr << endl;
  cerr <<" Error line number " <<lineno() << ", in file " << zzzfff->filename() 
       <<", before  token " <<zzzfff->YYText() << endl
       << s << endl;
   throw(ErrorCompile(s,lineno(),zzzfff->YYText() ));
}


#line 80 "lg.y"
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



#define	YYFINAL		342
#define	YYFLAG		-32768
#define	YYNTBASE	70

#define YYTRANSLATE(x) ((unsigned)(x) <= 299 ? yytranslate[x] : 110)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    30,     2,     2,     2,    24,    13,    32,    34,
    37,    22,    20,     5,    21,    36,    23,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,    69,    66,    16,
     6,    17,     2,     2,     2,     2,     2,     2,     2,     2,
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
   393,   397,   399,   403,   407,   411,   415,   419,   423,   427,
   431,   435,   439,   443,   447,   451,   455,   459,   463,   467,
   471,   475,   477,   479,   483,   489,   490,   492,   496,   498,
   502,   506,   512,   514,   518,   520,   523,   525,   529,   533,
   536,   538,   540,   542,   544,   546,   551,   556,   563,   567,
   571,   575,   580,   583,   586,   591,   595
};

static const short yyrhs[] = {    71,
    46,     0,    72,     0,    97,     0,    72,    97,     0,     0,
    75,     0,    75,     6,   102,     0,    57,    75,     0,    57,
    13,    75,     0,    78,    75,     0,    78,    13,    75,     0,
    35,    73,    38,     0,    73,     5,    75,     0,    73,     5,
    35,    73,    38,     0,    73,     5,    75,     6,   102,     0,
    73,     5,    57,    75,     0,    73,     5,    57,    13,    75,
     0,    73,     5,    78,    75,     0,    73,     5,    78,    13,
    75,     0,    75,     0,    74,     5,    75,     0,    42,     0,
    57,     0,    42,     0,    42,     6,   102,     0,    42,    34,
    77,    37,     0,    76,     5,    76,     0,   103,     0,    57,
    42,     0,    42,     6,   103,     0,    77,     5,   103,     0,
    77,     5,    75,     6,   103,     0,    55,     0,    55,    35,
    55,    38,     0,    55,    35,    55,     5,    55,    38,     0,
    55,    16,    55,    17,     0,    42,     0,    42,    35,   103,
    38,     0,    42,     6,   103,     0,    35,    74,    38,     0,
    35,    74,    38,    35,   103,    38,     0,    35,    74,    38,
     6,   103,     0,    42,    34,   103,    37,     0,    35,    74,
    38,    34,   103,    37,     0,    57,     0,    57,    16,    55,
    17,     0,    80,     0,    82,     5,    80,     0,    79,     0,
    83,     5,    79,     0,    81,    83,     0,    81,    35,    55,
    38,    82,     0,    42,    34,    77,    37,     0,    85,     0,
    86,     5,    85,     0,     0,    78,    88,    76,    66,     0,
    43,    86,    66,     0,    84,    66,     0,    56,    42,     6,
   100,    66,     0,     0,    56,    78,    42,    34,    73,    37,
    89,    67,    72,    68,     0,     0,    56,    42,    34,    73,
    37,    90,     6,   102,    66,     0,    67,     0,    68,     0,
    50,     0,    51,     0,     0,    78,    96,    76,     0,    66,
     0,    47,    45,     0,    48,    45,     0,   100,    66,     0,
    87,     0,    93,    34,   100,    66,   100,    66,   100,    37,
    97,     0,    93,    34,    95,    66,   100,    66,   100,    37,
    97,     0,    94,    34,   100,    37,    97,     0,     3,    34,
   100,    37,    97,     0,     3,    34,   100,    37,    97,     4,
    97,     0,    91,    72,    92,     0,    63,    42,    99,     0,
    63,    42,    35,   106,    38,    66,     0,    52,    66,     0,
    53,    66,     0,    54,   100,    66,     0,    34,    42,     6,
   100,     5,   100,    37,     0,    98,    97,     0,   102,     0,
   100,     5,   100,     0,    21,     0,    20,     0,    30,     0,
    28,     0,    29,     0,   103,     0,   103,     6,   102,     0,
   103,    58,   102,     0,   103,    59,   102,     0,   103,    60,
   102,     0,   103,    61,   102,     0,   107,     0,   103,    22,
   103,     0,   103,    25,   103,     0,   103,    26,   103,     0,
   103,    23,   103,     0,   103,    24,   103,     0,   103,    20,
   103,     0,   103,    21,   103,     0,   103,     8,   103,     0,
   103,     9,   103,     0,   103,    13,   103,     0,   103,    12,
   103,     0,   103,    11,   103,     0,   103,    10,   103,     0,
   103,    16,   103,     0,   103,    18,   103,     0,   103,    17,
   103,     0,   103,    19,   103,     0,   103,    14,   103,     0,
   103,    15,   103,     0,   103,     0,    69,     0,   103,    69,
   103,     0,   103,    69,   103,    69,   103,     0,     0,    57,
     0,    75,     6,   103,     0,   104,     0,   105,     5,    57,
     0,   105,     5,   104,     0,   105,     5,    75,     6,   103,
     0,   102,     0,   106,     5,   102,     0,   108,     0,   101,
   108,     0,   109,     0,   109,    31,   107,     0,   109,    33,
   107,     0,   109,    32,     0,    42,     0,    39,     0,    40,
     0,    41,     0,    45,     0,   109,    34,   105,    37,     0,
   109,    35,   104,    38,     0,   109,    35,   104,     5,   104,
    38,     0,   109,    35,    38,     0,   109,    36,    42,     0,
    57,    36,    42,     0,    57,    34,   105,    37,     0,   109,
    28,     0,   109,    29,     0,    55,    34,   100,    37,     0,
    34,   100,    37,     0,    35,   106,    38,     0
};

#endif

#if YYDEBUG != 0
static const short yyrline[] = { 0,
   198,   233,   236,   237,   240,   241,   242,   243,   244,   245,
   246,   247,   248,   249,   250,   251,   252,   253,   254,   257,
   258,   261,   261,   263,   264,   265,   266,   272,   274,   275,
   276,   277,   280,   281,   282,   283,   289,   291,   292,   293,
   294,   295,   297,   299,   303,   304,   309,   310,   312,   313,
   315,   316,   319,   323,   324,   327,   327,   328,   329,   330,
   331,   336,   341,   343,   348,   349,   351,   352,   354,   356,
   358,   359,   360,   361,   362,   363,   364,   368,   369,   370,
   371,   373,   375,   378,   382,   386,   394,   400,   405,   407,
   411,   413,   414,   415,   416,   419,   421,   422,   423,   424,
   425,   428,   430,   431,   432,   433,   434,   435,   436,   437,
   438,   439,   440,   441,   442,   443,   444,   445,   446,   447,
   448,   452,   454,   455,   456,   459,   460,   461,   462,   463,
   464,   465,   468,   469,   472,   474,   477,   478,   479,   480,
   483,   485,   486,   487,   488,   489,   490,   491,   492,   493,
   494,   495,   496,   497,   498,   507,   508
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
"'}'","':'","start","input","instructions","list_of_id_args","list_of_id1","id",
"list_of_dcls","parameters_list","type_of_dcl","ID_space","ID_array_space","fespace",
"spaceIDa","spaceIDb","spaceIDs","fespace_def","fespace_def_list","declaration",
"@1","@2","@3","begin","end","for_loop","while_loop","declaration_for","@4",
"instruction","bornes","border_expr","Expr","unop","no_comma_expr","no_set_expr",
"sub_script_expr","parameters","array","unary_expr","pow_expr","primary", NULL
};
#endif

static const short yyr1[] = {     0,
    70,    71,    72,    72,    73,    73,    73,    73,    73,    73,
    73,    73,    73,    73,    73,    73,    73,    73,    73,    74,
    74,    75,    75,    76,    76,    76,    76,    77,    77,    77,
    77,    77,    78,    78,    78,    78,    79,    79,    79,    79,
    79,    79,    80,    80,    81,    81,    82,    82,    83,    83,
    84,    84,    85,    86,    86,    88,    87,    87,    87,    87,
    89,    87,    90,    87,    91,    92,    93,    94,    96,    95,
    97,    97,    97,    97,    97,    97,    97,    97,    97,    97,
    97,    97,    97,    97,    97,    97,    98,    99,   100,   100,
   101,   101,   101,   101,   101,   102,   102,   102,   102,   102,
   102,   103,   103,   103,   103,   103,   103,   103,   103,   103,
   103,   103,   103,   103,   103,   103,   103,   103,   103,   103,
   103,   104,   104,   104,   104,   105,   105,   105,   105,   105,
   105,   105,   106,   106,   107,   107,   108,   108,   108,   108,
   109,   109,   109,   109,   109,   109,   109,   109,   109,   109,
   109,   109,   109,   109,   109,   109,   109
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
     3,     1,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     1,     1,     3,     5,     0,     1,     3,     1,     3,
     3,     5,     1,     3,     1,     2,     1,     3,     3,     2,
     1,     1,     1,     1,     1,     4,     4,     6,     3,     3,
     3,     4,     2,     2,     4,     3,     3
};

static const short yydefact[] = {     0,
     0,    92,    91,    94,    95,    93,     0,     0,   142,   143,
   144,   141,     0,   145,     0,     0,    67,    68,     0,     0,
     0,    33,     0,    45,     0,    71,    65,     0,     2,    56,
     0,     0,    75,     0,     0,     0,     3,     0,     0,    89,
    96,   102,   135,   137,     0,     0,     0,     0,   133,     0,
     0,    54,     0,    72,    73,    84,    85,     0,     0,     0,
     0,     0,    33,     0,     0,   126,     0,     0,     1,     4,
     0,     0,    37,    49,    51,    59,     0,     0,     0,     0,
    74,   136,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,   153,   154,     0,   140,
     0,   126,     0,     0,     0,   156,     0,   157,     0,     0,
    58,    86,     0,     0,     0,     0,     5,     0,     0,   141,
   127,   123,     0,   122,   129,     0,   151,     0,     0,     0,
    82,    24,     0,    22,     0,    23,     0,    20,     0,     0,
     0,    66,    81,    69,     0,     0,     0,    90,    97,   110,
   111,   115,   114,   113,   112,   120,   121,   116,   118,   117,
   119,   108,   109,   103,   106,   107,   104,   105,    98,    99,
   100,   101,   138,   139,     0,   149,     0,   150,     0,   134,
   141,     0,     0,    28,    55,    36,   155,     0,    34,     0,
     5,    23,     0,     6,     0,     5,    46,     0,     0,     0,
   152,     0,     0,    88,     0,     0,     0,    57,     0,     0,
    40,    39,     0,     0,    50,     0,     0,     0,     0,   146,
     0,   147,    79,     0,    29,     0,    53,     0,    60,     0,
     0,     8,     0,    63,     0,     0,    10,     0,   128,   124,
   130,     0,   131,     0,     0,    25,     0,    27,     0,     0,
    47,    52,    21,     0,     0,    38,    70,     0,     0,    78,
     0,     0,    30,    23,     0,    31,    35,    12,     9,     5,
    23,    13,     0,     0,     7,    11,    61,     0,     0,     0,
    83,    26,     0,     0,     0,    42,     0,     0,     0,   148,
    80,     0,     0,     0,    16,     0,     0,    18,     0,     0,
   125,   132,     0,     0,     0,    48,    41,     0,     0,    32,
    14,    17,    15,    19,     0,     0,    90,     0,    43,     0,
     0,    64,     0,    87,     0,    77,    76,    62,    44,     0,
     0,     0
};

static const short yydefgoto[] = {   340,
    28,    29,   203,   147,   204,   143,   193,    30,    74,   261,
    31,   262,    75,    32,    52,    53,    33,    71,   310,   284,
    34,   153,    35,    36,   155,   226,    37,   140,   141,    38,
    39,    40,    41,   135,   136,    50,    42,    43,    44
};

static const short yypact[] = {   392,
    14,-32768,-32768,-32768,-32768,-32768,   577,   577,-32768,-32768,
-32768,-32768,    -1,-32768,    48,    89,-32768,-32768,   103,   121,
   577,   161,    74,   214,    17,-32768,-32768,   145,   392,-32768,
    82,   175,-32768,   392,   154,   177,-32768,     0,   194,-32768,
   473,-32768,-32768,   240,   577,   197,   -24,    81,-32768,     6,
   203,-32768,     2,-32768,-32768,-32768,-32768,     8,   190,   577,
   191,    76,   163,   210,   198,   483,   216,   192,-32768,-32768,
   223,   200,   100,-32768,   261,-32768,   310,   607,   577,   577,
-32768,-32768,   577,   577,   577,   577,   577,   577,   577,   577,
   577,   577,   577,   577,   577,   577,   577,   577,   577,   577,
   577,   577,   577,   577,   577,   577,-32768,-32768,   577,-32768,
   577,   483,    22,   236,   122,-32768,   577,-32768,   637,    -1,
-32768,-32768,   262,   123,    49,   577,   150,   247,   265,   277,
   176,-32768,   278,     7,-32768,   125,-32768,   244,   577,   392,
-32768,   166,    29,-32768,   249,-32768,    50,-32768,   577,   577,
   180,-32768,-32768,-32768,   227,    30,   126,-32768,-32768,   856,
   856,   871,   871,   884,   884,   728,   728,   275,   275,   275,
   275,   286,   286,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,   127,-32768,    73,-32768,   392,-32768,
   283,   225,   131,   839,-32768,-32768,-32768,   259,-32768,    31,
   150,    27,   152,   309,    33,   150,-32768,   577,   577,   515,
-32768,   311,    75,-32768,   577,   637,   223,-32768,   182,   142,
   159,   839,   717,   142,-32768,   223,   577,   577,   392,-32768,
   546,-32768,   312,   577,-32768,   667,-32768,   282,-32768,    84,
   142,-32768,   183,-32768,   577,   142,-32768,   153,   839,   452,
   176,   315,-32768,   577,   256,-32768,   165,-32768,   142,   289,
-32768,   319,-32768,   577,   577,-32768,   320,    32,    42,-32768,
   288,   392,   839,   -24,   321,   839,-32768,-32768,-32768,   150,
    52,   322,    57,   323,-32768,-32768,-32768,   577,   577,   327,
-32768,-32768,    87,   577,   182,   839,   748,   577,   577,-32768,
-32768,   577,    95,   142,-32768,   577,   142,-32768,   577,   266,
   839,   839,   577,   300,   779,-32768,-32768,   169,   171,   839,
-32768,-32768,-32768,-32768,   269,   392,   299,   577,-32768,   392,
   392,-32768,   351,-32768,   809,-32768,-32768,-32768,-32768,   337,
   341,-32768
};

static const short yypgoto[] = {-32768,
-32768,   -32,  -197,    88,   -27,   -46,   130,   -20,   208,    79,
-32768,-32768,-32768,-32768,   248,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,   -28,-32768,-32768,    -7,
-32768,    -2,    54,  -105,   257,   231,    12,   336,-32768
};


#define	YYLAST		910


static const short yytable[] = {    48,
    70,    77,    64,   240,    80,    49,   120,   187,   248,    66,
   117,    67,    80,    58,    84,    85,    86,    87,    88,    89,
    90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
   100,   101,   102,   217,    80,    80,    80,   115,   133,   241,
    51,     2,     3,   118,   148,   246,    80,    45,    70,     4,
     5,     6,   124,   198,   220,     7,     8,   154,    68,   186,
     9,    10,    11,    12,   304,    81,    14,   121,   144,   307,
   156,   157,   158,   122,   144,   209,    46,   231,    47,   117,
   159,   126,   303,   146,   133,    80,   199,   221,   243,   146,
   132,   220,    54,   144,   218,   228,   239,   298,   144,   243,
   179,   180,   181,   182,   253,   149,   205,   299,   146,   127,
   232,   214,   255,   146,   190,    62,    72,   116,   200,   134,
   183,   278,   184,    73,   314,   271,    80,    80,    63,   210,
    80,   210,   321,    55,   150,   236,    49,   160,   161,   162,
   163,   164,   165,   166,   167,   168,   169,   170,   171,   172,
   173,   174,   175,   176,   177,   178,   243,   243,   189,   197,
   233,   211,   229,   230,   264,   134,   134,   237,    56,   236,
   258,   215,   194,    80,   242,    80,    59,   247,    59,   267,
   205,   -23,   252,   144,   201,   205,    57,    78,   244,   287,
    69,   144,   263,   265,    60,    61,   148,    61,   146,   216,
   270,   292,   222,   223,    63,   330,   202,   331,   275,    66,
    79,    67,   256,   279,   224,   282,   259,   280,   286,   268,
   269,    73,   283,   260,   144,   138,   139,     7,     8,    65,
    60,   148,     9,    10,    11,    12,   119,    63,    14,   281,
    76,   144,   285,   301,   123,   125,   290,    66,    46,    67,
    47,   128,   129,   305,   145,   308,   146,   137,    66,   205,
    67,   249,   250,   134,   142,   151,   235,   107,   108,   194,
   109,   110,   111,   112,   113,   114,   322,   188,   196,   324,
   206,   207,   -22,   208,   134,   212,   219,   273,   234,   276,
   318,   319,   227,   333,    96,    97,    98,    99,   100,   101,
   102,   336,   337,   323,    70,   327,   325,    98,    99,   100,
   101,   102,     1,   238,   245,   272,   254,   296,   297,   277,
   289,   291,   294,   295,   217,   300,   302,   306,   309,     2,
     3,   313,   326,   328,   332,   334,   341,     4,     5,     6,
   342,   311,   312,     7,     8,   257,   293,   315,     9,    10,
    11,    12,    13,     1,    14,   320,    15,    16,   225,    17,
    18,    19,    20,    21,    22,    23,    24,   195,   185,   213,
     2,     3,    25,   316,    82,    26,    27,   152,     4,     5,
     6,   335,     0,     0,     7,     8,     0,     0,     0,     9,
    10,    11,    12,    13,     1,    14,     0,    15,    16,     0,
    17,    18,    19,    20,    21,    22,    23,    24,     0,     0,
     0,     2,     3,    25,     0,     0,    26,    27,   338,     4,
     5,     6,     0,     0,     0,     7,     8,     0,     0,     0,
     9,    10,    11,    12,    13,     0,    14,     0,    15,    16,
     0,    17,    18,    19,    20,    21,    22,    23,    24,     0,
     0,     0,     0,     0,    25,     0,     0,    26,    27,    84,
    85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
    95,    96,    97,    98,    99,   100,   101,   102,    83,     0,
    84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
    94,    95,    96,    97,    98,    99,   100,   101,   102,     0,
     0,     0,     2,     3,     0,     0,     0,     0,     0,     0,
     4,     5,     6,     0,     0,     0,     7,     8,     0,     0,
   288,     9,    10,    11,   130,     0,     0,    14,     0,     0,
   103,   104,   105,   106,     2,     3,     0,    46,     0,   131,
     0,     0,     4,     5,     6,     0,     0,     0,     7,     8,
     0,   132,     0,     9,    10,    11,   130,     0,     0,    14,
     0,     0,     0,     0,     0,     2,     3,     0,     0,    46,
     0,   251,     0,     4,     5,     6,     0,     0,     0,     7,
     8,     0,     0,   132,     9,    10,    11,    12,     0,     0,
    14,     0,     0,     0,     0,     0,     2,     3,     0,     0,
    46,     0,    47,     0,     4,     5,     6,     0,     0,     0,
     7,     8,     0,     0,   132,     9,    10,    11,    12,     0,
     0,    14,     0,     0,     0,     0,     2,     3,     0,     0,
     0,    46,     0,    47,     4,     5,     6,     0,     0,     0,
     7,     8,     0,     0,     0,     9,    10,    11,    12,     0,
     0,    14,     0,     0,     0,     0,     2,     3,     0,     0,
     0,    22,     0,    47,     4,     5,     6,     0,     0,     0,
     7,     8,     0,     0,     0,     9,    10,    11,   191,     0,
     0,    14,     0,     0,     0,     0,     2,     3,     0,     0,
     0,    46,     0,   192,     4,     5,     6,     0,     0,     0,
     7,     8,     0,     0,     0,     9,    10,    11,   130,     0,
     0,    14,     0,     0,     0,     0,     0,     0,     0,     0,
     0,    46,     0,   274,    84,    85,    86,    87,    88,    89,
    90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
   100,   101,   102,    92,    93,    94,    95,    96,    97,    98,
    99,   100,   101,   102,   266,    84,    85,    86,    87,    88,
    89,    90,    91,    92,    93,    94,    95,    96,    97,    98,
    99,   100,   101,   102,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,   317,    84,    85,    86,    87,
    88,    89,    90,    91,    92,    93,    94,    95,    96,    97,
    98,    99,   100,   101,   102,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,   329,    84,    85,    86,    87,
    88,    89,    90,    91,    92,    93,    94,    95,    96,    97,
    98,    99,   100,   101,   102,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,   339,    84,    85,    86,    87,
    88,    89,    90,    91,    92,    93,    94,    95,    96,    97,
    98,    99,   100,   101,   102,    86,    87,    88,    89,    90,
    91,    92,    93,    94,    95,    96,    97,    98,    99,   100,
   101,   102,    88,    89,    90,    91,    92,    93,    94,    95,
    96,    97,    98,    99,   100,   101,   102,    90,    91,    92,
    93,    94,    95,    96,    97,    98,    99,   100,   101,   102
};

static const short yycheck[] = {     7,
    29,    34,    23,   201,     5,     8,     5,   113,   206,    34,
     5,    36,     5,    21,     8,     9,    10,    11,    12,    13,
    14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
    24,    25,    26,     5,     5,     5,     5,    45,    66,    13,
    42,    20,    21,    38,    72,    13,     5,    34,    77,    28,
    29,    30,    60,     5,     5,    34,    35,    78,    42,    38,
    39,    40,    41,    42,    13,    66,    45,    66,    42,    13,
    78,    79,    80,    66,    42,    69,    55,     5,    57,     5,
    83,     6,   280,    57,   112,     5,    38,    38,     5,    57,
    69,     5,    45,    42,    66,    66,    66,    66,    42,     5,
   103,   104,   105,   106,   210,     6,   127,    66,    57,    34,
    38,   140,    38,    57,   117,    42,    35,    37,   126,    66,
   109,    38,   111,    42,    38,   231,     5,     5,    55,     5,
     5,     5,    38,    45,    35,     5,   139,    84,    85,    86,
    87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
    97,    98,    99,   100,   101,   102,     5,     5,    37,    37,
   189,    37,    37,    37,     6,   112,   113,    37,    66,     5,
   217,     6,   119,     5,   202,     5,    16,   205,    16,   226,
   201,     6,   210,    42,    35,   206,    66,    34,    37,    37,
    46,    42,   220,    35,    34,    35,   224,    35,    57,    34,
   229,    37,   149,   150,    55,    37,    57,    37,   236,    34,
    34,    36,   215,   241,    35,   243,    35,    35,   246,   227,
   228,    42,   243,    42,    42,    34,    35,    34,    35,    16,
    34,   259,    39,    40,    41,    42,    34,    55,    45,    57,
    66,    42,   245,   272,    55,    55,   254,    34,    55,    36,
    57,    42,    55,   281,    55,   283,    57,    42,    34,   280,
    36,   208,   209,   210,    42,     5,    42,    28,    29,   216,
    31,    32,    33,    34,    35,    36,   304,    42,    17,   307,
    34,    17,     6,     6,   231,    42,    38,   234,     6,   236,
   298,   299,    66,   326,    20,    21,    22,    23,    24,    25,
    26,   330,   331,   306,   333,   313,   309,    22,    23,    24,
    25,    26,     3,    55,     6,     4,     6,   264,   265,    38,
     6,    66,    34,     5,     5,    38,     6,     6,     6,    20,
    21,     5,    67,    34,    66,    37,     0,    28,    29,    30,
     0,   288,   289,    34,    35,   216,   259,   294,    39,    40,
    41,    42,    43,     3,    45,   302,    47,    48,   151,    50,
    51,    52,    53,    54,    55,    56,    57,   120,   112,   139,
    20,    21,    63,   295,    39,    66,    67,    68,    28,    29,
    30,   328,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,
    40,    41,    42,    43,     3,    45,    -1,    47,    48,    -1,
    50,    51,    52,    53,    54,    55,    56,    57,    -1,    -1,
    -1,    20,    21,    63,    -1,    -1,    66,    67,    68,    28,
    29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
    39,    40,    41,    42,    43,    -1,    45,    -1,    47,    48,
    -1,    50,    51,    52,    53,    54,    55,    56,    57,    -1,
    -1,    -1,    -1,    -1,    63,    -1,    -1,    66,    67,     8,
     9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
    19,    20,    21,    22,    23,    24,    25,    26,     6,    -1,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
    18,    19,    20,    21,    22,    23,    24,    25,    26,    -1,
    -1,    -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    -1,
    28,    29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,
    69,    39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,
    58,    59,    60,    61,    20,    21,    -1,    55,    -1,    57,
    -1,    -1,    28,    29,    30,    -1,    -1,    -1,    34,    35,
    -1,    69,    -1,    39,    40,    41,    42,    -1,    -1,    45,
    -1,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    55,
    -1,    57,    -1,    28,    29,    30,    -1,    -1,    -1,    34,
    35,    -1,    -1,    69,    39,    40,    41,    42,    -1,    -1,
    45,    -1,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,
    55,    -1,    57,    -1,    28,    29,    30,    -1,    -1,    -1,
    34,    35,    -1,    -1,    69,    39,    40,    41,    42,    -1,
    -1,    45,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,
    -1,    55,    -1,    57,    28,    29,    30,    -1,    -1,    -1,
    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,
    -1,    45,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,
    -1,    55,    -1,    57,    28,    29,    30,    -1,    -1,    -1,
    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,
    -1,    45,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,
    -1,    55,    -1,    57,    28,    29,    30,    -1,    -1,    -1,
    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,
    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    55,    -1,    57,     8,     9,    10,    11,    12,    13,
    14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
    24,    25,    26,    16,    17,    18,    19,    20,    21,    22,
    23,    24,    25,    26,    38,     8,     9,    10,    11,    12,
    13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
    23,    24,    25,    26,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    38,     8,     9,    10,    11,
    12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
    22,    23,    24,    25,    26,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    37,     8,     9,    10,    11,
    12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
    22,    23,    24,    25,    26,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    37,     8,     9,    10,    11,
    12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
    22,    23,    24,    25,    26,    10,    11,    12,    13,    14,
    15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
    25,    26,    12,    13,    14,    15,    16,    17,    18,    19,
    20,    21,    22,    23,    24,    25,    26,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    23,    24,    25,    26
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
#line 198 "lg.y"
{
	            
                        size_t sizestack = currentblock->size()+1024 ; //  before close 
                        yyvsp[-1].cinst+=currentblock->close(currentblock);
                        cout << " sizestack + 1024 =" << sizestack << "  ( " << sizestack-1024 <<" )\n" ;                         
                        int NbPtr = ShowAlloc("init execution "); // number of un delele ptr
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
                        NbPtr = ShowAlloc("end execution -- ") - NbPtr;
                        
                        if (NbPtr) { cout << " ######## We forget of deleting   " << NbPtr << " Nb pointer  " << endl;}
  return 0;;
    break;}
case 3:
#line 236 "lg.y"
{yyval.cinst=yyvsp[0].cexp;;;;
    break;}
case 4:
#line 237 "lg.y"
{ yyval.cinst= (yyvsp[-1].cinst+=yyvsp[0].cexp) ;
    break;}
case 5:
#line 240 "lg.y"
{ yyval.clist_id=new ListOfId();;
    break;}
case 6:
#line 241 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str));
    break;}
case 7:
#line 242 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 8:
#line 243 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>()));
    break;}
case 9:
#line 244 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true));
    break;}
case 10:
#line 245 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 11:
#line 246 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 12:
#line 247 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 13:
#line 248 "lg.y"
{ yyval.clist_id = yyvsp[-2].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str)) ;
    break;}
case 14:
#line 249 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 15:
#line 250 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 16:
#line 251 "lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>())) ;
    break;}
case 17:
#line 252 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true)) ;
    break;}
case 18:
#line 253 "lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 19:
#line 254 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 20:
#line 257 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 21:
#line 258 "lg.y"
{ yyval.clist_id=yyvsp[-2].clist_id  ; yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 24:
#line 263 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[0].str,dcltype);
    break;}
case 25:
#line 264 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-2].str,dcltype,yyvsp[0].cexp);
    break;}
case 26:
#line 265 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-3].str,dcltype,yyvsp[-1].args);
    break;}
case 27:
#line 266 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 28:
#line 273 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 29:
#line 274 "lg.y"
{yyval.args=Find(yyvsp[-1].str);
    break;}
case 30:
#line 275 "lg.y"
{ yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);
    break;}
case 31:
#line 276 "lg.y"
{ yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 32:
#line 277 "lg.y"
{ yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp));
    break;}
case 34:
#line 281 "lg.y"
{yyval.type=TypeArray(yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 35:
#line 282 "lg.y"
{yyval.type=TypeArray(yyvsp[-5].type,yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 36:
#line 283 "lg.y"
{yyval.type=TypeTemplate(yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 37:
#line 290 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[0].str,currentblock,fespacetype,fespacecomplex); ;
    break;}
case 38:
#line 291 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;
    break;}
case 39:
#line 292 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-2].str,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;
    break;}
case 40:
#line 293 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-1].clist_id,currentblock,fespacetype,fespacecomplex) ;
    break;}
case 41:
#line 294 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;
    break;}
case 42:
#line 295 "lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-3].clist_id,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;
    break;}
case 43:
#line 298 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;
    break;}
case 44:
#line 299 "lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;
    break;}
case 45:
#line 303 "lg.y"
{fespacecomplex=false;  fespacetype = Find(yyvsp[0].str);;
    break;}
case 46:
#line 304 "lg.y"
{
             if (yyvsp[-1].type != typevarreal && yyvsp[-1].type != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=(yyvsp[-1].type==typevarcomplex);
             fespacetype = Find(yyvsp[-3].str);;
    break;}
case 47:
#line 309 "lg.y"
{  yyval.cexp = yyvsp[0].cexp  ;
    break;}
case 48:
#line 310 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 49:
#line 312 "lg.y"
{  yyval.cexp = yyvsp[0].cexp  ;
    break;}
case 50:
#line 313 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 51:
#line 315 "lg.y"
{ yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;
    break;}
case 52:
#line 316 "lg.y"
{ yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;
    break;}
case 53:
#line 321 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariableFES,size_t>(yyvsp[-3].str,atype<pfes*>(),yyvsp[-1].args,dimFESpaceImage(yyvsp[-1].args));
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
                      yyvsp[-1].routine=new Routine(yyvsp[-5].type,yyvsp[-4].type->right(),yyvsp[-3].str,yyvsp[-1].clist_id,currentblock);;
    break;}
case 62:
#line 337 "lg.y"
{ currentblock=yyvsp[-5].routine->Set(yyvsp[-1].cinst);
                       currentblock->Add(yyvsp[-7].str,"(",yyvsp[-5].routine);
                       kkembtype--;
                       yyval.cexp=0 ;
    break;}
case 63:
#line 342 "lg.y"
{currentblock = new Block(currentblock); yyvsp[-4].type->SetArgs(yyvsp[-1].clist_id);;
    break;}
case 64:
#line 344 "lg.y"
{  yyval.cinst=currentblock->close(currentblock);
                         yyval.cexp=currentblock->NewID(yyvsp[-8].type,yyvsp[-7].str,yyvsp[-1].cexp,*yyvsp[-5].clist_id);
    break;}
case 65:
#line 348 "lg.y"
{  currentblock = new Block(currentblock);
    break;}
case 66:
#line 349 "lg.y"
{  yyval.cexp=currentblock->close(currentblock);
    break;}
case 67:
#line 351 "lg.y"
{inloopcount++;;
    break;}
case 68:
#line 352 "lg.y"
{inloopcount++;
    break;}
case 69:
#line 355 "lg.y"
{dcltype=yyvsp[0].type;currentblock = new Block(currentblock);
    break;}
case 70:
#line 356 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 71:
#line 358 "lg.y"
{yyval.cexp=0;;
    break;}
case 72:
#line 359 "lg.y"
{zzzfff->input(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 73:
#line 360 "lg.y"
{load(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 74:
#line 361 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 75:
#line 362 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 76:
#line 363 "lg.y"
{inloopcount--; yyval.cexp=For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 77:
#line 365 "lg.y"
{inloopcount--; 
                yyval.cexp=C_F0(For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp),currentblock->close(currentblock));
    break;}
case 78:
#line 368 "lg.y"
{inloopcount--;yyval.cexp=While(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 79:
#line 369 "lg.y"
{yyval.cexp=FIf(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 80:
#line 370 "lg.y"
{yyval.cexp=FIf(yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 81:
#line 371 "lg.y"
{ 
                      yyval.cexp=C_F0(new E_block(yyvsp[-1].cinst,yyvsp[0].cexp),atype<void>()) ;
    break;}
case 82:
#line 373 "lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-1].str,C_F0(TheOperators,"[border]",yyvsp[0].args));
    break;}
case 83:
#line 375 "lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-4].str,C_F0(TheOperators,"[border]",yyvsp[-2].args));
    break;}
case 84:
#line 378 "lg.y"
{
                    if(inloopcount) 
                      yyval.cexp= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;
    break;}
case 85:
#line 382 "lg.y"
{ 
                    if(inloopcount)
                        yyval.cexp= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");
    break;}
case 86:
#line 386 "lg.y"
{ 
                    if (kkembtype>=0)
                      yyval.cexp= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo(yyvsp[-1].cexp)) ,atype<void>());
                     else lgerror(" return not in routine ") ;
    break;}
case 87:
#line 394 "lg.y"
{ 
   currentblock = new Block(currentblock);
   yyval.args = currentblock->NewVar<LocalVariable>(yyvsp[-5].str,atype<double*>());
   yyval.args+= yyvsp[-3].cexp;
   yyval.args+= yyvsp[-1].cexp ;
    break;}
case 88:
#line 400 "lg.y"
{   
   yyval.args = (yyvsp[-1].args += yyvsp[0].cexp);
   currentblock->close(currentblock);
    break;}
case 90:
#line 407 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 97:
#line 421 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 98:
#line 422 "lg.y"
{yyval.cexp=C_F0(TheOperators,"+=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 99:
#line 423 "lg.y"
{yyval.cexp=C_F0(TheOperators,"-=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 100:
#line 424 "lg.y"
{yyval.cexp=C_F0(TheOperators,"*=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 101:
#line 425 "lg.y"
{yyval.cexp=C_F0(TheOperators,"/=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 103:
#line 430 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 104:
#line 431 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 105:
#line 432 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 106:
#line 433 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 107:
#line 434 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 108:
#line 435 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 109:
#line 436 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 110:
#line 437 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 111:
#line 438 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 112:
#line 439 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 113:
#line 440 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 114:
#line 441 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 115:
#line 442 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 116:
#line 443 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 117:
#line 444 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 118:
#line 445 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 119:
#line 446 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 120:
#line 447 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 121:
#line 448 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 122:
#line 453 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 123:
#line 454 "lg.y"
{yyval.cexp=C_F0(TheOperators,":");
    break;}
case 124:
#line 455 "lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 125:
#line 456 "lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 126:
#line 459 "lg.y"
{yyval.args=0;
    break;}
case 127:
#line 460 "lg.y"
{yyval.args=Find(yyvsp[0].str);
    break;}
case 128:
#line 461 "lg.y"
{ yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);
    break;}
case 129:
#line 462 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 130:
#line 463 "lg.y"
{ yyval.args = (yyvsp[-2].args += Find(yyvsp[0].str)) ;
    break;}
case 131:
#line 464 "lg.y"
{ yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 132:
#line 465 "lg.y"
{ yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 133:
#line 468 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 134:
#line 469 "lg.y"
{yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 136:
#line 474 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[0].cexp);
    break;}
case 138:
#line 478 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 139:
#line 479 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 140:
#line 480 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 141:
#line 484 "lg.y"
{yyval.cexp=Find(yyvsp[0].str);;
    break;}
case 142:
#line 485 "lg.y"
{yyval.cexp= CConstant(yyvsp[0].lnum);
    break;}
case 143:
#line 486 "lg.y"
{yyval.cexp= CConstant(yyvsp[0].dnum);
    break;}
case 144:
#line 487 "lg.y"
{yyval.cexp= CConstant(complex<double>(0,yyvsp[0].dnum));
    break;}
case 145:
#line 488 "lg.y"
{yyval.cexp= CConstant<const char *>(yyvsp[0].str);
    break;}
case 146:
#line 489 "lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].args);;
    break;}
case 147:
#line 490 "lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].cexp);
    break;}
case 148:
#line 491 "lg.y"
{yyval.cexp=C_F0(yyvsp[-5].cexp,yyvsp[-4].oper,yyvsp[-3].cexp,yyvsp[-1].cexp);
    break;}
case 149:
#line 492 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,"[]");
    break;}
case 150:
#line 493 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].str) ;;
    break;}
case 151:
#line 494 "lg.y"
{ yyval.cexp=C_F0(Find(yyvsp[-2].str),yyvsp[0].str) ;;
    break;}
case 152:
#line 495 "lg.y"
{ yyval.cexp=C_F0(Find(yyvsp[-3].str),yyvsp[-2].oper,yyvsp[-1].args) ;;
    break;}
case 153:
#line 496 "lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 154:
#line 497 "lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 155:
#line 498 "lg.y"
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
case 156:
#line 507 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 157:
#line 508 "lg.y"
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
#line 513 "lg.y"
 


#include <fstream>
using namespace std;
// bool lgdebug;
 bool lexdebug;
void ForDebug();
void ForDebug()
{
  int i=0;
  i++;
}
extern void ShowAlloc(const char *s);
extern void ShowNbAlloc(const char *s);
void init_lgfem() ;
void init_lgmesh() ;
void init_algo();
bool withrgraphique = false;
const char * StrVersionNumber();

int mymain (int  argc, char **argv)
{
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
  int ok;
  
  currentblock=0;
  currentblock = new Block(currentblock);  
  try {
    retvalue=yyparse(); //  compile
    if(retvalue==0)  
      if(currentblock) 
	retvalue=1,cerr <<  "Error:a block is not close" << endl;   
      else 
	cerr <<  "Bien: On a fini Normalement" << endl; 
  }

  catch (Error & e) 
    {
      retvalue=e.errcode();
      cerr << "error " << e.what() 
	   << "\n code = "<<  retvalue << endl;
    }
#ifdef PARALLELE
  end_parallele();
#endif
  //  currentblock->close(currentblock).eval(thestack);
  fingraphique();
  return retvalue;
}


 
