
/*  A Bison parser, made from /Users/hecht/work/FreeFem++/src/lg.y
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

#line 1 "/Users/hecht/work/FreeFem++/src/lg.y"
 

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
#ifdef EIGENVALUE
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


#line 78 "/Users/hecht/work/FreeFem++/src/lg.y"
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



#define	YYFINAL		329
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
   107,   111,   117,   119,   124,   131,   133,   138,   142,   146,
   153,   159,   164,   171,   173,   175,   179,   181,   185,   188,
   194,   199,   201,   205,   206,   211,   215,   218,   224,   225,
   236,   237,   247,   249,   251,   253,   255,   256,   260,   262,
   265,   268,   271,   273,   283,   293,   299,   305,   313,   317,
   321,   328,   331,   334,   338,   346,   349,   351,   355,   357,
   359,   361,   363,   365,   367,   371,   375,   379,   383,   387,
   389,   393,   397,   401,   405,   409,   413,   417,   421,   425,
   429,   433,   437,   441,   445,   449,   453,   457,   461,   465,
   467,   469,   473,   479,   480,   482,   486,   488,   492,   496,
   502,   504,   508,   510,   513,   515,   519,   523,   526,   528,
   530,   532,   534,   536,   541,   546,   553,   557,   561,   564,
   567,   572,   576
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
    42,     0,    42,    35,   103,    38,     0,    42,     6,   103,
     0,    35,    74,    38,     0,    35,    74,    38,    35,   103,
    38,     0,    35,    74,    38,     6,   103,     0,    42,    34,
   103,    37,     0,    35,    74,    38,    34,   103,    37,     0,
    57,     0,    80,     0,    82,     5,    80,     0,    79,     0,
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
   109,    28,     0,   109,    29,     0,    55,    34,   100,    37,
     0,    34,   100,    37,     0,    35,   106,    38,     0
};

#endif

#if YYDEBUG != 0
static const short yyrline[] = { 0,
   197,   226,   229,   230,   233,   234,   235,   236,   237,   238,
   239,   240,   241,   242,   243,   244,   245,   246,   247,   250,
   251,   254,   254,   256,   257,   258,   259,   265,   267,   268,
   269,   270,   273,   274,   275,   281,   283,   284,   285,   286,
   287,   289,   291,   295,   298,   299,   301,   302,   304,   305,
   308,   312,   313,   316,   316,   317,   318,   319,   320,   325,
   330,   332,   337,   338,   340,   341,   343,   345,   347,   348,
   349,   350,   351,   352,   353,   357,   358,   359,   360,   362,
   364,   367,   371,   375,   383,   389,   394,   396,   400,   402,
   403,   404,   405,   408,   410,   411,   412,   413,   414,   417,
   419,   420,   421,   422,   423,   424,   425,   426,   427,   428,
   429,   430,   431,   432,   433,   434,   435,   436,   437,   441,
   443,   444,   445,   448,   449,   450,   451,   452,   453,   454,
   457,   458,   461,   463,   466,   467,   468,   469,   472,   474,
   475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
   485,   494,   495
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
    77,    77,    78,    78,    78,    79,    79,    79,    79,    79,
    79,    80,    80,    81,    82,    82,    83,    83,    84,    84,
    85,    86,    86,    88,    87,    87,    87,    87,    89,    87,
    90,    87,    91,    92,    93,    94,    96,    95,    97,    97,
    97,    97,    97,    97,    97,    97,    97,    97,    97,    97,
    97,    97,    97,    97,    98,    99,   100,   100,   101,   101,
   101,   101,   101,   102,   102,   102,   102,   102,   102,   103,
   103,   103,   103,   103,   103,   103,   103,   103,   103,   103,
   103,   103,   103,   103,   103,   103,   103,   103,   103,   104,
   104,   104,   104,   105,   105,   105,   105,   105,   105,   105,
   106,   106,   107,   107,   108,   108,   108,   108,   109,   109,
   109,   109,   109,   109,   109,   109,   109,   109,   109,   109,
   109,   109,   109
};

static const short yyr2[] = {     0,
     2,     1,     1,     2,     0,     1,     3,     2,     3,     2,
     3,     3,     3,     5,     5,     4,     5,     4,     5,     1,
     3,     1,     1,     1,     3,     4,     3,     1,     2,     3,
     3,     5,     1,     4,     6,     1,     4,     3,     3,     6,
     5,     4,     6,     1,     1,     3,     1,     3,     2,     5,
     4,     1,     3,     0,     4,     3,     2,     5,     0,    10,
     0,     9,     1,     1,     1,     1,     0,     3,     1,     2,
     2,     2,     1,     9,     9,     5,     5,     7,     3,     3,
     6,     2,     2,     3,     7,     2,     1,     3,     1,     1,
     1,     1,     1,     1,     3,     3,     3,     3,     3,     1,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
     1,     3,     5,     0,     1,     3,     1,     3,     3,     5,
     1,     3,     1,     2,     1,     3,     3,     2,     1,     1,
     1,     1,     1,     4,     4,     6,     3,     3,     2,     2,
     4,     3,     3
};

static const short yydefact[] = {     0,
     0,    90,    89,    92,    93,    91,     0,     0,   140,   141,
   142,   139,     0,   143,     0,     0,    65,    66,     0,     0,
     0,    33,     0,    44,     0,    69,    63,     0,     2,    54,
     0,     0,    73,     0,     0,     0,     3,     0,     0,    87,
    94,   100,   133,   135,     0,     0,     0,   131,     0,     0,
    52,     0,    70,    71,    82,    83,     0,     0,     0,     0,
    33,     0,     0,     1,     4,     0,     0,    36,    47,    49,
    57,     0,     0,     0,     0,    72,   134,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,   149,   150,     0,   138,     0,   124,     0,     0,     0,
   152,     0,   153,     0,     0,    56,    84,     0,     0,     0,
     5,     0,     0,     0,     0,    80,    24,     0,    22,     0,
    23,     0,    20,     0,     0,     0,    64,    79,    67,     0,
     0,     0,    88,    95,   108,   109,   113,   112,   111,   110,
   118,   119,   114,   116,   115,   117,   106,   107,   101,   104,
   105,   102,   103,    96,    97,    98,    99,   136,   137,   139,
   125,   121,     0,   120,   127,     0,   147,     0,   148,     0,
   132,   139,     0,     0,    28,    53,   151,     0,    34,     0,
     5,    23,     0,     6,     0,     5,     0,     0,    86,     0,
     0,     0,    55,     0,     0,    39,    38,     0,     0,    48,
     0,     0,     0,     0,     0,     0,     0,   144,     0,   145,
    77,     0,    29,     0,    51,     0,    58,     0,     0,     8,
     0,    61,     0,     0,    10,     0,     0,     0,    25,     0,
    27,     0,     0,    45,    50,    21,     0,     0,    37,    68,
     0,     0,    76,   126,   122,   128,     0,   129,     0,     0,
    30,     0,    31,    35,    12,     9,     5,    23,    13,     0,
     0,     7,    11,    59,     0,    81,    26,     0,     0,     0,
    41,     0,     0,     0,     0,     0,   146,    78,     0,     0,
     0,    16,     0,     0,    18,     0,     0,     0,     0,     0,
    46,    40,     0,     0,   123,   130,    32,    14,    17,    15,
    19,     0,     0,    88,     0,    42,     0,     0,    62,     0,
    85,     0,    75,    74,    60,    43,     0,     0,     0
};

static const short yydefgoto[] = {   327,
    28,    29,   193,   132,   194,   128,   184,    30,    69,   244,
    31,   245,    70,    32,    51,    52,    33,    66,   297,   271,
    34,   138,    35,    36,   140,   211,    37,   125,   126,    38,
    39,    40,    41,   175,   176,    49,    42,    43,    44
};

static const short yypact[] = {   417,
   -27,-32768,-32768,-32768,-32768,-32768,   292,   292,-32768,-32768,
-32768,-32768,    -5,-32768,    10,    16,-32768,-32768,    -3,    59,
   292,   131,   120,-32768,    86,-32768,-32768,   132,   417,-32768,
   155,    97,-32768,   417,   157,   165,-32768,    -2,   218,-32768,
   498,-32768,-32768,   661,   292,   177,    40,-32768,    34,   179,
-32768,     4,-32768,-32768,-32768,-32768,     5,   292,   124,    81,
   166,   172,   133,-32768,-32768,   173,   114,   145,-32768,   219,
-32768,   321,   646,   292,   292,-32768,-32768,   292,   292,   292,
   292,   292,   292,   292,   292,   292,   292,   292,   292,   292,
   292,   292,   292,   292,   292,   292,   292,   292,   292,   292,
   292,-32768,-32768,   292,-32768,   292,   153,   509,   181,    52,
-32768,   292,-32768,   586,    -5,-32768,-32768,    57,    43,   292,
   161,   191,   186,   292,   417,-32768,   142,     8,-32768,   193,
-32768,    45,-32768,   292,   292,   167,-32768,-32768,-32768,   168,
    29,    68,-32768,-32768,   833,   833,   848,   848,   861,   861,
   705,   705,   732,   732,   732,   732,    84,    84,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,   229,
   231,-32768,   232,     7,-32768,    83,-32768,    47,-32768,   417,
-32768,   233,   202,   112,   816,-32768,-32768,   190,-32768,    30,
   161,    -1,   113,   241,    23,   161,   243,    48,-32768,   292,
   586,   173,-32768,   184,    69,   154,   816,   694,    69,-32768,
   173,   292,   292,   417,   292,   292,   532,-32768,   563,-32768,
   246,   292,-32768,   616,-32768,   217,-32768,    53,    69,-32768,
   185,-32768,   292,    69,-32768,   116,   292,   195,-32768,   117,
-32768,    69,   222,-32768,   257,-32768,   292,   292,-32768,   259,
    37,    38,-32768,   816,   477,   231,   262,-32768,   237,   417,
   816,   264,   816,-32768,-32768,-32768,   161,    27,   265,    33,
   272,-32768,-32768,-32768,   274,-32768,-32768,    54,   292,   184,
   816,   725,   292,   292,   292,   292,-32768,-32768,   292,    55,
    69,-32768,   292,    69,-32768,   292,   216,   292,   250,   756,
-32768,-32768,   118,   122,   816,   816,   816,-32768,-32768,-32768,
-32768,   220,   417,   248,   292,-32768,   417,   417,-32768,   376,
-32768,   786,-32768,-32768,-32768,-32768,   287,   288,-32768
};

static const short yypgoto[] = {-32768,
-32768,   -32,  -185,    58,    12,   -41,    93,   -19,   159,    19,
-32768,-32768,-32768,-32768,   187,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,   -28,-32768,-32768,    -7,
-32768,     0,    50,  -103,-32768,   180,   -57,   266,-32768
};


#define	YYLAST		887


static const short yytable[] = {    47,
    65,    72,    75,    62,   178,   228,    45,    48,   115,    75,
   236,   229,   202,    57,    79,    80,    81,    82,    83,    84,
    85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
    95,    96,    97,    75,    75,   234,    50,   110,   112,   291,
   129,    75,    75,    65,    75,   294,   168,   188,   169,   205,
   118,   219,   112,   139,    53,   131,    75,   231,   205,   231,
    54,    75,    55,    76,   129,   141,   142,   143,   129,   116,
   117,   113,    75,   203,   129,   216,   111,   144,   133,   131,
   189,   290,   206,   131,   220,   238,   120,   217,   180,   131,
   265,   299,   308,   187,   213,   227,   199,   164,   165,   166,
   167,   195,   283,   284,   214,    93,    94,    95,    96,    97,
   129,   181,   190,   258,   121,   259,   224,   231,   173,   218,
   231,   224,    75,    48,    56,   131,    75,    63,   145,   146,
   147,   148,   149,   150,   151,   152,   153,   154,   155,   156,
   157,   158,   159,   160,   161,   162,   163,   200,   225,   232,
   134,   221,   274,   277,   317,   129,   174,   174,   318,   247,
   241,    60,    71,   185,    58,    59,   123,   124,   130,   250,
   131,   195,     2,     3,    61,   201,   195,    64,   119,   135,
     4,     5,     6,   207,   208,   253,     7,     8,   248,    67,
    73,     9,    10,    11,   170,   191,    68,    14,    74,   239,
    59,   209,   129,   230,   251,   252,   235,    46,    68,   171,
    58,   270,   114,   122,   127,    61,   246,   192,   242,   267,
   133,   172,   179,   136,   196,   243,   129,   197,   257,   275,
   204,   288,   272,   212,   -22,   262,   -23,   215,   222,    61,
   266,   268,   269,   223,   226,   273,   233,   195,   237,   260,
   185,     7,     8,   133,   264,   279,     9,    10,    11,    12,
   276,   280,    14,   202,   254,   255,   174,   286,   174,   289,
   293,   261,    46,   263,   287,   303,   304,   296,   298,   292,
   320,   295,   313,   315,   321,   319,   328,   329,   323,   324,
   314,    65,   310,   240,   210,   312,   281,   282,   301,   278,
     0,   186,   309,   198,    77,   311,     0,     0,     0,     0,
     0,     2,     3,     0,     0,     0,     0,     0,     0,     4,
     5,     6,     0,     1,     0,     7,     8,     0,   300,     0,
     9,    10,    11,    12,   305,   306,    14,     0,   307,     0,
     2,     3,     0,     0,     0,     0,    46,     0,     4,     5,
     6,     0,     0,     0,     7,     8,     0,     0,     0,     9,
    10,    11,    12,    13,   322,    14,     0,    15,    16,     0,
    17,    18,    19,    20,    21,    22,    23,    24,     1,     0,
     0,     0,     0,    25,     0,     0,    26,    27,   137,     0,
     0,     0,     0,     0,     0,     2,     3,     0,     0,     0,
     0,     0,     0,     4,     5,     6,     0,     0,     0,     7,
     8,     0,     0,     0,     9,    10,    11,    12,    13,     1,
    14,     0,    15,    16,     0,    17,    18,    19,    20,    21,
    22,    23,    24,     0,     0,     0,     2,     3,    25,     0,
     0,    26,    27,   325,     4,     5,     6,     0,     0,     0,
     7,     8,     0,     0,     0,     9,    10,    11,    12,    13,
     0,    14,     0,    15,    16,     0,    17,    18,    19,    20,
    21,    22,    23,    24,     0,     0,     0,     0,     0,    25,
     0,     0,    26,    27,    79,    80,    81,    82,    83,    84,
    85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
    95,    96,    97,    78,     0,    79,    80,    81,    82,    83,
    84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
    94,    95,    96,    97,     0,     0,     0,     0,     2,     3,
     0,     0,     0,     0,     0,     0,     4,     5,     6,     0,
     0,     0,     7,     8,     0,   285,   177,     9,    10,    11,
    12,     2,     3,    14,     0,    98,    99,   100,   101,     4,
     5,     6,     0,    46,     0,     7,     8,     0,     0,     0,
     9,    10,    11,   170,     0,     0,    14,   172,     0,     0,
     0,     0,     2,     3,     0,     0,    46,     0,   256,     0,
     4,     5,     6,     0,     0,     0,     7,     8,     0,     0,
   172,     9,    10,    11,    12,     2,     3,    14,     0,     0,
     0,     0,     0,     4,     5,     6,     0,    46,     0,     7,
     8,     0,     0,     0,     9,    10,    11,   182,     0,     0,
    14,   172,     0,     0,     0,     2,     3,     0,     0,     0,
    46,     0,   183,     4,     5,     6,     0,     0,     0,     7,
     8,     0,     0,     0,     9,    10,    11,   170,     0,     0,
    14,     0,     0,     0,     0,     2,     3,     0,     0,     0,
    46,     0,   131,     4,     5,     6,     0,     0,     0,     7,
     8,     0,     0,     0,     9,    10,    11,    12,   102,   103,
    14,   104,   105,   106,   107,   108,   109,     0,     0,     0,
    22,    79,    80,    81,    82,    83,    84,    85,    86,    87,
    88,    89,    90,    91,    92,    93,    94,    95,    96,    97,
    87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
    97,   249,    79,    80,    81,    82,    83,    84,    85,    86,
    87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
    97,    91,    92,    93,    94,    95,    96,    97,     0,     0,
     0,     0,   302,    79,    80,    81,    82,    83,    84,    85,
    86,    87,    88,    89,    90,    91,    92,    93,    94,    95,
    96,    97,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,   316,    79,    80,    81,    82,    83,    84,    85,
    86,    87,    88,    89,    90,    91,    92,    93,    94,    95,
    96,    97,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,   326,    79,    80,    81,    82,    83,    84,    85,
    86,    87,    88,    89,    90,    91,    92,    93,    94,    95,
    96,    97,    81,    82,    83,    84,    85,    86,    87,    88,
    89,    90,    91,    92,    93,    94,    95,    96,    97,    83,
    84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
    94,    95,    96,    97,    85,    86,    87,    88,    89,    90,
    91,    92,    93,    94,    95,    96,    97
};

static const short yycheck[] = {     7,
    29,    34,     5,    23,   108,   191,    34,     8,     5,     5,
   196,    13,     5,    21,     8,     9,    10,    11,    12,    13,
    14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
    24,    25,    26,     5,     5,    13,    42,    45,     5,    13,
    42,     5,     5,    72,     5,    13,   104,     5,   106,     5,
    58,     5,     5,    73,    45,    57,     5,     5,     5,     5,
    45,     5,    66,    66,    42,    73,    74,    75,    42,    66,
    66,    38,     5,    66,    42,    69,    37,    78,    67,    57,
    38,   267,    38,    57,    38,    38,     6,     5,    37,    57,
    38,    38,    38,    37,    66,    66,   125,    98,    99,   100,
   101,   121,    66,    66,    37,    22,    23,    24,    25,    26,
    42,   112,   120,   217,    34,   219,     5,     5,   107,    37,
     5,     5,     5,   124,    66,    57,     5,    42,    79,    80,
    81,    82,    83,    84,    85,    86,    87,    88,    89,    90,
    91,    92,    93,    94,    95,    96,    97,     6,    37,    37,
     6,   180,    37,    37,    37,    42,   107,   108,    37,     6,
   202,    42,    66,   114,    34,    35,    34,    35,    55,   211,
    57,   191,    20,    21,    55,    34,   196,    46,    55,    35,
    28,    29,    30,   134,   135,   214,    34,    35,    35,    35,
    34,    39,    40,    41,    42,    35,    42,    45,    34,   200,
    35,    35,    42,   192,   212,   213,   195,    55,    42,    57,
    34,   231,    34,    42,    42,    55,   205,    57,    35,    35,
   209,    69,    42,     5,    34,    42,    42,    42,   217,   237,
    38,   260,   233,    66,     6,   224,     6,     6,     6,    55,
   229,    57,   231,    42,    55,   234,     6,   267,     6,     4,
   201,    34,    35,   242,    38,    34,    39,    40,    41,    42,
    66,     5,    45,     5,   215,   216,   217,     6,   219,     6,
     6,   222,    55,   224,    38,   283,   284,     6,     5,   268,
   313,   270,    67,    34,    37,    66,     0,     0,   317,   318,
   298,   320,   293,   201,   136,   296,   247,   248,   280,   242,
    -1,   115,   291,   124,    39,   294,    -1,    -1,    -1,    -1,
    -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    -1,    28,
    29,    30,    -1,     3,    -1,    34,    35,    -1,   279,    -1,
    39,    40,    41,    42,   285,   286,    45,    -1,   289,    -1,
    20,    21,    -1,    -1,    -1,    -1,    55,    -1,    28,    29,
    30,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,
    40,    41,    42,    43,   315,    45,    -1,    47,    48,    -1,
    50,    51,    52,    53,    54,    55,    56,    57,     3,    -1,
    -1,    -1,    -1,    63,    -1,    -1,    66,    67,    68,    -1,
    -1,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    -1,
    -1,    -1,    -1,    28,    29,    30,    -1,    -1,    -1,    34,
    35,    -1,    -1,    -1,    39,    40,    41,    42,    43,     3,
    45,    -1,    47,    48,    -1,    50,    51,    52,    53,    54,
    55,    56,    57,    -1,    -1,    -1,    20,    21,    63,    -1,
    -1,    66,    67,    68,    28,    29,    30,    -1,    -1,    -1,
    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    43,
    -1,    45,    -1,    47,    48,    -1,    50,    51,    52,    53,
    54,    55,    56,    57,    -1,    -1,    -1,    -1,    -1,    63,
    -1,    -1,    66,    67,     8,     9,    10,    11,    12,    13,
    14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
    24,    25,    26,     6,    -1,     8,     9,    10,    11,    12,
    13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
    23,    24,    25,    26,    -1,    -1,    -1,    -1,    20,    21,
    -1,    -1,    -1,    -1,    -1,    -1,    28,    29,    30,    -1,
    -1,    -1,    34,    35,    -1,    69,    38,    39,    40,    41,
    42,    20,    21,    45,    -1,    58,    59,    60,    61,    28,
    29,    30,    -1,    55,    -1,    34,    35,    -1,    -1,    -1,
    39,    40,    41,    42,    -1,    -1,    45,    69,    -1,    -1,
    -1,    -1,    20,    21,    -1,    -1,    55,    -1,    57,    -1,
    28,    29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,
    69,    39,    40,    41,    42,    20,    21,    45,    -1,    -1,
    -1,    -1,    -1,    28,    29,    30,    -1,    55,    -1,    34,
    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
    45,    69,    -1,    -1,    -1,    20,    21,    -1,    -1,    -1,
    55,    -1,    57,    28,    29,    30,    -1,    -1,    -1,    34,
    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
    45,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    -1,
    55,    -1,    57,    28,    29,    30,    -1,    -1,    -1,    34,
    35,    -1,    -1,    -1,    39,    40,    41,    42,    28,    29,
    45,    31,    32,    33,    34,    35,    36,    -1,    -1,    -1,
    55,     8,     9,    10,    11,    12,    13,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
    16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
    26,    38,     8,     9,    10,    11,    12,    13,    14,    15,
    16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
    26,    20,    21,    22,    23,    24,    25,    26,    -1,    -1,
    -1,    -1,    38,     8,     9,    10,    11,    12,    13,    14,
    15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
    25,    26,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    37,     8,     9,    10,    11,    12,    13,    14,
    15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
    25,    26,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    37,     8,     9,    10,    11,    12,    13,    14,
    15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
    25,    26,    10,    11,    12,    13,    14,    15,    16,    17,
    18,    19,    20,    21,    22,    23,    24,    25,    26,    12,
    13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
    23,    24,    25,    26,    14,    15,    16,    17,    18,    19,
    20,    21,    22,    23,    24,    25,    26
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
#line 197 "/Users/hecht/work/FreeFem++/src/lg.y"
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
#line 229 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cinst=yyvsp[0].cexp;;;;
    break;}
case 4:
#line 230 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cinst= (yyvsp[-1].cinst+=yyvsp[0].cexp) ;
    break;}
case 5:
#line 233 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id=new ListOfId();;
    break;}
case 6:
#line 234 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str));
    break;}
case 7:
#line 235 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 8:
#line 236 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>()));
    break;}
case 9:
#line 237 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true));
    break;}
case 10:
#line 238 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 11:
#line 239 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 12:
#line 240 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 13:
#line 241 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = yyvsp[-2].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str)) ;
    break;}
case 14:
#line 242 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 15:
#line 243 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 16:
#line 244 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>())) ;
    break;}
case 17:
#line 245 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true)) ;
    break;}
case 18:
#line 246 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 19:
#line 247 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 20:
#line 250 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 21:
#line 251 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.clist_id=yyvsp[-2].clist_id  ; yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 24:
#line 256 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[0].str,dcltype);
    break;}
case 25:
#line 257 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-2].str,dcltype,yyvsp[0].cexp);
    break;}
case 26:
#line 258 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-3].str,dcltype,yyvsp[-1].args);
    break;}
case 27:
#line 259 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 28:
#line 266 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 29:
#line 267 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.args=Find(yyvsp[-1].str);
    break;}
case 30:
#line 268 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);
    break;}
case 31:
#line 269 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 32:
#line 270 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp));
    break;}
case 34:
#line 274 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.type=TypeArray(yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 35:
#line 275 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.type=TypeArray(yyvsp[-5].type,yyvsp[-3].type,yyvsp[-1].type);
    break;}
case 36:
#line 282 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[0].str,currentblock,fespacetype); ;
    break;}
case 37:
#line 283 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp); ;
    break;}
case 38:
#line 284 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-2].str,currentblock,fespacetype,yyvsp[0].cexp) ;
    break;}
case 39:
#line 285 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-1].clist_id,currentblock,fespacetype) ;
    break;}
case 40:
#line 286 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp) ;
    break;}
case 41:
#line 287 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp =  NewFEvariable(yyvsp[-3].clist_id,currentblock,fespacetype,yyvsp[0].cexp) ;
    break;}
case 42:
#line 290 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp); ;
    break;}
case 43:
#line 291 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp) ;
    break;}
case 44:
#line 295 "/Users/hecht/work/FreeFem++/src/lg.y"
{ fespacetype = Find(yyvsp[0].str);;
    break;}
case 45:
#line 298 "/Users/hecht/work/FreeFem++/src/lg.y"
{  yyval.cexp = yyvsp[0].cexp  ;
    break;}
case 46:
#line 299 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 47:
#line 301 "/Users/hecht/work/FreeFem++/src/lg.y"
{  yyval.cexp = yyvsp[0].cexp  ;
    break;}
case 48:
#line 302 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 49:
#line 304 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;
    break;}
case 50:
#line 305 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;
    break;}
case 51:
#line 310 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariableFES,size_t>(yyvsp[-3].str,atype<pfes*>(),yyvsp[-1].args,dimFESpaceImage(yyvsp[-1].args));
    break;}
case 53:
#line 313 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 54:
#line 316 "/Users/hecht/work/FreeFem++/src/lg.y"
{dcltype=yyvsp[0].type;
    break;}
case 55:
#line 316 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 56:
#line 317 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 57:
#line 318 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp=yyvsp[-1].cexp;
    break;}
case 58:
#line 319 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=currentblock->NewID(yyvsp[-4].type,yyvsp[-3].str,yyvsp[-1].cexp);;
    break;}
case 59:
#line 321 "/Users/hecht/work/FreeFem++/src/lg.y"
{   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = yyvsp[-4].type->right();
                      yyvsp[-1].routine=new Routine(yyvsp[-5].type,yyvsp[-4].type->right(),yyvsp[-3].str,yyvsp[-1].clist_id,currentblock);;
    break;}
case 60:
#line 326 "/Users/hecht/work/FreeFem++/src/lg.y"
{ currentblock=yyvsp[-5].routine->Set(yyvsp[-1].cinst);
                       currentblock->Add(yyvsp[-7].str,"(",yyvsp[-5].routine);
                       kkembtype--;
                       yyval.cexp=0 ;
    break;}
case 61:
#line 331 "/Users/hecht/work/FreeFem++/src/lg.y"
{currentblock = new Block(currentblock); yyvsp[-4].type->SetArgs(yyvsp[-1].clist_id);;
    break;}
case 62:
#line 333 "/Users/hecht/work/FreeFem++/src/lg.y"
{  yyval.cinst=currentblock->close(currentblock);
                         yyval.cexp=currentblock->NewID(yyvsp[-8].type,yyvsp[-7].str,yyvsp[-1].cexp,*yyvsp[-5].clist_id);
    break;}
case 63:
#line 337 "/Users/hecht/work/FreeFem++/src/lg.y"
{  currentblock = new Block(currentblock);
    break;}
case 64:
#line 338 "/Users/hecht/work/FreeFem++/src/lg.y"
{  yyval.cexp=currentblock->close(currentblock);
    break;}
case 65:
#line 340 "/Users/hecht/work/FreeFem++/src/lg.y"
{inloopcount++;;
    break;}
case 66:
#line 341 "/Users/hecht/work/FreeFem++/src/lg.y"
{inloopcount++;
    break;}
case 67:
#line 344 "/Users/hecht/work/FreeFem++/src/lg.y"
{dcltype=yyvsp[0].type;currentblock = new Block(currentblock);
    break;}
case 68:
#line 345 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 69:
#line 347 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=0;;
    break;}
case 70:
#line 348 "/Users/hecht/work/FreeFem++/src/lg.y"
{zzzfff->input(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 71:
#line 349 "/Users/hecht/work/FreeFem++/src/lg.y"
{load(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 72:
#line 350 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 73:
#line 351 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 74:
#line 352 "/Users/hecht/work/FreeFem++/src/lg.y"
{inloopcount--; yyval.cexp=For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 75:
#line 354 "/Users/hecht/work/FreeFem++/src/lg.y"
{inloopcount--; 
                yyval.cexp=C_F0(For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp),currentblock->close(currentblock));
    break;}
case 76:
#line 357 "/Users/hecht/work/FreeFem++/src/lg.y"
{inloopcount--;yyval.cexp=While(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 77:
#line 358 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=FIf(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 78:
#line 359 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=FIf(yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 79:
#line 360 "/Users/hecht/work/FreeFem++/src/lg.y"
{ 
                      yyval.cexp=C_F0(new E_block(yyvsp[-1].cinst,yyvsp[0].cexp),atype<void>()) ;
    break;}
case 80:
#line 362 "/Users/hecht/work/FreeFem++/src/lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-1].str,C_F0(TheOperators,"[border]",yyvsp[0].args));
    break;}
case 81:
#line 364 "/Users/hecht/work/FreeFem++/src/lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-4].str,C_F0(TheOperators,"[border]",yyvsp[-2].args));
    break;}
case 82:
#line 367 "/Users/hecht/work/FreeFem++/src/lg.y"
{
                    if(inloopcount) 
                      yyval.cexp= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;
    break;}
case 83:
#line 371 "/Users/hecht/work/FreeFem++/src/lg.y"
{ 
                    if(inloopcount)
                        yyval.cexp= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");
    break;}
case 84:
#line 375 "/Users/hecht/work/FreeFem++/src/lg.y"
{ 
                    if (kkembtype>=0)
                      yyval.cexp= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo(yyvsp[-1].cexp)) ,atype<void>());
                     else lgerror(" return not in routine ") ;
    break;}
case 85:
#line 383 "/Users/hecht/work/FreeFem++/src/lg.y"
{ 
   currentblock = new Block(currentblock);
   yyval.args = currentblock->NewVar<LocalVariable>(yyvsp[-5].str,atype<double*>());
   yyval.args+= yyvsp[-3].cexp;
   yyval.args+= yyvsp[-1].cexp ;
    break;}
case 86:
#line 389 "/Users/hecht/work/FreeFem++/src/lg.y"
{   
   yyval.args = (yyvsp[-1].args += yyvsp[0].cexp);
   currentblock->close(currentblock);
    break;}
case 88:
#line 396 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 95:
#line 410 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 96:
#line 411 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,"+=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 97:
#line 412 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,"-=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 98:
#line 413 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,"*=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 99:
#line 414 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,"/=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 101:
#line 419 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 102:
#line 420 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 103:
#line 421 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 104:
#line 422 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 105:
#line 423 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 106:
#line 424 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 107:
#line 425 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 108:
#line 426 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 109:
#line 427 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 110:
#line 428 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 111:
#line 429 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 112:
#line 430 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 113:
#line 431 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 114:
#line 432 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 115:
#line 433 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 116:
#line 434 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 117:
#line 435 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 118:
#line 436 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 119:
#line 437 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 120:
#line 442 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 121:
#line 443 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,":");
    break;}
case 122:
#line 444 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 123:
#line 445 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 124:
#line 448 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.args=0;
    break;}
case 125:
#line 449 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.args=Find(yyvsp[0].str);
    break;}
case 126:
#line 450 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);
    break;}
case 127:
#line 451 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 128:
#line 452 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.args = (yyvsp[-2].args += Find(yyvsp[0].str)) ;
    break;}
case 129:
#line 453 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 130:
#line 454 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 131:
#line 457 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 132:
#line 458 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 134:
#line 463 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[0].cexp);
    break;}
case 136:
#line 467 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 137:
#line 468 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 138:
#line 469 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 139:
#line 473 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=Find(yyvsp[0].str);;
    break;}
case 140:
#line 474 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp= CConstant(yyvsp[0].lnum);
    break;}
case 141:
#line 475 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp= CConstant(yyvsp[0].dnum);
    break;}
case 142:
#line 476 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp= CConstant(complex<double>(0,yyvsp[0].dnum));
    break;}
case 143:
#line 477 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp= CConstant<const char *>(yyvsp[0].str);
    break;}
case 144:
#line 478 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].args);;
    break;}
case 145:
#line 479 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].cexp);
    break;}
case 146:
#line 480 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(yyvsp[-5].cexp,yyvsp[-4].oper,yyvsp[-3].cexp,yyvsp[-1].cexp);
    break;}
case 147:
#line 481 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,"[]");
    break;}
case 148:
#line 482 "/Users/hecht/work/FreeFem++/src/lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].str) ;;
    break;}
case 149:
#line 483 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 150:
#line 484 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 151:
#line 485 "/Users/hecht/work/FreeFem++/src/lg.y"
{
             if (yyvsp[-3].type->right()->CastingFrom(yyvsp[-1].cexp.left()) ) 
                yyval.cexp=yyvsp[-3].type->right()->CastTo(yyvsp[-1].cexp)  ;
             else { yyval.cexp=C_F0(yyvsp[-1].cexp,yyvsp[-3].type->right()->name());
             if (!yyval.cexp.left()) { cerr << " no wait to change " << yyvsp[-1].cexp.left()->right()->name() << " in " << 
                                        yyvsp[-3].type->right()->name() << endl;
                                CompileError(" Error in type(exp) "); }
             }
            ;
    break;}
case 152:
#line 494 "/Users/hecht/work/FreeFem++/src/lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 153:
#line 495 "/Users/hecht/work/FreeFem++/src/lg.y"
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
#line 500 "/Users/hecht/work/FreeFem++/src/lg.y"
 


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
   
#ifdef EIGENVALUE
   init_eigenvalue();
#endif   
#ifdef PARALLELE
   init_lgparallele(); 
#endif 
#ifdef UMFPACK   
  cout << " UMFPACK ";  
#endif
 // callInitsFunct(); Pb opimisation 
   cout << endl;
  int ok;
  
  currentblock=0;
  currentblock = new Block(currentblock);  
  try {
    ok=yyparse(); //  compile
   if(ok==0)  
    if(currentblock) 
     cerr <<  "Error:a block is not close" << endl;   
    else 
     cerr <<  "Bien: On a fini Normalement" << endl; 
  }

  catch (Error & e) 
   {
     cerr << "error " << e.what() << endl;
   }
#ifdef PARALLELE
   end_parallele();
#endif
//  currentblock->close(currentblock).eval(thestack);
   fingraphique();
   return 0;
}


 
