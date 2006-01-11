/* A Bison parser, made from lg.y
   by GNU bison 1.35.  */

#define YYBISON 1  /* Identify Bison output.  */

#define yyparse lgparse
#define yylex lglex
#define yyerror lgerror
#define yylval lglval
#define yychar lgchar
#define yydebug lgdebug
#define yynerrs lgnerrs
# define	IF	257
# define	ELSE	258
# define	SET	259
# define	LTLT	260
# define	GTGT	261
# define	OR	262
# define	AND	263
# define	EQ	264
# define	NE	265
# define	LE	266
# define	GE	267
# define	DOTSTAR	268
# define	DOTSLASH	269
# define	UNARY	270
# define	PLUSPLUS	271
# define	MOINSMOINS	272
# define	LNUM	273
# define	DNUM	274
# define	CNUM	275
# define	ID	276
# define	FESPACEID	277
# define	IDPARAM	278
# define	STRING	279
# define	ENDOFFILE	280
# define	INCLUDE	281
# define	LOAD	282
# define	BIDON	283
# define	FOR	284
# define	WHILE	285
# define	BREAK	286
# define	CONTINUE	287
# define	RETURN	288
# define	TYPE	289
# define	FUNCTION	290
# define	FESPACE	291
# define	PLUSEQ	292
# define	MOINSEQ	293
# define	MULEQ	294
# define	DIVEQ	295
# define	ARROW	296
# define	BORDER	297
# define	CURVE	298
# define	SOLVE	299

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
#ifndef YYSTYPE
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
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
#ifndef YYDEBUG
# define YYDEBUG 1
#endif



#define	YYFINAL		347
#define	YYFLAG		-32768
#define	YYNTBASE	71

/* YYTRANSLATE(YYLEX) -- Bison token number corresponding to YYLEX. */
#define YYTRANSLATE(x) ((unsigned)(x) <= 299 ? yytranslate[x] : 112)

/* YYTRANSLATE[YYLEX] -- Bison token number corresponding to YYLEX. */
static const char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    30,     2,     2,     2,    24,    13,    32,
      34,    37,    22,    20,     5,    21,    36,    23,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    70,    66,
      16,     6,    17,    69,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    35,     2,    38,    31,    33,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    67,    11,    68,     2,     2,     2,     2,
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
       2,     2,     2,     2,     2,     2,     1,     3,     4,     7,
       8,     9,    10,    12,    14,    15,    18,    19,    25,    26,
      27,    28,    29,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
      56,    57,    58,    59,    60,    61,    62,    63,    64,    65
};

#if YYDEBUG
static const short yyprhs[] =
{
       0,     0,     3,     5,     7,    10,    11,    13,    17,    20,
      24,    27,    31,    35,    39,    45,    51,    56,    62,    67,
      73,    75,    79,    81,    83,    85,    89,    94,    98,   100,
     103,   107,   111,   117,   119,   124,   131,   136,   138,   143,
     147,   151,   158,   164,   169,   176,   178,   183,   185,   189,
     191,   195,   198,   204,   209,   211,   215,   216,   221,   225,
     228,   234,   235,   246,   247,   257,   259,   261,   263,   265,
     266,   270,   272,   275,   278,   281,   283,   293,   303,   309,
     315,   323,   327,   331,   338,   341,   344,   348,   356,   359,
     361,   365,   367,   369,   371,   373,   375,   377,   381,   385,
     389,   393,   397,   399,   405,   407,   411,   415,   419,   423,
     427,   431,   435,   439,   443,   447,   451,   455,   459,   463,
     467,   471,   475,   479,   483,   485,   487,   491,   497,   498,
     500,   504,   506,   510,   514,   520,   522,   526,   528,   531,
     533,   537,   541,   544,   546,   548,   550,   552,   554,   559,
     564,   571,   575,   579,   583,   588,   591,   594,   599,   603
};
static const short yyrhs[] =
{
      72,    46,     0,    73,     0,    98,     0,    73,    98,     0,
       0,    76,     0,    76,     6,   103,     0,    57,    76,     0,
      57,    13,    76,     0,    79,    76,     0,    79,    13,    76,
       0,    35,    74,    38,     0,    74,     5,    76,     0,    74,
       5,    35,    74,    38,     0,    74,     5,    76,     6,   103,
       0,    74,     5,    57,    76,     0,    74,     5,    57,    13,
      76,     0,    74,     5,    79,    76,     0,    74,     5,    79,
      13,    76,     0,    76,     0,    75,     5,    76,     0,    42,
       0,    57,     0,    42,     0,    42,     6,   103,     0,    42,
      34,    78,    37,     0,    77,     5,    77,     0,   104,     0,
      57,    42,     0,    42,     6,   104,     0,    78,     5,   104,
       0,    78,     5,    76,     6,   104,     0,    55,     0,    55,
      35,    55,    38,     0,    55,    35,    55,     5,    55,    38,
       0,    55,    16,    55,    17,     0,    42,     0,    42,    35,
     104,    38,     0,    42,     6,   104,     0,    35,    75,    38,
       0,    35,    75,    38,    35,   104,    38,     0,    35,    75,
      38,     6,   104,     0,    42,    34,   104,    37,     0,    35,
      75,    38,    34,   104,    37,     0,    57,     0,    57,    16,
      55,    17,     0,    81,     0,    83,     5,    81,     0,    80,
       0,    84,     5,    80,     0,    82,    84,     0,    82,    35,
      55,    38,    83,     0,    42,    34,    78,    37,     0,    86,
       0,    87,     5,    86,     0,     0,    79,    89,    77,    66,
       0,    43,    87,    66,     0,    85,    66,     0,    56,    42,
       6,   101,    66,     0,     0,    56,    79,    42,    34,    74,
      37,    90,    67,    73,    68,     0,     0,    56,    42,    34,
      74,    37,    91,     6,   103,    66,     0,    67,     0,    68,
       0,    50,     0,    51,     0,     0,    79,    97,    77,     0,
      66,     0,    47,    45,     0,    48,    45,     0,   101,    66,
       0,    88,     0,    94,    34,   101,    66,   101,    66,   101,
      37,    98,     0,    94,    34,    96,    66,   101,    66,   101,
      37,    98,     0,    95,    34,   101,    37,    98,     0,     3,
      34,   101,    37,    98,     0,     3,    34,   101,    37,    98,
       4,    98,     0,    92,    73,    93,     0,    63,    42,   100,
       0,    63,    42,    35,   108,    38,    66,     0,    52,    66,
       0,    53,    66,     0,    54,   101,    66,     0,    34,    42,
       6,   101,     5,   101,    37,     0,    99,    98,     0,   103,
       0,   101,     5,   101,     0,    21,     0,    20,     0,    30,
       0,    28,     0,    29,     0,   104,     0,   104,     6,   103,
       0,   104,    58,   103,     0,   104,    59,   103,     0,   104,
      60,   103,     0,   104,    61,   103,     0,   105,     0,   105,
      69,   104,    70,   104,     0,   109,     0,   105,    22,   105,
       0,   105,    25,   105,     0,   105,    26,   105,     0,   105,
      23,   105,     0,   105,    24,   105,     0,   105,    20,   105,
       0,   105,    21,   105,     0,   105,     8,   105,     0,   105,
       9,   105,     0,   105,    13,   105,     0,   105,    12,   105,
       0,   105,    11,   105,     0,   105,    10,   105,     0,   105,
      16,   105,     0,   105,    18,   105,     0,   105,    17,   105,
       0,   105,    19,   105,     0,   105,    14,   105,     0,   105,
      15,   105,     0,   104,     0,    70,     0,   104,    70,   104,
       0,   104,    70,   104,    70,   104,     0,     0,    57,     0,
      76,     6,   104,     0,   106,     0,   107,     5,    57,     0,
     107,     5,   106,     0,   107,     5,    76,     6,   104,     0,
     103,     0,   108,     5,   103,     0,   110,     0,   102,   110,
       0,   111,     0,   111,    31,   109,     0,   111,    33,   109,
       0,   111,    32,     0,    42,     0,    39,     0,    40,     0,
      41,     0,    45,     0,   111,    34,   107,    37,     0,   111,
      35,   106,    38,     0,   111,    35,   106,     5,   106,    38,
       0,   111,    35,    38,     0,   111,    36,    42,     0,    57,
      36,    42,     0,    57,    34,   107,    37,     0,   111,    28,
       0,   111,    29,     0,    55,    34,   101,    37,     0,    34,
     101,    37,     0,    35,   108,    38,     0
};

#endif

#if YYDEBUG
/* YYRLINE[YYN] -- source line where rule number YYN was defined. */
static const short yyrline[] =
{
       0,   196,   232,   235,   236,   239,   240,   241,   242,   243,
     244,   245,   246,   247,   248,   249,   250,   251,   252,   253,
     256,   257,   260,   260,   262,   263,   264,   266,   272,   274,
     275,   276,   277,   280,   281,   282,   283,   289,   291,   292,
     293,   294,   295,   297,   299,   303,   304,   309,   310,   312,
     313,   315,   316,   319,   324,   325,   328,   328,   329,   330,
     331,   332,   332,   347,   347,   356,   357,   359,   361,   365,
     365,   369,   370,   371,   372,   373,   374,   375,   379,   380,
     381,   382,   384,   386,   389,   393,   397,   405,   411,   416,
     418,   422,   424,   425,   426,   427,   430,   432,   433,   434,
     435,   436,   440,   442,   444,   446,   447,   448,   449,   450,
     451,   452,   453,   454,   455,   456,   457,   458,   459,   460,
     461,   462,   463,   464,   468,   470,   471,   472,   475,   476,
     477,   478,   479,   480,   481,   484,   485,   488,   490,   493,
     494,   495,   496,   499,   501,   502,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   523,   524
};
#endif


#if (YYDEBUG) || defined YYERROR_VERBOSE

/* YYTNAME[TOKEN_NUM] -- String name of the token TOKEN_NUM. */
static const char *const yytname[] =
{
  "$", "error", "$undefined.", "IF", "ELSE", "','", "'='", "SET", "LTLT", 
  "GTGT", "OR", "'|'", "AND", "'&'", "EQ", "NE", "'<'", "'>'", "LE", "GE", 
  "'+'", "'-'", "'*'", "'/'", "'%'", "DOTSTAR", "DOTSLASH", "UNARY", 
  "PLUSPLUS", "MOINSMOINS", "'!'", "'^'", "'\\''", "'_'", "'('", "'['", 
  "'.'", "')'", "']'", "LNUM", "DNUM", "CNUM", "ID", "FESPACEID", 
  "IDPARAM", "STRING", "ENDOFFILE", "INCLUDE", "LOAD", "BIDON", "FOR", 
  "WHILE", "BREAK", "CONTINUE", "RETURN", "TYPE", "FUNCTION", "FESPACE", 
  "PLUSEQ", "MOINSEQ", "MULEQ", "DIVEQ", "ARROW", "BORDER", "CURVE", 
  "SOLVE", "';'", "'{'", "'}'", "'?'", "':'", "start", "input", 
  "instructions", "list_of_id_args", "list_of_id1", "id", "list_of_dcls", 
  "parameters_list", "type_of_dcl", "ID_space", "ID_array_space", 
  "fespace", "spaceIDa", "spaceIDb", "spaceIDs", "fespace_def", 
  "fespace_def_list", "declaration", "@1", "@2", "@3", "begin", "end", 
  "for_loop", "while_loop", "declaration_for", "@4", "instruction", 
  "bornes", "border_expr", "Expr", "unop", "no_comma_expr", "no_set_expr", 
  "no_ternary_expr", "sub_script_expr", "parameters", "array", 
  "unary_expr", "pow_expr", "primary", 0
};
#endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives. */
static const short yyr1[] =
{
       0,    71,    72,    73,    73,    74,    74,    74,    74,    74,
      74,    74,    74,    74,    74,    74,    74,    74,    74,    74,
      75,    75,    76,    76,    77,    77,    77,    77,    78,    78,
      78,    78,    78,    79,    79,    79,    79,    80,    80,    80,
      80,    80,    80,    81,    81,    82,    82,    83,    83,    84,
      84,    85,    85,    86,    87,    87,    89,    88,    88,    88,
      88,    90,    88,    91,    88,    92,    93,    94,    95,    97,
      96,    98,    98,    98,    98,    98,    98,    98,    98,    98,
      98,    98,    98,    98,    98,    98,    98,    99,   100,   101,
     101,   102,   102,   102,   102,   102,   103,   103,   103,   103,
     103,   103,   104,   104,   105,   105,   105,   105,   105,   105,
     105,   105,   105,   105,   105,   105,   105,   105,   105,   105,
     105,   105,   105,   105,   106,   106,   106,   106,   107,   107,
     107,   107,   107,   107,   107,   108,   108,   109,   109,   110,
     110,   110,   110,   111,   111,   111,   111,   111,   111,   111,
     111,   111,   111,   111,   111,   111,   111,   111,   111,   111
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN. */
static const short yyr2[] =
{
       0,     2,     1,     1,     2,     0,     1,     3,     2,     3,
       2,     3,     3,     3,     5,     5,     4,     5,     4,     5,
       1,     3,     1,     1,     1,     3,     4,     3,     1,     2,
       3,     3,     5,     1,     4,     6,     4,     1,     4,     3,
       3,     6,     5,     4,     6,     1,     4,     1,     3,     1,
       3,     2,     5,     4,     1,     3,     0,     4,     3,     2,
       5,     0,    10,     0,     9,     1,     1,     1,     1,     0,
       3,     1,     2,     2,     2,     1,     9,     9,     5,     5,
       7,     3,     3,     6,     2,     2,     3,     7,     2,     1,
       3,     1,     1,     1,     1,     1,     1,     3,     3,     3,
       3,     3,     1,     5,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     1,     1,     3,     5,     0,     1,
       3,     1,     3,     3,     5,     1,     3,     1,     2,     1,
       3,     3,     2,     1,     1,     1,     1,     1,     4,     4,
       6,     3,     3,     3,     4,     2,     2,     4,     3,     3
};

/* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
   doesn't specify something else to do.  Zero means the default is an
   error. */
static const short yydefact[] =
{
       0,     0,    92,    91,    94,    95,    93,     0,     0,   144,
     145,   146,   143,     0,   147,     0,     0,    67,    68,     0,
       0,     0,    33,     0,    45,     0,    71,    65,     0,     2,
      56,     0,     0,    75,     0,     0,     0,     3,     0,     0,
      89,    96,   102,   104,   137,   139,     0,     0,     0,     0,
     135,     0,     0,    54,     0,    72,    73,    84,    85,     0,
       0,     0,     0,     0,    33,     0,     0,   128,     0,     0,
       1,     4,     0,     0,    37,    49,    51,    59,     0,     0,
       0,     0,    74,   138,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   155,
     156,     0,   142,     0,   128,     0,     0,     0,   158,     0,
     159,     0,     0,    58,    86,     0,     0,     0,     0,     5,
       0,     0,   143,   129,   125,     0,   124,   131,     0,   153,
       0,     0,     0,    82,    24,     0,    22,     0,    23,     0,
      20,     0,     0,     0,    66,    81,    69,     0,     0,     0,
      90,    97,    98,    99,   100,   101,   112,   113,   117,   116,
     115,   114,   122,   123,   118,   120,   119,   121,   110,   111,
     105,   108,   109,   106,   107,     0,   140,   141,     0,   151,
       0,   152,     0,   136,   143,     0,     0,    28,    55,    36,
     157,     0,    34,     0,     5,    23,     0,     6,     0,     5,
      46,     0,     0,     0,   154,     0,     0,    88,     0,     0,
       0,    57,     0,     0,    40,    39,     0,     0,    50,     0,
       0,     0,     0,     0,   148,     0,   149,    79,     0,    29,
       0,    53,     0,    60,     0,     0,     8,     0,    63,     0,
       0,    10,     0,   130,   126,   132,     0,   133,     0,     0,
      25,     0,    27,     0,     0,    47,    52,    21,     0,     0,
      38,    70,     0,     0,    78,   103,     0,     0,    30,    23,
       0,    31,    35,    12,     9,     5,    23,    13,     0,     0,
       7,    11,    61,     0,     0,     0,    83,    26,     0,     0,
       0,    42,     0,     0,     0,   150,    80,     0,     0,     0,
      16,     0,     0,    18,     0,     0,   127,   134,     0,     0,
       0,    48,    41,     0,     0,    32,    14,    17,    15,    19,
       0,     0,    90,     0,    43,     0,     0,    64,     0,    87,
       0,    77,    76,    62,    44,     0,     0,     0
};

static const short yydefgoto[] =
{
     345,    28,    29,   206,   149,   207,   145,   196,    30,    75,
     265,    31,   266,    76,    32,    53,    54,    33,    72,   315,
     289,    34,   155,    35,    36,   157,   229,    37,   142,   143,
      38,    39,    40,    41,    42,   137,   138,    51,    43,    44,
      45
};

static const short yypact[] =
{
     367,   -17,-32768,-32768,-32768,-32768,-32768,   503,   503,-32768,
  -32768,-32768,-32768,     6,-32768,     8,    74,-32768,-32768,    71,
      91,   503,   167,   148,   163,    66,-32768,-32768,    76,   367,
  -32768,    14,    96,-32768,   367,   158,   173,-32768,     2,   180,
  -32768,     5,    12,-32768,-32768,   247,   503,   178,   205,     7,
  -32768,    42,   183,-32768,     3,-32768,-32768,-32768,-32768,     4,
     103,   503,   127,   146,   100,   121,   174,   232,   191,   133,
  -32768,-32768,   200,   171,   112,-32768,   238,-32768,   285,   533,
     503,   503,-32768,-32768,   503,   503,   503,   503,   503,   503,
     503,   503,   503,   503,   503,   503,   503,   503,   503,   503,
     503,   503,   503,   503,   503,   503,   503,   503,   503,-32768,
  -32768,   503,-32768,   503,   232,   407,   203,    50,-32768,   503,
  -32768,   563,     6,-32768,-32768,   231,    73,    56,   503,   143,
     216,   237,   250,   140,-32768,   251,   189,-32768,    87,-32768,
     221,   503,   367,-32768,   147,    11,-32768,   226,-32768,    57,
  -32768,   503,   503,   152,-32768,-32768,-32768,   202,    13,   101,
  -32768,-32768,-32768,-32768,-32768,-32768,   660,   660,   675,   675,
     688,   688,   699,   699,   619,   619,   619,   619,   268,   268,
  -32768,-32768,-32768,-32768,-32768,   199,-32768,-32768,   107,-32768,
      88,-32768,   367,-32768,   278,   204,   108,-32768,-32768,-32768,
  -32768,   230,-32768,    35,   143,    33,   118,   280,    54,   143,
  -32768,   503,   503,   439,-32768,   289,    92,-32768,   503,   563,
     200,-32768,   176,     0,   136,-32768,   260,     0,-32768,   200,
     503,   503,   367,   503,-32768,   471,-32768,   296,   503,-32768,
     593,-32768,   263,-32768,    93,     0,-32768,   153,-32768,   503,
       0,-32768,   123,-32768,   233,   140,   298,-32768,   503,   252,
  -32768,   124,-32768,     0,   282,-32768,   312,-32768,   503,   503,
  -32768,   316,    36,    38,-32768,-32768,   284,   367,-32768,   205,
     317,-32768,-32768,-32768,-32768,   143,    58,   325,    63,   328,
  -32768,-32768,-32768,   503,   503,   338,-32768,-32768,    94,   503,
     176,-32768,   306,   503,   503,-32768,-32768,   503,   102,     0,
  -32768,   503,     0,-32768,   503,   283,-32768,-32768,   503,   311,
     320,-32768,-32768,   128,   129,-32768,-32768,-32768,-32768,-32768,
     292,   367,   322,   503,-32768,   367,   367,-32768,   326,-32768,
     327,-32768,-32768,-32768,-32768,   349,   362,-32768
};

static const short yypgoto[] =
{
  -32768,-32768,   -32,  -194,   109,   -54,   -93,   144,   -20,   222,
      84,-32768,-32768,-32768,-32768,   264,-32768,-32768,-32768,-32768,
  -32768,-32768,-32768,-32768,-32768,-32768,-32768,   -28,-32768,-32768,
      -7,-32768,    -2,   -63,   562,  -110,   271,   249,    30,   352,
  -32768
};


#define	YYLAST		725


static const short yytable[] =
{
      49,    71,    78,    65,   136,   190,    50,    81,   122,    81,
     244,    84,    81,   135,    59,   252,   220,    46,    81,   150,
      89,    90,    91,    92,    93,    94,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   117,
      81,    81,   146,    81,   118,   185,   245,   119,    52,    73,
      71,   136,   136,    55,   126,    81,    74,   148,   197,   156,
     135,   201,   223,    85,    86,    87,    88,   250,    82,   123,
     124,   309,   158,   159,   160,   146,   312,   221,    81,   231,
     120,   108,   161,   162,   163,   164,   165,   192,   225,   226,
     148,   308,   213,   235,   202,   224,   146,   119,   247,   223,
     146,   243,   303,   257,   304,   146,    81,   247,    69,   208,
     200,   148,   213,   240,   217,   148,    60,   193,   151,    56,
     148,   203,    70,   247,   214,   276,   236,   262,   247,   240,
     259,   283,   319,    81,    81,    62,   271,    57,   232,    50,
     326,   186,   268,   187,   234,   241,   -23,   152,   253,   254,
     136,   246,   128,   218,   251,   248,   197,    58,   125,   256,
     292,   297,    77,   130,   237,   335,   336,   140,   141,   267,
     275,   269,   136,   150,    67,   278,    68,   281,   204,    66,
     129,   219,   127,    60,   208,   146,   280,   227,   285,   208,
      63,   284,    79,   287,    74,   146,   291,    67,    64,    68,
     205,    61,    62,    64,   274,   301,   302,    80,    64,   150,
     286,   263,    61,   146,     7,     8,   260,   121,   264,     9,
      10,    11,    12,   272,   273,    14,   147,   288,   148,   131,
     316,   317,   310,   139,   313,    47,   320,    48,    67,    67,
      68,    68,   144,   153,   325,   191,   239,   290,   199,   306,
     209,   295,     2,     3,   210,   327,   -22,   211,   329,   212,
       4,     5,     6,   215,   222,   208,     7,     8,   230,   233,
     340,     9,    10,    11,   132,   109,   110,    14,   111,   112,
     113,   114,   115,   116,   238,   242,   249,    47,     1,   133,
     103,   104,   105,   106,   107,   258,   323,   324,   270,   338,
     277,   282,   134,   293,   294,     2,     3,   341,   342,   328,
      71,   332,   330,     4,     5,     6,   299,   300,   296,     7,
       8,   220,   305,   307,     9,    10,    11,    12,    13,     1,
      14,   311,    15,    16,   314,    17,    18,    19,    20,    21,
      22,    23,    24,   318,   322,   333,     2,     3,    25,   346,
     331,    26,    27,   154,     4,     5,     6,   334,   337,   339,
       7,     8,   347,   261,   344,     9,    10,    11,    12,    13,
       1,    14,   298,    15,    16,   228,    17,    18,    19,    20,
      21,    22,    23,    24,   321,   188,   198,     2,     3,    25,
     216,    83,    26,    27,   343,     4,     5,     6,     0,     0,
       0,     7,     8,     0,     0,     0,     9,    10,    11,    12,
      13,     0,    14,     0,    15,    16,     0,    17,    18,    19,
      20,    21,    22,    23,    24,     0,     0,     2,     3,     0,
      25,     0,     0,    26,    27,     4,     5,     6,     0,     0,
       0,     7,     8,     0,     0,   189,     9,    10,    11,    12,
       0,     0,    14,     0,     0,     0,     0,     0,     0,     2,
       3,     0,    47,     0,    48,     0,     0,     4,     5,     6,
       0,     0,     0,     7,     8,     0,     0,   134,     9,    10,
      11,   132,     0,     0,    14,     0,     0,     0,     0,     0,
       0,     2,     3,     0,    47,     0,   255,     0,     0,     4,
       5,     6,     0,     0,     0,     7,     8,     0,     0,   134,
       9,    10,    11,    12,     0,     0,    14,     0,     0,     0,
       0,     0,     0,     2,     3,     0,    47,     0,    48,     0,
       0,     4,     5,     6,     0,     0,     0,     7,     8,     0,
       0,   134,     9,    10,    11,    12,     0,     0,    14,     0,
       0,     0,     0,     2,     3,     0,     0,     0,    47,     0,
      48,     4,     5,     6,     0,     0,     0,     7,     8,     0,
       0,     0,     9,    10,    11,    12,     0,     0,    14,     0,
       0,     0,     0,     2,     3,     0,     0,     0,    22,     0,
      48,     4,     5,     6,     0,     0,     0,     7,     8,     0,
       0,     0,     9,    10,    11,   194,     0,     0,    14,     0,
       0,     0,     0,     2,     3,     0,     0,     0,    47,     0,
     195,     4,     5,     6,     0,     0,     0,     7,     8,     0,
       0,     0,     9,    10,    11,   132,     0,     0,    14,   101,
     102,   103,   104,   105,   106,   107,     0,     0,    47,     0,
     279,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
      91,    92,    93,    94,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,   105,   106,   107,    93,    94,    95,
      96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,    95,    96,    97,    98,    99,   100,   101,   102,
     103,   104,   105,   106,   107,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   107
};

static const short yycheck[] =
{
       7,    29,    34,    23,    67,   115,     8,     5,     5,     5,
     204,     6,     5,    67,    21,   209,     5,    34,     5,    73,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    46,
       5,     5,    42,     5,    37,   108,    13,     5,    42,    35,
      78,   114,   115,    45,    61,     5,    42,    57,   121,    79,
     114,     5,     5,    58,    59,    60,    61,    13,    66,    66,
      66,    13,    79,    80,    81,    42,    13,    66,     5,    66,
      38,    69,    84,    85,    86,    87,    88,    37,   151,   152,
      57,   285,     5,     5,    38,    38,    42,     5,     5,     5,
      42,    66,    66,   213,    66,    42,     5,     5,    42,   129,
      37,    57,     5,     5,   142,    57,    16,   119,     6,    45,
      57,   128,    46,     5,    37,   235,    38,   220,     5,     5,
      38,    38,    38,     5,     5,    35,   229,    66,    37,   141,
      38,   111,     6,   113,    37,    37,     6,    35,   211,   212,
     213,   205,     6,     6,   208,    37,   219,    66,    55,   213,
      37,    37,    66,    42,   192,    37,    37,    34,    35,   223,
     233,    35,   235,   227,    34,   238,    36,   240,    35,    16,
      34,    34,    55,    16,   204,    42,   240,    35,    35,   209,
      42,   245,    34,   247,    42,    42,   250,    34,    55,    36,
      57,    34,    35,    55,   232,   268,   269,    34,    55,   263,
      57,    35,    34,    42,    34,    35,   218,    34,    42,    39,
      40,    41,    42,   230,   231,    45,    55,   247,    57,    55,
     293,   294,   286,    42,   288,    55,   299,    57,    34,    34,
      36,    36,    42,     5,   307,    42,    42,   249,    17,   277,
      34,   258,    20,    21,    17,   309,     6,     6,   312,    70,
      28,    29,    30,    42,    38,   285,    34,    35,    66,    70,
     333,    39,    40,    41,    42,    28,    29,    45,    31,    32,
      33,    34,    35,    36,     6,    55,     6,    55,     3,    57,
      22,    23,    24,    25,    26,     6,   303,   304,    38,   331,
       4,    38,    70,    70,     6,    20,    21,   335,   336,   311,
     338,   318,   314,    28,    29,    30,    34,     5,    66,    34,
      35,     5,    38,     6,    39,    40,    41,    42,    43,     3,
      45,     6,    47,    48,     6,    50,    51,    52,    53,    54,
      55,    56,    57,     5,    38,    34,    20,    21,    63,     0,
      67,    66,    67,    68,    28,    29,    30,    37,    66,    37,
      34,    35,     0,   219,    37,    39,    40,    41,    42,    43,
       3,    45,   263,    47,    48,   153,    50,    51,    52,    53,
      54,    55,    56,    57,   300,   114,   122,    20,    21,    63,
     141,    39,    66,    67,    68,    28,    29,    30,    -1,    -1,
      -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      43,    -1,    45,    -1,    47,    48,    -1,    50,    51,    52,
      53,    54,    55,    56,    57,    -1,    -1,    20,    21,    -1,
      63,    -1,    -1,    66,    67,    28,    29,    30,    -1,    -1,
      -1,    34,    35,    -1,    -1,    38,    39,    40,    41,    42,
      -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    20,
      21,    -1,    55,    -1,    57,    -1,    -1,    28,    29,    30,
      -1,    -1,    -1,    34,    35,    -1,    -1,    70,    39,    40,
      41,    42,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    20,    21,    -1,    55,    -1,    57,    -1,    -1,    28,
      29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,    70,
      39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    20,    21,    -1,    55,    -1,    57,    -1,
      -1,    28,    29,    30,    -1,    -1,    -1,    34,    35,    -1,
      -1,    70,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    55,    -1,
      57,    28,    29,    30,    -1,    -1,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    55,    -1,
      57,    28,    29,    30,    -1,    -1,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    55,    -1,
      57,    28,    29,    30,    -1,    -1,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,    20,
      21,    22,    23,    24,    25,    26,    -1,    -1,    55,    -1,
      57,    89,    90,    91,    92,    93,    94,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,   107,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/sw/share/bison/bison.simple"

/* Skeleton output parser for bison,

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software
   Foundation, Inc.

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

/* This is the parser code that is written into each bison parser when
   the %semantic_parser declaration is not specified in the grammar.
   It was written by Richard Stallman by simplifying the hairy parser
   used when %semantic_parser is specified.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

#if ! defined (yyoverflow) || defined (YYERROR_VERBOSE)

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || defined (YYERROR_VERBOSE) */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (YYLTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
# if YYLSP_NEEDED
  YYLTYPE yyls;
# endif
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAX (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# if YYLSP_NEEDED
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE) + sizeof (YYLTYPE))	\
      + 2 * YYSTACK_GAP_MAX)
# else
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAX)
# endif

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAX;	\
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif


#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0
#define YYACCEPT	goto yyacceptlab
#define YYABORT 	goto yyabortlab
#define YYERROR		goto yyerrlab1
/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL		goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");			\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).

   When YYLLOC_DEFAULT is run, CURRENT is set the location of the
   first token.  By default, to implement support for ranges, extend
   its range to the last symbol.  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)       	\
   Current.last_line   = Rhs[N].last_line;	\
   Current.last_column = Rhs[N].last_column;
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#if YYPURE
# if YYLSP_NEEDED
#  ifdef YYLEX_PARAM
#   define YYLEX		yylex (&yylval, &yylloc, YYLEX_PARAM)
#  else
#   define YYLEX		yylex (&yylval, &yylloc)
#  endif
# else /* !YYLSP_NEEDED */
#  ifdef YYLEX_PARAM
#   define YYLEX		yylex (&yylval, YYLEX_PARAM)
#  else
#   define YYLEX		yylex (&yylval)
#  endif
# endif /* !YYLSP_NEEDED */
#else /* !YYPURE */
# define YYLEX			yylex ()
#endif /* !YYPURE */


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)
/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
#endif /* !YYDEBUG */

/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif

#ifdef YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif
#endif

#line 315 "/sw/share/bison/bison.simple"


/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
#  define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL
# else
#  define YYPARSE_PARAM_ARG YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
# endif
#else /* !YYPARSE_PARAM */
# define YYPARSE_PARAM_ARG
# define YYPARSE_PARAM_DECL
#endif /* !YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
# ifdef YYPARSE_PARAM
int yyparse (void *);
# else
int yyparse (void);
# endif
#endif

/* YY_DECL_VARIABLES -- depending whether we use a pure parser,
   variables are global, or local to YYPARSE.  */

#define YY_DECL_NON_LSP_VARIABLES			\
/* The lookahead symbol.  */				\
int yychar;						\
							\
/* The semantic value of the lookahead symbol. */	\
YYSTYPE yylval;						\
							\
/* Number of parse errors so far.  */			\
int yynerrs;

#if YYLSP_NEEDED
# define YY_DECL_VARIABLES			\
YY_DECL_NON_LSP_VARIABLES			\
						\
/* Location data for the lookahead symbol.  */	\
YYLTYPE yylloc;
#else
# define YY_DECL_VARIABLES			\
YY_DECL_NON_LSP_VARIABLES
#endif


/* If nonreentrant, generate the variables here. */

#if !YYPURE
YY_DECL_VARIABLES
#endif  /* !YYPURE */

int
yyparse (YYPARSE_PARAM_ARG)
     YYPARSE_PARAM_DECL
{
  /* If reentrant, generate the variables here. */
#if YYPURE
  YY_DECL_VARIABLES
#endif  /* !YYPURE */

  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yychar1 = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack. */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;

#if YYLSP_NEEDED
  /* The location stack.  */
  YYLTYPE yylsa[YYINITDEPTH];
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;
#endif

#if YYLSP_NEEDED
# define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
# define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  YYSIZE_T yystacksize = YYINITDEPTH;


  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
#if YYLSP_NEEDED
  YYLTYPE yyloc;
#endif

  /* When reducing, the number of symbols on the RHS of the reduced
     rule. */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;
#if YYLSP_NEEDED
  yylsp = yyls;
#endif
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  */
# if YYLSP_NEEDED
	YYLTYPE *yyls1 = yyls;
	/* This used to be a conditional around just the two extra args,
	   but that might be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yyls1, yysize * sizeof (*yylsp),
		    &yystacksize);
	yyls = yyls1;
# else
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);
# endif
	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);
# if YYLSP_NEEDED
	YYSTACK_RELOCATE (yyls);
# endif
# undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
#if YYLSP_NEEDED
      yylsp = yyls + yysize - 1;
#endif

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
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
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0)		/* This means end of input. */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more */

      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yychar1 = YYTRANSLATE (yychar);

#if YYDEBUG
     /* We have to keep this `#if YYDEBUG', since we use variables
	which are defined only if `YYDEBUG' is set.  */
      if (yydebug)
	{
	  YYFPRINTF (stderr, "Next token is %d (%s",
		     yychar, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise
	     meaning of a token, for further debugging info.  */
# ifdef YYPRINT
	  YYPRINT (stderr, yychar, yylval);
# endif
	  YYFPRINTF (stderr, ")\n");
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
  YYDPRINTF ((stderr, "Shifting token %d (%s), ",
	      yychar, yytname[yychar1]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#if YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to the semantic value of
     the lookahead token.  This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

#if YYLSP_NEEDED
  /* Similarly for the default location.  Let the user run additional
     commands if for instance locations are ranges.  */
  yyloc = yylsp[1-yylen];
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
#endif

#if YYDEBUG
  /* We have to keep this `#if YYDEBUG', since we use variables which
     are defined only if `YYDEBUG' is set.  */
  if (yydebug)
    {
      int yyi;

      YYFPRINTF (stderr, "Reducing via rule %d (line %d), ",
		 yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (yyi = yyprhs[yyn]; yyrhs[yyi] > 0; yyi++)
	YYFPRINTF (stderr, "%s ", yytname[yyrhs[yyi]]);
      YYFPRINTF (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif

  switch (yyn) {

case 1:
#line 196 "lg.y"
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
#line 235 "lg.y"
{yyval.cinst=yyvsp[0].cexp;;;;
    break;}
case 4:
#line 236 "lg.y"
{ yyval.cinst= (yyvsp[-1].cinst+=yyvsp[0].cexp) ;
    break;}
case 5:
#line 239 "lg.y"
{ yyval.clist_id=new ListOfId();;
    break;}
case 6:
#line 240 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str));
    break;}
case 7:
#line 241 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 8:
#line 242 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>()));
    break;}
case 9:
#line 243 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true));
    break;}
case 10:
#line 244 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 11:
#line 245 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 12:
#line 246 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 13:
#line 247 "lg.y"
{ yyval.clist_id = yyvsp[-2].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str)) ;
    break;}
case 14:
#line 248 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;
    break;}
case 15:
#line 249 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 16:
#line 250 "lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>())) ;
    break;}
case 17:
#line 251 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true)) ;
    break;}
case 18:
#line 252 "lg.y"
{ yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;
    break;}
case 19:
#line 253 "lg.y"
{ yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;
    break;}
case 20:
#line 256 "lg.y"
{ yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 21:
#line 257 "lg.y"
{ yyval.clist_id=yyvsp[-2].clist_id  ; yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;
    break;}
case 24:
#line 262 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[0].str,dcltype);
    break;}
case 25:
#line 263 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-2].str,dcltype,yyvsp[0].cexp);
    break;}
case 26:
#line 264 "lg.y"
{yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-3].str,dcltype,yyvsp[-1].args);
                                              yyvsp[-1].args.destroy();
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
     yyvsp[-1].args.destroy(); ;
    break;}
case 55:
#line 325 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 56:
#line 328 "lg.y"
{dcltype=yyvsp[0].type;
    break;}
case 57:
#line 328 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 58:
#line 329 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 59:
#line 330 "lg.y"
{ yyval.cexp=yyvsp[-1].cexp;
    break;}
case 60:
#line 331 "lg.y"
{yyval.cexp=currentblock->NewID(yyvsp[-4].type,yyvsp[-3].str,yyvsp[-1].cexp);;
    break;}
case 61:
#line 333 "lg.y"
{   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = yyvsp[-4].type->right();
                      routineinblock[kkembtype] = currentblock;
                      yyvsp[-1].routine=new Routine(yyvsp[-5].type,yyvsp[-4].type->right(),yyvsp[-3].str,yyvsp[-1].clist_id,currentblock);
                     // cout << " \n after new routine \n " << endl;                      
                      ;
    break;}
case 62:
#line 341 "lg.y"
{ currentblock=yyvsp[-5].routine->Set(yyvsp[-1].cinst);
                       currentblock->Add(yyvsp[-7].str,"(",yyvsp[-5].routine);
                       kkembtype--;
                       yyval.cexp=0;
                    
                        ;
    break;}
case 63:
#line 348 "lg.y"
{currentblock = new Block(currentblock); yyvsp[-4].type->SetArgs(yyvsp[-1].clist_id);;
    break;}
case 64:
#line 350 "lg.y"
{  yyval.cinst=currentblock->close(currentblock);
                         yyval.cexp=currentblock->NewID(yyvsp[-8].type,yyvsp[-7].str,yyvsp[-1].cexp,*yyvsp[-5].clist_id);
                         delete yyvsp[-5].clist_id; //  FH 23032005
                         ;
    break;}
case 65:
#line 356 "lg.y"
{  currentblock = new Block(currentblock);
    break;}
case 66:
#line 357 "lg.y"
{  yyval.cexp=currentblock->close(currentblock);
    break;}
case 67:
#line 359 "lg.y"
{ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;
    break;}
case 68:
#line 361 "lg.y"
{ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;
    break;}
case 69:
#line 366 "lg.y"
{dcltype=yyvsp[0].type;currentblock = new Block(currentblock);
    break;}
case 70:
#line 367 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 71:
#line 369 "lg.y"
{yyval.cexp=0;;
    break;}
case 72:
#line 370 "lg.y"
{zzzfff->input(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 73:
#line 371 "lg.y"
{load(yyvsp[0].str);yyval.cexp= 0; ;
    break;}
case 74:
#line 372 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 75:
#line 373 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 76:
#line 374 "lg.y"
{inloopcount--; yyval.cexp=For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 77:
#line 376 "lg.y"
{inloopcount--; 
                yyval.cexp=C_F0(For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp),currentblock->close(currentblock));
    break;}
case 78:
#line 379 "lg.y"
{inloopcount--;yyval.cexp=While(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 79:
#line 380 "lg.y"
{yyval.cexp=FIf(yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 80:
#line 381 "lg.y"
{yyval.cexp=FIf(yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 81:
#line 382 "lg.y"
{ 
                      yyval.cexp=C_F0(new E_block(yyvsp[-1].cinst,yyvsp[0].cexp),atype<void>()) ;
    break;}
case 82:
#line 384 "lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-1].str,C_F0(TheOperators,"[border]",yyvsp[0].args));
    break;}
case 83:
#line 386 "lg.y"
{
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-4].str,C_F0(TheOperators,"[border]",yyvsp[-2].args));
    break;}
case 84:
#line 389 "lg.y"
{
                    if(inloopcount) 
                      yyval.cexp= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;
    break;}
case 85:
#line 393 "lg.y"
{ 
                    if(inloopcount)
                        yyval.cexp= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");
    break;}
case 86:
#line 397 "lg.y"
{ 
                    if (kkembtype>=0)
                      yyval.cexp= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo(yyvsp[-1].cexp)) ,atype<void>());
                     else lgerror(" return not in routine ") ;
    break;}
case 87:
#line 405 "lg.y"
{ 
   currentblock = new Block(currentblock);
   yyval.args = currentblock->NewVar<LocalVariable>(yyvsp[-5].str,atype<double*>());
   yyval.args+= yyvsp[-3].cexp;
   yyval.args+= yyvsp[-1].cexp ;
    break;}
case 88:
#line 411 "lg.y"
{   
   yyval.args = (yyvsp[-1].args += yyvsp[0].cexp);
   currentblock->close(currentblock);
    break;}
case 90:
#line 418 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);;
    break;}
case 97:
#line 432 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 98:
#line 433 "lg.y"
{yyval.cexp=C_F0(TheOperators,"+=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 99:
#line 434 "lg.y"
{yyval.cexp=C_F0(TheOperators,"-=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 100:
#line 435 "lg.y"
{yyval.cexp=C_F0(TheOperators,"*=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 101:
#line 436 "lg.y"
{yyval.cexp=C_F0(TheOperators,"/=",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 103:
#line 442 "lg.y"
{yyval.cexp=C_F0(TheOperators,"?:",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 105:
#line 446 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 106:
#line 447 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 107:
#line 448 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 108:
#line 449 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 109:
#line 450 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 110:
#line 451 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 111:
#line 452 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 112:
#line 453 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 113:
#line 454 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 114:
#line 455 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 115:
#line 456 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 116:
#line 457 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 117:
#line 458 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 118:
#line 459 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 119:
#line 460 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 120:
#line 461 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 121:
#line 462 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 122:
#line 463 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 123:
#line 464 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 124:
#line 469 "lg.y"
{yyval.cexp=yyvsp[0].cexp;
    break;}
case 125:
#line 470 "lg.y"
{yyval.cexp=C_F0(TheOperators,":");
    break;}
case 126:
#line 471 "lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 127:
#line 472 "lg.y"
{yyval.cexp=C_F0(TheOperators,":",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 128:
#line 475 "lg.y"
{yyval.args=0;
    break;}
case 129:
#line 476 "lg.y"
{yyval.args=Find(yyvsp[0].str);
    break;}
case 130:
#line 477 "lg.y"
{ yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);
    break;}
case 131:
#line 478 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 132:
#line 479 "lg.y"
{ yyval.args = (yyvsp[-2].args += Find(yyvsp[0].str)) ;
    break;}
case 133:
#line 480 "lg.y"
{ yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 134:
#line 481 "lg.y"
{ yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp)) ;
    break;}
case 135:
#line 484 "lg.y"
{yyval.args=yyvsp[0].cexp;
    break;}
case 136:
#line 485 "lg.y"
{yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;
    break;}
case 138:
#line 490 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[0].cexp);
    break;}
case 140:
#line 494 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 141:
#line 495 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);
    break;}
case 142:
#line 496 "lg.y"
{yyval.cexp=C_F0(TheOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 143:
#line 500 "lg.y"
{yyval.cexp=Find(yyvsp[0].str);;
    break;}
case 144:
#line 501 "lg.y"
{yyval.cexp= CConstant(yyvsp[0].lnum);
    break;}
case 145:
#line 502 "lg.y"
{yyval.cexp= CConstant(yyvsp[0].dnum);
    break;}
case 146:
#line 503 "lg.y"
{yyval.cexp= CConstant(complex<double>(0,yyvsp[0].dnum));
    break;}
case 147:
#line 504 "lg.y"
{yyval.cexp= CConstant<const char *>(yyvsp[0].str);
    break;}
case 148:
#line 505 "lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].args);;
    break;}
case 149:
#line 506 "lg.y"
{yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].cexp);
    break;}
case 150:
#line 507 "lg.y"
{yyval.cexp=C_F0(yyvsp[-5].cexp,yyvsp[-4].oper,yyvsp[-3].cexp,yyvsp[-1].cexp);
    break;}
case 151:
#line 508 "lg.y"
{yyval.cexp=C_F0(yyvsp[-2].cexp,"[]");
    break;}
case 152:
#line 509 "lg.y"
{ yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].str) ;;
    break;}
case 153:
#line 510 "lg.y"
{ yyval.cexp=C_F0(Find(yyvsp[-2].str),yyvsp[0].str) ;;
    break;}
case 154:
#line 511 "lg.y"
{ yyval.cexp=C_F0(Find(yyvsp[-3].str),yyvsp[-2].oper,yyvsp[-1].args) ;;
    break;}
case 155:
#line 512 "lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 156:
#line 513 "lg.y"
{yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);
    break;}
case 157:
#line 514 "lg.y"
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
#line 523 "lg.y"
{yyval.cexp=yyvsp[-1].cexp;
    break;}
case 159:
#line 524 "lg.y"
{ yyval.cexp=C_F0(TheOperators,"[]",yyvsp[-1].args);
    break;}
}

#line 705 "/sw/share/bison/bison.simple"


  yyvsp -= yylen;
  yyssp -= yylen;
#if YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG
  if (yydebug)
    {
      short *yyssp1 = yyss - 1;
      YYFPRINTF (stderr, "state stack now");
      while (yyssp1 != yyssp)
	YYFPRINTF (stderr, " %d", *++yyssp1);
      YYFPRINTF (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;
#if YYLSP_NEEDED
  *++yylsp = yyloc;
#endif

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;

#ifdef YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  char *yymsg;
	  int yyx, yycount;

	  yycount = 0;
	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  for (yyx = yyn < 0 ? -yyn : 0;
	       yyx < (int) (sizeof (yytname) / sizeof (char *)); yyx++)
	    if (yycheck[yyx + yyn] == yyx)
	      yysize += yystrlen (yytname[yyx]) + 15, yycount++;
	  yysize += yystrlen ("parse error, unexpected ") + 1;
	  yysize += yystrlen (yytname[YYTRANSLATE (yychar)]);
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "parse error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[YYTRANSLATE (yychar)]);

	      if (yycount < 5)
		{
		  yycount = 0;
		  for (yyx = yyn < 0 ? -yyn : 0;
		       yyx < (int) (sizeof (yytname) / sizeof (char *));
		       yyx++)
		    if (yycheck[yyx + yyn] == yyx)
		      {
			const char *yyq = ! yycount ? ", expecting " : " or ";
			yyp = yystpcpy (yyp, yyq);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yycount++;
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("parse error; also virtual memory exhausted");
	}
      else
#endif /* defined (YYERROR_VERBOSE) */
	yyerror ("parse error");
    }
  goto yyerrlab1;


/*--------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action |
`--------------------------------------------------*/
yyerrlab1:
  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      /* return failure if at end of input */
      if (yychar == YYEOF)
	YYABORT;
      YYDPRINTF ((stderr, "Discarding token %d (%s).\n",
		  yychar, yytname[yychar1]));
      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */

  yyerrstatus = 3;		/* Each real token shifted decrements this */

  goto yyerrhandle;


/*-------------------------------------------------------------------.
| yyerrdefault -- current state does not do anything special for the |
| error token.                                                       |
`-------------------------------------------------------------------*/
yyerrdefault:
#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */

  /* If its default is to accept any token, ok.  Otherwise pop it.  */
  yyn = yydefact[yystate];
  if (yyn)
    goto yydefault;
#endif


/*---------------------------------------------------------------.
| yyerrpop -- pop the current state because it cannot handle the |
| error token                                                    |
`---------------------------------------------------------------*/
yyerrpop:
  if (yyssp == yyss)
    YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#if YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG
  if (yydebug)
    {
      short *yyssp1 = yyss - 1;
      YYFPRINTF (stderr, "Error: state stack now");
      while (yyssp1 != yyssp)
	YYFPRINTF (stderr, " %d", *++yyssp1);
      YYFPRINTF (stderr, "\n");
    }
#endif

/*--------------.
| yyerrhandle.  |
`--------------*/
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

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;
#if YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

/*---------------------------------------------.
| yyoverflowab -- parser overflow comes here.  |
`---------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}
#line 529 "lg.y"
 


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
  zzzfff->input(cc);
  GetEnvironment();     
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


 
