/* A Bison parser, made by GNU Bison 1.875b.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

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

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0

/* If NAME_PREFIX is specified substitute the variables and functions
   names.  */
#define yyparse lgparse
#define yylex   lglex
#define yyerror lgerror
#define yylval  lglval
#define yychar  lgchar
#define yydebug lgdebug
#define yynerrs lgnerrs


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     IF = 258,
     ELSE = 259,
     SET = 260,
     GTGT = 261,
     LTLT = 262,
     OR = 263,
     AND = 264,
     NE = 265,
     EQ = 266,
     GE = 267,
     LE = 268,
     DOTSLASH = 269,
     DOTSTAR = 270,
     MOINSMOINS = 271,
     PLUSPLUS = 272,
     UNARY = 273,
     LNUM = 274,
     DNUM = 275,
     CNUM = 276,
     ID = 277,
     FESPACEID = 278,
     IDPARAM = 279,
     STRING = 280,
     ENDOFFILE = 281,
     INCLUDE = 282,
     LOAD = 283,
     BIDON = 284,
     FOR = 285,
     WHILE = 286,
     BREAK = 287,
     CONTINUE = 288,
     RETURN = 289,
     TYPE = 290,
     FUNCTION = 291,
     FESPACE = 292,
     PLUSEQ = 293,
     MOINSEQ = 294,
     MULEQ = 295,
     DIVEQ = 296,
     ARROW = 297,
     BORDER = 298,
     CURVE = 299,
     SOLVE = 300
   };
#endif
#define IF 258
#define ELSE 259
#define SET 260
#define GTGT 261
#define LTLT 262
#define OR 263
#define AND 264
#define NE 265
#define EQ 266
#define GE 267
#define LE 268
#define DOTSLASH 269
#define DOTSTAR 270
#define MOINSMOINS 271
#define PLUSPLUS 272
#define UNARY 273
#define LNUM 274
#define DNUM 275
#define CNUM 276
#define ID 277
#define FESPACEID 278
#define IDPARAM 279
#define STRING 280
#define ENDOFFILE 281
#define INCLUDE 282
#define LOAD 283
#define BIDON 284
#define FOR 285
#define WHILE 286
#define BREAK 287
#define CONTINUE 288
#define RETURN 289
#define TYPE 290
#define FUNCTION 291
#define FESPACE 292
#define PLUSEQ 293
#define MOINSEQ 294
#define MULEQ 295
#define DIVEQ 296
#define ARROW 297
#define BORDER 298
#define CURVE 299
#define SOLVE 300




/* Copy the first part of user declarations.  */
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
double CPUcompileInit =0;
//class pfes;
C_F0  fespacetype;
bool fespacecomplex;

int ShowAlloc(char *s,size_t &);
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



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 80 "lg.y"
typedef union YYSTYPE { 
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
/* Line 191 of yacc.c.  */
#line 268 "lg.tab.cpp"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 280 "lg.tab.cpp"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

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
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAXIMUM)

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
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  70
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   959

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  70
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  41
/* YYNRULES -- Number of rules. */
#define YYNRULES  158
/* YYNRULES -- Number of states. */
#define YYNSTATES  342

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   300

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    27,     2,     2,     2,    24,    12,    32,
      34,    37,    22,    20,     5,    21,    36,    23,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    69,    66,
      16,     6,    17,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    35,     2,    38,    31,    33,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    67,    10,    68,     2,     2,     2,     2,
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
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       7,     8,     9,    11,    13,    14,    15,    18,    19,    25,
      26,    28,    29,    30,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short yyprhs[] =
{
       0,     0,     3,     6,     8,    10,    13,    14,    16,    20,
      23,    27,    30,    34,    38,    42,    48,    54,    59,    65,
      70,    76,    78,    82,    84,    86,    88,    92,    97,   101,
     103,   106,   110,   114,   120,   122,   127,   134,   139,   141,
     146,   150,   154,   161,   167,   172,   179,   181,   186,   188,
     192,   194,   198,   201,   207,   212,   214,   218,   219,   224,
     228,   231,   237,   238,   249,   250,   260,   262,   264,   266,
     268,   269,   273,   275,   278,   281,   284,   286,   296,   306,
     312,   318,   326,   330,   334,   341,   344,   347,   351,   359,
     362,   364,   368,   370,   372,   374,   376,   378,   380,   384,
     388,   392,   396,   400,   402,   406,   410,   414,   418,   422,
     426,   430,   434,   438,   442,   446,   450,   454,   458,   462,
     466,   470,   474,   478,   480,   482,   486,   492,   493,   495,
     499,   501,   505,   509,   515,   517,   521,   523,   526,   528,
     532,   536,   539,   541,   543,   545,   547,   549,   554,   559,
     566,   570,   574,   578,   583,   586,   589,   594,   598
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      71,     0,    -1,    72,    46,    -1,    73,    -1,    98,    -1,
      73,    98,    -1,    -1,    76,    -1,    76,     6,   103,    -1,
      57,    76,    -1,    57,    12,    76,    -1,    79,    76,    -1,
      79,    12,    76,    -1,    35,    74,    38,    -1,    74,     5,
      76,    -1,    74,     5,    35,    74,    38,    -1,    74,     5,
      76,     6,   103,    -1,    74,     5,    57,    76,    -1,    74,
       5,    57,    12,    76,    -1,    74,     5,    79,    76,    -1,
      74,     5,    79,    12,    76,    -1,    76,    -1,    75,     5,
      76,    -1,    42,    -1,    57,    -1,    42,    -1,    42,     6,
     103,    -1,    42,    34,    78,    37,    -1,    77,     5,    77,
      -1,   104,    -1,    57,    42,    -1,    42,     6,   104,    -1,
      78,     5,   104,    -1,    78,     5,    76,     6,   104,    -1,
      55,    -1,    55,    35,    55,    38,    -1,    55,    35,    55,
       5,    55,    38,    -1,    55,    16,    55,    17,    -1,    42,
      -1,    42,    35,   104,    38,    -1,    42,     6,   104,    -1,
      35,    75,    38,    -1,    35,    75,    38,    35,   104,    38,
      -1,    35,    75,    38,     6,   104,    -1,    42,    34,   104,
      37,    -1,    35,    75,    38,    34,   104,    37,    -1,    57,
      -1,    57,    16,    55,    17,    -1,    81,    -1,    83,     5,
      81,    -1,    80,    -1,    84,     5,    80,    -1,    82,    84,
      -1,    82,    35,    55,    38,    83,    -1,    42,    34,    78,
      37,    -1,    86,    -1,    87,     5,    86,    -1,    -1,    79,
      89,    77,    66,    -1,    43,    87,    66,    -1,    85,    66,
      -1,    56,    42,     6,   101,    66,    -1,    -1,    56,    79,
      42,    34,    74,    37,    90,    67,    73,    68,    -1,    -1,
      56,    42,    34,    74,    37,    91,     6,   103,    66,    -1,
      67,    -1,    68,    -1,    50,    -1,    51,    -1,    -1,    79,
      97,    77,    -1,    66,    -1,    47,    45,    -1,    48,    45,
      -1,   101,    66,    -1,    88,    -1,    94,    34,   101,    66,
     101,    66,   101,    37,    98,    -1,    94,    34,    96,    66,
     101,    66,   101,    37,    98,    -1,    95,    34,   101,    37,
      98,    -1,     3,    34,   101,    37,    98,    -1,     3,    34,
     101,    37,    98,     4,    98,    -1,    92,    73,    93,    -1,
      63,    42,   100,    -1,    63,    42,    35,   107,    38,    66,
      -1,    52,    66,    -1,    53,    66,    -1,    54,   101,    66,
      -1,    34,    42,     6,   101,     5,   101,    37,    -1,    99,
      98,    -1,   103,    -1,   101,     5,   101,    -1,    21,    -1,
      20,    -1,    27,    -1,    29,    -1,    28,    -1,   104,    -1,
     104,     6,   103,    -1,   104,    58,   103,    -1,   104,    59,
     103,    -1,   104,    60,   103,    -1,   104,    61,   103,    -1,
     108,    -1,   104,    22,   104,    -1,   104,    26,   104,    -1,
     104,    25,   104,    -1,   104,    23,   104,    -1,   104,    24,
     104,    -1,   104,    20,   104,    -1,   104,    21,   104,    -1,
     104,     9,   104,    -1,   104,     8,   104,    -1,   104,    12,
     104,    -1,   104,    13,   104,    -1,   104,    10,   104,    -1,
     104,    11,   104,    -1,   104,    16,   104,    -1,   104,    19,
     104,    -1,   104,    17,   104,    -1,   104,    18,   104,    -1,
     104,    15,   104,    -1,   104,    14,   104,    -1,   104,    -1,
      69,    -1,   104,    69,   104,    -1,   104,    69,   104,    69,
     104,    -1,    -1,    57,    -1,    76,     6,   104,    -1,   105,
      -1,   106,     5,    57,    -1,   106,     5,   105,    -1,   106,
       5,    76,     6,   104,    -1,   103,    -1,   107,     5,   103,
      -1,   109,    -1,   102,   109,    -1,   110,    -1,   110,    31,
     108,    -1,   110,    33,   108,    -1,   110,    32,    -1,    42,
      -1,    39,    -1,    40,    -1,    41,    -1,    45,    -1,   110,
      34,   106,    37,    -1,   110,    35,   105,    38,    -1,   110,
      35,   105,     5,   105,    38,    -1,   110,    35,    38,    -1,
     110,    36,    42,    -1,    57,    36,    42,    -1,    57,    34,
     106,    37,    -1,   110,    29,    -1,   110,    28,    -1,    55,
      34,   101,    37,    -1,    34,   101,    37,    -1,    35,   107,
      38,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   198,   198,   234,   237,   238,   241,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,   252,   253,   254,
     255,   258,   259,   262,   262,   264,   265,   266,   267,   274,
     275,   276,   277,   278,   281,   282,   283,   284,   291,   292,
     293,   294,   295,   296,   299,   300,   304,   305,   310,   311,
     313,   314,   316,   317,   321,   324,   325,   328,   328,   329,
     330,   331,   333,   332,   343,   342,   349,   350,   352,   353,
     356,   356,   359,   360,   361,   362,   363,   364,   365,   369,
     370,   371,   372,   374,   376,   379,   383,   387,   395,   401,
     407,   408,   413,   414,   415,   416,   417,   421,   422,   423,
     424,   425,   426,   430,   431,   432,   433,   434,   435,   436,
     437,   438,   439,   440,   441,   442,   443,   444,   445,   446,
     447,   448,   449,   454,   455,   456,   457,   460,   461,   462,
     463,   464,   465,   466,   469,   470,   474,   475,   478,   479,
     480,   481,   485,   486,   487,   488,   489,   490,   491,   492,
     493,   494,   495,   496,   497,   498,   499,   508,   509
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "IF", "ELSE", "','", "'='", "SET", "GTGT", 
  "LTLT", "'|'", "OR", "'&'", "AND", "NE", "EQ", "'<'", "'>'", "GE", "LE", 
  "'+'", "'-'", "'*'", "'/'", "'%'", "DOTSLASH", "DOTSTAR", "'!'", 
  "MOINSMOINS", "PLUSPLUS", "UNARY", "'^'", "'''", "'_'", "'('", "'['", 
  "'.'", "')'", "']'", "LNUM", "DNUM", "CNUM", "ID", "FESPACEID", 
  "IDPARAM", "STRING", "ENDOFFILE", "INCLUDE", "LOAD", "BIDON", "FOR", 
  "WHILE", "BREAK", "CONTINUE", "RETURN", "TYPE", "FUNCTION", "FESPACE", 
  "PLUSEQ", "MOINSEQ", "MULEQ", "DIVEQ", "ARROW", "BORDER", "CURVE", 
  "SOLVE", "';'", "'{'", "'}'", "':'", "$accept", "start", "input", 
  "instructions", "list_of_id_args", "list_of_id1", "id", "list_of_dcls", 
  "parameters_list", "type_of_dcl", "ID_space", "ID_array_space", 
  "fespace", "spaceIDa", "spaceIDb", "spaceIDs", "fespace_def", 
  "fespace_def_list", "declaration", "@1", "@2", "@3", "begin", "end", 
  "for_loop", "while_loop", "declaration_for", "@4", "instruction", 
  "bornes", "border_expr", "Expr", "unop", "no_comma_expr", "no_set_expr", 
  "sub_script_expr", "parameters", "array", "unary_expr", "pow_expr", 
  "primary", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
{
       0,   256,   257,   258,   259,    44,    61,   260,   261,   262,
     124,   263,    38,   264,   265,   266,    60,    62,   267,   268,
      43,    45,    42,    47,    37,   269,   270,    33,   271,   272,
     273,    94,    39,    95,    40,    91,    46,    41,    93,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,    59,   123,   125,    58
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    70,    71,    72,    73,    73,    74,    74,    74,    74,
      74,    74,    74,    74,    74,    74,    74,    74,    74,    74,
      74,    75,    75,    76,    76,    77,    77,    77,    77,    78,
      78,    78,    78,    78,    79,    79,    79,    79,    80,    80,
      80,    80,    80,    80,    81,    81,    82,    82,    83,    83,
      84,    84,    85,    85,    86,    87,    87,    89,    88,    88,
      88,    88,    90,    88,    91,    88,    92,    93,    94,    95,
      97,    96,    98,    98,    98,    98,    98,    98,    98,    98,
      98,    98,    98,    98,    98,    98,    98,    98,    99,   100,
     101,   101,   102,   102,   102,   102,   102,   103,   103,   103,
     103,   103,   103,   104,   104,   104,   104,   104,   104,   104,
     104,   104,   104,   104,   104,   104,   104,   104,   104,   104,
     104,   104,   104,   105,   105,   105,   105,   106,   106,   106,
     106,   106,   106,   106,   107,   107,   108,   108,   109,   109,
     109,   109,   110,   110,   110,   110,   110,   110,   110,   110,
     110,   110,   110,   110,   110,   110,   110,   110,   110
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     2,     1,     1,     2,     0,     1,     3,     2,
       3,     2,     3,     3,     3,     5,     5,     4,     5,     4,
       5,     1,     3,     1,     1,     1,     3,     4,     3,     1,
       2,     3,     3,     5,     1,     4,     6,     4,     1,     4,
       3,     3,     6,     5,     4,     6,     1,     4,     1,     3,
       1,     3,     2,     5,     4,     1,     3,     0,     4,     3,
       2,     5,     0,    10,     0,     9,     1,     1,     1,     1,
       0,     3,     1,     2,     2,     2,     1,     9,     9,     5,
       5,     7,     3,     3,     6,     2,     2,     3,     7,     2,
       1,     3,     1,     1,     1,     1,     1,     1,     3,     3,
       3,     3,     3,     1,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     1,     1,     3,     5,     0,     1,     3,
       1,     3,     3,     5,     1,     3,     1,     2,     1,     3,
       3,     2,     1,     1,     1,     1,     1,     4,     4,     6,
       3,     3,     3,     4,     2,     2,     4,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       0,     0,    93,    92,    94,    96,    95,     0,     0,   143,
     144,   145,   142,     0,   146,     0,     0,    68,    69,     0,
       0,     0,    34,     0,    46,     0,    72,    66,     0,     0,
       3,    57,     0,     0,    76,     0,     0,     0,     4,     0,
       0,    90,    97,   103,   136,   138,     0,     0,     0,     0,
     134,     0,     0,    55,     0,    73,    74,    85,    86,     0,
       0,     0,     0,     0,    34,     0,     0,   127,     0,     0,
       1,     2,     5,     0,     0,    38,    50,    52,    60,     0,
       0,     0,     0,    75,   137,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   155,
     154,     0,   141,     0,   127,     0,     0,     0,   157,     0,
     158,     0,     0,    59,    87,     0,     0,     0,     0,     6,
       0,     0,   142,   128,   124,     0,   123,   130,     0,   152,
       0,     0,     0,    83,    25,     0,    23,     0,    24,     0,
      21,     0,     0,     0,    67,    82,    70,     0,     0,     0,
      91,    98,   112,   111,   115,   116,   113,   114,   122,   121,
     117,   119,   120,   118,   109,   110,   104,   107,   108,   106,
     105,    99,   100,   101,   102,   139,   140,     0,   150,     0,
     151,     0,   135,   142,     0,     0,    29,    56,    37,   156,
       0,    35,     0,     6,    24,     0,     7,     0,     6,    47,
       0,     0,     0,   153,     0,     0,    89,     0,     0,     0,
      58,     0,     0,    41,    40,     0,     0,    51,     0,     0,
       0,     0,   147,     0,   148,    80,     0,    30,     0,    54,
       0,    61,     0,     0,     9,     0,    64,     0,     0,    11,
       0,   129,   125,   131,     0,   132,     0,     0,    26,     0,
      28,     0,     0,    48,    53,    22,     0,     0,    39,    71,
       0,     0,    79,     0,     0,    31,    24,     0,    32,    36,
      13,    10,     6,    24,    14,     0,     0,     8,    12,    62,
       0,     0,     0,    84,    27,     0,     0,     0,    43,     0,
       0,     0,   149,    81,     0,     0,     0,    17,     0,     0,
      19,     0,     0,   126,   133,     0,     0,     0,    49,    42,
       0,     0,    33,    15,    18,    16,    20,     0,     0,    91,
       0,    44,     0,     0,    65,     0,    88,     0,    78,    77,
      63,    45
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,    28,    29,    30,   205,   149,   206,   145,   195,    31,
      76,   263,    32,   264,    77,    33,    53,    54,    34,    73,
     312,   286,    35,   155,    36,    37,   157,   228,    38,   142,
     143,    39,    40,    41,    42,   137,   138,    51,    43,    44,
      45
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -198
static const short yypact[] =
{
     469,    11,  -198,  -198,  -198,  -198,  -198,    96,    96,  -198,
    -198,  -198,  -198,    90,  -198,    -4,    94,  -198,  -198,   134,
     138,    96,   216,   108,   169,   144,  -198,  -198,   157,   178,
     469,  -198,    68,   149,  -198,   469,   184,   187,  -198,     0,
     156,  -198,   550,  -198,  -198,   238,    96,   213,    28,    54,
    -198,    41,   218,  -198,     3,  -198,  -198,  -198,  -198,     4,
     198,    96,   203,    38,   196,   217,   205,   560,   222,   120,
    -198,  -198,  -198,   223,   206,    45,  -198,   270,  -198,   359,
     654,    96,    96,  -198,  -198,    96,    96,    96,    96,    96,
      96,    96,    96,    96,    96,    96,    96,    96,    96,    96,
      96,    96,    96,    96,    96,    96,    96,    96,    96,  -198,
    -198,    96,  -198,    96,   560,   199,   234,    82,  -198,    96,
    -198,   685,    90,  -198,  -198,   260,    97,    44,    96,   152,
     244,   263,   275,     6,  -198,   277,     7,  -198,   106,  -198,
     242,    96,   469,  -198,   150,     5,  -198,   247,  -198,    48,
    -198,    96,    96,   175,  -198,  -198,  -198,   220,    29,   107,
    -198,  -198,   905,   905,   920,   920,   933,   933,   308,   308,
     317,   317,   317,   317,   324,   324,  -198,  -198,  -198,  -198,
    -198,  -198,  -198,  -198,  -198,  -198,  -198,   109,  -198,    50,
    -198,   469,  -198,   285,   172,   110,   888,  -198,  -198,  -198,
     237,  -198,    30,   152,     1,   115,   292,    26,   152,  -198,
      96,    96,   592,  -198,   293,    51,  -198,    96,   685,   223,
    -198,   201,    91,    46,   888,   766,    91,  -198,   223,    96,
      96,   469,  -198,   623,  -198,   296,    96,  -198,   716,  -198,
     264,  -198,    52,    91,  -198,   200,  -198,    96,    91,  -198,
     121,   888,   529,     6,   299,  -198,    96,   243,  -198,   122,
    -198,    91,   278,  -198,   310,  -198,    96,    96,  -198,   312,
      31,    32,  -198,   281,   469,   888,    28,   315,   888,  -198,
    -198,  -198,   152,    35,   316,    36,   329,  -198,  -198,  -198,
      96,    96,   346,  -198,  -198,    56,    96,   201,   888,   797,
      96,    96,  -198,  -198,    96,    62,    91,  -198,    96,    91,
    -198,    96,   286,   888,   888,    96,   318,   828,  -198,  -198,
     123,   124,   888,  -198,  -198,  -198,  -198,   288,   469,   319,
      96,  -198,   469,   469,  -198,   414,  -198,   858,  -198,  -198,
    -198,  -198
};

/* YYPGOTO[NTERM-NUM].  */
static const short yypgoto[] =
{
    -198,  -198,  -198,   -33,  -197,    98,    75,  -120,   137,   -20,
     204,    64,  -198,  -198,  -198,  -198,   241,  -198,  -198,  -198,
    -198,  -198,  -198,  -198,  -198,  -198,  -198,  -198,   -29,  -198,
    -198,    -7,  -198,    -1,    78,  -111,   250,   224,   -48,   326,
    -198
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -25
static const short yytable[] =
{
      49,    72,    79,    65,   189,    82,   242,    50,   122,    82,
     219,   250,   -24,   243,    59,    86,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,    82,    82,    82,    82,   248,   117,
      67,    55,    68,   146,   128,    46,   119,   306,   309,   200,
      72,   151,   266,   222,   126,   233,   119,   245,   148,    82,
     156,   222,    67,   185,    68,   186,    83,   245,   146,   123,
     124,   220,   129,   158,   159,   160,   211,   146,   146,   120,
     152,   267,   201,   148,   161,   305,   223,    82,   234,   257,
     280,   118,   148,   148,   316,   230,   241,   300,   301,   260,
     323,   255,    82,    74,   181,   182,   183,   184,   269,   207,
      75,   212,    82,   216,   212,   238,     2,     3,   192,   191,
     245,   202,   273,     4,     5,     6,   245,   238,    82,    82,
       7,     8,    52,   146,   199,     9,    10,    11,    12,    56,
      50,    14,   135,   213,   231,   136,   232,   239,   148,   150,
      63,    47,   246,    48,   140,   141,   217,    70,   289,   294,
     332,   333,   235,    64,   162,   163,   164,   165,   166,   167,
     168,   169,   170,   171,   172,   173,   174,   175,   176,   177,
     178,   179,   180,   207,   218,    66,    69,   203,   207,   135,
       7,     8,   136,   136,   146,     9,    10,    11,    12,   196,
      57,    14,   272,    67,    58,    68,    67,    64,    68,   204,
     226,    47,    60,    48,   237,    78,   258,    75,    80,     2,
       3,    81,   270,   271,    71,   285,     4,     5,     6,   224,
     225,    62,    60,     7,     8,   282,   261,   188,     9,    10,
      11,    12,   146,   262,    14,   303,   287,    61,   146,   292,
      61,    62,   121,   125,    47,    64,    48,   283,   127,   130,
     131,   147,   207,   148,   139,   144,   109,   110,   134,   111,
     112,   113,   114,   115,   116,   153,   190,   198,   208,   244,
     209,   -23,   249,   210,   214,   221,   229,   254,   251,   252,
     136,   236,   240,   320,   321,   335,   196,   265,   247,   256,
     274,   150,   279,   338,   339,   291,    72,   325,   329,   293,
     327,   136,   296,   277,   275,   297,   278,   219,   281,   302,
     284,   304,   308,   288,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   311,   150,    98,    99,   100,
     101,   102,   103,   104,   298,   299,   100,   101,   102,   103,
     104,   315,   330,   328,   334,   259,   336,   227,   307,   295,
     310,   318,     1,   197,   187,   215,    84,     0,   313,   314,
       0,     0,     0,     0,   317,     0,     0,     0,     0,     2,
       3,   324,   322,     0,   326,     0,     4,     5,     6,     0,
       0,     0,     0,     7,     8,     0,     0,     0,     9,    10,
      11,    12,    13,     0,    14,     0,    15,    16,   337,    17,
      18,    19,    20,    21,    22,    23,    24,     1,     0,     0,
       0,     0,    25,     0,     0,    26,    27,   154,     0,     0,
       0,     0,     0,     0,     2,     3,     0,     0,     0,     0,
       0,     4,     5,     6,     0,     0,     0,     0,     7,     8,
       0,     0,     0,     9,    10,    11,    12,    13,     0,    14,
       0,    15,    16,     0,    17,    18,    19,    20,    21,    22,
      23,    24,     1,     0,     0,     0,     0,    25,     0,     0,
      26,    27,   340,     0,     0,     0,     0,     0,     0,     2,
       3,     0,     0,     0,     0,     0,     4,     5,     6,     0,
       0,     0,     0,     7,     8,     0,     0,     0,     9,    10,
      11,    12,    13,     0,    14,     0,    15,    16,     0,    17,
      18,    19,    20,    21,    22,    23,    24,     0,     0,     0,
       0,     0,    25,     0,     0,    26,    27,    86,    87,    88,
      89,    90,    91,    92,    93,    94,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,    85,     0,    86,    87,
      88,    89,    90,    91,    92,    93,    94,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,     0,     0,     0,
       2,     3,     0,     0,     0,     0,     0,     4,     5,     6,
       0,     0,     0,     0,     7,     8,     0,     0,   290,     9,
      10,    11,   132,     0,     0,    14,     0,     0,   105,   106,
     107,   108,     2,     3,     0,    47,     0,   133,     0,     4,
       5,     6,     0,     0,     0,     0,     7,     8,     0,   134,
       0,     9,    10,    11,   132,     0,     0,    14,     0,     0,
       0,     0,     0,     2,     3,     0,     0,    47,     0,   253,
       4,     5,     6,     0,     0,     0,     0,     7,     8,     0,
       0,   134,     9,    10,    11,    12,     0,     0,    14,     0,
       0,     0,     0,     0,     2,     3,     0,     0,    47,     0,
      48,     4,     5,     6,     0,     0,     0,     0,     7,     8,
       0,     0,   134,     9,    10,    11,    12,     0,     0,    14,
       0,     0,     0,     0,     0,     2,     3,     0,     0,    22,
       0,    48,     4,     5,     6,     0,     0,     0,     0,     7,
       8,     0,     0,     0,     9,    10,    11,   193,     0,     0,
      14,     0,     0,     0,     0,     0,     2,     3,     0,     0,
      47,     0,   194,     4,     5,     6,     0,     0,     0,     0,
       7,     8,     0,     0,     0,     9,    10,    11,   132,     0,
       0,    14,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    47,     0,   276,    86,    87,    88,    89,    90,    91,
      92,    93,    94,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   268,    86,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,    97,    98,    99,   100,
     101,   102,   103,   104,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   319,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   331,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   341,    86,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,    88,    89,    90,    91,    92,
      93,    94,    95,    96,    97,    98,    99,   100,   101,   102,
     103,   104,    90,    91,    92,    93,    94,    95,    96,    97,
      98,    99,   100,   101,   102,   103,   104,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104
};

static const short yycheck[] =
{
       7,    30,    35,    23,   115,     5,   203,     8,     5,     5,
       5,   208,     6,    12,    21,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,     5,     5,     5,     5,    12,    46,
      34,    45,    36,    42,     6,    34,     5,    12,    12,     5,
      79,     6,     6,     5,    61,     5,     5,     5,    57,     5,
      80,     5,    34,   111,    36,   113,    66,     5,    42,    66,
      66,    66,    34,    80,    81,    82,    69,    42,    42,    38,
      35,    35,    38,    57,    85,   282,    38,     5,    38,    38,
      38,    37,    57,    57,    38,    66,    66,    66,    66,   219,
      38,   212,     5,    35,   105,   106,   107,   108,   228,   129,
      42,     5,     5,   142,     5,     5,    20,    21,   119,    37,
       5,   128,   233,    27,    28,    29,     5,     5,     5,     5,
      34,    35,    42,    42,    37,    39,    40,    41,    42,    45,
     141,    45,    67,    37,    37,    67,    37,    37,    57,    74,
      42,    55,    37,    57,    34,    35,     6,     0,    37,    37,
      37,    37,   191,    55,    86,    87,    88,    89,    90,    91,
      92,    93,    94,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   203,    34,    16,    42,    35,   208,   114,
      34,    35,   114,   115,    42,    39,    40,    41,    42,   121,
      66,    45,   231,    34,    66,    36,    34,    55,    36,    57,
      35,    55,    16,    57,    42,    66,   217,    42,    34,    20,
      21,    34,   229,   230,    46,   245,    27,    28,    29,   151,
     152,    35,    16,    34,    35,    35,    35,    38,    39,    40,
      41,    42,    42,    42,    45,   274,   247,    34,    42,   256,
      34,    35,    34,    55,    55,    55,    57,    57,    55,    42,
      55,    55,   282,    57,    42,    42,    28,    29,    69,    31,
      32,    33,    34,    35,    36,     5,    42,    17,    34,   204,
      17,     6,   207,     6,    42,    38,    66,   212,   210,   211,
     212,     6,    55,   300,   301,   328,   218,   222,     6,     6,
       4,   226,    38,   332,   333,     6,   335,   308,   315,    66,
     311,   233,    34,   238,   236,     5,   238,     5,   243,    38,
     245,     6,     6,   248,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,     6,   261,    20,    21,    22,
      23,    24,    25,    26,   266,   267,    22,    23,    24,    25,
      26,     5,    34,    67,    66,   218,    37,   153,   283,   261,
     285,   297,     3,   122,   114,   141,    40,    -1,   290,   291,
      -1,    -1,    -1,    -1,   296,    -1,    -1,    -1,    -1,    20,
      21,   306,   304,    -1,   309,    -1,    27,    28,    29,    -1,
      -1,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,
      41,    42,    43,    -1,    45,    -1,    47,    48,   330,    50,
      51,    52,    53,    54,    55,    56,    57,     3,    -1,    -1,
      -1,    -1,    63,    -1,    -1,    66,    67,    68,    -1,    -1,
      -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    -1,
      -1,    27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,
      -1,    -1,    -1,    39,    40,    41,    42,    43,    -1,    45,
      -1,    47,    48,    -1,    50,    51,    52,    53,    54,    55,
      56,    57,     3,    -1,    -1,    -1,    -1,    63,    -1,    -1,
      66,    67,    68,    -1,    -1,    -1,    -1,    -1,    -1,    20,
      21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,    -1,
      -1,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,
      41,    42,    43,    -1,    45,    -1,    47,    48,    -1,    50,
      51,    52,    53,    54,    55,    56,    57,    -1,    -1,    -1,
      -1,    -1,    63,    -1,    -1,    66,    67,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,     6,    -1,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    -1,    -1,    -1,
      20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,
      -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    69,    39,
      40,    41,    42,    -1,    -1,    45,    -1,    -1,    58,    59,
      60,    61,    20,    21,    -1,    55,    -1,    57,    -1,    27,
      28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,    69,
      -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,
      -1,    -1,    -1,    20,    21,    -1,    -1,    55,    -1,    57,
      27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,
      -1,    69,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    55,    -1,
      57,    27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,
      -1,    -1,    69,    39,    40,    41,    42,    -1,    -1,    45,
      -1,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,    55,
      -1,    57,    27,    28,    29,    -1,    -1,    -1,    -1,    34,
      35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,
      55,    -1,    57,    27,    28,    29,    -1,    -1,    -1,    -1,
      34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,
      -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    55,    -1,    57,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    38,     8,     9,    10,    11,    12,
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

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     3,    20,    21,    27,    28,    29,    34,    35,    39,
      40,    41,    42,    43,    45,    47,    48,    50,    51,    52,
      53,    54,    55,    56,    57,    63,    66,    67,    71,    72,
      73,    79,    82,    85,    88,    92,    94,    95,    98,   101,
     102,   103,   104,   108,   109,   110,    34,    55,    57,   101,
     103,   107,    42,    86,    87,    45,    45,    66,    66,   101,
      16,    34,    35,    42,    55,    79,    16,    34,    36,    42,
       0,    46,    98,    89,    35,    42,    80,    84,    66,    73,
      34,    34,     5,    66,   109,     6,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    58,    59,    60,    61,    28,
      29,    31,    32,    33,    34,    35,    36,   101,    37,     5,
      38,    34,     5,    66,    66,    55,   101,    55,     6,    34,
      42,    55,    42,    57,    69,    76,   104,   105,   106,    42,
      34,    35,    99,   100,    42,    77,    42,    55,    57,    75,
      76,     6,    35,     5,    68,    93,    79,    96,   101,   101,
     101,   103,   104,   104,   104,   104,   104,   104,   104,   104,
     104,   104,   104,   104,   104,   104,   104,   104,   104,   104,
     104,   103,   103,   103,   103,   108,   108,   106,    38,   105,
      42,    37,   103,    42,    57,    78,   104,    86,    17,    37,
       5,    38,   101,    35,    57,    74,    76,    79,    34,    17,
       6,    69,     5,    37,    42,   107,    98,     6,    34,     5,
      66,    38,     5,    38,   104,   104,    35,    80,    97,    66,
      66,    37,    37,     5,    38,    98,     6,    42,     5,    37,
      55,    66,    74,    12,    76,     5,    37,     6,    12,    76,
      74,   104,   104,    57,    76,   105,     6,    38,   103,    78,
      77,    35,    42,    81,    83,    76,     6,    35,    38,    77,
     101,   101,    98,   105,     4,   104,    57,    76,   104,    38,
      38,    76,    35,    57,    76,    79,    91,   103,    76,    37,
      69,     6,   101,    66,    37,    75,    34,     5,   104,   104,
      66,    66,    38,    98,     6,    74,    12,    76,     6,    12,
      76,     6,    90,   104,   104,     5,    38,   104,    81,    38,
     101,   101,   104,    38,    76,   103,    76,   103,    67,   101,
      34,    37,    37,    37,    66,    73,    37,   104,    98,    98,
      68,    37
};

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
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
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
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)         \
  Current.first_line   = Rhs[1].first_line;      \
  Current.first_column = Rhs[1].first_column;    \
  Current.last_line    = Rhs[N].last_line;       \
  Current.last_column  = Rhs[N].last_column;
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

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

# define YYDSYMPRINT(Args)			\
do {						\
  if (yydebug)					\
    yysymprint Args;				\
} while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Token, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (cinluded).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short *bottom, short *top)
#else
static void
yy_stack_print (bottom, top)
    short *bottom;
    short *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
# define YYDSYMPRINTF(Title, Token, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
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



#if YYERROR_VERBOSE

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

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
    }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
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

  if (yyss + yystacksize - 1 <= yyssp)
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
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
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
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YYDSYMPRINTF ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


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

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 198 "lg.y"
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
  return 0;;}
    break;

  case 4:
#line 237 "lg.y"
    {yyval.cinst=yyvsp[0].cexp;;;;}
    break;

  case 5:
#line 238 "lg.y"
    { yyval.cinst= (yyvsp[-1].cinst+=yyvsp[0].cexp) ;}
    break;

  case 6:
#line 241 "lg.y"
    { yyval.clist_id=new ListOfId();;}
    break;

  case 7:
#line 242 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str));}
    break;

  case 8:
#line 243 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;}
    break;

  case 9:
#line 244 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>()));}
    break;

  case 10:
#line 245 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true));}
    break;

  case 11:
#line 246 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;}
    break;

  case 12:
#line 247 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;}
    break;

  case 13:
#line 248 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;}
    break;

  case 14:
#line 249 "lg.y"
    { yyval.clist_id = yyvsp[-2].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str)) ;}
    break;

  case 15:
#line 250 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;}
    break;

  case 16:
#line 251 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;}
    break;

  case 17:
#line 252 "lg.y"
    { yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>())) ;}
    break;

  case 18:
#line 253 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true)) ;}
    break;

  case 19:
#line 254 "lg.y"
    { yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;}
    break;

  case 20:
#line 255 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;}
    break;

  case 21:
#line 258 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;}
    break;

  case 22:
#line 259 "lg.y"
    { yyval.clist_id=yyvsp[-2].clist_id  ; yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;}
    break;

  case 25:
#line 264 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[0].str,dcltype);}
    break;

  case 26:
#line 265 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-2].str,dcltype,yyvsp[0].cexp);}
    break;

  case 27:
#line 266 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-3].str,dcltype,yyvsp[-1].args);}
    break;

  case 28:
#line 267 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 29:
#line 274 "lg.y"
    {yyval.args=yyvsp[0].cexp;}
    break;

  case 30:
#line 275 "lg.y"
    {yyval.args=Find(yyvsp[-1].str);}
    break;

  case 31:
#line 276 "lg.y"
    { yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);}
    break;

  case 32:
#line 277 "lg.y"
    { yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;}
    break;

  case 33:
#line 278 "lg.y"
    { yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp));}
    break;

  case 35:
#line 282 "lg.y"
    {yyval.type=TypeArray(yyvsp[-3].type,yyvsp[-1].type);}
    break;

  case 36:
#line 283 "lg.y"
    {yyval.type=TypeArray(yyvsp[-5].type,yyvsp[-3].type,yyvsp[-1].type);}
    break;

  case 37:
#line 284 "lg.y"
    {yyval.type=TypeTemplate(yyvsp[-3].type,yyvsp[-1].type);}
    break;

  case 38:
#line 291 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[0].str,currentblock,fespacetype,fespacecomplex); ;}
    break;

  case 39:
#line 292 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;}
    break;

  case 40:
#line 293 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[-2].str,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;}
    break;

  case 41:
#line 294 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[-1].clist_id,currentblock,fespacetype,fespacecomplex) ;}
    break;

  case 42:
#line 295 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;}
    break;

  case 43:
#line 296 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[-3].clist_id,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;}
    break;

  case 44:
#line 299 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;}
    break;

  case 45:
#line 300 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;}
    break;

  case 46:
#line 304 "lg.y"
    {fespacecomplex=false;  fespacetype = Find(yyvsp[0].str);;}
    break;

  case 47:
#line 305 "lg.y"
    {
             if (yyvsp[-1].type != typevarreal && yyvsp[-1].type != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=(yyvsp[-1].type==typevarcomplex);
             fespacetype = Find(yyvsp[-3].str);;}
    break;

  case 48:
#line 310 "lg.y"
    {  yyval.cexp = yyvsp[0].cexp  ;}
    break;

  case 49:
#line 311 "lg.y"
    { yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;}
    break;

  case 50:
#line 313 "lg.y"
    {  yyval.cexp = yyvsp[0].cexp  ;}
    break;

  case 51:
#line 314 "lg.y"
    { yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;}
    break;

  case 52:
#line 316 "lg.y"
    { yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;}
    break;

  case 53:
#line 317 "lg.y"
    { yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;}
    break;

  case 54:
#line 322 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariableFES,size_t>(yyvsp[-3].str,atype<pfes*>(),yyvsp[-1].args,dimFESpaceImage(yyvsp[-1].args));}
    break;

  case 56:
#line 325 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 57:
#line 328 "lg.y"
    {dcltype=yyvsp[0].type;}
    break;

  case 58:
#line 328 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 59:
#line 329 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 60:
#line 330 "lg.y"
    { yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 61:
#line 331 "lg.y"
    {yyval.cexp=currentblock->NewID(yyvsp[-4].type,yyvsp[-3].str,yyvsp[-1].cexp);;}
    break;

  case 62:
#line 333 "lg.y"
    {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = yyvsp[-4].type->right();
                      yyvsp[-1].routine=new Routine(yyvsp[-5].type,yyvsp[-4].type->right(),yyvsp[-3].str,yyvsp[-1].clist_id,currentblock);;}
    break;

  case 63:
#line 338 "lg.y"
    { currentblock=yyvsp[-5].routine->Set(yyvsp[-1].cinst);
                       currentblock->Add(yyvsp[-7].str,"(",yyvsp[-5].routine);
                       kkembtype--;
                       yyval.cexp=0 ;}
    break;

  case 64:
#line 343 "lg.y"
    {currentblock = new Block(currentblock); yyvsp[-4].type->SetArgs(yyvsp[-1].clist_id);;}
    break;

  case 65:
#line 345 "lg.y"
    {  yyval.cinst=currentblock->close(currentblock);
                         yyval.cexp=currentblock->NewID(yyvsp[-8].type,yyvsp[-7].str,yyvsp[-1].cexp,*yyvsp[-5].clist_id);}
    break;

  case 66:
#line 349 "lg.y"
    {  currentblock = new Block(currentblock);}
    break;

  case 67:
#line 350 "lg.y"
    {  yyval.cexp=currentblock->close(currentblock);}
    break;

  case 68:
#line 352 "lg.y"
    {inloopcount++;;}
    break;

  case 69:
#line 353 "lg.y"
    {inloopcount++;}
    break;

  case 70:
#line 356 "lg.y"
    {dcltype=yyvsp[0].type;currentblock = new Block(currentblock);}
    break;

  case 71:
#line 357 "lg.y"
    {yyval.cexp=yyvsp[0].cexp;}
    break;

  case 72:
#line 359 "lg.y"
    {yyval.cexp=0;;}
    break;

  case 73:
#line 360 "lg.y"
    {zzzfff->input(yyvsp[0].str);yyval.cexp= 0; ;}
    break;

  case 74:
#line 361 "lg.y"
    {load(yyvsp[0].str);yyval.cexp= 0; ;}
    break;

  case 75:
#line 362 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 76:
#line 363 "lg.y"
    {yyval.cexp=yyvsp[0].cexp;}
    break;

  case 77:
#line 364 "lg.y"
    {inloopcount--; yyval.cexp=For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 78:
#line 366 "lg.y"
    {inloopcount--; 
                yyval.cexp=C_F0(For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp),currentblock->close(currentblock));}
    break;

  case 79:
#line 369 "lg.y"
    {inloopcount--;yyval.cexp=While(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 80:
#line 370 "lg.y"
    {yyval.cexp=FIf(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 81:
#line 371 "lg.y"
    {yyval.cexp=FIf(yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 82:
#line 372 "lg.y"
    { 
                      yyval.cexp=C_F0(new E_block(yyvsp[-1].cinst,yyvsp[0].cexp),atype<void>()) ;}
    break;

  case 83:
#line 374 "lg.y"
    {
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-1].str,C_F0(TheOperators,"[border]",yyvsp[0].args));}
    break;

  case 84:
#line 376 "lg.y"
    {
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-4].str,C_F0(TheOperators,"[border]",yyvsp[-2].args));}
    break;

  case 85:
#line 379 "lg.y"
    {
                    if(inloopcount) 
                      yyval.cexp= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;}
    break;

  case 86:
#line 383 "lg.y"
    { 
                    if(inloopcount)
                        yyval.cexp= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");}
    break;

  case 87:
#line 387 "lg.y"
    { 
                    if (kkembtype>=0)
                      yyval.cexp= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo(yyvsp[-1].cexp)) ,atype<void>());
                     else lgerror(" return not in routine ") ;}
    break;

  case 88:
#line 395 "lg.y"
    { 
   currentblock = new Block(currentblock);
   yyval.args = currentblock->NewVar<LocalVariable>(yyvsp[-5].str,atype<double*>());
   yyval.args+= yyvsp[-3].cexp;
   yyval.args+= yyvsp[-1].cexp ;}
    break;

  case 89:
#line 401 "lg.y"
    {   
   yyval.args = (yyvsp[-1].args += yyvsp[0].cexp);
   currentblock->close(currentblock);}
    break;

  case 91:
#line 408 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);;}
    break;

  case 98:
#line 422 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 99:
#line 423 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"+=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 100:
#line 424 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"-=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 101:
#line 425 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"*=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 102:
#line 426 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"/=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 104:
#line 431 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 105:
#line 432 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 106:
#line 433 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 107:
#line 434 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 108:
#line 435 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 109:
#line 436 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 110:
#line 437 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 111:
#line 438 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 112:
#line 439 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 113:
#line 440 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 114:
#line 441 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 115:
#line 442 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 116:
#line 443 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 117:
#line 444 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 118:
#line 445 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 119:
#line 446 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 120:
#line 447 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 121:
#line 448 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 122:
#line 449 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 123:
#line 454 "lg.y"
    {yyval.cexp=yyvsp[0].cexp;}
    break;

  case 124:
#line 455 "lg.y"
    {yyval.cexp=C_F0(TheOperators,":");}
    break;

  case 125:
#line 456 "lg.y"
    {yyval.cexp=C_F0(TheOperators,":",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 126:
#line 457 "lg.y"
    {yyval.cexp=C_F0(TheOperators,":",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 127:
#line 460 "lg.y"
    {yyval.args=0;}
    break;

  case 128:
#line 461 "lg.y"
    {yyval.args=Find(yyvsp[0].str);}
    break;

  case 129:
#line 462 "lg.y"
    { yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);}
    break;

  case 130:
#line 463 "lg.y"
    {yyval.args=yyvsp[0].cexp;}
    break;

  case 131:
#line 464 "lg.y"
    { yyval.args = (yyvsp[-2].args += Find(yyvsp[0].str)) ;}
    break;

  case 132:
#line 465 "lg.y"
    { yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;}
    break;

  case 133:
#line 466 "lg.y"
    { yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp)) ;}
    break;

  case 134:
#line 469 "lg.y"
    {yyval.args=yyvsp[0].cexp;}
    break;

  case 135:
#line 470 "lg.y"
    {yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;}
    break;

  case 137:
#line 475 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[0].cexp);}
    break;

  case 139:
#line 479 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 140:
#line 480 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 141:
#line 481 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[0].oper,yyvsp[-1].cexp);}
    break;

  case 142:
#line 485 "lg.y"
    {yyval.cexp=Find(yyvsp[0].str);;}
    break;

  case 143:
#line 486 "lg.y"
    {yyval.cexp= CConstant(yyvsp[0].lnum);}
    break;

  case 144:
#line 487 "lg.y"
    {yyval.cexp= CConstant(yyvsp[0].dnum);}
    break;

  case 145:
#line 488 "lg.y"
    {yyval.cexp= CConstant(complex<double>(0,yyvsp[0].dnum));}
    break;

  case 146:
#line 489 "lg.y"
    {yyval.cexp= CConstant<const char *>(yyvsp[0].str);}
    break;

  case 147:
#line 490 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].args);;}
    break;

  case 148:
#line 491 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].cexp);}
    break;

  case 149:
#line 492 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-5].cexp,yyvsp[-4].oper,yyvsp[-3].cexp,yyvsp[-1].cexp);}
    break;

  case 150:
#line 493 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-2].cexp,"[]");}
    break;

  case 151:
#line 494 "lg.y"
    { yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].str) ;;}
    break;

  case 152:
#line 495 "lg.y"
    { yyval.cexp=C_F0(Find(yyvsp[-2].str),yyvsp[0].str) ;;}
    break;

  case 153:
#line 496 "lg.y"
    { yyval.cexp=C_F0(Find(yyvsp[-3].str),yyvsp[-2].oper,yyvsp[-1].args) ;;}
    break;

  case 154:
#line 497 "lg.y"
    {yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);}
    break;

  case 155:
#line 498 "lg.y"
    {yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);}
    break;

  case 156:
#line 499 "lg.y"
    {
             if (yyvsp[-3].type->right()->CastingFrom(yyvsp[-1].cexp.left()) ) 
                yyval.cexp=yyvsp[-3].type->right()->CastTo(yyvsp[-1].cexp)  ;
             else { yyval.cexp=yyvsp[-3].type->right()->Find("<--",basicAC_F0_wa(yyvsp[-1].cexp));
             if (!yyval.cexp.left()) { cerr << " no wait to change " << yyvsp[-1].cexp.left()->right()->name() << " in " << 
                                        yyvsp[-3].type->right()->name() << endl;
                                CompileError(" Error in type(exp) "); }
             }
            ;}
    break;

  case 157:
#line 508 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 158:
#line 509 "lg.y"
    { yyval.cexp=C_F0(TheOperators,"[]",yyvsp[-1].args);}
    break;


    }

/* Line 999 of yacc.c.  */
#line 2367 "lg.tab.cpp"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      /* Return failure if at end of input.  */
      if (yychar == YYEOF)
        {
	  /* Pop the error token.  */
          YYPOPSTACK;
	  /* Pop the rest of the stack.  */
	  while (yyss < yyssp)
	    {
	      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
	      yydestruct (yystos[*yyssp], yyvsp);
	      YYPOPSTACK;
	    }
	  YYABORT;
        }

      YYDSYMPRINTF ("Error: discarding", yytoken, &yylval, &yylloc);
      yydestruct (yytoken, &yylval);
      yychar = YYEMPTY;

    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*----------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action.  |
`----------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
      yydestruct (yystos[yystate], yyvsp);
      yyvsp--;
      yystate = *--yyssp;

      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


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

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 514 "lg.y"
 


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
extern void ShowAlloc(const char *s, size_t lg);
//extern void ShowNbAlloc(const char *s);
void init_lgfem() ;
void init_lgmesh() ;
void init_algo();
bool withrgraphique = false;
string StrVersionNumber();

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


 

