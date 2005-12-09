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
#line 76 "lg.y"
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
#line 264 "lg.tab.cpp"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 276 "lg.tab.cpp"

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
#define YYFINAL  71
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   771

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  71
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  42
/* YYNRULES -- Number of rules. */
#define YYNRULES  160
/* YYNRULES -- Number of states. */
#define YYNSTATES  347

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
       2,     2,     2,     2,     2,     2,     2,     2,    70,    66,
      16,     6,    17,    69,     2,     2,     2,     2,     2,     2,
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
     388,   392,   396,   400,   402,   408,   410,   414,   418,   422,
     426,   430,   434,   438,   442,   446,   450,   454,   458,   462,
     466,   470,   474,   478,   482,   486,   488,   490,   494,   500,
     501,   503,   507,   509,   513,   517,   523,   525,   529,   531,
     534,   536,   540,   544,   547,   549,   551,   553,   555,   557,
     562,   567,   574,   578,   582,   586,   591,   594,   597,   602,
     606
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      72,     0,    -1,    73,    46,    -1,    74,    -1,    99,    -1,
      74,    99,    -1,    -1,    77,    -1,    77,     6,   104,    -1,
      57,    77,    -1,    57,    12,    77,    -1,    80,    77,    -1,
      80,    12,    77,    -1,    35,    75,    38,    -1,    75,     5,
      77,    -1,    75,     5,    35,    75,    38,    -1,    75,     5,
      77,     6,   104,    -1,    75,     5,    57,    77,    -1,    75,
       5,    57,    12,    77,    -1,    75,     5,    80,    77,    -1,
      75,     5,    80,    12,    77,    -1,    77,    -1,    76,     5,
      77,    -1,    42,    -1,    57,    -1,    42,    -1,    42,     6,
     104,    -1,    42,    34,    79,    37,    -1,    78,     5,    78,
      -1,   105,    -1,    57,    42,    -1,    42,     6,   105,    -1,
      79,     5,   105,    -1,    79,     5,    77,     6,   105,    -1,
      55,    -1,    55,    35,    55,    38,    -1,    55,    35,    55,
       5,    55,    38,    -1,    55,    16,    55,    17,    -1,    42,
      -1,    42,    35,   105,    38,    -1,    42,     6,   105,    -1,
      35,    76,    38,    -1,    35,    76,    38,    35,   105,    38,
      -1,    35,    76,    38,     6,   105,    -1,    42,    34,   105,
      37,    -1,    35,    76,    38,    34,   105,    37,    -1,    57,
      -1,    57,    16,    55,    17,    -1,    82,    -1,    84,     5,
      82,    -1,    81,    -1,    85,     5,    81,    -1,    83,    85,
      -1,    83,    35,    55,    38,    84,    -1,    42,    34,    79,
      37,    -1,    87,    -1,    88,     5,    87,    -1,    -1,    80,
      90,    78,    66,    -1,    43,    88,    66,    -1,    86,    66,
      -1,    56,    42,     6,   102,    66,    -1,    -1,    56,    80,
      42,    34,    75,    37,    91,    67,    74,    68,    -1,    -1,
      56,    42,    34,    75,    37,    92,     6,   104,    66,    -1,
      67,    -1,    68,    -1,    50,    -1,    51,    -1,    -1,    80,
      98,    78,    -1,    66,    -1,    47,    45,    -1,    48,    45,
      -1,   102,    66,    -1,    89,    -1,    95,    34,   102,    66,
     102,    66,   102,    37,    99,    -1,    95,    34,    97,    66,
     102,    66,   102,    37,    99,    -1,    96,    34,   102,    37,
      99,    -1,     3,    34,   102,    37,    99,    -1,     3,    34,
     102,    37,    99,     4,    99,    -1,    93,    74,    94,    -1,
      63,    42,   101,    -1,    63,    42,    35,   109,    38,    66,
      -1,    52,    66,    -1,    53,    66,    -1,    54,   102,    66,
      -1,    34,    42,     6,   102,     5,   102,    37,    -1,   100,
      99,    -1,   104,    -1,   102,     5,   102,    -1,    21,    -1,
      20,    -1,    27,    -1,    29,    -1,    28,    -1,   105,    -1,
     105,     6,   104,    -1,   105,    58,   104,    -1,   105,    59,
     104,    -1,   105,    60,   104,    -1,   105,    61,   104,    -1,
     106,    -1,   106,    69,   105,    70,   105,    -1,   110,    -1,
     106,    22,   106,    -1,   106,    26,   106,    -1,   106,    25,
     106,    -1,   106,    23,   106,    -1,   106,    24,   106,    -1,
     106,    20,   106,    -1,   106,    21,   106,    -1,   106,     9,
     106,    -1,   106,     8,   106,    -1,   106,    12,   106,    -1,
     106,    13,   106,    -1,   106,    10,   106,    -1,   106,    11,
     106,    -1,   106,    16,   106,    -1,   106,    19,   106,    -1,
     106,    17,   106,    -1,   106,    18,   106,    -1,   106,    15,
     106,    -1,   106,    14,   106,    -1,   105,    -1,    70,    -1,
     105,    70,   105,    -1,   105,    70,   105,    70,   105,    -1,
      -1,    57,    -1,    77,     6,   105,    -1,   107,    -1,   108,
       5,    57,    -1,   108,     5,   107,    -1,   108,     5,    77,
       6,   105,    -1,   104,    -1,   109,     5,   104,    -1,   111,
      -1,   103,   111,    -1,   112,    -1,   112,    31,   110,    -1,
     112,    33,   110,    -1,   112,    32,    -1,    42,    -1,    39,
      -1,    40,    -1,    41,    -1,    45,    -1,   112,    34,   108,
      37,    -1,   112,    35,   107,    38,    -1,   112,    35,   107,
       5,   107,    38,    -1,   112,    35,    38,    -1,   112,    36,
      42,    -1,    57,    36,    42,    -1,    57,    34,   108,    37,
      -1,   112,    29,    -1,   112,    28,    -1,    55,    34,   102,
      37,    -1,    34,   102,    37,    -1,    35,   109,    38,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   195,   195,   231,   234,   235,   238,   239,   240,   241,
     242,   243,   244,   245,   246,   247,   248,   249,   250,   251,
     252,   255,   256,   259,   259,   261,   262,   263,   265,   272,
     273,   274,   275,   276,   279,   280,   281,   282,   289,   290,
     291,   292,   293,   294,   297,   298,   302,   303,   308,   309,
     311,   312,   314,   315,   319,   323,   324,   327,   327,   328,
     329,   330,   332,   331,   347,   346,   355,   356,   358,   360,
     365,   365,   368,   369,   370,   371,   372,   373,   374,   378,
     379,   380,   381,   383,   385,   388,   392,   396,   404,   410,
     416,   417,   422,   423,   424,   425,   426,   430,   431,   432,
     433,   434,   435,   440,   441,   444,   445,   446,   447,   448,
     449,   450,   451,   452,   453,   454,   455,   456,   457,   458,
     459,   460,   461,   462,   463,   468,   469,   470,   471,   474,
     475,   476,   477,   478,   479,   480,   483,   484,   488,   489,
     492,   493,   494,   495,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   522,
     523
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
  "SOLVE", "';'", "'{'", "'}'", "'?'", "':'", "$accept", "start", "input", 
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
     295,   296,   297,   298,   299,   300,    59,   123,   125,    63,
      58
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    71,    72,    73,    74,    74,    75,    75,    75,    75,
      75,    75,    75,    75,    75,    75,    75,    75,    75,    75,
      75,    76,    76,    77,    77,    78,    78,    78,    78,    79,
      79,    79,    79,    79,    80,    80,    80,    80,    81,    81,
      81,    81,    81,    81,    82,    82,    83,    83,    84,    84,
      85,    85,    86,    86,    87,    88,    88,    90,    89,    89,
      89,    89,    91,    89,    92,    89,    93,    94,    95,    96,
      98,    97,    99,    99,    99,    99,    99,    99,    99,    99,
      99,    99,    99,    99,    99,    99,    99,    99,   100,   101,
     102,   102,   103,   103,   103,   103,   103,   104,   104,   104,
     104,   104,   104,   105,   105,   106,   106,   106,   106,   106,
     106,   106,   106,   106,   106,   106,   106,   106,   106,   106,
     106,   106,   106,   106,   106,   107,   107,   107,   107,   108,
     108,   108,   108,   108,   108,   108,   109,   109,   110,   110,
     111,   111,   111,   111,   112,   112,   112,   112,   112,   112,
     112,   112,   112,   112,   112,   112,   112,   112,   112,   112,
     112
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
       3,     3,     3,     1,     5,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     1,     1,     3,     5,     0,
       1,     3,     1,     3,     3,     5,     1,     3,     1,     2,
       1,     3,     3,     2,     1,     1,     1,     1,     1,     4,
       4,     6,     3,     3,     3,     4,     2,     2,     4,     3,
       3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       0,     0,    93,    92,    94,    96,    95,     0,     0,   145,
     146,   147,   144,     0,   148,     0,     0,    68,    69,     0,
       0,     0,    34,     0,    46,     0,    72,    66,     0,     0,
       3,    57,     0,     0,    76,     0,     0,     0,     4,     0,
       0,    90,    97,   103,   105,   138,   140,     0,     0,     0,
       0,   136,     0,     0,    55,     0,    73,    74,    85,    86,
       0,     0,     0,     0,     0,    34,     0,     0,   129,     0,
       0,     1,     2,     5,     0,     0,    38,    50,    52,    60,
       0,     0,     0,     0,    75,   139,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   157,   156,     0,   143,     0,   129,     0,     0,     0,
     159,     0,   160,     0,     0,    59,    87,     0,     0,     0,
       0,     6,     0,     0,   144,   130,   126,     0,   125,   132,
       0,   154,     0,     0,     0,    83,    25,     0,    23,     0,
      24,     0,    21,     0,     0,     0,    67,    82,    70,     0,
       0,     0,    91,    98,    99,   100,   101,   102,   114,   113,
     117,   118,   115,   116,   124,   123,   119,   121,   122,   120,
     111,   112,   106,   109,   110,   108,   107,     0,   141,   142,
       0,   152,     0,   153,     0,   137,   144,     0,     0,    29,
      56,    37,   158,     0,    35,     0,     6,    24,     0,     7,
       0,     6,    47,     0,     0,     0,   155,     0,     0,    89,
       0,     0,     0,    58,     0,     0,    41,    40,     0,     0,
      51,     0,     0,     0,     0,     0,   149,     0,   150,    80,
       0,    30,     0,    54,     0,    61,     0,     0,     9,     0,
      64,     0,     0,    11,     0,   131,   127,   133,     0,   134,
       0,     0,    26,     0,    28,     0,     0,    48,    53,    22,
       0,     0,    39,    71,     0,     0,    79,   104,     0,     0,
      31,    24,     0,    32,    36,    13,    10,     6,    24,    14,
       0,     0,     8,    12,    62,     0,     0,     0,    84,    27,
       0,     0,     0,    43,     0,     0,     0,   151,    81,     0,
       0,     0,    17,     0,     0,    19,     0,     0,   128,   135,
       0,     0,     0,    49,    42,     0,     0,    33,    15,    18,
      16,    20,     0,     0,    91,     0,    44,     0,     0,    65,
       0,    88,     0,    78,    77,    63,    45
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,    28,    29,    30,   208,   151,   209,   147,   198,    31,
      77,   267,    32,   268,    78,    33,    54,    55,    34,    74,
     317,   291,    35,   157,    36,    37,   159,   231,    38,   144,
     145,    39,    40,    41,    42,    43,   139,   140,    52,    44,
      45,    46
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -202
static const short yypact[] =
{
     408,   -28,  -202,  -202,  -202,  -202,  -202,   249,   249,  -202,
    -202,  -202,  -202,     0,  -202,   -32,     4,  -202,  -202,    -8,
      18,   249,   135,    12,   140,    36,  -202,  -202,   150,    75,
     408,  -202,   148,    73,  -202,   408,   144,   169,  -202,     2,
     173,  -202,    58,     8,  -202,  -202,   226,   249,   189,   113,
     100,  -202,    41,   193,  -202,     6,  -202,  -202,  -202,  -202,
       7,   105,   249,   129,    35,    27,   190,   176,   481,   192,
     165,  -202,  -202,  -202,   195,   124,   136,  -202,   234,  -202,
     298,   577,   249,   249,  -202,  -202,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,  -202,  -202,   249,  -202,   249,   481,   449,   200,   101,
    -202,   249,  -202,   608,     0,  -202,  -202,   246,   103,    42,
     249,   162,   217,   247,   259,   123,  -202,   260,   198,  -202,
     104,  -202,   229,   249,   408,  -202,   139,    30,  -202,   235,
    -202,    43,  -202,   249,   249,   153,  -202,  -202,  -202,   206,
      31,   109,  -202,  -202,  -202,  -202,  -202,  -202,   706,   706,
     721,   721,   734,   734,   745,   745,   665,   665,   665,   665,
     222,   222,  -202,  -202,  -202,  -202,  -202,   204,  -202,  -202,
     115,  -202,    47,  -202,   408,  -202,   274,   151,   121,  -202,
    -202,  -202,  -202,   227,  -202,    32,   162,     3,   125,   279,
      53,   162,  -202,   249,   249,   513,  -202,   280,    65,  -202,
     249,   608,   195,  -202,   154,    91,   137,  -202,   254,    91,
    -202,   195,   249,   249,   408,   249,  -202,   545,  -202,   283,
     249,  -202,   639,  -202,   255,  -202,    66,    91,  -202,   167,
    -202,   249,    91,  -202,   126,  -202,   225,   123,   290,  -202,
     249,   231,  -202,   127,  -202,    91,   269,  -202,   302,  -202,
     249,   249,  -202,   305,    33,    34,  -202,  -202,   277,   408,
    -202,   113,   306,  -202,  -202,  -202,  -202,   162,    70,   310,
      71,   314,  -202,  -202,  -202,   249,   249,   316,  -202,  -202,
      87,   249,   154,  -202,   284,   249,   249,  -202,  -202,   249,
      96,    91,  -202,   249,    91,  -202,   249,   256,  -202,  -202,
     249,   294,   287,  -202,  -202,   130,   131,  -202,  -202,  -202,
    -202,  -202,   263,   408,   293,   249,  -202,   408,   408,  -202,
     353,  -202,   297,  -202,  -202,  -202,  -202
};

/* YYPGOTO[NTERM-NUM].  */
static const short yypgoto[] =
{
    -202,  -202,  -202,   -33,  -201,    77,    -9,  -178,   110,   -20,
     180,    45,  -202,  -202,  -202,  -202,   212,  -202,  -202,  -202,
    -202,  -202,  -202,  -202,  -202,  -202,  -202,  -202,   -29,  -202,
    -202,    -7,  -202,     1,   -60,   606,  -113,   228,   214,    79,
     318,  -202
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -25
static const short yytable[] =
{
      50,    73,    80,    66,   192,   246,    47,    83,   138,    51,
     254,   124,    83,    56,    60,   247,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   222,    83,    83,    83,    83,
     119,   130,    53,    61,   264,   148,   121,   203,   225,    57,
     187,    73,   237,   273,    64,   128,   138,   138,    58,   137,
     150,   158,    63,   199,    86,   252,   152,    65,    84,   131,
     121,   249,   125,   126,   160,   161,   162,   110,    70,   122,
     204,   226,   311,   314,    59,   238,   310,   163,   164,   165,
     166,   167,   225,   227,   228,   148,   223,   233,   245,   305,
     306,   249,   259,   261,   285,    83,    83,   137,    83,   215,
     150,   210,   148,   148,    83,   219,    87,    88,    89,    90,
     215,    72,   195,   205,   278,   321,   242,   150,   150,   -24,
     249,   249,   242,   148,   328,    83,    83,   120,   194,    79,
     202,   216,   153,   270,    51,   220,   234,    68,   150,    69,
      71,    61,   236,   255,   256,   138,    67,    68,   243,    69,
     127,   199,   250,   294,   299,   239,   148,   337,   338,    62,
      63,   154,   271,   221,    68,   277,    69,   138,    81,   149,
     280,   150,   283,    75,   129,    68,   210,    69,   229,   265,
      76,   210,   188,   241,   189,    76,   266,   206,   248,   142,
     143,   253,   287,    82,   148,   276,   258,     7,     8,   148,
     303,   304,     9,    10,    11,    12,   269,    65,    14,   207,
     152,   262,    65,    62,   288,   274,   275,   123,    48,   290,
      49,   133,   132,   282,   141,   318,   319,   146,   286,   155,
     289,   322,   193,   293,   105,   106,   107,   108,   109,   327,
     308,   211,   292,   297,   111,   112,   152,   113,   114,   115,
     116,   117,   118,   201,   212,   -23,   213,   210,   214,     2,
       3,   217,   232,   224,   235,   342,     4,     5,     6,   312,
     240,   315,   244,     7,     8,   251,   260,   279,     9,    10,
      11,    12,   272,   284,    14,   295,   296,   298,   325,   326,
     340,     1,   329,   301,    48,   331,    49,   302,   343,   344,
     222,    73,   309,   334,   330,   307,   313,   332,     2,     3,
     316,   320,   324,   333,   336,     4,     5,     6,   335,   339,
     341,   263,     7,     8,   346,   230,   200,     9,    10,    11,
      12,    13,   300,    14,   190,    15,    16,   323,    17,    18,
      19,    20,    21,    22,    23,    24,     1,   218,    85,     0,
       0,    25,     0,     0,    26,    27,   156,     0,     0,     0,
       0,     0,     0,     2,     3,     0,     0,     0,     0,     0,
       4,     5,     6,     0,     0,     0,     0,     7,     8,     0,
       0,     0,     9,    10,    11,    12,    13,     0,    14,     0,
      15,    16,     0,    17,    18,    19,    20,    21,    22,    23,
      24,     1,     0,     0,     0,     0,    25,     0,     0,    26,
      27,   345,     0,     0,     0,     0,     0,     0,     2,     3,
       0,     0,     0,     0,     0,     4,     5,     6,     0,     0,
       0,     0,     7,     8,     0,     0,     0,     9,    10,    11,
      12,    13,     0,    14,     0,    15,    16,     0,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     2,
       3,    25,     0,     0,    26,    27,     4,     5,     6,     0,
       0,     0,     0,     7,     8,     0,     0,   191,     9,    10,
      11,    12,     0,     0,    14,     0,     0,     0,     0,     0,
       0,     2,     3,     0,    48,     0,    49,     0,     4,     5,
       6,     0,     0,     0,     0,     7,     8,     0,     0,   136,
       9,    10,    11,   134,     0,     0,    14,     0,     0,     0,
       0,     0,     0,     2,     3,     0,    48,     0,   135,     0,
       4,     5,     6,     0,     0,     0,     0,     7,     8,     0,
       0,   136,     9,    10,    11,   134,     0,     0,    14,     0,
       0,     0,     0,     0,     0,     2,     3,     0,    48,     0,
     257,     0,     4,     5,     6,     0,     0,     0,     0,     7,
       8,     0,     0,   136,     9,    10,    11,    12,     0,     0,
      14,     0,     0,     0,     0,     0,     0,     2,     3,     0,
      48,     0,    49,     0,     4,     5,     6,     0,     0,     0,
       0,     7,     8,     0,     0,   136,     9,    10,    11,    12,
       0,     0,    14,     0,     0,     0,     0,     0,     2,     3,
       0,     0,    22,     0,    49,     4,     5,     6,     0,     0,
       0,     0,     7,     8,     0,     0,     0,     9,    10,    11,
     196,     0,     0,    14,     0,     0,     0,     0,     0,     2,
       3,     0,     0,    48,     0,   197,     4,     5,     6,     0,
       0,     0,     0,     7,     8,     0,     0,     0,     9,    10,
      11,   134,     0,     0,    14,   103,   104,   105,   106,   107,
     108,   109,     0,     0,    48,     0,   281,   168,   169,   170,
     171,   172,   173,   174,   175,   176,   177,   178,   179,   180,
     181,   182,   183,   184,   185,   186,    93,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
     107,   108,   109,    95,    96,    97,    98,    99,   100,   101,
     102,   103,   104,   105,   106,   107,   108,   109,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,    99,   100,   101,   102,   103,   104,   105,   106,   107,
     108,   109
};

static const short yycheck[] =
{
       7,    30,    35,    23,   117,   206,    34,     5,    68,     8,
     211,     5,     5,    45,    21,    12,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,     5,     5,     5,     5,     5,
      47,     6,    42,    16,   222,    42,     5,     5,     5,    45,
     110,    80,     5,   231,    42,    62,   116,   117,    66,    68,
      57,    81,    35,   123,     6,    12,    75,    55,    66,    34,
       5,     5,    66,    66,    81,    82,    83,    69,    42,    38,
      38,    38,    12,    12,    66,    38,   287,    86,    87,    88,
      89,    90,     5,   153,   154,    42,    66,    66,    66,    66,
      66,     5,   215,    38,    38,     5,     5,   116,     5,     5,
      57,   131,    42,    42,     5,   144,    58,    59,    60,    61,
       5,    46,   121,   130,   237,    38,     5,    57,    57,     6,
       5,     5,     5,    42,    38,     5,     5,    37,    37,    66,
      37,    37,     6,     6,   143,     6,    37,    34,    57,    36,
       0,    16,    37,   213,   214,   215,    16,    34,    37,    36,
      55,   221,    37,    37,    37,   194,    42,    37,    37,    34,
      35,    35,    35,    34,    34,   235,    36,   237,    34,    55,
     240,    57,   242,    35,    55,    34,   206,    36,    35,    35,
      42,   211,   113,    42,   115,    42,    42,    35,   207,    34,
      35,   210,    35,    34,    42,   234,   215,    34,    35,    42,
     270,   271,    39,    40,    41,    42,   225,    55,    45,    57,
     229,   220,    55,    34,    57,   232,   233,    34,    55,   249,
      57,    55,    42,   242,    42,   295,   296,    42,   247,     5,
     249,   301,    42,   252,    22,    23,    24,    25,    26,   309,
     279,    34,   251,   260,    28,    29,   265,    31,    32,    33,
      34,    35,    36,    17,    17,     6,     6,   287,    70,    20,
      21,    42,    66,    38,    70,   335,    27,    28,    29,   288,
       6,   290,    55,    34,    35,     6,     6,     4,    39,    40,
      41,    42,    38,    38,    45,    70,     6,    66,   305,   306,
     333,     3,   311,    34,    55,   314,    57,     5,   337,   338,
       5,   340,     6,   320,   313,    38,     6,   316,    20,    21,
       6,     5,    38,    67,    37,    27,    28,    29,    34,    66,
      37,   221,    34,    35,    37,   155,   124,    39,    40,    41,
      42,    43,   265,    45,   116,    47,    48,   302,    50,    51,
      52,    53,    54,    55,    56,    57,     3,   143,    40,    -1,
      -1,    63,    -1,    -1,    66,    67,    68,    -1,    -1,    -1,
      -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    -1,    -1,
      27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    43,    -1,    45,    -1,
      47,    48,    -1,    50,    51,    52,    53,    54,    55,    56,
      57,     3,    -1,    -1,    -1,    -1,    63,    -1,    -1,    66,
      67,    68,    -1,    -1,    -1,    -1,    -1,    -1,    20,    21,
      -1,    -1,    -1,    -1,    -1,    27,    28,    29,    -1,    -1,
      -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
      42,    43,    -1,    45,    -1,    47,    48,    -1,    50,    51,
      52,    53,    54,    55,    56,    57,    -1,    -1,    -1,    20,
      21,    63,    -1,    -1,    66,    67,    27,    28,    29,    -1,
      -1,    -1,    -1,    34,    35,    -1,    -1,    38,    39,    40,
      41,    42,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    20,    21,    -1,    55,    -1,    57,    -1,    27,    28,
      29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    70,
      39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    20,    21,    -1,    55,    -1,    57,    -1,
      27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,
      -1,    70,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    -1,    20,    21,    -1,    55,    -1,
      57,    -1,    27,    28,    29,    -1,    -1,    -1,    -1,    34,
      35,    -1,    -1,    70,    39,    40,    41,    42,    -1,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    20,    21,    -1,
      55,    -1,    57,    -1,    27,    28,    29,    -1,    -1,    -1,
      -1,    34,    35,    -1,    -1,    70,    39,    40,    41,    42,
      -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,    20,    21,
      -1,    -1,    55,    -1,    57,    27,    28,    29,    -1,    -1,
      -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
      42,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,    20,
      21,    -1,    -1,    55,    -1,    57,    27,    28,    29,    -1,
      -1,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,
      41,    42,    -1,    -1,    45,    20,    21,    22,    23,    24,
      25,    26,    -1,    -1,    55,    -1,    57,    91,    92,    93,
      94,    95,    96,    97,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     3,    20,    21,    27,    28,    29,    34,    35,    39,
      40,    41,    42,    43,    45,    47,    48,    50,    51,    52,
      53,    54,    55,    56,    57,    63,    66,    67,    72,    73,
      74,    80,    83,    86,    89,    93,    95,    96,    99,   102,
     103,   104,   105,   106,   110,   111,   112,    34,    55,    57,
     102,   104,   109,    42,    87,    88,    45,    45,    66,    66,
     102,    16,    34,    35,    42,    55,    80,    16,    34,    36,
      42,     0,    46,    99,    90,    35,    42,    81,    85,    66,
      74,    34,    34,     5,    66,   111,     6,    58,    59,    60,
      61,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      69,    28,    29,    31,    32,    33,    34,    35,    36,   102,
      37,     5,    38,    34,     5,    66,    66,    55,   102,    55,
       6,    34,    42,    55,    42,    57,    70,    77,   105,   107,
     108,    42,    34,    35,   100,   101,    42,    78,    42,    55,
      57,    76,    77,     6,    35,     5,    68,    94,    80,    97,
     102,   102,   102,   104,   104,   104,   104,   104,   106,   106,
     106,   106,   106,   106,   106,   106,   106,   106,   106,   106,
     106,   106,   106,   106,   106,   106,   106,   105,   110,   110,
     108,    38,   107,    42,    37,   104,    42,    57,    79,   105,
      87,    17,    37,     5,    38,   102,    35,    57,    75,    77,
      80,    34,    17,     6,    70,     5,    37,    42,   109,    99,
       6,    34,     5,    66,    38,     5,    38,   105,   105,    35,
      81,    98,    66,    66,    37,    70,    37,     5,    38,    99,
       6,    42,     5,    37,    55,    66,    75,    12,    77,     5,
      37,     6,    12,    77,    75,   105,   105,    57,    77,   107,
       6,    38,   104,    79,    78,    35,    42,    82,    84,    77,
       6,    35,    38,    78,   102,   102,    99,   105,   107,     4,
     105,    57,    77,   105,    38,    38,    77,    35,    57,    77,
      80,    92,   104,    77,    37,    70,     6,   102,    66,    37,
      76,    34,     5,   105,   105,    66,    66,    38,    99,     6,
      75,    12,    77,     6,    12,    77,     6,    91,   105,   105,
       5,    38,   105,    82,    38,   102,   102,   105,    38,    77,
     104,    77,   104,    67,   102,    34,    37,    37,    37,    66,
      74,    37,   105,    99,    99,    68,    37
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
  return 0;;}
    break;

  case 4:
#line 234 "lg.y"
    {yyval.cinst=yyvsp[0].cexp;;;;}
    break;

  case 5:
#line 235 "lg.y"
    { yyval.cinst= (yyvsp[-1].cinst+=yyvsp[0].cexp) ;}
    break;

  case 6:
#line 238 "lg.y"
    { yyval.clist_id=new ListOfId();;}
    break;

  case 7:
#line 239 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str));}
    break;

  case 8:
#line 240 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;}
    break;

  case 9:
#line 241 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>()));}
    break;

  case 10:
#line 242 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true));}
    break;

  case 11:
#line 243 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;}
    break;

  case 12:
#line 244 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;}
    break;

  case 13:
#line 245 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;}
    break;

  case 14:
#line 246 "lg.y"
    { yyval.clist_id = yyvsp[-2].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str)) ;}
    break;

  case 15:
#line 247 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;}
    break;

  case 16:
#line 248 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;}
    break;

  case 17:
#line 249 "lg.y"
    { yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>())) ;}
    break;

  case 18:
#line 250 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true)) ;}
    break;

  case 19:
#line 251 "lg.y"
    { yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;}
    break;

  case 20:
#line 252 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;}
    break;

  case 21:
#line 255 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;}
    break;

  case 22:
#line 256 "lg.y"
    { yyval.clist_id=yyvsp[-2].clist_id  ; yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;}
    break;

  case 25:
#line 261 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[0].str,dcltype);}
    break;

  case 26:
#line 262 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-2].str,dcltype,yyvsp[0].cexp);}
    break;

  case 27:
#line 263 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-3].str,dcltype,yyvsp[-1].args);
                                              yyvsp[-1].args.destroy();}
    break;

  case 28:
#line 265 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 29:
#line 272 "lg.y"
    {yyval.args=yyvsp[0].cexp;}
    break;

  case 30:
#line 273 "lg.y"
    {yyval.args=Find(yyvsp[-1].str);}
    break;

  case 31:
#line 274 "lg.y"
    { yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);}
    break;

  case 32:
#line 275 "lg.y"
    { yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;}
    break;

  case 33:
#line 276 "lg.y"
    { yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp));}
    break;

  case 35:
#line 280 "lg.y"
    {yyval.type=TypeArray(yyvsp[-3].type,yyvsp[-1].type);}
    break;

  case 36:
#line 281 "lg.y"
    {yyval.type=TypeArray(yyvsp[-5].type,yyvsp[-3].type,yyvsp[-1].type);}
    break;

  case 37:
#line 282 "lg.y"
    {yyval.type=TypeTemplate(yyvsp[-3].type,yyvsp[-1].type);}
    break;

  case 38:
#line 289 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[0].str,currentblock,fespacetype,fespacecomplex); ;}
    break;

  case 39:
#line 290 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;}
    break;

  case 40:
#line 291 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[-2].str,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;}
    break;

  case 41:
#line 292 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[-1].clist_id,currentblock,fespacetype,fespacecomplex) ;}
    break;

  case 42:
#line 293 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;}
    break;

  case 43:
#line 294 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[-3].clist_id,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;}
    break;

  case 44:
#line 297 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;}
    break;

  case 45:
#line 298 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;}
    break;

  case 46:
#line 302 "lg.y"
    {fespacecomplex=false;  fespacetype = Find(yyvsp[0].str);;}
    break;

  case 47:
#line 303 "lg.y"
    {
             if (yyvsp[-1].type != typevarreal && yyvsp[-1].type != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=(yyvsp[-1].type==typevarcomplex);
             fespacetype = Find(yyvsp[-3].str);;}
    break;

  case 48:
#line 308 "lg.y"
    {  yyval.cexp = yyvsp[0].cexp  ;}
    break;

  case 49:
#line 309 "lg.y"
    { yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;}
    break;

  case 50:
#line 311 "lg.y"
    {  yyval.cexp = yyvsp[0].cexp  ;}
    break;

  case 51:
#line 312 "lg.y"
    { yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;}
    break;

  case 52:
#line 314 "lg.y"
    { yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;}
    break;

  case 53:
#line 315 "lg.y"
    { yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;}
    break;

  case 54:
#line 320 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariableFES,size_t>(yyvsp[-3].str,atype<pfes*>(),yyvsp[-1].args,dimFESpaceImage(yyvsp[-1].args));
     yyvsp[-1].args.destroy(); ;}
    break;

  case 56:
#line 324 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 57:
#line 327 "lg.y"
    {dcltype=yyvsp[0].type;}
    break;

  case 58:
#line 327 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 59:
#line 328 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 60:
#line 329 "lg.y"
    { yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 61:
#line 330 "lg.y"
    {yyval.cexp=currentblock->NewID(yyvsp[-4].type,yyvsp[-3].str,yyvsp[-1].cexp);;}
    break;

  case 62:
#line 332 "lg.y"
    {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = yyvsp[-4].type->right();
                      routineinblock[kkembtype] = currentblock;
                      yyvsp[-1].routine=new Routine(yyvsp[-5].type,yyvsp[-4].type->right(),yyvsp[-3].str,yyvsp[-1].clist_id,currentblock);
                     // cout << " \n after new routine \n " << endl;                      
                      ;}
    break;

  case 63:
#line 340 "lg.y"
    { currentblock=yyvsp[-5].routine->Set(yyvsp[-1].cinst);
                       currentblock->Add(yyvsp[-7].str,"(",yyvsp[-5].routine);
                       kkembtype--;
                       yyval.cexp=0;
                    
                        ;}
    break;

  case 64:
#line 347 "lg.y"
    {currentblock = new Block(currentblock); yyvsp[-4].type->SetArgs(yyvsp[-1].clist_id);;}
    break;

  case 65:
#line 349 "lg.y"
    {  yyval.cinst=currentblock->close(currentblock);
                         yyval.cexp=currentblock->NewID(yyvsp[-8].type,yyvsp[-7].str,yyvsp[-1].cexp,*yyvsp[-5].clist_id);
                         delete yyvsp[-5].clist_id; //  FH 23032005
                         ;}
    break;

  case 66:
#line 355 "lg.y"
    {  currentblock = new Block(currentblock);}
    break;

  case 67:
#line 356 "lg.y"
    {  yyval.cexp=currentblock->close(currentblock);}
    break;

  case 68:
#line 358 "lg.y"
    {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 69:
#line 360 "lg.y"
    {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 70:
#line 365 "lg.y"
    {dcltype=yyvsp[0].type;currentblock = new Block(currentblock);}
    break;

  case 71:
#line 366 "lg.y"
    {yyval.cexp=yyvsp[0].cexp;}
    break;

  case 72:
#line 368 "lg.y"
    {yyval.cexp=0;;}
    break;

  case 73:
#line 369 "lg.y"
    {zzzfff->input(yyvsp[0].str);yyval.cexp= 0; ;}
    break;

  case 74:
#line 370 "lg.y"
    {load(yyvsp[0].str);yyval.cexp= 0; ;}
    break;

  case 75:
#line 371 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 76:
#line 372 "lg.y"
    {yyval.cexp=yyvsp[0].cexp;}
    break;

  case 77:
#line 373 "lg.y"
    {inloopcount--; yyval.cexp=For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 78:
#line 375 "lg.y"
    {inloopcount--; 
                yyval.cexp=C_F0(For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp),currentblock->close(currentblock));}
    break;

  case 79:
#line 378 "lg.y"
    {inloopcount--;yyval.cexp=While(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 80:
#line 379 "lg.y"
    {yyval.cexp=FIf(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 81:
#line 380 "lg.y"
    {yyval.cexp=FIf(yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 82:
#line 381 "lg.y"
    { 
                      yyval.cexp=C_F0(new E_block(yyvsp[-1].cinst,yyvsp[0].cexp),atype<void>()) ;}
    break;

  case 83:
#line 383 "lg.y"
    {
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-1].str,C_F0(TheOperators,"[border]",yyvsp[0].args));}
    break;

  case 84:
#line 385 "lg.y"
    {
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-4].str,C_F0(TheOperators,"[border]",yyvsp[-2].args));}
    break;

  case 85:
#line 388 "lg.y"
    {
                    if(inloopcount) 
                      yyval.cexp= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;}
    break;

  case 86:
#line 392 "lg.y"
    { 
                    if(inloopcount)
                        yyval.cexp= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");}
    break;

  case 87:
#line 396 "lg.y"
    { 
                    if (kkembtype>=0)
                      yyval.cexp= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo(yyvsp[-1].cexp)) ,atype<void>());
                     else lgerror(" return not in routine ") ;}
    break;

  case 88:
#line 404 "lg.y"
    { 
   currentblock = new Block(currentblock);
   yyval.args = currentblock->NewVar<LocalVariable>(yyvsp[-5].str,atype<double*>());
   yyval.args+= yyvsp[-3].cexp;
   yyval.args+= yyvsp[-1].cexp ;}
    break;

  case 89:
#line 410 "lg.y"
    {   
   yyval.args = (yyvsp[-1].args += yyvsp[0].cexp);
   currentblock->close(currentblock);}
    break;

  case 91:
#line 417 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);;}
    break;

  case 98:
#line 431 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 99:
#line 432 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"+=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 100:
#line 433 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"-=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 101:
#line 434 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"*=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 102:
#line 435 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"/=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 104:
#line 441 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"?:",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 106:
#line 445 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 107:
#line 446 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 108:
#line 447 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 109:
#line 448 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 110:
#line 449 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 111:
#line 450 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 112:
#line 451 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 113:
#line 452 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 114:
#line 453 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 115:
#line 454 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 116:
#line 455 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 117:
#line 456 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 118:
#line 457 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 119:
#line 458 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 120:
#line 459 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 121:
#line 460 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 122:
#line 461 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 123:
#line 462 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 124:
#line 463 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 125:
#line 468 "lg.y"
    {yyval.cexp=yyvsp[0].cexp;}
    break;

  case 126:
#line 469 "lg.y"
    {yyval.cexp=C_F0(TheOperators,":");}
    break;

  case 127:
#line 470 "lg.y"
    {yyval.cexp=C_F0(TheOperators,":",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 128:
#line 471 "lg.y"
    {yyval.cexp=C_F0(TheOperators,":",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 129:
#line 474 "lg.y"
    {yyval.args=0;}
    break;

  case 130:
#line 475 "lg.y"
    {yyval.args=Find(yyvsp[0].str);}
    break;

  case 131:
#line 476 "lg.y"
    { yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);}
    break;

  case 132:
#line 477 "lg.y"
    {yyval.args=yyvsp[0].cexp;}
    break;

  case 133:
#line 478 "lg.y"
    { yyval.args = (yyvsp[-2].args += Find(yyvsp[0].str)) ;}
    break;

  case 134:
#line 479 "lg.y"
    { yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;}
    break;

  case 135:
#line 480 "lg.y"
    { yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp)) ;}
    break;

  case 136:
#line 483 "lg.y"
    {yyval.args=yyvsp[0].cexp;}
    break;

  case 137:
#line 484 "lg.y"
    {yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;}
    break;

  case 139:
#line 489 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[0].cexp);}
    break;

  case 141:
#line 493 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 142:
#line 494 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 143:
#line 495 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[0].oper,yyvsp[-1].cexp);}
    break;

  case 144:
#line 499 "lg.y"
    {yyval.cexp=Find(yyvsp[0].str);;}
    break;

  case 145:
#line 500 "lg.y"
    {yyval.cexp= CConstant(yyvsp[0].lnum);}
    break;

  case 146:
#line 501 "lg.y"
    {yyval.cexp= CConstant(yyvsp[0].dnum);}
    break;

  case 147:
#line 502 "lg.y"
    {yyval.cexp= CConstant(complex<double>(0,yyvsp[0].dnum));}
    break;

  case 148:
#line 503 "lg.y"
    {yyval.cexp= CConstant<const char *>(yyvsp[0].str);}
    break;

  case 149:
#line 504 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].args);;}
    break;

  case 150:
#line 505 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].cexp);}
    break;

  case 151:
#line 506 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-5].cexp,yyvsp[-4].oper,yyvsp[-3].cexp,yyvsp[-1].cexp);}
    break;

  case 152:
#line 507 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-2].cexp,"[]");}
    break;

  case 153:
#line 508 "lg.y"
    { yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].str) ;;}
    break;

  case 154:
#line 509 "lg.y"
    { yyval.cexp=C_F0(Find(yyvsp[-2].str),yyvsp[0].str) ;;}
    break;

  case 155:
#line 510 "lg.y"
    { yyval.cexp=C_F0(Find(yyvsp[-3].str),yyvsp[-2].oper,yyvsp[-1].args) ;;}
    break;

  case 156:
#line 511 "lg.y"
    {yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);}
    break;

  case 157:
#line 512 "lg.y"
    {yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);}
    break;

  case 158:
#line 513 "lg.y"
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

  case 159:
#line 522 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 160:
#line 523 "lg.y"
    { yyval.cexp=C_F0(TheOperators,"[]",yyvsp[-1].args);}
    break;


    }

/* Line 999 of yacc.c.  */
#line 2348 "lg.tab.cpp"

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


#line 528 "lg.y"
 


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


 

