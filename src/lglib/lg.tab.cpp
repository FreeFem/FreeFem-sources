/* A Bison parser, made by GNU Bison 2.0.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

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

/* Substitute the variable and function names.  */
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
     TRY = 290,
     CATCH = 291,
     THROW = 292,
     TYPE = 293,
     FUNCTION = 294,
     FESPACE = 295,
     FESPACE1 = 296,
     FESPACE3 = 297,
     PLUSEQ = 298,
     MOINSEQ = 299,
     MULEQ = 300,
     DIVEQ = 301,
     ARROW = 302,
     BORDER = 303,
     CURVE = 304,
     SOLVE = 305
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
#define TRY 290
#define CATCH 291
#define THROW 292
#define TYPE 293
#define FUNCTION 294
#define FESPACE 295
#define FESPACE1 296
#define FESPACE3 297
#define PLUSEQ 298
#define MOINSEQ 299
#define MULEQ 300
#define DIVEQ 301
#define ARROW 302
#define BORDER 303
#define CURVE 304
#define SOLVE 305




/* Copy the first part of user declarations.  */
#line 1 "lg.y"
 
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
  // for reset cout,cin  in windows  dll
#ifdef _WIN32
#include <ext/stdio_filebuf.h>
#include <iostream>
#include <cstdio>
#endif

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
#include "MeshPoint.hpp"

#include "lgfem.hpp" 
#include "lex.hpp"
#include "environment.hpp"

    extern FILE *ThePlotStream;
    
class Routine;
bool load(string s);

 template <class R,int d> class FE;
 template <class R,int d,int i> class FE_;

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
int fespacedim;

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
#line 118 "lg.y"
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
/* ListCatch * clist_Catchs;*/
} YYSTYPE;
/* Line 190 of yacc.c.  */
#line 316 "lg.tab.cpp"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 213 of yacc.c.  */
#line 328 "lg.tab.cpp"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

# ifndef YYFREE
#  define YYFREE free
# endif
# ifndef YYMALLOC
#  define YYMALLOC malloc
# endif

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   else
#    define YYSTACK_ALLOC alloca
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
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short int yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short int) + sizeof (YYSTYPE))			\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
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
   typedef short int yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  81
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   909

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  76
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  45
/* YYNRULES -- Number of rules. */
#define YYNRULES  182
/* YYNRULES -- Number of states. */
#define YYNSTATES  395

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   305

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
       2,     2,     2,     2,     2,     2,     2,     2,    75,    71,
      16,     6,    17,    74,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    35,     2,    38,    31,    33,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    72,    10,    73,     2,     2,     2,     2,
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
      65,    66,    67,    68,    69,    70
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short int yyprhs[] =
{
       0,     0,     3,     6,     8,    10,    13,    14,    16,    20,
      23,    27,    30,    34,    37,    41,    45,    49,    55,    61,
      66,    72,    77,    83,    88,    94,    96,   100,   102,   104,
     106,   108,   110,   114,   119,   123,   125,   128,   131,   134,
     138,   142,   148,   150,   155,   162,   167,   169,   174,   178,
     182,   189,   195,   200,   207,   209,   211,   213,   215,   220,
     222,   226,   228,   232,   235,   241,   246,   248,   252,   253,
     258,   262,   265,   271,   272,   283,   284,   294,   296,   298,
     300,   302,   303,   307,   309,   311,   314,   317,   323,   326,
     328,   338,   348,   354,   360,   368,   372,   376,   383,   386,
     389,   393,   401,   409,   412,   414,   418,   420,   422,   424,
     426,   428,   430,   434,   438,   442,   446,   450,   452,   458,
     460,   464,   468,   472,   476,   480,   484,   488,   492,   496,
     500,   504,   508,   512,   516,   520,   524,   528,   532,   536,
     538,   540,   544,   550,   551,   553,   555,   557,   561,   563,
     567,   571,   575,   579,   585,   587,   591,   593,   596,   598,
     602,   606,   609,   611,   613,   615,   617,   619,   624,   629,
     636,   640,   644,   648,   653,   657,   662,   666,   671,   674,
     677,   682,   686
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      77,     0,    -1,    78,    46,    -1,    79,    -1,   106,    -1,
      79,   106,    -1,    -1,    82,    -1,    82,     6,   112,    -1,
      60,    82,    -1,    60,    12,    82,    -1,    62,    82,    -1,
      62,    12,    82,    -1,    85,    82,    -1,    85,    12,    82,
      -1,    35,    80,    38,    -1,    80,     5,    82,    -1,    80,
       5,    35,    80,    38,    -1,    80,     5,    82,     6,   112,
      -1,    80,     5,    60,    82,    -1,    80,     5,    60,    12,
      82,    -1,    80,     5,    62,    82,    -1,    80,     5,    62,
      12,    82,    -1,    80,     5,    85,    82,    -1,    80,     5,
      85,    12,    82,    -1,    82,    -1,    81,     5,    82,    -1,
      42,    -1,    60,    -1,    62,    -1,    61,    -1,    42,    -1,
      42,     6,   112,    -1,    42,    34,    84,    37,    -1,    83,
       5,    83,    -1,   113,    -1,    60,    42,    -1,    61,    42,
      -1,    62,    42,    -1,    42,     6,   113,    -1,    84,     5,
     113,    -1,    84,     5,    82,     6,   113,    -1,    58,    -1,
      58,    35,    58,    38,    -1,    58,    35,    58,     5,    58,
      38,    -1,    58,    16,    58,    17,    -1,    42,    -1,    42,
      35,   113,    38,    -1,    42,     6,   113,    -1,    35,    81,
      38,    -1,    35,    81,    38,    35,   113,    38,    -1,    35,
      81,    38,     6,   113,    -1,    42,    34,   113,    37,    -1,
      35,    81,    38,    34,   113,    37,    -1,    60,    -1,    61,
      -1,    62,    -1,    88,    -1,    88,    16,    58,    17,    -1,
      87,    -1,    90,     5,    87,    -1,    86,    -1,    91,     5,
      86,    -1,    89,    91,    -1,    89,    35,    58,    38,    90,
      -1,    42,    34,    84,    37,    -1,    93,    -1,    94,     5,
      93,    -1,    -1,    85,    96,    83,    71,    -1,    43,    94,
      71,    -1,    92,    71,    -1,    59,    42,     6,   110,    71,
      -1,    -1,    59,    85,    42,    34,    80,    37,    97,    72,
      79,    73,    -1,    -1,    59,    42,    34,    80,    37,    98,
       6,   112,    71,    -1,    72,    -1,    73,    -1,    50,    -1,
      51,    -1,    -1,    85,   104,    83,    -1,    55,    -1,    71,
      -1,    47,    45,    -1,    48,    45,    -1,   105,    72,    79,
      73,   107,    -1,   110,    71,    -1,    95,    -1,   101,    34,
     110,    71,   110,    71,   110,    37,   106,    -1,   101,    34,
     103,    71,   110,    71,   110,    37,   106,    -1,   102,    34,
     110,    37,   106,    -1,     3,    34,   110,    37,   106,    -1,
       3,    34,   110,    37,   106,     4,   106,    -1,    99,    79,
     100,    -1,    68,    42,   109,    -1,    68,    42,    35,   117,
      38,    71,    -1,    52,    71,    -1,    53,    71,    -1,    54,
     110,    71,    -1,    56,    34,    36,    36,    36,    37,   106,
      -1,    34,    42,     6,   110,     5,   110,    37,    -1,   108,
     106,    -1,   112,    -1,   110,     5,   110,    -1,    21,    -1,
      20,    -1,    27,    -1,    29,    -1,    28,    -1,   113,    -1,
     113,     6,   112,    -1,   113,    63,   112,    -1,   113,    64,
     112,    -1,   113,    65,   112,    -1,   113,    66,   112,    -1,
     114,    -1,   114,    74,   113,    75,   113,    -1,   118,    -1,
     114,    22,   114,    -1,   114,    26,   114,    -1,   114,    25,
     114,    -1,   114,    23,   114,    -1,   114,    24,   114,    -1,
     114,    20,   114,    -1,   114,    21,   114,    -1,   114,     9,
     114,    -1,   114,     8,   114,    -1,   114,    12,   114,    -1,
     114,    13,   114,    -1,   114,    10,   114,    -1,   114,    11,
     114,    -1,   114,    16,   114,    -1,   114,    19,   114,    -1,
     114,    17,   114,    -1,   114,    18,   114,    -1,   114,    15,
     114,    -1,   114,    14,   114,    -1,   113,    -1,    75,    -1,
     113,    75,   113,    -1,   113,    75,   113,    75,   113,    -1,
      -1,    60,    -1,    61,    -1,    62,    -1,    82,     6,   113,
      -1,   115,    -1,   116,     5,    60,    -1,   116,     5,    61,
      -1,   116,     5,    62,    -1,   116,     5,   115,    -1,   116,
       5,    82,     6,   113,    -1,   112,    -1,   117,     5,   112,
      -1,   119,    -1,   111,   119,    -1,   120,    -1,   120,    31,
     118,    -1,   120,    33,   118,    -1,   120,    32,    -1,    42,
      -1,    39,    -1,    40,    -1,    41,    -1,    45,    -1,   120,
      34,   116,    37,    -1,   120,    35,   115,    38,    -1,   120,
      35,   115,     5,   115,    38,    -1,   120,    35,    38,    -1,
     120,    36,    42,    -1,    60,    36,    42,    -1,    60,    34,
     116,    37,    -1,    61,    36,    42,    -1,    61,    34,   116,
      37,    -1,    62,    36,    42,    -1,    62,    34,   116,    37,
      -1,   120,    29,    -1,   120,    28,    -1,    58,    34,   110,
      37,    -1,    34,   110,    37,    -1,    35,   117,    38,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short int yyrline[] =
{
       0,   249,   249,   289,   292,   293,   296,   297,   298,   299,
     300,   301,   302,   303,   304,   305,   306,   307,   308,   309,
     310,   311,   312,   313,   314,   317,   318,   321,   321,   321,
     321,   323,   324,   325,   327,   334,   335,   336,   337,   338,
     339,   340,   343,   344,   345,   346,   353,   354,   355,   356,
     357,   358,   361,   362,   366,   366,   366,   367,   368,   373,
     374,   376,   377,   379,   380,   384,   387,   388,   391,   391,
     392,   393,   394,   396,   395,   411,   410,   419,   420,   422,
     424,   429,   429,   432,   434,   435,   436,   437,   438,   439,
     440,   441,   445,   446,   447,   448,   450,   452,   455,   459,
     463,   470,   473,   479,   485,   486,   491,   492,   493,   494,
     495,   499,   500,   501,   502,   503,   504,   509,   510,   513,
     514,   515,   516,   517,   518,   519,   520,   521,   522,   523,
     524,   525,   526,   527,   528,   529,   530,   531,   532,   537,
     538,   539,   540,   543,   544,   545,   546,   547,   548,   549,
     550,   551,   552,   553,   556,   557,   561,   562,   565,   566,
     567,   568,   572,   573,   574,   575,   576,   577,   578,   579,
     580,   581,   582,   583,   584,   585,   586,   587,   588,   589,
     590,   599,   600
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "IF", "ELSE", "','", "'='", "SET",
  "GTGT", "LTLT", "'|'", "OR", "'&'", "AND", "NE", "EQ", "'<'", "'>'",
  "GE", "LE", "'+'", "'-'", "'*'", "'/'", "'%'", "DOTSLASH", "DOTSTAR",
  "'!'", "MOINSMOINS", "PLUSPLUS", "UNARY", "'^'", "'''", "'_'", "'('",
  "'['", "'.'", "')'", "']'", "LNUM", "DNUM", "CNUM", "ID", "FESPACEID",
  "IDPARAM", "STRING", "ENDOFFILE", "INCLUDE", "LOAD", "BIDON", "FOR",
  "WHILE", "BREAK", "CONTINUE", "RETURN", "TRY", "CATCH", "THROW", "TYPE",
  "FUNCTION", "FESPACE", "FESPACE1", "FESPACE3", "PLUSEQ", "MOINSEQ",
  "MULEQ", "DIVEQ", "ARROW", "BORDER", "CURVE", "SOLVE", "';'", "'{'",
  "'}'", "'?'", "':'", "$accept", "start", "input", "instructions",
  "list_of_id_args", "list_of_id1", "id", "list_of_dcls",
  "parameters_list", "type_of_dcl", "ID_space", "ID_array_space",
  "fespace123", "fespace", "spaceIDa", "spaceIDb", "spaceIDs",
  "fespace_def", "fespace_def_list", "declaration", "@1", "@2", "@3",
  "begin", "end", "for_loop", "while_loop", "declaration_for", "@4", "try",
  "instruction", "catchs", "bornes", "border_expr", "Expr", "unop",
  "no_comma_expr", "no_set_expr", "no_ternary_expr", "sub_script_expr",
  "parameters", "array", "unary_expr", "pow_expr", "primary", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short int yytoknum[] =
{
       0,   256,   257,   258,   259,    44,    61,   260,   261,   262,
     124,   263,    38,   264,   265,   266,    60,    62,   267,   268,
      43,    45,    42,    47,    37,   269,   270,    33,   271,   272,
     273,    94,    39,    95,    40,    91,    46,    41,    93,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,    59,   123,   125,    63,    58
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    76,    77,    78,    79,    79,    80,    80,    80,    80,
      80,    80,    80,    80,    80,    80,    80,    80,    80,    80,
      80,    80,    80,    80,    80,    81,    81,    82,    82,    82,
      82,    83,    83,    83,    83,    84,    84,    84,    84,    84,
      84,    84,    85,    85,    85,    85,    86,    86,    86,    86,
      86,    86,    87,    87,    88,    88,    88,    89,    89,    90,
      90,    91,    91,    92,    92,    93,    94,    94,    96,    95,
      95,    95,    95,    97,    95,    98,    95,    99,   100,   101,
     102,   104,   103,   105,   106,   106,   106,   106,   106,   106,
     106,   106,   106,   106,   106,   106,   106,   106,   106,   106,
     106,   107,   108,   109,   110,   110,   111,   111,   111,   111,
     111,   112,   112,   112,   112,   112,   112,   113,   113,   114,
     114,   114,   114,   114,   114,   114,   114,   114,   114,   114,
     114,   114,   114,   114,   114,   114,   114,   114,   114,   115,
     115,   115,   115,   116,   116,   116,   116,   116,   116,   116,
     116,   116,   116,   116,   117,   117,   118,   118,   119,   119,
     119,   119,   120,   120,   120,   120,   120,   120,   120,   120,
     120,   120,   120,   120,   120,   120,   120,   120,   120,   120,
     120,   120,   120
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     2,     1,     1,     2,     0,     1,     3,     2,
       3,     2,     3,     2,     3,     3,     3,     5,     5,     4,
       5,     4,     5,     4,     5,     1,     3,     1,     1,     1,
       1,     1,     3,     4,     3,     1,     2,     2,     2,     3,
       3,     5,     1,     4,     6,     4,     1,     4,     3,     3,
       6,     5,     4,     6,     1,     1,     1,     1,     4,     1,
       3,     1,     3,     2,     5,     4,     1,     3,     0,     4,
       3,     2,     5,     0,    10,     0,     9,     1,     1,     1,
       1,     0,     3,     1,     1,     2,     2,     5,     2,     1,
       9,     9,     5,     5,     7,     3,     3,     6,     2,     2,
       3,     7,     7,     2,     1,     3,     1,     1,     1,     1,
       1,     1,     3,     3,     3,     3,     3,     1,     5,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       1,     3,     5,     0,     1,     1,     1,     3,     1,     3,
       3,     3,     3,     5,     1,     3,     1,     2,     1,     3,
       3,     2,     1,     1,     1,     1,     1,     4,     4,     6,
       3,     3,     3,     4,     3,     4,     3,     4,     2,     2,
       4,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       0,     0,   107,   106,   108,   110,   109,     0,     0,   163,
     164,   165,   162,     0,   166,     0,     0,    79,    80,     0,
       0,     0,    83,    42,     0,    54,    55,    56,     0,    84,
      77,     0,     0,     3,    68,    57,     0,     0,    89,     0,
       0,     0,     0,     4,     0,     0,   104,   111,   117,   119,
     156,   158,     0,     0,     0,     0,     0,     0,   154,     0,
       0,    66,     0,    85,    86,    98,    99,     0,     0,     0,
       0,     0,    42,     0,   143,     0,   143,     0,   143,     0,
       0,     1,     2,     5,     0,     0,     0,    46,    61,    63,
      71,     0,     0,     0,     0,     0,    88,   157,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   179,   178,     0,   161,     0,   143,     0,
       0,     0,   181,     0,   182,     0,     0,    70,   100,     0,
       0,     0,     0,     6,     0,   162,   144,   145,   146,   140,
       0,   139,   148,     0,   172,     0,   174,     0,   176,     0,
       0,     0,    96,    31,     0,     0,    27,     0,    28,    30,
      29,     0,    25,     0,     0,     0,    78,    95,    81,     0,
       0,     0,     0,   105,   112,   113,   114,   115,   116,   128,
     127,   131,   132,   129,   130,   138,   137,   133,   135,   136,
     134,   125,   126,   120,   123,   124,   122,   121,     0,   159,
     160,     0,   170,     0,   171,     0,   155,   162,     0,     0,
       0,     0,    35,    67,    45,   180,     0,    43,     0,     6,
      28,    29,     0,     7,     0,     6,     0,     0,     0,   173,
     175,   177,     0,     0,   103,     0,     0,     0,    69,    58,
       0,     0,    49,    48,     0,     0,    62,     0,     0,     0,
       0,     0,     0,   167,     0,   168,    93,     0,    36,    37,
      38,     0,    65,     0,    72,     0,     0,     9,     0,    11,
       0,    75,     0,     0,    13,     0,   147,   141,   149,   150,
     151,     0,   152,     0,     0,    32,     0,    34,     0,     0,
      59,    64,    26,     0,     0,    47,    82,     0,     0,    92,
       0,    87,   118,     0,     0,    39,    28,    30,    29,     0,
      40,    44,    15,    10,    12,     6,    28,    29,    16,     0,
       0,     8,    14,    73,     0,     0,     0,    97,    33,     0,
       0,     0,    51,     0,     0,     0,     0,   169,    94,     0,
       0,     0,    19,     0,    21,     0,     0,    23,     0,     0,
     142,   153,     0,     0,     0,    60,    50,     0,     0,     0,
      41,    17,    20,    22,    18,    24,     0,     0,   105,     0,
      52,     0,     0,     0,    76,     0,   102,     0,    91,    90,
       0,    74,    53,     0,   101
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short int yydefgoto[] =
{
      -1,    31,    32,    33,   232,   171,   150,   164,   221,    34,
      88,   300,    35,    36,   301,    89,    37,    61,    62,    38,
      84,   359,   330,    39,   177,    40,    41,   179,   257,    42,
      43,   311,   161,   162,    44,    45,    46,    47,    48,   152,
     153,    59,    49,    50,    51
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -220
static const short int yypact[] =
{
     500,   -20,  -220,  -220,  -220,  -220,  -220,   675,   675,  -220,
    -220,  -220,  -220,    28,  -220,   101,   105,  -220,  -220,     5,
     145,   675,  -220,   248,   205,    14,    99,   222,   183,  -220,
    -220,   176,   190,   500,  -220,   230,    26,   179,  -220,   500,
     221,   227,   216,  -220,     3,   251,  -220,    48,     9,  -220,
    -220,   604,   675,   250,    14,    99,   222,   116,  -220,    34,
     255,  -220,     6,  -220,  -220,  -220,  -220,     7,   236,   675,
     240,    58,    40,   266,   103,   273,   103,   274,   103,   275,
      10,  -220,  -220,  -220,   277,   270,   272,    56,  -220,   319,
    -220,   338,   711,   675,   500,   675,  -220,  -220,   675,   675,
     675,   675,   675,   675,   675,   675,   675,   675,   675,   675,
     675,   675,   675,   675,   675,   675,   675,   675,   675,   675,
     675,   675,   675,  -220,  -220,   675,  -220,   675,   103,   546,
     287,   129,  -220,   675,  -220,   747,    28,  -220,  -220,   314,
     134,    42,   675,   260,   301,   330,   181,   226,   234,  -220,
     331,   278,  -220,   146,  -220,   149,  -220,   152,  -220,   305,
     675,   500,  -220,   203,     8,   333,  -220,   317,  -220,  -220,
    -220,    44,  -220,   675,   675,    85,  -220,  -220,  -220,   286,
      32,   167,   392,  -220,  -220,  -220,  -220,  -220,  -220,   855,
     855,   870,   870,   883,   883,   576,   576,   658,   658,   658,
     658,   320,   320,  -220,  -220,  -220,  -220,  -220,   288,  -220,
    -220,   185,  -220,    47,  -220,   500,  -220,   355,    31,   215,
     263,   186,  -220,  -220,  -220,  -220,   304,  -220,    36,   260,
      39,    98,   192,   358,   113,   260,   675,   675,   589,  -220,
    -220,  -220,   362,    50,  -220,   675,   747,   277,  -220,  -220,
     161,   212,   112,  -220,   332,   212,  -220,   277,   675,   675,
     500,   313,   675,  -220,   632,  -220,   367,   675,  -220,  -220,
    -220,   783,  -220,   336,  -220,    64,   212,  -220,   212,  -220,
     265,  -220,   675,   212,  -220,   194,  -220,   300,   181,   226,
     234,   370,  -220,   675,   311,  -220,   201,  -220,   212,   350,
    -220,   382,  -220,   675,   675,  -220,   389,    37,    38,  -220,
     368,  -220,  -220,   363,   500,  -220,    14,    99,   222,   397,
    -220,  -220,  -220,  -220,  -220,   260,   150,   158,   398,   182,
     399,  -220,  -220,  -220,   675,   675,   402,  -220,  -220,    79,
     675,   161,  -220,   376,   675,   675,   372,  -220,  -220,   675,
     114,   212,  -220,   212,  -220,   675,   212,  -220,   675,   343,
    -220,  -220,   675,   383,   379,  -220,  -220,   202,   208,   386,
    -220,  -220,  -220,  -220,  -220,  -220,   347,   500,   387,   675,
    -220,   500,   500,   393,  -220,   446,  -220,   388,  -220,  -220,
     391,  -220,  -220,   500,  -220
};

/* YYPGOTO[NTERM-NUM].  */
static const short int yypgoto[] =
{
    -220,  -220,  -220,   -37,  -219,   125,   -50,  -131,   184,   -21,
     261,    97,  -220,  -220,  -220,  -220,  -220,   312,  -220,  -220,
    -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,  -220,
     -33,  -220,  -220,  -220,    -6,  -220,    -4,   -69,   743,  -123,
     -38,   281,   142,   410,  -220
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -31
static const short int yytable[] =
{
      83,    57,    91,    73,    58,   151,   213,   151,    95,   151,
     275,   136,    95,   247,    52,    67,   285,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,   119,   120,   121,   172,    95,   155,   133,
     157,    95,    95,    95,   159,   160,   131,   226,    74,   251,
      75,   276,   264,   208,    98,   133,    68,   182,    83,   151,
     151,    86,   173,   140,   142,    74,   222,    75,    87,   280,
      60,   178,   134,   268,    96,    70,    65,   137,   138,   248,
     227,   166,   252,   122,   251,   265,   180,   181,   294,   183,
     211,   174,   143,   233,   184,   185,   186,   187,   188,   168,
     169,   170,   322,   259,   253,   254,   350,   274,   344,   345,
     278,    99,   100,   101,   102,   292,   297,   363,   303,   280,
     255,    95,   234,     2,     3,   283,   306,    87,   244,   216,
       4,     5,     6,    76,    95,    77,   228,     7,     8,    95,
     166,   313,     9,    10,    11,   145,    63,   304,    14,    83,
      64,   238,   371,   132,   238,   166,    58,   238,   168,   169,
     170,    53,   351,   146,   147,   148,   215,   286,   287,   151,
     353,   225,    95,   168,   169,   170,    81,   222,   149,   233,
     277,   279,   266,   239,   284,   233,   240,   -28,   291,   241,
     238,   271,   166,   312,   356,   151,   298,   280,   315,   280,
     166,   302,   320,   299,   260,   172,   271,    95,   234,   245,
     168,   169,   170,    95,   234,    74,    66,    75,   168,   169,
     170,   319,   263,   272,   166,    80,   323,   309,   324,   281,
     328,   333,   -30,   332,   342,   343,    82,   246,   338,   381,
     -29,   295,   168,   169,   170,   382,    85,    71,   172,    76,
      90,    77,   307,   308,   166,    92,    78,   269,    79,   329,
      76,    93,    77,    72,    68,   360,   361,   209,    78,   210,
      79,   364,   168,   169,   170,   233,   352,   354,   331,   357,
     370,   348,    69,    70,    69,     7,     8,   336,    94,   135,
       9,    10,    11,    12,   139,   229,    14,    78,   141,    79,
     325,   372,   166,   373,   234,   270,   375,   166,   144,    53,
     387,    54,    55,    56,   166,   154,   156,   158,    72,   163,
     230,   169,   231,    72,   175,   326,   169,   327,   165,   214,
     167,   224,   168,   169,   170,   235,   -27,   236,   367,   368,
     385,     1,   117,   118,   119,   120,   121,   242,   388,   389,
     249,   374,    83,   237,   376,   250,   378,   258,     2,     3,
     394,   267,   273,   262,   282,     4,     5,     6,   293,   310,
     305,   314,     7,     8,   321,   334,   335,     9,    10,    11,
      12,    13,   337,    14,   340,    15,    16,   341,    17,    18,
      19,    20,    21,    22,   247,     1,    23,    24,    25,    26,
      27,   347,   346,   349,   355,   358,    28,   362,   369,    29,
      30,   176,     2,     3,   366,   377,   380,   379,   384,     4,
       5,     6,   383,   339,   386,   392,     7,     8,   393,   390,
     296,     9,    10,    11,    12,    13,   256,    14,   365,    15,
      16,   243,    17,    18,    19,    20,    21,    22,   223,     1,
      23,    24,    25,    26,    27,    97,     0,     0,     0,     0,
      28,     0,     0,    29,    30,   261,     2,     3,     0,     0,
       0,     0,     0,     4,     5,     6,     0,     0,     0,     0,
       7,     8,     0,     0,     0,     9,    10,    11,    12,    13,
       0,    14,     0,    15,    16,     0,    17,    18,    19,    20,
      21,    22,     0,     1,    23,    24,    25,    26,    27,     0,
       0,     0,     0,     0,    28,     0,     0,    29,    30,   391,
       2,     3,     0,     0,     0,     0,     0,     4,     5,     6,
       0,     0,     0,     0,     7,     8,     0,     0,     0,     9,
      10,    11,    12,    13,     0,    14,     0,    15,    16,     0,
      17,    18,    19,    20,    21,    22,     0,     0,    23,    24,
      25,    26,    27,     0,     0,     0,     2,     3,    28,     0,
       0,    29,    30,     4,     5,     6,     0,     0,     0,     0,
       7,     8,     0,     0,   212,     9,    10,    11,    12,     0,
       0,    14,   111,   112,   113,   114,   115,   116,   117,   118,
     119,   120,   121,     0,    53,     0,    54,    55,    56,     2,
       3,     0,     0,     0,     0,     0,     4,     5,     6,     0,
       0,   149,     0,     7,     8,     0,     0,     0,     9,    10,
      11,   145,   123,   124,    14,   125,   126,   127,   128,   129,
     130,     0,     0,     0,     0,     0,     0,    53,     0,   288,
     289,   290,     2,     3,     0,     0,     0,     0,     0,     4,
       5,     6,     0,     0,   149,     0,     7,     8,     0,     0,
       0,     9,    10,    11,    12,     0,     0,    14,   115,   116,
     117,   118,   119,   120,   121,     0,     0,     0,     0,     0,
      53,     0,    54,    55,    56,     2,     3,     0,     0,     0,
       0,     0,     4,     5,     6,     0,     0,   149,     0,     7,
       8,     0,     0,     0,     9,    10,    11,    12,     0,     0,
      14,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     2,     3,    53,     0,    54,    55,    56,     4,     5,
       6,     0,     0,     0,     0,     7,     8,     0,     0,     0,
       9,    10,    11,    12,     0,     0,    14,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     2,     3,    23,
       0,    54,    55,    56,     4,     5,     6,     0,     0,     0,
       0,     7,     8,     0,     0,     0,     9,    10,    11,   217,
       0,     0,    14,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     2,     3,    53,     0,   218,   219,   220,
       4,     5,     6,     0,     0,     0,     0,     7,     8,     0,
       0,     0,     9,    10,    11,   145,     0,     0,    14,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    53,     0,   316,   317,   318,   189,   190,   191,   192,
     193,   194,   195,   196,   197,   198,   199,   200,   201,   202,
     203,   204,   205,   206,   207,   105,   106,   107,   108,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     120,   121,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   109,   110,   111,
     112,   113,   114,   115,   116,   117,   118,   119,   120,   121
};

static const short int yycheck[] =
{
      33,     7,    39,    24,     8,    74,   129,    76,     5,    78,
     229,     5,     5,     5,    34,    21,   235,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    86,     5,    76,     5,
      78,     5,     5,     5,    34,    35,    52,     5,    34,     5,
      36,    12,     5,   122,     6,     5,    16,    94,    91,   128,
     129,    35,     6,    69,     6,    34,   135,    36,    42,     5,
      42,    92,    38,    42,    71,    35,    71,    71,    71,    71,
      38,    42,    38,    74,     5,    38,    92,    93,    38,    95,
     128,    35,    34,   143,    98,    99,   100,   101,   102,    60,
      61,    62,    38,    71,   173,   174,   325,    71,    71,    71,
      12,    63,    64,    65,    66,   238,   247,    38,     6,     5,
      35,     5,   143,    20,    21,    12,   257,    42,   161,   133,
      27,    28,    29,    34,     5,    36,   142,    34,    35,     5,
      42,   264,    39,    40,    41,    42,    45,    35,    45,   182,
      45,     5,    38,    37,     5,    42,   160,     5,    60,    61,
      62,    58,    12,    60,    61,    62,    37,   236,   237,   238,
      12,    37,     5,    60,    61,    62,     0,   246,    75,   229,
     230,   231,   215,    37,   234,   235,    37,     6,   238,    37,
       5,     5,    42,   262,    12,   264,    35,     5,   267,     5,
      42,   251,   271,    42,    37,   255,     5,     5,   229,     6,
      60,    61,    62,     5,   235,    34,    71,    36,    60,    61,
      62,   271,    37,    37,    42,    42,   276,   260,   278,    37,
     280,    37,     6,   283,   303,   304,    46,    34,    37,    37,
       6,   245,    60,    61,    62,    37,    16,    42,   298,    34,
      71,    36,   258,   259,    42,    34,    34,    42,    36,   280,
      34,    34,    36,    58,    16,   334,   335,   125,    34,   127,
      36,   340,    60,    61,    62,   325,   326,   327,   282,   329,
     349,   314,    34,    35,    34,    34,    35,   293,    72,    34,
      39,    40,    41,    42,    58,    35,    45,    34,    58,    36,
      35,   351,    42,   353,   325,    42,   356,    42,    42,    58,
     379,    60,    61,    62,    42,    42,    42,    42,    58,    42,
      60,    61,    62,    58,     5,    60,    61,    62,    58,    42,
      58,    17,    60,    61,    62,    34,     6,     6,   344,   345,
     377,     3,    22,    23,    24,    25,    26,    42,   381,   382,
      17,   355,   385,    75,   358,    38,   362,    71,    20,    21,
     393,     6,    58,    75,     6,    27,    28,    29,     6,    56,
      38,     4,    34,    35,    38,    75,     6,    39,    40,    41,
      42,    43,    71,    45,    34,    47,    48,     5,    50,    51,
      52,    53,    54,    55,     5,     3,    58,    59,    60,    61,
      62,    38,    34,     6,     6,     6,    68,     5,    36,    71,
      72,    73,    20,    21,    38,    72,    37,    34,    71,    27,
      28,    29,    36,   298,    37,    37,    34,    35,    37,    36,
     246,    39,    40,    41,    42,    43,   175,    45,   341,    47,
      48,   160,    50,    51,    52,    53,    54,    55,   136,     3,
      58,    59,    60,    61,    62,    45,    -1,    -1,    -1,    -1,
      68,    -1,    -1,    71,    72,    73,    20,    21,    -1,    -1,
      -1,    -1,    -1,    27,    28,    29,    -1,    -1,    -1,    -1,
      34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    43,
      -1,    45,    -1,    47,    48,    -1,    50,    51,    52,    53,
      54,    55,    -1,     3,    58,    59,    60,    61,    62,    -1,
      -1,    -1,    -1,    -1,    68,    -1,    -1,    71,    72,    73,
      20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,
      -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,
      40,    41,    42,    43,    -1,    45,    -1,    47,    48,    -1,
      50,    51,    52,    53,    54,    55,    -1,    -1,    58,    59,
      60,    61,    62,    -1,    -1,    -1,    20,    21,    68,    -1,
      -1,    71,    72,    27,    28,    29,    -1,    -1,    -1,    -1,
      34,    35,    -1,    -1,    38,    39,    40,    41,    42,    -1,
      -1,    45,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    -1,    58,    -1,    60,    61,    62,    20,
      21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,    -1,
      -1,    75,    -1,    34,    35,    -1,    -1,    -1,    39,    40,
      41,    42,    28,    29,    45,    31,    32,    33,    34,    35,
      36,    -1,    -1,    -1,    -1,    -1,    -1,    58,    -1,    60,
      61,    62,    20,    21,    -1,    -1,    -1,    -1,    -1,    27,
      28,    29,    -1,    -1,    75,    -1,    34,    35,    -1,    -1,
      -1,    39,    40,    41,    42,    -1,    -1,    45,    20,    21,
      22,    23,    24,    25,    26,    -1,    -1,    -1,    -1,    -1,
      58,    -1,    60,    61,    62,    20,    21,    -1,    -1,    -1,
      -1,    -1,    27,    28,    29,    -1,    -1,    75,    -1,    34,
      35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    20,    21,    58,    -1,    60,    61,    62,    27,    28,
      29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    20,    21,    58,
      -1,    60,    61,    62,    27,    28,    29,    -1,    -1,    -1,
      -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    20,    21,    58,    -1,    60,    61,    62,
      27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    58,    -1,    60,    61,    62,   103,   104,   105,   106,
     107,   108,   109,   110,   111,   112,   113,   114,   115,   116,
     117,   118,   119,   120,   121,    10,    11,    12,    13,    14,
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
      53,    54,    55,    58,    59,    60,    61,    62,    68,    71,
      72,    77,    78,    79,    85,    88,    89,    92,    95,    99,
     101,   102,   105,   106,   110,   111,   112,   113,   114,   118,
     119,   120,    34,    58,    60,    61,    62,   110,   112,   117,
      42,    93,    94,    45,    45,    71,    71,   110,    16,    34,
      35,    42,    58,    85,    34,    36,    34,    36,    34,    36,
      42,     0,    46,   106,    96,    16,    35,    42,    86,    91,
      71,    79,    34,    34,    72,     5,    71,   119,     6,    63,
      64,    65,    66,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    74,    28,    29,    31,    32,    33,    34,    35,
      36,   110,    37,     5,    38,    34,     5,    71,    71,    58,
     110,    58,     6,    34,    42,    42,    60,    61,    62,    75,
      82,   113,   115,   116,    42,   116,    42,   116,    42,    34,
      35,   108,   109,    42,    83,    58,    42,    58,    60,    61,
      62,    81,    82,     6,    35,     5,    73,   100,    85,   103,
     110,   110,    79,   110,   112,   112,   112,   112,   112,   114,
     114,   114,   114,   114,   114,   114,   114,   114,   114,   114,
     114,   114,   114,   114,   114,   114,   114,   114,   113,   118,
     118,   116,    38,   115,    42,    37,   112,    42,    60,    61,
      62,    84,   113,    93,    17,    37,     5,    38,   110,    35,
      60,    62,    80,    82,    85,    34,     6,    75,     5,    37,
      37,    37,    42,   117,   106,     6,    34,     5,    71,    17,
      38,     5,    38,   113,   113,    35,    86,   104,    71,    71,
      37,    73,    75,    37,     5,    38,   106,     6,    42,    42,
      42,     5,    37,    58,    71,    80,    12,    82,    12,    82,
       5,    37,     6,    12,    82,    80,   113,   113,    60,    61,
      62,    82,   115,     6,    38,   112,    84,    83,    35,    42,
      87,    90,    82,     6,    35,    38,    83,   110,   110,   106,
      56,   107,   113,   115,     4,   113,    60,    61,    62,    82,
     113,    38,    38,    82,    82,    35,    60,    62,    82,    85,
      98,   112,    82,    37,    75,     6,   110,    71,    37,    81,
      34,     5,   113,   113,    71,    71,    34,    38,   106,     6,
      80,    12,    82,    12,    82,     6,    12,    82,     6,    97,
     113,   113,     5,    38,   113,    87,    38,   110,   110,    36,
     113,    38,    82,    82,   112,    82,   112,    72,   110,    34,
      37,    37,    37,    36,    71,    79,    37,   113,   106,   106,
      36,    73,    37,    37,   106
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
#define YYERROR		goto yyerrorlab


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


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (N)								\
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (0)
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
              (Loc).first_line, (Loc).first_column,	\
              (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
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

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Type, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short int *bottom, short int *top)
#else
static void
yy_stack_print (bottom, top)
    short int *bottom;
    short int *top;
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
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
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
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);


# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

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



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
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
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short int yyssa[YYINITDEPTH];
  short int *yyss = yyssa;
  register short int *yyssp;

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


  yyvsp[0] = yylval;

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
	short int *yyss1 = yyss;


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
	short int *yyss1 = yyss;
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
/* Read a look-ahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to look-ahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
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
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
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

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

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
#line 249 "lg.y"
    {
		        const char *  magicffglut="#!ffglutdata...\n";
                        if(ThePlotStream) fwrite(magicffglut,strlen(magicffglut),1,ThePlotStream);	            
                        size_t sizestack = currentblock->size()+1024 ; //  before close 
                        (yyvsp[-1].cinst)+=currentblock->close(currentblock);
                        if(verbosity) cout << " sizestack + 1024 =" << sizestack << "  ( " << sizestack-1024 <<" )\n" ;   
                        size_t lg0,lg1;                       
                        int NbPtr = ShowAlloc("init execution ",lg0); // number of un delele ptr
                        if(verbosity) cout << endl;  
                        { Stack stack = newStack(sizestack);
                        double CPUcompile= CPUtime();
                        try {                  
                          (yyvsp[-1].cinst).eval(stack);}
                        catch ( E_exception & e)  {
                          cerr << e.what() << " ,  mpirank " << mpirank << endl;
                          return 1; }
                        catch( Error & err) {
                          cerr << err.what() << endl;
			  cerr << " err code " << err.errcode() << " ,  mpirank " << mpirank << endl;
                          return err.errcode();
                        }
                         catch( ...) { cerr << "Strange catch exception ???\n"; 
                          cerr << " at exec line  " << TheCurrentLine << " ,  mpirank " << mpirank << endl;
                          return 1; 
                         }

                        if(verbosity)  cout << "times: compile "<< CPUcompile-CPUcompileInit <<"s, execution " 
			    <<  CPUtime()-CPUcompile  <<"s,  mpirank:" << mpirank << endl;
                        deleteStack(stack);
                        //debugstack.clear() 
                        } 
                        fingraphique();
			if(ThePlotStream) {pclose(ThePlotStream); ThePlotStream=0;}
                        NbPtr = ShowAlloc("end execution -- ",lg1) - NbPtr;
                        
			    if (NbPtr) { cout << " ######## We forget of deleting   " << NbPtr 
			                      << " Nb pointer,   " <<  lg1-lg0 << "Bytes " << " ,  mpirank " << mpirank <<endl;}
  return 0;;}
    break;

  case 4:
#line 292 "lg.y"
    {(yyval.cinst)=(yyvsp[0].cexp);;;;}
    break;

  case 5:
#line 293 "lg.y"
    { (yyval.cinst)= ((yyvsp[-1].cinst)+=(yyvsp[0].cexp)) ;}
    break;

  case 6:
#line 296 "lg.y"
    { (yyval.clist_id)=new ListOfId();;}
    break;

  case 7:
#line 297 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
    break;

  case 8:
#line 298 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str),(yyvsp[0].cexp))) ;}
    break;

  case 9:
#line 299 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,2> **>()));}
    break;

  case 10:
#line 300 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,2> **>(),true));}
    break;

  case 11:
#line 301 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,3> **>()));}
    break;

  case 12:
#line 302 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,3> **>(),true));}
    break;

  case 13:
#line 303 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-1].type)->right())) ;}
    break;

  case 14:
#line 304 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-2].type),true)) ;}
    break;

  case 15:
#line 305 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id))) ;}
    break;

  case 16:
#line 306 "lg.y"
    { (yyval.clist_id) = (yyvsp[-2].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str))) ;}
    break;

  case 17:
#line 307 "lg.y"
    { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id))) ;}
    break;

  case 18:
#line 308 "lg.y"
    { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str),(yyvsp[0].cexp))) ;}
    break;

  case 19:
#line 309 "lg.y"
    { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,2> **>())) ;}
    break;

  case 20:
#line 310 "lg.y"
    { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,2> **>(),true)) ;}
    break;

  case 21:
#line 311 "lg.y"
    { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,3> **>())) ;}
    break;

  case 22:
#line 312 "lg.y"
    { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,3> **>(),true)) ;}
    break;

  case 23:
#line 313 "lg.y"
    { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-1].type)->right())) ;}
    break;

  case 24:
#line 314 "lg.y"
    { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-2].type),true)) ;}
    break;

  case 25:
#line 317 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str))); ;}
    break;

  case 26:
#line 318 "lg.y"
    { (yyval.clist_id)=(yyvsp[-2].clist_id)  ; (yyval.clist_id)->push_back(UnId((yyvsp[0].str))); ;}
    break;

  case 31:
#line 323 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[0].str),dcltype);}
    break;

  case 32:
#line 324 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[-2].str),dcltype,(yyvsp[0].cexp));}
    break;

  case 33:
#line 325 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[-3].str),dcltype,(yyvsp[-1].args));
                                              (yyvsp[-1].args).destroy();}
    break;

  case 34:
#line 327 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 35:
#line 334 "lg.y"
    {(yyval.args)=(yyvsp[0].cexp);}
    break;

  case 36:
#line 335 "lg.y"
    {(yyval.args)=Find((yyvsp[-1].str));}
    break;

  case 37:
#line 336 "lg.y"
    {(yyval.args)=Find((yyvsp[-1].str));}
    break;

  case 38:
#line 337 "lg.y"
    {(yyval.args)=Find((yyvsp[-1].str));}
    break;

  case 39:
#line 338 "lg.y"
    { (yyval.args)=make_pair<const char *,const C_F0>((yyvsp[-2].str),(yyvsp[0].cexp));}
    break;

  case 40:
#line 339 "lg.y"
    { (yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp)) ;}
    break;

  case 41:
#line 340 "lg.y"
    { (yyval.args)= ((yyvsp[-4].args)+= make_pair<const char *,const C_F0>((yyvsp[-2].str),(yyvsp[0].cexp)));}
    break;

  case 43:
#line 344 "lg.y"
    {(yyval.type)=TypeArray((yyvsp[-3].type),(yyvsp[-1].type));}
    break;

  case 44:
#line 345 "lg.y"
    {(yyval.type)=TypeArray((yyvsp[-5].type),(yyvsp[-3].type),(yyvsp[-1].type));}
    break;

  case 45:
#line 346 "lg.y"
    {(yyval.type)=TypeTemplate((yyvsp[-3].type),(yyvsp[-1].type));}
    break;

  case 46:
#line 353 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[0].str),currentblock,fespacetype,fespacecomplex,fespacedim); ;}
    break;

  case 47:
#line 354 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[-3].str),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 48:
#line 355 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[-2].str),currentblock,fespacetype,(yyvsp[0].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 49:
#line 356 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[-1].clist_id),currentblock,fespacetype,fespacecomplex,fespacedim) ;}
    break;

  case 50:
#line 357 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[-4].clist_id),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 51:
#line 358 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[-3].clist_id),currentblock,fespacetype,(yyvsp[0].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 52:
#line 361 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[-3].str),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 53:
#line 362 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[-4].clist_id),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 54:
#line 366 "lg.y"
    { fespacedim=2;}
    break;

  case 55:
#line 366 "lg.y"
    { fespacedim=1;}
    break;

  case 56:
#line 366 "lg.y"
    { fespacedim=3;}
    break;

  case 57:
#line 367 "lg.y"
    {fespacecomplex=false;  fespacetype = Find((yyvsp[0].str));;}
    break;

  case 58:
#line 368 "lg.y"
    {
             if ((yyvsp[-1].type) != typevarreal && (yyvsp[-1].type) != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=((yyvsp[-1].type)==typevarcomplex);
             fespacetype = Find((yyvsp[-3].str));;}
    break;

  case 59:
#line 373 "lg.y"
    {  (yyval.cexp) = (yyvsp[0].cexp)  ;}
    break;

  case 60:
#line 374 "lg.y"
    { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));;}
    break;

  case 61:
#line 376 "lg.y"
    {  (yyval.cexp) = (yyvsp[0].cexp)  ;}
    break;

  case 62:
#line 377 "lg.y"
    { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));;}
    break;

  case 63:
#line 379 "lg.y"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[0].cexp);}
    break;

  case 64:
#line 380 "lg.y"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[0].cexp);}
    break;

  case 65:
#line 384 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,size_t>((yyvsp[-3].str),typeFESpace((yyvsp[-1].args)),(yyvsp[-1].args),dimFESpaceImage((yyvsp[-1].args)));
     (yyvsp[-1].args).destroy(); ;}
    break;

  case 67:
#line 388 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 68:
#line 391 "lg.y"
    {dcltype=(yyvsp[0].type);}
    break;

  case 69:
#line 391 "lg.y"
    {(yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 70:
#line 392 "lg.y"
    {(yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 71:
#line 393 "lg.y"
    { (yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 72:
#line 394 "lg.y"
    {(yyval.cexp)=currentblock->NewID((yyvsp[-4].type),(yyvsp[-3].str),(yyvsp[-1].cexp));;}
    break;

  case 73:
#line 396 "lg.y"
    {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = (yyvsp[-4].type)->right();
                      routineinblock[kkembtype] = currentblock;
                      (yyvsp[-1].routine)=new Routine((yyvsp[-5].type),(yyvsp[-4].type)->right(),(yyvsp[-3].str),(yyvsp[-1].clist_id),currentblock);
                     // cout << " \n after new routine \n " << endl;                      
                      ;}
    break;

  case 74:
#line 404 "lg.y"
    { currentblock=(yyvsp[-5].routine)->Set((yyvsp[-1].cinst));
                       currentblock->Add((yyvsp[-7].str),"(",(yyvsp[-5].routine));
                       kkembtype--;
                       (yyval.cexp)=0;
                    
                        ;}
    break;

  case 75:
#line 411 "lg.y"
    {Block::open(currentblock); (yyvsp[-4].type)->SetArgs((yyvsp[-1].clist_id));;}
    break;

  case 76:
#line 413 "lg.y"
    {  (yyval.cinst)=currentblock->close(currentblock);
                         (yyval.cexp)=currentblock->NewID((yyvsp[-8].type),(yyvsp[-7].str),(yyvsp[-1].cexp),*(yyvsp[-5].clist_id));
                         delete (yyvsp[-5].clist_id); //  FH 23032005
                         ;}
    break;

  case 77:
#line 419 "lg.y"
    {  Block::open(currentblock);}
    break;

  case 78:
#line 420 "lg.y"
    {  (yyval.cexp)=currentblock->close(currentblock);}
    break;

  case 79:
#line 422 "lg.y"
    {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 80:
#line 424 "lg.y"
    {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 81:
#line 429 "lg.y"
    {dcltype=(yyvsp[0].type); Block::open(currentblock);  ;}
    break;

  case 82:
#line 430 "lg.y"
    {(yyval.cexp)=(yyvsp[0].cexp);}
    break;

  case 83:
#line 432 "lg.y"
    { Block::open(currentblock) ;}
    break;

  case 84:
#line 434 "lg.y"
    {(yyval.cexp)=0;;}
    break;

  case 85:
#line 435 "lg.y"
    {zzzfff->input((yyvsp[0].str));(yyval.cexp)= 0; ;}
    break;

  case 86:
#line 436 "lg.y"
    {load((yyvsp[0].str));(yyval.cexp)= 0; ;}
    break;

  case 87:
#line 437 "lg.y"
    {(yyval.cexp)=Try((yyvsp[-2].cinst),(yyvsp[0].cexp),currentblock->close(currentblock));;}
    break;

  case 88:
#line 438 "lg.y"
    {(yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 89:
#line 439 "lg.y"
    {(yyval.cexp)=(yyvsp[0].cexp);}
    break;

  case 90:
#line 440 "lg.y"
    {inloopcount--; (yyval.cexp)=For((yyvsp[-6].cexp),(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 91:
#line 442 "lg.y"
    {inloopcount--; 
                (yyval.cexp)=C_F0(For((yyvsp[-6].cexp),(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp)),currentblock->close(currentblock));}
    break;

  case 92:
#line 445 "lg.y"
    {inloopcount--;(yyval.cexp)=While((yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 93:
#line 446 "lg.y"
    {(yyval.cexp)=FIf((yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 94:
#line 447 "lg.y"
    {(yyval.cexp)=FIf((yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 95:
#line 448 "lg.y"
    { 
                      (yyval.cexp)=C_F0(new E_block((yyvsp[-1].cinst),(yyvsp[0].cexp)),atype<void>()) ;}
    break;

  case 96:
#line 450 "lg.y"
    {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[-1].str),C_F0(TheOperators,"[border]",(yyvsp[0].args)));}
    break;

  case 97:
#line 452 "lg.y"
    {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[-4].str),C_F0(TheOperators,"[border]",(yyvsp[-2].args)));}
    break;

  case 98:
#line 455 "lg.y"
    {
                    if(inloopcount) 
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;}
    break;

  case 99:
#line 459 "lg.y"
    { 
                    if(inloopcount)
                        (yyval.cexp)= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");}
    break;

  case 100:
#line 463 "lg.y"
    { 
                    if (kkembtype>=0)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo((yyvsp[-1].cexp))) ,atype<void>());
                     else lgerror(" return not in routine ") ;}
    break;

  case 101:
#line 470 "lg.y"
    {(yyval.cexp) =  (yyvsp[0].cexp); ;}
    break;

  case 102:
#line 473 "lg.y"
    { 
   Block::open(currentblock);
   (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[-5].str),atype<double*>());
   (yyval.args)+= (yyvsp[-3].cexp);
   (yyval.args)+= (yyvsp[-1].cexp) ;}
    break;

  case 103:
#line 479 "lg.y"
    {   
   (yyval.args) = ((yyvsp[-1].args) += (yyvsp[0].cexp));
   currentblock->close(currentblock);}
    break;

  case 105:
#line 486 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));;}
    break;

  case 112:
#line 500 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 113:
#line 501 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"+=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 114:
#line 502 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"-=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 115:
#line 503 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"*=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 116:
#line 504 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"/=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 118:
#line 510 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"?:",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 120:
#line 514 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 121:
#line 515 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 122:
#line 516 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 123:
#line 517 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 124:
#line 518 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 125:
#line 519 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 126:
#line 520 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 127:
#line 521 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 128:
#line 522 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 129:
#line 523 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 130:
#line 524 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 131:
#line 525 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 132:
#line 526 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 133:
#line 527 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 134:
#line 528 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 135:
#line 529 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 136:
#line 530 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 137:
#line 531 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 138:
#line 532 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 139:
#line 537 "lg.y"
    {(yyval.cexp)=(yyvsp[0].cexp);}
    break;

  case 140:
#line 538 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,":");}
    break;

  case 141:
#line 539 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 142:
#line 540 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 143:
#line 543 "lg.y"
    {(yyval.args)=0;}
    break;

  case 144:
#line 544 "lg.y"
    {(yyval.args)=Find((yyvsp[0].str));}
    break;

  case 145:
#line 545 "lg.y"
    {(yyval.args)=Find((yyvsp[0].str));}
    break;

  case 146:
#line 546 "lg.y"
    {(yyval.args)=Find((yyvsp[0].str));}
    break;

  case 147:
#line 547 "lg.y"
    { (yyval.args)=make_pair<const char *,const C_F0>((yyvsp[-2].str),(yyvsp[0].cexp));}
    break;

  case 148:
#line 548 "lg.y"
    {(yyval.args)=(yyvsp[0].cexp);}
    break;

  case 149:
#line 549 "lg.y"
    { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str))) ;}
    break;

  case 150:
#line 550 "lg.y"
    { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str))) ;}
    break;

  case 151:
#line 551 "lg.y"
    { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str))) ;}
    break;

  case 152:
#line 552 "lg.y"
    { (yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp)) ;}
    break;

  case 153:
#line 553 "lg.y"
    { (yyval.args)= ((yyvsp[-4].args)+= make_pair<const char *,const C_F0>((yyvsp[-2].str),(yyvsp[0].cexp))) ;}
    break;

  case 154:
#line 556 "lg.y"
    {(yyval.args)=(yyvsp[0].cexp);}
    break;

  case 155:
#line 557 "lg.y"
    {(yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp)) ;}
    break;

  case 157:
#line 562 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[0].cexp));}
    break;

  case 159:
#line 566 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 160:
#line 567 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 161:
#line 568 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
    break;

  case 162:
#line 572 "lg.y"
    {(yyval.cexp)=Find((yyvsp[0].str));;}
    break;

  case 163:
#line 573 "lg.y"
    {(yyval.cexp)= CConstant((yyvsp[0].lnum));}
    break;

  case 164:
#line 574 "lg.y"
    {(yyval.cexp)= CConstant((yyvsp[0].dnum));}
    break;

  case 165:
#line 575 "lg.y"
    {(yyval.cexp)= CConstant(complex<double>(0,(yyvsp[0].dnum)));}
    break;

  case 166:
#line 576 "lg.y"
    {(yyval.cexp)= CConstant<const char *>((yyvsp[0].str));}
    break;

  case 167:
#line 577 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-3].cexp),(yyvsp[-2].oper),(yyvsp[-1].args));;}
    break;

  case 168:
#line 578 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-3].cexp),(yyvsp[-2].oper),(yyvsp[-1].cexp));}
    break;

  case 169:
#line 579 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-5].cexp),(yyvsp[-4].oper),(yyvsp[-3].cexp),(yyvsp[-1].cexp));}
    break;

  case 170:
#line 580 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),"[]");}
    break;

  case 171:
#line 581 "lg.y"
    { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].str)) ;;}
    break;

  case 172:
#line 582 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;;}
    break;

  case 173:
#line 583 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;;}
    break;

  case 174:
#line 584 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;;}
    break;

  case 175:
#line 585 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;;}
    break;

  case 176:
#line 586 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;;}
    break;

  case 177:
#line 587 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;;}
    break;

  case 178:
#line 588 "lg.y"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
    break;

  case 179:
#line 589 "lg.y"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
    break;

  case 180:
#line 590 "lg.y"
    {
             if ((yyvsp[-3].type)->right()->CastingFrom((yyvsp[-1].cexp).left()) ) 
                (yyval.cexp)=(yyvsp[-3].type)->right()->CastTo((yyvsp[-1].cexp))  ;
             else { (yyval.cexp)=(yyvsp[-3].type)->right()->Find("<--",basicAC_F0_wa((yyvsp[-1].cexp)));
             if (!(yyval.cexp).left()) { cerr << " no wait to change " << (yyvsp[-1].cexp).left()->right()->name() << " in " << 
                                        (yyvsp[-3].type)->right()->name() << endl;
                                CompileError(" Error in type(exp) "); }
             }
            ;}
    break;

  case 181:
#line 599 "lg.y"
    {(yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 182:
#line 600 "lg.y"
    { (yyval.cexp)=C_F0(TheOperators,"[]",(yyvsp[-1].args));}
    break;


    }

/* Line 1037 of yacc.c.  */
#line 2595 "lg.tab.cpp"

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
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* If at end of input, pop the error token,
	     then the rest of the stack, then return failure.  */
	  if (yychar == YYEOF)
	     for (;;)
	       {

		 YYPOPSTACK;
		 if (yyssp == yyss)
		   YYABORT;
		 yydestruct ("Error: popping",
                             yystos[*yyssp], yyvsp);
	       }
        }
      else
	{
	  yydestruct ("Error: discarding", yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

#ifdef __GNUC__
  /* Pacify GCC when the user code never invokes YYERROR and the label
     yyerrorlab therefore never appears in user code.  */
  if (0)
     goto yyerrorlab;
#endif

yyvsp -= yylen;
  yyssp -= yylen;
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
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


      yydestruct ("Error: popping", yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token. */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

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
  yydestruct ("Error: discarding lookahead",
              yytoken, &yylval);
  yychar = YYEMPTY;
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


#line 605 "lg.y"
 


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
	{retvalue=1; if(!mpirank) cerr <<  "Error:a block is not close" << endl; }  
      else {
	  if( verbosity  ) {
	      cerr << " CodeAlloc : nb ptr  "<< CodeAlloc::nb << ",  size :"  <<  CodeAlloc::lg << " mpirank: " <<mpirank << endl;
	      if(!mpirank) cerr <<  "Bien: On a fini Normalement" << endl; }
	}
  }

  catch (Error & e) 
    {
      retvalue=e.errcode();
      cerr << "error " << e.what() 
	   << "\n code = "<<  retvalue << " mpirank: " <<mpirank  << endl;
    }
  catch(std::ios_base::failure & e)
    {
     cerr << "std  catch io failure \n what : " << e.what() << endl;; 
     cerr << " at exec line  " << TheCurrentLine << " mpirank: " <<mpirank  << endl; 
    }
  catch(std::exception & e)
    {
     cerr << "std  catch exception \n what : " << e.what() << endl;; 
     cerr << " at exec line  " << TheCurrentLine << " mpirank: " <<mpirank  << endl; 
    
    }
  catch(...)
   {
     cerr << "Strange catch exception ???\n"; 
     cerr << " at exec line  " << TheCurrentLine << " mpirank: " <<mpirank << endl; 
    }
  return retvalue; 
}
static void SetcppIo()
{

#ifdef _WIN32XXXX
  freopen("conin$", "r", stdin);
  freopen("conout$", "w", stdout);
  using namespace __gnu_cxx;
  //  stdio_filebuf<char> * ccout = new stdio_filebuf<char>(stdout, std::ios_base::out);
  static  stdio_filebuf<char> ccout(stdout, std::ios_base::out);
  static  stdio_filebuf<char> ccin(stdin, std::ios_base::in);
   //stdio_filebuf<char> *ccin= new stdio_filebuf<char>(stdin, std::ios_base::in);
   
   cout.rdbuf(&ccout);
   cin.rdbuf(&ccin);
   cerr.rdbuf(&ccout);
   cout << " -- SetcppIo --" << endl; 
#endif
   ios::sync_with_stdio();
}
// pour l'environement.
extern const char *  prognamearg;
int mainff (int  argc, char **argv)
{
  if(argc)  
    prognamearg=argv[0];

    int vvold=verbosity; 
    if(mpirank !=0) verbosity=0;
  SetcppIo();
  GetEnvironment();   
    vvold=verbosity; 
    if(mpirank !=0) verbosity=0; 
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
   //#ifdef HAVE_LIBUMFPACK   
     //if(verbosity)  cout << " UMFPACK ";  
   // #endif
 // callInitsFunct(); Pb opimisation 
  if(verbosity)  cout << endl;
  zzzfff->input(cc);
  EnvironmentLoad(); // just before compile
  verbosity=vvold; 
    
  retvalue= Compile(); 
      
#ifdef PARALLELE
  end_parallele();
#endif
  //  currentblock->close(currentblock).eval(thestack);
  fingraphique();
  if(ThePlotStream) {pclose(ThePlotStream); ThePlotStream=0;}  
  Destroylex( zzzfff);
  
   // ClearMem();
  return retvalue;
}


 

