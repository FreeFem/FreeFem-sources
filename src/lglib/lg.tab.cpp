/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

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
/* Tokens.  */
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

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 118 "lg.y"
{ 
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
/* Line 193 of yacc.c.  */
#line 336 "lg.tab.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 349 "lg.tab.cpp"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
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
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  81
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1023

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  76
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  45
/* YYNRULES -- Number of rules.  */
#define YYNRULES  184
/* YYNRULES -- Number of states.  */
#define YYNSTATES  399

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   305

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
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
static const yytype_uint16 yyprhs[] =
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
     462,   468,   470,   474,   478,   482,   486,   490,   494,   498,
     502,   506,   510,   514,   518,   522,   526,   530,   534,   538,
     542,   546,   548,   550,   554,   560,   561,   563,   565,   567,
     571,   573,   577,   581,   585,   589,   595,   597,   601,   603,
     606,   608,   612,   616,   619,   621,   623,   625,   627,   629,
     634,   639,   646,   650,   654,   658,   663,   667,   672,   676,
     681,   684,   687,   692,   696
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
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
     114,    -1,   114,    74,   114,    75,   114,    -1,   114,    75,
     114,    -1,   114,    75,   114,    75,   114,    -1,   118,    -1,
     114,    22,   114,    -1,   114,    26,   114,    -1,   114,    25,
     114,    -1,   114,    23,   114,    -1,   114,    24,   114,    -1,
     114,    20,   114,    -1,   114,    21,   114,    -1,   114,     9,
     114,    -1,   114,     8,   114,    -1,   114,    12,   114,    -1,
     114,    13,   114,    -1,   114,    10,   114,    -1,   114,    11,
     114,    -1,   114,    16,   114,    -1,   114,    19,   114,    -1,
     114,    17,   114,    -1,   114,    18,   114,    -1,   114,    15,
     114,    -1,   114,    14,   114,    -1,   114,    -1,    75,    -1,
     114,    75,   114,    -1,   114,    75,   114,    75,   114,    -1,
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
static const yytype_uint16 yyrline[] =
{
       0,   249,   249,   289,   292,   293,   296,   297,   298,   299,
     300,   301,   302,   303,   304,   305,   306,   307,   308,   309,
     310,   311,   312,   313,   314,   317,   318,   321,   321,   321,
     321,   323,   324,   325,   327,   334,   335,   336,   337,   338,
     339,   340,   343,   344,   345,   346,   353,   354,   355,   356,
     357,   358,   361,   362,   366,   366,   366,   367,   368,   373,
     374,   376,   377,   379,   380,   384,   387,   388,   391,   391,
     392,   393,   394,   396,   395,   412,   411,   420,   421,   423,
     425,   430,   430,   433,   435,   436,   437,   438,   439,   440,
     441,   442,   446,   447,   448,   449,   451,   453,   456,   460,
     464,   471,   474,   480,   486,   487,   492,   493,   494,   495,
     496,   500,   501,   502,   503,   504,   505,   510,   511,   512,
     513,   517,   518,   519,   520,   521,   522,   523,   524,   525,
     526,   527,   528,   529,   530,   531,   532,   533,   534,   535,
     536,   541,   542,   543,   544,   548,   549,   550,   551,   552,
     553,   554,   555,   556,   557,   558,   561,   562,   566,   567,
     570,   571,   572,   573,   577,   578,   579,   580,   581,   582,
     583,   584,   585,   586,   587,   588,   589,   590,   591,   592,
     593,   594,   595,   604,   605
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
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
static const yytype_uint16 yytoknum[] =
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
static const yytype_uint8 yyr1[] =
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
     111,   112,   112,   112,   112,   112,   112,   113,   113,   113,
     113,   114,   114,   114,   114,   114,   114,   114,   114,   114,
     114,   114,   114,   114,   114,   114,   114,   114,   114,   114,
     114,   115,   115,   115,   115,   116,   116,   116,   116,   116,
     116,   116,   116,   116,   116,   116,   117,   117,   118,   118,
     119,   119,   119,   119,   120,   120,   120,   120,   120,   120,
     120,   120,   120,   120,   120,   120,   120,   120,   120,   120,
     120,   120,   120,   120,   120
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
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
       1,     1,     3,     3,     3,     3,     3,     1,     5,     3,
       5,     1,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     1,     1,     3,     5,     0,     1,     1,     1,     3,
       1,     3,     3,     3,     3,     5,     1,     3,     1,     2,
       1,     3,     3,     2,     1,     1,     1,     1,     1,     4,
       4,     6,     3,     3,     3,     4,     3,     4,     3,     4,
       2,     2,     4,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,   107,   106,   108,   110,   109,     0,     0,   165,
     166,   167,   164,     0,   168,     0,     0,    79,    80,     0,
       0,     0,    83,    42,     0,    54,    55,    56,     0,    84,
      77,     0,     0,     3,    68,    57,     0,     0,    89,     0,
       0,     0,     0,     4,     0,     0,   104,   111,   117,   121,
     158,   160,     0,     0,     0,     0,     0,     0,   156,     0,
       0,    66,     0,    85,    86,    98,    99,     0,     0,     0,
       0,     0,    42,     0,   145,     0,   145,     0,   145,     0,
       0,     1,     2,     5,     0,     0,     0,    46,    61,    63,
      71,     0,     0,     0,     0,     0,    88,   159,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   181,   180,     0,   163,     0,   145,
       0,     0,     0,   183,     0,   184,     0,     0,    70,   100,
       0,     0,     0,     0,     6,     0,   164,   146,   147,   148,
     142,     0,   141,   150,     0,   174,     0,   176,     0,   178,
       0,     0,     0,    96,    31,     0,     0,    27,     0,    28,
      30,    29,     0,    25,     0,     0,     0,    78,    95,    81,
       0,     0,     0,     0,   105,   112,   113,   114,   115,   116,
     130,   129,   133,   134,   131,   132,   140,   139,   135,   137,
     138,   136,   127,   128,   122,   125,   126,   124,   123,     0,
     119,   161,   162,     0,   172,     0,   173,     0,   157,   164,
       0,     0,     0,     0,    35,    67,    45,   182,     0,    43,
       0,     6,    28,    29,     0,     7,     0,     6,     0,     0,
       0,   175,   177,   179,     0,     0,   103,     0,     0,     0,
      69,    58,     0,     0,    49,    48,     0,     0,    62,     0,
       0,     0,     0,     0,     0,     0,   169,     0,   170,    93,
       0,    36,    37,    38,     0,    65,     0,    72,     0,     0,
       9,     0,    11,     0,    75,     0,     0,    13,     0,   149,
     143,   151,   152,   153,     0,   154,     0,     0,    32,     0,
      34,     0,     0,    59,    64,    26,     0,     0,    47,    82,
       0,     0,    92,     0,    87,   118,   120,     0,     0,    39,
      28,    30,    29,     0,    40,    44,    15,    10,    12,     6,
      28,    29,    16,     0,     0,     8,    14,    73,     0,     0,
       0,    97,    33,     0,     0,     0,    51,     0,     0,     0,
       0,   171,    94,     0,     0,     0,    19,     0,    21,     0,
       0,    23,     0,     0,   144,   155,     0,     0,     0,    60,
      50,     0,     0,     0,    41,    17,    20,    22,    18,    24,
       0,     0,   105,     0,    52,     0,     0,     0,    76,     0,
     102,     0,    91,    90,     0,    74,    53,     0,   101
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    31,    32,    33,   234,   172,   151,   165,   223,    34,
      88,   303,    35,    36,   304,    89,    37,    61,    62,    38,
      84,   363,   334,    39,   178,    40,    41,   180,   259,    42,
      43,   314,   162,   163,    44,    45,    46,    47,    48,   153,
     154,    59,    49,    50,    51
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -225
static const yytype_int16 yypact[] =
{
     554,   -22,  -225,  -225,  -225,  -225,  -225,   793,   793,  -225,
    -225,  -225,  -225,     8,  -225,    79,   101,  -225,  -225,   -18,
     -10,   793,  -225,   210,   127,   125,   128,   252,   131,  -225,
    -225,   163,   135,   554,  -225,   187,   190,   181,  -225,   554,
     225,   232,   211,  -225,     3,   910,  -225,    37,    16,  -225,
    -225,   683,   793,   240,   125,   128,   252,    42,  -225,    17,
     244,  -225,     4,  -225,  -225,  -225,  -225,     5,   229,   793,
     238,   177,    35,   265,   664,   292,   664,   295,   664,   298,
     258,  -225,  -225,  -225,   299,   293,     7,    53,  -225,   349,
    -225,   342,   829,   793,   554,   793,  -225,  -225,   793,   793,
     793,   793,   793,   793,   793,   793,   793,   793,   793,   793,
     793,   793,   793,   793,   793,   793,   793,   793,   793,   793,
     793,   793,   793,   793,  -225,  -225,   793,  -225,   793,   664,
      93,   319,    99,  -225,   793,  -225,   865,     8,  -225,  -225,
     348,   102,    18,   793,   215,   332,   361,   224,   234,   235,
    -225,   362,   409,  -225,   104,  -225,   121,  -225,   142,  -225,
     330,   793,   554,  -225,   227,     6,   356,  -225,   336,  -225,
    -225,  -225,    43,  -225,   793,   793,   192,  -225,  -225,  -225,
     304,     9,   155,   440,  -225,  -225,  -225,  -225,  -225,  -225,
     982,   982,   997,   997,   433,   433,   737,   737,   280,   280,
     280,   280,   324,   324,  -225,  -225,  -225,  -225,  -225,   619,
     638,  -225,  -225,   175,  -225,    47,  -225,   554,  -225,   372,
      83,   296,   297,   184,  -225,  -225,  -225,  -225,   321,  -225,
      11,   215,    50,   133,   185,   374,   136,   215,   793,   793,
     707,  -225,  -225,  -225,   380,    55,  -225,   793,   865,   299,
    -225,  -225,   207,   186,   137,  -225,   350,   186,  -225,   299,
     793,   793,   554,   335,   793,   793,  -225,   750,  -225,   394,
     793,  -225,  -225,  -225,   901,  -225,   367,  -225,    61,   186,
    -225,   186,  -225,   237,  -225,   793,   186,  -225,   194,  -225,
     657,   224,   234,   235,   393,  -225,   793,   337,  -225,   202,
    -225,   186,   373,  -225,   401,  -225,   793,   793,  -225,   404,
      12,    13,  -225,   377,  -225,   965,   965,   378,   554,  -225,
     125,   128,   252,   406,  -225,  -225,  -225,  -225,  -225,   215,
     140,   144,   430,   176,   431,  -225,  -225,  -225,   793,   793,
     434,  -225,  -225,    68,   793,   207,  -225,   400,   793,   793,
     405,  -225,  -225,   793,   111,   186,  -225,   186,  -225,   793,
     186,  -225,   793,   368,   965,  -225,   793,   408,   407,  -225,
    -225,   214,   219,   426,  -225,  -225,  -225,  -225,  -225,  -225,
     392,   554,   427,   793,  -225,   554,   554,   429,  -225,   500,
    -225,   435,  -225,  -225,   439,  -225,  -225,   554,  -225
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -225,  -225,  -225,   -37,  -224,   165,   -66,  -195,   230,   -21,
     301,   141,  -225,  -225,  -225,  -225,  -225,   352,  -225,  -225,
    -225,  -225,  -225,  -225,  -225,  -225,  -225,  -225,  -225,  -225,
     -33,  -225,  -225,  -225,    -6,  -225,    -4,  -130,   206,  -125,
     -57,   343,   231,   451,  -225
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -31
static const yytype_int16 yytable[] =
{
      83,    57,    91,    73,    58,   215,   224,   278,    95,   137,
      95,   249,    52,   288,    95,    67,    95,    95,    95,   156,
     173,   158,   134,   228,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,   113,   114,   115,   116,   117,   118,
     119,   120,   121,    98,   255,   256,   132,    95,   253,   167,
      60,    68,   267,    65,   300,   135,   229,   183,    83,   174,
     134,    66,   279,   141,   309,   168,   283,   169,   170,   171,
      70,   179,   213,   253,    96,   138,   139,   250,   235,   133,
     261,   254,   277,   348,   349,   268,   181,   182,   175,   184,
     122,   123,   167,   297,   185,   186,   187,   188,   189,   326,
      99,   100,   101,   102,    95,   354,   367,    95,   289,   240,
     169,   170,   171,     2,     3,   295,   283,    74,   224,    75,
       4,     5,     6,   236,    63,   271,   240,     7,     8,   246,
     218,   214,     9,    10,    11,    12,   217,   230,    14,   227,
     319,   241,   317,   306,   324,   281,    64,   240,   286,   375,
      83,    53,   355,    54,    55,    56,   357,    58,   242,    74,
      95,    75,    76,    81,    77,   235,   280,   282,   150,    71,
     287,   235,   307,    80,   294,   167,   346,   347,   167,   243,
     240,    82,   167,   143,   269,    72,   167,   305,   360,   274,
     283,   173,   262,   169,   170,   171,   169,   170,   171,   283,
     169,   170,   171,    85,   169,   170,   171,   274,   323,   365,
     236,   144,   266,   327,   368,   328,   236,   332,   167,    95,
     336,   275,   284,   374,    95,    86,    68,   257,   167,   312,
     -28,   337,    87,   247,    87,   173,   169,   170,   171,   342,
     -30,   -29,   301,   298,    69,    70,   169,   170,   171,   302,
     231,   385,    90,   391,   310,   311,   386,   167,    74,    92,
      75,   248,   333,   235,   356,   358,    93,   361,    76,    78,
      77,    79,   329,    72,    69,   232,   170,   233,   136,   167,
     152,   335,   152,    94,   152,   352,    78,   140,    79,   376,
     340,   377,   160,   161,   379,    72,   142,   330,   170,   331,
     115,   116,   117,   118,   119,   120,   121,   145,   236,   190,
     191,   192,   193,   194,   195,   196,   197,   198,   199,   200,
     201,   202,   203,   204,   205,   206,   207,   208,   209,   210,
      76,    78,    77,    79,   155,   152,   152,   157,   272,   273,
     159,   164,   371,   372,   389,     1,   117,   118,   119,   120,
     121,   166,   392,   393,   176,   378,    83,   211,   380,   212,
     382,   216,     2,     3,   398,   226,   237,   -27,   238,     4,
       5,     6,   244,   251,   252,   260,     7,     8,   270,   276,
     285,     9,    10,    11,    12,    13,   296,    14,   308,    15,
      16,   313,    17,    18,    19,    20,    21,    22,   318,   339,
      23,    24,    25,    26,    27,   325,   345,   344,   341,   249,
      28,   350,   353,    29,    30,   177,   351,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,   119,   120,   121,   359,   362,   370,   366,
     381,   373,   383,     1,   384,   290,   152,   109,   110,   111,
     112,   113,   114,   115,   116,   117,   118,   119,   120,   121,
       2,     3,   387,   388,   390,   394,   343,     4,     5,     6,
     315,   316,   396,   152,     7,     8,   397,   258,   299,     9,
      10,    11,    12,    13,   239,    14,   369,    15,    16,   225,
      17,    18,    19,    20,    21,    22,    97,     0,    23,    24,
      25,    26,    27,     1,   245,     0,     0,     0,    28,     0,
       0,    29,    30,   263,     0,     0,     0,     0,     0,     0,
       2,     3,     0,     0,     0,     0,     0,     4,     5,     6,
       0,     0,     0,     0,     7,     8,     0,     0,     0,     9,
      10,    11,    12,    13,   364,    14,     0,    15,    16,     0,
      17,    18,    19,    20,    21,    22,     0,     1,    23,    24,
      25,    26,    27,     0,     0,     0,     0,     0,    28,     0,
       0,    29,    30,   395,     2,     3,     0,     0,     0,     0,
       0,     4,     5,     6,     0,     0,     0,     0,     7,     8,
       0,     0,     0,     9,    10,    11,    12,    13,     0,    14,
       0,    15,    16,     0,    17,    18,    19,    20,    21,    22,
       0,     0,    23,    24,    25,    26,    27,     0,     0,     0,
       0,     0,    28,     0,     0,    29,    30,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,   119,   120,   121,   103,   104,   105,   106,
     107,   108,   109,   110,   111,   112,   113,   114,   115,   116,
     117,   118,   119,   120,   121,   103,   104,   105,   106,   107,
     108,   109,   110,   111,   112,   113,   114,   115,   116,   117,
     118,   119,   120,   121,     2,     3,     0,     0,     0,     0,
       0,     4,     5,     6,   264,     0,     0,     0,     7,     8,
       0,     0,     0,     9,    10,    11,   146,     0,     0,    14,
       0,   124,   125,   265,   126,   127,   128,   129,   130,   131,
       0,     0,    53,     0,   147,   148,   149,     2,     3,     0,
       0,     0,   338,     0,     4,     5,     6,     0,     0,   150,
       0,     7,     8,     0,     0,     0,     9,    10,    11,   146,
       0,     0,    14,   111,   112,   113,   114,   115,   116,   117,
     118,   119,   120,   121,     0,    53,     0,   291,   292,   293,
       2,     3,     0,     0,     0,     0,     0,     4,     5,     6,
       0,     0,   150,     0,     7,     8,     0,     0,     0,     9,
      10,    11,    12,     0,     0,    14,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    53,     0,
      54,    55,    56,     2,     3,     0,     0,     0,     0,     0,
       4,     5,     6,     0,     0,   150,     0,     7,     8,     0,
       0,     0,     9,    10,    11,    12,     0,     0,    14,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     2,
       3,    53,     0,    54,    55,    56,     4,     5,     6,     0,
       0,     0,     0,     7,     8,     0,     0,     0,     9,    10,
      11,    12,     0,     0,    14,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     2,     3,    23,     0,    54,
      55,    56,     4,     5,     6,     0,     0,     0,     0,     7,
       8,     0,     0,     0,     9,    10,    11,   219,     0,     0,
      14,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     2,     3,    53,     0,   220,   221,   222,     4,     5,
       6,     0,     0,     0,     0,     7,     8,     0,     0,     0,
       9,    10,    11,   146,     7,     8,    14,     0,     0,     9,
      10,    11,    12,     0,     0,    14,     0,     0,     0,    53,
       0,   320,   321,   322,     0,     0,     0,     0,    53,     0,
      54,    55,    56,   103,   104,   105,   106,   107,   108,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     120,   121,   105,   106,   107,   108,   109,   110,   111,   112,
     113,   114,   115,   116,   117,   118,   119,   120,   121,   107,
     108,   109,   110,   111,   112,   113,   114,   115,   116,   117,
     118,   119,   120,   121
};

static const yytype_int16 yycheck[] =
{
      33,     7,    39,    24,     8,   130,   136,   231,     5,     5,
       5,     5,    34,   237,     5,    21,     5,     5,     5,    76,
      86,    78,     5,     5,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,     6,   174,   175,    52,     5,     5,    42,
      42,    16,     5,    71,   249,    38,    38,    94,    91,     6,
       5,    71,    12,    69,   259,    58,     5,    60,    61,    62,
      35,    92,   129,     5,    71,    71,    71,    71,   144,    37,
      71,    38,    71,    71,    71,    38,    92,    93,    35,    95,
      74,    75,    42,    38,    98,    99,   100,   101,   102,    38,
      63,    64,    65,    66,     5,   329,    38,     5,   238,     5,
      60,    61,    62,    20,    21,   240,     5,    34,   248,    36,
      27,    28,    29,   144,    45,    42,     5,    34,    35,   162,
     134,    38,    39,    40,    41,    42,    37,   143,    45,    37,
     270,    37,   267,     6,   274,    12,    45,     5,    12,    38,
     183,    58,    12,    60,    61,    62,    12,   161,    37,    34,
       5,    36,    34,     0,    36,   231,   232,   233,    75,    42,
     236,   237,    35,    42,   240,    42,   306,   307,    42,    37,
       5,    46,    42,     6,   217,    58,    42,   253,    12,     5,
       5,   257,    37,    60,    61,    62,    60,    61,    62,     5,
      60,    61,    62,    16,    60,    61,    62,     5,   274,   339,
     231,    34,    37,   279,   344,   281,   237,   283,    42,     5,
     286,    37,    37,   353,     5,    35,    16,    35,    42,   262,
       6,    37,    42,     6,    42,   301,    60,    61,    62,    37,
       6,     6,    35,   247,    34,    35,    60,    61,    62,    42,
      35,    37,    71,   383,   260,   261,    37,    42,    34,    34,
      36,    34,   283,   329,   330,   331,    34,   333,    34,    34,
      36,    36,    35,    58,    34,    60,    61,    62,    34,    42,
      74,   285,    76,    72,    78,   318,    34,    58,    36,   355,
     296,   357,    34,    35,   360,    58,    58,    60,    61,    62,
      20,    21,    22,    23,    24,    25,    26,    42,   329,   103,
     104,   105,   106,   107,   108,   109,   110,   111,   112,   113,
     114,   115,   116,   117,   118,   119,   120,   121,   122,   123,
      34,    34,    36,    36,    42,   129,   130,    42,    42,    42,
      42,    42,   348,   349,   381,     3,    22,    23,    24,    25,
      26,    58,   385,   386,     5,   359,   389,   126,   362,   128,
     366,    42,    20,    21,   397,    17,    34,     6,     6,    27,
      28,    29,    42,    17,    38,    71,    34,    35,     6,    58,
       6,    39,    40,    41,    42,    43,     6,    45,    38,    47,
      48,    56,    50,    51,    52,    53,    54,    55,     4,     6,
      58,    59,    60,    61,    62,    38,     5,    34,    71,     5,
      68,    34,     6,    71,    72,    73,    38,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,     6,     6,    38,     5,
      72,    36,    34,     3,    37,   239,   240,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      20,    21,    36,    71,    37,    36,   301,    27,    28,    29,
     264,   265,    37,   267,    34,    35,    37,   176,   248,    39,
      40,    41,    42,    43,    75,    45,   345,    47,    48,   137,
      50,    51,    52,    53,    54,    55,    45,    -1,    58,    59,
      60,    61,    62,     3,   161,    -1,    -1,    -1,    68,    -1,
      -1,    71,    72,    73,    -1,    -1,    -1,    -1,    -1,    -1,
      20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,
      -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,
      40,    41,    42,    43,   338,    45,    -1,    47,    48,    -1,
      50,    51,    52,    53,    54,    55,    -1,     3,    58,    59,
      60,    61,    62,    -1,    -1,    -1,    -1,    -1,    68,    -1,
      -1,    71,    72,    73,    20,    21,    -1,    -1,    -1,    -1,
      -1,    27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,
      -1,    -1,    -1,    39,    40,    41,    42,    43,    -1,    45,
      -1,    47,    48,    -1,    50,    51,    52,    53,    54,    55,
      -1,    -1,    58,    59,    60,    61,    62,    -1,    -1,    -1,
      -1,    -1,    68,    -1,    -1,    71,    72,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    20,    21,    -1,    -1,    -1,    -1,
      -1,    27,    28,    29,    75,    -1,    -1,    -1,    34,    35,
      -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,
      -1,    28,    29,    75,    31,    32,    33,    34,    35,    36,
      -1,    -1,    58,    -1,    60,    61,    62,    20,    21,    -1,
      -1,    -1,    75,    -1,    27,    28,    29,    -1,    -1,    75,
      -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      -1,    -1,    45,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    -1,    58,    -1,    60,    61,    62,
      20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,
      -1,    -1,    75,    -1,    34,    35,    -1,    -1,    -1,    39,
      40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    58,    -1,
      60,    61,    62,    20,    21,    -1,    -1,    -1,    -1,    -1,
      27,    28,    29,    -1,    -1,    75,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    20,
      21,    58,    -1,    60,    61,    62,    27,    28,    29,    -1,
      -1,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,
      41,    42,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    20,    21,    58,    -1,    60,
      61,    62,    27,    28,    29,    -1,    -1,    -1,    -1,    34,
      35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    20,    21,    58,    -1,    60,    61,    62,    27,    28,
      29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    34,    35,    45,    -1,    -1,    39,
      40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,    58,
      -1,    60,    61,    62,    -1,    -1,    -1,    -1,    58,    -1,
      60,    61,    62,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
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
      25,    26,    74,    75,    28,    29,    31,    32,    33,    34,
      35,    36,   110,    37,     5,    38,    34,     5,    71,    71,
      58,   110,    58,     6,    34,    42,    42,    60,    61,    62,
      75,    82,   114,   115,   116,    42,   116,    42,   116,    42,
      34,    35,   108,   109,    42,    83,    58,    42,    58,    60,
      61,    62,    81,    82,     6,    35,     5,    73,   100,    85,
     103,   110,   110,    79,   110,   112,   112,   112,   112,   112,
     114,   114,   114,   114,   114,   114,   114,   114,   114,   114,
     114,   114,   114,   114,   114,   114,   114,   114,   114,   114,
     114,   118,   118,   116,    38,   115,    42,    37,   112,    42,
      60,    61,    62,    84,   113,    93,    17,    37,     5,    38,
     110,    35,    60,    62,    80,    82,    85,    34,     6,    75,
       5,    37,    37,    37,    42,   117,   106,     6,    34,     5,
      71,    17,    38,     5,    38,   113,   113,    35,    86,   104,
      71,    71,    37,    73,    75,    75,    37,     5,    38,   106,
       6,    42,    42,    42,     5,    37,    58,    71,    80,    12,
      82,    12,    82,     5,    37,     6,    12,    82,    80,   113,
     114,    60,    61,    62,    82,   115,     6,    38,   112,    84,
      83,    35,    42,    87,    90,    82,     6,    35,    38,    83,
     110,   110,   106,    56,   107,   114,   114,   115,     4,   113,
      60,    61,    62,    82,   113,    38,    38,    82,    82,    35,
      60,    62,    82,    85,    98,   112,    82,    37,    75,     6,
     110,    71,    37,    81,    34,     5,   113,   113,    71,    71,
      34,    38,   106,     6,    80,    12,    82,    12,    82,     6,
      12,    82,     6,    97,   114,   113,     5,    38,   113,    87,
      38,   110,   110,    36,   113,    38,    82,    82,   112,    82,
     112,    72,   110,    34,    37,    37,    37,    36,    71,    79,
      37,   113,   106,   106,    36,    73,    37,    37,   106
};

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
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
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
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
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
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

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
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
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
  YYUSE (yyvaluep);

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
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
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
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

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
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
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

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

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

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

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
		        const char *  magicffglut="#!ffglutdata2..\n";
                        if(ThePlotStream) fwrite(magicffglut,strlen(magicffglut),1,ThePlotStream);	            
                        size_t sizestack = currentblock->size()+1024 ; //  before close 
                        (yyvsp[(1) - (2)].cinst)+=currentblock->close(currentblock);
                        if(verbosity) cout << " sizestack + 1024 =" << sizestack << "  ( " << sizestack-1024 <<" )\n" ;   
                        size_t lg0,lg1;                       
                        int NbPtr = ShowAlloc("init execution ",lg0); // number of un delele ptr
                        if(verbosity) cout << endl;  
                        { Stack stack = newStack(sizestack);
                        double CPUcompile= CPUtime();
                        try {                  
                          (yyvsp[(1) - (2)].cinst).eval(stack);}
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
    {(yyval.cinst)=(yyvsp[(1) - (1)].cexp);;;;}
    break;

  case 5:
#line 293 "lg.y"
    { (yyval.cinst)= ((yyvsp[(1) - (2)].cinst)+=(yyvsp[(2) - (2)].cexp)) ;}
    break;

  case 6:
#line 296 "lg.y"
    { (yyval.clist_id)=new ListOfId();;}
    break;

  case 7:
#line 297 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str)));}
    break;

  case 8:
#line 298 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (3)].str),(yyvsp[(3) - (3)].cexp))) ;}
    break;

  case 9:
#line 299 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,2> **>()));}
    break;

  case 10:
#line 300 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,2> **>(),true));}
    break;

  case 11:
#line 301 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,3> **>()));}
    break;

  case 12:
#line 302 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,3> **>(),true));}
    break;

  case 13:
#line 303 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),C_F0(),(yyvsp[(1) - (2)].type)->right())) ;}
    break;

  case 14:
#line 304 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),C_F0(),(yyvsp[(1) - (3)].type),true)) ;}
    break;

  case 15:
#line 305 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (3)].clist_id))) ;}
    break;

  case 16:
#line 306 "lg.y"
    { (yyval.clist_id) = (yyvsp[(1) - (3)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str))) ;}
    break;

  case 17:
#line 307 "lg.y"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (5)].clist_id))) ;}
    break;

  case 18:
#line 308 "lg.y"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (5)].str),(yyvsp[(5) - (5)].cexp))) ;}
    break;

  case 19:
#line 309 "lg.y"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,2> **>())) ;}
    break;

  case 20:
#line 310 "lg.y"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,2> **>(),true)) ;}
    break;

  case 21:
#line 311 "lg.y"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,3> **>())) ;}
    break;

  case 22:
#line 312 "lg.y"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,3> **>(),true)) ;}
    break;

  case 23:
#line 313 "lg.y"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),C_F0(),(yyvsp[(3) - (4)].type)->right())) ;}
    break;

  case 24:
#line 314 "lg.y"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),C_F0(),(yyvsp[(3) - (5)].type),true)) ;}
    break;

  case 25:
#line 317 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str))); ;}
    break;

  case 26:
#line 318 "lg.y"
    { (yyval.clist_id)=(yyvsp[(1) - (3)].clist_id)  ; (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str))); ;}
    break;

  case 31:
#line 323 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (1)].str),dcltype);}
    break;

  case 32:
#line 324 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (3)].str),dcltype,(yyvsp[(3) - (3)].cexp));}
    break;

  case 33:
#line 325 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (4)].str),dcltype,(yyvsp[(3) - (4)].args));
                                              (yyvsp[(3) - (4)].args).destroy();}
    break;

  case 34:
#line 327 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 35:
#line 334 "lg.y"
    {(yyval.args)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 36:
#line 335 "lg.y"
    {(yyval.args)=Find((yyvsp[(1) - (2)].str));}
    break;

  case 37:
#line 336 "lg.y"
    {(yyval.args)=Find((yyvsp[(1) - (2)].str));}
    break;

  case 38:
#line 337 "lg.y"
    {(yyval.args)=Find((yyvsp[(1) - (2)].str));}
    break;

  case 39:
#line 338 "lg.y"
    { (yyval.args)=make_pair<const char *,const C_F0>((yyvsp[(1) - (3)].str),(yyvsp[(3) - (3)].cexp));}
    break;

  case 40:
#line 339 "lg.y"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp)) ;}
    break;

  case 41:
#line 340 "lg.y"
    { (yyval.args)= ((yyvsp[(1) - (5)].args)+= make_pair<const char *,const C_F0>((yyvsp[(3) - (5)].str),(yyvsp[(5) - (5)].cexp)));}
    break;

  case 43:
#line 344 "lg.y"
    {(yyval.type)=TypeArray((yyvsp[(1) - (4)].type),(yyvsp[(3) - (4)].type));}
    break;

  case 44:
#line 345 "lg.y"
    {(yyval.type)=TypeArray((yyvsp[(1) - (6)].type),(yyvsp[(3) - (6)].type),(yyvsp[(5) - (6)].type));}
    break;

  case 45:
#line 346 "lg.y"
    {(yyval.type)=TypeTemplate((yyvsp[(1) - (4)].type),(yyvsp[(3) - (4)].type));}
    break;

  case 46:
#line 353 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(1) - (1)].str),currentblock,fespacetype,fespacecomplex,fespacedim); ;}
    break;

  case 47:
#line 354 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[(1) - (4)].str),currentblock,fespacetype,(yyvsp[(3) - (4)].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 48:
#line 355 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(1) - (3)].str),currentblock,fespacetype,(yyvsp[(3) - (3)].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 49:
#line 356 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(2) - (3)].clist_id),currentblock,fespacetype,fespacecomplex,fespacedim) ;}
    break;

  case 50:
#line 357 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[(2) - (6)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (6)].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 51:
#line 358 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(2) - (5)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (5)].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 52:
#line 361 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[(1) - (4)].str),currentblock,fespacetype,(yyvsp[(3) - (4)].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 53:
#line 362 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[(2) - (6)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (6)].cexp),fespacecomplex,fespacedim) ;}
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
    {fespacecomplex=false;  fespacetype = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 58:
#line 368 "lg.y"
    {
             if ((yyvsp[(3) - (4)].type) != typevarreal && (yyvsp[(3) - (4)].type) != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=((yyvsp[(3) - (4)].type)==typevarcomplex);
             fespacetype = Find((yyvsp[(1) - (4)].str));;}
    break;

  case 59:
#line 373 "lg.y"
    {  (yyval.cexp) = (yyvsp[(1) - (1)].cexp)  ;}
    break;

  case 60:
#line 374 "lg.y"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 61:
#line 376 "lg.y"
    {  (yyval.cexp) = (yyvsp[(1) - (1)].cexp)  ;}
    break;

  case 62:
#line 377 "lg.y"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 63:
#line 379 "lg.y"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[(2) - (2)].cexp);}
    break;

  case 64:
#line 380 "lg.y"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[(5) - (5)].cexp);}
    break;

  case 65:
#line 384 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,size_t>((yyvsp[(1) - (4)].str),typeFESpace((yyvsp[(3) - (4)].args)),(yyvsp[(3) - (4)].args),dimFESpaceImage((yyvsp[(3) - (4)].args)));
     (yyvsp[(3) - (4)].args).destroy(); ;}
    break;

  case 67:
#line 388 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 68:
#line 391 "lg.y"
    {dcltype=(yyvsp[(1) - (1)].type);}
    break;

  case 69:
#line 391 "lg.y"
    {(yyval.cexp)=(yyvsp[(3) - (4)].cexp);}
    break;

  case 70:
#line 392 "lg.y"
    {(yyval.cexp)=(yyvsp[(2) - (3)].cexp);}
    break;

  case 71:
#line 393 "lg.y"
    { (yyval.cexp)=(yyvsp[(1) - (2)].cexp);}
    break;

  case 72:
#line 394 "lg.y"
    {(yyval.cexp)=currentblock->NewID((yyvsp[(1) - (5)].type),(yyvsp[(2) - (5)].str),(yyvsp[(4) - (5)].cexp));;}
    break;

  case 73:
#line 396 "lg.y"
    {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = (yyvsp[(2) - (6)].type)->right();
                      routineinblock[kkembtype] = currentblock;
                      (yyvsp[(5) - (6)].routine)=new Routine((yyvsp[(1) - (6)].type),(yyvsp[(2) - (6)].type)->right(),(yyvsp[(3) - (6)].str),(yyvsp[(5) - (6)].clist_id),currentblock);
		      // routineinblock[kkembtype]->Add($3,"(",$<routine>5); //pas recursif pour l'instanat test  FH 27 dec 2008
                     // cout << " \n after new routine \n " << endl;                      
                      ;}
    break;

  case 74:
#line 405 "lg.y"
    { currentblock=(yyvsp[(5) - (10)].routine)->Set((yyvsp[(9) - (10)].cinst));
                       currentblock->Add((yyvsp[(3) - (10)].str),"(",(yyvsp[(5) - (10)].routine)); //pas recursif pour l'instant test  FH 27 dec 2008
                       kkembtype--;
                       (yyval.cexp)=0;
                    
                        ;}
    break;

  case 75:
#line 412 "lg.y"
    {Block::open(currentblock); (yyvsp[(1) - (5)].type)->SetArgs((yyvsp[(4) - (5)].clist_id));;}
    break;

  case 76:
#line 414 "lg.y"
    {  (yyval.cinst)=currentblock->close(currentblock);
                         (yyval.cexp)=currentblock->NewID((yyvsp[(1) - (9)].type),(yyvsp[(2) - (9)].str),(yyvsp[(8) - (9)].cexp),*(yyvsp[(4) - (9)].clist_id));
                         delete (yyvsp[(4) - (9)].clist_id); //  FH 23032005
                         ;}
    break;

  case 77:
#line 420 "lg.y"
    {  Block::open(currentblock);}
    break;

  case 78:
#line 421 "lg.y"
    {  (yyval.cexp)=currentblock->close(currentblock);}
    break;

  case 79:
#line 423 "lg.y"
    {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 80:
#line 425 "lg.y"
    {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 81:
#line 430 "lg.y"
    {dcltype=(yyvsp[(1) - (1)].type); Block::open(currentblock);  ;}
    break;

  case 82:
#line 431 "lg.y"
    {(yyval.cexp)=(yyvsp[(3) - (3)].cexp);}
    break;

  case 83:
#line 433 "lg.y"
    { Block::open(currentblock) ;}
    break;

  case 84:
#line 435 "lg.y"
    {(yyval.cexp)=0;;}
    break;

  case 85:
#line 436 "lg.y"
    {zzzfff->input((yyvsp[(2) - (2)].str));(yyval.cexp)= 0; ;}
    break;

  case 86:
#line 437 "lg.y"
    {load((yyvsp[(2) - (2)].str));(yyval.cexp)= 0; ;}
    break;

  case 87:
#line 438 "lg.y"
    {(yyval.cexp)=Try((yyvsp[(3) - (5)].cinst),(yyvsp[(5) - (5)].cexp),currentblock->close(currentblock));;}
    break;

  case 88:
#line 439 "lg.y"
    {(yyval.cexp)=(yyvsp[(1) - (2)].cexp);}
    break;

  case 89:
#line 440 "lg.y"
    {(yyval.cexp)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 90:
#line 441 "lg.y"
    {inloopcount--; (yyval.cexp)=For((yyvsp[(3) - (9)].cexp),(yyvsp[(5) - (9)].cexp),(yyvsp[(7) - (9)].cexp),(yyvsp[(9) - (9)].cexp));}
    break;

  case 91:
#line 443 "lg.y"
    {inloopcount--; 
                (yyval.cexp)=C_F0(For((yyvsp[(3) - (9)].cexp),(yyvsp[(5) - (9)].cexp),(yyvsp[(7) - (9)].cexp),(yyvsp[(9) - (9)].cexp)),currentblock->close(currentblock));}
    break;

  case 92:
#line 446 "lg.y"
    {inloopcount--;(yyval.cexp)=While((yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 93:
#line 447 "lg.y"
    {(yyval.cexp)=FIf((yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 94:
#line 448 "lg.y"
    {(yyval.cexp)=FIf((yyvsp[(3) - (7)].cexp),(yyvsp[(5) - (7)].cexp),(yyvsp[(7) - (7)].cexp));}
    break;

  case 95:
#line 449 "lg.y"
    { 
                      (yyval.cexp)=C_F0(new E_block((yyvsp[(2) - (3)].cinst),(yyvsp[(3) - (3)].cexp)),atype<void>()) ;}
    break;

  case 96:
#line 451 "lg.y"
    {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[(2) - (3)].str),C_F0(TheOperators,"[border]",(yyvsp[(3) - (3)].args)));}
    break;

  case 97:
#line 453 "lg.y"
    {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[(2) - (6)].str),C_F0(TheOperators,"[border]",(yyvsp[(4) - (6)].args)));}
    break;

  case 98:
#line 456 "lg.y"
    {
                    if(inloopcount) 
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;}
    break;

  case 99:
#line 460 "lg.y"
    { 
                    if(inloopcount)
                        (yyval.cexp)= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");}
    break;

  case 100:
#line 464 "lg.y"
    { 
                    if (kkembtype>=0)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_return,(rettype[kkembtype]->CastTo((yyvsp[(2) - (3)].cexp))).OnReturn()) ,atype<void>());
                     else lgerror(" return not in routine ") ;}
    break;

  case 101:
#line 471 "lg.y"
    {(yyval.cexp) =  (yyvsp[(7) - (7)].cexp); ;}
    break;

  case 102:
#line 474 "lg.y"
    { 
   Block::open(currentblock);
   (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[(2) - (7)].str),atype<double*>());
   (yyval.args)+= (yyvsp[(4) - (7)].cexp);
   (yyval.args)+= (yyvsp[(6) - (7)].cexp) ;}
    break;

  case 103:
#line 480 "lg.y"
    {   
   (yyval.args) = ((yyvsp[(1) - (2)].args) += (yyvsp[(2) - (2)].cexp));
   currentblock->close(currentblock);}
    break;

  case 105:
#line 487 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 112:
#line 501 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 113:
#line 502 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"+=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 114:
#line 503 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"-=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 115:
#line 504 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"*=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 116:
#line 505 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"/=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 118:
#line 511 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"?:",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 119:
#line 512 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 120:
#line 513 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 122:
#line 518 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 123:
#line 519 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 124:
#line 520 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 125:
#line 521 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 126:
#line 522 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 127:
#line 523 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 128:
#line 524 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 129:
#line 525 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 130:
#line 526 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 131:
#line 527 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 132:
#line 528 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 133:
#line 529 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 134:
#line 530 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 135:
#line 531 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 136:
#line 532 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 137:
#line 533 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 138:
#line 534 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 139:
#line 535 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 140:
#line 536 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 141:
#line 541 "lg.y"
    {(yyval.cexp)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 142:
#line 542 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,":");}
    break;

  case 143:
#line 543 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 144:
#line 544 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 145:
#line 548 "lg.y"
    {(yyval.args)=0;}
    break;

  case 146:
#line 549 "lg.y"
    {(yyval.args)=Find((yyvsp[(1) - (1)].str));}
    break;

  case 147:
#line 550 "lg.y"
    {(yyval.args)=Find((yyvsp[(1) - (1)].str));}
    break;

  case 148:
#line 551 "lg.y"
    {(yyval.args)=Find((yyvsp[(1) - (1)].str));}
    break;

  case 149:
#line 552 "lg.y"
    { (yyval.args)=make_pair<const char *,const C_F0>((yyvsp[(1) - (3)].str),(yyvsp[(3) - (3)].cexp));}
    break;

  case 150:
#line 553 "lg.y"
    {(yyval.args)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 151:
#line 554 "lg.y"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str))) ;}
    break;

  case 152:
#line 555 "lg.y"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str))) ;}
    break;

  case 153:
#line 556 "lg.y"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str))) ;}
    break;

  case 154:
#line 557 "lg.y"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp)) ;}
    break;

  case 155:
#line 558 "lg.y"
    { (yyval.args)= ((yyvsp[(1) - (5)].args)+= make_pair<const char *,const C_F0>((yyvsp[(3) - (5)].str),(yyvsp[(5) - (5)].cexp))) ;}
    break;

  case 156:
#line 561 "lg.y"
    {(yyval.args)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 157:
#line 562 "lg.y"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp)) ;}
    break;

  case 159:
#line 567 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(1) - (2)].oper),(yyvsp[(2) - (2)].cexp));}
    break;

  case 161:
#line 571 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 162:
#line 572 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 163:
#line 573 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));}
    break;

  case 164:
#line 577 "lg.y"
    {(yyval.cexp)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 165:
#line 578 "lg.y"
    {(yyval.cexp)= CConstant((yyvsp[(1) - (1)].lnum));}
    break;

  case 166:
#line 579 "lg.y"
    {(yyval.cexp)= CConstant((yyvsp[(1) - (1)].dnum));}
    break;

  case 167:
#line 580 "lg.y"
    {(yyval.cexp)= CConstant(complex<double>(0,(yyvsp[(1) - (1)].dnum)));}
    break;

  case 168:
#line 581 "lg.y"
    {(yyval.cexp)= CConstant<const char *>((yyvsp[(1) - (1)].str));}
    break;

  case 169:
#line 582 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (4)].cexp),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args));;}
    break;

  case 170:
#line 583 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (4)].cexp),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].cexp));}
    break;

  case 171:
#line 584 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (6)].cexp),(yyvsp[(2) - (6)].oper),(yyvsp[(3) - (6)].cexp),(yyvsp[(5) - (6)].cexp));}
    break;

  case 172:
#line 585 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),"[]");}
    break;

  case 173:
#line 586 "lg.y"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 174:
#line 587 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 175:
#line 588 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 176:
#line 589 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 177:
#line 590 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 178:
#line 591 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 179:
#line 592 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 180:
#line 593 "lg.y"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));}
    break;

  case 181:
#line 594 "lg.y"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));}
    break;

  case 182:
#line 595 "lg.y"
    {
             if ((yyvsp[(1) - (4)].type)->right()->CastingFrom((yyvsp[(3) - (4)].cexp).left()) ) 
                (yyval.cexp)=(yyvsp[(1) - (4)].type)->right()->CastTo((yyvsp[(3) - (4)].cexp))  ;
             else { (yyval.cexp)=(yyvsp[(1) - (4)].type)->right()->Find("<--",basicAC_F0_wa((yyvsp[(3) - (4)].cexp)));
             if (!(yyval.cexp).left()) { cerr << " no wait to change " << (yyvsp[(3) - (4)].cexp).left()->right()->name() << " in " << 
                                        (yyvsp[(1) - (4)].type)->right()->name() << endl;
                                CompileError(" Error in type(exp) "); }
             }
            ;}
    break;

  case 183:
#line 604 "lg.y"
    {(yyval.cexp)=(yyvsp[(2) - (3)].cexp);}
    break;

  case 184:
#line 605 "lg.y"
    { (yyval.cexp)=C_F0(TheOperators,"[]",(yyvsp[(2) - (3)].args));}
    break;


/* Line 1267 of yacc.c.  */
#line 2937 "lg.tab.cpp"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
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
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
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

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
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


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
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
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 610 "lg.y"
 


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


 

