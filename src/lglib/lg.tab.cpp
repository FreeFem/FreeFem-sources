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
     DOTMULEQ = 302,
     DOTDIVEQ = 303,
     ARROW = 304,
     BORDER = 305,
     CURVE = 306,
     SOLVE = 307
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
#define DOTMULEQ 302
#define DOTDIVEQ 303
#define ARROW 304
#define BORDER 305
#define CURVE 306
#define SOLVE 307




/* Copy the first part of user declarations.  */
#line 3 "lg.ypp"

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
#include "InitFunct.hpp"
#ifdef __MWERKS__
#ifdef __INTEL__
#include <malloc.h>
#else
#include <alloca.h>
#endif
#endif
#include "RNM.hpp"

#include "AFunction.hpp"
//  to reserve space to graphical pointer function
#include "rgraph.hpp"
#include "fem.hpp"
#include "FESpacen.hpp"
#include "FESpace.hpp"
#include "MeshPoint.hpp"

#include "lgfem.hpp"
#include "lex.hpp"
#include "environment.hpp"
extern long storageused();
    extern FILE *ThePlotStream;
    extern KN<String> *pkarg;

class Routine;
bool load(string s);

 template <class R,int d> class FE;
 template <class R,int d,int i> class FE_;

extern mylex *zzzfff;
// modif FH for window to have 1 dll  for mpi and none mpi ..
extern  void (*initparallele)(int &, char **&);
extern  void (*init_lgparallele)();
// extern  void (*end_parallele)();
//
#ifdef HAVE_LIBARPACK
  void init_eigenvalue();
#endif

  aType dcltype;
const int nbembtype=10;
aType rettype[nbembtype];
Block * routineinblock[nbembtype]; // Add FH july 2005 pb clean on return
int kkembtype=-1;
int inloopcount=0;

/// <<currentblock>> Block class from [[file:../fflib/AFunction.hpp::Block]]

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
extern int UnShowAlloc;
int ShowAlloc(const char *s,size_t &);
// <<yylex>> Connection from grammar to lexer object zzzfff [[file:../fflib/lex.hpp::zzzfff]] of class mylex
// [[file:../fflib/lex.hpp::class mylex]]. Method mylex::scan() is implemented at [[file:../fflib/lex.cpp::mylex_scan]]

inline int yylex()  {return zzzfff->scan();}
inline int lineno() {return zzzfff->lineno();}

extern bool withrgraphique;

/// <<fingraphique>>

inline void fingraphique()
 { if(withrgraphique)
   { withrgraphique=false;
    rattente(1);
    closegraphique();
  }}

void lgerror (const char* s) ;


 // mpi ptr to function ...
void (*initparallele)(int &argc, char **& argv)=0 ;
void (*init_lgparallele)()=0;
//void (*end_parallele)()=0;

// Add dec 2014
#include <vector>
typedef void (*AtEnd)();
vector<AtEnd> AtFFEnd;
void ff_finalize()
{
    for (vector<AtEnd>::const_reverse_iterator i=AtFFEnd.rbegin(); i !=AtFFEnd.rend(); ++ i)
    (**i)();
    AtFFEnd.clear();
}
void ff_atend(AtEnd f)
{
    AtFFEnd.push_back(f);
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

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 156 "lg.ypp"
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
#line 373 "lg.tab.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 386 "lg.tab.cpp"

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
#define YYLAST   1093

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  78
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  47
/* YYNRULES -- Number of rules.  */
#define YYNRULES  196
/* YYNRULES -- Number of states.  */
#define YYNSTATES  430

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   307

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
       2,     2,     2,     2,     2,     2,     2,     2,    76,    73,
      16,     6,    17,    77,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    35,     2,    38,    31,    33,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    74,    10,    75,     2,     2,     2,     2,
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
      65,    66,    67,    68,    69,    70,    71,    72
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
     138,   142,   148,   150,   155,   163,   170,   180,   185,   193,
     203,   205,   210,   214,   218,   225,   231,   236,   243,   245,
     247,   249,   251,   256,   258,   262,   264,   268,   271,   277,
     282,   284,   288,   289,   294,   298,   301,   307,   308,   319,
     320,   330,   332,   334,   336,   338,   339,   343,   345,   347,
     351,   357,   359,   362,   365,   371,   374,   376,   377,   386,
     396,   406,   412,   418,   426,   430,   434,   441,   444,   447,
     451,   459,   467,   477,   480,   482,   486,   488,   490,   492,
     494,   496,   498,   502,   506,   510,   514,   518,   522,   526,
     528,   534,   538,   544,   546,   550,   554,   558,   562,   566,
     570,   574,   578,   582,   586,   590,   594,   598,   602,   606,
     610,   614,   618,   622,   624,   626,   630,   636,   637,   639,
     641,   643,   647,   649,   653,   657,   661,   665,   671,   673,
     677,   679,   682,   684,   688,   692,   695,   697,   699,   701,
     703,   705,   710,   715,   722,   726,   730,   734,   739,   743,
     748,   752,   757,   760,   763,   768,   772
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      79,     0,    -1,    80,    46,    -1,    81,    -1,   109,    -1,
      81,   109,    -1,    -1,    84,    -1,    84,     6,   116,    -1,
      60,    84,    -1,    60,    12,    84,    -1,    62,    84,    -1,
      62,    12,    84,    -1,    87,    84,    -1,    87,    12,    84,
      -1,    35,    82,    38,    -1,    82,     5,    84,    -1,    82,
       5,    35,    82,    38,    -1,    82,     5,    84,     6,   116,
      -1,    82,     5,    60,    84,    -1,    82,     5,    60,    12,
      84,    -1,    82,     5,    62,    84,    -1,    82,     5,    62,
      12,    84,    -1,    82,     5,    87,    84,    -1,    82,     5,
      87,    12,    84,    -1,    84,    -1,    83,     5,    84,    -1,
      42,    -1,    60,    -1,    62,    -1,    61,    -1,    42,    -1,
      42,     6,   116,    -1,    42,    34,    86,    37,    -1,    85,
       5,    85,    -1,   117,    -1,    60,    42,    -1,    61,    42,
      -1,    62,    42,    -1,    42,     6,   117,    -1,    86,     5,
     117,    -1,    86,     5,    84,     6,   117,    -1,    58,    -1,
      58,    35,    58,    38,    -1,    58,    35,    58,    38,    35,
      58,    38,    -1,    58,    35,    58,     5,    58,    38,    -1,
      58,    35,    58,     5,    58,    38,    35,    58,    38,    -1,
      58,    16,    58,    17,    -1,    58,    16,    58,    17,    35,
      58,    38,    -1,    58,    16,    58,    17,    35,    58,     5,
      58,    38,    -1,    42,    -1,    42,    35,   117,    38,    -1,
      42,     6,   117,    -1,    35,    83,    38,    -1,    35,    83,
      38,    35,   117,    38,    -1,    35,    83,    38,     6,   117,
      -1,    42,    34,   117,    37,    -1,    35,    83,    38,    34,
     117,    37,    -1,    60,    -1,    61,    -1,    62,    -1,    90,
      -1,    90,    16,    58,    17,    -1,    89,    -1,    92,     5,
      89,    -1,    88,    -1,    93,     5,    88,    -1,    91,    93,
      -1,    91,    35,    58,    38,    92,    -1,    42,    34,    86,
      37,    -1,    95,    -1,    96,     5,    95,    -1,    -1,    87,
      98,    85,    73,    -1,    43,    96,    73,    -1,    94,    73,
      -1,    59,    42,     6,   114,    73,    -1,    -1,    59,    87,
      42,    34,    82,    37,    99,    74,    81,    75,    -1,    -1,
      59,    42,    34,    82,    37,   100,     6,   116,    73,    -1,
      74,    -1,    75,    -1,    50,    -1,    51,    -1,    -1,    87,
     106,    85,    -1,    55,    -1,    84,    -1,    84,     5,    84,
      -1,    84,     5,    84,     5,    84,    -1,    73,    -1,    47,
      45,    -1,    48,    45,    -1,   107,    74,    81,    75,   111,
      -1,   114,    73,    -1,    97,    -1,    -1,   103,    35,   108,
      76,   124,    38,   110,   109,    -1,   103,    34,   114,    73,
     114,    73,   114,    37,   109,    -1,   103,    34,   105,    73,
     114,    73,   114,    37,   109,    -1,   104,    34,   114,    37,
     109,    -1,     3,    34,   114,    37,   109,    -1,     3,    34,
     114,    37,   109,     4,   109,    -1,   101,    81,   102,    -1,
      70,    42,   113,    -1,    70,    42,    35,   121,    38,    73,
      -1,    52,    73,    -1,    53,    73,    -1,    54,   114,    73,
      -1,    56,    34,    36,    36,    36,    37,   109,    -1,    34,
      42,     6,   114,     5,   114,    37,    -1,    34,    42,     6,
     114,     5,   114,    73,    42,    37,    -1,   112,   109,    -1,
     116,    -1,   114,     5,   114,    -1,    21,    -1,    20,    -1,
      27,    -1,    29,    -1,    28,    -1,   117,    -1,   117,     6,
     116,    -1,   117,    63,   116,    -1,   117,    64,   116,    -1,
     117,    65,   116,    -1,   117,    66,   116,    -1,   117,    67,
     116,    -1,   117,    68,   116,    -1,   118,    -1,   118,    77,
     118,    76,   118,    -1,   118,    76,   118,    -1,   118,    76,
     118,    76,   118,    -1,   122,    -1,   118,    22,   118,    -1,
     118,    26,   118,    -1,   118,    25,   118,    -1,   118,    23,
     118,    -1,   118,    24,   118,    -1,   118,    20,   118,    -1,
     118,    21,   118,    -1,   118,     9,   118,    -1,   118,     8,
     118,    -1,   118,    12,   118,    -1,   118,    13,   118,    -1,
     118,    10,   118,    -1,   118,    11,   118,    -1,   118,    16,
     118,    -1,   118,    19,   118,    -1,   118,    17,   118,    -1,
     118,    18,   118,    -1,   118,    15,   118,    -1,   118,    14,
     118,    -1,   118,    -1,    76,    -1,   118,    76,   118,    -1,
     118,    76,   118,    76,   118,    -1,    -1,    60,    -1,    61,
      -1,    62,    -1,    84,     6,   117,    -1,   119,    -1,   120,
       5,    60,    -1,   120,     5,    61,    -1,   120,     5,    62,
      -1,   120,     5,   119,    -1,   120,     5,    84,     6,   117,
      -1,   116,    -1,   121,     5,   116,    -1,   123,    -1,   115,
     123,    -1,   124,    -1,   124,    31,   122,    -1,   124,    33,
     122,    -1,   124,    32,    -1,    42,    -1,    39,    -1,    40,
      -1,    41,    -1,    45,    -1,   124,    34,   120,    37,    -1,
     124,    35,   119,    38,    -1,   124,    35,   119,     5,   119,
      38,    -1,   124,    35,    38,    -1,   124,    36,    42,    -1,
      60,    36,    42,    -1,    60,    34,   120,    37,    -1,    61,
      36,    42,    -1,    61,    34,   120,    37,    -1,    62,    36,
      42,    -1,    62,    34,   120,    37,    -1,   124,    29,    -1,
     124,    28,    -1,    58,    34,   114,    37,    -1,    34,   114,
      37,    -1,    35,   121,    38,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   290,   290,   344,   349,   350,   355,   356,   357,   358,
     359,   360,   361,   362,   363,   364,   365,   366,   367,   368,
     369,   370,   371,   372,   373,   376,   377,   382,   382,   382,
     382,   386,   387,   388,   390,   397,   398,   399,   400,   401,
     402,   403,   408,   409,   410,   411,   412,   413,   414,   415,
     423,   424,   425,   426,   427,   428,   431,   432,   438,   438,
     438,   439,   440,   445,   446,   448,   449,   451,   452,   458,
     461,   462,   467,   467,   468,   469,   470,   472,   471,   488,
     487,   496,   497,   499,   501,   506,   506,   509,   511,   512,
     513,   516,   517,   518,   519,   520,   521,   523,   522,   525,
     526,   530,   531,   532,   533,   535,   537,   540,   544,   548,
     555,   558,   565,   575,   583,   584,   590,   591,   592,   593,
     594,   600,   601,   602,   603,   604,   605,   606,   607,   613,
     614,   615,   616,   622,   623,   624,   625,   626,   627,   628,
     629,   630,   631,   632,   633,   634,   635,   636,   637,   638,
     639,   640,   641,   648,   649,   650,   651,   657,   658,   659,
     660,   663,   665,   666,   667,   668,   669,   672,   675,   676,
     681,   682,   687,   688,   689,   690,   700,   702,   703,   704,
     705,   710,   712,   713,   714,   715,   716,   717,   718,   719,
     720,   721,   722,   723,   724,   733,   734
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
  "MULEQ", "DIVEQ", "DOTMULEQ", "DOTDIVEQ", "ARROW", "BORDER", "CURVE",
  "SOLVE", "';'", "'{'", "'}'", "':'", "'?'", "$accept", "start", "input",
  "instructions", "list_of_id_args", "list_of_id1", "id", "list_of_dcls",
  "parameters_list", "type_of_dcl", "ID_space", "ID_array_space",
  "fespace123", "fespace", "spaceIDa", "spaceIDb", "spaceIDs",
  "fespace_def", "fespace_def_list", "declaration", "@1", "@2", "@3",
  "begin", "end", "for_loop", "while_loop", "declaration_for", "@4", "try",
  "IDfor", "instruction", "@5", "catchs", "bornes", "border_expr", "Expr",
  "unop", "no_comma_expr", "no_set_expr", "no_ternary_expr",
  "sub_script_expr", "parameters", "array", "unary_expr", "pow_expr",
  "primary", 0
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
     305,   306,   307,    59,   123,   125,    58,    63
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    78,    79,    80,    81,    81,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    82,    82,    82,    82,    82,
      82,    82,    82,    82,    82,    83,    83,    84,    84,    84,
      84,    85,    85,    85,    85,    86,    86,    86,    86,    86,
      86,    86,    87,    87,    87,    87,    87,    87,    87,    87,
      88,    88,    88,    88,    88,    88,    89,    89,    90,    90,
      90,    91,    91,    92,    92,    93,    93,    94,    94,    95,
      96,    96,    98,    97,    97,    97,    97,    99,    97,   100,
      97,   101,   102,   103,   104,   106,   105,   107,   108,   108,
     108,   109,   109,   109,   109,   109,   109,   110,   109,   109,
     109,   109,   109,   109,   109,   109,   109,   109,   109,   109,
     111,   112,   112,   113,   114,   114,   115,   115,   115,   115,
     115,   116,   116,   116,   116,   116,   116,   116,   116,   117,
     117,   117,   117,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   119,   119,   119,   119,   120,   120,   120,
     120,   120,   120,   120,   120,   120,   120,   120,   121,   121,
     122,   122,   123,   123,   123,   123,   124,   124,   124,   124,
     124,   124,   124,   124,   124,   124,   124,   124,   124,   124,
     124,   124,   124,   124,   124,   124,   124
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     1,     1,     2,     0,     1,     3,     2,
       3,     2,     3,     2,     3,     3,     3,     5,     5,     4,
       5,     4,     5,     4,     5,     1,     3,     1,     1,     1,
       1,     1,     3,     4,     3,     1,     2,     2,     2,     3,
       3,     5,     1,     4,     7,     6,     9,     4,     7,     9,
       1,     4,     3,     3,     6,     5,     4,     6,     1,     1,
       1,     1,     4,     1,     3,     1,     3,     2,     5,     4,
       1,     3,     0,     4,     3,     2,     5,     0,    10,     0,
       9,     1,     1,     1,     1,     0,     3,     1,     1,     3,
       5,     1,     2,     2,     5,     2,     1,     0,     8,     9,
       9,     5,     5,     7,     3,     3,     6,     2,     2,     3,
       7,     7,     9,     2,     1,     3,     1,     1,     1,     1,
       1,     1,     3,     3,     3,     3,     3,     3,     3,     1,
       5,     3,     5,     1,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     1,     1,     3,     5,     0,     1,     1,
       1,     3,     1,     3,     3,     3,     3,     5,     1,     3,
       1,     2,     1,     3,     3,     2,     1,     1,     1,     1,
       1,     4,     4,     6,     3,     3,     3,     4,     3,     4,
       3,     4,     2,     2,     4,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,   117,   116,   118,   120,   119,     0,     0,   177,
     178,   179,   176,     0,   180,     0,     0,    83,    84,     0,
       0,     0,    87,    42,     0,    58,    59,    60,     0,    91,
      81,     0,     0,     3,    72,    61,     0,     0,    96,     0,
       0,     0,     0,     4,     0,     0,   114,   121,   129,   133,
     170,   172,     0,     0,     0,     0,     0,     0,   168,     0,
       0,    70,     0,    92,    93,   107,   108,     0,     0,     0,
       0,     0,    42,     0,   157,     0,   157,     0,   157,     0,
       0,     1,     2,     5,     0,     0,     0,    50,    65,    67,
      75,     0,     0,     0,     0,     0,     0,    95,   171,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   193,   192,     0,
     175,     0,   157,     0,     0,     0,   195,     0,   196,     0,
       0,    74,   109,     0,     0,     0,     0,     6,     0,   176,
     158,   159,   160,   154,     0,   153,   162,     0,   186,     0,
     188,     0,   190,     0,     0,     0,   105,    31,     0,     0,
      27,     0,    28,    30,    29,     0,    25,     0,     0,     0,
      82,   104,    85,     0,     0,    88,     0,     0,     0,   115,
     122,   123,   124,   125,   126,   127,   128,   142,   141,   145,
     146,   143,   144,   152,   151,   147,   149,   150,   148,   139,
     140,   134,   137,   138,   136,   135,   131,     0,   173,   174,
       0,   184,     0,   185,     0,   169,   176,     0,     0,     0,
       0,    35,    71,    47,   194,     0,    43,     0,     6,    28,
      29,     0,     7,     0,     6,     0,     0,     0,   187,   189,
     191,     0,     0,   113,     0,     0,     0,    73,    62,     0,
       0,    53,    52,     0,     0,    66,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   181,     0,   182,   102,     0,
      36,    37,    38,     0,    69,     0,     0,     0,    76,     0,
       0,     9,     0,    11,     0,    79,     0,     0,    13,     0,
     161,   155,   163,   164,   165,     0,   166,     0,     0,    32,
       0,    34,     0,     0,    63,    68,    26,     0,     0,    51,
      86,     0,     0,    89,     0,   101,     0,    94,   132,   130,
       0,     0,    39,    28,    30,    29,     0,    40,     0,    45,
       0,    15,    10,    12,     6,    28,    29,    16,     0,     0,
       8,    14,    77,     0,     0,     0,   106,    33,     0,     0,
       0,    55,     0,     0,     0,     0,    97,     0,   183,   103,
       0,     0,    48,     0,    44,     0,     0,    19,     0,    21,
       0,     0,    23,     0,     0,   156,   167,     0,     0,     0,
      64,    54,     0,     0,    90,     0,     0,    41,     0,     0,
      17,    20,    22,    18,    24,     0,     0,   115,     0,    56,
       0,     0,    98,     0,    49,    46,    80,     0,   111,     0,
       0,   100,    99,     0,    78,     0,    57,     0,   112,   110
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    31,    32,    33,   241,   175,   154,   168,   230,    34,
      88,   314,    35,    36,   315,    89,    37,    61,    62,    38,
      84,   384,   349,    39,   181,    40,    41,   183,   266,    42,
     186,    43,   395,   327,   165,   166,    44,    45,    46,    47,
      48,   156,   157,    59,    49,    50,    51
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -244
static const yytype_int16 yypact[] =
{
     542,    10,  -244,  -244,  -244,  -244,  -244,   654,   654,  -244,
    -244,  -244,  -244,   -14,  -244,    -2,    82,  -244,  -244,    77,
      88,   654,  -244,   180,   -16,   120,   156,   259,   163,  -244,
    -244,   176,   164,   542,  -244,   202,     5,   107,  -244,   542,
      57,   205,   147,  -244,     1,   714,  -244,    81,   316,  -244,
    -244,   674,   654,   211,   120,   156,   259,    89,  -244,    13,
     222,  -244,     3,  -244,  -244,  -244,  -244,     4,   175,   654,
     196,   135,    67,   217,   766,   227,   766,   239,   766,   247,
     241,  -244,  -244,  -244,   252,   206,   225,    18,  -244,   292,
    -244,   362,   895,   -12,   654,   542,   654,  -244,  -244,   654,
     654,   654,   654,   654,   654,   654,   654,   654,   654,   654,
     654,   654,   654,   654,   654,   654,   654,   654,   654,   654,
     654,   654,   654,   654,   654,   654,   654,  -244,  -244,   654,
    -244,   654,   766,   723,   269,   131,  -244,   654,  -244,   931,
     -14,  -244,  -244,   299,   134,    14,   654,   242,   284,   337,
      32,   224,   229,  -244,   338,   609,  -244,   138,  -244,   152,
    -244,   160,  -244,   275,   654,   542,  -244,   145,     6,   331,
    -244,   313,  -244,  -244,  -244,    16,  -244,   654,   654,   100,
    -244,  -244,  -244,   276,     7,   347,   278,   165,   422,  -244,
    -244,  -244,  -244,  -244,  -244,  -244,  -244,  1039,  1039,  1054,
    1054,  1067,  1067,   796,   796,   585,   585,   585,   585,   283,
     283,  -244,  -244,  -244,  -244,  -244,   628,   647,  -244,  -244,
     190,  -244,    24,  -244,   542,  -244,   349,   198,   232,   254,
     199,  -244,  -244,   321,  -244,   301,   325,     9,   242,    -5,
     122,   204,   355,   125,   242,   654,   654,   809,  -244,  -244,
    -244,   357,    27,  -244,   654,   931,   252,  -244,  -244,   117,
     -12,    58,  -244,   326,   -12,  -244,   252,   654,   654,   -12,
     714,   542,   310,   654,   654,  -244,   852,  -244,   363,   654,
    -244,  -244,  -244,   967,  -244,   317,   330,   322,  -244,    29,
     -12,  -244,   -12,  -244,   311,  -244,   654,   -12,  -244,   214,
    -244,   716,    32,   224,   229,   368,  -244,   654,   312,  -244,
     215,  -244,   -12,   352,  -244,   382,  -244,   654,   654,  -244,
     383,    11,    12,   390,   285,  -244,   364,  -244,  1022,  1022,
     361,   542,  -244,   120,   156,   259,   394,  -244,    31,   371,
     370,  -244,  -244,  -244,   242,   146,   151,   405,   187,   412,
    -244,  -244,  -244,   654,   654,   414,  -244,  -244,    34,   654,
     117,  -244,   388,   654,   654,   -12,  -244,   391,  -244,  -244,
     654,   372,  -244,   373,  -244,    40,   -12,  -244,   -12,  -244,
     654,   -12,  -244,   654,   354,  1022,  -244,   654,   395,   396,
    -244,  -244,   218,   220,  -244,   542,   398,  -244,   400,   401,
    -244,  -244,  -244,  -244,  -244,   367,   542,   -32,   654,  -244,
     542,   542,  -244,   408,  -244,  -244,  -244,   482,  -244,   399,
     409,  -244,  -244,   410,  -244,   411,  -244,   542,  -244,  -244
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -244,  -244,  -244,   -36,  -213,   133,   -66,  -243,   197,   -22,
     274,    94,  -244,  -244,  -244,  -244,  -244,   315,  -244,  -244,
    -244,  -244,  -244,  -244,  -244,  -244,  -244,  -244,  -244,  -244,
    -244,   -33,  -244,  -244,  -244,  -244,    -6,  -244,    -4,  -117,
      -3,  -123,   -43,   294,   216,   415,   189
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -31
static const yytype_int16 yytable[] =
{
      83,    57,    73,    91,    58,   418,    96,   290,   140,    96,
     222,   256,    96,   311,    96,    67,    96,    96,   137,   235,
     176,   260,   231,   320,   177,   289,    71,   185,    60,   276,
     170,   299,   137,   159,   294,   161,   371,   170,   -28,   260,
      86,   419,    72,    63,    52,   294,   135,    87,   172,   173,
     174,   138,   236,   178,   261,   172,   173,   174,    83,   188,
     262,   263,   277,   144,   317,   308,    74,   341,    75,   372,
     182,   155,   388,   155,    97,   155,   141,   142,   400,   257,
     268,   242,   288,    68,   363,   364,   184,    99,   187,   220,
     189,    92,    93,   318,    96,   190,   191,   192,   193,   194,
     195,   196,    70,   197,   198,   199,   200,   201,   202,   203,
     204,   205,   206,   207,   208,   209,   210,   211,   212,   213,
     214,   215,   216,   217,   306,   243,   136,    64,   300,   155,
     155,   375,   253,   225,   292,   264,    96,   297,   231,    96,
     237,   146,    87,   247,   100,   101,   102,   103,   104,   105,
      65,   254,   312,   330,    74,    83,    75,   247,   376,   313,
      58,    66,   332,   378,   170,   247,   337,   170,   224,   147,
      96,   234,   242,   291,   293,   248,    81,   298,   242,   255,
      90,   305,   172,   173,   174,   172,   173,   174,   170,   249,
      76,   278,    77,   170,   316,   247,    68,   250,   176,   381,
     361,   362,   271,   323,   283,    80,   172,   173,   174,   294,
      82,   172,   173,   174,    69,    70,   243,   336,    85,   294,
     283,    95,   243,    96,   342,    96,   343,   275,   347,   170,
     -30,   351,    74,   143,    75,   -29,   284,   386,   325,    94,
     280,   295,   389,   301,   155,    69,   176,   172,   173,   174,
     309,   352,   357,   397,   145,   410,   139,   411,    76,   148,
      77,   321,   322,    78,   169,    79,    76,   170,    77,   158,
     328,   329,   348,   155,   281,   163,   164,   238,   242,   377,
     379,   160,   382,   171,   170,   172,   173,   174,    78,   162,
      79,   420,   350,    78,   167,    79,   282,   179,   369,   394,
      72,   355,   239,   173,   240,   120,   121,   122,   123,   124,
     401,   223,   402,   127,   128,   404,   233,   251,   244,   132,
     133,   134,   243,   366,   106,   107,   108,   109,   110,   111,
     112,   113,   114,   115,   116,   117,   118,   119,   120,   121,
     122,   123,   124,   -27,   245,   218,   344,   219,   258,   267,
     385,   259,   269,   170,   270,   279,   285,   392,   393,   286,
     287,   296,   412,   307,   319,     1,   326,   331,   339,    72,
     417,   345,   173,   346,   354,   338,   403,   421,   422,   405,
     340,   407,     2,     3,    83,   356,   359,   360,   256,     4,
       5,     6,   125,   126,   429,   365,     7,     8,   367,   368,
     370,     9,    10,    11,    12,    13,   373,    14,   374,    15,
      16,   380,    17,    18,    19,    20,    21,    22,   383,   387,
      23,    24,    25,    26,    27,     1,   391,   396,   406,   408,
     398,   399,    28,   409,   413,    29,    30,   180,   414,   415,
     416,   425,     2,     3,   423,   358,   426,   427,   428,     4,
       5,     6,   310,   265,   390,   232,     7,     8,   252,   324,
      98,     9,    10,    11,    12,    13,     0,    14,     0,    15,
      16,     0,    17,    18,    19,    20,    21,    22,     0,     0,
      23,    24,    25,    26,    27,     1,     0,     0,     0,     0,
       0,     0,    28,     0,     0,    29,    30,   272,     0,     0,
       0,     0,     2,     3,     0,     0,     0,     0,     0,     4,
       5,     6,     0,     0,     0,     0,     7,     8,     0,     0,
       0,     9,    10,    11,    12,    13,     0,    14,     0,    15,
      16,     0,    17,    18,    19,    20,    21,    22,     0,     0,
      23,    24,    25,    26,    27,     1,     0,     0,     0,     0,
       0,     0,    28,     0,     0,    29,    30,   424,     0,     0,
       0,     0,     2,     3,     0,     0,     0,     0,     0,     4,
       5,     6,     0,     0,     0,     0,     7,     8,     0,     0,
       0,     9,    10,    11,    12,    13,     0,    14,     0,    15,
      16,     0,    17,    18,    19,    20,    21,    22,     0,     0,
      23,    24,    25,    26,    27,   118,   119,   120,   121,   122,
     123,   124,    28,     0,     0,    29,    30,   106,   107,   108,
     109,   110,   111,   112,   113,   114,   115,   116,   117,   118,
     119,   120,   121,   122,   123,   124,   106,   107,   108,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     120,   121,   122,   123,   124,   106,   107,   108,   109,   110,
     111,   112,   113,   114,   115,   116,   117,   118,   119,   120,
     121,   122,   123,   124,     2,     3,     0,     0,     0,     0,
       0,     4,     5,     6,     0,   246,     0,     0,     7,     8,
       0,     0,     0,     9,    10,    11,    12,     0,     0,    14,
       0,     0,   127,   128,   273,   129,   130,   131,   132,   133,
     134,     0,    53,     0,    54,    55,    56,     0,     0,     0,
       0,     0,     0,   274,   106,   107,   108,   109,   110,   111,
     112,   113,   114,   115,   116,   117,   118,   119,   120,   121,
     122,   123,   124,     2,     3,     0,     0,     0,     7,     8,
       4,     5,     6,     9,    10,    11,    12,     7,     8,    14,
       0,   221,     9,    10,    11,    12,     0,     0,    14,     0,
       0,     0,    53,     0,    54,    55,    56,     0,     0,     0,
       0,    53,     0,    54,    55,    56,     2,     3,     0,     0,
       0,     0,   353,     4,     5,     6,     0,     0,     0,   153,
       7,     8,     0,     0,     0,     9,    10,    11,   149,     0,
       0,    14,   114,   115,   116,   117,   118,   119,   120,   121,
     122,   123,   124,     0,    53,     0,   150,   151,   152,     2,
       3,     0,     0,     0,     0,     0,     4,     5,     6,     0,
       0,     0,   153,     7,     8,     0,     0,     0,     9,    10,
      11,   149,     0,     0,    14,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    53,     0,   302,
     303,   304,     2,     3,     0,     0,     0,     0,     0,     4,
       5,     6,     0,     0,     0,   153,     7,     8,     0,     0,
       0,     9,    10,    11,    12,     0,     0,    14,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      53,     0,    54,    55,    56,     2,     3,     0,     0,     0,
       0,     0,     4,     5,     6,     0,     0,     0,   153,     7,
       8,     0,     0,     0,     9,    10,    11,    12,     0,     0,
      14,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     2,     3,    23,     0,    54,    55,    56,     4,     5,
       6,     0,     0,     0,     0,     7,     8,     0,     0,     0,
       9,    10,    11,   226,     0,     0,    14,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     2,     3,    53,
       0,   227,   228,   229,     4,     5,     6,     0,     0,     0,
       0,     7,     8,     0,     0,     0,     9,    10,    11,   149,
       0,     0,    14,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    53,     0,   333,   334,   335,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,   119,   120,   121,   122,   123,   124,   108,
     109,   110,   111,   112,   113,   114,   115,   116,   117,   118,
     119,   120,   121,   122,   123,   124,   110,   111,   112,   113,
     114,   115,   116,   117,   118,   119,   120,   121,   122,   123,
     124,   112,   113,   114,   115,   116,   117,   118,   119,   120,
     121,   122,   123,   124
};

static const yytype_int16 yycheck[] =
{
      33,     7,    24,    39,     8,    37,     5,    12,     5,     5,
     133,     5,     5,   256,     5,    21,     5,     5,     5,     5,
      86,     5,   139,   266,     6,   238,    42,    93,    42,     5,
      42,   244,     5,    76,     5,    78,     5,    42,     6,     5,
      35,    73,    58,    45,    34,     5,    52,    42,    60,    61,
      62,    38,    38,    35,    38,    60,    61,    62,    91,    95,
     177,   178,    38,    69,     6,    38,    34,    38,    36,    38,
      92,    74,    38,    76,    73,    78,    73,    73,    38,    73,
      73,   147,    73,    16,    73,    73,    92,     6,    94,   132,
      96,    34,    35,    35,     5,    99,   100,   101,   102,   103,
     104,   105,    35,   106,   107,   108,   109,   110,   111,   112,
     113,   114,   115,   116,   117,   118,   119,   120,   121,   122,
     123,   124,   125,   126,   247,   147,    37,    45,   245,   132,
     133,   344,   165,   137,    12,    35,     5,    12,   255,     5,
     146,     6,    42,     5,    63,    64,    65,    66,    67,    68,
      73,     6,    35,   276,    34,   188,    36,     5,    12,    42,
     164,    73,   279,    12,    42,     5,   283,    42,    37,    34,
       5,    37,   238,   239,   240,    37,     0,   243,   244,    34,
      73,   247,    60,    61,    62,    60,    61,    62,    42,    37,
      34,   224,    36,    42,   260,     5,    16,    37,   264,    12,
     317,   318,    37,   269,     5,    42,    60,    61,    62,     5,
      46,    60,    61,    62,    34,    35,   238,   283,    16,     5,
       5,    74,   244,     5,   290,     5,   292,    37,   294,    42,
       6,   297,    34,    58,    36,     6,    37,   354,   271,    34,
      42,    37,   359,   246,   247,    34,   312,    60,    61,    62,
     254,    37,    37,   370,    58,    37,    34,    37,    34,    42,
      36,   267,   268,    34,    58,    36,    34,    42,    36,    42,
     273,   274,   294,   276,    42,    34,    35,    35,   344,   345,
     346,    42,   348,    58,    42,    60,    61,    62,    34,    42,
      36,   408,   296,    34,    42,    36,    42,     5,   331,   365,
      58,   307,    60,    61,    62,    22,    23,    24,    25,    26,
     376,    42,   378,    28,    29,   381,    17,    42,    34,    34,
      35,    36,   344,    38,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,     6,     6,   129,    35,   131,    17,    73,
     353,    38,     5,    42,    76,     6,    35,   363,   364,    58,
      35,     6,   395,     6,    38,     3,    56,     4,    38,    58,
     406,    60,    61,    62,     6,    58,   380,   410,   411,   383,
      58,   387,    20,    21,   417,    73,    34,     5,     5,    27,
      28,    29,    76,    77,   427,     5,    34,    35,    34,    38,
       6,    39,    40,    41,    42,    43,    35,    45,    38,    47,
      48,     6,    50,    51,    52,    53,    54,    55,     6,     5,
      58,    59,    60,    61,    62,     3,    38,    36,    74,    34,
      58,    58,    70,    37,    36,    73,    74,    75,    38,    38,
      73,    42,    20,    21,    36,   312,    37,    37,    37,    27,
      28,    29,   255,   179,   360,   140,    34,    35,   164,   270,
      45,    39,    40,    41,    42,    43,    -1,    45,    -1,    47,
      48,    -1,    50,    51,    52,    53,    54,    55,    -1,    -1,
      58,    59,    60,    61,    62,     3,    -1,    -1,    -1,    -1,
      -1,    -1,    70,    -1,    -1,    73,    74,    75,    -1,    -1,
      -1,    -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    27,
      28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,
      -1,    39,    40,    41,    42,    43,    -1,    45,    -1,    47,
      48,    -1,    50,    51,    52,    53,    54,    55,    -1,    -1,
      58,    59,    60,    61,    62,     3,    -1,    -1,    -1,    -1,
      -1,    -1,    70,    -1,    -1,    73,    74,    75,    -1,    -1,
      -1,    -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    27,
      28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,
      -1,    39,    40,    41,    42,    43,    -1,    45,    -1,    47,
      48,    -1,    50,    51,    52,    53,    54,    55,    -1,    -1,
      58,    59,    60,    61,    62,    20,    21,    22,    23,    24,
      25,    26,    70,    -1,    -1,    73,    74,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    20,    21,    -1,    -1,    -1,    -1,
      -1,    27,    28,    29,    -1,    76,    -1,    -1,    34,    35,
      -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,
      -1,    -1,    28,    29,    76,    31,    32,    33,    34,    35,
      36,    -1,    58,    -1,    60,    61,    62,    -1,    -1,    -1,
      -1,    -1,    -1,    76,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    20,    21,    -1,    -1,    -1,    34,    35,
      27,    28,    29,    39,    40,    41,    42,    34,    35,    45,
      -1,    38,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    58,    -1,    60,    61,    62,    -1,    -1,    -1,
      -1,    58,    -1,    60,    61,    62,    20,    21,    -1,    -1,
      -1,    -1,    76,    27,    28,    29,    -1,    -1,    -1,    76,
      34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,
      -1,    45,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    -1,    58,    -1,    60,    61,    62,    20,
      21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,    -1,
      -1,    -1,    76,    34,    35,    -1,    -1,    -1,    39,    40,
      41,    42,    -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    58,    -1,    60,
      61,    62,    20,    21,    -1,    -1,    -1,    -1,    -1,    27,
      28,    29,    -1,    -1,    -1,    76,    34,    35,    -1,    -1,
      -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      58,    -1,    60,    61,    62,    20,    21,    -1,    -1,    -1,
      -1,    -1,    27,    28,    29,    -1,    -1,    -1,    76,    34,
      35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    20,    21,    58,    -1,    60,    61,    62,    27,    28,
      29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    20,    21,    58,
      -1,    60,    61,    62,    27,    28,    29,    -1,    -1,    -1,
      -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      -1,    -1,    45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    58,    -1,    60,    61,    62,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    20,    21,    27,    28,    29,    34,    35,    39,
      40,    41,    42,    43,    45,    47,    48,    50,    51,    52,
      53,    54,    55,    58,    59,    60,    61,    62,    70,    73,
      74,    79,    80,    81,    87,    90,    91,    94,    97,   101,
     103,   104,   107,   109,   114,   115,   116,   117,   118,   122,
     123,   124,    34,    58,    60,    61,    62,   114,   116,   121,
      42,    95,    96,    45,    45,    73,    73,   114,    16,    34,
      35,    42,    58,    87,    34,    36,    34,    36,    34,    36,
      42,     0,    46,   109,    98,    16,    35,    42,    88,    93,
      73,    81,    34,    35,    34,    74,     5,    73,   123,     6,
      63,    64,    65,    66,    67,    68,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    76,    77,    28,    29,    31,
      32,    33,    34,    35,    36,   114,    37,     5,    38,    34,
       5,    73,    73,    58,   114,    58,     6,    34,    42,    42,
      60,    61,    62,    76,    84,   118,   119,   120,    42,   120,
      42,   120,    42,    34,    35,   112,   113,    42,    85,    58,
      42,    58,    60,    61,    62,    83,    84,     6,    35,     5,
      75,   102,    87,   105,   114,    84,   108,   114,    81,   114,
     116,   116,   116,   116,   116,   116,   116,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   118,   118,   118,   118,
     118,   118,   118,   118,   118,   118,   118,   118,   122,   122,
     120,    38,   119,    42,    37,   116,    42,    60,    61,    62,
      86,   117,    95,    17,    37,     5,    38,   114,    35,    60,
      62,    82,    84,    87,    34,     6,    76,     5,    37,    37,
      37,    42,   121,   109,     6,    34,     5,    73,    17,    38,
       5,    38,   117,   117,    35,    88,   106,    73,    73,     5,
      76,    37,    75,    76,    76,    37,     5,    38,   109,     6,
      42,    42,    42,     5,    37,    35,    58,    35,    73,    82,
      12,    84,    12,    84,     5,    37,     6,    12,    84,    82,
     117,   118,    60,    61,    62,    84,   119,     6,    38,   116,
      86,    85,    35,    42,    89,    92,    84,     6,    35,    38,
      85,   114,   114,    84,   124,   109,    56,   111,   118,   118,
     119,     4,   117,    60,    61,    62,    84,   117,    58,    38,
      58,    38,    84,    84,    35,    60,    62,    84,    87,   100,
     116,    84,    37,    76,     6,   114,    73,    37,    83,    34,
       5,   117,   117,    73,    73,     5,    38,    34,    38,   109,
       6,     5,    38,    35,    38,    82,    12,    84,    12,    84,
       6,    12,    84,     6,    99,   118,   117,     5,    38,   117,
      89,    38,   114,   114,    84,   110,    36,   117,    58,    58,
      38,    84,    84,   116,    84,   116,    74,   114,    34,    37,
      37,    37,   109,    36,    38,    38,    73,    81,    37,    73,
     117,   109,   109,    36,    75,    42,    37,    37,    37,   109
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
#line 290 "lg.ypp"
    {
                        if(  ffapi::ff_justcompile) exit(0);
    // clean FH  mach 2014
		        const char *  magicffglut="#!ffglutdata3.1\n";// for complex and vector 3d plot
			//FFCS: divert stream to FFCS
                        if(ThePlotStream) ffapi::fwriteinit(magicffglut,strlen(magicffglut),1,ThePlotStream);
                        size_t sizestack = currentblock->size()+1024 ; //  before close
                        (yyvsp[(1) - (2)].cinst)+=currentblock->close(currentblock);
                        if(verbosity>2 || ( (mpirank==0) && verbosity)) cout << " sizestack + 1024 =" << sizestack << "  ( " << sizestack-1024 <<" )\n" ;
                        size_t lg0,lg1;
                        ffapi::ifchtmpdir(); // change  to tmp after compile FH sep 2015 ...
                        int NbPtr = ShowAlloc("init execution ",lg0); // number of un delele ptr
                        debugstack= new vector<pair<const E_Routine*,int> >;
                        size_t stu0=storageused(); // get Storage usage
			UnShowAlloc =0;// add FH for parallee
                        if(verbosity>2  || ( (mpirank==0) && verbosity) ) cout << endl;
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
			//FFCS: divert stream to FFCS
			if(ThePlotStream) {ffapi::ff_pclose(ThePlotStream); ThePlotStream=0;}
			UnShowAlloc =1;
                        if(debugstack) delete debugstack;
                        NbPtr = ShowAlloc("end execution -- ",lg1) - NbPtr;
                        long stu1 =storageused()-stu0    ;


			    if (verbosity && (NbPtr || (stu1>100000) )) { cout << " ######## We forget of deleting   " << NbPtr
			                      << " Nb pointer,   " <<  lg1-lg0 << "Bytes " << " ,  mpirank " << mpirank << ", memory leak ="<< stu1 <<  endl;}
  return 0;;}
    break;

  case 4:
#line 349 "lg.ypp"
    {(yyval.cinst)=(yyvsp[(1) - (1)].cexp);;;;}
    break;

  case 5:
#line 350 "lg.ypp"
    { (yyval.cinst)= ((yyvsp[(1) - (2)].cinst)+=(yyvsp[(2) - (2)].cexp)) ;}
    break;

  case 6:
#line 355 "lg.ypp"
    { (yyval.clist_id)=new ListOfId();;}
    break;

  case 7:
#line 356 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str)));}
    break;

  case 8:
#line 357 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (3)].str),(yyvsp[(3) - (3)].cexp))) ;}
    break;

  case 9:
#line 358 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,2> **>()));}
    break;

  case 10:
#line 359 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,2> **>(),true));}
    break;

  case 11:
#line 360 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,3> **>()));}
    break;

  case 12:
#line 361 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,3> **>(),true));}
    break;

  case 13:
#line 362 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),C_F0(),(yyvsp[(1) - (2)].type)->right())) ;}
    break;

  case 14:
#line 363 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),C_F0(),(yyvsp[(1) - (3)].type),true)) ;}
    break;

  case 15:
#line 364 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (3)].clist_id))) ;}
    break;

  case 16:
#line 365 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (3)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str))) ;}
    break;

  case 17:
#line 366 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (5)].clist_id))) ;}
    break;

  case 18:
#line 367 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (5)].str),(yyvsp[(5) - (5)].cexp))) ;}
    break;

  case 19:
#line 368 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,2> **>())) ;}
    break;

  case 20:
#line 369 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,2> **>(),true)) ;}
    break;

  case 21:
#line 370 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,3> **>())) ;}
    break;

  case 22:
#line 371 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,3> **>(),true)) ;}
    break;

  case 23:
#line 372 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),C_F0(),(yyvsp[(3) - (4)].type)->right())) ;}
    break;

  case 24:
#line 373 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),C_F0(),(yyvsp[(3) - (5)].type),true)) ;}
    break;

  case 25:
#line 376 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str))); ;}
    break;

  case 26:
#line 377 "lg.ypp"
    { (yyval.clist_id)=(yyvsp[(1) - (3)].clist_id)  ; (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str))); ;}
    break;

  case 31:
#line 386 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (1)].str),dcltype);}
    break;

  case 32:
#line 387 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (3)].str),dcltype,(yyvsp[(3) - (3)].cexp));}
    break;

  case 33:
#line 388 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (4)].str),dcltype,(yyvsp[(3) - (4)].args));
                                              (yyvsp[(3) - (4)].args).destroy();}
    break;

  case 34:
#line 390 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 35:
#line 397 "lg.ypp"
    {(yyval.args)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 36:
#line 398 "lg.ypp"
    {(yyval.args)=Find((yyvsp[(1) - (2)].str));}
    break;

  case 37:
#line 399 "lg.ypp"
    {(yyval.args)=Find((yyvsp[(1) - (2)].str));}
    break;

  case 38:
#line 400 "lg.ypp"
    {(yyval.args)=Find((yyvsp[(1) - (2)].str));}
    break;

  case 39:
#line 401 "lg.ypp"
    { (yyval.args)=make_pair<const char *,const C_F0>((yyvsp[(1) - (3)].str),(yyvsp[(3) - (3)].cexp));}
    break;

  case 40:
#line 402 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp)) ;}
    break;

  case 41:
#line 403 "lg.ypp"
    { (yyval.args)= ((yyvsp[(1) - (5)].args)+= make_pair<const char *,const C_F0>((yyvsp[(3) - (5)].str),(yyvsp[(5) - (5)].cexp)));}
    break;

  case 43:
#line 409 "lg.ypp"
    {(yyval.type)=TypeArray((yyvsp[(1) - (4)].type),(yyvsp[(3) - (4)].type));}
    break;

  case 44:
#line 410 "lg.ypp"
    {(yyval.type)=TypeArray(TypeArray((yyvsp[(1) - (7)].type),(yyvsp[(3) - (7)].type)),(yyvsp[(6) - (7)].type));}
    break;

  case 45:
#line 411 "lg.ypp"
    {(yyval.type)=TypeArray((yyvsp[(1) - (6)].type),(yyvsp[(3) - (6)].type),(yyvsp[(5) - (6)].type));}
    break;

  case 46:
#line 412 "lg.ypp"
    {(yyval.type)=TypeArray(TypeArray((yyvsp[(1) - (9)].type),(yyvsp[(3) - (9)].type),(yyvsp[(5) - (9)].type)),(yyvsp[(8) - (9)].type));}
    break;

  case 47:
#line 413 "lg.ypp"
    {(yyval.type)=TypeTemplate((yyvsp[(1) - (4)].type),(yyvsp[(3) - (4)].type));}
    break;

  case 48:
#line 414 "lg.ypp"
    {(yyval.type)=TypeArray(TypeTemplate((yyvsp[(1) - (7)].type),(yyvsp[(3) - (7)].type)),(yyvsp[(6) - (7)].type));}
    break;

  case 49:
#line 415 "lg.ypp"
    {(yyval.type)=TypeArray(TypeTemplate((yyvsp[(1) - (9)].type),(yyvsp[(3) - (9)].type)),(yyvsp[(6) - (9)].type),(yyvsp[(8) - (9)].type));}
    break;

  case 50:
#line 423 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(1) - (1)].str),currentblock,fespacetype,fespacecomplex,fespacedim); ;}
    break;

  case 51:
#line 424 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(1) - (4)].str),currentblock,fespacetype,(yyvsp[(3) - (4)].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 52:
#line 425 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(1) - (3)].str),currentblock,fespacetype,(yyvsp[(3) - (3)].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 53:
#line 426 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(2) - (3)].clist_id),currentblock,fespacetype,fespacecomplex,fespacedim) ;}
    break;

  case 54:
#line 427 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(2) - (6)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (6)].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 55:
#line 428 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(2) - (5)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (5)].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 56:
#line 431 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(1) - (4)].str),currentblock,fespacetype,(yyvsp[(3) - (4)].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 57:
#line 432 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(2) - (6)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (6)].cexp),fespacecomplex,fespacedim) ;}
    break;

  case 58:
#line 438 "lg.ypp"
    { fespacedim=2;}
    break;

  case 59:
#line 438 "lg.ypp"
    { fespacedim=1;}
    break;

  case 60:
#line 438 "lg.ypp"
    { fespacedim=3;}
    break;

  case 61:
#line 439 "lg.ypp"
    {fespacecomplex=false;  fespacetype = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 62:
#line 440 "lg.ypp"
    {
             if ((yyvsp[(3) - (4)].type) != typevarreal && (yyvsp[(3) - (4)].type) != typevarcomplex) lgerror (" type of finite element <real> or <complex>");
             fespacecomplex=((yyvsp[(3) - (4)].type)==typevarcomplex);
             fespacetype = Find((yyvsp[(1) - (4)].str));;}
    break;

  case 63:
#line 445 "lg.ypp"
    {  (yyval.cexp) = (yyvsp[(1) - (1)].cexp)  ;}
    break;

  case 64:
#line 446 "lg.ypp"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 65:
#line 448 "lg.ypp"
    {  (yyval.cexp) = (yyvsp[(1) - (1)].cexp)  ;}
    break;

  case 66:
#line 449 "lg.ypp"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 67:
#line 451 "lg.ypp"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[(2) - (2)].cexp);}
    break;

  case 68:
#line 452 "lg.ypp"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[(5) - (5)].cexp);}
    break;

  case 69:
#line 458 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,size_t>((yyvsp[(1) - (4)].str),typeFESpace((yyvsp[(3) - (4)].args)),(yyvsp[(3) - (4)].args),dimFESpaceImage((yyvsp[(3) - (4)].args)));
     (yyvsp[(3) - (4)].args).destroy(); ;}
    break;

  case 71:
#line 462 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 72:
#line 467 "lg.ypp"
    {dcltype=(yyvsp[(1) - (1)].type);}
    break;

  case 73:
#line 467 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(3) - (4)].cexp);}
    break;

  case 74:
#line 468 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(2) - (3)].cexp);}
    break;

  case 75:
#line 469 "lg.ypp"
    { (yyval.cexp)=(yyvsp[(1) - (2)].cexp);}
    break;

  case 76:
#line 470 "lg.ypp"
    {(yyval.cexp)=currentblock->NewID((yyvsp[(1) - (5)].type),(yyvsp[(2) - (5)].str),(yyvsp[(4) - (5)].cexp));;}
    break;

  case 77:
#line 472 "lg.ypp"
    {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = (yyvsp[(2) - (6)].type)->right();
                      routineinblock[kkembtype] = currentblock;
                      (yyvsp[(5) - (6)].routine)=new Routine((yyvsp[(1) - (6)].type),(yyvsp[(2) - (6)].type)->right(),(yyvsp[(3) - (6)].str),(yyvsp[(5) - (6)].clist_id),currentblock);
		      // routineinblock[kkembtype]->Add($3,"(",$<routine>5); //pas recursif pour l'instanat test  FH 27 dec 2008
                     // cout << " \n after new routine \n " << endl;
                      ;}
    break;

  case 78:
#line 481 "lg.ypp"
    { currentblock=(yyvsp[(5) - (10)].routine)->Set((yyvsp[(9) - (10)].cinst));
                       currentblock->Add((yyvsp[(3) - (10)].str),"(",(yyvsp[(5) - (10)].routine)); //pas recursif pour l'instant test  FH 27 dec 2008
                       kkembtype--;
                       (yyval.cexp)=0;

                        ;}
    break;

  case 79:
#line 488 "lg.ypp"
    {Block::open(currentblock); (yyvsp[(1) - (5)].type)->SetArgs((yyvsp[(4) - (5)].clist_id));;}
    break;

  case 80:
#line 490 "lg.ypp"
    {  (yyval.cinst)=currentblock->close(currentblock);
                         (yyval.cexp)=currentblock->NewID((yyvsp[(1) - (9)].type),(yyvsp[(2) - (9)].str),(yyvsp[(8) - (9)].cexp),*(yyvsp[(4) - (9)].clist_id));
                         delete (yyvsp[(4) - (9)].clist_id); //  FH 23032005
                         ;}
    break;

  case 81:
#line 496 "lg.ypp"
    {  Block::open(currentblock);}
    break;

  case 82:
#line 497 "lg.ypp"
    {  (yyval.cexp)=currentblock->close(currentblock);}
    break;

  case 83:
#line 499 "lg.ypp"
    {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 84:
#line 501 "lg.ypp"
    {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 85:
#line 506 "lg.ypp"
    {dcltype=(yyvsp[(1) - (1)].type); Block::open(currentblock);  ;}
    break;

  case 86:
#line 507 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(3) - (3)].cexp);}
    break;

  case 87:
#line 509 "lg.ypp"
    { Block::open(currentblock) ;}
    break;

  case 88:
#line 511 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str)));Block::open(currentblock); ;}
    break;

  case 89:
#line 512 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (3)].str)));(yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str)));Block::open(currentblock); ;}
    break;

  case 90:
#line 513 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (5)].str)));(yyval.clist_id)->push_back(UnId((yyvsp[(3) - (5)].str)));(yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str)));Block::open(currentblock); ;}
    break;

  case 91:
#line 516 "lg.ypp"
    {(yyval.cexp)=0;;}
    break;

  case 92:
#line 517 "lg.ypp"
    {zzzfff->input((yyvsp[(2) - (2)].str));(yyval.cexp)= 0; ;}
    break;

  case 93:
#line 518 "lg.ypp"
    {load((yyvsp[(2) - (2)].str));(yyval.cexp)= 0; ;}
    break;

  case 94:
#line 519 "lg.ypp"
    {(yyval.cexp)=Try((yyvsp[(3) - (5)].cinst),(yyvsp[(5) - (5)].cexp),currentblock->close(currentblock));;}
    break;

  case 95:
#line 520 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(1) - (2)].cexp);}
    break;

  case 96:
#line 521 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 97:
#line 523 "lg.ypp"
    {(yyvsp[(5) - (6)].cexp)=ForAll(currentblock,(yyvsp[(3) - (6)].clist_id),(yyvsp[(5) - (6)].cexp));;}
    break;

  case 98:
#line 524 "lg.ypp"
    {inloopcount--; (yyval.cexp)=ForAll((yyvsp[(5) - (8)].cexp),(yyvsp[(8) - (8)].cexp),currentblock->close(currentblock));}
    break;

  case 99:
#line 525 "lg.ypp"
    {inloopcount--; (yyval.cexp)=For((yyvsp[(3) - (9)].cexp),(yyvsp[(5) - (9)].cexp),(yyvsp[(7) - (9)].cexp),(yyvsp[(9) - (9)].cexp));}
    break;

  case 100:
#line 527 "lg.ypp"
    {inloopcount--;
                (yyval.cexp)=C_F0(For((yyvsp[(3) - (9)].cexp),(yyvsp[(5) - (9)].cexp),(yyvsp[(7) - (9)].cexp),(yyvsp[(9) - (9)].cexp)),currentblock->close(currentblock));}
    break;

  case 101:
#line 530 "lg.ypp"
    {inloopcount--;(yyval.cexp)=While((yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 102:
#line 531 "lg.ypp"
    {(yyval.cexp)=FIf((yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 103:
#line 532 "lg.ypp"
    {(yyval.cexp)=FIf((yyvsp[(3) - (7)].cexp),(yyvsp[(5) - (7)].cexp),(yyvsp[(7) - (7)].cexp));}
    break;

  case 104:
#line 533 "lg.ypp"
    {
                      (yyval.cexp)=C_F0(new E_block((yyvsp[(2) - (3)].cinst),(yyvsp[(3) - (3)].cexp)),atype<void>()) ;}
    break;

  case 105:
#line 535 "lg.ypp"
    { /* <<BORDER_ID>> */
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[(2) - (3)].str),C_F0(TheOperators,"[border]",(yyvsp[(3) - (3)].args)));}
    break;

  case 106:
#line 537 "lg.ypp"
    {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[(2) - (6)].str),C_F0(TheOperators,"[border]",(yyvsp[(4) - (6)].args)));}
    break;

  case 107:
#line 540 "lg.ypp"
    {
                    if(inloopcount)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_break),atype<void>());
                    else lgerror("break not in loop") ;}
    break;

  case 108:
#line 544 "lg.ypp"
    {
                    if(inloopcount)
                        (yyval.cexp)= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");}
    break;

  case 109:
#line 548 "lg.ypp"
    {
                    if (kkembtype>=0)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_return,(rettype[kkembtype]->CastTo((yyvsp[(2) - (3)].cexp))).OnReturn()) ,atype<void>());
                     else lgerror(" return not in routine ") ;}
    break;

  case 110:
#line 555 "lg.ypp"
    {(yyval.cexp) =  (yyvsp[(7) - (7)].cexp); ;}
    break;

  case 111:
#line 558 "lg.ypp"
    {
   Block::open(currentblock);
   (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[(2) - (7)].str),atype<double*>());
   (yyval.args)+= (yyvsp[(4) - (7)].cexp);
   (yyval.args)+= (yyvsp[(6) - (7)].cexp);
   (yyval.args)+= currentblock->NewVar<LocalVariable>("IndexBorder",atype<long*>());;}
    break;

  case 112:
#line 565 "lg.ypp"
    {
    Block::open(currentblock);
    (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[(2) - (9)].str),atype<double*>());
    (yyval.args)+= (yyvsp[(4) - (9)].cexp);
    (yyval.args)+= (yyvsp[(6) - (9)].cexp);
    (yyval.args)+= currentblock->NewVar<LocalVariable>((yyvsp[(8) - (9)].str),atype<long*>());;}
    break;

  case 113:
#line 575 "lg.ypp"
    {
   (yyval.args) = ((yyvsp[(1) - (2)].args) += (yyvsp[(2) - (2)].cexp));
   currentblock->close(currentblock);}
    break;

  case 115:
#line 584 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 122:
#line 601 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 123:
#line 602 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"+=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 124:
#line 603 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"-=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 125:
#line 604 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"*=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 126:
#line 605 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"/=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 127:
#line 606 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,".*=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 128:
#line 607 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"./=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 130:
#line 614 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"?:",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 131:
#line 615 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 132:
#line 616 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 134:
#line 623 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 135:
#line 624 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 136:
#line 625 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 137:
#line 626 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 138:
#line 627 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 139:
#line 628 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 140:
#line 629 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 141:
#line 630 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 142:
#line 631 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 143:
#line 632 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 144:
#line 633 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 145:
#line 634 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 146:
#line 635 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 147:
#line 636 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 148:
#line 637 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 149:
#line 638 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 150:
#line 639 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 151:
#line 640 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 152:
#line 641 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 153:
#line 648 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 154:
#line 649 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,":");}
    break;

  case 155:
#line 650 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 156:
#line 651 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));}
    break;

  case 157:
#line 657 "lg.ypp"
    {(yyval.args)=0;}
    break;

  case 158:
#line 658 "lg.ypp"
    {(yyval.args)=Find((yyvsp[(1) - (1)].str));}
    break;

  case 159:
#line 659 "lg.ypp"
    {(yyval.args)=Find((yyvsp[(1) - (1)].str));}
    break;

  case 160:
#line 660 "lg.ypp"
    {(yyval.args)=Find((yyvsp[(1) - (1)].str));}
    break;

  case 161:
#line 663 "lg.ypp"
    { (yyval.args)=make_pair<const char *,const C_F0>((yyvsp[(1) - (3)].str),(yyvsp[(3) - (3)].cexp));}
    break;

  case 162:
#line 665 "lg.ypp"
    {(yyval.args)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 163:
#line 666 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str))) ;}
    break;

  case 164:
#line 667 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str))) ;}
    break;

  case 165:
#line 668 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str))) ;}
    break;

  case 166:
#line 669 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp)) ;}
    break;

  case 167:
#line 672 "lg.ypp"
    { (yyval.args)= ((yyvsp[(1) - (5)].args)+= make_pair<const char *,const C_F0>((yyvsp[(3) - (5)].str),(yyvsp[(5) - (5)].cexp))) ;}
    break;

  case 168:
#line 675 "lg.ypp"
    {(yyval.args)=(yyvsp[(1) - (1)].cexp);}
    break;

  case 169:
#line 676 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp)) ;}
    break;

  case 171:
#line 682 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(1) - (2)].oper),(yyvsp[(2) - (2)].cexp));}
    break;

  case 173:
#line 688 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 174:
#line 689 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));}
    break;

  case 175:
#line 690 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));}
    break;

  case 176:
#line 700 "lg.ypp"
    {(yyval.cexp)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 177:
#line 702 "lg.ypp"
    {(yyval.cexp)= CConstant((yyvsp[(1) - (1)].lnum));}
    break;

  case 178:
#line 703 "lg.ypp"
    {(yyval.cexp)= CConstant((yyvsp[(1) - (1)].dnum));}
    break;

  case 179:
#line 704 "lg.ypp"
    {(yyval.cexp)= CConstant(complex<double>(0,(yyvsp[(1) - (1)].dnum)));}
    break;

  case 180:
#line 705 "lg.ypp"
    {(yyval.cexp)= CConstant<const char *>((yyvsp[(1) - (1)].str));}
    break;

  case 181:
#line 710 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (4)].cexp),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args));;}
    break;

  case 182:
#line 712 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (4)].cexp),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].cexp));}
    break;

  case 183:
#line 713 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (6)].cexp),(yyvsp[(2) - (6)].oper),(yyvsp[(3) - (6)].cexp),(yyvsp[(5) - (6)].cexp));}
    break;

  case 184:
#line 714 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),"[]");}
    break;

  case 185:
#line 715 "lg.ypp"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 186:
#line 716 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 187:
#line 717 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 188:
#line 718 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 189:
#line 719 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 190:
#line 720 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 191:
#line 721 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 192:
#line 722 "lg.ypp"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));}
    break;

  case 193:
#line 723 "lg.ypp"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));}
    break;

  case 194:
#line 724 "lg.ypp"
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

  case 195:
#line 733 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(2) - (3)].cexp);}
    break;

  case 196:
#line 734 "lg.ypp"
    { (yyval.cexp)=C_F0(TheOperators,"[]",(yyvsp[(2) - (3)].args));}
    break;


/* Line 1267 of yacc.c.  */
#line 3088 "lg.tab.cpp"
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


#line 739 "lg.ypp"



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

/// <<Compile>> Called by mainff(). Activates the bison parser by calling yyparse()
int Compile()
{

  // see [[YYSTYPE]] [[yylval]] [[lglval]]
  extern   YYSTYPE *plglval;  // modif FH

  /// plglval is allocated at [[file:../fflib/global.cpp::plglval]]
  plglval = &lglval;

  int retvalue=0;

  currentblock=0;
  Block::open(currentblock);
  try {
    UnShowAlloc =0;

    retvalue=yyparse(); // grammar analysis starting from [[start_symbol]]

    if(retvalue==0){
      if(currentblock)
        {retvalue=1; if(!mpirank) cerr <<  "Error:a block is not close" << endl; }
      else {
        if( verbosity  ) {
	      UnShowAlloc =1;
	      cerr << " CodeAlloc : nb ptr  "<< CodeAlloc::nb << ",  size :"  <<  CodeAlloc::lg
              << " mpirank: " <<mpirank <<  endl    ;
              extern   long npichon2d, npichon3d;
              extern   long npichon2d1, npichon3d1;
              if( npichon2d || npichon3d ) cout << " WARNING NUMBER bad SearchMethod cas in 2d: "
                 <<npichon2d << " int 3d "<< npichon3d << "(essai d2: " <<npichon2d1  <<" 3d: " << npichon3d1 <<" )"<< endl;
	      if(!mpirank) cerr <<  "Ok: Normal End" << endl;
	    }
	}
    }
  }

  catch (Error & e)
    {
      retvalue=e.errcode();
      if(mpirank ==0)
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
extern  bool echo_edp;

/// Called by mymain() and calls Compile() to run the FF language parser
int mainff (int  argc, char **argv)
{

  if(argc)
    prognamearg=argv[0];

 //   int vvold=verbosity;
    if(mpirank !=0) verbosity=0;

  // ALH - 14/10/8 - This breaks FFCS output redirection
#ifndef ENABLE_FFCS
  SetcppIo();
#endif

  GetEnvironment();
//    vvold=verbosity;
    if(mpirank !=0) verbosity=0;
  //  size_t lg000;
 // ShowAlloc("begin main ",lg000);
  int retvalue=0;
   ff_atend(fingraphique);
   if (initparallele)initparallele(argc,argv);

  CPUcompileInit= CPUtime();
  withrgraphique = false;
   atexit(ForDebug);
//  AllFunctions::maptype  xlocal;
//  local=&xlocal;
  lexdebug = false;
  lgdebug = false;

  char *  cc= new char [1024];
  //  istream * ccin=0;
  if ( ! (getprog(cc,argc,argv) >0)  )
    {
      cout << "-- FreeFem++ v" << StrVersionNumber() << " (error parameter!)\n"  ;
      if(ThePlotStream) {ffapi::ff_pclose(ThePlotStream); ThePlotStream=0;}
      return 1;
    }

  if(verbosity && (mpirank==0)) {
      cout << "-- FreeFem++ v" << StrVersionNumber() << endl;
      if(verbosity>1) cout << "   file :" << cc << " " << " verbosity= " << verbosity << endl;
  }

    KN<String> karg(argc);
    for(int i=0;i< argc;++i)
	karg[i]=argv[i];
    pkarg= &karg;

    /// <<zzzfff>>
    zzzfff = Newlex(cout,echo_edp);


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
 //  zzzfff->Add("forall",FORALL);
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
//   Init_map_type();
   if(verbosity>2 || ( (mpirank==0) && verbosity ) ) cout << " Load: ";
   callInitsFunct() ; //  init for dynamique libs ...
  // init_lgfem() ;
   init_lgmesh() ;
   init_lgmesh3() ;
   init_algo();

#ifdef HAVE_LIBARPACK
   init_eigenvalue();
#endif

   if(init_lgparallele)  init_lgparallele();
  //  callInitsFunct() ; //  init for dynamique libs ...

   if(verbosity>2 || ((mpirank==0)&& verbosity)  )  cout << endl;
  zzzfff->input(cc); // [[file:../fflib/lex.cpp::void mylex input]]
  EnvironmentLoad(); // just before compile
  //verbosity=vvold;

  retvalue= Compile(); // [[Compile]]
   // cout << " xxxxx " <<  retvalue << " " << ThePlotStream << endl;

  //if(end_parallele) end_parallele();
  ff_finalize();
  //  currentblock->close(currentblock).eval(thestack);
 // fingraphique();
  // FFCS: divert stream to FFCS
  if(ThePlotStream){
    ffapi::ff_pclose(ThePlotStream);
    ThePlotStream=0;
  }
  Destroylex( zzzfff);
  delete [] cc;
   // ClearMem();
  return retvalue;
}

/* FFCS: emacs configuration for this file */

/*!
 * Local Variables:
 * mode:antlr
 * ispell-local-dictionary:"british"
 * coding:utf-8
 * End:
 */

