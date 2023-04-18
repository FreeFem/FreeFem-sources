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
     FESPACES = 298,
     FESPACEL = 299,
     VGFESPACE = 300,
     GFESPACE = 301,
     PLUSEQ = 302,
     MOINSEQ = 303,
     MULEQ = 304,
     DIVEQ = 305,
     DOTMULEQ = 306,
     DOTDIVEQ = 307,
     ARROW = 308,
     BORDER = 309,
     SOLVE = 310
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
#define FESPACES 298
#define FESPACEL 299
#define VGFESPACE 300
#define GFESPACE 301
#define PLUSEQ 302
#define MOINSEQ 303
#define MULEQ 304
#define DIVEQ 305
#define DOTMULEQ 306
#define DOTDIVEQ 307
#define ARROW 308
#define BORDER 309
#define SOLVE 310




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
/*
  bug of clang if optimisation !!!!!
*/
#ifdef   __clang__
#pragma clang optimize off
#endif

#include <config.h>
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

#include <csignal>
void signalCPUHandler( int signum ) {
    ff_finalize();
    std::cout << "Cputime limit exceeded:  (" << signum << ") received.\n";
    
    exit(24);
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
#line 169 "lg.ypp"
{
 double dnum;

 /* <<YYSTYPE_lnum>> */
 long   lnum;// to read long long number !!!! FH dec 2022

 /* <<YYSTYPE_str>> */
 char * str;
 char oper[8];

 /* <<YYSTYPE_cexp>> [[file:../fflib/AFunction.hpp::CC_F0]] */
 CC_F0 cexp;

 Routine   *routine;

 /* <<YYSTYPE_args>> [[file:~/ff/src/fflib/AFunction.hpp::AC_F0]] */
 AC_F0 args;

 /* <<YYSTYPE_type>> refers to [[file:~/ff/src/fflib/AnyType.hpp::aType]] */
 aType type;

 /* <<YYSTYPE_cinst>> refers to [[file:~/ff/src/fflib/AFunction.hpp::CListOfInst]] */
 CListOfInst cinst;

 Block * block;

 /* <<YYSTYPE_clist_id>> [[file:~/ff/src/fflib/AFunction.hpp::ListOfId]] */
 ListOfId *clist_id;

/* ListCatch * clist_Catchs;*/

 vectorOfInst * endb;
}
/* Line 193 of yacc.c.  */
#line 411 "lg.tab.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 424 "lg.tab.cpp"

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
#define YYFINAL  100
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1274

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  81
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  50
/* YYNRULES -- Number of rules.  */
#define YYNRULES  252
/* YYNRULES -- Number of states.  */
#define YYNSTATES  507

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   310

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
       2,     2,     2,     2,     2,     2,     2,     2,    79,    76,
      16,     6,    17,    80,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    35,     2,    38,    31,    33,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    77,    10,    78,     2,     2,     2,     2,
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
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     6,     8,    10,    13,    14,    16,    20,
      23,    27,    30,    34,    37,    41,    44,    48,    51,    55,
      59,    65,    70,    74,    80,    87,    95,   102,   108,   113,
     119,   124,   130,   135,   141,   146,   152,   157,   163,   165,
     169,   171,   173,   175,   177,   179,   181,   183,   187,   192,
     196,   198,   201,   204,   207,   210,   213,   215,   217,   219,
     221,   223,   227,   231,   233,   238,   246,   253,   263,   268,
     276,   286,   288,   293,   297,   301,   308,   314,   319,   326,
     328,   330,   332,   334,   336,   338,   340,   342,   347,   349,
     353,   355,   359,   362,   368,   373,   377,   379,   383,   384,
     389,   393,   396,   402,   403,   414,   415,   425,   427,   429,
     431,   433,   434,   438,   440,   442,   446,   452,   454,   457,
     460,   466,   469,   471,   472,   481,   491,   501,   507,   513,
     521,   525,   529,   536,   539,   542,   546,   554,   562,   572,
     575,   577,   581,   583,   585,   587,   589,   591,   593,   597,
     601,   605,   609,   613,   617,   621,   623,   629,   633,   639,
     641,   645,   649,   653,   657,   661,   665,   669,   673,   677,
     681,   685,   689,   693,   697,   701,   705,   709,   713,   717,
     719,   721,   725,   731,   735,   736,   738,   740,   742,   744,
     746,   748,   752,   754,   758,   762,   766,   770,   774,   778,
     782,   788,   790,   794,   796,   798,   800,   802,   804,   808,
     812,   816,   820,   824,   828,   832,   836,   840,   844,   846,
     849,   851,   855,   859,   861,   864,   866,   868,   870,   872,
     874,   879,   884,   891,   895,   899,   903,   908,   912,   917,
     921,   926,   930,   935,   939,   944,   947,   950,   955,   960,
     964,   968,   972
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
      82,     0,    -1,    83,    46,    -1,    84,    -1,   112,    -1,
      84,   112,    -1,    -1,    87,    -1,    87,     6,   119,    -1,
      60,    87,    -1,    60,    12,    87,    -1,    62,    87,    -1,
      62,    12,    87,    -1,    63,    87,    -1,    63,    12,    87,
      -1,    64,    87,    -1,    64,    12,    87,    -1,    90,    87,
      -1,    90,    12,    87,    -1,    35,    85,    38,    -1,    16,
      35,    85,    38,    17,    -1,    16,    35,    85,    38,    -1,
      85,     5,    87,    -1,    85,     5,    35,    85,    38,    -1,
      85,     5,    35,    85,    38,    17,    -1,    85,     5,    16,
      35,    85,    38,    17,    -1,    85,     5,    16,    35,    85,
      38,    -1,    85,     5,    87,     6,   119,    -1,    85,     5,
      60,    87,    -1,    85,     5,    60,    12,    87,    -1,    85,
       5,    62,    87,    -1,    85,     5,    62,    12,    87,    -1,
      85,     5,    63,    87,    -1,    85,     5,    63,    12,    87,
      -1,    85,     5,    64,    87,    -1,    85,     5,    64,    12,
      87,    -1,    85,     5,    90,    87,    -1,    85,     5,    90,
      12,    87,    -1,    87,    -1,    86,     5,    87,    -1,    42,
      -1,    60,    -1,    62,    -1,    63,    -1,    64,    -1,    61,
      -1,    42,    -1,    42,     6,   119,    -1,    42,    34,    89,
      37,    -1,    88,     5,    88,    -1,   120,    -1,    60,    42,
      -1,    61,    42,    -1,    62,    42,    -1,    63,    42,    -1,
      64,    42,    -1,    60,    -1,    61,    -1,    62,    -1,    63,
      -1,    64,    -1,    42,     6,   120,    -1,    89,     5,    89,
      -1,    58,    -1,    58,    35,    58,    38,    -1,    58,    35,
      58,    38,    35,    58,    38,    -1,    58,    35,    58,     5,
      58,    38,    -1,    58,    35,    58,     5,    58,    38,    35,
      58,    38,    -1,    58,    16,    58,    17,    -1,    58,    16,
      58,    17,    35,    58,    38,    -1,    58,    16,    58,    17,
      35,    58,     5,    58,    38,    -1,    42,    -1,    42,    35,
     120,    38,    -1,    42,     6,   120,    -1,    35,    86,    38,
      -1,    35,    86,    38,    35,   120,    38,    -1,    35,    86,
      38,     6,   120,    -1,    42,    34,   120,    37,    -1,    35,
      86,    38,    34,   120,    37,    -1,    60,    -1,    61,    -1,
      62,    -1,    63,    -1,    64,    -1,    65,    -1,    66,    -1,
      93,    -1,    93,    16,    58,    17,    -1,    92,    -1,    95,
       5,    92,    -1,    91,    -1,    96,     5,    91,    -1,    94,
      96,    -1,    94,    35,    58,    38,    95,    -1,    42,    34,
      89,    37,    -1,    42,     6,    89,    -1,    98,    -1,    99,
       5,    98,    -1,    -1,    90,   101,    88,    76,    -1,    43,
      99,    76,    -1,    97,    76,    -1,    59,    42,     6,   117,
      76,    -1,    -1,    59,    90,    42,    34,    85,    37,   102,
      77,    84,    78,    -1,    -1,    59,    42,    34,    85,    37,
     103,     6,   119,    76,    -1,    77,    -1,    78,    -1,    50,
      -1,    51,    -1,    -1,    90,   109,    88,    -1,    55,    -1,
      87,    -1,    87,     5,    87,    -1,    87,     5,    87,     5,
      87,    -1,    76,    -1,    47,    45,    -1,    48,    45,    -1,
     110,    77,    84,    78,   114,    -1,   117,    76,    -1,   100,
      -1,    -1,   106,    35,   111,    79,   130,    38,   113,   112,
      -1,   106,    34,   117,    76,   117,    76,   117,    37,   112,
      -1,   106,    34,   108,    76,   117,    76,   117,    37,   112,
      -1,   107,    34,   117,    37,   112,    -1,     3,    34,   117,
      37,   112,    -1,     3,    34,   117,    37,   112,     4,   112,
      -1,   104,    84,   105,    -1,    74,    42,   116,    -1,    74,
      42,    35,   125,    38,    76,    -1,    52,    76,    -1,    53,
      76,    -1,    54,   117,    76,    -1,    56,    34,    36,    36,
      36,    37,   112,    -1,    34,    42,     6,   117,     5,   117,
      37,    -1,    34,    42,     6,   117,     5,   117,    76,    42,
      37,    -1,   115,   112,    -1,   119,    -1,   117,     5,   117,
      -1,    21,    -1,    20,    -1,    27,    -1,    29,    -1,    28,
      -1,   120,    -1,   120,     6,   119,    -1,   120,    67,   119,
      -1,   120,    68,   119,    -1,   120,    69,   119,    -1,   120,
      70,   119,    -1,   120,    71,   119,    -1,   120,    72,   119,
      -1,   121,    -1,   121,    80,   121,    79,   121,    -1,   121,
      79,   121,    -1,   121,    79,   121,    79,   121,    -1,   127,
      -1,   121,    22,   121,    -1,   121,    26,   121,    -1,   121,
      25,   121,    -1,   121,    23,   121,    -1,   121,    24,   121,
      -1,   121,    20,   121,    -1,   121,    21,   121,    -1,   121,
       9,   121,    -1,   121,     8,   121,    -1,   121,    12,   121,
      -1,   121,    13,   121,    -1,   121,    10,   121,    -1,   121,
      11,   121,    -1,   121,    16,   121,    -1,   121,    19,   121,
      -1,   121,    17,   121,    -1,   121,    18,   121,    -1,   121,
      15,   121,    -1,   121,    14,   121,    -1,   121,    -1,    79,
      -1,   121,    79,   121,    -1,   121,    79,   121,    79,   121,
      -1,   122,     5,   124,    -1,    -1,    60,    -1,    61,    -1,
      62,    -1,    63,    -1,    64,    -1,    65,    -1,    87,     6,
     120,    -1,   122,    -1,   124,     5,    60,    -1,   124,     5,
      61,    -1,   124,     5,    62,    -1,   124,     5,    63,    -1,
     124,     5,    64,    -1,   124,     5,    65,    -1,   124,     5,
     122,    -1,   124,     5,    87,     6,   120,    -1,   119,    -1,
     125,     5,   119,    -1,    60,    -1,    61,    -1,    62,    -1,
      63,    -1,    64,    -1,   126,     5,    60,    -1,   126,     5,
      61,    -1,   126,     5,    62,    -1,   126,     5,    63,    -1,
     126,     5,    64,    -1,   126,    22,    60,    -1,   126,    22,
      61,    -1,   126,    22,    62,    -1,   126,    22,    63,    -1,
     126,    22,    64,    -1,   128,    -1,   118,   128,    -1,   129,
      -1,   129,    31,   127,    -1,   129,    33,   127,    -1,   130,
      -1,   130,    32,    -1,    42,    -1,    39,    -1,    40,    -1,
      41,    -1,    45,    -1,   130,    34,   124,    37,    -1,   130,
      35,   122,    38,    -1,   130,    35,   122,     5,   122,    38,
      -1,   130,    35,    38,    -1,   130,    36,    42,    -1,    60,
      36,    42,    -1,    60,    34,   124,    37,    -1,    61,    36,
      42,    -1,    61,    34,   124,    37,    -1,    62,    36,    42,
      -1,    62,    34,   124,    37,    -1,    63,    36,    42,    -1,
      63,    34,   124,    37,    -1,    64,    36,    42,    -1,    64,
      34,   124,    37,    -1,   130,    29,    -1,   130,    28,    -1,
      58,    34,   122,    37,    -1,    58,    34,   123,    37,    -1,
      34,   117,    37,    -1,    35,   125,    38,    -1,    16,   126,
      17,    -1,   126,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   338,   338,   408,   412,   413,   419,   420,   421,   422,
     423,   424,   425,   426,   427,   428,   429,   430,   431,   432,
     433,   434,   435,   436,   437,   438,   439,   440,   441,   442,
     443,   444,   445,   446,   447,   448,   449,   450,   453,   454,
     459,   459,   459,   459,   459,   459,   462,   463,   464,   465,
     471,   472,   473,   474,   475,   476,   477,   478,   479,   480,
     481,   482,   483,   488,   489,   490,   491,   492,   493,   494,
     495,   499,   500,   501,   502,   503,   504,   508,   509,   514,
     515,   516,   517,   518,   519,   520,   523,   524,   529,   530,
     532,   533,   535,   536,   539,   542,   546,   547,   550,   550,
     551,   552,   553,   555,   554,   571,   570,   580,   581,   585,
     587,   591,   591,   594,   596,   597,   598,   600,   601,   602,
     603,   604,   605,   607,   606,   612,   613,   617,   618,   619,
     620,   625,   627,   630,   634,   638,   645,   648,   656,   664,
     671,   672,   676,   677,   678,   679,   680,   684,   685,   686,
     687,   688,   689,   690,   691,   696,   697,   698,   699,   703,
     704,   705,   706,   707,   708,   709,   710,   711,   712,   713,
     714,   715,   716,   717,   718,   719,   720,   721,   722,   726,
     727,   728,   729,   734,   738,   739,   740,   741,   742,   743,
     744,   745,   746,   747,   748,   749,   750,   751,   752,   753,
     756,   759,   760,   763,   764,   765,   766,   767,   768,   769,
     770,   771,   772,   773,   774,   775,   776,   777,   781,   782,
     786,   787,   788,   792,   793,   801,   805,   806,   807,   808,
     813,   815,   816,   817,   818,   819,   820,   821,   822,   823,
     824,   825,   826,   827,   828,   829,   830,   831,   840,   848,
     849,   850,   851
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
  "FUNCTION", "FESPACE", "FESPACE1", "FESPACE3", "FESPACES", "FESPACEL",
  "VGFESPACE", "GFESPACE", "PLUSEQ", "MOINSEQ", "MULEQ", "DIVEQ",
  "DOTMULEQ", "DOTDIVEQ", "ARROW", "BORDER", "SOLVE", "';'", "'{'", "'}'",
  "':'", "'?'", "$accept", "start", "input", "instructions",
  "list_of_id_args", "list_of_id1", "id", "list_of_dcls",
  "parameters_list", "type_of_dcl", "ID_space", "ID_array_space",
  "fespace123", "fespace", "spaceIDa", "spaceIDb", "spaceIDs",
  "fespace_def", "fespace_def_list", "declaration", "@1", "@2", "@3",
  "begin", "end", "for_loop", "while_loop", "declaration_for", "@4", "try",
  "IDfor", "instruction", "@5", "catchs", "bornes", "border_expr", "Expr",
  "unop", "no_comma_expr", "no_set_expr", "no_ternary_expr",
  "sub_script_expr", "parameterstype", "parameters", "array", "FEarray",
  "unary_expr", "pow_expr", "primaryp", "primary", 0
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
     305,   306,   307,   308,   309,   310,    59,   123,   125,    58,
      63
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    81,    82,    83,    84,    84,    85,    85,    85,    85,
      85,    85,    85,    85,    85,    85,    85,    85,    85,    85,
      85,    85,    85,    85,    85,    85,    85,    85,    85,    85,
      85,    85,    85,    85,    85,    85,    85,    85,    86,    86,
      87,    87,    87,    87,    87,    87,    88,    88,    88,    88,
      89,    89,    89,    89,    89,    89,    89,    89,    89,    89,
      89,    89,    89,    90,    90,    90,    90,    90,    90,    90,
      90,    91,    91,    91,    91,    91,    91,    92,    92,    93,
      93,    93,    93,    93,    93,    93,    94,    94,    95,    95,
      96,    96,    97,    97,    98,    98,    99,    99,   101,   100,
     100,   100,   100,   102,   100,   103,   100,   104,   105,   106,
     107,   109,   108,   110,   111,   111,   111,   112,   112,   112,
     112,   112,   112,   113,   112,   112,   112,   112,   112,   112,
     112,   112,   112,   112,   112,   112,   114,   115,   115,   116,
     117,   117,   118,   118,   118,   118,   118,   119,   119,   119,
     119,   119,   119,   119,   119,   120,   120,   120,   120,   121,
     121,   121,   121,   121,   121,   121,   121,   121,   121,   121,
     121,   121,   121,   121,   121,   121,   121,   121,   121,   122,
     122,   122,   122,   123,   124,   124,   124,   124,   124,   124,
     124,   124,   124,   124,   124,   124,   124,   124,   124,   124,
     124,   125,   125,   126,   126,   126,   126,   126,   126,   126,
     126,   126,   126,   126,   126,   126,   126,   126,   127,   127,
     128,   128,   128,   129,   129,   130,   130,   130,   130,   130,
     130,   130,   130,   130,   130,   130,   130,   130,   130,   130,
     130,   130,   130,   130,   130,   130,   130,   130,   130,   130,
     130,   130,   130
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     1,     1,     2,     0,     1,     3,     2,
       3,     2,     3,     2,     3,     2,     3,     2,     3,     3,
       5,     4,     3,     5,     6,     7,     6,     5,     4,     5,
       4,     5,     4,     5,     4,     5,     4,     5,     1,     3,
       1,     1,     1,     1,     1,     1,     1,     3,     4,     3,
       1,     2,     2,     2,     2,     2,     1,     1,     1,     1,
       1,     3,     3,     1,     4,     7,     6,     9,     4,     7,
       9,     1,     4,     3,     3,     6,     5,     4,     6,     1,
       1,     1,     1,     1,     1,     1,     1,     4,     1,     3,
       1,     3,     2,     5,     4,     3,     1,     3,     0,     4,
       3,     2,     5,     0,    10,     0,     9,     1,     1,     1,
       1,     0,     3,     1,     1,     3,     5,     1,     2,     2,
       5,     2,     1,     0,     8,     9,     9,     5,     5,     7,
       3,     3,     6,     2,     2,     3,     7,     7,     9,     2,
       1,     3,     1,     1,     1,     1,     1,     1,     3,     3,
       3,     3,     3,     3,     3,     1,     5,     3,     5,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       1,     3,     5,     3,     0,     1,     1,     1,     1,     1,
       1,     3,     1,     3,     3,     3,     3,     3,     3,     3,
       5,     1,     3,     1,     1,     1,     1,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     2,
       1,     3,     3,     1,     2,     1,     1,     1,     1,     1,
       4,     4,     6,     3,     3,     3,     4,     3,     4,     3,
       4,     3,     4,     3,     4,     2,     2,     4,     4,     3,
       3,     3,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,   143,   142,   144,   146,   145,     0,     0,
     226,   227,   228,   225,     0,   229,     0,     0,   109,   110,
       0,     0,     0,   113,    63,     0,   203,   204,   205,   206,
     207,    84,    85,     0,   117,   107,     0,     0,     3,    98,
      86,     0,     0,   122,     0,     0,     0,     0,     4,     0,
       0,   140,   147,   155,   252,   159,   218,   220,   223,     0,
     203,   204,   205,   206,   207,     0,     0,   203,   204,   205,
     206,   207,     0,   201,     0,     0,    96,     0,   118,   119,
     133,   134,     0,     0,     0,     0,     0,    63,     0,   184,
       0,   184,     0,   184,     0,   184,     0,   184,     0,     0,
       1,     2,     5,     0,     0,     0,    71,    90,    92,   101,
       0,     0,     0,     0,     0,     0,   121,   219,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     246,   245,   224,   184,     0,     0,     0,   251,   249,     0,
     250,     0,     0,     0,   100,   135,     0,   180,   179,     0,
       0,     0,     0,     6,     0,   225,   203,   204,   205,   206,
     207,   190,     0,   192,     0,   235,     0,   237,     0,   239,
       0,   241,     0,   243,     0,     0,     0,   131,    46,     0,
       0,    40,     0,    41,    45,    42,    43,    44,     0,    38,
       0,     0,     0,   108,   130,   111,     0,     0,   114,     0,
       0,     0,   141,   148,   149,   150,   151,   152,   153,   154,
     168,   167,   171,   172,   169,   170,   178,   177,   173,   175,
     176,   174,   165,   166,   160,   163,   164,   162,   161,   157,
       0,   208,   209,   210,   211,   212,   213,   214,   215,   216,
     217,   221,   222,     0,   233,     0,   234,     0,   202,   225,
     203,   204,   205,   206,   207,    95,    50,     0,    97,    68,
       0,   184,   247,   248,     0,    64,     0,     0,     6,    41,
      42,    43,    44,     0,     7,     0,     6,     0,     0,   236,
     238,   240,   242,   244,     0,     0,   139,     0,     0,     0,
      99,    87,     0,     0,    74,    73,     0,     0,    91,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   230,     0,
     231,   128,     0,    51,    52,    53,    54,    55,     0,    94,
       0,   181,   183,     0,     0,   102,     6,     0,     0,     9,
       0,    11,     0,    13,     0,    15,     0,   105,     0,     0,
      17,     0,   191,   203,   204,   205,   206,   207,   198,     0,
     199,     0,     0,    47,     0,    49,     0,     0,    88,    93,
      39,     0,     0,    72,   112,     0,     0,   115,     0,   127,
       0,   120,   158,   156,     0,     0,    61,    62,     0,     0,
      66,     0,     0,    19,    10,    12,    14,    16,     0,     6,
      41,    42,    43,    44,    22,     0,     0,     8,    18,   103,
       0,     0,   132,    48,     0,     0,     0,    76,     0,     0,
       0,     0,   123,     0,   232,   129,     0,    69,   182,     0,
      65,    21,     6,     0,     0,    28,     0,    30,     0,    32,
       0,    34,     0,     0,    36,     0,     0,   200,     0,     0,
       0,    89,    75,     0,     0,   116,     0,     0,     0,     0,
      20,     0,    23,    29,    31,    33,    35,    27,    37,     0,
       0,   141,     0,    77,     0,     0,   124,     0,    70,    67,
      26,    24,   106,     0,   137,     0,     0,   126,   125,     0,
      25,   104,     0,    78,     0,   138,   136
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    36,    37,    38,   293,   208,   182,   199,   275,    39,
     107,   378,    40,    41,   379,   108,    42,    76,    77,    43,
     103,   456,   416,    44,   214,    45,    46,   216,   319,    47,
     219,    48,   466,   391,   196,   197,    49,    50,    51,    52,
      53,   183,   170,   184,    74,    54,    55,    56,    57,    58
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -284
static const yytype_int16 yypact[] =
{
     697,    77,   642,  -284,  -284,  -284,  -284,  -284,  1020,  1020,
    -284,  -284,  -284,  -284,    -6,  -284,    31,    69,  -284,  -284,
     -15,    73,  1020,  -284,   285,     0,   402,   442,   523,   621,
     685,  -284,  -284,   121,  -284,  -284,   197,   140,   697,  -284,
     184,    -5,   134,  -284,   697,   209,   179,   149,  -284,    16,
    1060,  -284,    99,   119,   102,  -284,  -284,    70,   806,  1020,
    -284,  -284,  -284,  -284,  -284,   245,   200,   139,   218,   246,
     275,   276,    34,  -284,     2,    11,  -284,    17,  -284,  -284,
    -284,  -284,    18,   191,   971,   208,    78,   268,   230,   816,
     253,   816,   265,   816,   292,   816,   293,   816,   294,   287,
    -284,  -284,  -284,   295,   260,   907,   154,  -284,   234,  -284,
     440,  1069,   332,  1020,   697,  1020,  -284,  -284,  1020,  1020,
    1020,  1020,  1020,  1020,  1020,  1020,  1020,  1020,  1020,  1020,
    1020,  1020,  1020,  1020,  1020,  1020,  1020,  1020,  1020,  1020,
    1020,  1020,  1020,  1020,  1020,  1020,  1080,  1205,  1020,  1020,
    -284,  -284,  -284,   816,   869,   301,    50,  -284,  -284,  1020,
    -284,  1118,  1118,    -6,  -284,  -284,   329,  -284,   511,    52,
     312,    15,  1020,   959,   328,   360,   117,   142,   227,   311,
     353,  -284,   362,  -284,    76,  -284,   147,  -284,   157,  -284,
     167,  -284,   172,  -284,   330,  1020,   697,  -284,   182,    19,
     359,  -284,   340,  -284,  -284,  -284,  -284,  -284,    24,  -284,
    1020,  1020,   252,  -284,  -284,  -284,   304,    20,   365,   302,
     175,   568,  -284,  -284,  -284,  -284,  -284,  -284,  -284,  -284,
    1192,  1192,  1207,  1207,  1220,  1220,  1148,  1148,  1238,  1238,
    1238,  1238,  1248,  1248,  -284,  -284,  -284,  -284,  -284,   767,
     786,  -284,  -284,  -284,  -284,  -284,  -284,  -284,  -284,  -284,
    -284,  -284,  -284,   177,  -284,    26,  -284,   697,  -284,   377,
      14,    33,    49,   159,   169,  -284,  -284,   186,  -284,   349,
    1020,   816,  -284,  -284,   327,   356,    21,   385,   959,   195,
     229,   236,   264,   203,   382,   290,   959,  1020,   918,  -284,
    -284,  -284,  -284,  -284,   415,    27,  -284,  1020,  1118,   295,
    -284,  -284,   273,   332,   155,  -284,   386,   332,  -284,   295,
    1020,  1020,   332,  1060,   697,   371,  1020,  1020,  -284,   971,
    -284,   425,  1020,  -284,  -284,  -284,  -284,  -284,  1118,  -284,
     372,   805,   426,   394,   387,  -284,   959,    28,   332,  -284,
     332,  -284,   332,  -284,   332,  -284,  1009,  -284,  1020,   332,
    -284,   210,  -284,   428,   436,   538,   542,   549,  -284,   429,
    -284,  1020,   364,  -284,   223,  -284,   332,   414,  -284,   445,
    -284,  1020,  1020,  -284,   448,    22,    23,   449,  1219,  -284,
     423,  -284,  1175,  1175,   433,   697,  -284,  -284,    30,  1020,
     451,   458,    39,  -284,  -284,  -284,  -284,  -284,   462,   959,
     447,   500,   704,   857,   457,   880,   507,  -284,  -284,  -284,
    1020,   510,  -284,  -284,    41,  1020,   273,  -284,   502,  1020,
    1020,   332,  -284,   505,  -284,  -284,   480,  -284,  1175,   487,
    -284,   529,   959,    42,   332,  -284,   332,  -284,   332,  -284,
     332,  -284,  1020,   332,  -284,  1020,   472,  -284,  1020,   516,
     514,  -284,  -284,   231,   233,  -284,   697,   520,   528,   530,
    -284,    44,   550,  -284,  -284,  -284,  -284,  -284,  -284,   493,
     697,   -33,  1020,  -284,   697,   697,  -284,   537,  -284,  -284,
     560,  -284,  -284,   633,  -284,   539,   543,  -284,  -284,   545,
    -284,  -284,   554,  -284,   697,  -284,  -284
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -284,  -284,  -284,   -41,  -283,   211,   -71,  -209,  -153,   -23,
     380,   168,  -284,  -284,  -284,  -284,  -284,   430,  -284,  -284,
    -284,  -284,  -284,  -284,  -284,  -284,  -284,  -284,  -284,  -284,
    -284,   -38,  -284,  -284,  -284,  -284,    -7,  -284,    -3,  -151,
     272,   -76,  -284,   -79,   405,   602,   181,   555,  -284,   283
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -198
static const yytype_int16 yytable[] =
{
     102,    72,    88,   110,   494,   347,    73,   159,   169,   277,
     276,   276,   186,   361,   188,    82,   190,   161,   192,   -56,
     284,   115,   163,   115,   309,   115,   115,   115,   115,   313,
     105,   329,   159,   356,   209,   436,    75,   106,   -57,   115,
     160,   218,    86,   495,   356,   162,   313,   356,    89,   356,
      90,   -56,   156,   285,   -58,   115,   333,   281,    87,   315,
     316,    80,   314,   402,   330,   372,   403,    91,   437,    92,
     -57,   158,   102,   221,   263,   334,    78,   441,   265,   459,
     472,   298,   490,    93,   172,    94,   -58,   267,   215,   282,
     -56,   335,   116,   164,   165,   310,   321,   345,   429,   430,
     375,   148,   294,   149,   217,   118,   220,   146,   222,   -57,
     384,    59,   173,   299,    79,   223,   224,   225,   226,   227,
     228,   229,  -185,   -41,   147,   -58,   443,   125,   126,   127,
     128,   129,   130,   131,   132,   133,   134,   135,   136,   137,
     138,   139,   140,   141,   142,   143,   362,  -186,   -45,    81,
     295,    89,   298,    90,  -185,   374,   268,   276,   306,   471,
     210,   381,   298,    99,   -59,   286,   119,   120,   121,   122,
     123,   124,   298,    89,   -60,    90,    91,   298,    92,  -186,
     115,   396,   298,   102,   300,   397,   101,   276,   307,   211,
     382,   338,    73,    95,   301,    96,   -59,   100,   144,   145,
     104,   336,   342,    97,   302,    98,   -60,   348,   356,   303,
     109,   337,   324,   113,   328,   356,   308,   294,   349,   351,
     353,   355,   370,   339,   360,   294,   114,   369,   338,   331,
     427,   428,  -187,   -42,    84,   -59,   115,   201,   115,   212,
     357,   350,   380,   111,   112,   -60,   209,   419,   352,   166,
     146,   387,    91,   394,    92,   203,   204,   205,   206,   207,
     423,    93,   157,    94,  -187,   295,   171,   147,   484,   457,
     485,   201,   174,   295,   460,   294,   354,   404,   201,   405,
      93,   406,    94,   407,    83,   414,   389,   317,   418,   203,
     204,   205,   206,   207,   106,   185,   203,   204,   205,   206,
     207,    83,   359,    85,   373,   209,   201,   187,   376,    95,
      97,    96,    98,   385,   386,   377,  -188,   -43,   200,    84,
      85,   194,   195,   295,   203,   204,   205,   206,   207,   261,
     262,   496,   201,   415,   189,   191,   193,   198,   294,   445,
     447,   449,   451,   266,   454,    95,   279,    96,  -188,   283,
     203,   204,   205,   206,   207,   417,   168,   435,  -189,   -44,
     465,   168,   296,   168,   421,   168,   -40,   168,   297,   168,
     322,   294,   304,   473,   201,   474,   311,   475,   312,   476,
     320,   323,   478,   332,   340,   343,   295,    97,   358,    98,
    -189,   344,   203,   204,   205,   206,   207,   230,   231,   232,
     233,   234,   235,   236,   237,   238,   239,   240,   241,   242,
     243,   244,   245,   246,   247,   248,   249,   250,   -79,   295,
     346,   371,   463,   464,   383,   168,   168,   390,   486,   395,
     398,   298,   400,  -193,   -41,   420,    89,   -79,    90,   493,
     422,  -194,   -45,     1,   -79,   401,   497,   498,   425,   477,
     426,   481,   479,   309,   431,   102,     2,   433,   -80,   444,
       3,     4,    89,   452,    90,  -193,   506,     5,     6,     7,
      91,   434,    92,  -194,     8,     9,    91,   -80,    92,    10,
      11,    12,    13,    14,   -80,    15,   439,    16,    17,   201,
      18,    19,    20,    21,    22,    23,   440,   442,    24,    25,
      26,    27,    28,    29,    30,    31,    32,   203,   204,   205,
     206,   207,   446,   455,    33,   458,    34,    35,   213,   125,
     126,   127,   128,   129,   130,   131,   132,   133,   134,   135,
     136,   137,   138,   139,   140,   141,   142,   143,   468,   -81,
     462,   467,   201,  -195,   -42,   469,   470,  -196,   -43,   480,
     482,   483,   341,   168,  -197,   -44,   487,    93,   -81,    94,
     203,   204,   205,   206,   207,   -81,   488,   491,   489,   492,
     168,     1,    93,   499,    94,  -195,    95,   500,    96,  -196,
     503,   502,   504,    97,     2,    98,  -197,   424,     3,     4,
     280,   505,   318,   278,   461,     5,     6,     7,   392,   393,
     305,   168,     8,     9,    65,   117,   388,    10,    11,    12,
      13,    14,     0,    15,     0,    16,    17,     0,    18,    19,
      20,    21,    22,    23,     0,     0,    24,    25,    26,    27,
      28,    29,    30,    31,    32,     0,     1,   -82,     0,     0,
       0,     0,    33,     0,    34,    35,   325,     0,     0,     2,
       0,     0,     0,     3,     4,    95,   -82,    96,     0,     0,
       5,     6,     7,   -82,     0,     0,     0,     8,     9,     0,
       0,   438,    10,    11,    12,    13,    14,     0,    15,     0,
      16,    17,     0,    18,    19,    20,    21,    22,    23,     0,
       0,    24,    25,    26,    27,    28,    29,    30,    31,    32,
       1,   -83,    60,    61,    62,    63,    64,    33,     0,    34,
      35,   501,     0,     2,     0,     0,   448,     3,     4,    97,
     -83,    98,     0,     0,     5,     6,     7,   -83,     0,     0,
       0,     8,     9,     0,     0,     0,    10,    11,    12,    13,
      14,     0,    15,     0,    16,    17,   201,    18,    19,    20,
      21,    22,    23,     0,     0,    24,    25,    26,    27,    28,
      29,    30,    31,    32,   203,   204,   205,   206,   207,     0,
       0,    33,     0,    34,    35,   125,   126,   127,   128,   129,
     130,   131,   132,   133,   134,   135,   136,   137,   138,   139,
     140,   141,   142,   143,   125,   126,   127,   128,   129,   130,
     131,   132,   133,   134,   135,   136,   137,   138,   139,   140,
     141,   142,   143,   125,   126,   127,   128,   129,   130,   131,
     132,   133,   134,   135,   136,   137,   138,   139,   140,   141,
     142,   143,     2,     0,   150,   151,     3,     4,   152,     0,
     153,   154,   155,     5,     6,     7,   326,     0,     0,     0,
       8,     9,     0,     0,     0,    10,    11,    12,   175,     0,
       0,    15,     0,     0,     0,   327,     0,     0,     0,   450,
       0,     0,     0,     0,    66,     0,   176,   177,   178,   179,
     180,   181,     0,     0,   399,     2,     0,     0,     0,     3,
       4,     0,   453,     0,     0,   167,     5,     6,     7,   201,
       0,     0,     0,     8,     9,     0,     0,   264,    10,    11,
      12,    13,     0,     0,    15,     0,     0,   203,   204,   205,
     206,   207,   201,     0,     0,     0,     0,    66,     0,    67,
      68,    69,    70,    71,     2,     0,     0,     0,     3,     4,
     203,   204,   205,   206,   207,     5,     6,     7,   167,   201,
       0,     0,     8,     9,     0,     0,     0,    10,    11,    12,
     175,     0,     0,    15,     0,   202,     0,   203,   204,   205,
     206,   207,     0,     0,     0,   287,    66,     0,   363,   364,
     365,   366,   367,   368,     0,     0,     0,     2,     0,     0,
       0,     3,     4,     0,   288,     0,     0,   167,     5,     6,
       7,   201,     0,     0,     0,     8,     9,     0,     0,     0,
      10,    11,    12,    13,     0,     0,    15,    87,     0,   289,
     204,   290,   291,   292,     0,   408,     0,     0,     0,    66,
       0,    67,    68,    69,    70,    71,     2,     0,     0,     0,
       3,     4,     0,     0,   409,     0,     0,     5,     6,     7,
     167,   201,     0,     0,     8,     9,     0,     0,     0,    10,
      11,    12,    13,     0,     0,    15,     0,    87,     0,   410,
     204,   411,   412,   413,     0,     0,     2,     0,    66,     0,
      67,    68,    69,    70,    71,     2,     0,     0,     0,     3,
       4,     0,     0,     0,     8,     9,     5,     6,     7,    10,
      11,    12,    13,     8,     9,    15,     0,     0,    10,    11,
      12,    13,     0,     0,    15,     0,     0,     0,    66,     0,
      67,    68,    69,    70,    71,     0,     0,    24,     0,    67,
      68,    69,    70,    71,     2,     0,     0,     0,     3,     4,
     251,   252,   253,   254,   255,     5,     6,     7,     0,     0,
       0,     0,     8,     9,     0,     0,     0,    10,    11,    12,
     269,     0,     0,    15,   133,   134,   135,   136,   137,   138,
     139,   140,   141,   142,   143,     0,    66,     0,   270,   271,
     272,   273,   274,   125,   126,   127,   128,   129,   130,   131,
     132,   133,   134,   135,   136,   137,   138,   139,   140,   141,
     142,   143,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   129,
     130,   131,   132,   133,   134,   135,   136,   137,   138,   139,
     140,   141,   142,   143,   131,   132,   133,   134,   135,   136,
     137,   138,   139,   140,   141,   142,   143,   150,   151,     0,
       0,     0,     0,   153,   154,   155,     0,   432,   137,   138,
     139,   140,   141,   142,   143,   256,   257,   258,   259,   260,
     139,   140,   141,   142,   143
};

static const yytype_int16 yycheck[] =
{
      38,     8,    25,    44,    37,   288,     9,     5,    84,   162,
     161,   162,    91,   296,    93,    22,    95,     6,    97,     5,
       5,     5,     5,     5,     5,     5,     5,     5,     5,     5,
      35,     5,     5,     5,   105,     5,    42,    42,     5,     5,
      38,   112,    42,    76,     5,    34,     5,     5,    34,     5,
      36,    37,    59,    38,     5,     5,    42,     5,    58,   210,
     211,    76,    38,   346,    38,    38,    38,    34,    38,    36,
      37,    37,   110,   114,   153,    42,    45,    38,   154,    38,
      38,     5,    38,    34,     6,    36,    37,    37,   111,    37,
      76,    42,    76,    76,    76,    76,    76,    76,    76,    76,
     309,    31,   173,    33,   111,     6,   113,     5,   115,    76,
     319,    34,    34,    37,    45,   118,   119,   120,   121,   122,
     123,   124,     5,     6,    22,    76,   409,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,   297,     5,     6,    76,
     173,    34,     5,    36,    37,   308,   159,   308,   196,   442,
       6,     6,     5,    42,     5,   172,    67,    68,    69,    70,
      71,    72,     5,    34,     5,    36,    34,     5,    36,    37,
       5,   332,     5,   221,    37,   338,    46,   338,     6,    35,
      35,     5,   195,    34,    37,    36,    37,     0,    79,    80,
      16,    42,   281,    34,    37,    36,    37,    12,     5,    37,
      76,    42,    37,    34,    37,     5,    34,   288,   289,   290,
     291,   292,   298,    37,   295,   296,    77,   298,     5,   267,
     381,   382,     5,     6,    34,    76,     5,    42,     5,     5,
      37,    12,   313,    34,    35,    76,   317,    37,    12,    58,
       5,   322,    34,   329,    36,    60,    61,    62,    63,    64,
      37,    34,    17,    36,    37,   288,    58,    22,    37,   420,
      37,    42,    42,   296,   425,   346,    12,   348,    42,   350,
      34,   352,    36,   354,    16,   356,   324,    35,   359,    60,
      61,    62,    63,    64,    42,    42,    60,    61,    62,    63,
      64,    16,    12,    35,   307,   376,    42,    42,    35,    34,
      34,    36,    36,   320,   321,    42,     5,     6,    58,    34,
      35,    34,    35,   346,    60,    61,    62,    63,    64,   148,
     149,   482,    42,   356,    42,    42,    42,    42,   409,   410,
     411,   412,   413,    42,   415,    34,    17,    36,    37,    37,
      60,    61,    62,    63,    64,   358,    84,   395,     5,     6,
     431,    89,    34,    91,   371,    93,     6,    95,     6,    97,
       5,   442,    42,   444,    42,   446,    17,   448,    38,   450,
      76,    79,   453,     6,    35,    58,   409,    34,     6,    36,
      37,    35,    60,    61,    62,    63,    64,   125,   126,   127,
     128,   129,   130,   131,   132,   133,   134,   135,   136,   137,
     138,   139,   140,   141,   142,   143,   144,   145,    16,   442,
      35,     6,   429,   430,    38,   153,   154,    56,   466,     4,
      58,     5,    38,     5,     6,     6,    34,    35,    36,   480,
      76,     5,     6,     3,    42,    58,   484,   485,    34,   452,
       5,   458,   455,     5,     5,   493,    16,    34,    16,    12,
      20,    21,    34,     6,    36,    37,   504,    27,    28,    29,
      34,    38,    36,    37,    34,    35,    34,    35,    36,    39,
      40,    41,    42,    43,    42,    45,    35,    47,    48,    42,
      50,    51,    52,    53,    54,    55,    38,    35,    58,    59,
      60,    61,    62,    63,    64,    65,    66,    60,    61,    62,
      63,    64,    12,     6,    74,     5,    76,    77,    78,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    58,    16,
      38,    36,    42,     5,     6,    58,    17,     5,     6,    77,
      34,    37,   280,   281,     5,     6,    36,    34,    35,    36,
      60,    61,    62,    63,    64,    42,    38,    17,    38,    76,
     298,     3,    34,    36,    36,    37,    34,    17,    36,    37,
      37,    42,    37,    34,    16,    36,    37,   376,    20,    21,
      79,    37,   212,   163,   426,    27,    28,    29,   326,   327,
     195,   329,    34,    35,     2,    50,   323,    39,    40,    41,
      42,    43,    -1,    45,    -1,    47,    48,    -1,    50,    51,
      52,    53,    54,    55,    -1,    -1,    58,    59,    60,    61,
      62,    63,    64,    65,    66,    -1,     3,    16,    -1,    -1,
      -1,    -1,    74,    -1,    76,    77,    78,    -1,    -1,    16,
      -1,    -1,    -1,    20,    21,    34,    35,    36,    -1,    -1,
      27,    28,    29,    42,    -1,    -1,    -1,    34,    35,    -1,
      -1,   399,    39,    40,    41,    42,    43,    -1,    45,    -1,
      47,    48,    -1,    50,    51,    52,    53,    54,    55,    -1,
      -1,    58,    59,    60,    61,    62,    63,    64,    65,    66,
       3,    16,    60,    61,    62,    63,    64,    74,    -1,    76,
      77,    78,    -1,    16,    -1,    -1,    12,    20,    21,    34,
      35,    36,    -1,    -1,    27,    28,    29,    42,    -1,    -1,
      -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      43,    -1,    45,    -1,    47,    48,    42,    50,    51,    52,
      53,    54,    55,    -1,    -1,    58,    59,    60,    61,    62,
      63,    64,    65,    66,    60,    61,    62,    63,    64,    -1,
      -1,    74,    -1,    76,    77,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    16,    -1,    28,    29,    20,    21,    32,    -1,
      34,    35,    36,    27,    28,    29,    79,    -1,    -1,    -1,
      34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,
      -1,    45,    -1,    -1,    -1,    79,    -1,    -1,    -1,    12,
      -1,    -1,    -1,    -1,    58,    -1,    60,    61,    62,    63,
      64,    65,    -1,    -1,    79,    16,    -1,    -1,    -1,    20,
      21,    -1,    12,    -1,    -1,    79,    27,    28,    29,    42,
      -1,    -1,    -1,    34,    35,    -1,    -1,    38,    39,    40,
      41,    42,    -1,    -1,    45,    -1,    -1,    60,    61,    62,
      63,    64,    42,    -1,    -1,    -1,    -1,    58,    -1,    60,
      61,    62,    63,    64,    16,    -1,    -1,    -1,    20,    21,
      60,    61,    62,    63,    64,    27,    28,    29,    79,    42,
      -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
      42,    -1,    -1,    45,    -1,    58,    -1,    60,    61,    62,
      63,    64,    -1,    -1,    -1,    16,    58,    -1,    60,    61,
      62,    63,    64,    65,    -1,    -1,    -1,    16,    -1,    -1,
      -1,    20,    21,    -1,    35,    -1,    -1,    79,    27,    28,
      29,    42,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    -1,    -1,    45,    58,    -1,    60,
      61,    62,    63,    64,    -1,    16,    -1,    -1,    -1,    58,
      -1,    60,    61,    62,    63,    64,    16,    -1,    -1,    -1,
      20,    21,    -1,    -1,    35,    -1,    -1,    27,    28,    29,
      79,    42,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,
      40,    41,    42,    -1,    -1,    45,    -1,    58,    -1,    60,
      61,    62,    63,    64,    -1,    -1,    16,    -1,    58,    -1,
      60,    61,    62,    63,    64,    16,    -1,    -1,    -1,    20,
      21,    -1,    -1,    -1,    34,    35,    27,    28,    29,    39,
      40,    41,    42,    34,    35,    45,    -1,    -1,    39,    40,
      41,    42,    -1,    -1,    45,    -1,    -1,    -1,    58,    -1,
      60,    61,    62,    63,    64,    -1,    -1,    58,    -1,    60,
      61,    62,    63,    64,    16,    -1,    -1,    -1,    20,    21,
      60,    61,    62,    63,    64,    27,    28,    29,    -1,    -1,
      -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
      42,    -1,    -1,    45,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    -1,    58,    -1,    60,    61,
      62,    63,    64,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    28,    29,    -1,
      -1,    -1,    -1,    34,    35,    36,    -1,    38,    20,    21,
      22,    23,    24,    25,    26,    60,    61,    62,    63,    64,
      22,    23,    24,    25,    26
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    16,    20,    21,    27,    28,    29,    34,    35,
      39,    40,    41,    42,    43,    45,    47,    48,    50,    51,
      52,    53,    54,    55,    58,    59,    60,    61,    62,    63,
      64,    65,    66,    74,    76,    77,    82,    83,    84,    90,
      93,    94,    97,   100,   104,   106,   107,   110,   112,   117,
     118,   119,   120,   121,   126,   127,   128,   129,   130,    34,
      60,    61,    62,    63,    64,   126,    58,    60,    61,    62,
      63,    64,   117,   119,   125,    42,    98,    99,    45,    45,
      76,    76,   117,    16,    34,    35,    42,    58,    90,    34,
      36,    34,    36,    34,    36,    34,    36,    34,    36,    42,
       0,    46,   112,   101,    16,    35,    42,    91,    96,    76,
      84,    34,    35,    34,    77,     5,    76,   128,     6,    67,
      68,    69,    70,    71,    72,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    79,    80,     5,    22,    31,    33,
      28,    29,    32,    34,    35,    36,   117,    17,    37,     5,
      38,     6,    34,     5,    76,    76,    58,    79,   121,   122,
     123,    58,     6,    34,    42,    42,    60,    61,    62,    63,
      64,    65,    87,   122,   124,    42,   124,    42,   124,    42,
     124,    42,   124,    42,    34,    35,   115,   116,    42,    88,
      58,    42,    58,    60,    61,    62,    63,    64,    86,    87,
       6,    35,     5,    78,   105,    90,   108,   117,    87,   111,
     117,    84,   117,   119,   119,   119,   119,   119,   119,   119,
     121,   121,   121,   121,   121,   121,   121,   121,   121,   121,
     121,   121,   121,   121,   121,   121,   121,   121,   121,   121,
     121,    60,    61,    62,    63,    64,    60,    61,    62,    63,
      64,   127,   127,   124,    38,   122,    42,    37,   119,    42,
      60,    61,    62,    63,    64,    89,   120,    89,    98,    17,
      79,     5,    37,    37,     5,    38,   117,    16,    35,    60,
      62,    63,    64,    85,    87,    90,    34,     6,     5,    37,
      37,    37,    37,    37,    42,   125,   112,     6,    34,     5,
      76,    17,    38,     5,    38,   120,   120,    35,    91,   109,
      76,    76,     5,    79,    37,    78,    79,    79,    37,     5,
      38,   112,     6,    42,    42,    42,    42,    42,     5,    37,
      35,   121,   124,    58,    35,    76,    35,    85,    12,    87,
      12,    87,    12,    87,    12,    87,     5,    37,     6,    12,
      87,    85,   120,    60,    61,    62,    63,    64,    65,    87,
     122,     6,    38,   119,    89,    88,    35,    42,    92,    95,
      87,     6,    35,    38,    88,   117,   117,    87,   130,   112,
      56,   114,   121,   121,   122,     4,   120,    89,    58,    79,
      38,    58,    85,    38,    87,    87,    87,    87,    16,    35,
      60,    62,    63,    64,    87,    90,   103,   119,    87,    37,
       6,   117,    76,    37,    86,    34,     5,   120,   120,    76,
      76,     5,    38,    34,    38,   112,     5,    38,   121,    35,
      38,    38,    35,    85,    12,    87,    12,    87,    12,    87,
      12,    87,     6,    12,    87,     6,   102,   120,     5,    38,
     120,    92,    38,   117,   117,    87,   113,    36,    58,    58,
      17,    85,    38,    87,    87,    87,    87,   119,    87,   119,
      77,   117,    34,    37,    37,    37,   112,    36,    38,    38,
      38,    17,    76,    84,    37,    76,   120,   112,   112,    36,
      17,    78,    42,    37,    37,    37,   112
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
#line 338 "lg.ypp"
    {
                        if(  ffapi::ff_justcompile) exit(0);
    // clean FH  mach 2014
		        const char *  magicffglut="#!ffglutdata4.1\n";// for complex and vector 3d plot
			//FFCS: divert stream to FFCS
                        if(ThePlotStream) ffapi::fwriteinit(magicffglut,strlen(magicffglut),1,ThePlotStream);

                        // <<sizestack_set>>
                        size_t sizestack = currentblock->size()+1024 ; //  before close

                        // <<close_final_block>>
                       // $1+=currentblock->close(currentblock);
                       (yyvsp[(1) - (2)].cinst).setclose(Block::snewclose(currentblock));// Sep 2016 FH
                        if(verbosity>2 || ( (mpirank==0) && verbosity)) cout << " sizestack + 1024 =" << sizestack << "  ( " << sizestack-1024 <<" )\n" ;
                        size_t lg0,lg1;
                        ffapi::ifchtmpdir(); // change  to tmp after compile FH sep 2015 ...
                        int NbPtr = ShowAlloc("init execution ",lg0); // number of un delele ptr
                        debugstack= new vector<pair<const E_Routine*,int> >;
                        size_t stu0=storageused(); // get Storage usage
			UnShowAlloc =0;// add FH for parallee
                        if(verbosity>2  || ( (mpirank==0) && verbosity) ) cout << endl;
                        {

                            // <<create_global_FF_stack>> calls [[file:../fflib/ffstack.hpp::newStack]]

                            Stack stack = newStack(sizestack);

                        double CPUcompile= CPUtime();
                        try {

                          // <<evaluate_parsed_FF_script>> calls [[file:../fflib/AFunction.hpp::CListOfInst::eval]]
                          (yyvsp[(1) - (2)].cinst).eval(stack);
                        }
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


                        // <<delete_global_FF_stack>>

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
#line 412 "lg.ypp"
    {(yyval.cinst) = (yyvsp[(1) - (1)].cexp);;}
    break;

  case 5:
#line 413 "lg.ypp"
    {(yyval.cinst) = ((yyvsp[(1) - (2)].cinst)+=(yyvsp[(2) - (2)].cexp));;}
    break;

  case 6:
#line 419 "lg.ypp"
    { (yyval.clist_id) = new ListOfId();;}
    break;

  case 7:
#line 420 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str)));;}
    break;

  case 8:
#line 421 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (3)].str),(yyvsp[(3) - (3)].cexp)));;}
    break;

  case 9:
#line 422 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,2> **>()));;}
    break;

  case 10:
#line 423 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,2> **>(),true));;}
    break;

  case 11:
#line 424 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,3> **>()));;}
    break;

  case 12:
#line 425 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,3> **>(),true));;}
    break;

  case 13:
#line 426 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,4> **>()));;}
    break;

  case 14:
#line 427 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,4> **>(),true));;}
    break;

  case 15:
#line 428 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,5> **>()));;}
    break;

  case 16:
#line 429 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,5> **>(),true));;}
    break;

  case 17:
#line 430 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),C_F0(),(yyvsp[(1) - (2)].type)->right()));;}
    break;

  case 18:
#line 431 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),C_F0(),(yyvsp[(1) - (3)].type),true));;}
    break;

  case 19:
#line 432 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (3)].clist_id)));;}
    break;

  case 20:
#line 433 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (5)].clist_id),true,true));;}
    break;

  case 21:
#line 434 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (4)].clist_id),true,false));;}
    break;

  case 22:
#line 435 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (3)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str)));;}
    break;

  case 23:
#line 436 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (5)].clist_id)));;}
    break;

  case 24:
#line 437 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (6)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (6)].clist_id),false,true));;}
    break;

  case 25:
#line 438 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (7)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (7)].clist_id),true,true));;}
    break;

  case 26:
#line 439 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (6)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (6)].clist_id),true,false));;}
    break;

  case 27:
#line 440 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (5)].str),(yyvsp[(5) - (5)].cexp)));;}
    break;

  case 28:
#line 441 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,2> **>()));;}
    break;

  case 29:
#line 442 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,2> **>(),true));;}
    break;

  case 30:
#line 443 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,3> **>()));;}
    break;

  case 31:
#line 444 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,3> **>(),true));;}
    break;

  case 32:
#line 445 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,4> **>()));;}
    break;

  case 33:
#line 446 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,4> **>(),true));;}
    break;

  case 34:
#line 447 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,5> **>()));;}
    break;

  case 35:
#line 448 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,5> **>(),true));;}
    break;

  case 36:
#line 449 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),C_F0(),(yyvsp[(3) - (4)].type)->right()));;}
    break;

  case 37:
#line 450 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),C_F0(),(yyvsp[(3) - (5)].type),true));;}
    break;

  case 38:
#line 453 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str)));;}
    break;

  case 39:
#line 454 "lg.ypp"
    { (yyval.clist_id)=(yyvsp[(1) - (3)].clist_id)  ; (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str)));;}
    break;

  case 46:
#line 462 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (1)].str),dcltype);;}
    break;

  case 47:
#line 463 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (3)].str),dcltype,(yyvsp[(3) - (3)].cexp));;}
    break;

  case 48:
#line 464 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (4)].str),dcltype,(yyvsp[(3) - (4)].args));(yyvsp[(3) - (4)].args).destroy();;}
    break;

  case 49:
#line 465 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 50:
#line 471 "lg.ypp"
    { (yyval.args)=(yyvsp[(1) - (1)].cexp);;}
    break;

  case 51:
#line 472 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (2)].str));;}
    break;

  case 52:
#line 473 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (2)].str));;}
    break;

  case 53:
#line 474 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (2)].str));;}
    break;

  case 54:
#line 475 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (2)].str));;}
    break;

  case 55:
#line 476 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (2)].str));;}
    break;

  case 56:
#line 477 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 57:
#line 478 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 58:
#line 479 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 59:
#line 480 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 60:
#line 481 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 61:
#line 482 "lg.ypp"
    { (yyval.args)=make_pair<const char *,const C_F0>((const char *) (yyvsp[(1) - (3)].str),(C_F0) (yyvsp[(3) - (3)].cexp));;}
    break;

  case 62:
#line 483 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].args));;}
    break;

  case 64:
#line 489 "lg.ypp"
    {(yyval.type)=TypeArray((yyvsp[(1) - (4)].type),(yyvsp[(3) - (4)].type));;}
    break;

  case 65:
#line 490 "lg.ypp"
    {(yyval.type)=TypeArray(TypeArray((yyvsp[(1) - (7)].type),(yyvsp[(3) - (7)].type)),(yyvsp[(6) - (7)].type));;}
    break;

  case 66:
#line 491 "lg.ypp"
    {(yyval.type)=TypeArray((yyvsp[(1) - (6)].type),(yyvsp[(3) - (6)].type),(yyvsp[(5) - (6)].type));;}
    break;

  case 67:
#line 492 "lg.ypp"
    {(yyval.type)=TypeArray(TypeArray((yyvsp[(1) - (9)].type),(yyvsp[(3) - (9)].type),(yyvsp[(5) - (9)].type)),(yyvsp[(8) - (9)].type));;}
    break;

  case 68:
#line 493 "lg.ypp"
    {(yyval.type)=TypeTemplate((yyvsp[(1) - (4)].type),(yyvsp[(3) - (4)].type));;}
    break;

  case 69:
#line 494 "lg.ypp"
    {(yyval.type)=TypeArray(TypeTemplate((yyvsp[(1) - (7)].type),(yyvsp[(3) - (7)].type)),(yyvsp[(6) - (7)].type));;}
    break;

  case 70:
#line 495 "lg.ypp"
    {(yyval.type)=TypeArray(TypeTemplate((yyvsp[(1) - (9)].type),(yyvsp[(3) - (9)].type)),(yyvsp[(6) - (9)].type),(yyvsp[(8) - (9)].type));;}
    break;

  case 71:
#line 499 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(1) - (1)].str),currentblock,fespacetype,fespacecomplex,fespacedim); ;}
    break;

  case 72:
#line 500 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(1) - (4)].str),currentblock,fespacetype,(yyvsp[(3) - (4)].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 73:
#line 501 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(1) - (3)].str),currentblock,fespacetype,(yyvsp[(3) - (3)].cexp),fespacecomplex,fespacedim);;}
    break;

  case 74:
#line 502 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(2) - (3)].clist_id),currentblock,fespacetype,fespacecomplex,fespacedim);;}
    break;

  case 75:
#line 503 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(2) - (6)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (6)].cexp),fespacecomplex,fespacedim);;}
    break;

  case 76:
#line 504 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(2) - (5)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (5)].cexp),fespacecomplex,fespacedim);;}
    break;

  case 77:
#line 508 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(1) - (4)].str),currentblock,fespacetype,(yyvsp[(3) - (4)].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 78:
#line 509 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(2) - (6)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (6)].cexp),fespacecomplex,fespacedim);;}
    break;

  case 79:
#line 514 "lg.ypp"
    { fespacedim=2;;}
    break;

  case 80:
#line 515 "lg.ypp"
    { fespacedim=1;;}
    break;

  case 81:
#line 516 "lg.ypp"
    { fespacedim=3;;}
    break;

  case 82:
#line 517 "lg.ypp"
    { fespacedim=4;;}
    break;

  case 83:
#line 518 "lg.ypp"
    { fespacedim=5;;}
    break;

  case 84:
#line 519 "lg.ypp"
    { fespacedim=6;;}
    break;

  case 85:
#line 520 "lg.ypp"
    { fespacedim=7;;}
    break;

  case 86:
#line 523 "lg.ypp"
    {fespacecomplex=false;  fespacetype = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 87:
#line 524 "lg.ypp"
    {
             if ((yyvsp[(3) - (4)].type) != typevarreal && (yyvsp[(3) - (4)].type) != typevarcomplex) lgerror (" type of finite element <real> or <complex>");
             fespacecomplex=((yyvsp[(3) - (4)].type)==typevarcomplex);
             fespacetype = Find((yyvsp[(1) - (4)].str));;}
    break;

  case 88:
#line 529 "lg.ypp"
    {  (yyval.cexp) = (yyvsp[(1) - (1)].cexp);;}
    break;

  case 89:
#line 530 "lg.ypp"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 90:
#line 532 "lg.ypp"
    {  (yyval.cexp) = (yyvsp[(1) - (1)].cexp);;}
    break;

  case 91:
#line 533 "lg.ypp"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 92:
#line 535 "lg.ypp"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[(2) - (2)].cexp);;}
    break;

  case 93:
#line 536 "lg.ypp"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[(5) - (5)].cexp);;}
    break;

  case 94:
#line 539 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,KN<size_t>>((yyvsp[(1) - (4)].str),typeFESpace((yyvsp[(3) - (4)].args)),(yyvsp[(3) - (4)].args),dimFESpaceImage((yyvsp[(3) - (4)].args)));
     (yyvsp[(3) - (4)].args).destroy(); ;}
    break;

  case 95:
#line 542 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,KN<size_t>>((yyvsp[(1) - (3)].str),typeFESpace((yyvsp[(3) - (3)].args)),(yyvsp[(3) - (3)].args),dimFESpaceImage((yyvsp[(3) - (3)].args)));
     (yyvsp[(3) - (3)].args).destroy(); ;}
    break;

  case 97:
#line 547 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 98:
#line 550 "lg.ypp"
    {dcltype=(yyvsp[(1) - (1)].type);;}
    break;

  case 99:
#line 550 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(3) - (4)].cexp);;}
    break;

  case 100:
#line 551 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(2) - (3)].cexp);;}
    break;

  case 101:
#line 552 "lg.ypp"
    { (yyval.cexp)=(yyvsp[(1) - (2)].cexp);;}
    break;

  case 102:
#line 553 "lg.ypp"
    {(yyval.cexp)=currentblock->NewID((yyvsp[(1) - (5)].type),(yyvsp[(2) - (5)].str),(yyvsp[(4) - (5)].cexp));;}
    break;

  case 103:
#line 555 "lg.ypp"
    {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = (yyvsp[(2) - (6)].type)->right();
                      routineinblock[kkembtype] = currentblock;
                      (yyvsp[(5) - (6)].routine)=new Routine((yyvsp[(1) - (6)].type),(yyvsp[(2) - (6)].type)->right(),(yyvsp[(3) - (6)].str),(yyvsp[(5) - (6)].clist_id),currentblock);
		      // routineinblock[kkembtype]->Add($3,"(",$<routine>5); //pas recursif pour l'instanat test  FH 27 dec 2008
                     // cout << " \n after new routine \n " << endl;
                      ;}
    break;

  case 104:
#line 564 "lg.ypp"
    { currentblock=(yyvsp[(5) - (10)].routine)->Set((yyvsp[(9) - (10)].cinst));
                       currentblock->Add((yyvsp[(3) - (10)].str),"(",(yyvsp[(5) - (10)].routine)); //pas recursif pour l'instant test  FH 27 dec 2008
                       kkembtype--;
                       (yyval.cexp)=0;

                        ;}
    break;

  case 105:
#line 571 "lg.ypp"
    {Block::open(currentblock); (yyvsp[(1) - (5)].type)->SetArgs((yyvsp[(4) - (5)].clist_id));;}
    break;

  case 106:
#line 573 "lg.ypp"
    {  //$<cinst>$=currentblock->close(currentblock);
                         (yyval.cinst).setclose(Block::snewclose(currentblock));// Sep 2016 FH.
                         (yyval.cexp)=currentblock->NewID((yyvsp[(1) - (9)].type),(yyvsp[(2) - (9)].str),(yyvsp[(8) - (9)].cexp),*(yyvsp[(4) - (9)].clist_id));
                         delete (yyvsp[(4) - (9)].clist_id); //  FH 23032005
                         ;}
    break;

  case 107:
#line 580 "lg.ypp"
    {  Block::open(currentblock);;}
    break;

  case 108:
#line 581 "lg.ypp"
    { (yyval.endb)=Block::snewclose(currentblock);
//  $$=currentblock->close(currentblock);
;}
    break;

  case 109:
#line 585 "lg.ypp"
    {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 110:
#line 587 "lg.ypp"
    {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 111:
#line 591 "lg.ypp"
    {dcltype=(yyvsp[(1) - (1)].type); Block::open(currentblock);  ;}
    break;

  case 112:
#line 592 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(3) - (3)].cexp);;}
    break;

  case 113:
#line 594 "lg.ypp"
    { Block::open(currentblock);;}
    break;

  case 114:
#line 596 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str)));Block::open(currentblock); ;}
    break;

  case 115:
#line 597 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (3)].str)));(yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str)));Block::open(currentblock); ;}
    break;

  case 116:
#line 598 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (5)].str)));(yyval.clist_id)->push_back(UnId((yyvsp[(3) - (5)].str)));(yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str)));Block::open(currentblock); ;}
    break;

  case 117:
#line 600 "lg.ypp"
    {(yyval.cexp)=0;;}
    break;

  case 118:
#line 601 "lg.ypp"
    {zzzfff->input((yyvsp[(2) - (2)].str));(yyval.cexp)= 0; ;}
    break;

  case 119:
#line 602 "lg.ypp"
    {load((yyvsp[(2) - (2)].str));(yyval.cexp)= 0; ;}
    break;

  case 120:
#line 603 "lg.ypp"
    {(yyval.cexp)=Try((yyvsp[(3) - (5)].cinst),currentblock->close(currentblock,(yyvsp[(5) - (5)].cexp)));;}
    break;

  case 121:
#line 604 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(1) - (2)].cexp);;}
    break;

  case 122:
#line 605 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(1) - (1)].cexp);;}
    break;

  case 123:
#line 607 "lg.ypp"
    {(yyvsp[(5) - (6)].cexp)=ForAll(currentblock,(yyvsp[(3) - (6)].clist_id),(yyvsp[(5) - (6)].cexp));;}
    break;

  case 124:
#line 608 "lg.ypp"
    {
                    inloopcount--;
                    (yyval.cexp)=Block::close(currentblock,C_F0(ForAll((yyvsp[(5) - (8)].cexp),(yyvsp[(8) - (8)].cexp))));
                 ;}
    break;

  case 125:
#line 612 "lg.ypp"
    {inloopcount--; (yyval.cexp)=For((yyvsp[(3) - (9)].cexp),(yyvsp[(5) - (9)].cexp),(yyvsp[(7) - (9)].cexp),(yyvsp[(9) - (9)].cexp));;}
    break;

  case 126:
#line 614 "lg.ypp"
    {inloopcount--;
                 (yyval.cexp)=Block::close(currentblock,C_F0(For((yyvsp[(3) - (9)].cexp),(yyvsp[(5) - (9)].cexp),(yyvsp[(7) - (9)].cexp),(yyvsp[(9) - (9)].cexp))));
                ;}
    break;

  case 127:
#line 617 "lg.ypp"
    {inloopcount--;(yyval.cexp)=While((yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 128:
#line 618 "lg.ypp"
    {(yyval.cexp)=FIf((yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 129:
#line 619 "lg.ypp"
    {(yyval.cexp)=FIf((yyvsp[(3) - (7)].cexp),(yyvsp[(5) - (7)].cexp),(yyvsp[(7) - (7)].cexp));;}
    break;

  case 130:
#line 620 "lg.ypp"
    { /* [[begin:]] [[end:]] */
             (yyvsp[(2) - (3)].cinst).setclose((yyvsp[(3) - (3)].endb));
             (yyval.cexp)=(yyvsp[(2) - (3)].cinst);
                    //  $$=C_F0(new E_block($2,$3),atype<void>());
         ;}
    break;

  case 131:
#line 625 "lg.ypp"
    { /* <<BORDER_ID>> */
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[(2) - (3)].str),C_F0(TheOperators,"[border]",(yyvsp[(3) - (3)].args)));;}
    break;

  case 132:
#line 627 "lg.ypp"
    {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[(2) - (6)].str),C_F0(TheOperators,"[border]",(yyvsp[(4) - (6)].args)));;}
    break;

  case 133:
#line 630 "lg.ypp"
    {
                    if(inloopcount)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_break),atype<void>());
                    else lgerror("break not in loop");;}
    break;

  case 134:
#line 634 "lg.ypp"
    {
                    if(inloopcount)
                        (yyval.cexp)= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");;}
    break;

  case 135:
#line 638 "lg.ypp"
    {
                    if (kkembtype>=0)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_return,(rettype[kkembtype]->CastTo((yyvsp[(2) - (3)].cexp))).OnReturn()) ,atype<void>());
                     else lgerror(" return not in routine ");;}
    break;

  case 136:
#line 645 "lg.ypp"
    {(yyval.cexp) =  (yyvsp[(7) - (7)].cexp); ;}
    break;

  case 137:
#line 648 "lg.ypp"
    {
   Block::open(currentblock);
   (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[(2) - (7)].str),atype<double*>());
   (yyval.args)+= (yyvsp[(4) - (7)].cexp);
   (yyval.args)+= (yyvsp[(6) - (7)].cexp);
   (yyval.args)+= currentblock->NewVar<LocalVariable>("IndexBorder",atype<long*>());;}
    break;

  case 138:
#line 656 "lg.ypp"
    {
    Block::open(currentblock);
    (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[(2) - (9)].str),atype<double*>());
    (yyval.args)+= (yyvsp[(4) - (9)].cexp);
    (yyval.args)+= (yyvsp[(6) - (9)].cexp);
    (yyval.args)+= currentblock->NewVar<LocalVariable>((yyvsp[(8) - (9)].str),atype<long*>());;}
    break;

  case 139:
#line 664 "lg.ypp"
    {
    //currentblock->close(currentblock;);
   (yyval.args) = ((yyvsp[(1) - (2)].args) += currentblock->close(currentblock,(yyvsp[(2) - (2)].cexp)));
   ;}
    break;

  case 141:
#line 672 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 148:
#line 685 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 149:
#line 686 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"+=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 150:
#line 687 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"-=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 151:
#line 688 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"*=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 152:
#line 689 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"/=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 153:
#line 690 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,".*=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 154:
#line 691 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"./=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 156:
#line 697 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"?:",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 157:
#line 698 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 158:
#line 699 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 160:
#line 704 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 161:
#line 705 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 162:
#line 706 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 163:
#line 707 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 164:
#line 708 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 165:
#line 709 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 166:
#line 710 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 167:
#line 711 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 168:
#line 712 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 169:
#line 713 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 170:
#line 714 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 171:
#line 715 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 172:
#line 716 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 173:
#line 717 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 174:
#line 718 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 175:
#line 719 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 176:
#line 720 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 177:
#line 721 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 178:
#line 722 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 179:
#line 726 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(1) - (1)].cexp);;}
    break;

  case 180:
#line 727 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,":");;}
    break;

  case 181:
#line 728 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 182:
#line 729 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 183:
#line 734 "lg.ypp"
    {
      (yyval.args) = (yyvsp[(1) - (3)].cexp);
      (yyval.args) += (yyvsp[(3) - (3)].args); ;}
    break;

  case 184:
#line 738 "lg.ypp"
    {(yyval.args) = 0;;}
    break;

  case 185:
#line 739 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 186:
#line 740 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 187:
#line 741 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 188:
#line 742 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 189:
#line 743 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 190:
#line 744 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 191:
#line 745 "lg.ypp"
    {(yyval.args) = make_pair<const char *,const C_F0>((const char *) (yyvsp[(1) - (3)].str),(C_F0) (yyvsp[(3) - (3)].cexp));;}
    break;

  case 192:
#line 746 "lg.ypp"
    {(yyval.args) = (yyvsp[(1) - (1)].cexp);;}
    break;

  case 193:
#line 747 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 194:
#line 748 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 195:
#line 749 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 196:
#line 750 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 197:
#line 751 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 198:
#line 752 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 199:
#line 753 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp));;}
    break;

  case 200:
#line 756 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (5)].args)+= make_pair<const char *,const C_F0>((const char *)(yyvsp[(3) - (5)].str),(C_F0) (yyvsp[(5) - (5)].cexp)));;}
    break;

  case 201:
#line 759 "lg.ypp"
    {(yyval.args)=(yyvsp[(1) - (1)].cexp);;}
    break;

  case 202:
#line 760 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp));;}
    break;

  case 203:
#line 763 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));}
    break;

  case 204:
#line 764 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 205:
#line 765 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 206:
#line 766 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 207:
#line 767 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 208:
#line 768 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 209:
#line 769 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 210:
#line 770 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 211:
#line 771 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 212:
#line 772 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 213:
#line 773 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 214:
#line 774 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 215:
#line 775 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 216:
#line 776 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 217:
#line 777 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 219:
#line 782 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(1) - (2)].oper),(yyvsp[(2) - (2)].cexp));;}
    break;

  case 221:
#line 787 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 222:
#line 788 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 224:
#line 793 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));;}
    break;

  case 225:
#line 801 "lg.ypp"
    {(yyval.cexp)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 226:
#line 805 "lg.ypp"
    {(yyval.cexp)= CConstant((yyvsp[(1) - (1)].lnum));;}
    break;

  case 227:
#line 806 "lg.ypp"
    {(yyval.cexp)= CConstant((yyvsp[(1) - (1)].dnum));;}
    break;

  case 228:
#line 807 "lg.ypp"
    {(yyval.cexp)= CConstant(complex<double>(0,(yyvsp[(1) - (1)].dnum)));;}
    break;

  case 229:
#line 808 "lg.ypp"
    {(yyval.cexp)= CConstant<const char *>((yyvsp[(1) - (1)].str));;}
    break;

  case 230:
#line 813 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (4)].cexp),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args));;}
    break;

  case 231:
#line 815 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (4)].cexp),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].cexp));;}
    break;

  case 232:
#line 816 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (6)].cexp),(yyvsp[(2) - (6)].oper),(yyvsp[(3) - (6)].cexp),(yyvsp[(5) - (6)].cexp));;}
    break;

  case 233:
#line 817 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),"[]");;}
    break;

  case 234:
#line 818 "lg.ypp"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 235:
#line 819 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 236:
#line 820 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 237:
#line 821 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 238:
#line 822 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 239:
#line 823 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 240:
#line 824 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 241:
#line 825 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 242:
#line 826 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 243:
#line 827 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 244:
#line 828 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 245:
#line 829 "lg.ypp"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));;}
    break;

  case 246:
#line 830 "lg.ypp"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));;}
    break;

  case 247:
#line 831 "lg.ypp"
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

  case 248:
#line 840 "lg.ypp"
    {
           { (yyval.cexp)=(yyvsp[(1) - (4)].type)->right()->Find("<--",basicAC_F0_wa((yyvsp[(3) - (4)].args)));
           if (!(yyval.cexp).left()) { cerr << " no wait to change (args) in " <<
                                      (yyvsp[(1) - (4)].type)->right()->name() << endl;
                              CompileError(" Error in type(exp) "); }
           }
          ;}
    break;

  case 249:
#line 848 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(2) - (3)].cexp);;}
    break;

  case 250:
#line 849 "lg.ypp"
    { (yyval.cexp)=C_F0(TheOperators,"[]",(yyvsp[(2) - (3)].args));;}
    break;

  case 251:
#line 850 "lg.ypp"
    { (yyval.cexp)=C_F0(TheOperators,"<>",(yyvsp[(2) - (3)].args));;}
    break;

  case 252:
#line 851 "lg.ypp"
    { (yyval.cexp)=C_F0(TheOperators,"<>",(yyvsp[(1) - (1)].args));;}
    break;


/* Line 1267 of yacc.c.  */
#line 3535 "lg.tab.cpp"
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


#line 855 "lg.ypp"



#include <fstream>
using namespace std;
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

  // <<initialize_currentblock>>

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

/// <<mainff>> Called by [[file:mymain.cpp::mymain]] and calls [[Compile]] to run the FF language parser

int mainff (int  argc, char **argv)
{
#ifndef _WIN32
  	signal(SIGXCPU, signalCPUHandler);
#endif	
  if(argc)
    prognamearg=argv[0];

 //   int vvold=verbosity;
  if(mpirank !=0) verbosity=0;

  // ALH - 14/10/8 - This breaks FFCS output redirection
#ifndef ENABLE_FFCS
  SetcppIo();
#endif

  GetEnvironment(); // [[file:~/ff/src/fflib/environment.cpp::GetEnvironment]]
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
  if ( ! (getprog(cc,argc,argv) >0)  ) // [[file:~/ff/src/Graphics/getprog-unix.hpp::getprog]]
    {
      cout << "-- FreeFem++ v" << StrVersionNumber() << " (error parameter!)\n"  ;
      if(ThePlotStream) {ffapi::ff_pclose(ThePlotStream); ThePlotStream=0;}
      return 1;
    }

  if(verbosity && (mpirank==0)) {
      cout << "-- FreeFem++ v" << StrVersionNumber() << endl;
      cout << "   file : " << cc ;
      if(verbosity>1) cout << " " << " verbosity= " << verbosity ;
      cout  << endl;
  }

    KN<String> karg(argc);
    for(int i=0;i< argc;++i)
	karg[i]=argv[i];
    pkarg= &karg;

    /// <<zzzfff>>
    zzzfff = Newlex(cout,echo_edp,pkarg);


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
  zzzfff->input(cc); // [[file:../fflib/lex.cpp::mylex_input_filename]]
  EnvironmentLoad(); // just before compile [[file:~/ff/src/fflib/environment.cpp::EnvironmentLoad]]

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

