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
     PLUSEQ = 299,
     MOINSEQ = 300,
     MULEQ = 301,
     DIVEQ = 302,
     DOTMULEQ = 303,
     DOTDIVEQ = 304,
     ARROW = 305,
     BORDER = 306,
     CURVE = 307,
     SOLVE = 308
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
#define PLUSEQ 299
#define MOINSEQ 300
#define MULEQ 301
#define DIVEQ 302
#define DOTMULEQ 303
#define DOTDIVEQ 304
#define ARROW 305
#define BORDER 306
#define CURVE 307
#define SOLVE 308




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

#ifndef _WIN32
#include <csignal>

void signalCPUHandler( int signum ) {
    ff_finalize();
    std::cout << "Cputime limit exceeded:  (" << signum << ") received.\n";

    exit(24);
}
#endif


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
#line 163 "lg.ypp"
{
 double dnum;

 /* <<YYSTYPE_lnum>> */
 long lnum;

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
#line 401 "lg.tab.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 414 "lg.tab.cpp"

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
#define YYFINAL  86
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1130

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  79
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  48
/* YYNRULES -- Number of rules.  */
#define YYNRULES  208
/* YYNRULES -- Number of states.  */
#define YYNSTATES  452

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   308

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
       2,     2,     2,     2,     2,     2,     2,     2,    77,    74,
      16,     6,    17,    78,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    35,     2,    38,    31,    33,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    75,    10,    76,     2,     2,     2,     2,
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
      65,    66,    67,    68,    69,    70,    71,    72,    73
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     6,     8,    10,    13,    14,    16,    20,
      23,    27,    30,    34,    37,    41,    44,    48,    52,    56,
      62,    68,    73,    79,    84,    90,    95,   101,   106,   112,
     114,   118,   120,   122,   124,   126,   128,   130,   134,   139,
     143,   145,   148,   151,   154,   157,   161,   165,   171,   173,
     178,   186,   193,   203,   208,   216,   226,   228,   233,   237,
     241,   248,   254,   259,   266,   268,   270,   272,   274,   276,
     281,   283,   287,   289,   293,   296,   302,   307,   309,   313,
     314,   319,   323,   326,   332,   333,   344,   345,   355,   357,
     359,   361,   363,   364,   368,   370,   372,   376,   382,   384,
     387,   390,   396,   399,   401,   402,   411,   421,   431,   437,
     443,   451,   455,   459,   466,   469,   472,   476,   484,   492,
     502,   505,   507,   511,   513,   515,   517,   519,   521,   523,
     527,   531,   535,   539,   543,   547,   551,   553,   559,   563,
     569,   571,   575,   579,   583,   587,   591,   595,   599,   603,
     607,   611,   615,   619,   623,   627,   631,   635,   639,   643,
     647,   649,   651,   655,   661,   662,   664,   666,   668,   670,
     674,   676,   680,   684,   688,   692,   696,   702,   704,   708,
     710,   713,   715,   719,   723,   725,   728,   730,   732,   734,
     736,   738,   743,   748,   755,   759,   763,   767,   772,   776,
     781,   785,   790,   794,   799,   802,   805,   810,   814
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      80,     0,    -1,    81,    46,    -1,    82,    -1,   110,    -1,
      82,   110,    -1,    -1,    85,    -1,    85,     6,   117,    -1,
      60,    85,    -1,    60,    12,    85,    -1,    62,    85,    -1,
      62,    12,    85,    -1,    63,    85,    -1,    63,    12,    85,
      -1,    88,    85,    -1,    88,    12,    85,    -1,    35,    83,
      38,    -1,    83,     5,    85,    -1,    83,     5,    35,    83,
      38,    -1,    83,     5,    85,     6,   117,    -1,    83,     5,
      60,    85,    -1,    83,     5,    60,    12,    85,    -1,    83,
       5,    62,    85,    -1,    83,     5,    62,    12,    85,    -1,
      83,     5,    63,    85,    -1,    83,     5,    63,    12,    85,
      -1,    83,     5,    88,    85,    -1,    83,     5,    88,    12,
      85,    -1,    85,    -1,    84,     5,    85,    -1,    42,    -1,
      60,    -1,    62,    -1,    63,    -1,    61,    -1,    42,    -1,
      42,     6,   117,    -1,    42,    34,    87,    37,    -1,    86,
       5,    86,    -1,   118,    -1,    60,    42,    -1,    61,    42,
      -1,    62,    42,    -1,    63,    42,    -1,    42,     6,   118,
      -1,    87,     5,   118,    -1,    87,     5,    85,     6,   118,
      -1,    58,    -1,    58,    35,    58,    38,    -1,    58,    35,
      58,    38,    35,    58,    38,    -1,    58,    35,    58,     5,
      58,    38,    -1,    58,    35,    58,     5,    58,    38,    35,
      58,    38,    -1,    58,    16,    58,    17,    -1,    58,    16,
      58,    17,    35,    58,    38,    -1,    58,    16,    58,    17,
      35,    58,     5,    58,    38,    -1,    42,    -1,    42,    35,
     118,    38,    -1,    42,     6,   118,    -1,    35,    84,    38,
      -1,    35,    84,    38,    35,   118,    38,    -1,    35,    84,
      38,     6,   118,    -1,    42,    34,   118,    37,    -1,    35,
      84,    38,    34,   118,    37,    -1,    60,    -1,    61,    -1,
      62,    -1,    63,    -1,    91,    -1,    91,    16,    58,    17,
      -1,    90,    -1,    93,     5,    90,    -1,    89,    -1,    94,
       5,    89,    -1,    92,    94,    -1,    92,    35,    58,    38,
      93,    -1,    42,    34,    87,    37,    -1,    96,    -1,    97,
       5,    96,    -1,    -1,    88,    99,    86,    74,    -1,    43,
      97,    74,    -1,    95,    74,    -1,    59,    42,     6,   115,
      74,    -1,    -1,    59,    88,    42,    34,    83,    37,   100,
      75,    82,    76,    -1,    -1,    59,    42,    34,    83,    37,
     101,     6,   117,    74,    -1,    75,    -1,    76,    -1,    50,
      -1,    51,    -1,    -1,    88,   107,    86,    -1,    55,    -1,
      85,    -1,    85,     5,    85,    -1,    85,     5,    85,     5,
      85,    -1,    74,    -1,    47,    45,    -1,    48,    45,    -1,
     108,    75,    82,    76,   112,    -1,   115,    74,    -1,    98,
      -1,    -1,   104,    35,   109,    77,   126,    38,   111,   110,
      -1,   104,    34,   115,    74,   115,    74,   115,    37,   110,
      -1,   104,    34,   106,    74,   115,    74,   115,    37,   110,
      -1,   105,    34,   115,    37,   110,    -1,     3,    34,   115,
      37,   110,    -1,     3,    34,   115,    37,   110,     4,   110,
      -1,   102,    82,   103,    -1,    71,    42,   114,    -1,    71,
      42,    35,   122,    38,    74,    -1,    52,    74,    -1,    53,
      74,    -1,    54,   115,    74,    -1,    56,    34,    36,    36,
      36,    37,   110,    -1,    34,    42,     6,   115,     5,   115,
      37,    -1,    34,    42,     6,   115,     5,   115,    74,    42,
      37,    -1,   113,   110,    -1,   117,    -1,   115,     5,   115,
      -1,    21,    -1,    20,    -1,    27,    -1,    29,    -1,    28,
      -1,   118,    -1,   118,     6,   117,    -1,   118,    64,   117,
      -1,   118,    65,   117,    -1,   118,    66,   117,    -1,   118,
      67,   117,    -1,   118,    68,   117,    -1,   118,    69,   117,
      -1,   119,    -1,   119,    78,   119,    77,   119,    -1,   119,
      77,   119,    -1,   119,    77,   119,    77,   119,    -1,   123,
      -1,   119,    22,   119,    -1,   119,    26,   119,    -1,   119,
      25,   119,    -1,   119,    23,   119,    -1,   119,    24,   119,
      -1,   119,    20,   119,    -1,   119,    21,   119,    -1,   119,
       9,   119,    -1,   119,     8,   119,    -1,   119,    12,   119,
      -1,   119,    13,   119,    -1,   119,    10,   119,    -1,   119,
      11,   119,    -1,   119,    16,   119,    -1,   119,    19,   119,
      -1,   119,    17,   119,    -1,   119,    18,   119,    -1,   119,
      15,   119,    -1,   119,    14,   119,    -1,   119,    -1,    77,
      -1,   119,    77,   119,    -1,   119,    77,   119,    77,   119,
      -1,    -1,    60,    -1,    61,    -1,    62,    -1,    63,    -1,
      85,     6,   118,    -1,   120,    -1,   121,     5,    60,    -1,
     121,     5,    61,    -1,   121,     5,    62,    -1,   121,     5,
      63,    -1,   121,     5,   120,    -1,   121,     5,    85,     6,
     118,    -1,   117,    -1,   122,     5,   117,    -1,   124,    -1,
     116,   124,    -1,   125,    -1,   125,    31,   123,    -1,   125,
      33,   123,    -1,   126,    -1,   126,    32,    -1,    42,    -1,
      39,    -1,    40,    -1,    41,    -1,    45,    -1,   126,    34,
     121,    37,    -1,   126,    35,   120,    38,    -1,   126,    35,
     120,     5,   120,    38,    -1,   126,    35,    38,    -1,   126,
      36,    42,    -1,    60,    36,    42,    -1,    60,    34,   121,
      37,    -1,    61,    36,    42,    -1,    61,    34,   121,    37,
      -1,    62,    36,    42,    -1,    62,    34,   121,    37,    -1,
      63,    36,    42,    -1,    63,    34,   121,    37,    -1,   126,
      29,    -1,   126,    28,    -1,    58,    34,   115,    37,    -1,
      34,   115,    37,    -1,    35,   122,    38,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   329,   329,   399,   403,   404,   410,   411,   412,   413,
     414,   415,   416,   417,   418,   419,   420,   421,   422,   423,
     424,   425,   426,   427,   428,   429,   430,   431,   432,   435,
     436,   441,   441,   441,   441,   441,   444,   445,   446,   447,
     453,   454,   455,   456,   457,   458,   459,   460,   463,   464,
     465,   466,   467,   468,   469,   470,   474,   475,   476,   477,
     478,   479,   483,   484,   489,   490,   491,   492,   495,   496,
     501,   502,   504,   505,   507,   508,   512,   515,   516,   519,
     519,   520,   521,   522,   524,   523,   540,   539,   549,   550,
     554,   556,   560,   560,   563,   565,   566,   567,   569,   570,
     571,   572,   573,   574,   576,   575,   581,   582,   586,   587,
     588,   589,   594,   596,   599,   603,   607,   614,   617,   625,
     633,   640,   641,   645,   646,   647,   648,   649,   653,   654,
     655,   656,   657,   658,   659,   660,   665,   666,   667,   668,
     672,   673,   674,   675,   676,   677,   678,   679,   680,   681,
     682,   683,   684,   685,   686,   687,   688,   689,   690,   691,
     695,   696,   697,   698,   703,   704,   705,   706,   707,   708,
     709,   710,   711,   712,   713,   714,   717,   720,   721,   725,
     726,   730,   731,   732,   736,   737,   745,   749,   750,   751,
     752,   757,   759,   760,   761,   762,   763,   764,   765,   766,
     767,   768,   769,   770,   771,   772,   773,   782,   783
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
  "FUNCTION", "FESPACE", "FESPACE1", "FESPACE3", "FESPACES", "PLUSEQ",
  "MOINSEQ", "MULEQ", "DIVEQ", "DOTMULEQ", "DOTDIVEQ", "ARROW", "BORDER",
  "CURVE", "SOLVE", "';'", "'{'", "'}'", "':'", "'?'", "$accept", "start",
  "input", "instructions", "list_of_id_args", "list_of_id1", "id",
  "list_of_dcls", "parameters_list", "type_of_dcl", "ID_space",
  "ID_array_space", "fespace123", "fespace", "spaceIDa", "spaceIDb",
  "spaceIDs", "fespace_def", "fespace_def_list", "declaration", "@1", "@2",
  "@3", "begin", "end", "for_loop", "while_loop", "declaration_for", "@4",
  "try", "IDfor", "instruction", "@5", "catchs", "bornes", "border_expr",
  "Expr", "unop", "no_comma_expr", "no_set_expr", "no_ternary_expr",
  "sub_script_expr", "parameters", "array", "unary_expr", "pow_expr",
  "primaryp", "primary", 0
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
     305,   306,   307,   308,    59,   123,   125,    58,    63
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    79,    80,    81,    82,    82,    83,    83,    83,    83,
      83,    83,    83,    83,    83,    83,    83,    83,    83,    83,
      83,    83,    83,    83,    83,    83,    83,    83,    83,    84,
      84,    85,    85,    85,    85,    85,    86,    86,    86,    86,
      87,    87,    87,    87,    87,    87,    87,    87,    88,    88,
      88,    88,    88,    88,    88,    88,    89,    89,    89,    89,
      89,    89,    90,    90,    91,    91,    91,    91,    92,    92,
      93,    93,    94,    94,    95,    95,    96,    97,    97,    99,
      98,    98,    98,    98,   100,    98,   101,    98,   102,   103,
     104,   105,   107,   106,   108,   109,   109,   109,   110,   110,
     110,   110,   110,   110,   111,   110,   110,   110,   110,   110,
     110,   110,   110,   110,   110,   110,   110,   112,   113,   113,
     114,   115,   115,   116,   116,   116,   116,   116,   117,   117,
     117,   117,   117,   117,   117,   117,   118,   118,   118,   118,
     119,   119,   119,   119,   119,   119,   119,   119,   119,   119,
     119,   119,   119,   119,   119,   119,   119,   119,   119,   119,
     120,   120,   120,   120,   121,   121,   121,   121,   121,   121,
     121,   121,   121,   121,   121,   121,   121,   122,   122,   123,
     123,   124,   124,   124,   125,   125,   126,   126,   126,   126,
     126,   126,   126,   126,   126,   126,   126,   126,   126,   126,
     126,   126,   126,   126,   126,   126,   126,   126,   126
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     1,     1,     2,     0,     1,     3,     2,
       3,     2,     3,     2,     3,     2,     3,     3,     3,     5,
       5,     4,     5,     4,     5,     4,     5,     4,     5,     1,
       3,     1,     1,     1,     1,     1,     1,     3,     4,     3,
       1,     2,     2,     2,     2,     3,     3,     5,     1,     4,
       7,     6,     9,     4,     7,     9,     1,     4,     3,     3,
       6,     5,     4,     6,     1,     1,     1,     1,     1,     4,
       1,     3,     1,     3,     2,     5,     4,     1,     3,     0,
       4,     3,     2,     5,     0,    10,     0,     9,     1,     1,
       1,     1,     0,     3,     1,     1,     3,     5,     1,     2,
       2,     5,     2,     1,     0,     8,     9,     9,     5,     5,
       7,     3,     3,     6,     2,     2,     3,     7,     7,     9,
       2,     1,     3,     1,     1,     1,     1,     1,     1,     3,
       3,     3,     3,     3,     3,     3,     1,     5,     3,     5,
       1,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       1,     1,     3,     5,     0,     1,     1,     1,     1,     3,
       1,     3,     3,     3,     3,     3,     5,     1,     3,     1,
       2,     1,     3,     3,     1,     2,     1,     1,     1,     1,
       1,     4,     4,     6,     3,     3,     3,     4,     3,     4,
       3,     4,     3,     4,     2,     2,     4,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,   124,   123,   125,   127,   126,     0,     0,   187,
     188,   189,   186,     0,   190,     0,     0,    90,    91,     0,
       0,     0,    94,    48,     0,    64,    65,    66,    67,     0,
      98,    88,     0,     0,     3,    79,    68,     0,     0,   103,
       0,     0,     0,     0,     4,     0,     0,   121,   128,   136,
     140,   179,   181,   184,     0,     0,     0,     0,     0,     0,
       0,   177,     0,     0,    77,     0,    99,   100,   114,   115,
       0,     0,     0,     0,     0,    48,     0,   164,     0,   164,
       0,   164,     0,   164,     0,     0,     1,     2,     5,     0,
       0,     0,    56,    72,    74,    82,     0,     0,     0,     0,
       0,     0,   102,   180,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   205,   204,   185,   164,     0,     0,
       0,   207,     0,   208,     0,     0,    81,   116,     0,     0,
       0,     0,     6,     0,   186,   165,   166,   167,   168,   161,
       0,   160,   170,     0,   196,     0,   198,     0,   200,     0,
     202,     0,     0,     0,   112,    36,     0,     0,    31,     0,
      32,    35,    33,    34,     0,    29,     0,     0,     0,    89,
     111,    92,     0,     0,    95,     0,     0,     0,   122,   129,
     130,   131,   132,   133,   134,   135,   149,   148,   152,   153,
     150,   151,   159,   158,   154,   156,   157,   155,   146,   147,
     141,   144,   145,   143,   142,   138,     0,   182,   183,     0,
     194,     0,   195,     0,   178,   186,     0,     0,     0,     0,
       0,    40,    78,    53,   206,     0,    49,     0,     6,    32,
      33,    34,     0,     7,     0,     6,     0,     0,     0,   197,
     199,   201,   203,     0,     0,   120,     0,     0,     0,    80,
      69,     0,     0,    59,    58,     0,     0,    73,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   191,     0,   192,
     109,     0,    41,    42,    43,    44,     0,    76,     0,     0,
       0,    83,     0,     0,     9,     0,    11,     0,    13,     0,
      86,     0,     0,    15,     0,   169,   162,   171,   172,   173,
     174,     0,   175,     0,     0,    37,     0,    39,     0,     0,
      70,    75,    30,     0,     0,    57,    93,     0,     0,    96,
       0,   108,     0,   101,   139,   137,     0,     0,    45,    32,
      35,    33,    34,     0,    46,     0,    51,     0,    17,    10,
      12,    14,     6,    32,    33,    34,    18,     0,     0,     8,
      16,    84,     0,     0,     0,   113,    38,     0,     0,     0,
      61,     0,     0,     0,     0,   104,     0,   193,   110,     0,
       0,    54,     0,    50,     0,     0,    21,     0,    23,     0,
      25,     0,     0,    27,     0,     0,   163,   176,     0,     0,
       0,    71,    60,     0,     0,    97,     0,     0,    47,     0,
       0,    19,    22,    24,    26,    20,    28,     0,     0,   122,
       0,    62,     0,     0,   105,     0,    55,    52,    87,     0,
     118,     0,     0,   107,   106,     0,    85,     0,    63,     0,
     119,   117
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    32,    33,    34,   252,   184,   160,   176,   240,    35,
      93,   330,    36,    37,   331,    94,    38,    64,    65,    39,
      89,   405,   368,    40,   190,    41,    42,   192,   278,    43,
     195,    44,   416,   343,   173,   174,    45,    46,    47,    48,
      49,   162,   163,    62,    50,    51,    52,    53
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -232
static const yytype_int16 yypact[] =
{
     581,   -15,  -232,  -232,  -232,  -232,  -232,   694,   694,  -232,
    -232,  -232,  -232,    -5,  -232,    49,    52,  -232,  -232,   -32,
     -24,   694,  -232,   162,    -9,   186,   208,   227,   242,    30,
    -232,  -232,    59,    61,   581,  -232,   120,    57,    80,  -232,
     581,   211,   123,    86,  -232,     3,   950,  -232,    82,   649,
    -232,  -232,   278,   388,   694,   133,   186,   208,   227,   242,
      33,  -232,     2,   137,  -232,     5,  -232,  -232,  -232,  -232,
       7,   125,   694,   163,   171,   188,   166,   783,   215,   783,
     226,   783,   237,   783,   238,   230,  -232,  -232,  -232,   246,
     245,   312,    17,  -232,   295,  -232,   398,   915,   300,   694,
     581,   694,  -232,  -232,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,  -232,  -232,  -232,   783,   190,   272,
      34,  -232,   694,  -232,   959,    -5,  -232,  -232,   302,    36,
       6,   694,   283,   290,   321,    62,   203,   241,   263,  -232,
     323,   668,  -232,   127,  -232,   132,  -232,   136,  -232,   138,
    -232,   297,   694,   581,  -232,   250,     8,   315,  -232,   317,
    -232,  -232,  -232,  -232,     9,  -232,   694,   694,    98,  -232,
    -232,  -232,   279,    11,   366,   306,   139,   459,  -232,  -232,
    -232,  -232,  -232,  -232,  -232,  -232,  1076,  1076,  1091,  1091,
    1104,  1104,   857,   857,   625,   625,   625,   625,   342,   342,
    -232,  -232,  -232,  -232,  -232,   687,   757,  -232,  -232,   147,
    -232,    20,  -232,   581,  -232,   373,   224,   270,   344,   394,
     150,  -232,  -232,   349,  -232,   327,   354,    13,   283,    -7,
     130,   153,   151,   384,   229,   283,   694,   694,   827,  -232,
    -232,  -232,  -232,   387,    22,  -232,   694,   959,   246,  -232,
    -232,   118,   300,    51,  -232,   356,   300,  -232,   246,   694,
     694,   300,   950,   581,   339,   694,   694,  -232,   871,  -232,
     392,   694,  -232,  -232,  -232,  -232,  1003,  -232,   352,   374,
     353,  -232,    26,   300,  -232,   300,  -232,   300,  -232,   346,
    -232,   694,   300,  -232,   157,  -232,   776,    62,   203,   241,
     263,   407,  -232,   694,   340,  -232,   161,  -232,   300,   395,
    -232,   416,  -232,   694,   694,  -232,   426,    15,    16,   429,
     447,  -232,   401,  -232,  1059,  1059,   404,   581,  -232,   186,
     208,   227,   242,   441,  -232,    27,   419,   417,  -232,  -232,
    -232,  -232,   283,   260,   274,   289,   457,   296,   458,  -232,
    -232,  -232,   694,   694,   460,  -232,  -232,    29,   694,   118,
    -232,   428,   694,   694,   300,  -232,   431,  -232,  -232,   694,
     410,  -232,   412,  -232,    31,   300,  -232,   300,  -232,   300,
    -232,   694,   300,  -232,   694,   396,  1059,  -232,   694,   443,
     452,  -232,  -232,   165,   169,  -232,   581,   442,  -232,   446,
     453,  -232,  -232,  -232,  -232,  -232,  -232,   418,   581,   -31,
     694,  -232,   581,   581,  -232,   454,  -232,  -232,  -232,   520,
    -232,   455,   466,  -232,  -232,   468,  -232,   471,  -232,   581,
    -232,  -232
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -232,  -232,  -232,   -37,  -231,   167,   -69,  -217,   248,   -22,
     308,   145,  -232,  -232,  -232,  -232,  -232,   371,  -232,  -232,
    -232,  -232,  -232,  -232,  -232,  -232,  -232,  -232,  -232,  -232,
    -232,   -34,  -232,  -232,  -232,  -232,    -6,  -232,    -4,    14,
      -3,  -129,   -53,   355,   -87,   479,  -232,   244
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -36
static const yytype_int16 yytable[] =
{
      88,    60,    76,    96,    61,   303,   440,   142,   101,   231,
     145,   245,   101,   268,   272,    70,   101,   302,   101,    54,
     101,   101,   185,   186,   314,   288,   165,   142,   167,   194,
     169,   309,   390,    74,   272,   178,   309,    63,   101,   101,
     143,   101,    68,   441,   246,   227,   228,   273,   140,    75,
      69,   327,   187,   180,   181,   182,   183,   333,   289,    86,
     324,   336,    88,   197,   358,   391,   149,   409,   -32,   421,
     141,   233,    85,   244,   161,   191,   161,   102,   161,   146,
     161,   147,   269,   253,   229,   280,   334,   301,   104,   382,
     383,   193,    91,   196,    66,   198,    77,    67,    78,    92,
     199,   200,   201,   202,   203,   204,   205,    87,   206,   207,
     208,   209,   210,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   322,
     254,   394,   258,   276,   161,   161,    90,   258,   234,   265,
      92,   258,   305,   258,   101,   247,   105,   106,   107,   108,
     109,   110,   258,   328,    95,   296,   309,    99,   241,   346,
     329,   100,   309,    88,   259,   307,   296,    72,    61,   260,
     101,   144,   178,   261,   101,   262,   283,   151,    71,   253,
     304,   306,   308,   148,   287,   313,   253,   297,   310,   321,
     180,   181,   182,   183,   371,   178,    72,    73,   376,   290,
     274,   275,   432,   332,    71,   152,   433,   185,   153,   -35,
       2,     3,   339,   180,   181,   182,   183,     4,     5,     6,
      77,   150,    78,    73,     7,     8,   254,   353,   230,     9,
      10,    11,    12,   254,   359,    14,   360,    79,   361,    80,
     366,   312,    79,   370,    80,    97,    98,   -33,    55,   341,
      56,    57,    58,    59,   316,   161,   266,   164,    77,   185,
      78,    81,   325,    82,   171,   172,   292,   159,   166,   -34,
     315,   178,   395,   337,   338,    81,    83,    82,    84,   168,
     170,   241,   344,   345,   267,   161,   397,   367,   175,   180,
     181,   182,   183,   253,   396,   398,   400,    83,   403,    84,
     188,   399,   178,   177,    79,   348,    80,   369,   402,   132,
     354,   133,   293,   388,   232,   415,   178,   374,   248,   243,
     180,   181,   182,   183,   255,   178,   422,   -31,   423,   256,
     424,   178,   270,   426,   180,   181,   182,   183,   178,   263,
     254,    75,   178,   249,   181,   250,   251,   380,   381,   180,
     181,   182,   183,   279,   178,   271,   180,   181,   182,   183,
     180,   181,   182,   183,   125,   126,   127,   128,   129,   406,
     179,   281,   180,   181,   182,   183,   413,   414,    81,   291,
      82,   362,   434,   282,   298,   299,   294,   407,   178,   300,
     311,   439,   410,   323,   335,   342,   347,   425,   443,   444,
     427,     1,   429,   418,    75,    88,   363,   181,   364,   365,
     355,   357,   356,   373,   375,   451,   134,   135,     2,     3,
     136,   379,   137,   138,   139,     4,     5,     6,    83,   378,
      84,   268,     7,     8,   384,   386,   295,     9,    10,    11,
      12,    13,   387,    14,   442,    15,    16,   389,    17,    18,
      19,    20,    21,    22,   392,   393,    23,    24,    25,    26,
      27,    28,     1,   401,   404,   408,   412,   417,   419,    29,
     420,   428,    30,    31,   189,   134,   135,   430,   435,     2,
       3,   137,   138,   139,   436,   385,     4,     5,     6,   431,
     445,   437,   438,     7,     8,   377,   277,   447,     9,    10,
      11,    12,    13,   448,    14,   449,    15,    16,   450,    17,
      18,    19,    20,    21,    22,   326,   242,    23,    24,    25,
      26,    27,    28,     1,   411,   103,   340,   264,     0,     0,
      29,     0,     0,    30,    31,   284,     0,     0,     0,     0,
       2,     3,     0,     0,     0,     0,     0,     4,     5,     6,
       0,     0,     0,     0,     7,     8,     0,     0,     0,     9,
      10,    11,    12,    13,     0,    14,     0,    15,    16,     0,
      17,    18,    19,    20,    21,    22,     0,     0,    23,    24,
      25,    26,    27,    28,     1,     0,     0,     0,     0,     0,
       0,    29,     0,     0,    30,    31,   446,     0,     0,     0,
       0,     2,     3,     0,     0,     0,     0,     0,     4,     5,
       6,     0,     0,     0,     0,     7,     8,     0,     0,     0,
       9,    10,    11,    12,    13,     0,    14,     0,    15,    16,
       0,    17,    18,    19,    20,    21,    22,     0,     0,    23,
      24,    25,    26,    27,    28,   123,   124,   125,   126,   127,
     128,   129,    29,     0,     0,    30,    31,   111,   112,   113,
     114,   115,   116,   117,   118,   119,   120,   121,   122,   123,
     124,   125,   126,   127,   128,   129,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   111,   112,   113,   114,   115,
     116,   117,   118,   119,   120,   121,   122,   123,   124,   125,
     126,   127,   128,   129,     2,     3,     0,     0,     0,     0,
       0,     4,     5,     6,     0,     0,   130,   131,     7,     8,
       0,     0,     0,     9,    10,    11,    12,     0,     0,    14,
       0,     0,     0,     0,     0,   257,     0,     0,     0,     0,
       0,     0,    55,     0,    56,    57,    58,    59,     0,     0,
       0,     0,     0,     0,   285,   111,   112,   113,   114,   115,
     116,   117,   118,   119,   120,   121,   122,   123,   124,   125,
     126,   127,   128,   129,   111,   112,   113,   114,   115,   116,
     117,   118,   119,   120,   121,   122,   123,   124,   125,   126,
     127,   128,   129,     2,     3,     0,     0,     0,     0,     0,
       4,     5,     6,     0,     0,     0,     0,     7,     8,     0,
       0,     0,     9,    10,    11,   154,     0,     0,    14,     0,
       0,     0,     0,     0,   286,     0,     0,     0,     0,     0,
       0,    55,     0,   155,   156,   157,   158,     2,     3,     0,
       0,     0,     0,   372,     4,     5,     6,     0,     0,     0,
     159,     7,     8,     0,     0,     0,     9,    10,    11,   154,
       0,     0,    14,   119,   120,   121,   122,   123,   124,   125,
     126,   127,   128,   129,     0,    55,     0,   317,   318,   319,
     320,     2,     3,     0,     0,     0,     0,     0,     4,     5,
       6,     0,     0,     0,   159,     7,     8,     0,     0,     0,
       9,    10,    11,    12,     0,     0,    14,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    55,
       0,    56,    57,    58,    59,     2,     3,     0,     0,     0,
       0,     0,     4,     5,     6,     0,     0,     0,   159,     7,
       8,     0,     0,     0,     9,    10,    11,    12,     0,     0,
      14,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    23,     0,    56,    57,    58,    59,     2,
       3,     0,     0,     0,     7,     8,     4,     5,     6,     9,
      10,    11,    12,     7,     8,    14,     0,     0,     9,    10,
      11,   235,     0,     0,    14,     0,     0,     0,    55,     0,
      56,    57,    58,    59,     0,     0,     0,    55,     0,   236,
     237,   238,   239,     2,     3,     0,     0,     0,     0,     0,
       4,     5,     6,     0,     0,     0,     0,     7,     8,     0,
       0,     0,     9,    10,    11,   154,     0,     0,    14,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    55,     0,   349,   350,   351,   352,   111,   112,   113,
     114,   115,   116,   117,   118,   119,   120,   121,   122,   123,
     124,   125,   126,   127,   128,   129,   113,   114,   115,   116,
     117,   118,   119,   120,   121,   122,   123,   124,   125,   126,
     127,   128,   129,   115,   116,   117,   118,   119,   120,   121,
     122,   123,   124,   125,   126,   127,   128,   129,   117,   118,
     119,   120,   121,   122,   123,   124,   125,   126,   127,   128,
     129
};

static const yytype_int16 yycheck[] =
{
      34,     7,    24,    40,     8,    12,    37,     5,     5,   138,
       5,     5,     5,     5,     5,    21,     5,   248,     5,    34,
       5,     5,    91,     6,   255,     5,    79,     5,    81,    98,
      83,     5,     5,    42,     5,    42,     5,    42,     5,     5,
      38,     5,    74,    74,    38,   132,   133,    38,    54,    58,
      74,   268,    35,    60,    61,    62,    63,     6,    38,     0,
      38,   278,    96,   100,    38,    38,    72,    38,     6,    38,
      37,    37,    42,    37,    77,    97,    79,    74,    81,    74,
      83,    74,    74,   152,   137,    74,    35,    74,     6,    74,
      74,    97,    35,    99,    45,   101,    34,    45,    36,    42,
     104,   105,   106,   107,   108,   109,   110,    46,   111,   112,
     113,   114,   115,   116,   117,   118,   119,   120,   121,   122,
     123,   124,   125,   126,   127,   128,   129,   130,   131,   258,
     152,   362,     5,    35,   137,   138,    16,     5,   142,   173,
      42,     5,    12,     5,     5,   151,    64,    65,    66,    67,
      68,    69,     5,    35,    74,     5,     5,    34,   144,   288,
      42,    75,     5,   197,    37,    12,     5,    34,   172,    37,
       5,    34,    42,    37,     5,    37,    37,     6,    16,   248,
     249,   250,   251,    58,    37,   254,   255,    37,    37,   258,
      60,    61,    62,    63,    37,    42,    34,    35,    37,   233,
     186,   187,    37,   272,    16,    34,    37,   276,    42,     6,
      20,    21,   281,    60,    61,    62,    63,    27,    28,    29,
      34,    58,    36,    35,    34,    35,   248,   296,    38,    39,
      40,    41,    42,   255,   303,    45,   305,    34,   307,    36,
     309,    12,    34,   312,    36,    34,    35,     6,    58,   283,
      60,    61,    62,    63,   257,   258,     6,    42,    34,   328,
      36,    34,   266,    36,    34,    35,    42,    77,    42,     6,
     256,    42,    12,   279,   280,    34,    34,    36,    36,    42,
      42,   267,   285,   286,    34,   288,    12,   309,    42,    60,
      61,    62,    63,   362,   363,   364,   365,    34,   367,    36,
       5,    12,    42,    58,    34,   291,    36,   311,    12,    31,
     296,    33,    42,   347,    42,   384,    42,   323,    35,    17,
      60,    61,    62,    63,    34,    42,   395,     6,   397,     6,
     399,    42,    17,   402,    60,    61,    62,    63,    42,    42,
     362,    58,    42,    60,    61,    62,    63,   333,   334,    60,
      61,    62,    63,    74,    42,    38,    60,    61,    62,    63,
      60,    61,    62,    63,    22,    23,    24,    25,    26,   372,
      58,     5,    60,    61,    62,    63,   382,   383,    34,     6,
      36,    35,   416,    77,    35,    58,    42,   373,    42,    35,
       6,   428,   378,     6,    38,    56,     4,   401,   432,   433,
     404,     3,   408,   389,    58,   439,    60,    61,    62,    63,
      58,    58,    38,     6,    74,   449,    28,    29,    20,    21,
      32,     5,    34,    35,    36,    27,    28,    29,    34,    34,
      36,     5,    34,    35,     5,    34,    42,    39,    40,    41,
      42,    43,    38,    45,   430,    47,    48,     6,    50,    51,
      52,    53,    54,    55,    35,    38,    58,    59,    60,    61,
      62,    63,     3,     6,     6,     5,    38,    36,    58,    71,
      58,    75,    74,    75,    76,    28,    29,    34,    36,    20,
      21,    34,    35,    36,    38,    38,    27,    28,    29,    37,
      36,    38,    74,    34,    35,   328,   188,    42,    39,    40,
      41,    42,    43,    37,    45,    37,    47,    48,    37,    50,
      51,    52,    53,    54,    55,   267,   145,    58,    59,    60,
      61,    62,    63,     3,   379,    46,   282,   172,    -1,    -1,
      71,    -1,    -1,    74,    75,    76,    -1,    -1,    -1,    -1,
      20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,
      -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,
      40,    41,    42,    43,    -1,    45,    -1,    47,    48,    -1,
      50,    51,    52,    53,    54,    55,    -1,    -1,    58,    59,
      60,    61,    62,    63,     3,    -1,    -1,    -1,    -1,    -1,
      -1,    71,    -1,    -1,    74,    75,    76,    -1,    -1,    -1,
      -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,
      29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    43,    -1,    45,    -1,    47,    48,
      -1,    50,    51,    52,    53,    54,    55,    -1,    -1,    58,
      59,    60,    61,    62,    63,    20,    21,    22,    23,    24,
      25,    26,    71,    -1,    -1,    74,    75,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    20,    21,    -1,    -1,    -1,    -1,
      -1,    27,    28,    29,    -1,    -1,    77,    78,    34,    35,
      -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,
      -1,    -1,    -1,    -1,    -1,    77,    -1,    -1,    -1,    -1,
      -1,    -1,    58,    -1,    60,    61,    62,    63,    -1,    -1,
      -1,    -1,    -1,    -1,    77,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    20,    21,    -1,    -1,    -1,    -1,    -1,
      27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    77,    -1,    -1,    -1,    -1,    -1,
      -1,    58,    -1,    60,    61,    62,    63,    20,    21,    -1,
      -1,    -1,    -1,    77,    27,    28,    29,    -1,    -1,    -1,
      77,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      -1,    -1,    45,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    -1,    58,    -1,    60,    61,    62,
      63,    20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,
      29,    -1,    -1,    -1,    77,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    58,
      -1,    60,    61,    62,    63,    20,    21,    -1,    -1,    -1,
      -1,    -1,    27,    28,    29,    -1,    -1,    -1,    77,    34,
      35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
      45,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    58,    -1,    60,    61,    62,    63,    20,
      21,    -1,    -1,    -1,    34,    35,    27,    28,    29,    39,
      40,    41,    42,    34,    35,    45,    -1,    -1,    39,    40,
      41,    42,    -1,    -1,    45,    -1,    -1,    -1,    58,    -1,
      60,    61,    62,    63,    -1,    -1,    -1,    58,    -1,    60,
      61,    62,    63,    20,    21,    -1,    -1,    -1,    -1,    -1,
      27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    58,    -1,    60,    61,    62,    63,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    20,    21,    27,    28,    29,    34,    35,    39,
      40,    41,    42,    43,    45,    47,    48,    50,    51,    52,
      53,    54,    55,    58,    59,    60,    61,    62,    63,    71,
      74,    75,    80,    81,    82,    88,    91,    92,    95,    98,
     102,   104,   105,   108,   110,   115,   116,   117,   118,   119,
     123,   124,   125,   126,    34,    58,    60,    61,    62,    63,
     115,   117,   122,    42,    96,    97,    45,    45,    74,    74,
     115,    16,    34,    35,    42,    58,    88,    34,    36,    34,
      36,    34,    36,    34,    36,    42,     0,    46,   110,    99,
      16,    35,    42,    89,    94,    74,    82,    34,    35,    34,
      75,     5,    74,   124,     6,    64,    65,    66,    67,    68,
      69,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      77,    78,    31,    33,    28,    29,    32,    34,    35,    36,
     115,    37,     5,    38,    34,     5,    74,    74,    58,   115,
      58,     6,    34,    42,    42,    60,    61,    62,    63,    77,
      85,   119,   120,   121,    42,   121,    42,   121,    42,   121,
      42,    34,    35,   113,   114,    42,    86,    58,    42,    58,
      60,    61,    62,    63,    84,    85,     6,    35,     5,    76,
     103,    88,   106,   115,    85,   109,   115,    82,   115,   117,
     117,   117,   117,   117,   117,   117,   119,   119,   119,   119,
     119,   119,   119,   119,   119,   119,   119,   119,   119,   119,
     119,   119,   119,   119,   119,   119,   119,   123,   123,   121,
      38,   120,    42,    37,   117,    42,    60,    61,    62,    63,
      87,   118,    96,    17,    37,     5,    38,   115,    35,    60,
      62,    63,    83,    85,    88,    34,     6,    77,     5,    37,
      37,    37,    37,    42,   122,   110,     6,    34,     5,    74,
      17,    38,     5,    38,   118,   118,    35,    89,   107,    74,
      74,     5,    77,    37,    76,    77,    77,    37,     5,    38,
     110,     6,    42,    42,    42,    42,     5,    37,    35,    58,
      35,    74,    83,    12,    85,    12,    85,    12,    85,     5,
      37,     6,    12,    85,    83,   118,   119,    60,    61,    62,
      63,    85,   120,     6,    38,   117,    87,    86,    35,    42,
      90,    93,    85,     6,    35,    38,    86,   115,   115,    85,
     126,   110,    56,   112,   119,   119,   120,     4,   118,    60,
      61,    62,    63,    85,   118,    58,    38,    58,    38,    85,
      85,    85,    35,    60,    62,    63,    85,    88,   101,   117,
      85,    37,    77,     6,   115,    74,    37,    84,    34,     5,
     118,   118,    74,    74,     5,    38,    34,    38,   110,     6,
       5,    38,    35,    38,    83,    12,    85,    12,    85,    12,
      85,     6,    12,    85,     6,   100,   119,   118,     5,    38,
     118,    90,    38,   115,   115,    85,   111,    36,   118,    58,
      58,    38,    85,    85,    85,   117,    85,   117,    75,   115,
      34,    37,    37,    37,   110,    36,    38,    38,    74,    82,
      37,    74,   118,   110,   110,    36,    76,    42,    37,    37,
      37,   110
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
#line 329 "lg.ypp"
    {
                        if(  ffapi::ff_justcompile) exit(0);
    // clean FH  mach 2014
		        const char *  magicffglut="#!ffglutdata4.0\n";// for complex and vector 3d plot
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
#line 403 "lg.ypp"
    {(yyval.cinst) = (yyvsp[(1) - (1)].cexp);;}
    break;

  case 5:
#line 404 "lg.ypp"
    {(yyval.cinst) = ((yyvsp[(1) - (2)].cinst)+=(yyvsp[(2) - (2)].cexp));;}
    break;

  case 6:
#line 410 "lg.ypp"
    { (yyval.clist_id) = new ListOfId();;}
    break;

  case 7:
#line 411 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str)));;}
    break;

  case 8:
#line 412 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (3)].str),(yyvsp[(3) - (3)].cexp)));;}
    break;

  case 9:
#line 413 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,2> **>()));;}
    break;

  case 10:
#line 414 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,2> **>(),true));;}
    break;

  case 11:
#line 415 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,3> **>()));;}
    break;

  case 12:
#line 416 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,3> **>(),true));;}
    break;

  case 13:
#line 417 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),Find((yyvsp[(1) - (2)].str)),atype<FE<double,4> **>()));;}
    break;

  case 14:
#line 418 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),Find((yyvsp[(1) - (3)].str)),atype<FE<double,4> **>(),true));;}
    break;

  case 15:
#line 419 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (2)].str),C_F0(),(yyvsp[(1) - (2)].type)->right()));;}
    break;

  case 16:
#line 420 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str),C_F0(),(yyvsp[(1) - (3)].type),true));;}
    break;

  case 17:
#line 421 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(2) - (3)].clist_id)));;}
    break;

  case 18:
#line 422 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (3)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str)));;}
    break;

  case 19:
#line 423 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (5)].clist_id)));;}
    break;

  case 20:
#line 424 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (5)].str),(yyvsp[(5) - (5)].cexp)));;}
    break;

  case 21:
#line 425 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,2> **>()));;}
    break;

  case 22:
#line 426 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,2> **>(),true));;}
    break;

  case 23:
#line 427 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,3> **>()));;}
    break;

  case 24:
#line 428 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,3> **>(),true));;}
    break;

  case 25:
#line 429 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),Find((yyvsp[(3) - (4)].str)),atype<FE<double,4> **>()));;}
    break;

  case 26:
#line 430 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),Find((yyvsp[(3) - (5)].str)),atype<FE<double,4> **>(),true));;}
    break;

  case 27:
#line 431 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (4)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(4) - (4)].str),C_F0(),(yyvsp[(3) - (4)].type)->right()));;}
    break;

  case 28:
#line 432 "lg.ypp"
    { (yyval.clist_id) = (yyvsp[(1) - (5)].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str),C_F0(),(yyvsp[(3) - (5)].type),true));;}
    break;

  case 29:
#line 435 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str)));;}
    break;

  case 30:
#line 436 "lg.ypp"
    { (yyval.clist_id)=(yyvsp[(1) - (3)].clist_id)  ; (yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str)));;}
    break;

  case 36:
#line 444 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (1)].str),dcltype);;}
    break;

  case 37:
#line 445 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (3)].str),dcltype,(yyvsp[(3) - (3)].cexp));;}
    break;

  case 38:
#line 446 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[(1) - (4)].str),dcltype,(yyvsp[(3) - (4)].args));(yyvsp[(3) - (4)].args).destroy();;}
    break;

  case 39:
#line 447 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 40:
#line 453 "lg.ypp"
    { (yyval.args)=(yyvsp[(1) - (1)].cexp);;}
    break;

  case 41:
#line 454 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (2)].str));;}
    break;

  case 42:
#line 455 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (2)].str));;}
    break;

  case 43:
#line 456 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (2)].str));;}
    break;

  case 44:
#line 457 "lg.ypp"
    { (yyval.args)=Find((yyvsp[(1) - (2)].str));;}
    break;

  case 45:
#line 458 "lg.ypp"
    { (yyval.args)=make_pair<const char *,const C_F0>((const char *) (yyvsp[(1) - (3)].str),(C_F0) (yyvsp[(3) - (3)].cexp));;}
    break;

  case 46:
#line 459 "lg.ypp"
    { (yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp));;}
    break;

  case 47:
#line 460 "lg.ypp"
    { (yyval.args)= ((yyvsp[(1) - (5)].args)+= make_pair<const char *,const C_F0>((const char *) (yyvsp[(3) - (5)].str),(C_F0) (yyvsp[(5) - (5)].cexp)));;}
    break;

  case 49:
#line 464 "lg.ypp"
    {(yyval.type)=TypeArray((yyvsp[(1) - (4)].type),(yyvsp[(3) - (4)].type));;}
    break;

  case 50:
#line 465 "lg.ypp"
    {(yyval.type)=TypeArray(TypeArray((yyvsp[(1) - (7)].type),(yyvsp[(3) - (7)].type)),(yyvsp[(6) - (7)].type));;}
    break;

  case 51:
#line 466 "lg.ypp"
    {(yyval.type)=TypeArray((yyvsp[(1) - (6)].type),(yyvsp[(3) - (6)].type),(yyvsp[(5) - (6)].type));;}
    break;

  case 52:
#line 467 "lg.ypp"
    {(yyval.type)=TypeArray(TypeArray((yyvsp[(1) - (9)].type),(yyvsp[(3) - (9)].type),(yyvsp[(5) - (9)].type)),(yyvsp[(8) - (9)].type));;}
    break;

  case 53:
#line 468 "lg.ypp"
    {(yyval.type)=TypeTemplate((yyvsp[(1) - (4)].type),(yyvsp[(3) - (4)].type));;}
    break;

  case 54:
#line 469 "lg.ypp"
    {(yyval.type)=TypeArray(TypeTemplate((yyvsp[(1) - (7)].type),(yyvsp[(3) - (7)].type)),(yyvsp[(6) - (7)].type));;}
    break;

  case 55:
#line 470 "lg.ypp"
    {(yyval.type)=TypeArray(TypeTemplate((yyvsp[(1) - (9)].type),(yyvsp[(3) - (9)].type)),(yyvsp[(6) - (9)].type),(yyvsp[(8) - (9)].type));;}
    break;

  case 56:
#line 474 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(1) - (1)].str),currentblock,fespacetype,fespacecomplex,fespacedim); ;}
    break;

  case 57:
#line 475 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(1) - (4)].str),currentblock,fespacetype,(yyvsp[(3) - (4)].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 58:
#line 476 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(1) - (3)].str),currentblock,fespacetype,(yyvsp[(3) - (3)].cexp),fespacecomplex,fespacedim);;}
    break;

  case 59:
#line 477 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(2) - (3)].clist_id),currentblock,fespacetype,fespacecomplex,fespacedim);;}
    break;

  case 60:
#line 478 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(2) - (6)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (6)].cexp),fespacecomplex,fespacedim);;}
    break;

  case 61:
#line 479 "lg.ypp"
    { (yyval.cexp) =  NewFEvariable((yyvsp[(2) - (5)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (5)].cexp),fespacecomplex,fespacedim);;}
    break;

  case 62:
#line 483 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(1) - (4)].str),currentblock,fespacetype,(yyvsp[(3) - (4)].cexp),fespacecomplex,fespacedim); ;}
    break;

  case 63:
#line 484 "lg.ypp"
    { (yyval.cexp) =  NewFEarray((yyvsp[(2) - (6)].clist_id),currentblock,fespacetype,(yyvsp[(5) - (6)].cexp),fespacecomplex,fespacedim);;}
    break;

  case 64:
#line 489 "lg.ypp"
    { fespacedim=2;;}
    break;

  case 65:
#line 490 "lg.ypp"
    { fespacedim=1;;}
    break;

  case 66:
#line 491 "lg.ypp"
    { fespacedim=3;;}
    break;

  case 67:
#line 492 "lg.ypp"
    { fespacedim=4;;}
    break;

  case 68:
#line 495 "lg.ypp"
    {fespacecomplex=false;  fespacetype = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 69:
#line 496 "lg.ypp"
    {
             if ((yyvsp[(3) - (4)].type) != typevarreal && (yyvsp[(3) - (4)].type) != typevarcomplex) lgerror (" type of finite element <real> or <complex>");
             fespacecomplex=((yyvsp[(3) - (4)].type)==typevarcomplex);
             fespacetype = Find((yyvsp[(1) - (4)].str));;}
    break;

  case 70:
#line 501 "lg.ypp"
    {  (yyval.cexp) = (yyvsp[(1) - (1)].cexp);;}
    break;

  case 71:
#line 502 "lg.ypp"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 72:
#line 504 "lg.ypp"
    {  (yyval.cexp) = (yyvsp[(1) - (1)].cexp);;}
    break;

  case 73:
#line 505 "lg.ypp"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 74:
#line 507 "lg.ypp"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[(2) - (2)].cexp);;}
    break;

  case 75:
#line 508 "lg.ypp"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[(5) - (5)].cexp);;}
    break;

  case 76:
#line 512 "lg.ypp"
    {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,size_t>((yyvsp[(1) - (4)].str),typeFESpace((yyvsp[(3) - (4)].args)),(yyvsp[(3) - (4)].args),dimFESpaceImage((yyvsp[(3) - (4)].args)));
     (yyvsp[(3) - (4)].args).destroy(); ;}
    break;

  case 78:
#line 516 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 79:
#line 519 "lg.ypp"
    {dcltype=(yyvsp[(1) - (1)].type);;}
    break;

  case 80:
#line 519 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(3) - (4)].cexp);;}
    break;

  case 81:
#line 520 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(2) - (3)].cexp);;}
    break;

  case 82:
#line 521 "lg.ypp"
    { (yyval.cexp)=(yyvsp[(1) - (2)].cexp);;}
    break;

  case 83:
#line 522 "lg.ypp"
    {(yyval.cexp)=currentblock->NewID((yyvsp[(1) - (5)].type),(yyvsp[(2) - (5)].str),(yyvsp[(4) - (5)].cexp));;}
    break;

  case 84:
#line 524 "lg.ypp"
    {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = (yyvsp[(2) - (6)].type)->right();
                      routineinblock[kkembtype] = currentblock;
                      (yyvsp[(5) - (6)].routine)=new Routine((yyvsp[(1) - (6)].type),(yyvsp[(2) - (6)].type)->right(),(yyvsp[(3) - (6)].str),(yyvsp[(5) - (6)].clist_id),currentblock);
		      // routineinblock[kkembtype]->Add($3,"(",$<routine>5); //pas recursif pour l'instanat test  FH 27 dec 2008
                     // cout << " \n after new routine \n " << endl;
                      ;}
    break;

  case 85:
#line 533 "lg.ypp"
    { currentblock=(yyvsp[(5) - (10)].routine)->Set((yyvsp[(9) - (10)].cinst));
                       currentblock->Add((yyvsp[(3) - (10)].str),"(",(yyvsp[(5) - (10)].routine)); //pas recursif pour l'instant test  FH 27 dec 2008
                       kkembtype--;
                       (yyval.cexp)=0;

                        ;}
    break;

  case 86:
#line 540 "lg.ypp"
    {Block::open(currentblock); (yyvsp[(1) - (5)].type)->SetArgs((yyvsp[(4) - (5)].clist_id));;}
    break;

  case 87:
#line 542 "lg.ypp"
    {  //$<cinst>$=currentblock->close(currentblock);
                         (yyval.cinst).setclose(Block::snewclose(currentblock));// Sep 2016 FH.
                         (yyval.cexp)=currentblock->NewID((yyvsp[(1) - (9)].type),(yyvsp[(2) - (9)].str),(yyvsp[(8) - (9)].cexp),*(yyvsp[(4) - (9)].clist_id));
                         delete (yyvsp[(4) - (9)].clist_id); //  FH 23032005
                         ;}
    break;

  case 88:
#line 549 "lg.ypp"
    {  Block::open(currentblock);;}
    break;

  case 89:
#line 550 "lg.ypp"
    { (yyval.endb)=Block::snewclose(currentblock);
//  $$=currentblock->close(currentblock);
;}
    break;

  case 90:
#line 554 "lg.ypp"
    {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 91:
#line 556 "lg.ypp"
    {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 92:
#line 560 "lg.ypp"
    {dcltype=(yyvsp[(1) - (1)].type); Block::open(currentblock);  ;}
    break;

  case 93:
#line 561 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(3) - (3)].cexp);;}
    break;

  case 94:
#line 563 "lg.ypp"
    { Block::open(currentblock);;}
    break;

  case 95:
#line 565 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (1)].str)));Block::open(currentblock); ;}
    break;

  case 96:
#line 566 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (3)].str)));(yyval.clist_id)->push_back(UnId((yyvsp[(3) - (3)].str)));Block::open(currentblock); ;}
    break;

  case 97:
#line 567 "lg.ypp"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[(1) - (5)].str)));(yyval.clist_id)->push_back(UnId((yyvsp[(3) - (5)].str)));(yyval.clist_id)->push_back(UnId((yyvsp[(5) - (5)].str)));Block::open(currentblock); ;}
    break;

  case 98:
#line 569 "lg.ypp"
    {(yyval.cexp)=0;;}
    break;

  case 99:
#line 570 "lg.ypp"
    {zzzfff->input((yyvsp[(2) - (2)].str));(yyval.cexp)= 0; ;}
    break;

  case 100:
#line 571 "lg.ypp"
    {load((yyvsp[(2) - (2)].str));(yyval.cexp)= 0; ;}
    break;

  case 101:
#line 572 "lg.ypp"
    {(yyval.cexp)=Try((yyvsp[(3) - (5)].cinst),currentblock->close(currentblock,(yyvsp[(5) - (5)].cexp)));;}
    break;

  case 102:
#line 573 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(1) - (2)].cexp);;}
    break;

  case 103:
#line 574 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(1) - (1)].cexp);;}
    break;

  case 104:
#line 576 "lg.ypp"
    {(yyvsp[(5) - (6)].cexp)=ForAll(currentblock,(yyvsp[(3) - (6)].clist_id),(yyvsp[(5) - (6)].cexp));;}
    break;

  case 105:
#line 577 "lg.ypp"
    {
                    inloopcount--;
                    (yyval.cexp)=Block::close(currentblock,C_F0(ForAll((yyvsp[(5) - (8)].cexp),(yyvsp[(8) - (8)].cexp))));
                 ;}
    break;

  case 106:
#line 581 "lg.ypp"
    {inloopcount--; (yyval.cexp)=For((yyvsp[(3) - (9)].cexp),(yyvsp[(5) - (9)].cexp),(yyvsp[(7) - (9)].cexp),(yyvsp[(9) - (9)].cexp));;}
    break;

  case 107:
#line 583 "lg.ypp"
    {inloopcount--;
                 (yyval.cexp)=Block::close(currentblock,C_F0(For((yyvsp[(3) - (9)].cexp),(yyvsp[(5) - (9)].cexp),(yyvsp[(7) - (9)].cexp),(yyvsp[(9) - (9)].cexp))));
                ;}
    break;

  case 108:
#line 586 "lg.ypp"
    {inloopcount--;(yyval.cexp)=While((yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 109:
#line 587 "lg.ypp"
    {(yyval.cexp)=FIf((yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 110:
#line 588 "lg.ypp"
    {(yyval.cexp)=FIf((yyvsp[(3) - (7)].cexp),(yyvsp[(5) - (7)].cexp),(yyvsp[(7) - (7)].cexp));;}
    break;

  case 111:
#line 589 "lg.ypp"
    { /* [[begin:]] [[end:]] */
             (yyvsp[(2) - (3)].cinst).setclose((yyvsp[(3) - (3)].endb));
             (yyval.cexp)=(yyvsp[(2) - (3)].cinst);
                    //  $$=C_F0(new E_block($2,$3),atype<void>());
         ;}
    break;

  case 112:
#line 594 "lg.ypp"
    { /* <<BORDER_ID>> */
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[(2) - (3)].str),C_F0(TheOperators,"[border]",(yyvsp[(3) - (3)].args)));;}
    break;

  case 113:
#line 596 "lg.ypp"
    {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[(2) - (6)].str),C_F0(TheOperators,"[border]",(yyvsp[(4) - (6)].args)));;}
    break;

  case 114:
#line 599 "lg.ypp"
    {
                    if(inloopcount)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_break),atype<void>());
                    else lgerror("break not in loop");;}
    break;

  case 115:
#line 603 "lg.ypp"
    {
                    if(inloopcount)
                        (yyval.cexp)= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");;}
    break;

  case 116:
#line 607 "lg.ypp"
    {
                    if (kkembtype>=0)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_return,(rettype[kkembtype]->CastTo((yyvsp[(2) - (3)].cexp))).OnReturn()) ,atype<void>());
                     else lgerror(" return not in routine ");;}
    break;

  case 117:
#line 614 "lg.ypp"
    {(yyval.cexp) =  (yyvsp[(7) - (7)].cexp); ;}
    break;

  case 118:
#line 617 "lg.ypp"
    {
   Block::open(currentblock);
   (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[(2) - (7)].str),atype<double*>());
   (yyval.args)+= (yyvsp[(4) - (7)].cexp);
   (yyval.args)+= (yyvsp[(6) - (7)].cexp);
   (yyval.args)+= currentblock->NewVar<LocalVariable>("IndexBorder",atype<long*>());;}
    break;

  case 119:
#line 625 "lg.ypp"
    {
    Block::open(currentblock);
    (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[(2) - (9)].str),atype<double*>());
    (yyval.args)+= (yyvsp[(4) - (9)].cexp);
    (yyval.args)+= (yyvsp[(6) - (9)].cexp);
    (yyval.args)+= currentblock->NewVar<LocalVariable>((yyvsp[(8) - (9)].str),atype<long*>());;}
    break;

  case 120:
#line 633 "lg.ypp"
    {
    //currentblock->close(currentblock;);
   (yyval.args) = ((yyvsp[(1) - (2)].args) += currentblock->close(currentblock,(yyvsp[(2) - (2)].cexp)));
   ;}
    break;

  case 122:
#line 641 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 129:
#line 654 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 130:
#line 655 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"+=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 131:
#line 656 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"-=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 132:
#line 657 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"*=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 133:
#line 658 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"/=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 134:
#line 659 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,".*=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 135:
#line 660 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"./=",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 137:
#line 666 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"?:",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 138:
#line 667 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 139:
#line 668 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 141:
#line 673 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 142:
#line 674 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 143:
#line 675 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 144:
#line 676 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 145:
#line 677 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 146:
#line 678 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 147:
#line 679 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 148:
#line 680 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 149:
#line 681 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 150:
#line 682 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 151:
#line 683 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 152:
#line 684 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 153:
#line 685 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 154:
#line 686 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 155:
#line 687 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 156:
#line 688 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 157:
#line 689 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 158:
#line 690 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 159:
#line 691 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 160:
#line 695 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(1) - (1)].cexp);;}
    break;

  case 161:
#line 696 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,":");;}
    break;

  case 162:
#line 697 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 163:
#line 698 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[(1) - (5)].cexp),(yyvsp[(3) - (5)].cexp),(yyvsp[(5) - (5)].cexp));;}
    break;

  case 164:
#line 703 "lg.ypp"
    {(yyval.args) = 0;;}
    break;

  case 165:
#line 704 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 166:
#line 705 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 167:
#line 706 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 168:
#line 707 "lg.ypp"
    {(yyval.args) = Find((yyvsp[(1) - (1)].str));;}
    break;

  case 169:
#line 708 "lg.ypp"
    {(yyval.args) = make_pair<const char *,const C_F0>((const char *) (yyvsp[(1) - (3)].str),(C_F0) (yyvsp[(3) - (3)].cexp));;}
    break;

  case 170:
#line 709 "lg.ypp"
    {(yyval.args) = (yyvsp[(1) - (1)].cexp);;}
    break;

  case 171:
#line 710 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 172:
#line 711 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 173:
#line 712 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 174:
#line 713 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += Find((yyvsp[(3) - (3)].str)));;}
    break;

  case 175:
#line 714 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp));;}
    break;

  case 176:
#line 717 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (5)].args)+= make_pair<const char *,const C_F0>((const char *)(yyvsp[(3) - (5)].str),(C_F0) (yyvsp[(5) - (5)].cexp)));;}
    break;

  case 177:
#line 720 "lg.ypp"
    {(yyval.args)=(yyvsp[(1) - (1)].cexp);;}
    break;

  case 178:
#line 721 "lg.ypp"
    {(yyval.args) = ((yyvsp[(1) - (3)].args) += (yyvsp[(3) - (3)].cexp));;}
    break;

  case 180:
#line 726 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(1) - (2)].oper),(yyvsp[(2) - (2)].cexp));;}
    break;

  case 182:
#line 731 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 183:
#line 732 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (3)].oper),(yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].cexp));;}
    break;

  case 185:
#line 737 "lg.ypp"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));;}
    break;

  case 186:
#line 745 "lg.ypp"
    {(yyval.cexp)=Find((yyvsp[(1) - (1)].str));;}
    break;

  case 187:
#line 749 "lg.ypp"
    {(yyval.cexp)= CConstant((yyvsp[(1) - (1)].lnum));;}
    break;

  case 188:
#line 750 "lg.ypp"
    {(yyval.cexp)= CConstant((yyvsp[(1) - (1)].dnum));;}
    break;

  case 189:
#line 751 "lg.ypp"
    {(yyval.cexp)= CConstant(complex<double>(0,(yyvsp[(1) - (1)].dnum)));;}
    break;

  case 190:
#line 752 "lg.ypp"
    {(yyval.cexp)= CConstant<const char *>((yyvsp[(1) - (1)].str));;}
    break;

  case 191:
#line 757 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (4)].cexp),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args));;}
    break;

  case 192:
#line 759 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (4)].cexp),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].cexp));;}
    break;

  case 193:
#line 760 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (6)].cexp),(yyvsp[(2) - (6)].oper),(yyvsp[(3) - (6)].cexp),(yyvsp[(5) - (6)].cexp));;}
    break;

  case 194:
#line 761 "lg.ypp"
    {(yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),"[]");;}
    break;

  case 195:
#line 762 "lg.ypp"
    { (yyval.cexp)=C_F0((yyvsp[(1) - (3)].cexp),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 196:
#line 763 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 197:
#line 764 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 198:
#line 765 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 199:
#line 766 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 200:
#line 767 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 201:
#line 768 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 202:
#line 769 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (3)].str)),(yyvsp[(3) - (3)].str)) ;;}
    break;

  case 203:
#line 770 "lg.ypp"
    { (yyval.cexp)=C_F0(Find((yyvsp[(1) - (4)].str)),(yyvsp[(2) - (4)].oper),(yyvsp[(3) - (4)].args)) ;;}
    break;

  case 204:
#line 771 "lg.ypp"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));;}
    break;

  case 205:
#line 772 "lg.ypp"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[(2) - (2)].oper),(yyvsp[(1) - (2)].cexp));;}
    break;

  case 206:
#line 773 "lg.ypp"
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

  case 207:
#line 782 "lg.ypp"
    {(yyval.cexp)=(yyvsp[(2) - (3)].cexp);;}
    break;

  case 208:
#line 783 "lg.ypp"
    { (yyval.cexp)=C_F0(TheOperators,"[]",(yyvsp[(2) - (3)].args));;}
    break;


/* Line 1267 of yacc.c.  */
#line 3219 "lg.tab.cpp"
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


#line 788 "lg.ypp"



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

#include "ffServer.hpp"
extern ffServer *graphicServer;

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
      if(verbosity>1) cout << "   file :" << cc << " " << " verbosity= " << verbosity << endl;
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
  if (graphicServer) {
    // graphicServer->stop();
    delete graphicServer;
  }
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
