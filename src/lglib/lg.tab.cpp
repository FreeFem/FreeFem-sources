/* A Bison parser, made by GNU Bison 3.5.1.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2020 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

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

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.5.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         lgparse
#define yylex           lglex
#define yyerror         lgerror
#define yydebug         lgdebug
#define yynerrs         lgnerrs
#define yylval          lglval
#define yychar          lgchar

/* First part of user prologue.  */
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


#line 241 "lg.tab.cpp"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Use api.header.include to #include this header
   instead of duplicating it here.  */
#ifndef YY_LG_LG_TAB_HPP_INCLUDED
# define YY_LG_LG_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int lgdebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    IF = 258,
    ELSE = 259,
    SET = 260,
    LTLT = 261,
    GTGT = 262,
    OR = 263,
    AND = 264,
    EQ = 265,
    NE = 266,
    LE = 267,
    GE = 268,
    DOTSTAR = 269,
    DOTSLASH = 270,
    UNARY = 271,
    PLUSPLUS = 272,
    MOINSMOINS = 273,
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
    FESPACEDS = 299,
    FESPACEL = 300,
    VGFESPACE = 301,
    GFESPACE = 302,
    PLUSEQ = 303,
    MOINSEQ = 304,
    MULEQ = 305,
    DIVEQ = 306,
    DOTMULEQ = 307,
    DOTDIVEQ = 308,
    ARROW = 309,
    BORDER = 310,
    SOLVE = 311
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 170 "lg.ypp"

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

#line 384 "lg.tab.cpp"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE lglval;

int lgparse (void);

#endif /* !YY_LG_LG_TAB_HPP_INCLUDED  */



#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))

/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                            \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
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
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
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
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  112
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1519

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  82
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  50
/* YYNRULES -- Number of rules.  */
#define YYNRULES  265
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  526

#define YYUNDEFTOK  2
#define YYMAXUTOK   311


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    30,     2,     2,     2,    24,    13,    32,
      34,    37,    22,    20,     5,    21,    36,    23,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    80,    77,
      16,     6,    17,    81,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    35,     2,    38,    31,    33,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    78,    11,    79,     2,     2,     2,     2,
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
       7,     8,     9,    10,    12,    14,    15,    18,    19,    25,
      26,    27,    28,    29,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   340,   340,   410,   414,   415,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     457,   458,   463,   463,   463,   463,   463,   463,   463,   466,
     467,   468,   469,   475,   476,   477,   478,   479,   480,   481,
     482,   483,   484,   485,   486,   487,   488,   489,   494,   495,
     496,   497,   498,   499,   500,   501,   505,   506,   507,   508,
     509,   510,   514,   515,   520,   521,   522,   523,   524,   525,
     526,   527,   530,   531,   536,   537,   539,   540,   542,   543,
     546,   549,   553,   554,   557,   557,   558,   559,   560,   562,
     561,   578,   577,   587,   588,   592,   594,   598,   598,   601,
     603,   604,   605,   607,   608,   609,   610,   611,   612,   614,
     613,   619,   620,   624,   625,   626,   627,   632,   634,   637,
     641,   645,   652,   655,   663,   671,   678,   679,   683,   684,
     685,   686,   687,   691,   692,   693,   694,   695,   696,   697,
     698,   703,   704,   705,   706,   710,   711,   712,   713,   714,
     715,   716,   717,   718,   719,   720,   721,   722,   723,   724,
     725,   726,   727,   728,   729,   733,   734,   735,   736,   741,
     745,   746,   747,   748,   749,   750,   751,   752,   753,   754,
     755,   756,   757,   758,   759,   760,   761,   762,   765,   768,
     769,   772,   773,   774,   775,   776,   777,   778,   779,   780,
     781,   782,   783,   784,   785,   786,   787,   788,   789,   793,
     794,   798,   799,   800,   804,   805,   813,   817,   818,   819,
     820,   825,   827,   828,   829,   830,   831,   832,   833,   834,
     835,   836,   837,   838,   839,   840,   841,   842,   843,   844,
     845,   854,   862,   863,   864,   865
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "IF", "ELSE", "','", "'='", "SET",
  "LTLT", "GTGT", "OR", "'|'", "AND", "'&'", "EQ", "NE", "'<'", "'>'",
  "LE", "GE", "'+'", "'-'", "'*'", "'/'", "'%'", "DOTSTAR", "DOTSLASH",
  "UNARY", "PLUSPLUS", "MOINSMOINS", "'!'", "'^'", "'\\''", "'_'", "'('",
  "'['", "'.'", "')'", "']'", "LNUM", "DNUM", "CNUM", "ID", "FESPACEID",
  "IDPARAM", "STRING", "ENDOFFILE", "INCLUDE", "LOAD", "BIDON", "FOR",
  "WHILE", "BREAK", "CONTINUE", "RETURN", "TRY", "CATCH", "THROW", "TYPE",
  "FUNCTION", "FESPACE", "FESPACE1", "FESPACE3", "FESPACES", "FESPACEDS",
  "FESPACEL", "VGFESPACE", "GFESPACE", "PLUSEQ", "MOINSEQ", "MULEQ",
  "DIVEQ", "DOTMULEQ", "DOTDIVEQ", "ARROW", "BORDER", "SOLVE", "';'",
  "'{'", "'}'", "':'", "'?'", "$accept", "start", "input", "instructions",
  "list_of_id_args", "list_of_id1", "id", "list_of_dcls",
  "parameters_list", "type_of_dcl", "ID_space", "ID_array_space",
  "fespace123", "fespace", "spaceIDa", "spaceIDb", "spaceIDs",
  "fespace_def", "fespace_def_list", "declaration", "$@1", "$@2", "$@3",
  "begin", "end", "for_loop", "while_loop", "declaration_for", "$@4",
  "try", "IDfor", "instruction", "$@5", "catchs", "bornes", "border_expr",
  "Expr", "unop", "no_comma_expr", "no_set_expr", "no_ternary_expr",
  "sub_script_expr", "parameterstype", "parameters", "array", "FEarray",
  "unary_expr", "pow_expr", "primaryp", "primary", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_int16 yytoknum[] =
{
       0,   256,   257,   258,   259,    44,    61,   260,   261,   262,
     263,   124,   264,    38,   265,   266,    60,    62,   267,   268,
      43,    45,    42,    47,    37,   269,   270,   271,   272,   273,
      33,    94,    39,    95,    40,    91,    46,    41,    93,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,    59,   123,   125,
      58,    63
};
# endif

#define YYPACT_NINF (-296)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-206)

#define yytable_value_is_error(Yyn) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     756,   -21,   700,  -296,  -296,  -296,  -296,  -296,  1194,  1194,
    -296,  -296,  -296,  -296,   492,  -296,   -16,    23,  -296,  -296,
      67,    80,  1194,  -296,   129,   854,   971,  1025,  1078,  1131,
    1184,  1234,  -296,  -296,   492,  -296,  -296,   149,   126,   756,
    -296,   171,  1343,   113,  -296,   756,    38,   169,   121,  -296,
      21,  1303,  -296,   111,    24,    73,  -296,  -296,   -15,   869,
    1194,  -296,  -296,  -296,  -296,  -296,  -296,   319,   190,    35,
     243,   262,   284,   306,   312,    17,  -296,    19,  -296,  -296,
    -296,  -296,  -296,  -296,  -296,   141,  -296,    22,  -296,  -296,
    -296,  -296,    25,   181,  1141,   223,   192,   191,   492,   982,
     492,   982,   492,   982,   492,   982,   492,   982,   492,   982,
     492,   117,  -296,  -296,  -296,   492,   227,  1373,   115,  -296,
     302,  -296,   530,  1244,   492,  1194,   756,  1194,  -296,  -296,
    1194,  1194,  1194,  1194,  1194,  1194,  1194,  1194,  1194,  1194,
    1194,  1194,  1194,  1194,  1194,  1194,  1194,  1194,  1194,  1194,
    1194,  1194,  1194,  1194,  1194,  1194,  1194,  1194,  1365,  1454,
    1194,  1194,  -296,  -296,  -296,   982,  1035,   492,    56,  -296,
    -296,  1194,  -296,  1294,  1294,   492,  -296,  -296,   297,  -296,
     376,    60,   288,    26,  1194,  1311,   287,   320,   206,   228,
     268,   373,   398,   502,  -296,   337,  -296,    72,  -296,   157,
    -296,   165,  -296,   168,  -296,   173,  -296,   193,  -296,   492,
    1194,   756,  -296,   220,    46,   333,   316,    53,  -296,  1194,
    1194,  1349,  -296,  -296,  -296,   281,    54,   355,   325,   201,
     583,  -296,  -296,  -296,  -296,  -296,  -296,  -296,  -296,  1448,
    1448,  1463,  1463,  1476,  1476,  1487,  1487,   804,   804,   804,
     804,  1357,  1357,  -296,  -296,  -296,  -296,  -296,   827,   846,
    -296,  -296,  -296,  -296,  -296,  -296,  -296,  -296,  -296,  -296,
    -296,  -296,  -296,  -296,   208,  -296,    62,  -296,   756,  -296,
      77,   186,   195,   887,   894,   932,   377,  -296,  -296,   239,
    -296,   326,  1194,   982,  -296,  -296,   323,   378,    55,   384,
    1311,   154,   269,   303,   310,   364,   247,   402,   430,  1311,
    1194,  1088,  -296,  -296,  -296,  -296,  -296,  -296,   406,    69,
    -296,  1194,  1294,   492,  -296,  -296,  1359,   492,   150,  -296,
     383,   492,  -296,   492,  1194,  1194,   492,  1303,   756,   366,
    1194,  1194,  -296,  1141,  -296,   419,  -296,  -296,  -296,  -296,
    -296,  -296,  1194,  1294,  -296,   372,   865,   428,   399,   380,
    -296,  1311,    70,   492,  -296,   492,  -296,   492,  -296,   492,
    -296,   492,  -296,  1335,  -296,  1194,   492,  -296,   256,  -296,
     647,   659,   905,  1047,  1100,  1153,  -296,   433,  -296,  1194,
     363,  -296,   273,  -296,   492,   411,  -296,   445,  -296,  1194,
    1194,  -296,   447,    57,    58,   452,   745,  -296,   412,  -296,
    1431,  1431,   421,   756,  -296,  -296,    74,  1194,   426,   424,
      84,  -296,  -296,  -296,  -296,  -296,  -296,   429,  1311,   441,
     614,   627,   638,   457,   919,   465,  -296,  -296,  -296,  1194,
     468,  -296,  -296,   110,  1194,  1359,  -296,   439,  1194,  1194,
     492,  -296,   443,  -296,  -296,   423,  -296,  1431,   427,  -296,
     470,  1311,   131,   492,  -296,   492,  -296,   492,  -296,   492,
    -296,  1194,   492,  -296,  1194,   410,  -296,  1194,   455,   459,
    -296,  -296,   298,   301,  -296,   756,   461,   460,   471,  -296,
     138,   482,  -296,  -296,  -296,  -296,  -296,  -296,   434,   756,
     -25,  1194,  -296,   756,   756,  -296,   474,  -296,  -296,   518,
    -296,  -296,   691,  -296,   492,   500,  -296,  -296,   505,  -296,
    -296,   506,  -296,   756,  -296,  -296
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int16 yydefact[] =
{
       0,     0,     0,   149,   148,   151,   152,   150,     0,     0,
     237,   238,   239,   236,     0,   240,     0,     0,   115,   116,
       0,     0,     0,   119,    68,     0,   211,   212,   213,   214,
     216,   215,    89,    90,     0,   123,   113,     0,     0,     3,
     104,    92,     0,     0,   128,     0,     0,     0,     0,     4,
       0,     0,   146,   153,   161,   265,   165,   229,   231,   234,
       0,   211,   212,   213,   214,   216,   215,     0,     0,   211,
     212,   213,   214,   216,   215,     0,   209,     0,    42,    43,
      47,    44,    45,    48,    46,     0,   102,     0,   124,   125,
     139,   140,     0,     0,     0,     0,    68,     0,     0,   190,
       0,   190,     0,   190,     0,   190,     0,   190,     0,   190,
       0,     0,     1,     2,     5,     0,     0,     0,    76,    96,
      98,   107,     0,     0,     0,     0,     0,     0,   127,   230,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   258,   259,   235,   190,     0,     0,     0,   264,
     262,     0,   263,     0,     0,     0,   106,   141,     0,   186,
     185,     0,     0,     0,     0,     6,     0,   236,   211,   212,
     213,   214,   216,   215,   197,     0,   199,     0,   246,     0,
     248,     0,   250,     0,   252,     0,   256,     0,   254,     0,
       0,     0,   137,    49,     0,     0,     0,     0,    40,     0,
       0,     0,   114,   136,   117,     0,     0,   120,     0,     0,
       0,   147,   154,   155,   156,   157,   158,   159,   160,   173,
     174,   178,   177,   176,   175,   183,   184,   179,   181,   180,
     182,   171,   172,   166,   169,   170,   167,   168,   163,     0,
     217,   218,   219,   220,   222,   221,   223,   224,   225,   226,
     228,   227,   232,   233,     0,   244,     0,   245,     0,   210,
     211,   212,   213,   214,   216,   215,     0,   101,    53,     0,
     103,    73,     0,   190,   260,   261,     0,    69,     0,     0,
       6,    43,    44,    45,    48,    46,     0,     7,     0,     6,
       0,     0,   247,   249,   251,   253,   257,   255,     0,     0,
     145,     0,     0,     0,   105,    93,     0,     0,    79,    78,
       0,     0,    97,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   241,     0,   242,   134,    54,    55,    56,    57,
      64,    58,     0,     0,   100,     0,   187,   189,     0,     0,
     108,     6,     0,     0,     9,     0,    11,     0,    13,     0,
      15,     0,    17,     0,   111,     0,     0,    19,     0,   198,
     211,   212,   213,   214,   216,   215,   206,     0,   207,     0,
       0,    50,     0,    52,     0,     0,    94,    99,    41,     0,
       0,    77,   118,     0,     0,   121,     0,   133,     0,   126,
     164,   162,     0,     0,    66,    67,     0,     0,    71,     0,
       0,    21,    10,    12,    14,    16,    18,     0,     6,    43,
      44,    45,    46,    24,     0,     0,     8,    20,   109,     0,
       0,   138,    51,     0,     0,     0,    81,     0,     0,     0,
       0,   129,     0,   243,   135,     0,    74,   188,     0,    70,
      23,     6,     0,     0,    30,     0,    32,     0,    34,     0,
      36,     0,     0,    38,     0,     0,   208,     0,     0,     0,
      95,    80,     0,     0,   122,     0,     0,     0,     0,    22,
       0,    25,    31,    33,    35,    37,    29,    39,     0,     0,
     147,     0,    82,     0,     0,   130,     0,    75,    72,    28,
      26,   112,     0,   143,     0,     0,   132,   131,     0,    27,
     110,     0,    83,     0,   144,   142
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -296,  -296,  -296,   -39,  -295,   151,   -14,  -253,  -167,   -17,
     327,    99,  -296,  -296,  -296,  -296,  -296,   374,  -296,  -296,
    -296,  -296,  -296,  -296,  -296,  -296,  -296,  -296,  -296,  -296,
    -296,   -37,  -296,  -296,  -296,  -296,    -7,  -296,    -6,  -164,
     375,   -90,  -296,   -84,   351,   545,   139,   511,  -296,   226
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    37,    38,    39,   306,   217,   195,   214,   287,    40,
     119,   396,    41,    42,   397,   120,    43,    86,    87,    44,
     115,   475,   435,    45,   223,    46,    47,   225,   333,    48,
     228,    49,   485,   409,   211,   212,    50,    51,    52,    53,
      54,   196,   182,   197,    77,    55,    56,    57,    58,    59
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      85,    75,   114,    76,   181,   362,   122,   289,    98,   288,
     288,    97,   513,    60,   378,    92,   160,   199,   161,   201,
     111,   203,   127,   205,   171,   207,   127,   175,   118,    88,
     127,   296,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   323,   514,   168,   170,   329,   330,   172,   327,   127,
     127,   127,   127,   127,   297,   293,   420,   343,    89,    99,
     393,   100,   123,   124,   171,   373,   276,   311,   158,   455,
     402,   274,   -59,   -43,   186,   114,   198,   230,   200,   373,
     202,   328,   204,   278,   206,   159,   208,   294,   128,   176,
     344,   213,   177,   218,   156,   157,   224,   390,   421,   312,
     227,    99,   456,   100,   -59,   327,   226,   130,   229,    78,
     231,   219,   460,   324,   232,   233,   234,   235,   236,   237,
     238,   335,   360,   462,   448,   449,   373,    79,    80,    81,
      82,    83,    84,   373,    90,    93,   379,   173,   478,   112,
     220,   209,   210,   277,   -59,   392,   399,    91,   288,   286,
     286,    85,   311,    94,    95,   279,   490,   363,   308,   491,
     311,   307,   113,   311,   320,   174,   509,   298,   311,   131,
     132,   133,   134,   135,   136,   400,   415,   116,   414,   288,
     121,   -60,   -47,   114,   313,   318,    78,   184,   311,   126,
     -61,   -44,   314,   125,    76,   315,   127,   118,    93,   357,
     316,  -191,   -43,   311,    79,    80,    81,    82,    83,    84,
     101,   388,   102,   -60,    94,   185,   321,    95,    78,   103,
     317,   104,   -61,  -192,   -47,   446,   447,    78,   338,   178,
      99,   345,   100,  -191,   353,   342,    79,    80,    81,    82,
      83,    84,   373,   412,   322,    79,    80,    81,    82,    83,
      84,   373,   101,   -60,   102,  -192,   346,   347,   348,   349,
     350,   351,   -61,  -193,   -44,   476,   354,   101,   353,   102,
     479,   183,   365,   308,   374,   215,   307,   364,   366,   368,
     370,   372,   308,   438,   377,   307,   103,   387,   104,   272,
     273,   407,   103,   127,   104,  -193,   127,   221,   286,   213,
     442,    78,   395,   398,   291,   391,   367,   218,   105,   213,
     106,   309,   405,   369,   158,   295,   -42,   403,   404,    79,
      80,    81,    82,    83,    84,   503,   169,   515,   504,   286,
     107,   159,   108,   310,   308,    78,   109,   307,   110,   422,
     325,   423,    78,   424,   326,   425,   434,   426,   334,   433,
     336,   355,   437,    79,    80,    81,    82,    83,    84,   436,
      79,    80,    81,    82,    83,    84,   454,   371,  -194,   -45,
     218,   358,   440,   352,   137,   138,   139,   140,   141,   142,
     143,   144,   145,   146,   147,   148,   149,   150,   151,   152,
     153,   154,   155,  -196,   -48,   337,    78,   105,   375,   106,
    -194,   308,   389,   359,   307,   464,   466,   468,   470,   361,
     473,   401,   408,   413,    79,    80,    81,    82,    83,    84,
     416,   395,   107,   311,   108,  -196,   484,   418,   419,   439,
     441,   482,   483,   376,   308,   444,   452,   307,   505,   492,
     445,   493,   323,   494,   463,   495,   292,   450,   497,   453,
     512,   458,   459,   471,   461,   496,   516,   517,   498,   180,
     500,   474,    78,   477,   180,   114,   180,   481,   180,   486,
     180,   487,   180,    78,   180,   488,   525,   489,   499,   501,
      79,    80,    81,    82,    83,    84,   502,   506,   507,   510,
     521,    79,    80,    81,    82,    83,    84,  -195,   -46,   508,
     518,   511,   239,   240,   241,   242,   243,   244,   245,   246,
     247,   248,   249,   250,   251,   252,   253,   254,   255,   256,
     257,   258,   259,     1,    78,   519,   109,   522,   110,  -195,
     180,   180,   523,   524,   480,   443,     2,    67,   332,   290,
       3,     4,    79,    80,    81,    82,    83,    84,     5,     6,
       7,   319,   129,   406,     8,     9,     0,     0,     0,    10,
      11,    12,    13,    14,     0,    15,     0,    16,    17,     0,
      18,    19,    20,    21,    22,    23,     1,     0,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,     0,     2,
       0,     0,     0,     3,     4,    34,     0,    35,    36,   222,
       0,     5,     6,     7,     0,     0,     0,     8,     9,     0,
       0,     0,    10,    11,    12,    13,    14,   465,    15,     0,
      16,    17,     0,    18,    19,    20,    21,    22,    23,     0,
     467,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,   469,  -200,   -43,     0,     0,    78,     0,    34,     0,
      35,    36,   339,     0,  -201,   -47,     0,   356,   180,    78,
       0,     0,     0,     0,    79,    80,    81,    82,    83,    84,
      78,    99,     0,   100,  -200,     0,   180,    79,    80,    81,
      82,    83,    84,   101,     1,   102,  -201,     0,    79,    80,
      81,    82,    83,    84,     0,     0,     0,     2,     0,     0,
       0,     3,     4,     0,     0,   410,   411,     0,   180,     5,
       6,     7,     0,     0,     0,     8,     9,     0,     0,     0,
      10,    11,    12,    13,    14,     0,    15,     0,    16,    17,
       0,    18,    19,    20,    21,    22,    23,     0,     0,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,     1,
      61,    62,    63,    64,    65,    66,    34,     0,    35,    36,
     520,     0,     2,   162,   163,     0,     3,     4,     0,   165,
     166,   167,     0,   451,     5,     6,     7,     0,     0,     0,
       8,     9,   457,     0,     0,    10,    11,    12,    13,    14,
       0,    15,     0,    16,    17,     0,    18,    19,    20,    21,
      22,    23,     0,     0,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,   149,   150,   151,   152,   153,   154,
     155,    34,     0,    35,    36,   137,   138,   139,   140,   141,
     142,   143,   144,   145,   146,   147,   148,   149,   150,   151,
     152,   153,   154,   155,   137,   138,   139,   140,   141,   142,
     143,   144,   145,   146,   147,   148,   149,   150,   151,   152,
     153,   154,   155,   137,   138,   139,   140,   141,   142,   143,
     144,   145,   146,   147,   148,   149,   150,   151,   152,   153,
     154,   155,   -62,   -45,     0,     0,    78,   162,   163,   -65,
     -48,   164,     0,   165,   166,   167,     0,   340,     0,     0,
    -202,   -44,    96,     0,    79,    80,    81,    82,    83,    84,
       0,   105,     0,   106,   -62,     0,   341,     0,   107,    78,
     108,   -65,   472,     0,     0,     0,    78,   -63,   -46,   103,
       0,   104,  -202,     0,     0,   417,     0,    79,    80,    81,
      82,    83,    84,     0,    79,    80,    81,    82,    83,    84,
       0,    78,     0,     0,   -62,     0,   109,     0,   110,   -63,
       0,   -65,     0,     0,    78,     0,     0,     0,     0,    79,
      80,    81,    82,    83,    84,     0,     0,   -84,     0,     0,
       0,     0,    79,    80,    81,    82,    83,    84,     2,     0,
       0,     0,     3,     4,     0,    99,   -84,   100,     0,   -63,
       5,     6,     7,   -84,     0,     0,     8,     9,     0,     0,
       0,    10,    11,    12,   187,     0,     0,    15,     0,     0,
       0,   -84,   -84,   -84,   -84,   -84,   -84,     0,     0,     0,
      68,   -85,   188,   189,   190,   191,   192,   193,   194,     0,
       0,     2,  -203,   -45,     0,     3,     4,     0,     0,   101,
     -85,   102,   179,     5,     6,     7,     0,   -85,     0,     8,
       9,     0,     0,   275,    10,    11,    12,    13,     0,     0,
      15,   105,     0,   106,  -203,   -85,   -85,   -85,   -85,   -85,
     -85,     0,     0,    68,   -86,    69,    70,    71,    72,    73,
      74,     0,     0,     0,     2,  -205,   -48,     0,     3,     4,
       0,     0,   103,   -86,   104,   179,     5,     6,     7,     0,
     -86,     0,     8,     9,     0,     0,     0,    10,    11,    12,
     187,     0,     0,    15,   107,     0,   108,  -205,   -86,   -86,
     -86,   -86,   -86,   -86,     0,     0,    68,   -87,   380,   381,
     382,   383,   384,   385,   386,     0,     0,     2,  -204,   -46,
       0,     3,     4,     0,     0,   105,   -87,   106,   179,     5,
       6,     7,     0,   -87,     0,     8,     9,     0,     0,     0,
      10,    11,    12,    13,     0,     0,    15,   109,     0,   110,
    -204,   -87,   -87,   -87,   -87,   -87,   -87,     0,     0,    68,
     -91,    69,    70,    71,    72,    73,    74,     0,     0,     0,
       2,     0,     0,     0,     3,     4,     0,     0,   107,   -91,
     108,   179,     5,     6,     7,     0,   -91,     0,     8,     9,
       0,     0,     0,    10,    11,    12,    13,     0,     0,    15,
       0,     0,     0,     0,   -91,   -91,   -91,   -91,   -91,   -91,
     -88,     0,    68,     0,    69,    70,    71,    72,    73,    74,
       2,     0,     0,     0,     3,     4,     0,     0,   109,   -88,
     110,     0,     5,     6,     7,     0,   -88,     0,     8,     9,
       0,     0,     0,    10,    11,    12,    13,     0,     0,    15,
       0,     0,     0,     0,   -88,   -88,   -88,   -88,   -88,   -88,
       0,     0,    24,     0,    69,    70,    71,    72,    73,    74,
       2,     0,     0,     0,     3,     4,     0,     0,     0,     2,
       0,     0,     5,     6,     7,     0,     0,   299,     8,     9,
       0,     0,     0,    10,    11,    12,   187,     8,     9,    15,
       0,     0,    10,    11,    12,    13,   300,     0,    15,     0,
       0,   427,    68,    78,   280,   281,   282,   283,   284,   285,
       0,    68,     0,    69,    70,    71,    72,    73,    74,    96,
     428,   301,    80,   302,   303,   304,   305,    78,   117,   151,
     152,   153,   154,   155,   331,    78,     0,     0,     0,     0,
       0,    78,     0,    96,   394,   429,    80,   430,   431,    83,
     432,    78,     0,    79,    80,    81,    82,    83,    84,    79,
      80,    81,    82,    83,    84,    78,     0,     0,     0,    79,
      80,    81,    82,    83,    84,   260,   261,   262,   263,   264,
     265,   216,     0,    79,    80,    81,    82,    83,    84,   137,
     138,   139,   140,   141,   142,   143,   144,   145,   146,   147,
     148,   149,   150,   151,   152,   153,   154,   155,   139,   140,
     141,   142,   143,   144,   145,   146,   147,   148,   149,   150,
     151,   152,   153,   154,   155,   141,   142,   143,   144,   145,
     146,   147,   148,   149,   150,   151,   152,   153,   154,   155,
     143,   144,   145,   146,   147,   148,   149,   150,   151,   152,
     153,   154,   155,   145,   146,   147,   148,   149,   150,   151,
     152,   153,   154,   155,   266,   267,   268,   269,   270,   271
};

static const yytype_int16 yycheck[] =
{
      14,     8,    39,     9,    94,   300,    45,   174,    25,   173,
     174,    25,    37,    34,   309,    22,    31,   101,    33,   103,
      34,   105,     5,   107,     5,   109,     5,     5,    42,    45,
       5,     5,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,     5,    77,    60,    37,   219,   220,    38,     5,     5,
       5,     5,     5,     5,    38,     5,   361,     5,    45,    34,
     323,    36,    34,    35,     5,     5,   166,     5,     5,     5,
     333,   165,     5,     6,    98,   122,   100,   126,   102,     5,
     104,    38,   106,    37,   108,    22,   110,    37,    77,    77,
      38,   115,    77,   117,    80,    81,   123,    38,    38,    37,
     124,    34,    38,    36,    37,     5,   123,     6,   125,    42,
     127,     6,    38,    77,   130,   131,   132,   133,   134,   135,
     136,    77,    77,   428,    77,    77,     5,    60,    61,    62,
      63,    64,    65,     5,    77,    16,   310,     6,    38,     0,
      35,    34,    35,   167,    77,   322,     6,    77,   322,   173,
     174,   175,     5,    34,    35,   171,   461,    13,   185,    38,
       5,   185,    46,     5,   211,    34,    38,   184,     5,    68,
      69,    70,    71,    72,    73,    35,   353,    16,   352,   353,
      77,     5,     6,   230,    37,   209,    42,     6,     5,    78,
       5,     6,    37,    34,   210,    37,     5,   221,    16,   293,
      37,     5,     6,     5,    60,    61,    62,    63,    64,    65,
      34,   311,    36,    37,    34,    34,     6,    35,    42,    34,
      37,    36,    37,     5,     6,   399,   400,    42,    37,    58,
      34,   278,    36,    37,     5,    37,    60,    61,    62,    63,
      64,    65,     5,   343,    34,    60,    61,    62,    63,    64,
      65,     5,    34,    77,    36,    37,   280,   281,   282,   283,
     284,   285,    77,     5,     6,   439,    37,    34,     5,    36,
     444,    58,    13,   300,    37,    58,   300,   301,   302,   303,
     304,   305,   309,    37,   308,   309,    34,   311,    36,   160,
     161,   338,    34,     5,    36,    37,     5,     5,   322,   323,
      37,    42,   326,   327,    17,   321,    13,   331,    34,   333,
      36,    34,   336,    13,     5,    37,     6,   334,   335,    60,
      61,    62,    63,    64,    65,    37,    17,   501,    37,   353,
      34,    22,    36,     6,   361,    42,    34,   361,    36,   363,
      17,   365,    42,   367,    38,   369,   373,   371,    77,   373,
       5,    35,   376,    60,    61,    62,    63,    64,    65,   375,
      60,    61,    62,    63,    64,    65,   413,    13,     5,     6,
     394,    58,   389,     6,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,     5,     6,    80,    42,    34,     6,    36,
      37,   428,     6,    35,   428,   429,   430,   431,   432,    35,
     434,    38,    56,     4,    60,    61,    62,    63,    64,    65,
      58,   445,    34,     5,    36,    37,   450,    38,    58,     6,
      77,   448,   449,    13,   461,    34,    34,   461,   485,   463,
       5,   465,     5,   467,    13,   469,    80,     5,   472,    38,
     499,    35,    38,     6,    35,   471,   503,   504,   474,    94,
     477,     6,    42,     5,    99,   512,   101,    38,   103,    36,
     105,    58,   107,    42,   109,    58,   523,    17,    78,    34,
      60,    61,    62,    63,    64,    65,    37,    36,    38,    17,
     514,    60,    61,    62,    63,    64,    65,     5,     6,    38,
      36,    77,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,     3,    42,    17,    34,    37,    36,    37,
     165,   166,    37,    37,   445,   394,    16,     2,   221,   175,
      20,    21,    60,    61,    62,    63,    64,    65,    28,    29,
      30,   210,    51,   337,    34,    35,    -1,    -1,    -1,    39,
      40,    41,    42,    43,    -1,    45,    -1,    47,    48,    -1,
      50,    51,    52,    53,    54,    55,     3,    -1,    58,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    -1,    16,
      -1,    -1,    -1,    20,    21,    75,    -1,    77,    78,    79,
      -1,    28,    29,    30,    -1,    -1,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    43,    13,    45,    -1,
      47,    48,    -1,    50,    51,    52,    53,    54,    55,    -1,
      13,    58,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    13,     5,     6,    -1,    -1,    42,    -1,    75,    -1,
      77,    78,    79,    -1,     5,     6,    -1,   292,   293,    42,
      -1,    -1,    -1,    -1,    60,    61,    62,    63,    64,    65,
      42,    34,    -1,    36,    37,    -1,   311,    60,    61,    62,
      63,    64,    65,    34,     3,    36,    37,    -1,    60,    61,
      62,    63,    64,    65,    -1,    -1,    -1,    16,    -1,    -1,
      -1,    20,    21,    -1,    -1,   340,   341,    -1,   343,    28,
      29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    43,    -1,    45,    -1,    47,    48,
      -1,    50,    51,    52,    53,    54,    55,    -1,    -1,    58,
      59,    60,    61,    62,    63,    64,    65,    66,    67,     3,
      60,    61,    62,    63,    64,    65,    75,    -1,    77,    78,
      79,    -1,    16,    28,    29,    -1,    20,    21,    -1,    34,
      35,    36,    -1,    38,    28,    29,    30,    -1,    -1,    -1,
      34,    35,   417,    -1,    -1,    39,    40,    41,    42,    43,
      -1,    45,    -1,    47,    48,    -1,    50,    51,    52,    53,
      54,    55,    -1,    -1,    58,    59,    60,    61,    62,    63,
      64,    65,    66,    67,    20,    21,    22,    23,    24,    25,
      26,    75,    -1,    77,    78,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,     5,     6,    -1,    -1,    42,    28,    29,     5,
       6,    32,    -1,    34,    35,    36,    -1,    80,    -1,    -1,
       5,     6,    58,    -1,    60,    61,    62,    63,    64,    65,
      -1,    34,    -1,    36,    37,    -1,    80,    -1,    34,    42,
      36,    37,    13,    -1,    -1,    -1,    42,     5,     6,    34,
      -1,    36,    37,    -1,    -1,    80,    -1,    60,    61,    62,
      63,    64,    65,    -1,    60,    61,    62,    63,    64,    65,
      -1,    42,    -1,    -1,    77,    -1,    34,    -1,    36,    37,
      -1,    77,    -1,    -1,    42,    -1,    -1,    -1,    -1,    60,
      61,    62,    63,    64,    65,    -1,    -1,    16,    -1,    -1,
      -1,    -1,    60,    61,    62,    63,    64,    65,    16,    -1,
      -1,    -1,    20,    21,    -1,    34,    35,    36,    -1,    77,
      28,    29,    30,    42,    -1,    -1,    34,    35,    -1,    -1,
      -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,
      -1,    60,    61,    62,    63,    64,    65,    -1,    -1,    -1,
      58,    16,    60,    61,    62,    63,    64,    65,    66,    -1,
      -1,    16,     5,     6,    -1,    20,    21,    -1,    -1,    34,
      35,    36,    80,    28,    29,    30,    -1,    42,    -1,    34,
      35,    -1,    -1,    38,    39,    40,    41,    42,    -1,    -1,
      45,    34,    -1,    36,    37,    60,    61,    62,    63,    64,
      65,    -1,    -1,    58,    16,    60,    61,    62,    63,    64,
      65,    -1,    -1,    -1,    16,     5,     6,    -1,    20,    21,
      -1,    -1,    34,    35,    36,    80,    28,    29,    30,    -1,
      42,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
      42,    -1,    -1,    45,    34,    -1,    36,    37,    60,    61,
      62,    63,    64,    65,    -1,    -1,    58,    16,    60,    61,
      62,    63,    64,    65,    66,    -1,    -1,    16,     5,     6,
      -1,    20,    21,    -1,    -1,    34,    35,    36,    80,    28,
      29,    30,    -1,    42,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    -1,    -1,    45,    34,    -1,    36,
      37,    60,    61,    62,    63,    64,    65,    -1,    -1,    58,
      16,    60,    61,    62,    63,    64,    65,    -1,    -1,    -1,
      16,    -1,    -1,    -1,    20,    21,    -1,    -1,    34,    35,
      36,    80,    28,    29,    30,    -1,    42,    -1,    34,    35,
      -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,
      -1,    -1,    -1,    -1,    60,    61,    62,    63,    64,    65,
      16,    -1,    58,    -1,    60,    61,    62,    63,    64,    65,
      16,    -1,    -1,    -1,    20,    21,    -1,    -1,    34,    35,
      36,    -1,    28,    29,    30,    -1,    42,    -1,    34,    35,
      -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,
      -1,    -1,    -1,    -1,    60,    61,    62,    63,    64,    65,
      -1,    -1,    58,    -1,    60,    61,    62,    63,    64,    65,
      16,    -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    16,
      -1,    -1,    28,    29,    30,    -1,    -1,    16,    34,    35,
      -1,    -1,    -1,    39,    40,    41,    42,    34,    35,    45,
      -1,    -1,    39,    40,    41,    42,    35,    -1,    45,    -1,
      -1,    16,    58,    42,    60,    61,    62,    63,    64,    65,
      -1,    58,    -1,    60,    61,    62,    63,    64,    65,    58,
      35,    60,    61,    62,    63,    64,    65,    42,    35,    22,
      23,    24,    25,    26,    35,    42,    -1,    -1,    -1,    -1,
      -1,    42,    -1,    58,    35,    60,    61,    62,    63,    64,
      65,    42,    -1,    60,    61,    62,    63,    64,    65,    60,
      61,    62,    63,    64,    65,    42,    -1,    -1,    -1,    60,
      61,    62,    63,    64,    65,    60,    61,    62,    63,    64,
      65,    58,    -1,    60,    61,    62,    63,    64,    65,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    60,    61,    62,    63,    64,    65
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    16,    20,    21,    28,    29,    30,    34,    35,
      39,    40,    41,    42,    43,    45,    47,    48,    50,    51,
      52,    53,    54,    55,    58,    59,    60,    61,    62,    63,
      64,    65,    66,    67,    75,    77,    78,    83,    84,    85,
      91,    94,    95,    98,   101,   105,   107,   108,   111,   113,
     118,   119,   120,   121,   122,   127,   128,   129,   130,   131,
      34,    60,    61,    62,    63,    64,    65,   127,    58,    60,
      61,    62,    63,    64,    65,   118,   120,   126,    42,    60,
      61,    62,    63,    64,    65,    88,    99,   100,    45,    45,
      77,    77,   118,    16,    34,    35,    58,    88,    91,    34,
      36,    34,    36,    34,    36,    34,    36,    34,    36,    34,
      36,    88,     0,    46,   113,   102,    16,    35,    88,    92,
      97,    77,    85,    34,    35,    34,    78,     5,    77,   129,
       6,    68,    69,    70,    71,    72,    73,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    80,    81,     5,    22,
      31,    33,    28,    29,    32,    34,    35,    36,   118,    17,
      37,     5,    38,     6,    34,     5,    77,    77,    58,    80,
     122,   123,   124,    58,     6,    34,    88,    42,    60,    61,
      62,    63,    64,    65,    66,    88,   123,   125,    88,   125,
      88,   125,    88,   125,    88,   125,    88,   125,    88,    34,
      35,   116,   117,    88,    89,    58,    58,    87,    88,     6,
      35,     5,    79,   106,    91,   109,   118,    88,   112,   118,
      85,   118,   120,   120,   120,   120,   120,   120,   120,   122,
     122,   122,   122,   122,   122,   122,   122,   122,   122,   122,
     122,   122,   122,   122,   122,   122,   122,   122,   122,   122,
      60,    61,    62,    63,    64,    65,    60,    61,    62,    63,
      64,    65,   128,   128,   125,    38,   123,    88,    37,   120,
      60,    61,    62,    63,    64,    65,    88,    90,   121,    90,
      99,    17,    80,     5,    37,    37,     5,    38,   118,    16,
      35,    60,    62,    63,    64,    65,    86,    88,    91,    34,
       6,     5,    37,    37,    37,    37,    37,    37,    88,   126,
     113,     6,    34,     5,    77,    17,    38,     5,    38,   121,
     121,    35,    92,   110,    77,    77,     5,    80,    37,    79,
      80,    80,    37,     5,    38,   113,    88,    88,    88,    88,
      88,    88,     6,     5,    37,    35,   122,   125,    58,    35,
      77,    35,    86,    13,    88,    13,    88,    13,    88,    13,
      88,    13,    88,     5,    37,     6,    13,    88,    86,   121,
      60,    61,    62,    63,    64,    65,    66,    88,   123,     6,
      38,   120,    90,    89,    35,    88,    93,    96,    88,     6,
      35,    38,    89,   118,   118,    88,   131,   113,    56,   115,
     122,   122,   123,     4,   121,    90,    58,    80,    38,    58,
      86,    38,    88,    88,    88,    88,    88,    16,    35,    60,
      62,    63,    65,    88,    91,   104,   120,    88,    37,     6,
     118,    77,    37,    87,    34,     5,   121,   121,    77,    77,
       5,    38,    34,    38,   113,     5,    38,   122,    35,    38,
      38,    35,    86,    13,    88,    13,    88,    13,    88,    13,
      88,     6,    13,    88,     6,   103,   121,     5,    38,   121,
      93,    38,   118,   118,    88,   114,    36,    58,    58,    17,
      86,    38,    88,    88,    88,    88,   120,    88,   120,    78,
     118,    34,    37,    37,    37,   113,    36,    38,    38,    38,
      17,    77,    85,    37,    77,   121,   113,   113,    36,    17,
      79,    88,    37,    37,    37,   113
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    82,    83,    84,    85,    85,    86,    86,    86,    86,
      86,    86,    86,    86,    86,    86,    86,    86,    86,    86,
      86,    86,    86,    86,    86,    86,    86,    86,    86,    86,
      86,    86,    86,    86,    86,    86,    86,    86,    86,    86,
      87,    87,    88,    88,    88,    88,    88,    88,    88,    89,
      89,    89,    89,    90,    90,    90,    90,    90,    90,    90,
      90,    90,    90,    90,    90,    90,    90,    90,    91,    91,
      91,    91,    91,    91,    91,    91,    92,    92,    92,    92,
      92,    92,    93,    93,    94,    94,    94,    94,    94,    94,
      94,    94,    95,    95,    96,    96,    97,    97,    98,    98,
      99,    99,   100,   100,   102,   101,   101,   101,   101,   103,
     101,   104,   101,   105,   106,   107,   108,   110,   109,   111,
     112,   112,   112,   113,   113,   113,   113,   113,   113,   114,
     113,   113,   113,   113,   113,   113,   113,   113,   113,   113,
     113,   113,   115,   116,   116,   117,   118,   118,   119,   119,
     119,   119,   119,   120,   120,   120,   120,   120,   120,   120,
     120,   121,   121,   121,   121,   122,   122,   122,   122,   122,
     122,   122,   122,   122,   122,   122,   122,   122,   122,   122,
     122,   122,   122,   122,   122,   123,   123,   123,   123,   124,
     125,   125,   125,   125,   125,   125,   125,   125,   125,   125,
     125,   125,   125,   125,   125,   125,   125,   125,   125,   126,
     126,   127,   127,   127,   127,   127,   127,   127,   127,   127,
     127,   127,   127,   127,   127,   127,   127,   127,   127,   128,
     128,   129,   129,   129,   130,   130,   131,   131,   131,   131,
     131,   131,   131,   131,   131,   131,   131,   131,   131,   131,
     131,   131,   131,   131,   131,   131,   131,   131,   131,   131,
     131,   131,   131,   131,   131,   131
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     2,     1,     1,     2,     0,     1,     3,     2,
       3,     2,     3,     2,     3,     2,     3,     2,     3,     2,
       3,     3,     5,     4,     3,     5,     6,     7,     6,     5,
       4,     5,     4,     5,     4,     5,     4,     5,     4,     5,
       1,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       3,     4,     3,     1,     2,     2,     2,     2,     2,     1,
       1,     1,     1,     1,     2,     1,     3,     3,     1,     4,
       7,     6,     9,     4,     7,     9,     1,     4,     3,     3,
       6,     5,     4,     6,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     4,     1,     3,     1,     3,     2,     5,
       4,     3,     1,     3,     0,     4,     3,     2,     5,     0,
      10,     0,     9,     1,     1,     1,     1,     0,     3,     1,
       1,     3,     5,     1,     2,     2,     5,     2,     1,     0,
       8,     9,     9,     5,     5,     7,     3,     3,     6,     2,
       2,     3,     7,     7,     9,     2,     1,     3,     1,     1,
       1,     1,     1,     1,     3,     3,     3,     3,     3,     3,
       3,     1,     5,     3,     5,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     1,     1,     3,     5,     3,
       0,     1,     1,     1,     1,     1,     1,     1,     3,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     5,     1,
       3,     1,     1,     1,     1,     1,     1,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       2,     1,     3,     3,     1,     2,     1,     1,     1,     1,
       1,     4,     4,     6,     3,     3,     3,     4,     3,     4,
       3,     4,     3,     4,     3,     4,     3,     4,     2,     2,
       4,     4,     3,     3,     3,     1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YYUSE (yyoutput);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyo, yytoknum[yytype], *yyvaluep);
# endif
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyo, yytype, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[+yyssp[yyi + 1 - yynrhs]],
                       &yyvsp[(yyi + 1) - (yynrhs)]
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
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
#ifndef YYINITDEPTH
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
#   define yystrlen(S) (YY_CAST (YYPTRDIFF_T, strlen (S)))
#  else
/* Return the length of YYSTR.  */
static YYPTRDIFF_T
yystrlen (const char *yystr)
{
  YYPTRDIFF_T yylen;
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
static char *
yystpcpy (char *yydest, const char *yysrc)
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
static YYPTRDIFF_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYPTRDIFF_T yyn = 0;
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
            else
              goto append;

          append:
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

  if (yyres)
    return yystpcpy (yyres, yystr) - yyres;
  else
    return yystrlen (yystr);
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYPTRDIFF_T *yymsg_alloc, char **yymsg,
                yy_state_t *yyssp, int yytoken)
{
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat: reported tokens (one for the "unexpected",
     one per "expected"). */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Actual size of YYARG. */
  int yycount = 0;
  /* Cumulated lengths of YYARG.  */
  YYPTRDIFF_T yysize = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[+*yyssp];
      YYPTRDIFF_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
      yysize = yysize0;
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYPTRDIFF_T yysize1
                    = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
                    yysize = yysize1;
                  else
                    return 2;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
    default: /* Avoid compiler warnings. */
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    /* Don't count the "%s"s in the final size, but reserve room for
       the terminator.  */
    YYPTRDIFF_T yysize1 = yysize + (yystrlen (yyformat) - 2 * yycount) + 1;
    if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
      yysize = yysize1;
    else
      return 2;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          ++yyp;
          ++yyformat;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss;
    yy_state_t *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYPTRDIFF_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYPTRDIFF_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    goto yyexhaustedlab;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
# undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
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
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
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
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

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
#line 340 "lg.ypp"
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
                       (yyvsp[-1].cinst).setclose(Block::snewclose(currentblock));// Sep 2016 FH
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
                          (yyvsp[-1].cinst).eval(stack);
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


			    if (verbosity && (NbPtr || (stu1>100000) )) { cout << " ######## unfreed pointers   " << NbPtr
			                      << " Nb pointer,   " <<  lg1-lg0 << "Bytes " << " ,  mpirank " << mpirank << ", memory leak ="<< stu1 <<  endl;}
  return 0;}
#line 2189 "lg.tab.cpp"
    break;

  case 4:
#line 414 "lg.ypp"
                                                                          {(yyval.cinst) = (yyvsp[0].cexp);}
#line 2195 "lg.tab.cpp"
    break;

  case 5:
#line 415 "lg.ypp"
                                                                          {(yyval.cinst) = ((yyvsp[-1].cinst)+=(yyvsp[0].cexp));}
#line 2201 "lg.tab.cpp"
    break;

  case 6:
#line 421 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId();}
#line 2207 "lg.tab.cpp"
    break;

  case 7:
#line 422 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
#line 2213 "lg.tab.cpp"
    break;

  case 8:
#line 423 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str),(yyvsp[0].cexp)));}
#line 2219 "lg.tab.cpp"
    break;

  case 9:
#line 424 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,2> **>()));}
#line 2225 "lg.tab.cpp"
    break;

  case 10:
#line 425 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,2> **>(),true));}
#line 2231 "lg.tab.cpp"
    break;

  case 11:
#line 426 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,3> **>()));}
#line 2237 "lg.tab.cpp"
    break;

  case 12:
#line 427 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,3> **>(),true));}
#line 2243 "lg.tab.cpp"
    break;

  case 13:
#line 428 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,4> **>()));}
#line 2249 "lg.tab.cpp"
    break;

  case 14:
#line 429 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,4> **>(),true));}
#line 2255 "lg.tab.cpp"
    break;

  case 15:
#line 430 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,40> **>())); }
#line 2261 "lg.tab.cpp"
    break;

  case 16:
#line 431 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,40> **>(),true)); }
#line 2267 "lg.tab.cpp"
    break;

  case 17:
#line 432 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,5> **>()));}
#line 2273 "lg.tab.cpp"
    break;

  case 18:
#line 433 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,5> **>(),true));}
#line 2279 "lg.tab.cpp"
    break;

  case 19:
#line 434 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-1].type)->right()));}
#line 2285 "lg.tab.cpp"
    break;

  case 20:
#line 435 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-2].type),true));}
#line 2291 "lg.tab.cpp"
    break;

  case 21:
#line 436 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id)));}
#line 2297 "lg.tab.cpp"
    break;

  case 22:
#line 437 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-2].clist_id),true,true));}
#line 2303 "lg.tab.cpp"
    break;

  case 23:
#line 438 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id),true,false));}
#line 2309 "lg.tab.cpp"
    break;

  case 24:
#line 439 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-2].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
#line 2315 "lg.tab.cpp"
    break;

  case 25:
#line 440 "lg.ypp"
                                                { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id)));}
#line 2321 "lg.tab.cpp"
    break;

  case 26:
#line 441 "lg.ypp"
                                                    { (yyval.clist_id) = (yyvsp[-5].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-2].clist_id),false,true));}
#line 2327 "lg.tab.cpp"
    break;

  case 27:
#line 442 "lg.ypp"
                                                        { (yyval.clist_id) = (yyvsp[-6].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-2].clist_id),true,true));}
#line 2333 "lg.tab.cpp"
    break;

  case 28:
#line 443 "lg.ypp"
                                                     { (yyval.clist_id) = (yyvsp[-5].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id),true,false));}
#line 2339 "lg.tab.cpp"
    break;

  case 29:
#line 444 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str),(yyvsp[0].cexp)));}
#line 2345 "lg.tab.cpp"
    break;

  case 30:
#line 445 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,2> **>()));}
#line 2351 "lg.tab.cpp"
    break;

  case 31:
#line 446 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,2> **>(),true));}
#line 2357 "lg.tab.cpp"
    break;

  case 32:
#line 447 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,3> **>()));}
#line 2363 "lg.tab.cpp"
    break;

  case 33:
#line 448 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,3> **>(),true));}
#line 2369 "lg.tab.cpp"
    break;

  case 34:
#line 449 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,4> **>()));}
#line 2375 "lg.tab.cpp"
    break;

  case 35:
#line 450 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,4> **>(),true));}
#line 2381 "lg.tab.cpp"
    break;

  case 36:
#line 451 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,5> **>()));}
#line 2387 "lg.tab.cpp"
    break;

  case 37:
#line 452 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,5> **>(),true));}
#line 2393 "lg.tab.cpp"
    break;

  case 38:
#line 453 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-1].type)->right()));}
#line 2399 "lg.tab.cpp"
    break;

  case 39:
#line 454 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-2].type),true));}
#line 2405 "lg.tab.cpp"
    break;

  case 40:
#line 457 "lg.ypp"
                                      { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
#line 2411 "lg.tab.cpp"
    break;

  case 41:
#line 458 "lg.ypp"
                                      { (yyval.clist_id)=(yyvsp[-2].clist_id)  ; (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
#line 2417 "lg.tab.cpp"
    break;

  case 49:
#line 466 "lg.ypp"
                                    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[0].str),dcltype);}
#line 2423 "lg.tab.cpp"
    break;

  case 50:
#line 467 "lg.ypp"
                                    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[-2].str),dcltype,(yyvsp[0].cexp));}
#line 2429 "lg.tab.cpp"
    break;

  case 51:
#line 468 "lg.ypp"
                                    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[-3].str),dcltype,(yyvsp[-1].args));(yyvsp[-1].args).destroy();}
#line 2435 "lg.tab.cpp"
    break;

  case 52:
#line 469 "lg.ypp"
                                    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2441 "lg.tab.cpp"
    break;

  case 53:
#line 475 "lg.ypp"
                                                 { (yyval.args)=(yyvsp[0].cexp);}
#line 2447 "lg.tab.cpp"
    break;

  case 54:
#line 476 "lg.ypp"
                                                 { (yyval.args)=Find((yyvsp[-1].str));}
#line 2453 "lg.tab.cpp"
    break;

  case 55:
#line 477 "lg.ypp"
                                                 { (yyval.args)=Find((yyvsp[-1].str));}
#line 2459 "lg.tab.cpp"
    break;

  case 56:
#line 478 "lg.ypp"
                                                 { (yyval.args)=Find((yyvsp[-1].str));}
#line 2465 "lg.tab.cpp"
    break;

  case 57:
#line 479 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[-1].str));}
#line 2471 "lg.tab.cpp"
    break;

  case 58:
#line 480 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[-1].str));}
#line 2477 "lg.tab.cpp"
    break;

  case 59:
#line 481 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2483 "lg.tab.cpp"
    break;

  case 60:
#line 482 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2489 "lg.tab.cpp"
    break;

  case 61:
#line 483 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2495 "lg.tab.cpp"
    break;

  case 62:
#line 484 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2501 "lg.tab.cpp"
    break;

  case 63:
#line 485 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2507 "lg.tab.cpp"
    break;

  case 64:
#line 486 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[-1].str));}
#line 2513 "lg.tab.cpp"
    break;

  case 65:
#line 487 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2519 "lg.tab.cpp"
    break;

  case 66:
#line 488 "lg.ypp"
                                                 { (yyval.args)=make_pair<const char *,const C_F0>((const char *) (yyvsp[-2].str),(C_F0) (yyvsp[0].cexp));}
#line 2525 "lg.tab.cpp"
    break;

  case 67:
#line 489 "lg.ypp"
                                                 { (yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].args));}
#line 2531 "lg.tab.cpp"
    break;

  case 69:
#line 495 "lg.ypp"
                                              {(yyval.type)=TypeArray((yyvsp[-3].type),(yyvsp[-1].type));}
#line 2537 "lg.tab.cpp"
    break;

  case 70:
#line 496 "lg.ypp"
                                              {(yyval.type)=TypeArray(TypeArray((yyvsp[-6].type),(yyvsp[-4].type)),(yyvsp[-1].type));}
#line 2543 "lg.tab.cpp"
    break;

  case 71:
#line 497 "lg.ypp"
                                              {(yyval.type)=TypeArray((yyvsp[-5].type),(yyvsp[-3].type),(yyvsp[-1].type));}
#line 2549 "lg.tab.cpp"
    break;

  case 72:
#line 498 "lg.ypp"
                                              {(yyval.type)=TypeArray(TypeArray((yyvsp[-8].type),(yyvsp[-6].type),(yyvsp[-4].type)),(yyvsp[-1].type));}
#line 2555 "lg.tab.cpp"
    break;

  case 73:
#line 499 "lg.ypp"
                                              {(yyval.type)=TypeTemplate((yyvsp[-3].type),(yyvsp[-1].type));}
#line 2561 "lg.tab.cpp"
    break;

  case 74:
#line 500 "lg.ypp"
                                              {(yyval.type)=TypeArray(TypeTemplate((yyvsp[-6].type),(yyvsp[-4].type)),(yyvsp[-1].type));}
#line 2567 "lg.tab.cpp"
    break;

  case 75:
#line 501 "lg.ypp"
                                              {(yyval.type)=TypeArray(TypeTemplate((yyvsp[-8].type),(yyvsp[-6].type)),(yyvsp[-3].type),(yyvsp[-1].type));}
#line 2573 "lg.tab.cpp"
    break;

  case 76:
#line 505 "lg.ypp"
                                              { (yyval.cexp) =  NewFEvariable((yyvsp[0].str),currentblock,fespacetype,fespacecomplex,fespacedim); }
#line 2579 "lg.tab.cpp"
    break;

  case 77:
#line 506 "lg.ypp"
                                              { (yyval.cexp) =  NewFEarray((yyvsp[-3].str),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim); }
#line 2585 "lg.tab.cpp"
    break;

  case 78:
#line 507 "lg.ypp"
                                              { (yyval.cexp) =  NewFEvariable((yyvsp[-2].str),currentblock,fespacetype,(yyvsp[0].cexp),fespacecomplex,fespacedim);}
#line 2591 "lg.tab.cpp"
    break;

  case 79:
#line 508 "lg.ypp"
                                              { (yyval.cexp) =  NewFEvariable((yyvsp[-1].clist_id),currentblock,fespacetype,fespacecomplex,fespacedim);}
#line 2597 "lg.tab.cpp"
    break;

  case 80:
#line 509 "lg.ypp"
                                              { (yyval.cexp) =  NewFEarray((yyvsp[-4].clist_id),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim);}
#line 2603 "lg.tab.cpp"
    break;

  case 81:
#line 510 "lg.ypp"
                                              { (yyval.cexp) =  NewFEvariable((yyvsp[-3].clist_id),currentblock,fespacetype,(yyvsp[0].cexp),fespacecomplex,fespacedim);}
#line 2609 "lg.tab.cpp"
    break;

  case 82:
#line 514 "lg.ypp"
                                     { (yyval.cexp) =  NewFEarray((yyvsp[-3].str),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim); }
#line 2615 "lg.tab.cpp"
    break;

  case 83:
#line 515 "lg.ypp"
                                              { (yyval.cexp) =  NewFEarray((yyvsp[-4].clist_id),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim);}
#line 2621 "lg.tab.cpp"
    break;

  case 84:
#line 520 "lg.ypp"
              { fespacedim=2;}
#line 2627 "lg.tab.cpp"
    break;

  case 85:
#line 521 "lg.ypp"
               { fespacedim=1;}
#line 2633 "lg.tab.cpp"
    break;

  case 86:
#line 522 "lg.ypp"
               { fespacedim=3;}
#line 2639 "lg.tab.cpp"
    break;

  case 87:
#line 523 "lg.ypp"
               { fespacedim=4;}
#line 2645 "lg.tab.cpp"
    break;

  case 88:
#line 524 "lg.ypp"
               { fespacedim=5;}
#line 2651 "lg.tab.cpp"
    break;

  case 89:
#line 525 "lg.ypp"
                { fespacedim=6;}
#line 2657 "lg.tab.cpp"
    break;

  case 90:
#line 526 "lg.ypp"
               { fespacedim=7;}
#line 2663 "lg.tab.cpp"
    break;

  case 91:
#line 527 "lg.ypp"
                { fespacedim = 40;}
#line 2669 "lg.tab.cpp"
    break;

  case 92:
#line 530 "lg.ypp"
                     {fespacecomplex=false;  fespacetype = Find((yyvsp[0].str));}
#line 2675 "lg.tab.cpp"
    break;

  case 93:
#line 531 "lg.ypp"
                                  {
             if ((yyvsp[-1].type) != typevarreal && (yyvsp[-1].type) != typevarcomplex) lgerror (" type of finite element <real> or <complex>");
             fespacecomplex=((yyvsp[-1].type)==typevarcomplex);
             fespacetype = Find((yyvsp[-3].str));}
#line 2684 "lg.tab.cpp"
    break;

  case 94:
#line 536 "lg.ypp"
                                {  (yyval.cexp) = (yyvsp[0].cexp);}
#line 2690 "lg.tab.cpp"
    break;

  case 95:
#line 537 "lg.ypp"
                                             { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2696 "lg.tab.cpp"
    break;

  case 96:
#line 539 "lg.ypp"
                          {  (yyval.cexp) = (yyvsp[0].cexp);}
#line 2702 "lg.tab.cpp"
    break;

  case 97:
#line 540 "lg.ypp"
                                       { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2708 "lg.tab.cpp"
    break;

  case 98:
#line 542 "lg.ypp"
                                                { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[0].cexp);}
#line 2714 "lg.tab.cpp"
    break;

  case 99:
#line 543 "lg.ypp"
                                                { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[0].cexp);}
#line 2720 "lg.tab.cpp"
    break;

  case 100:
#line 546 "lg.ypp"
                               {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,KN<size_t>>((yyvsp[-3].str),typeFESpace((yyvsp[-1].args)),(yyvsp[-1].args),dimFESpaceImage((yyvsp[-1].args)));
     (yyvsp[-1].args).destroy(); }
#line 2727 "lg.tab.cpp"
    break;

  case 101:
#line 549 "lg.ypp"
                            {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,KN<size_t>>((yyvsp[-2].str),typeFESpace((yyvsp[0].args)),(yyvsp[0].args),dimFESpaceImage((yyvsp[0].args)));
     (yyvsp[0].args).destroy(); }
#line 2734 "lg.tab.cpp"
    break;

  case 103:
#line 554 "lg.ypp"
                                                    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2740 "lg.tab.cpp"
    break;

  case 104:
#line 557 "lg.ypp"
                           {dcltype=(yyvsp[0].type);}
#line 2746 "lg.tab.cpp"
    break;

  case 105:
#line 557 "lg.ypp"
                                                          {(yyval.cexp)=(yyvsp[-1].cexp);}
#line 2752 "lg.tab.cpp"
    break;

  case 106:
#line 558 "lg.ypp"
                                                  {(yyval.cexp)=(yyvsp[-1].cexp);}
#line 2758 "lg.tab.cpp"
    break;

  case 107:
#line 559 "lg.ypp"
                           { (yyval.cexp)=(yyvsp[-1].cexp);}
#line 2764 "lg.tab.cpp"
    break;

  case 108:
#line 560 "lg.ypp"
                                        {(yyval.cexp)=currentblock->NewID((yyvsp[-4].type),(yyvsp[-3].str),(yyvsp[-1].cexp));}
#line 2770 "lg.tab.cpp"
    break;

  case 109:
#line 562 "lg.ypp"
                   {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = (yyvsp[-4].type)->right();
                      routineinblock[kkembtype] = currentblock;
                      (yyvsp[-1].routine)=new Routine((yyvsp[-5].type),(yyvsp[-4].type)->right(),(yyvsp[-3].str),(yyvsp[-1].clist_id),currentblock);
		      // routineinblock[kkembtype]->Add($3,"(",$<routine>5); //pas recursif pour l'instanat test  FH 27 dec 2008
                     // cout << " \n after new routine \n " << endl;
                      }
#line 2783 "lg.tab.cpp"
    break;

  case 110:
#line 571 "lg.ypp"
                     { currentblock=(yyvsp[-5].routine)->Set((yyvsp[-1].cinst));
                       currentblock->Add((yyvsp[-7].str),"(",(yyvsp[-5].routine)); //pas recursif pour l'instant test  FH 27 dec 2008
                       kkembtype--;
                       (yyval.cexp)=0;

                        }
#line 2794 "lg.tab.cpp"
    break;

  case 111:
#line 578 "lg.ypp"
                      {Block::open(currentblock); (yyvsp[-4].type)->SetArgs((yyvsp[-1].clist_id));}
#line 2800 "lg.tab.cpp"
    break;

  case 112:
#line 580 "lg.ypp"
                      {  //$<cinst>$=currentblock->close(currentblock);
                         (yyval.cinst).setclose(Block::snewclose(currentblock));// Sep 2016 FH.
                         (yyval.cexp)=currentblock->NewID((yyvsp[-8].type),(yyvsp[-7].str),(yyvsp[-1].cexp),*(yyvsp[-5].clist_id));
                         delete (yyvsp[-5].clist_id); //  FH 23032005
                         }
#line 2810 "lg.tab.cpp"
    break;

  case 113:
#line 587 "lg.ypp"
            {  Block::open(currentblock);}
#line 2816 "lg.tab.cpp"
    break;

  case 114:
#line 588 "lg.ypp"
            { (yyval.endb)=Block::snewclose(currentblock);
//  $$=currentblock->close(currentblock);
}
#line 2824 "lg.tab.cpp"
    break;

  case 115:
#line 592 "lg.ypp"
               {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;}
#line 2831 "lg.tab.cpp"
    break;

  case 116:
#line 594 "lg.ypp"
                   {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;}
#line 2838 "lg.tab.cpp"
    break;

  case 117:
#line 598 "lg.ypp"
                {dcltype=(yyvsp[0].type); Block::open(currentblock);  }
#line 2844 "lg.tab.cpp"
    break;

  case 118:
#line 599 "lg.ypp"
                            {(yyval.cexp)=(yyvsp[0].cexp);}
#line 2850 "lg.tab.cpp"
    break;

  case 119:
#line 601 "lg.ypp"
         { Block::open(currentblock);}
#line 2856 "lg.tab.cpp"
    break;

  case 120:
#line 603 "lg.ypp"
               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));Block::open(currentblock); }
#line 2862 "lg.tab.cpp"
    break;

  case 121:
#line 604 "lg.ypp"
                       { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str)));(yyval.clist_id)->push_back(UnId((yyvsp[0].str)));Block::open(currentblock); }
#line 2868 "lg.tab.cpp"
    break;

  case 122:
#line 605 "lg.ypp"
                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-4].str)));(yyval.clist_id)->push_back(UnId((yyvsp[-2].str)));(yyval.clist_id)->push_back(UnId((yyvsp[0].str)));Block::open(currentblock); }
#line 2874 "lg.tab.cpp"
    break;

  case 123:
#line 607 "lg.ypp"
                   {(yyval.cexp)=0;}
#line 2880 "lg.tab.cpp"
    break;

  case 124:
#line 608 "lg.ypp"
                            {zzzfff->input((yyvsp[0].str));(yyval.cexp)= 0; }
#line 2886 "lg.tab.cpp"
    break;

  case 125:
#line 609 "lg.ypp"
                         {load((yyvsp[0].str));(yyval.cexp)= 0; }
#line 2892 "lg.tab.cpp"
    break;

  case 126:
#line 610 "lg.ypp"
                                             {(yyval.cexp)=Try((yyvsp[-2].cinst),currentblock->close(currentblock,(yyvsp[0].cexp)));}
#line 2898 "lg.tab.cpp"
    break;

  case 127:
#line 611 "lg.ypp"
                     {(yyval.cexp)=(yyvsp[-1].cexp);}
#line 2904 "lg.tab.cpp"
    break;

  case 128:
#line 612 "lg.ypp"
                         {(yyval.cexp)=(yyvsp[0].cexp);}
#line 2910 "lg.tab.cpp"
    break;

  case 129:
#line 614 "lg.ypp"
                {(yyvsp[-1].cexp)=ForAll(currentblock,(yyvsp[-3].clist_id),(yyvsp[-1].cexp));}
#line 2916 "lg.tab.cpp"
    break;

  case 130:
#line 615 "lg.ypp"
                         {
                    inloopcount--;
                    (yyval.cexp)=Block::close(currentblock,C_F0(ForAll((yyvsp[-3].cexp),(yyvsp[0].cexp))));
                 }
#line 2925 "lg.tab.cpp"
    break;

  case 131:
#line 619 "lg.ypp"
                                                                 {inloopcount--; (yyval.cexp)=For((yyvsp[-6].cexp),(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2931 "lg.tab.cpp"
    break;

  case 132:
#line 621 "lg.ypp"
                {inloopcount--;
                 (yyval.cexp)=Block::close(currentblock,C_F0(For((yyvsp[-6].cexp),(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp))));
                }
#line 2939 "lg.tab.cpp"
    break;

  case 133:
#line 624 "lg.ypp"
                                                {inloopcount--;(yyval.cexp)=While((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2945 "lg.tab.cpp"
    break;

  case 134:
#line 625 "lg.ypp"
                                           {(yyval.cexp)=FIf((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2951 "lg.tab.cpp"
    break;

  case 135:
#line 626 "lg.ypp"
                                                            {(yyval.cexp)=FIf((yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2957 "lg.tab.cpp"
    break;

  case 136:
#line 627 "lg.ypp"
                                    { /* [[begin:]] [[end:]] */
             (yyvsp[-1].cinst).setclose((yyvsp[0].endb));
             (yyval.cexp)=(yyvsp[-1].cinst);
                    //  $$=C_F0(new E_block($2,$3),atype<void>());
         }
#line 2967 "lg.tab.cpp"
    break;

  case 137:
#line 632 "lg.ypp"
                                     { /* <<BORDER_ID>> */
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[-1].str),C_F0(TheOperators,"[border]",(yyvsp[0].args)));}
#line 2974 "lg.tab.cpp"
    break;

  case 138:
#line 634 "lg.ypp"
                                           {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[-4].str),C_F0(TheOperators,"[border]",(yyvsp[-2].args)));}
#line 2981 "lg.tab.cpp"
    break;

  case 139:
#line 637 "lg.ypp"
                      {
                    if(inloopcount)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_break),atype<void>());
                    else lgerror("break not in loop");}
#line 2990 "lg.tab.cpp"
    break;

  case 140:
#line 641 "lg.ypp"
                         {
                    if(inloopcount)
                        (yyval.cexp)= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");}
#line 2999 "lg.tab.cpp"
    break;

  case 141:
#line 645 "lg.ypp"
                             {
                    if (kkembtype>=0)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_return,(rettype[kkembtype]->CastTo((yyvsp[-1].cexp))).OnReturn()) ,atype<void>());
                     else lgerror(" return not in routine ");}
#line 3008 "lg.tab.cpp"
    break;

  case 142:
#line 652 "lg.ypp"
                                         {(yyval.cexp) =  (yyvsp[0].cexp); }
#line 3014 "lg.tab.cpp"
    break;

  case 143:
#line 655 "lg.ypp"
                                     {
   Block::open(currentblock);
   (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[-5].str),atype<double*>());
   (yyval.args)+= (yyvsp[-3].cexp);
   (yyval.args)+= (yyvsp[-1].cexp);
   (yyval.args)+= currentblock->NewVar<LocalVariable>("IndexBorder",atype<long*>());}
#line 3025 "lg.tab.cpp"
    break;

  case 144:
#line 663 "lg.ypp"
                                            {
    Block::open(currentblock);
    (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[-7].str),atype<double*>());
    (yyval.args)+= (yyvsp[-5].cexp);
    (yyval.args)+= (yyvsp[-3].cexp);
    (yyval.args)+= currentblock->NewVar<LocalVariable>((yyvsp[-1].str),atype<long*>());}
#line 3036 "lg.tab.cpp"
    break;

  case 145:
#line 671 "lg.ypp"
                                  {
    //currentblock->close(currentblock;);
   (yyval.args) = ((yyvsp[-1].args) += currentblock->close(currentblock,(yyvsp[0].cexp)));
   }
#line 3045 "lg.tab.cpp"
    break;

  case 147:
#line 679 "lg.ypp"
                  {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3051 "lg.tab.cpp"
    break;

  case 154:
#line 692 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3057 "lg.tab.cpp"
    break;

  case 155:
#line 693 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"+=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3063 "lg.tab.cpp"
    break;

  case 156:
#line 694 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"-=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3069 "lg.tab.cpp"
    break;

  case 157:
#line 695 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"*=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3075 "lg.tab.cpp"
    break;

  case 158:
#line 696 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"/=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3081 "lg.tab.cpp"
    break;

  case 159:
#line 697 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,".*=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3087 "lg.tab.cpp"
    break;

  case 160:
#line 698 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"./=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3093 "lg.tab.cpp"
    break;

  case 162:
#line 704 "lg.ypp"
                                                            {(yyval.cexp)=C_F0(TheOperators,"?:",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3099 "lg.tab.cpp"
    break;

  case 163:
#line 705 "lg.ypp"
                                                            {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3105 "lg.tab.cpp"
    break;

  case 164:
#line 706 "lg.ypp"
                                                            {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3111 "lg.tab.cpp"
    break;

  case 166:
#line 711 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3117 "lg.tab.cpp"
    break;

  case 167:
#line 712 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3123 "lg.tab.cpp"
    break;

  case 168:
#line 713 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3129 "lg.tab.cpp"
    break;

  case 169:
#line 714 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3135 "lg.tab.cpp"
    break;

  case 170:
#line 715 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3141 "lg.tab.cpp"
    break;

  case 171:
#line 716 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3147 "lg.tab.cpp"
    break;

  case 172:
#line 717 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3153 "lg.tab.cpp"
    break;

  case 173:
#line 718 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3159 "lg.tab.cpp"
    break;

  case 174:
#line 719 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3165 "lg.tab.cpp"
    break;

  case 175:
#line 720 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3171 "lg.tab.cpp"
    break;

  case 176:
#line 721 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3177 "lg.tab.cpp"
    break;

  case 177:
#line 722 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3183 "lg.tab.cpp"
    break;

  case 178:
#line 723 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3189 "lg.tab.cpp"
    break;

  case 179:
#line 724 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3195 "lg.tab.cpp"
    break;

  case 180:
#line 725 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3201 "lg.tab.cpp"
    break;

  case 181:
#line 726 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3207 "lg.tab.cpp"
    break;

  case 182:
#line 727 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3213 "lg.tab.cpp"
    break;

  case 183:
#line 728 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3219 "lg.tab.cpp"
    break;

  case 184:
#line 729 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3225 "lg.tab.cpp"
    break;

  case 185:
#line 733 "lg.ypp"
                                                                    {(yyval.cexp)=(yyvsp[0].cexp);}
#line 3231 "lg.tab.cpp"
    break;

  case 186:
#line 734 "lg.ypp"
                                                                {(yyval.cexp)=C_F0(TheOperators,":");}
#line 3237 "lg.tab.cpp"
    break;

  case 187:
#line 735 "lg.ypp"
                                                                {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3243 "lg.tab.cpp"
    break;

  case 188:
#line 736 "lg.ypp"
                                                                {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3249 "lg.tab.cpp"
    break;

  case 189:
#line 741 "lg.ypp"
                                 {
      (yyval.args) = (yyvsp[-2].cexp);
      (yyval.args) += (yyvsp[0].args); }
#line 3257 "lg.tab.cpp"
    break;

  case 190:
#line 745 "lg.ypp"
                                                            {(yyval.args) = 0;}
#line 3263 "lg.tab.cpp"
    break;

  case 191:
#line 746 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3269 "lg.tab.cpp"
    break;

  case 192:
#line 747 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3275 "lg.tab.cpp"
    break;

  case 193:
#line 748 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3281 "lg.tab.cpp"
    break;

  case 194:
#line 749 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3287 "lg.tab.cpp"
    break;

  case 195:
#line 750 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3293 "lg.tab.cpp"
    break;

  case 196:
#line 751 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3299 "lg.tab.cpp"
    break;

  case 197:
#line 752 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3305 "lg.tab.cpp"
    break;

  case 198:
#line 753 "lg.ypp"
                                                            {(yyval.args) = make_pair<const char *,const C_F0>((const char *) (yyvsp[-2].str),(C_F0) (yyvsp[0].cexp));}
#line 3311 "lg.tab.cpp"
    break;

  case 199:
#line 754 "lg.ypp"
                                                            {(yyval.args) = (yyvsp[0].cexp);}
#line 3317 "lg.tab.cpp"
    break;

  case 200:
#line 755 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3323 "lg.tab.cpp"
    break;

  case 201:
#line 756 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3329 "lg.tab.cpp"
    break;

  case 202:
#line 757 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3335 "lg.tab.cpp"
    break;

  case 203:
#line 758 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3341 "lg.tab.cpp"
    break;

  case 204:
#line 759 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3347 "lg.tab.cpp"
    break;

  case 205:
#line 760 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3353 "lg.tab.cpp"
    break;

  case 206:
#line 761 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3359 "lg.tab.cpp"
    break;

  case 207:
#line 762 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp));}
#line 3365 "lg.tab.cpp"
    break;

  case 208:
#line 765 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-4].args)+= make_pair<const char *,const C_F0>((const char *)(yyvsp[-2].str),(C_F0) (yyvsp[0].cexp)));}
#line 3371 "lg.tab.cpp"
    break;

  case 209:
#line 768 "lg.ypp"
                                 {(yyval.args)=(yyvsp[0].cexp);}
#line 3377 "lg.tab.cpp"
    break;

  case 210:
#line 769 "lg.ypp"
                                 {(yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp));}
#line 3383 "lg.tab.cpp"
    break;

  case 211:
#line 772 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3389 "lg.tab.cpp"
    break;

  case 212:
#line 773 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3395 "lg.tab.cpp"
    break;

  case 213:
#line 774 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3401 "lg.tab.cpp"
    break;

  case 214:
#line 775 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3407 "lg.tab.cpp"
    break;

  case 215:
#line 776 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3413 "lg.tab.cpp"
    break;

  case 216:
#line 777 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3419 "lg.tab.cpp"
    break;

  case 217:
#line 778 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3425 "lg.tab.cpp"
    break;

  case 218:
#line 779 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3431 "lg.tab.cpp"
    break;

  case 219:
#line 780 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3437 "lg.tab.cpp"
    break;

  case 220:
#line 781 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3443 "lg.tab.cpp"
    break;

  case 221:
#line 782 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3449 "lg.tab.cpp"
    break;

  case 222:
#line 783 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3455 "lg.tab.cpp"
    break;

  case 223:
#line 784 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3461 "lg.tab.cpp"
    break;

  case 224:
#line 785 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3467 "lg.tab.cpp"
    break;

  case 225:
#line 786 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3473 "lg.tab.cpp"
    break;

  case 226:
#line 787 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3479 "lg.tab.cpp"
    break;

  case 227:
#line 788 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3485 "lg.tab.cpp"
    break;

  case 228:
#line 789 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3491 "lg.tab.cpp"
    break;

  case 230:
#line 794 "lg.ypp"
                              {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[0].cexp));}
#line 3497 "lg.tab.cpp"
    break;

  case 232:
#line 799 "lg.ypp"
                                    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3503 "lg.tab.cpp"
    break;

  case 233:
#line 800 "lg.ypp"
                                    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3509 "lg.tab.cpp"
    break;

  case 235:
#line 805 "lg.ypp"
                                   {(yyval.cexp)=C_F0(TheOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
#line 3515 "lg.tab.cpp"
    break;

  case 236:
#line 813 "lg.ypp"
                        {(yyval.cexp)=Find((yyvsp[0].str));}
#line 3521 "lg.tab.cpp"
    break;

  case 237:
#line 817 "lg.ypp"
                        {(yyval.cexp)= CConstant((yyvsp[0].lnum));}
#line 3527 "lg.tab.cpp"
    break;

  case 238:
#line 818 "lg.ypp"
                        {(yyval.cexp)= CConstant((yyvsp[0].dnum));}
#line 3533 "lg.tab.cpp"
    break;

  case 239:
#line 819 "lg.ypp"
                        {(yyval.cexp)= CConstant(complex<double>(0,(yyvsp[0].dnum)));}
#line 3539 "lg.tab.cpp"
    break;

  case 240:
#line 820 "lg.ypp"
                  {(yyval.cexp)= CConstant<const char *>((yyvsp[0].str));}
#line 3545 "lg.tab.cpp"
    break;

  case 241:
#line 825 "lg.ypp"
                                                            {(yyval.cexp)=C_F0((yyvsp[-3].cexp),(yyvsp[-2].oper),(yyvsp[-1].args));}
#line 3551 "lg.tab.cpp"
    break;

  case 242:
#line 827 "lg.ypp"
                                              {(yyval.cexp)=C_F0((yyvsp[-3].cexp),(yyvsp[-2].oper),(yyvsp[-1].cexp));}
#line 3557 "lg.tab.cpp"
    break;

  case 243:
#line 828 "lg.ypp"
                                                                {(yyval.cexp)=C_F0((yyvsp[-5].cexp),(yyvsp[-4].oper),(yyvsp[-3].cexp),(yyvsp[-1].cexp));}
#line 3563 "lg.tab.cpp"
    break;

  case 244:
#line 829 "lg.ypp"
                                   {(yyval.cexp)=C_F0((yyvsp[-2].cexp),"[]");}
#line 3569 "lg.tab.cpp"
    break;

  case 245:
#line 830 "lg.ypp"
                                 { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].str)) ;}
#line 3575 "lg.tab.cpp"
    break;

  case 246:
#line 831 "lg.ypp"
                                 { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3581 "lg.tab.cpp"
    break;

  case 247:
#line 832 "lg.ypp"
                                          { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3587 "lg.tab.cpp"
    break;

  case 248:
#line 833 "lg.ypp"
                                  { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3593 "lg.tab.cpp"
    break;

  case 249:
#line 834 "lg.ypp"
                                           { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3599 "lg.tab.cpp"
    break;

  case 250:
#line 835 "lg.ypp"
                                  { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3605 "lg.tab.cpp"
    break;

  case 251:
#line 836 "lg.ypp"
                                           { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3611 "lg.tab.cpp"
    break;

  case 252:
#line 837 "lg.ypp"
                                  { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3617 "lg.tab.cpp"
    break;

  case 253:
#line 838 "lg.ypp"
                                           { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3623 "lg.tab.cpp"
    break;

  case 254:
#line 839 "lg.ypp"
                                  { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3629 "lg.tab.cpp"
    break;

  case 255:
#line 840 "lg.ypp"
                                           { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3635 "lg.tab.cpp"
    break;

  case 256:
#line 841 "lg.ypp"
                                   { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3641 "lg.tab.cpp"
    break;

  case 257:
#line 842 "lg.ypp"
                                            { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3647 "lg.tab.cpp"
    break;

  case 258:
#line 843 "lg.ypp"
                                 {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
#line 3653 "lg.tab.cpp"
    break;

  case 259:
#line 844 "lg.ypp"
                                 {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
#line 3659 "lg.tab.cpp"
    break;

  case 260:
#line 845 "lg.ypp"
                                         {
             if ((yyvsp[-3].type)->right()->CastingFrom((yyvsp[-1].cexp).left()) )
                (yyval.cexp)=(yyvsp[-3].type)->right()->CastTo((yyvsp[-1].cexp))  ;
             else { (yyval.cexp)=(yyvsp[-3].type)->right()->Find("<--",basicAC_F0_wa((yyvsp[-1].cexp)));
             if (!(yyval.cexp).left()) { cerr << " no wait to change " << (yyvsp[-1].cexp).left()->right()->name() << " in " <<
                                        (yyvsp[-3].type)->right()->name() << endl;
                                CompileError(" Error in type(exp) "); }
             }
            }
#line 3673 "lg.tab.cpp"
    break;

  case 261:
#line 854 "lg.ypp"
                                        {
           { (yyval.cexp)=(yyvsp[-3].type)->right()->Find("<--",basicAC_F0_wa((yyvsp[-1].args)));
           if (!(yyval.cexp).left()) { cerr << " no wait to change (args) in " <<
                                      (yyvsp[-3].type)->right()->name() << endl;
                              CompileError(" Error in type(exp) "); }
           }
          }
#line 3685 "lg.tab.cpp"
    break;

  case 262:
#line 862 "lg.ypp"
                        {(yyval.cexp)=(yyvsp[-1].cexp);}
#line 3691 "lg.tab.cpp"
    break;

  case 263:
#line 863 "lg.ypp"
                          { (yyval.cexp)=C_F0(TheOperators,"[]",(yyvsp[-1].args));}
#line 3697 "lg.tab.cpp"
    break;

  case 264:
#line 864 "lg.ypp"
                           { (yyval.cexp)=C_F0(TheOperators,"<>",(yyvsp[-1].args));}
#line 3703 "lg.tab.cpp"
    break;

  case 265:
#line 865 "lg.ypp"
                   { (yyval.cexp)=C_F0(TheOperators,"<>",(yyvsp[0].args));}
#line 3709 "lg.tab.cpp"
    break;


#line 3713 "lg.tab.cpp"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = YY_CAST (char *, YYSTACK_ALLOC (YY_CAST (YYSIZE_T, yymsg_alloc)));
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
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

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;

  /* Do not reclaim the symbols of the rule whose action triggered
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
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
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

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


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


#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif


/*-----------------------------------------------------.
| yyreturn -- parsing is finished, return the result.  |
`-----------------------------------------------------*/
yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[+*yyssp], yyvsp);
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
  return yyresult;
}
#line 869 "lg.ypp"



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
void msh3_Load_Init( ); //

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
      if( typeofscript ) cout << " " << " Markdown = " << typeofscript;
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
   msh3_Load_Init(); // Add msh3 lib !!! 
     //  callInitsFunct() ; //  init for dynamique libs ...

   if(verbosity>2 || ((mpirank==0)&& verbosity)  )  cout << endl;
  zzzfff->input(cc,typeofscript); // [[file:../fflib/lex.cpp::mylex_input_filename]]
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
