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
    FESPACED3 = 298,
    FESPACES = 299,
    FESPACEDS = 300,
    FESPACEL = 301,
    FESPACEDL = 302,
    VGFESPACE = 303,
    GFESPACE = 304,
    PLUSEQ = 305,
    MOINSEQ = 306,
    MULEQ = 307,
    DIVEQ = 308,
    DOTMULEQ = 309,
    DOTDIVEQ = 310,
    ARROW = 311,
    BORDER = 312,
    SOLVE = 313
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

#line 386 "lg.tab.cpp"

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
#define YYFINAL  124
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1863

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  84
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  50
/* YYNRULES -- Number of rules.  */
#define YYNRULES  291
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  564

#define YYUNDEFTOK  2
#define YYMAXUTOK   313


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
       2,     2,     2,     2,     2,     2,     2,     2,    82,    79,
      16,     6,    17,    83,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    35,     2,    38,    31,    33,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    80,    11,    81,     2,     2,     2,     2,
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
      75,    76,    77,    78
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   342,   342,   412,   416,   417,   423,   424,   425,   426,
     427,   428,   429,   430,   431,   432,   433,   434,   435,   436,
     437,   438,   439,   440,   441,   442,   443,   444,   445,   446,
     447,   448,   449,   450,   451,   452,   453,   454,   455,   456,
     457,   458,   459,   460,   463,   464,   469,   469,   469,   469,
     469,   469,   469,   469,   469,   472,   473,   474,   475,   481,
     482,   483,   484,   485,   486,   487,   488,   489,   490,   491,
     492,   493,   494,   495,   496,   497,   498,   499,   504,   505,
     506,   507,   508,   509,   510,   511,   515,   516,   517,   518,
     519,   520,   524,   525,   530,   531,   532,   533,   534,   535,
     536,   537,   538,   539,   542,   543,   548,   549,   551,   552,
     554,   555,   558,   561,   565,   566,   569,   569,   570,   571,
     572,   574,   573,   590,   589,   599,   600,   604,   606,   610,
     610,   613,   615,   616,   617,   619,   620,   621,   622,   623,
     624,   626,   625,   631,   632,   636,   637,   638,   639,   644,
     646,   649,   653,   657,   664,   667,   675,   683,   690,   691,
     695,   696,   697,   698,   699,   703,   704,   705,   706,   707,
     708,   709,   710,   715,   716,   717,   718,   722,   723,   724,
     725,   726,   727,   728,   729,   730,   731,   732,   733,   734,
     735,   736,   737,   738,   739,   740,   741,   745,   746,   747,
     748,   753,   757,   758,   759,   760,   761,   762,   763,   764,
     765,   766,   767,   768,   769,   770,   771,   772,   773,   774,
     775,   776,   777,   778,   781,   784,   785,   788,   789,   790,
     791,   792,   793,   794,   795,   796,   797,   798,   799,   800,
     801,   802,   803,   804,   805,   806,   807,   808,   809,   810,
     811,   815,   816,   820,   821,   822,   826,   827,   835,   839,
     840,   841,   842,   847,   849,   850,   851,   852,   853,   854,
     855,   856,   857,   858,   859,   860,   861,   862,   863,   864,
     865,   866,   867,   868,   869,   870,   871,   880,   888,   889,
     890,   891
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
  "FUNCTION", "FESPACE", "FESPACE1", "FESPACE3", "FESPACED3", "FESPACES",
  "FESPACEDS", "FESPACEL", "FESPACEDL", "VGFESPACE", "GFESPACE", "PLUSEQ",
  "MOINSEQ", "MULEQ", "DIVEQ", "DOTMULEQ", "DOTDIVEQ", "ARROW", "BORDER",
  "SOLVE", "';'", "'{'", "'}'", "':'", "'?'", "$accept", "start", "input",
  "instructions", "list_of_id_args", "list_of_id1", "id", "list_of_dcls",
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
     305,   306,   307,   308,   309,   310,   311,   312,   313,    59,
     123,   125,    58,    63
};
# endif

#define YYPACT_NINF (-318)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-222)

#define yytable_value_is_error(Yyn) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     651,   -20,   972,  -318,  -318,  -318,  -318,  -318,  1367,  1367,
    -318,  -318,  -318,  -318,  1721,  -318,   -13,     3,  -318,  -318,
     -51,    30,  1367,  -318,   308,  1694,  1356,  1408,  1539,  1548,
    1582,  1591,  1625,  1634,  -318,  -318,  1721,  -318,  -318,    82,
     105,   651,  -318,   145,   203,    95,  -318,   651,   101,   143,
     115,  -318,    11,  1505,  -318,    75,   412,   236,  -318,  -318,
     213,   849,  1367,  -318,  -318,  -318,  -318,  -318,  -318,  -318,
    -318,   356,   171,   214,   262,   315,   324,   346,   350,   358,
     381,    13,  -318,    15,  -318,  -318,  -318,  -318,  -318,  -318,
    -318,  -318,  -318,    20,  -318,    28,  -318,  -318,  -318,  -318,
      31,   150,  1203,   182,   313,    99,  1721,  1038,  1721,  1038,
    1721,  1038,  1721,  1038,  1721,  1038,  1721,  1038,  1721,  1038,
    1721,  1038,  1721,   304,  -318,  -318,  -318,  1721,   204,  1704,
     151,  -318,   270,  -318,   489,  1419,  1721,  1367,   651,  1367,
    -318,  -318,  1367,  1367,  1367,  1367,  1367,  1367,  1367,  1367,
    1367,  1367,  1367,  1367,  1367,  1367,  1367,  1367,  1367,  1367,
    1367,  1367,  1367,  1367,  1367,  1367,  1367,  1367,  1367,  1367,
    1216,  1712,  1367,  1367,  -318,  -318,  -318,  1038,  1093,  1721,
      92,  -318,  -318,  1367,  -318,  1471,  1471,  1721,  -318,  -318,
     288,  -318,   724,    94,   247,    47,  1367,  1667,   274,   326,
     157,   223,   246,   267,   294,   328,   361,   365,  -318,   329,
    -318,   129,  -318,   147,  -318,   148,  -318,   159,  -318,   164,
    -318,   174,  -318,   177,  -318,   187,  -318,  1721,  1367,   651,
    -318,   153,    33,   309,   314,    51,  -318,  1367,  1367,   754,
    -318,  -318,  -318,   289,    35,   341,   271,   234,   543,  -318,
    -318,  -318,  -318,  -318,  -318,  -318,  -318,  1798,  1798,  1813,
    1813,  1826,  1826,  1837,  1837,   701,   701,   701,   701,   383,
     383,  -318,  -318,  -318,  -318,  -318,   743,   762,  -318,  -318,
    -318,  -318,  -318,  -318,  -318,  -318,  -318,  -318,  -318,  -318,
    -318,  -318,  -318,  -318,  -318,  -318,   248,  -318,    79,  -318,
     651,  -318,   170,   830,   864,   874,   908,   942,   952,   986,
     366,  -318,  -318,   249,  -318,   340,  1367,  1038,  -318,  -318,
     318,   353,    37,   368,  1667,   155,   399,  1027,  1082,  1137,
    1192,  1233,   269,   394,  1247,  1667,  1367,  1148,  -318,  -318,
    -318,  -318,  -318,  -318,  -318,  -318,   405,    81,  -318,  1367,
    1471,  1721,  -318,  -318,   766,  1721,   175,  -318,   410,  1721,
    -318,  1721,  1367,  1367,  1721,  1505,   651,   363,  1367,  1367,
    -318,  1203,  -318,   439,  -318,  -318,  -318,  -318,  -318,  -318,
    -318,  -318,  1367,  1471,  -318,   392,   837,   446,   419,   409,
    -318,  1667,    86,  1721,  -318,  1721,  -318,  1721,  -318,  1721,
    -318,  1721,  -318,  1721,  -318,  1721,  -318,  1677,  -318,  1367,
    1721,  -318,   272,  -318,   408,   434,   441,   804,   868,   881,
     915,   927,  -318,   466,  -318,  1367,   397,  -318,   282,  -318,
    1721,   447,  -318,   468,  -318,  1367,  1367,  -318,   478,    41,
      44,   481,   764,  -318,   454,  -318,  1781,  1781,   452,   651,
    -318,  -318,    88,  1367,   462,   460,   106,  -318,  -318,  -318,
    -318,  -318,  -318,  -318,  -318,   464,  1667,  1259,  1273,  1289,
    1303,   495,  1315,   500,  -318,  -318,  -318,  1367,   507,  -318,
    -318,   116,  1367,   766,  -318,   475,  1367,  1367,  1721,  -318,
     479,  -318,  -318,   456,  -318,  1781,   458,  -318,   503,  1667,
     120,  1721,  -318,  1721,  -318,  1721,  -318,  1721,  -318,  1367,
    1721,  -318,  1367,   445,  -318,  1367,   487,   490,  -318,  -318,
     317,   320,  -318,   651,   497,   488,   522,  -318,   122,   518,
    -318,  -318,  -318,  -318,  -318,  -318,   482,   651,   -28,  1367,
    -318,   651,   651,  -318,   509,  -318,  -318,   545,  -318,  -318,
     597,  -318,  1721,   528,  -318,  -318,   530,  -318,  -318,   537,
    -318,   651,  -318,  -318
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int16 yydefact[] =
{
       0,     0,     0,   161,   160,   163,   164,   162,     0,     0,
     259,   260,   261,   258,     0,   262,     0,     0,   127,   128,
       0,     0,     0,   131,    78,     0,   227,   228,   229,   234,
     230,   233,   231,   232,    99,   100,     0,   135,   125,     0,
       0,     3,   116,   104,     0,     0,   140,     0,     0,     0,
       0,     4,     0,     0,   158,   165,   173,   291,   177,   251,
     253,   256,     0,   227,   228,   229,   234,   230,   233,   231,
     232,     0,     0,   227,   228,   229,   234,   230,   233,   231,
     232,     0,   225,     0,    46,    47,    51,    48,    53,    49,
      52,    50,    54,     0,   114,     0,   136,   137,   151,   152,
       0,     0,     0,     0,    78,     0,     0,   202,     0,   202,
       0,   202,     0,   202,     0,   202,     0,   202,     0,   202,
       0,   202,     0,     0,     1,     2,     5,     0,     0,     0,
      86,   108,   110,   119,     0,     0,     0,     0,     0,     0,
     139,   252,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   284,   285,   257,   202,     0,     0,
       0,   290,   288,     0,   289,     0,     0,     0,   118,   153,
       0,   198,   197,     0,     0,     0,     0,     6,     0,   258,
     227,   228,   229,   234,   230,   233,   231,   232,   211,     0,
     213,     0,   268,     0,   270,     0,   272,     0,   280,     0,
     274,     0,   278,     0,   276,     0,   282,     0,     0,     0,
     149,    55,     0,     0,     0,     0,    44,     0,     0,     0,
     126,   148,   129,     0,     0,   132,     0,     0,     0,   159,
     166,   167,   168,   169,   170,   171,   172,   185,   186,   190,
     189,   188,   187,   195,   196,   191,   193,   192,   194,   183,
     184,   178,   181,   182,   179,   180,   175,     0,   235,   236,
     237,   242,   238,   241,   239,   240,   243,   244,   245,   249,
     246,   248,   247,   250,   254,   255,     0,   266,     0,   267,
       0,   226,   227,   228,   229,   234,   230,   233,   231,   232,
       0,   113,    59,     0,   115,    83,     0,   202,   286,   287,
       0,    79,     0,     0,     6,    47,    48,    53,    49,    52,
      50,    54,     0,     7,     0,     6,     0,     0,   269,   271,
     273,   281,   275,   279,   277,   283,     0,     0,   157,     0,
       0,     0,   117,   105,     0,     0,    89,    88,     0,     0,
     109,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     263,     0,   264,   146,    60,    61,    62,    72,    63,    70,
      64,    74,     0,     0,   112,     0,   199,   201,     0,     0,
     120,     6,     0,     0,     9,     0,    11,     0,    13,     0,
      15,     0,    17,     0,    19,     0,    21,     0,   123,     0,
       0,    23,     0,   212,   227,   228,   229,   234,   230,   233,
     231,   232,   222,     0,   223,     0,     0,    56,     0,    58,
       0,     0,   106,   111,    45,     0,     0,    87,   130,     0,
       0,   133,     0,   145,     0,   138,   176,   174,     0,     0,
      76,    77,     0,     0,    81,     0,     0,    25,    10,    12,
      14,    16,    18,    20,    22,     0,     6,    47,    48,    49,
      50,    28,     0,     0,     8,    24,   121,     0,     0,   150,
      57,     0,     0,     0,    91,     0,     0,     0,     0,   141,
       0,   265,   147,     0,    84,   200,     0,    80,    27,     6,
       0,     0,    34,     0,    36,     0,    38,     0,    40,     0,
       0,    42,     0,     0,   224,     0,     0,     0,   107,    90,
       0,     0,   134,     0,     0,     0,     0,    26,     0,    29,
      35,    37,    39,    41,    33,    43,     0,     0,   159,     0,
      92,     0,     0,   142,     0,    85,    82,    32,    30,   124,
       0,   155,     0,     0,   144,   143,     0,    31,   122,     0,
      93,     0,   156,   154
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -318,  -318,  -318,   -37,  -311,   146,   -14,  -317,  -183,   -17,
     336,    96,  -318,  -318,  -318,  -318,  -318,   393,  -318,  -318,
    -318,  -318,  -318,  -318,  -318,  -318,  -318,  -318,  -318,  -318,
    -318,   -39,  -318,  -318,  -318,  -318,    -7,  -318,    -5,  -180,
     -90,   -95,  -318,   -74,   359,   579,    26,   536,  -318,   227
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    39,    40,    41,   332,   235,   209,   232,   311,    42,
     131,   432,    43,    44,   433,   132,    45,    94,    95,    46,
     127,   513,   473,    47,   241,    48,    49,   243,   361,    50,
     246,    51,   523,   445,   229,   230,    52,    53,    54,    55,
      56,   210,   194,   211,    83,    57,    58,    59,    60,    61
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      93,    81,   126,   313,    82,   312,   312,   193,   106,   551,
     134,   105,   192,   392,    62,   100,   139,   192,   139,   192,
     183,   192,   123,   192,   412,   192,   185,   192,    98,   192,
     130,   192,    96,   187,   429,   213,   139,   215,   351,   217,
     139,   219,   139,   221,   438,   223,   139,   225,    97,   139,
     182,   552,   320,   184,   186,   180,   355,   357,   358,   257,
     258,   259,   260,   261,   262,   263,   264,   265,   266,   267,
     268,   269,   270,   271,   272,   273,   274,   275,   276,   277,
     456,   142,   124,   298,   371,   321,   183,   192,   192,   356,
     140,   407,   198,   493,   212,   126,   214,   139,   216,   317,
     218,   248,   220,   296,   222,   196,   224,   188,   226,    99,
     189,   407,   352,   231,   363,   236,   390,   372,   242,   426,
     486,   355,   245,   487,   457,   407,   494,   407,   244,   300,
     247,   318,   249,   197,   337,   135,   136,   250,   251,   252,
     253,   254,   255,   256,   498,   143,   144,   145,   146,   147,
     148,   125,   337,   337,   516,   500,   413,   237,   529,   349,
     547,   128,  -203,   -47,   337,   299,   338,   428,   393,   337,
     312,   310,   310,    93,   133,   -65,   -47,   137,   301,   337,
     334,   435,   337,   333,   339,   340,   238,   350,   528,   322,
     348,   107,   337,   108,  -203,   138,   341,    84,   294,   295,
     451,   342,   450,   312,   107,   102,   108,   -65,   190,   126,
     436,   343,    84,   346,   344,    85,    86,    87,    88,    89,
      90,    91,    92,    82,   345,   130,   386,   192,  -204,   -51,
      85,    86,    87,    88,    89,    90,    91,    92,   129,   139,
     195,   170,   424,   387,   172,    84,   173,   192,   107,   -65,
     108,  -205,   -48,   337,   383,   484,   485,   109,   171,   110,
    -204,   373,   233,    85,    86,    87,    88,    89,    90,    91,
      92,   366,  -210,   -53,   407,   239,   448,   407,   446,   447,
     111,   192,   112,  -205,   319,   370,   384,   383,   374,   375,
     376,   377,   378,   379,   380,   381,   109,   514,   110,  -206,
     -49,   113,   517,   114,  -210,   315,   408,   334,   335,   476,
     333,   394,   396,   398,   400,   402,   404,   406,   334,   480,
     411,   333,   139,   423,   101,   139,   353,   443,   115,   101,
     116,  -206,   -46,  -209,   -52,   336,   310,   231,   227,   228,
     431,   434,   102,   103,   427,   236,   364,   231,   103,   111,
     441,   112,   354,   365,   541,   439,   440,   542,   113,   553,
     114,   170,   117,   495,   118,  -209,  -207,   -50,   362,   310,
    -208,   -54,   382,   181,   334,   385,   388,   333,   171,   458,
     115,   459,   116,   460,   117,   461,   118,   462,   389,   463,
     472,   464,   119,   471,   120,   119,   475,   120,  -207,   121,
     409,   122,  -208,   391,   474,   163,   164,   165,   166,   167,
     492,   425,   395,  -214,   -47,   121,   236,   122,   478,   444,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,   163,   164,   165,   166,   167,  -215,
     -51,    84,   107,   449,   108,  -214,  -216,   -48,   437,   334,
     452,   337,   333,   502,   504,   506,   508,   454,   511,    85,
      86,    87,    88,    89,    90,    91,    92,   455,   109,   431,
     110,  -215,   477,   483,   522,   111,   479,   112,  -216,   520,
     521,   482,   334,   351,   543,   333,   488,   530,   490,   531,
     491,   532,     1,   533,   168,   169,   535,   496,   497,   499,
     550,   509,   554,   555,   534,     2,   512,   536,   538,     3,
       4,   126,   515,   519,   525,   524,   526,     5,     6,     7,
     527,   539,   563,     8,     9,   537,   545,   540,    10,    11,
      12,    13,    14,   544,    15,   548,    16,    17,   559,    18,
      19,    20,    21,    22,    23,   556,     1,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,     2,
     546,   549,   557,     3,     4,   560,    36,   561,    37,    38,
     240,     5,     6,     7,   562,   360,   481,     8,     9,   518,
     314,    71,    10,    11,    12,    13,    14,   347,    15,   141,
      16,    17,   442,    18,    19,    20,    21,    22,    23,     0,
       1,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      33,    34,    35,     2,     0,     0,     0,     3,     4,     0,
      36,     0,    37,    38,   367,     5,     6,     7,     0,     0,
       0,     8,     9,     0,     0,     0,    10,    11,    12,    13,
      14,     0,    15,     0,    16,    17,     0,    18,    19,    20,
      21,    22,    23,     0,     1,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,     2,     0,     0,
       0,     3,     4,     0,    36,     0,    37,    38,   558,     5,
       6,     7,     0,     0,     0,     8,     9,     0,     0,     0,
      10,    11,    12,    13,    14,     0,    15,     0,    16,    17,
       0,    18,    19,    20,    21,    22,    23,     0,     0,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,   161,   162,   163,   164,   165,   166,   167,    36,     0,
      37,    38,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   160,   161,   162,   163,   164,   165,   166,
     167,   149,   150,   151,   152,   153,   154,   155,   156,   157,
     158,   159,   160,   161,   162,   163,   164,   165,   166,   167,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,   163,   164,   165,   166,   167,   359,
       0,     0,   174,   175,     0,     0,    84,     0,   177,   178,
     179,   430,   489,     0,     0,     0,   316,     0,    84,  -221,
     -53,     0,     0,     0,    85,    86,    87,    88,    89,    90,
      91,    92,     0,     0,     0,   368,    85,    86,    87,    88,
      89,    90,    91,    92,     0,   -66,   -51,     0,   113,     0,
     114,  -221,     0,     0,   369,   149,   150,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   160,   161,   162,   163,
     164,   165,   166,   167,   109,     0,   110,   -66,     0,   -67,
     -48,     0,    84,  -217,   -49,     0,     0,   174,   175,   -73,
     -53,   176,     0,   177,   178,   179,  -220,   -52,     0,     0,
      85,    86,    87,    88,    89,    90,    91,    92,   111,     0,
     112,   -67,   115,     0,   116,  -217,    84,     0,   113,   -66,
     114,   -73,     0,   -68,   -49,   117,    84,   118,  -220,   453,
    -218,   -50,     0,     0,    85,    86,    87,    88,    89,    90,
      91,    92,  -219,   -54,    85,    86,    87,    88,    89,    90,
      91,    92,   115,   -67,   116,   -68,     0,   -71,   -52,   119,
      84,   120,  -218,   -73,     0,     0,     0,   -69,   -50,     0,
       0,   121,     0,   122,  -219,     0,     0,     0,    85,    86,
      87,    88,    89,    90,    91,    92,   117,     0,   118,   -71,
       0,     0,     0,     0,    84,     0,   119,   -68,   120,   -69,
       0,   -75,   -54,     0,    84,     0,     0,     0,     0,     0,
       0,     0,    85,    86,    87,    88,    89,    90,    91,    92,
       0,     0,    85,    86,    87,    88,    89,    90,    91,    92,
     121,   -71,   122,   -75,     0,     0,     0,     0,    84,     0,
       0,   -69,    63,    64,    65,    66,    67,    68,    69,    70,
     397,     0,     0,     0,     0,     0,    85,    86,    87,    88,
      89,    90,    91,    92,     2,     0,     0,     0,     3,     4,
       0,     0,     0,     0,     0,   -75,     5,     6,     7,    84,
       0,     0,     8,     9,     0,     0,     0,    10,    11,    12,
     199,     0,     0,    15,     0,     0,     0,    85,    86,    87,
      88,    89,    90,    91,    92,   399,    72,     0,   200,   201,
     202,   203,   204,   205,   206,   207,   208,     0,     0,     2,
       0,     0,     0,     3,     4,     0,     0,     0,     0,     0,
     191,     5,     6,     7,    84,     0,     0,     8,     9,     0,
       0,   297,    10,    11,    12,    13,     0,     0,    15,     0,
       0,     0,    85,    86,    87,    88,    89,    90,    91,    92,
     401,    72,     0,    73,    74,    75,    76,    77,    78,    79,
      80,     0,     0,     0,     2,     0,     0,     0,     3,     4,
       0,     0,     0,     0,     0,   191,     5,     6,     7,    84,
       0,     0,     8,     9,     0,     0,     0,    10,    11,    12,
     199,     0,     0,    15,     0,     0,     0,    85,    86,    87,
      88,    89,    90,    91,    92,   403,    72,     0,   414,   415,
     416,   417,   418,   419,   420,   421,   422,     0,     0,     2,
       0,     0,     0,     3,     4,     0,     0,     0,     0,     0,
     191,     5,     6,     7,    84,     0,     0,     8,     9,     0,
       0,     0,    10,    11,    12,    13,   405,     0,    15,     0,
       0,     0,    85,    86,    87,    88,    89,    90,    91,    92,
     410,    72,     0,    73,    74,    75,    76,    77,    78,    79,
      80,     0,   501,     0,     0,    84,   278,   279,   280,   281,
     282,   283,   284,   285,     0,   191,   503,     0,     0,    84,
       0,     0,     0,    85,    86,    87,    88,    89,    90,    91,
      92,    84,   505,     0,     0,     0,     0,    85,    86,    87,
      88,    89,    90,    91,    92,    84,   507,     0,     0,    85,
      86,    87,    88,    89,    90,    91,    92,     0,   510,     0,
       0,    84,     0,    85,    86,    87,    88,    89,    90,    91,
      92,     0,     0,     0,     0,    84,     0,     0,     0,    85,
      86,    87,    88,    89,    90,    91,    92,    84,     0,     0,
       0,     0,     0,    85,    86,    87,    88,    89,    90,    91,
      92,     0,   -94,     0,     0,    85,    86,    87,    88,    89,
      90,    91,    92,     2,     0,     0,     0,     3,     4,     0,
     107,   -94,   108,     0,     0,     5,     6,     7,   -94,     0,
       0,     8,     9,     0,     0,     0,    10,    11,    12,    13,
       0,     0,    15,     0,     0,     0,   -94,   -94,   -94,   -94,
     -94,   -94,   -94,   -94,   -95,    72,     0,    73,    74,    75,
      76,    77,    78,    79,    80,     2,     0,     0,     0,     3,
       4,     0,   109,   -95,   110,     0,     0,     5,     6,     7,
     -95,     0,     0,     8,     9,     0,     0,     0,    10,    11,
      12,    13,     0,     0,    15,     0,     0,     0,   -95,   -95,
     -95,   -95,   -95,   -95,   -95,   -95,     0,    24,     0,    73,
      74,    75,    76,    77,    78,    79,    80,     2,     0,     0,
       0,     3,     4,     0,     0,     0,     0,     0,     0,     5,
       6,     7,     0,     0,     0,     8,     9,     0,     0,     0,
      10,    11,    12,   199,     0,     0,    15,     0,     0,     0,
       0,     2,     0,     0,     0,     0,     0,     0,     0,    72,
       0,   302,   303,   304,   305,   306,   307,   308,   309,     8,
       9,     0,     0,     0,    10,    11,    12,    13,     0,     0,
      15,     0,     0,     0,     0,   -96,     0,     0,     0,     0,
       0,     0,     0,    72,  -102,    73,    74,    75,    76,    77,
      78,    79,    80,   111,   -96,   112,     0,     0,     0,     0,
       0,   -96,   113,  -102,   114,     0,     0,     0,     0,     0,
    -102,     0,     0,     0,     0,     0,     0,     0,   -97,   -96,
     -96,   -96,   -96,   -96,   -96,   -96,   -96,  -101,  -102,  -102,
    -102,  -102,  -102,  -102,  -102,  -102,   115,   -97,   116,     0,
       0,     0,     0,     0,   -97,   117,  -101,   118,     0,     0,
       0,     0,     0,  -101,     0,     0,     0,     0,     0,     0,
       0,   -98,   -97,   -97,   -97,   -97,   -97,   -97,   -97,   -97,
    -103,  -101,  -101,  -101,  -101,  -101,  -101,  -101,  -101,   119,
     -98,   120,     0,     0,     0,     0,     0,   -98,   121,  -103,
     122,     0,     0,     0,     0,     0,  -103,     0,     0,     0,
       0,     0,     0,   323,     0,   -98,   -98,   -98,   -98,   -98,
     -98,   -98,   -98,   465,  -103,  -103,  -103,  -103,  -103,  -103,
    -103,  -103,   324,     0,     0,     0,     0,     0,     0,    84,
       0,     0,   466,     0,     0,     0,     0,     0,     0,    84,
       0,     0,     0,     0,     0,   104,     0,   325,    86,   326,
     327,   328,   329,   330,   331,   104,    84,   467,    86,   468,
      88,   469,    90,   470,    92,     0,    84,     0,     0,     0,
       0,     0,   104,     0,    85,    86,    87,    88,    89,    90,
      91,    92,   234,    84,    85,    86,    87,    88,    89,    90,
      91,    92,   286,   287,   288,   289,   290,   291,   292,   293,
       0,    85,    86,    87,    88,    89,    90,    91,    92,   149,
     150,   151,   152,   153,   154,   155,   156,   157,   158,   159,
     160,   161,   162,   163,   164,   165,   166,   167,   151,   152,
     153,   154,   155,   156,   157,   158,   159,   160,   161,   162,
     163,   164,   165,   166,   167,   153,   154,   155,   156,   157,
     158,   159,   160,   161,   162,   163,   164,   165,   166,   167,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   157,   158,   159,   160,   161,   162,   163,
     164,   165,   166,   167
};

static const yytype_int16 yycheck[] =
{
      14,     8,    41,   186,     9,   185,   186,   102,    25,    37,
      47,    25,   102,   324,    34,    22,     5,   107,     5,   109,
       5,   111,    36,   113,   335,   115,     6,   117,    79,   119,
      44,   121,    45,     5,   351,   109,     5,   111,     5,   113,
       5,   115,     5,   117,   361,   119,     5,   121,    45,     5,
      37,    79,     5,    38,    34,    62,     5,   237,   238,   149,
     150,   151,   152,   153,   154,   155,   156,   157,   158,   159,
     160,   161,   162,   163,   164,   165,   166,   167,   168,   169,
     391,     6,     0,   178,     5,    38,     5,   177,   178,    38,
      79,     5,   106,     5,   108,   134,   110,     5,   112,     5,
     114,   138,   116,   177,   118,     6,   120,    79,   122,    79,
      79,     5,    79,   127,    79,   129,    79,    38,   135,    38,
      79,     5,   136,    79,    38,     5,    38,     5,   135,    37,
     137,    37,   139,    34,     5,    34,    35,   142,   143,   144,
     145,   146,   147,   148,    38,    70,    71,    72,    73,    74,
      75,    46,     5,     5,    38,   466,   336,     6,    38,     6,
      38,    16,     5,     6,     5,   179,    37,   350,    13,     5,
     350,   185,   186,   187,    79,     5,     6,    34,   183,     5,
     197,     6,     5,   197,    37,    37,    35,    34,   499,   196,
     229,    34,     5,    36,    37,    80,    37,    42,   172,   173,
     383,    37,   382,   383,    34,    34,    36,    37,    58,   248,
      35,    37,    42,   227,    37,    60,    61,    62,    63,    64,
      65,    66,    67,   228,    37,   239,   316,   317,     5,     6,
      60,    61,    62,    63,    64,    65,    66,    67,    35,     5,
      58,     5,   337,   317,    31,    42,    33,   337,    34,    79,
      36,     5,     6,     5,     5,   435,   436,    34,    22,    36,
      37,   300,    58,    60,    61,    62,    63,    64,    65,    66,
      67,    37,     5,     6,     5,     5,   371,     5,   368,   369,
      34,   371,    36,    37,    37,    37,    37,     5,   302,   303,
     304,   305,   306,   307,   308,   309,    34,   477,    36,     5,
       6,    34,   482,    36,    37,    17,    37,   324,    34,    37,
     324,   325,   326,   327,   328,   329,   330,   331,   335,    37,
     334,   335,     5,   337,    16,     5,    17,   366,    34,    16,
      36,    37,     6,     5,     6,     6,   350,   351,    34,    35,
     354,   355,    34,    35,   349,   359,     5,   361,    35,    34,
     364,    36,    38,    82,    37,   362,   363,    37,    34,   539,
      36,     5,    34,   453,    36,    37,     5,     6,    79,   383,
       5,     6,     6,    17,   391,    35,    58,   391,    22,   393,
      34,   395,    36,   397,    34,   399,    36,   401,    35,   403,
     407,   405,    34,   407,    36,    34,   410,    36,    37,    34,
       6,    36,    37,    35,   409,    22,    23,    24,    25,    26,
     449,     6,    13,     5,     6,    34,   430,    36,   425,    56,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,     5,
       6,    42,    34,     4,    36,    37,     5,     6,    38,   466,
      58,     5,   466,   467,   468,   469,   470,    38,   472,    60,
      61,    62,    63,    64,    65,    66,    67,    58,    34,   483,
      36,    37,     6,     5,   488,    34,    79,    36,    37,   486,
     487,    34,   499,     5,   523,   499,     5,   501,    34,   503,
      38,   505,     3,   507,    82,    83,   510,    35,    38,    35,
     537,     6,   541,   542,   509,    16,     6,   512,   515,    20,
      21,   550,     5,    38,    58,    36,    58,    28,    29,    30,
      17,    34,   561,    34,    35,    80,    38,    37,    39,    40,
      41,    42,    43,    36,    45,    17,    47,    48,   552,    50,
      51,    52,    53,    54,    55,    36,     3,    58,    59,    60,
      61,    62,    63,    64,    65,    66,    67,    68,    69,    16,
      38,    79,    17,    20,    21,    37,    77,    37,    79,    80,
      81,    28,    29,    30,    37,   239,   430,    34,    35,   483,
     187,     2,    39,    40,    41,    42,    43,   228,    45,    53,
      47,    48,   365,    50,    51,    52,    53,    54,    55,    -1,
       3,    58,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    69,    16,    -1,    -1,    -1,    20,    21,    -1,
      77,    -1,    79,    80,    81,    28,    29,    30,    -1,    -1,
      -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      43,    -1,    45,    -1,    47,    48,    -1,    50,    51,    52,
      53,    54,    55,    -1,     3,    58,    59,    60,    61,    62,
      63,    64,    65,    66,    67,    68,    69,    16,    -1,    -1,
      -1,    20,    21,    -1,    77,    -1,    79,    80,    81,    28,
      29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    43,    -1,    45,    -1,    47,    48,
      -1,    50,    51,    52,    53,    54,    55,    -1,    -1,    58,
      59,    60,    61,    62,    63,    64,    65,    66,    67,    68,
      69,    20,    21,    22,    23,    24,    25,    26,    77,    -1,
      79,    80,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    35,
      -1,    -1,    28,    29,    -1,    -1,    42,    -1,    34,    35,
      36,    35,    38,    -1,    -1,    -1,    82,    -1,    42,     5,
       6,    -1,    -1,    -1,    60,    61,    62,    63,    64,    65,
      66,    67,    -1,    -1,    -1,    82,    60,    61,    62,    63,
      64,    65,    66,    67,    -1,     5,     6,    -1,    34,    -1,
      36,    37,    -1,    -1,    82,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    34,    -1,    36,    37,    -1,     5,
       6,    -1,    42,     5,     6,    -1,    -1,    28,    29,     5,
       6,    32,    -1,    34,    35,    36,     5,     6,    -1,    -1,
      60,    61,    62,    63,    64,    65,    66,    67,    34,    -1,
      36,    37,    34,    -1,    36,    37,    42,    -1,    34,    79,
      36,    37,    -1,     5,     6,    34,    42,    36,    37,    82,
       5,     6,    -1,    -1,    60,    61,    62,    63,    64,    65,
      66,    67,     5,     6,    60,    61,    62,    63,    64,    65,
      66,    67,    34,    79,    36,    37,    -1,     5,     6,    34,
      42,    36,    37,    79,    -1,    -1,    -1,     5,     6,    -1,
      -1,    34,    -1,    36,    37,    -1,    -1,    -1,    60,    61,
      62,    63,    64,    65,    66,    67,    34,    -1,    36,    37,
      -1,    -1,    -1,    -1,    42,    -1,    34,    79,    36,    37,
      -1,     5,     6,    -1,    42,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    60,    61,    62,    63,    64,    65,    66,    67,
      -1,    -1,    60,    61,    62,    63,    64,    65,    66,    67,
      34,    79,    36,    37,    -1,    -1,    -1,    -1,    42,    -1,
      -1,    79,    60,    61,    62,    63,    64,    65,    66,    67,
      13,    -1,    -1,    -1,    -1,    -1,    60,    61,    62,    63,
      64,    65,    66,    67,    16,    -1,    -1,    -1,    20,    21,
      -1,    -1,    -1,    -1,    -1,    79,    28,    29,    30,    42,
      -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
      42,    -1,    -1,    45,    -1,    -1,    -1,    60,    61,    62,
      63,    64,    65,    66,    67,    13,    58,    -1,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    -1,    -1,    16,
      -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    -1,    -1,
      82,    28,    29,    30,    42,    -1,    -1,    34,    35,    -1,
      -1,    38,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      -1,    -1,    60,    61,    62,    63,    64,    65,    66,    67,
      13,    58,    -1,    60,    61,    62,    63,    64,    65,    66,
      67,    -1,    -1,    -1,    16,    -1,    -1,    -1,    20,    21,
      -1,    -1,    -1,    -1,    -1,    82,    28,    29,    30,    42,
      -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
      42,    -1,    -1,    45,    -1,    -1,    -1,    60,    61,    62,
      63,    64,    65,    66,    67,    13,    58,    -1,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    -1,    -1,    16,
      -1,    -1,    -1,    20,    21,    -1,    -1,    -1,    -1,    -1,
      82,    28,    29,    30,    42,    -1,    -1,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    13,    -1,    45,    -1,
      -1,    -1,    60,    61,    62,    63,    64,    65,    66,    67,
      13,    58,    -1,    60,    61,    62,    63,    64,    65,    66,
      67,    -1,    13,    -1,    -1,    42,    60,    61,    62,    63,
      64,    65,    66,    67,    -1,    82,    13,    -1,    -1,    42,
      -1,    -1,    -1,    60,    61,    62,    63,    64,    65,    66,
      67,    42,    13,    -1,    -1,    -1,    -1,    60,    61,    62,
      63,    64,    65,    66,    67,    42,    13,    -1,    -1,    60,
      61,    62,    63,    64,    65,    66,    67,    -1,    13,    -1,
      -1,    42,    -1,    60,    61,    62,    63,    64,    65,    66,
      67,    -1,    -1,    -1,    -1,    42,    -1,    -1,    -1,    60,
      61,    62,    63,    64,    65,    66,    67,    42,    -1,    -1,
      -1,    -1,    -1,    60,    61,    62,    63,    64,    65,    66,
      67,    -1,    16,    -1,    -1,    60,    61,    62,    63,    64,
      65,    66,    67,    16,    -1,    -1,    -1,    20,    21,    -1,
      34,    35,    36,    -1,    -1,    28,    29,    30,    42,    -1,
      -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      -1,    -1,    45,    -1,    -1,    -1,    60,    61,    62,    63,
      64,    65,    66,    67,    16,    58,    -1,    60,    61,    62,
      63,    64,    65,    66,    67,    16,    -1,    -1,    -1,    20,
      21,    -1,    34,    35,    36,    -1,    -1,    28,    29,    30,
      42,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,
      41,    42,    -1,    -1,    45,    -1,    -1,    -1,    60,    61,
      62,    63,    64,    65,    66,    67,    -1,    58,    -1,    60,
      61,    62,    63,    64,    65,    66,    67,    16,    -1,    -1,
      -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    -1,    28,
      29,    30,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,
      -1,    16,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    58,
      -1,    60,    61,    62,    63,    64,    65,    66,    67,    34,
      35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
      45,    -1,    -1,    -1,    -1,    16,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    58,    16,    60,    61,    62,    63,    64,
      65,    66,    67,    34,    35,    36,    -1,    -1,    -1,    -1,
      -1,    42,    34,    35,    36,    -1,    -1,    -1,    -1,    -1,
      42,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    16,    60,
      61,    62,    63,    64,    65,    66,    67,    16,    60,    61,
      62,    63,    64,    65,    66,    67,    34,    35,    36,    -1,
      -1,    -1,    -1,    -1,    42,    34,    35,    36,    -1,    -1,
      -1,    -1,    -1,    42,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    16,    60,    61,    62,    63,    64,    65,    66,    67,
      16,    60,    61,    62,    63,    64,    65,    66,    67,    34,
      35,    36,    -1,    -1,    -1,    -1,    -1,    42,    34,    35,
      36,    -1,    -1,    -1,    -1,    -1,    42,    -1,    -1,    -1,
      -1,    -1,    -1,    16,    -1,    60,    61,    62,    63,    64,
      65,    66,    67,    16,    60,    61,    62,    63,    64,    65,
      66,    67,    35,    -1,    -1,    -1,    -1,    -1,    -1,    42,
      -1,    -1,    35,    -1,    -1,    -1,    -1,    -1,    -1,    42,
      -1,    -1,    -1,    -1,    -1,    58,    -1,    60,    61,    62,
      63,    64,    65,    66,    67,    58,    42,    60,    61,    62,
      63,    64,    65,    66,    67,    -1,    42,    -1,    -1,    -1,
      -1,    -1,    58,    -1,    60,    61,    62,    63,    64,    65,
      66,    67,    58,    42,    60,    61,    62,    63,    64,    65,
      66,    67,    60,    61,    62,    63,    64,    65,    66,    67,
      -1,    60,    61,    62,    63,    64,    65,    66,    67,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    16,    20,    21,    28,    29,    30,    34,    35,
      39,    40,    41,    42,    43,    45,    47,    48,    50,    51,
      52,    53,    54,    55,    58,    59,    60,    61,    62,    63,
      64,    65,    66,    67,    68,    69,    77,    79,    80,    85,
      86,    87,    93,    96,    97,   100,   103,   107,   109,   110,
     113,   115,   120,   121,   122,   123,   124,   129,   130,   131,
     132,   133,    34,    60,    61,    62,    63,    64,    65,    66,
      67,   129,    58,    60,    61,    62,    63,    64,    65,    66,
      67,   120,   122,   128,    42,    60,    61,    62,    63,    64,
      65,    66,    67,    90,   101,   102,    45,    45,    79,    79,
     120,    16,    34,    35,    58,    90,    93,    34,    36,    34,
      36,    34,    36,    34,    36,    34,    36,    34,    36,    34,
      36,    34,    36,    90,     0,    46,   115,   104,    16,    35,
      90,    94,    99,    79,    87,    34,    35,    34,    80,     5,
      79,   131,     6,    70,    71,    72,    73,    74,    75,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    82,    83,
       5,    22,    31,    33,    28,    29,    32,    34,    35,    36,
     120,    17,    37,     5,    38,     6,    34,     5,    79,    79,
      58,    82,   124,   125,   126,    58,     6,    34,    90,    42,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    90,
     125,   127,    90,   127,    90,   127,    90,   127,    90,   127,
      90,   127,    90,   127,    90,   127,    90,    34,    35,   118,
     119,    90,    91,    58,    58,    89,    90,     6,    35,     5,
      81,   108,    93,   111,   120,    90,   114,   120,    87,   120,
     122,   122,   122,   122,   122,   122,   122,   124,   124,   124,
     124,   124,   124,   124,   124,   124,   124,   124,   124,   124,
     124,   124,   124,   124,   124,   124,   124,   124,    60,    61,
      62,    63,    64,    65,    66,    67,    60,    61,    62,    63,
      64,    65,    66,    67,   130,   130,   127,    38,   125,    90,
      37,   122,    60,    61,    62,    63,    64,    65,    66,    67,
      90,    92,   123,    92,   101,    17,    82,     5,    37,    37,
       5,    38,   120,    16,    35,    60,    62,    63,    64,    65,
      66,    67,    88,    90,    93,    34,     6,     5,    37,    37,
      37,    37,    37,    37,    37,    37,    90,   128,   115,     6,
      34,     5,    79,    17,    38,     5,    38,   123,   123,    35,
      94,   112,    79,    79,     5,    82,    37,    81,    82,    82,
      37,     5,    38,   115,    90,    90,    90,    90,    90,    90,
      90,    90,     6,     5,    37,    35,   124,   127,    58,    35,
      79,    35,    88,    13,    90,    13,    90,    13,    90,    13,
      90,    13,    90,    13,    90,    13,    90,     5,    37,     6,
      13,    90,    88,   123,    60,    61,    62,    63,    64,    65,
      66,    67,    68,    90,   125,     6,    38,   122,    92,    91,
      35,    90,    95,    98,    90,     6,    35,    38,    91,   120,
     120,    90,   133,   115,    56,   117,   124,   124,   125,     4,
     123,    92,    58,    82,    38,    58,    88,    38,    90,    90,
      90,    90,    90,    90,    90,    16,    35,    60,    62,    64,
      66,    90,    93,   106,   122,    90,    37,     6,   120,    79,
      37,    89,    34,     5,   123,   123,    79,    79,     5,    38,
      34,    38,   115,     5,    38,   124,    35,    38,    38,    35,
      88,    13,    90,    13,    90,    13,    90,    13,    90,     6,
      13,    90,     6,   105,   123,     5,    38,   123,    95,    38,
     120,   120,    90,   116,    36,    58,    58,    17,    88,    38,
      90,    90,    90,    90,   122,    90,   122,    80,   120,    34,
      37,    37,    37,   115,    36,    38,    38,    38,    17,    79,
      87,    37,    79,   123,   115,   115,    36,    17,    81,    90,
      37,    37,    37,   115
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    84,    85,    86,    87,    87,    88,    88,    88,    88,
      88,    88,    88,    88,    88,    88,    88,    88,    88,    88,
      88,    88,    88,    88,    88,    88,    88,    88,    88,    88,
      88,    88,    88,    88,    88,    88,    88,    88,    88,    88,
      88,    88,    88,    88,    89,    89,    90,    90,    90,    90,
      90,    90,    90,    90,    90,    91,    91,    91,    91,    92,
      92,    92,    92,    92,    92,    92,    92,    92,    92,    92,
      92,    92,    92,    92,    92,    92,    92,    92,    93,    93,
      93,    93,    93,    93,    93,    93,    94,    94,    94,    94,
      94,    94,    95,    95,    96,    96,    96,    96,    96,    96,
      96,    96,    96,    96,    97,    97,    98,    98,    99,    99,
     100,   100,   101,   101,   102,   102,   104,   103,   103,   103,
     103,   105,   103,   106,   103,   107,   108,   109,   110,   112,
     111,   113,   114,   114,   114,   115,   115,   115,   115,   115,
     115,   116,   115,   115,   115,   115,   115,   115,   115,   115,
     115,   115,   115,   115,   117,   118,   118,   119,   120,   120,
     121,   121,   121,   121,   121,   122,   122,   122,   122,   122,
     122,   122,   122,   123,   123,   123,   123,   124,   124,   124,
     124,   124,   124,   124,   124,   124,   124,   124,   124,   124,
     124,   124,   124,   124,   124,   124,   124,   125,   125,   125,
     125,   126,   127,   127,   127,   127,   127,   127,   127,   127,
     127,   127,   127,   127,   127,   127,   127,   127,   127,   127,
     127,   127,   127,   127,   127,   128,   128,   129,   129,   129,
     129,   129,   129,   129,   129,   129,   129,   129,   129,   129,
     129,   129,   129,   129,   129,   129,   129,   129,   129,   129,
     129,   130,   130,   131,   131,   131,   132,   132,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133,   133,   133,   133,   133,   133,   133,   133,   133,
     133,   133
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     2,     1,     1,     2,     0,     1,     3,     2,
       3,     2,     3,     2,     3,     2,     3,     2,     3,     2,
       3,     2,     3,     2,     3,     3,     5,     4,     3,     5,
       6,     7,     6,     5,     4,     5,     4,     5,     4,     5,
       4,     5,     4,     5,     1,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     4,     3,     1,
       2,     2,     2,     2,     2,     1,     1,     1,     1,     1,
       2,     1,     2,     1,     2,     1,     3,     3,     1,     4,
       7,     6,     9,     4,     7,     9,     1,     4,     3,     3,
       6,     5,     4,     6,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     4,     1,     3,     1,     3,
       2,     5,     4,     3,     1,     3,     0,     4,     3,     2,
       5,     0,    10,     0,     9,     1,     1,     1,     1,     0,
       3,     1,     1,     3,     5,     1,     2,     2,     5,     2,
       1,     0,     8,     9,     9,     5,     5,     7,     3,     3,
       6,     2,     2,     3,     7,     7,     9,     2,     1,     3,
       1,     1,     1,     1,     1,     1,     3,     3,     3,     3,
       3,     3,     3,     1,     5,     3,     5,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     1,     3,
       5,     3,     0,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     3,     1,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     5,     1,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     1,     2,     1,     3,     3,     1,     2,     1,     1,
       1,     1,     1,     4,     4,     6,     3,     3,     3,     4,
       3,     4,     3,     4,     3,     4,     3,     4,     3,     4,
       3,     4,     3,     4,     2,     2,     4,     4,     3,     3,
       3,     1
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
#line 342 "lg.ypp"
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
#line 2282 "lg.tab.cpp"
    break;

  case 4:
#line 416 "lg.ypp"
                                                                          {(yyval.cinst) = (yyvsp[0].cexp);}
#line 2288 "lg.tab.cpp"
    break;

  case 5:
#line 417 "lg.ypp"
                                                                          {(yyval.cinst) = ((yyvsp[-1].cinst)+=(yyvsp[0].cexp));}
#line 2294 "lg.tab.cpp"
    break;

  case 6:
#line 423 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId();}
#line 2300 "lg.tab.cpp"
    break;

  case 7:
#line 424 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
#line 2306 "lg.tab.cpp"
    break;

  case 8:
#line 425 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str),(yyvsp[0].cexp)));}
#line 2312 "lg.tab.cpp"
    break;

  case 9:
#line 426 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,2> **>()));}
#line 2318 "lg.tab.cpp"
    break;

  case 10:
#line 427 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,2> **>(),true));}
#line 2324 "lg.tab.cpp"
    break;

  case 11:
#line 428 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,3> **>()));}
#line 2330 "lg.tab.cpp"
    break;

  case 12:
#line 429 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,3> **>(),true));}
#line 2336 "lg.tab.cpp"
    break;

  case 13:
#line 430 "lg.ypp"
                                                { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,30> **>()));}
#line 2342 "lg.tab.cpp"
    break;

  case 14:
#line 431 "lg.ypp"
                                                { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,30> **>(),true));}
#line 2348 "lg.tab.cpp"
    break;

  case 15:
#line 432 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,4> **>()));}
#line 2354 "lg.tab.cpp"
    break;

  case 16:
#line 433 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,4> **>(),true));}
#line 2360 "lg.tab.cpp"
    break;

  case 17:
#line 434 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,40> **>())); }
#line 2366 "lg.tab.cpp"
    break;

  case 18:
#line 435 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,40> **>(),true)); }
#line 2372 "lg.tab.cpp"
    break;

  case 19:
#line 436 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,5> **>()));}
#line 2378 "lg.tab.cpp"
    break;

  case 20:
#line 437 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,5> **>(),true));}
#line 2384 "lg.tab.cpp"
    break;

  case 21:
#line 438 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,50> **>()));}
#line 2390 "lg.tab.cpp"
    break;

  case 22:
#line 439 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,50> **>(),true));}
#line 2396 "lg.tab.cpp"
    break;

  case 23:
#line 440 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-1].type)->right()));}
#line 2402 "lg.tab.cpp"
    break;

  case 24:
#line 441 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-2].type),true));}
#line 2408 "lg.tab.cpp"
    break;

  case 25:
#line 442 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id)));}
#line 2414 "lg.tab.cpp"
    break;

  case 26:
#line 443 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-2].clist_id),true,true));}
#line 2420 "lg.tab.cpp"
    break;

  case 27:
#line 444 "lg.ypp"
                                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id),true,false));}
#line 2426 "lg.tab.cpp"
    break;

  case 28:
#line 445 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-2].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
#line 2432 "lg.tab.cpp"
    break;

  case 29:
#line 446 "lg.ypp"
                                                { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id)));}
#line 2438 "lg.tab.cpp"
    break;

  case 30:
#line 447 "lg.ypp"
                                                    { (yyval.clist_id) = (yyvsp[-5].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-2].clist_id),false,true));}
#line 2444 "lg.tab.cpp"
    break;

  case 31:
#line 448 "lg.ypp"
                                                        { (yyval.clist_id) = (yyvsp[-6].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-2].clist_id),true,true));}
#line 2450 "lg.tab.cpp"
    break;

  case 32:
#line 449 "lg.ypp"
                                                     { (yyval.clist_id) = (yyvsp[-5].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id),true,false));}
#line 2456 "lg.tab.cpp"
    break;

  case 33:
#line 450 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str),(yyvsp[0].cexp)));}
#line 2462 "lg.tab.cpp"
    break;

  case 34:
#line 451 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,2> **>()));}
#line 2468 "lg.tab.cpp"
    break;

  case 35:
#line 452 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,2> **>(),true));}
#line 2474 "lg.tab.cpp"
    break;

  case 36:
#line 453 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,3> **>()));}
#line 2480 "lg.tab.cpp"
    break;

  case 37:
#line 454 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,3> **>(),true));}
#line 2486 "lg.tab.cpp"
    break;

  case 38:
#line 455 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,4> **>()));}
#line 2492 "lg.tab.cpp"
    break;

  case 39:
#line 456 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,4> **>(),true));}
#line 2498 "lg.tab.cpp"
    break;

  case 40:
#line 457 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double,5> **>()));}
#line 2504 "lg.tab.cpp"
    break;

  case 41:
#line 458 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double,5> **>(),true));}
#line 2510 "lg.tab.cpp"
    break;

  case 42:
#line 459 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-1].type)->right()));}
#line 2516 "lg.tab.cpp"
    break;

  case 43:
#line 460 "lg.ypp"
                                               { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-2].type),true));}
#line 2522 "lg.tab.cpp"
    break;

  case 44:
#line 463 "lg.ypp"
                                      { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
#line 2528 "lg.tab.cpp"
    break;

  case 45:
#line 464 "lg.ypp"
                                      { (yyval.clist_id)=(yyvsp[-2].clist_id)  ; (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
#line 2534 "lg.tab.cpp"
    break;

  case 55:
#line 472 "lg.ypp"
                                    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[0].str),dcltype);}
#line 2540 "lg.tab.cpp"
    break;

  case 56:
#line 473 "lg.ypp"
                                    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[-2].str),dcltype,(yyvsp[0].cexp));}
#line 2546 "lg.tab.cpp"
    break;

  case 57:
#line 474 "lg.ypp"
                                    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[-3].str),dcltype,(yyvsp[-1].args));(yyvsp[-1].args).destroy();}
#line 2552 "lg.tab.cpp"
    break;

  case 58:
#line 475 "lg.ypp"
                                    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2558 "lg.tab.cpp"
    break;

  case 59:
#line 481 "lg.ypp"
                                                 { (yyval.args)=(yyvsp[0].cexp);}
#line 2564 "lg.tab.cpp"
    break;

  case 60:
#line 482 "lg.ypp"
                                                 { (yyval.args)=Find((yyvsp[-1].str));}
#line 2570 "lg.tab.cpp"
    break;

  case 61:
#line 483 "lg.ypp"
                                                 { (yyval.args)=Find((yyvsp[-1].str));}
#line 2576 "lg.tab.cpp"
    break;

  case 62:
#line 484 "lg.ypp"
                                                 { (yyval.args)=Find((yyvsp[-1].str));}
#line 2582 "lg.tab.cpp"
    break;

  case 63:
#line 485 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[-1].str));}
#line 2588 "lg.tab.cpp"
    break;

  case 64:
#line 486 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[-1].str));}
#line 2594 "lg.tab.cpp"
    break;

  case 65:
#line 487 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2600 "lg.tab.cpp"
    break;

  case 66:
#line 488 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2606 "lg.tab.cpp"
    break;

  case 67:
#line 489 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2612 "lg.tab.cpp"
    break;

  case 68:
#line 490 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2618 "lg.tab.cpp"
    break;

  case 69:
#line 491 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2624 "lg.tab.cpp"
    break;

  case 70:
#line 492 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[-1].str));}
#line 2630 "lg.tab.cpp"
    break;

  case 71:
#line 493 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2636 "lg.tab.cpp"
    break;

  case 72:
#line 494 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[-1].str));}
#line 2642 "lg.tab.cpp"
    break;

  case 73:
#line 495 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2648 "lg.tab.cpp"
    break;

  case 74:
#line 496 "lg.ypp"
                                                 { (yyval.args)=Find((yyvsp[-1].str));}
#line 2654 "lg.tab.cpp"
    break;

  case 75:
#line 497 "lg.ypp"
                                           { (yyval.args)=Find((yyvsp[0].str));}
#line 2660 "lg.tab.cpp"
    break;

  case 76:
#line 498 "lg.ypp"
                                                 { (yyval.args)=make_pair<const char *,const C_F0>((const char *) (yyvsp[-2].str),(C_F0) (yyvsp[0].cexp));}
#line 2666 "lg.tab.cpp"
    break;

  case 77:
#line 499 "lg.ypp"
                                                 { (yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].args));}
#line 2672 "lg.tab.cpp"
    break;

  case 79:
#line 505 "lg.ypp"
                                              {(yyval.type)=TypeArray((yyvsp[-3].type),(yyvsp[-1].type));}
#line 2678 "lg.tab.cpp"
    break;

  case 80:
#line 506 "lg.ypp"
                                              {(yyval.type)=TypeArray(TypeArray((yyvsp[-6].type),(yyvsp[-4].type)),(yyvsp[-1].type));}
#line 2684 "lg.tab.cpp"
    break;

  case 81:
#line 507 "lg.ypp"
                                              {(yyval.type)=TypeArray((yyvsp[-5].type),(yyvsp[-3].type),(yyvsp[-1].type));}
#line 2690 "lg.tab.cpp"
    break;

  case 82:
#line 508 "lg.ypp"
                                              {(yyval.type)=TypeArray(TypeArray((yyvsp[-8].type),(yyvsp[-6].type),(yyvsp[-4].type)),(yyvsp[-1].type));}
#line 2696 "lg.tab.cpp"
    break;

  case 83:
#line 509 "lg.ypp"
                                              {(yyval.type)=TypeTemplate((yyvsp[-3].type),(yyvsp[-1].type));}
#line 2702 "lg.tab.cpp"
    break;

  case 84:
#line 510 "lg.ypp"
                                              {(yyval.type)=TypeArray(TypeTemplate((yyvsp[-6].type),(yyvsp[-4].type)),(yyvsp[-1].type));}
#line 2708 "lg.tab.cpp"
    break;

  case 85:
#line 511 "lg.ypp"
                                              {(yyval.type)=TypeArray(TypeTemplate((yyvsp[-8].type),(yyvsp[-6].type)),(yyvsp[-3].type),(yyvsp[-1].type));}
#line 2714 "lg.tab.cpp"
    break;

  case 86:
#line 515 "lg.ypp"
                                              { (yyval.cexp) =  NewFEvariable((yyvsp[0].str),currentblock,fespacetype,fespacecomplex,fespacedim); }
#line 2720 "lg.tab.cpp"
    break;

  case 87:
#line 516 "lg.ypp"
                                              { (yyval.cexp) =  NewFEarray((yyvsp[-3].str),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim); }
#line 2726 "lg.tab.cpp"
    break;

  case 88:
#line 517 "lg.ypp"
                                              { (yyval.cexp) =  NewFEvariable((yyvsp[-2].str),currentblock,fespacetype,(yyvsp[0].cexp),fespacecomplex,fespacedim);}
#line 2732 "lg.tab.cpp"
    break;

  case 89:
#line 518 "lg.ypp"
                                              { (yyval.cexp) =  NewFEvariable((yyvsp[-1].clist_id),currentblock,fespacetype,fespacecomplex,fespacedim);}
#line 2738 "lg.tab.cpp"
    break;

  case 90:
#line 519 "lg.ypp"
                                              { (yyval.cexp) =  NewFEarray((yyvsp[-4].clist_id),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim);}
#line 2744 "lg.tab.cpp"
    break;

  case 91:
#line 520 "lg.ypp"
                                              { (yyval.cexp) =  NewFEvariable((yyvsp[-3].clist_id),currentblock,fespacetype,(yyvsp[0].cexp),fespacecomplex,fespacedim);}
#line 2750 "lg.tab.cpp"
    break;

  case 92:
#line 524 "lg.ypp"
                                     { (yyval.cexp) =  NewFEarray((yyvsp[-3].str),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim); }
#line 2756 "lg.tab.cpp"
    break;

  case 93:
#line 525 "lg.ypp"
                                              { (yyval.cexp) =  NewFEarray((yyvsp[-4].clist_id),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex,fespacedim);}
#line 2762 "lg.tab.cpp"
    break;

  case 94:
#line 530 "lg.ypp"
              { fespacedim=2;}
#line 2768 "lg.tab.cpp"
    break;

  case 95:
#line 531 "lg.ypp"
               { fespacedim=1;}
#line 2774 "lg.tab.cpp"
    break;

  case 96:
#line 532 "lg.ypp"
               { fespacedim=3;}
#line 2780 "lg.tab.cpp"
    break;

  case 97:
#line 533 "lg.ypp"
               { fespacedim=4;}
#line 2786 "lg.tab.cpp"
    break;

  case 98:
#line 534 "lg.ypp"
               { fespacedim=5;}
#line 2792 "lg.tab.cpp"
    break;

  case 99:
#line 535 "lg.ypp"
                { fespacedim=6;}
#line 2798 "lg.tab.cpp"
    break;

  case 100:
#line 536 "lg.ypp"
               { fespacedim=7;}
#line 2804 "lg.tab.cpp"
    break;

  case 101:
#line 537 "lg.ypp"
                { fespacedim = 40;}
#line 2810 "lg.tab.cpp"
    break;

  case 102:
#line 538 "lg.ypp"
                { fespacedim = 30;}
#line 2816 "lg.tab.cpp"
    break;

  case 103:
#line 539 "lg.ypp"
                { fespacedim = 50;}
#line 2822 "lg.tab.cpp"
    break;

  case 104:
#line 542 "lg.ypp"
                     {fespacecomplex=false;  fespacetype = Find((yyvsp[0].str));}
#line 2828 "lg.tab.cpp"
    break;

  case 105:
#line 543 "lg.ypp"
                                  {
             if ((yyvsp[-1].type) != typevarreal && (yyvsp[-1].type) != typevarcomplex) lgerror (" type of finite element <real> or <complex>");
             fespacecomplex=((yyvsp[-1].type)==typevarcomplex);
             fespacetype = Find((yyvsp[-3].str));}
#line 2837 "lg.tab.cpp"
    break;

  case 106:
#line 548 "lg.ypp"
                                {  (yyval.cexp) = (yyvsp[0].cexp);}
#line 2843 "lg.tab.cpp"
    break;

  case 107:
#line 549 "lg.ypp"
                                             { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2849 "lg.tab.cpp"
    break;

  case 108:
#line 551 "lg.ypp"
                          {  (yyval.cexp) = (yyvsp[0].cexp);}
#line 2855 "lg.tab.cpp"
    break;

  case 109:
#line 552 "lg.ypp"
                                       { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2861 "lg.tab.cpp"
    break;

  case 110:
#line 554 "lg.ypp"
                                                { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[0].cexp);}
#line 2867 "lg.tab.cpp"
    break;

  case 111:
#line 555 "lg.ypp"
                                                { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[0].cexp);}
#line 2873 "lg.tab.cpp"
    break;

  case 112:
#line 558 "lg.ypp"
                               {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,KN<size_t>>((yyvsp[-3].str),typeFESpace((yyvsp[-1].args)),(yyvsp[-1].args),dimFESpaceImage((yyvsp[-1].args)));
     (yyvsp[-1].args).destroy(); }
#line 2880 "lg.tab.cpp"
    break;

  case 113:
#line 561 "lg.ypp"
                            {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,KN<size_t>>((yyvsp[-2].str),typeFESpace((yyvsp[0].args)),(yyvsp[0].args),dimFESpaceImage((yyvsp[0].args)));
     (yyvsp[0].args).destroy(); }
#line 2887 "lg.tab.cpp"
    break;

  case 115:
#line 566 "lg.ypp"
                                                    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 2893 "lg.tab.cpp"
    break;

  case 116:
#line 569 "lg.ypp"
                           {dcltype=(yyvsp[0].type);}
#line 2899 "lg.tab.cpp"
    break;

  case 117:
#line 569 "lg.ypp"
                                                          {(yyval.cexp)=(yyvsp[-1].cexp);}
#line 2905 "lg.tab.cpp"
    break;

  case 118:
#line 570 "lg.ypp"
                                                  {(yyval.cexp)=(yyvsp[-1].cexp);}
#line 2911 "lg.tab.cpp"
    break;

  case 119:
#line 571 "lg.ypp"
                           { (yyval.cexp)=(yyvsp[-1].cexp);}
#line 2917 "lg.tab.cpp"
    break;

  case 120:
#line 572 "lg.ypp"
                                        {(yyval.cexp)=currentblock->NewID((yyvsp[-4].type),(yyvsp[-3].str),(yyvsp[-1].cexp));}
#line 2923 "lg.tab.cpp"
    break;

  case 121:
#line 574 "lg.ypp"
                   {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = (yyvsp[-4].type)->right();
                      routineinblock[kkembtype] = currentblock;
                      (yyvsp[-1].routine)=new Routine((yyvsp[-5].type),(yyvsp[-4].type)->right(),(yyvsp[-3].str),(yyvsp[-1].clist_id),currentblock);
		      // routineinblock[kkembtype]->Add($3,"(",$<routine>5); //pas recursif pour l'instanat test  FH 27 dec 2008
                     // cout << " \n after new routine \n " << endl;
                      }
#line 2936 "lg.tab.cpp"
    break;

  case 122:
#line 583 "lg.ypp"
                     { currentblock=(yyvsp[-5].routine)->Set((yyvsp[-1].cinst));
                       currentblock->Add((yyvsp[-7].str),"(",(yyvsp[-5].routine)); //pas recursif pour l'instant test  FH 27 dec 2008
                       kkembtype--;
                       (yyval.cexp)=0;

                        }
#line 2947 "lg.tab.cpp"
    break;

  case 123:
#line 590 "lg.ypp"
                      {Block::open(currentblock); (yyvsp[-4].type)->SetArgs((yyvsp[-1].clist_id));}
#line 2953 "lg.tab.cpp"
    break;

  case 124:
#line 592 "lg.ypp"
                      {  //$<cinst>$=currentblock->close(currentblock);
                         (yyval.cinst).setclose(Block::snewclose(currentblock));// Sep 2016 FH.
                         (yyval.cexp)=currentblock->NewID((yyvsp[-8].type),(yyvsp[-7].str),(yyvsp[-1].cexp),*(yyvsp[-5].clist_id));
                         delete (yyvsp[-5].clist_id); //  FH 23032005
                         }
#line 2963 "lg.tab.cpp"
    break;

  case 125:
#line 599 "lg.ypp"
            {  Block::open(currentblock);}
#line 2969 "lg.tab.cpp"
    break;

  case 126:
#line 600 "lg.ypp"
            { (yyval.endb)=Block::snewclose(currentblock);
//  $$=currentblock->close(currentblock);
}
#line 2977 "lg.tab.cpp"
    break;

  case 127:
#line 604 "lg.ypp"
               {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;}
#line 2984 "lg.tab.cpp"
    break;

  case 128:
#line 606 "lg.ypp"
                   {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;}
#line 2991 "lg.tab.cpp"
    break;

  case 129:
#line 610 "lg.ypp"
                {dcltype=(yyvsp[0].type); Block::open(currentblock);  }
#line 2997 "lg.tab.cpp"
    break;

  case 130:
#line 611 "lg.ypp"
                            {(yyval.cexp)=(yyvsp[0].cexp);}
#line 3003 "lg.tab.cpp"
    break;

  case 131:
#line 613 "lg.ypp"
         { Block::open(currentblock);}
#line 3009 "lg.tab.cpp"
    break;

  case 132:
#line 615 "lg.ypp"
               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));Block::open(currentblock); }
#line 3015 "lg.tab.cpp"
    break;

  case 133:
#line 616 "lg.ypp"
                       { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str)));(yyval.clist_id)->push_back(UnId((yyvsp[0].str)));Block::open(currentblock); }
#line 3021 "lg.tab.cpp"
    break;

  case 134:
#line 617 "lg.ypp"
                               { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-4].str)));(yyval.clist_id)->push_back(UnId((yyvsp[-2].str)));(yyval.clist_id)->push_back(UnId((yyvsp[0].str)));Block::open(currentblock); }
#line 3027 "lg.tab.cpp"
    break;

  case 135:
#line 619 "lg.ypp"
                   {(yyval.cexp)=0;}
#line 3033 "lg.tab.cpp"
    break;

  case 136:
#line 620 "lg.ypp"
                            {zzzfff->input((yyvsp[0].str));(yyval.cexp)= 0; }
#line 3039 "lg.tab.cpp"
    break;

  case 137:
#line 621 "lg.ypp"
                         {load((yyvsp[0].str));(yyval.cexp)= 0; }
#line 3045 "lg.tab.cpp"
    break;

  case 138:
#line 622 "lg.ypp"
                                             {(yyval.cexp)=Try((yyvsp[-2].cinst),currentblock->close(currentblock,(yyvsp[0].cexp)));}
#line 3051 "lg.tab.cpp"
    break;

  case 139:
#line 623 "lg.ypp"
                     {(yyval.cexp)=(yyvsp[-1].cexp);}
#line 3057 "lg.tab.cpp"
    break;

  case 140:
#line 624 "lg.ypp"
                         {(yyval.cexp)=(yyvsp[0].cexp);}
#line 3063 "lg.tab.cpp"
    break;

  case 141:
#line 626 "lg.ypp"
                {(yyvsp[-1].cexp)=ForAll(currentblock,(yyvsp[-3].clist_id),(yyvsp[-1].cexp));}
#line 3069 "lg.tab.cpp"
    break;

  case 142:
#line 627 "lg.ypp"
                         {
                    inloopcount--;
                    (yyval.cexp)=Block::close(currentblock,C_F0(ForAll((yyvsp[-3].cexp),(yyvsp[0].cexp))));
                 }
#line 3078 "lg.tab.cpp"
    break;

  case 143:
#line 631 "lg.ypp"
                                                                 {inloopcount--; (yyval.cexp)=For((yyvsp[-6].cexp),(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3084 "lg.tab.cpp"
    break;

  case 144:
#line 633 "lg.ypp"
                {inloopcount--;
                 (yyval.cexp)=Block::close(currentblock,C_F0(For((yyvsp[-6].cexp),(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp))));
                }
#line 3092 "lg.tab.cpp"
    break;

  case 145:
#line 636 "lg.ypp"
                                                {inloopcount--;(yyval.cexp)=While((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3098 "lg.tab.cpp"
    break;

  case 146:
#line 637 "lg.ypp"
                                           {(yyval.cexp)=FIf((yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3104 "lg.tab.cpp"
    break;

  case 147:
#line 638 "lg.ypp"
                                                            {(yyval.cexp)=FIf((yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3110 "lg.tab.cpp"
    break;

  case 148:
#line 639 "lg.ypp"
                                    { /* [[begin:]] [[end:]] */
             (yyvsp[-1].cinst).setclose((yyvsp[0].endb));
             (yyval.cexp)=(yyvsp[-1].cinst);
                    //  $$=C_F0(new E_block($2,$3),atype<void>());
         }
#line 3120 "lg.tab.cpp"
    break;

  case 149:
#line 644 "lg.ypp"
                                     { /* <<BORDER_ID>> */
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[-1].str),C_F0(TheOperators,"[border]",(yyvsp[0].args)));}
#line 3127 "lg.tab.cpp"
    break;

  case 150:
#line 646 "lg.ypp"
                                           {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[-4].str),C_F0(TheOperators,"[border]",(yyvsp[-2].args)));}
#line 3134 "lg.tab.cpp"
    break;

  case 151:
#line 649 "lg.ypp"
                      {
                    if(inloopcount)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_break),atype<void>());
                    else lgerror("break not in loop");}
#line 3143 "lg.tab.cpp"
    break;

  case 152:
#line 653 "lg.ypp"
                         {
                    if(inloopcount)
                        (yyval.cexp)= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");}
#line 3152 "lg.tab.cpp"
    break;

  case 153:
#line 657 "lg.ypp"
                             {
                    if (kkembtype>=0)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_return,(rettype[kkembtype]->CastTo((yyvsp[-1].cexp))).OnReturn()) ,atype<void>());
                     else lgerror(" return not in routine ");}
#line 3161 "lg.tab.cpp"
    break;

  case 154:
#line 664 "lg.ypp"
                                         {(yyval.cexp) =  (yyvsp[0].cexp); }
#line 3167 "lg.tab.cpp"
    break;

  case 155:
#line 667 "lg.ypp"
                                     {
   Block::open(currentblock);
   (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[-5].str),atype<double*>());
   (yyval.args)+= (yyvsp[-3].cexp);
   (yyval.args)+= (yyvsp[-1].cexp);
   (yyval.args)+= currentblock->NewVar<LocalVariable>("IndexBorder",atype<long*>());}
#line 3178 "lg.tab.cpp"
    break;

  case 156:
#line 675 "lg.ypp"
                                            {
    Block::open(currentblock);
    (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[-7].str),atype<double*>());
    (yyval.args)+= (yyvsp[-5].cexp);
    (yyval.args)+= (yyvsp[-3].cexp);
    (yyval.args)+= currentblock->NewVar<LocalVariable>((yyvsp[-1].str),atype<long*>());}
#line 3189 "lg.tab.cpp"
    break;

  case 157:
#line 683 "lg.ypp"
                                  {
    //currentblock->close(currentblock;);
   (yyval.args) = ((yyvsp[-1].args) += currentblock->close(currentblock,(yyvsp[0].cexp)));
   }
#line 3198 "lg.tab.cpp"
    break;

  case 159:
#line 691 "lg.ypp"
                  {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3204 "lg.tab.cpp"
    break;

  case 166:
#line 704 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3210 "lg.tab.cpp"
    break;

  case 167:
#line 705 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"+=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3216 "lg.tab.cpp"
    break;

  case 168:
#line 706 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"-=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3222 "lg.tab.cpp"
    break;

  case 169:
#line 707 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"*=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3228 "lg.tab.cpp"
    break;

  case 170:
#line 708 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"/=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3234 "lg.tab.cpp"
    break;

  case 171:
#line 709 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,".*=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3240 "lg.tab.cpp"
    break;

  case 172:
#line 710 "lg.ypp"
                                       {(yyval.cexp)=C_F0(TheOperators,"./=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3246 "lg.tab.cpp"
    break;

  case 174:
#line 716 "lg.ypp"
                                                            {(yyval.cexp)=C_F0(TheOperators,"?:",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3252 "lg.tab.cpp"
    break;

  case 175:
#line 717 "lg.ypp"
                                                            {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3258 "lg.tab.cpp"
    break;

  case 176:
#line 718 "lg.ypp"
                                                            {(yyval.cexp)=C_F0(TheOperators,"::",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3264 "lg.tab.cpp"
    break;

  case 178:
#line 723 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3270 "lg.tab.cpp"
    break;

  case 179:
#line 724 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3276 "lg.tab.cpp"
    break;

  case 180:
#line 725 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3282 "lg.tab.cpp"
    break;

  case 181:
#line 726 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3288 "lg.tab.cpp"
    break;

  case 182:
#line 727 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3294 "lg.tab.cpp"
    break;

  case 183:
#line 728 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3300 "lg.tab.cpp"
    break;

  case 184:
#line 729 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3306 "lg.tab.cpp"
    break;

  case 185:
#line 730 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3312 "lg.tab.cpp"
    break;

  case 186:
#line 731 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3318 "lg.tab.cpp"
    break;

  case 187:
#line 732 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3324 "lg.tab.cpp"
    break;

  case 188:
#line 733 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3330 "lg.tab.cpp"
    break;

  case 189:
#line 734 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3336 "lg.tab.cpp"
    break;

  case 190:
#line 735 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3342 "lg.tab.cpp"
    break;

  case 191:
#line 736 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3348 "lg.tab.cpp"
    break;

  case 192:
#line 737 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3354 "lg.tab.cpp"
    break;

  case 193:
#line 738 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3360 "lg.tab.cpp"
    break;

  case 194:
#line 739 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3366 "lg.tab.cpp"
    break;

  case 195:
#line 740 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3372 "lg.tab.cpp"
    break;

  case 196:
#line 741 "lg.ypp"
                                             {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3378 "lg.tab.cpp"
    break;

  case 197:
#line 745 "lg.ypp"
                                                                    {(yyval.cexp)=(yyvsp[0].cexp);}
#line 3384 "lg.tab.cpp"
    break;

  case 198:
#line 746 "lg.ypp"
                                                                {(yyval.cexp)=C_F0(TheOperators,":");}
#line 3390 "lg.tab.cpp"
    break;

  case 199:
#line 747 "lg.ypp"
                                                                {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3396 "lg.tab.cpp"
    break;

  case 200:
#line 748 "lg.ypp"
                                                                {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3402 "lg.tab.cpp"
    break;

  case 201:
#line 753 "lg.ypp"
                                 {
      (yyval.args) = (yyvsp[-2].cexp);
      (yyval.args) += (yyvsp[0].args); }
#line 3410 "lg.tab.cpp"
    break;

  case 202:
#line 757 "lg.ypp"
                                                            {(yyval.args) = 0;}
#line 3416 "lg.tab.cpp"
    break;

  case 203:
#line 758 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3422 "lg.tab.cpp"
    break;

  case 204:
#line 759 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3428 "lg.tab.cpp"
    break;

  case 205:
#line 760 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3434 "lg.tab.cpp"
    break;

  case 206:
#line 761 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3440 "lg.tab.cpp"
    break;

  case 207:
#line 762 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3446 "lg.tab.cpp"
    break;

  case 208:
#line 763 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3452 "lg.tab.cpp"
    break;

  case 209:
#line 764 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3458 "lg.tab.cpp"
    break;

  case 210:
#line 765 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3464 "lg.tab.cpp"
    break;

  case 211:
#line 766 "lg.ypp"
                                                            {(yyval.args) = Find((yyvsp[0].str));}
#line 3470 "lg.tab.cpp"
    break;

  case 212:
#line 767 "lg.ypp"
                                                            {(yyval.args) = make_pair<const char *,const C_F0>((const char *) (yyvsp[-2].str),(C_F0) (yyvsp[0].cexp));}
#line 3476 "lg.tab.cpp"
    break;

  case 213:
#line 768 "lg.ypp"
                                                            {(yyval.args) = (yyvsp[0].cexp);}
#line 3482 "lg.tab.cpp"
    break;

  case 214:
#line 769 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3488 "lg.tab.cpp"
    break;

  case 215:
#line 770 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3494 "lg.tab.cpp"
    break;

  case 216:
#line 771 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3500 "lg.tab.cpp"
    break;

  case 217:
#line 772 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3506 "lg.tab.cpp"
    break;

  case 218:
#line 773 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3512 "lg.tab.cpp"
    break;

  case 219:
#line 774 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3518 "lg.tab.cpp"
    break;

  case 220:
#line 775 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3524 "lg.tab.cpp"
    break;

  case 221:
#line 776 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3530 "lg.tab.cpp"
    break;

  case 222:
#line 777 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3536 "lg.tab.cpp"
    break;

  case 223:
#line 778 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp));}
#line 3542 "lg.tab.cpp"
    break;

  case 224:
#line 781 "lg.ypp"
                                                            {(yyval.args) = ((yyvsp[-4].args)+= make_pair<const char *,const C_F0>((const char *)(yyvsp[-2].str),(C_F0) (yyvsp[0].cexp)));}
#line 3548 "lg.tab.cpp"
    break;

  case 225:
#line 784 "lg.ypp"
                                 {(yyval.args)=(yyvsp[0].cexp);}
#line 3554 "lg.tab.cpp"
    break;

  case 226:
#line 785 "lg.ypp"
                                 {(yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp));}
#line 3560 "lg.tab.cpp"
    break;

  case 227:
#line 788 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3566 "lg.tab.cpp"
    break;

  case 228:
#line 789 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3572 "lg.tab.cpp"
    break;

  case 229:
#line 790 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3578 "lg.tab.cpp"
    break;

  case 230:
#line 791 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3584 "lg.tab.cpp"
    break;

  case 231:
#line 792 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3590 "lg.tab.cpp"
    break;

  case 232:
#line 793 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3596 "lg.tab.cpp"
    break;

  case 233:
#line 794 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3602 "lg.tab.cpp"
    break;

  case 234:
#line 795 "lg.ypp"
                                   { (yyval.args)=Find((yyvsp[0].str));}
#line 3608 "lg.tab.cpp"
    break;

  case 235:
#line 796 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3614 "lg.tab.cpp"
    break;

  case 236:
#line 797 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3620 "lg.tab.cpp"
    break;

  case 237:
#line 798 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3626 "lg.tab.cpp"
    break;

  case 238:
#line 799 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3632 "lg.tab.cpp"
    break;

  case 239:
#line 800 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3638 "lg.tab.cpp"
    break;

  case 240:
#line 801 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3644 "lg.tab.cpp"
    break;

  case 241:
#line 802 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3650 "lg.tab.cpp"
    break;

  case 242:
#line 803 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3656 "lg.tab.cpp"
    break;

  case 243:
#line 804 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3662 "lg.tab.cpp"
    break;

  case 244:
#line 805 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3668 "lg.tab.cpp"
    break;

  case 245:
#line 806 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3674 "lg.tab.cpp"
    break;

  case 246:
#line 807 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3680 "lg.tab.cpp"
    break;

  case 247:
#line 808 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3686 "lg.tab.cpp"
    break;

  case 248:
#line 809 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3692 "lg.tab.cpp"
    break;

  case 249:
#line 810 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3698 "lg.tab.cpp"
    break;

  case 250:
#line 811 "lg.ypp"
                                   { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str)));}
#line 3704 "lg.tab.cpp"
    break;

  case 252:
#line 816 "lg.ypp"
                              {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[0].cexp));}
#line 3710 "lg.tab.cpp"
    break;

  case 254:
#line 821 "lg.ypp"
                                    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3716 "lg.tab.cpp"
    break;

  case 255:
#line 822 "lg.ypp"
                                    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
#line 3722 "lg.tab.cpp"
    break;

  case 257:
#line 827 "lg.ypp"
                                   {(yyval.cexp)=C_F0(TheOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
#line 3728 "lg.tab.cpp"
    break;

  case 258:
#line 835 "lg.ypp"
                        {(yyval.cexp)=Find((yyvsp[0].str));}
#line 3734 "lg.tab.cpp"
    break;

  case 259:
#line 839 "lg.ypp"
                        {(yyval.cexp)= CConstant((yyvsp[0].lnum));}
#line 3740 "lg.tab.cpp"
    break;

  case 260:
#line 840 "lg.ypp"
                        {(yyval.cexp)= CConstant((yyvsp[0].dnum));}
#line 3746 "lg.tab.cpp"
    break;

  case 261:
#line 841 "lg.ypp"
                        {(yyval.cexp)= CConstant(complex<double>(0,(yyvsp[0].dnum)));}
#line 3752 "lg.tab.cpp"
    break;

  case 262:
#line 842 "lg.ypp"
                  {(yyval.cexp)= CConstant<const char *>((yyvsp[0].str));}
#line 3758 "lg.tab.cpp"
    break;

  case 263:
#line 847 "lg.ypp"
                                                            {(yyval.cexp)=C_F0((yyvsp[-3].cexp),(yyvsp[-2].oper),(yyvsp[-1].args));}
#line 3764 "lg.tab.cpp"
    break;

  case 264:
#line 849 "lg.ypp"
                                              {(yyval.cexp)=C_F0((yyvsp[-3].cexp),(yyvsp[-2].oper),(yyvsp[-1].cexp));}
#line 3770 "lg.tab.cpp"
    break;

  case 265:
#line 850 "lg.ypp"
                                                                {(yyval.cexp)=C_F0((yyvsp[-5].cexp),(yyvsp[-4].oper),(yyvsp[-3].cexp),(yyvsp[-1].cexp));}
#line 3776 "lg.tab.cpp"
    break;

  case 266:
#line 851 "lg.ypp"
                                   {(yyval.cexp)=C_F0((yyvsp[-2].cexp),"[]");}
#line 3782 "lg.tab.cpp"
    break;

  case 267:
#line 852 "lg.ypp"
                                 { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].str)) ;}
#line 3788 "lg.tab.cpp"
    break;

  case 268:
#line 853 "lg.ypp"
                                 { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3794 "lg.tab.cpp"
    break;

  case 269:
#line 854 "lg.ypp"
                                          { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3800 "lg.tab.cpp"
    break;

  case 270:
#line 855 "lg.ypp"
                                  { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3806 "lg.tab.cpp"
    break;

  case 271:
#line 856 "lg.ypp"
                                           { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3812 "lg.tab.cpp"
    break;

  case 272:
#line 857 "lg.ypp"
                                  { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3818 "lg.tab.cpp"
    break;

  case 273:
#line 858 "lg.ypp"
                                           { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3824 "lg.tab.cpp"
    break;

  case 274:
#line 859 "lg.ypp"
                                  { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3830 "lg.tab.cpp"
    break;

  case 275:
#line 860 "lg.ypp"
                                           { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3836 "lg.tab.cpp"
    break;

  case 276:
#line 861 "lg.ypp"
                                  { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3842 "lg.tab.cpp"
    break;

  case 277:
#line 862 "lg.ypp"
                                           { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3848 "lg.tab.cpp"
    break;

  case 278:
#line 863 "lg.ypp"
                                   { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3854 "lg.tab.cpp"
    break;

  case 279:
#line 864 "lg.ypp"
                                            { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3860 "lg.tab.cpp"
    break;

  case 280:
#line 865 "lg.ypp"
                                   { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3866 "lg.tab.cpp"
    break;

  case 281:
#line 866 "lg.ypp"
                                            { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3872 "lg.tab.cpp"
    break;

  case 282:
#line 867 "lg.ypp"
                                   { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;}
#line 3878 "lg.tab.cpp"
    break;

  case 283:
#line 868 "lg.ypp"
                                            { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;}
#line 3884 "lg.tab.cpp"
    break;

  case 284:
#line 869 "lg.ypp"
                                 {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
#line 3890 "lg.tab.cpp"
    break;

  case 285:
#line 870 "lg.ypp"
                                 {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
#line 3896 "lg.tab.cpp"
    break;

  case 286:
#line 871 "lg.ypp"
                                         {
             if ((yyvsp[-3].type)->right()->CastingFrom((yyvsp[-1].cexp).left()) )
                (yyval.cexp)=(yyvsp[-3].type)->right()->CastTo((yyvsp[-1].cexp))  ;
             else { (yyval.cexp)=(yyvsp[-3].type)->right()->Find("<--",basicAC_F0_wa((yyvsp[-1].cexp)));
             if (!(yyval.cexp).left()) { cerr << " no wait to change " << (yyvsp[-1].cexp).left()->right()->name() << " in " <<
                                        (yyvsp[-3].type)->right()->name() << endl;
                                CompileError(" Error in type(exp) "); }
             }
            }
#line 3910 "lg.tab.cpp"
    break;

  case 287:
#line 880 "lg.ypp"
                                        {
           { (yyval.cexp)=(yyvsp[-3].type)->right()->Find("<--",basicAC_F0_wa((yyvsp[-1].args)));
           if (!(yyval.cexp).left()) { cerr << " no wait to change (args) in " <<
                                      (yyvsp[-3].type)->right()->name() << endl;
                              CompileError(" Error in type(exp) "); }
           }
          }
#line 3922 "lg.tab.cpp"
    break;

  case 288:
#line 888 "lg.ypp"
                        {(yyval.cexp)=(yyvsp[-1].cexp);}
#line 3928 "lg.tab.cpp"
    break;

  case 289:
#line 889 "lg.ypp"
                          { (yyval.cexp)=C_F0(TheOperators,"[]",(yyvsp[-1].args));}
#line 3934 "lg.tab.cpp"
    break;

  case 290:
#line 890 "lg.ypp"
                           { (yyval.cexp)=C_F0(TheOperators,"<>",(yyvsp[-1].args));}
#line 3940 "lg.tab.cpp"
    break;

  case 291:
#line 891 "lg.ypp"
                   { (yyval.cexp)=C_F0(TheOperators,"<>",(yyvsp[0].args));}
#line 3946 "lg.tab.cpp"
    break;


#line 3950 "lg.tab.cpp"

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
#line 895 "lg.ypp"



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
