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
     PLUSEQ = 296,
     MOINSEQ = 297,
     MULEQ = 298,
     DIVEQ = 299,
     ARROW = 300,
     BORDER = 301,
     CURVE = 302,
     SOLVE = 303
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
#define PLUSEQ 296
#define MOINSEQ 297
#define MULEQ 298
#define DIVEQ 299
#define ARROW 300
#define BORDER 301
#define CURVE 302
#define SOLVE 303




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
#line 105 "lg.y"
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
#line 299 "lg.tab.cpp"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 213 of yacc.c.  */
#line 311 "lg.tab.cpp"

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
#define YYFINAL  73
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   821

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  74
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  44
/* YYNRULES -- Number of rules. */
#define YYNRULES  163
/* YYNRULES -- Number of states. */
#define YYNSTATES  360

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   303

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
       2,     2,     2,     2,     2,     2,     2,     2,    73,    69,
      16,     6,    17,    72,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    35,     2,    38,    31,    33,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    70,    10,    71,     2,     2,     2,     2,
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
      65,    66,    67,    68
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short int yyprhs[] =
{
       0,     0,     3,     6,     8,    10,    13,    14,    16,    20,
      23,    27,    30,    34,    38,    42,    48,    54,    59,    65,
      70,    76,    78,    82,    84,    86,    88,    92,    97,   101,
     103,   106,   110,   114,   120,   122,   127,   134,   139,   141,
     146,   150,   154,   161,   167,   172,   179,   181,   186,   188,
     192,   194,   198,   201,   207,   212,   214,   218,   219,   224,
     228,   231,   237,   238,   249,   250,   260,   262,   264,   266,
     268,   269,   273,   275,   277,   280,   283,   289,   292,   294,
     304,   314,   320,   326,   334,   338,   342,   349,   352,   355,
     359,   367,   375,   378,   380,   384,   386,   388,   390,   392,
     394,   396,   400,   404,   408,   412,   416,   418,   424,   426,
     430,   434,   438,   442,   446,   450,   454,   458,   462,   466,
     470,   474,   478,   482,   486,   490,   494,   498,   502,   504,
     506,   510,   516,   517,   519,   523,   525,   529,   533,   539,
     541,   545,   547,   550,   552,   556,   560,   563,   565,   567,
     569,   571,   573,   578,   583,   590,   594,   598,   602,   607,
     610,   613,   618,   622
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      75,     0,    -1,    76,    46,    -1,    77,    -1,   103,    -1,
      77,   103,    -1,    -1,    80,    -1,    80,     6,   109,    -1,
      60,    80,    -1,    60,    12,    80,    -1,    83,    80,    -1,
      83,    12,    80,    -1,    35,    78,    38,    -1,    78,     5,
      80,    -1,    78,     5,    35,    78,    38,    -1,    78,     5,
      80,     6,   109,    -1,    78,     5,    60,    80,    -1,    78,
       5,    60,    12,    80,    -1,    78,     5,    83,    80,    -1,
      78,     5,    83,    12,    80,    -1,    80,    -1,    79,     5,
      80,    -1,    42,    -1,    60,    -1,    42,    -1,    42,     6,
     109,    -1,    42,    34,    82,    37,    -1,    81,     5,    81,
      -1,   110,    -1,    60,    42,    -1,    42,     6,   110,    -1,
      82,     5,   110,    -1,    82,     5,    80,     6,   110,    -1,
      58,    -1,    58,    35,    58,    38,    -1,    58,    35,    58,
       5,    58,    38,    -1,    58,    16,    58,    17,    -1,    42,
      -1,    42,    35,   110,    38,    -1,    42,     6,   110,    -1,
      35,    79,    38,    -1,    35,    79,    38,    35,   110,    38,
      -1,    35,    79,    38,     6,   110,    -1,    42,    34,   110,
      37,    -1,    35,    79,    38,    34,   110,    37,    -1,    60,
      -1,    60,    16,    58,    17,    -1,    85,    -1,    87,     5,
      85,    -1,    84,    -1,    88,     5,    84,    -1,    86,    88,
      -1,    86,    35,    58,    38,    87,    -1,    42,    34,    82,
      37,    -1,    90,    -1,    91,     5,    90,    -1,    -1,    83,
      93,    81,    69,    -1,    43,    91,    69,    -1,    89,    69,
      -1,    59,    42,     6,   107,    69,    -1,    -1,    59,    83,
      42,    34,    78,    37,    94,    70,    77,    71,    -1,    -1,
      59,    42,    34,    78,    37,    95,     6,   109,    69,    -1,
      70,    -1,    71,    -1,    50,    -1,    51,    -1,    -1,    83,
     101,    81,    -1,    55,    -1,    69,    -1,    47,    45,    -1,
      48,    45,    -1,   102,    70,    77,    71,   104,    -1,   107,
      69,    -1,    92,    -1,    98,    34,   107,    69,   107,    69,
     107,    37,   103,    -1,    98,    34,   100,    69,   107,    69,
     107,    37,   103,    -1,    99,    34,   107,    37,   103,    -1,
       3,    34,   107,    37,   103,    -1,     3,    34,   107,    37,
     103,     4,   103,    -1,    96,    77,    97,    -1,    66,    42,
     106,    -1,    66,    42,    35,   114,    38,    69,    -1,    52,
      69,    -1,    53,    69,    -1,    54,   107,    69,    -1,    56,
      34,    36,    36,    36,    37,   103,    -1,    34,    42,     6,
     107,     5,   107,    37,    -1,   105,   103,    -1,   109,    -1,
     107,     5,   107,    -1,    21,    -1,    20,    -1,    27,    -1,
      29,    -1,    28,    -1,   110,    -1,   110,     6,   109,    -1,
     110,    61,   109,    -1,   110,    62,   109,    -1,   110,    63,
     109,    -1,   110,    64,   109,    -1,   111,    -1,   111,    72,
     110,    73,   110,    -1,   115,    -1,   111,    22,   111,    -1,
     111,    26,   111,    -1,   111,    25,   111,    -1,   111,    23,
     111,    -1,   111,    24,   111,    -1,   111,    20,   111,    -1,
     111,    21,   111,    -1,   111,     9,   111,    -1,   111,     8,
     111,    -1,   111,    12,   111,    -1,   111,    13,   111,    -1,
     111,    10,   111,    -1,   111,    11,   111,    -1,   111,    16,
     111,    -1,   111,    19,   111,    -1,   111,    17,   111,    -1,
     111,    18,   111,    -1,   111,    15,   111,    -1,   111,    14,
     111,    -1,   110,    -1,    73,    -1,   110,    73,   110,    -1,
     110,    73,   110,    73,   110,    -1,    -1,    60,    -1,    80,
       6,   110,    -1,   112,    -1,   113,     5,    60,    -1,   113,
       5,   112,    -1,   113,     5,    80,     6,   110,    -1,   109,
      -1,   114,     5,   109,    -1,   116,    -1,   108,   116,    -1,
     117,    -1,   117,    31,   115,    -1,   117,    33,   115,    -1,
     117,    32,    -1,    42,    -1,    39,    -1,    40,    -1,    41,
      -1,    45,    -1,   117,    34,   113,    37,    -1,   117,    35,
     112,    38,    -1,   117,    35,   112,     5,   112,    38,    -1,
     117,    35,    38,    -1,   117,    36,    42,    -1,    60,    36,
      42,    -1,    60,    34,   113,    37,    -1,   117,    29,    -1,
     117,    28,    -1,    58,    34,   107,    37,    -1,    34,   107,
      37,    -1,    35,   114,    38,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short int yyrline[] =
{
       0,   231,   231,   267,   270,   271,   274,   275,   276,   277,
     278,   279,   280,   281,   282,   283,   284,   285,   286,   287,
     288,   291,   292,   295,   295,   297,   298,   299,   301,   308,
     309,   310,   311,   312,   315,   316,   317,   318,   325,   326,
     327,   328,   329,   330,   333,   334,   338,   339,   344,   345,
     347,   348,   350,   351,   355,   359,   360,   363,   363,   364,
     365,   366,   368,   367,   383,   382,   391,   392,   394,   396,
     401,   401,   404,   406,   407,   408,   409,   410,   411,   412,
     413,   417,   418,   419,   420,   422,   424,   427,   431,   435,
     442,   445,   451,   457,   458,   463,   464,   465,   466,   467,
     471,   472,   473,   474,   475,   476,   481,   482,   485,   486,
     487,   488,   489,   490,   491,   492,   493,   494,   495,   496,
     497,   498,   499,   500,   501,   502,   503,   504,   509,   510,
     511,   512,   515,   516,   517,   518,   519,   520,   521,   524,
     525,   529,   530,   533,   534,   535,   536,   540,   541,   542,
     543,   544,   545,   546,   547,   548,   549,   550,   551,   552,
     553,   554,   563,   564
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
  "FUNCTION", "FESPACE", "PLUSEQ", "MOINSEQ", "MULEQ", "DIVEQ", "ARROW",
  "BORDER", "CURVE", "SOLVE", "';'", "'{'", "'}'", "'?'", "':'", "$accept",
  "start", "input", "instructions", "list_of_id_args", "list_of_id1", "id",
  "list_of_dcls", "parameters_list", "type_of_dcl", "ID_space",
  "ID_array_space", "fespace", "spaceIDa", "spaceIDb", "spaceIDs",
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
     295,   296,   297,   298,   299,   300,   301,   302,   303,    59,
     123,   125,    63,    58
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    74,    75,    76,    77,    77,    78,    78,    78,    78,
      78,    78,    78,    78,    78,    78,    78,    78,    78,    78,
      78,    79,    79,    80,    80,    81,    81,    81,    81,    82,
      82,    82,    82,    82,    83,    83,    83,    83,    84,    84,
      84,    84,    84,    84,    85,    85,    86,    86,    87,    87,
      88,    88,    89,    89,    90,    91,    91,    93,    92,    92,
      92,    92,    94,    92,    95,    92,    96,    97,    98,    99,
     101,   100,   102,   103,   103,   103,   103,   103,   103,   103,
     103,   103,   103,   103,   103,   103,   103,   103,   103,   103,
     104,   105,   106,   107,   107,   108,   108,   108,   108,   108,
     109,   109,   109,   109,   109,   109,   110,   110,   111,   111,
     111,   111,   111,   111,   111,   111,   111,   111,   111,   111,
     111,   111,   111,   111,   111,   111,   111,   111,   112,   112,
     112,   112,   113,   113,   113,   113,   113,   113,   113,   114,
     114,   115,   115,   116,   116,   116,   116,   117,   117,   117,
     117,   117,   117,   117,   117,   117,   117,   117,   117,   117,
     117,   117,   117,   117
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
       0,     3,     1,     1,     2,     2,     5,     2,     1,     9,
       9,     5,     5,     7,     3,     3,     6,     2,     2,     3,
       7,     7,     2,     1,     3,     1,     1,     1,     1,     1,
       1,     3,     3,     3,     3,     3,     1,     5,     1,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     1,
       3,     5,     0,     1,     3,     1,     3,     3,     5,     1,
       3,     1,     2,     1,     3,     3,     2,     1,     1,     1,
       1,     1,     4,     4,     6,     3,     3,     3,     4,     2,
       2,     4,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       0,     0,    96,    95,    97,    99,    98,     0,     0,   148,
     149,   150,   147,     0,   151,     0,     0,    68,    69,     0,
       0,     0,    72,    34,     0,    46,     0,    73,    66,     0,
       0,     3,    57,     0,     0,    78,     0,     0,     0,     0,
       4,     0,     0,    93,   100,   106,   108,   141,   143,     0,
       0,     0,     0,   139,     0,     0,    55,     0,    74,    75,
      87,    88,     0,     0,     0,     0,     0,    34,     0,     0,
     132,     0,     0,     1,     2,     5,     0,     0,    38,    50,
      52,    60,     0,     0,     0,     0,     0,    77,   142,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   160,   159,     0,   146,     0,   132,
       0,     0,     0,   162,     0,   163,     0,     0,    59,    89,
       0,     0,     0,     0,     6,     0,     0,   147,   133,   129,
       0,   128,   135,     0,   157,     0,     0,     0,    85,    25,
       0,    23,     0,    24,     0,    21,     0,     0,     0,    67,
      84,    70,     0,     0,     0,     0,    94,   101,   102,   103,
     104,   105,   117,   116,   120,   121,   118,   119,   127,   126,
     122,   124,   125,   123,   114,   115,   109,   112,   113,   111,
     110,     0,   144,   145,     0,   155,     0,   156,     0,   140,
     147,     0,     0,    29,    56,    37,   161,     0,    35,     0,
       6,    24,     0,     7,     0,     6,    47,     0,     0,     0,
     158,     0,     0,    92,     0,     0,     0,    58,     0,     0,
      41,    40,     0,     0,    51,     0,     0,     0,     0,     0,
       0,   152,     0,   153,    82,     0,    30,     0,    54,     0,
      61,     0,     0,     9,     0,    64,     0,     0,    11,     0,
     134,   130,   136,     0,   137,     0,     0,    26,     0,    28,
       0,     0,    48,    53,    22,     0,     0,    39,    71,     0,
       0,    81,     0,    76,   107,     0,     0,    31,    24,     0,
      32,    36,    13,    10,     6,    24,    14,     0,     0,     8,
      12,    62,     0,     0,     0,    86,    27,     0,     0,     0,
      43,     0,     0,     0,     0,   154,    83,     0,     0,     0,
      17,     0,     0,    19,     0,     0,   131,   138,     0,     0,
       0,    49,    42,     0,     0,     0,    33,    15,    18,    16,
      20,     0,     0,    94,     0,    44,     0,     0,     0,    65,
       0,    91,     0,    80,    79,     0,    63,    45,     0,    90
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short int yydefgoto[] =
{
      -1,    29,    30,    31,   212,   154,   213,   150,   202,    32,
      79,   272,    33,   273,    80,    34,    56,    57,    35,    76,
     325,   298,    36,   160,    37,    38,   162,   235,    39,    40,
     283,   147,   148,    41,    42,    43,    44,    45,   142,   143,
      54,    46,    47,    48
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -200
static const short int yypact[] =
{
     464,    -9,  -200,  -200,  -200,  -200,  -200,   616,   616,  -200,
    -200,  -200,  -200,     6,  -200,   -12,    -7,  -200,  -200,   -24,
      22,   616,  -200,   157,   -29,    94,    18,  -200,  -200,   125,
      56,   464,  -200,   129,    45,  -200,   464,    83,    88,    61,
    -200,     0,   155,  -200,    35,   277,  -200,  -200,   229,   616,
     135,   151,    17,  -200,     9,   148,  -200,     3,  -200,  -200,
    -200,  -200,     4,    93,   616,   105,   113,   122,   159,   156,
     535,   174,   175,  -200,  -200,  -200,   176,   -18,   109,  -200,
     212,  -200,   305,   643,   616,   464,   616,  -200,  -200,   616,
     616,   616,   616,   616,   616,   616,   616,   616,   616,   616,
     616,   616,   616,   616,   616,   616,   616,   616,   616,   616,
     616,   616,   616,   616,  -200,  -200,   616,  -200,   616,   535,
     508,   177,    89,  -200,   616,  -200,   670,     6,  -200,  -200,
     223,    95,    21,   616,   114,   195,   224,   236,   112,  -200,
     237,   172,  -200,    96,  -200,   205,   616,   464,  -200,   118,
       5,  -200,   213,  -200,    25,  -200,   616,   616,   133,  -200,
    -200,  -200,   180,     7,    98,   358,  -200,  -200,  -200,  -200,
    -200,  -200,   767,   767,   782,   782,   795,   795,   727,   727,
     202,   202,   202,   202,   211,   211,  -200,  -200,  -200,  -200,
    -200,   181,  -200,  -200,    99,  -200,    26,  -200,   464,  -200,
     244,   124,   100,  -200,  -200,  -200,  -200,   198,  -200,    12,
     114,     8,   102,   260,    11,   114,  -200,   616,   616,   562,
    -200,   261,    27,  -200,   616,   670,   176,  -200,   164,   117,
     115,  -200,   230,   117,  -200,   176,   616,   616,   464,   215,
     616,  -200,   589,  -200,   265,   616,  -200,   697,  -200,   238,
    -200,    29,   117,  -200,   144,  -200,   616,   117,  -200,   103,
    -200,   201,   112,   269,  -200,   616,   208,  -200,   104,  -200,
     117,   245,  -200,   273,  -200,   616,   616,  -200,   276,    13,
      14,  -200,   248,  -200,  -200,   246,   464,  -200,   151,   298,
    -200,  -200,  -200,  -200,   114,    15,   299,    24,   304,  -200,
    -200,  -200,   616,   616,   278,  -200,  -200,    41,   616,   164,
    -200,   274,   616,   616,   275,  -200,  -200,   616,    85,   117,
    -200,   616,   117,  -200,   616,   243,  -200,  -200,   616,   280,
     281,  -200,  -200,   106,   108,   285,  -200,  -200,  -200,  -200,
    -200,   254,   464,   287,   616,  -200,   464,   464,   292,  -200,
     411,  -200,   293,  -200,  -200,   294,  -200,  -200,   464,  -200
};

/* YYPGOTO[NTERM-NUM].  */
static const short int yypgoto[] =
{
    -200,  -200,  -200,   -33,  -199,    59,   -49,  -191,   110,   -22,
     178,    28,  -200,  -200,  -200,  -200,   214,  -200,  -200,  -200,
    -200,  -200,  -200,  -200,  -200,  -200,  -200,  -200,  -200,   -31,
    -200,  -200,  -200,    -6,  -200,    -4,   -64,   664,  -113,   219,
     196,   -79,   301,  -200
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -25
static const short int yytable[] =
{
      75,    52,    68,    82,    53,    86,   141,   196,   127,    86,
     226,   251,    86,    66,   124,    62,   259,    86,    86,    86,
     252,   140,    86,   257,   151,    49,   207,   319,   155,    67,
     229,   242,   124,    58,   254,   269,   322,   192,    59,   193,
     152,    89,   153,   122,   278,    60,   229,   125,    55,   191,
     151,    75,   165,   151,   123,   141,   141,   151,   131,   208,
      72,   161,   203,   230,   243,   266,   151,   292,   153,    87,
     140,   153,   128,   129,   227,   153,   237,   163,   164,   329,
     166,   250,   312,   313,   153,   167,   168,   169,   170,   171,
     254,    61,   231,   232,    86,   318,    90,    91,    92,    93,
      86,   219,    74,    86,   219,   247,   264,   254,   254,   247,
      69,    86,   214,    86,    81,   156,   223,    83,   -24,   133,
     199,   275,    84,   337,   224,    73,   198,   209,    70,   285,
      71,    85,   206,   220,    75,   238,   241,   248,    63,   255,
     301,   306,    53,   346,   157,   347,    70,   134,    71,   210,
     276,   130,   225,   260,   261,   141,   151,    65,    70,   151,
      71,   203,   253,   132,    77,   258,   246,   244,   233,    64,
     263,    78,    67,    63,   211,    78,   284,   153,   141,   294,
     274,   287,   126,   290,   155,    70,   151,    71,   214,     7,
       8,    64,    65,   214,     9,    10,    11,    12,   289,   270,
      14,   135,    67,   293,   295,   296,   271,   281,   300,   145,
     146,   310,   311,    50,   136,    51,   144,   158,   149,   197,
     267,   155,   106,   107,   108,   109,   110,   111,   112,   215,
     279,   280,   297,   108,   109,   110,   111,   112,   326,   327,
     205,   216,   -23,   217,   330,   218,   320,   221,   323,   236,
     245,   228,   299,   336,   240,   316,   249,   114,   115,   304,
     116,   117,   118,   119,   120,   121,   256,   265,   277,   286,
     338,   282,   214,   340,   302,   303,   291,   305,   309,   308,
     352,   226,   314,   328,   315,    94,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,   317,   321,   333,   334,     1,   350,
     324,   335,   332,   342,   344,   353,   354,   339,   345,    75,
     341,   348,   343,   349,   351,     2,     3,   359,   355,   307,
     357,   358,     4,     5,     6,   268,   234,   331,   194,     7,
       8,   204,   222,    88,     9,    10,    11,    12,    13,   113,
      14,     0,    15,    16,     0,    17,    18,    19,    20,    21,
      22,     1,     0,    23,    24,    25,     0,     0,     0,     0,
       0,    26,     0,     0,    27,    28,   159,     0,     2,     3,
       0,     0,     0,     0,     0,     4,     5,     6,     0,     0,
       0,     0,     7,     8,     0,     0,     0,     9,    10,    11,
      12,    13,     0,    14,     0,    15,    16,     0,    17,    18,
      19,    20,    21,    22,     1,     0,    23,    24,    25,     0,
       0,     0,     0,     0,    26,     0,     0,    27,    28,   239,
       0,     2,     3,     0,     0,     0,     0,     0,     4,     5,
       6,     0,     0,     0,     0,     7,     8,     0,     0,     0,
       9,    10,    11,    12,    13,     0,    14,     0,    15,    16,
       0,    17,    18,    19,    20,    21,    22,     1,     0,    23,
      24,    25,     0,     0,     0,     0,     0,    26,     0,     0,
      27,    28,   356,     0,     2,     3,     0,     0,     0,     0,
       0,     4,     5,     6,     0,     0,     0,     0,     7,     8,
       0,     0,     0,     9,    10,    11,    12,    13,     0,    14,
       0,    15,    16,     0,    17,    18,    19,    20,    21,    22,
       0,     0,    23,    24,    25,     0,     0,     0,     2,     3,
      26,     0,     0,    27,    28,     4,     5,     6,     0,     0,
       0,     0,     7,     8,     0,     0,   195,     9,    10,    11,
      12,     0,     0,    14,     0,     2,     3,     0,     0,     0,
       0,     0,     4,     5,     6,     0,    50,     0,    51,     7,
       8,     0,     0,     0,     9,    10,    11,   137,     0,     0,
      14,   139,     2,     3,     0,     0,     0,     0,     0,     4,
       5,     6,     0,    50,     0,   138,     7,     8,     0,     0,
       0,     9,    10,    11,   137,     0,     0,    14,   139,     2,
       3,     0,     0,     0,     0,     0,     4,     5,     6,     0,
      50,     0,   262,     7,     8,     0,     0,     0,     9,    10,
      11,    12,     0,     0,    14,   139,     2,     3,     0,     0,
       0,     0,     0,     4,     5,     6,     0,    50,     0,    51,
       7,     8,     0,     0,     0,     9,    10,    11,    12,     0,
       0,    14,   139,     2,     3,     0,     0,     0,     0,     0,
       4,     5,     6,     0,    50,     0,    51,     7,     8,     0,
       0,     0,     9,    10,    11,    12,     0,     0,    14,     0,
       2,     3,     0,     0,     0,     0,     0,     4,     5,     6,
       0,    23,     0,    51,     7,     8,     0,     0,     0,     9,
      10,    11,   200,     0,     0,    14,     0,     2,     3,     0,
       0,     0,     0,     0,     4,     5,     6,     0,    50,     0,
     201,     7,     8,     0,     0,     0,     9,    10,    11,   137,
       0,     0,    14,   102,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,     0,    50,     0,   288,   172,   173,
     174,   175,   176,   177,   178,   179,   180,   181,   182,   183,
     184,   185,   186,   187,   188,   189,   190,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,    98,    99,   100,   101,   102,   103,
     104,   105,   106,   107,   108,   109,   110,   111,   112,   100,
     101,   102,   103,   104,   105,   106,   107,   108,   109,   110,
     111,   112
};

static const short int yycheck[] =
{
      31,     7,    24,    36,     8,     5,    70,   120,     5,     5,
       5,   210,     5,    42,     5,    21,   215,     5,     5,     5,
      12,    70,     5,    12,    42,    34,     5,    12,    77,    58,
       5,     5,     5,    45,     5,   226,    12,   116,    45,   118,
      58,     6,    60,    49,   235,    69,     5,    38,    42,   113,
      42,    82,    85,    42,    37,   119,   120,    42,    64,    38,
      42,    83,   126,    38,    38,    38,    42,    38,    60,    69,
     119,    60,    69,    69,    69,    60,    69,    83,    84,    38,
      86,    69,    69,    69,    60,    89,    90,    91,    92,    93,
       5,    69,   156,   157,     5,   294,    61,    62,    63,    64,
       5,     5,    46,     5,     5,     5,   219,     5,     5,     5,
      16,     5,   134,     5,    69,     6,   147,    34,     6,     6,
     124,     6,    34,    38,     6,     0,    37,   133,    34,   242,
      36,    70,    37,    37,   165,    37,    37,    37,    16,    37,
      37,    37,   146,    37,    35,    37,    34,    34,    36,    35,
      35,    58,    34,   217,   218,   219,    42,    35,    34,    42,
      36,   225,   211,    58,    35,   214,    42,   198,    35,    34,
     219,    42,    58,    16,    60,    42,   240,    60,   242,    35,
     229,   245,    34,   247,   233,    34,    42,    36,   210,    34,
      35,    34,    35,   215,    39,    40,    41,    42,   247,    35,
      45,    42,    58,   252,    60,   254,    42,   238,   257,    34,
      35,   275,   276,    58,    58,    60,    42,     5,    42,    42,
     224,   270,    20,    21,    22,    23,    24,    25,    26,    34,
     236,   237,   254,    22,    23,    24,    25,    26,   302,   303,
      17,    17,     6,     6,   308,    73,   295,    42,   297,    69,
       6,    38,   256,   317,    73,   286,    58,    28,    29,   265,
      31,    32,    33,    34,    35,    36,     6,     6,    38,     4,
     319,    56,   294,   322,    73,     6,    38,    69,     5,    34,
     344,     5,    34,     5,    38,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,     6,     6,   312,   313,     3,   342,
       6,    36,    38,    70,    34,   346,   347,   321,    37,   350,
     324,    36,   328,    69,    37,    20,    21,   358,    36,   270,
      37,    37,    27,    28,    29,   225,   158,   309,   119,    34,
      35,   127,   146,    42,    39,    40,    41,    42,    43,    72,
      45,    -1,    47,    48,    -1,    50,    51,    52,    53,    54,
      55,     3,    -1,    58,    59,    60,    -1,    -1,    -1,    -1,
      -1,    66,    -1,    -1,    69,    70,    71,    -1,    20,    21,
      -1,    -1,    -1,    -1,    -1,    27,    28,    29,    -1,    -1,
      -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
      42,    43,    -1,    45,    -1,    47,    48,    -1,    50,    51,
      52,    53,    54,    55,     3,    -1,    58,    59,    60,    -1,
      -1,    -1,    -1,    -1,    66,    -1,    -1,    69,    70,    71,
      -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,
      29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    43,    -1,    45,    -1,    47,    48,
      -1,    50,    51,    52,    53,    54,    55,     3,    -1,    58,
      59,    60,    -1,    -1,    -1,    -1,    -1,    66,    -1,    -1,
      69,    70,    71,    -1,    20,    21,    -1,    -1,    -1,    -1,
      -1,    27,    28,    29,    -1,    -1,    -1,    -1,    34,    35,
      -1,    -1,    -1,    39,    40,    41,    42,    43,    -1,    45,
      -1,    47,    48,    -1,    50,    51,    52,    53,    54,    55,
      -1,    -1,    58,    59,    60,    -1,    -1,    -1,    20,    21,
      66,    -1,    -1,    69,    70,    27,    28,    29,    -1,    -1,
      -1,    -1,    34,    35,    -1,    -1,    38,    39,    40,    41,
      42,    -1,    -1,    45,    -1,    20,    21,    -1,    -1,    -1,
      -1,    -1,    27,    28,    29,    -1,    58,    -1,    60,    34,
      35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
      45,    73,    20,    21,    -1,    -1,    -1,    -1,    -1,    27,
      28,    29,    -1,    58,    -1,    60,    34,    35,    -1,    -1,
      -1,    39,    40,    41,    42,    -1,    -1,    45,    73,    20,
      21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,    -1,
      58,    -1,    60,    34,    35,    -1,    -1,    -1,    39,    40,
      41,    42,    -1,    -1,    45,    73,    20,    21,    -1,    -1,
      -1,    -1,    -1,    27,    28,    29,    -1,    58,    -1,    60,
      34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,
      -1,    45,    73,    20,    21,    -1,    -1,    -1,    -1,    -1,
      27,    28,    29,    -1,    58,    -1,    60,    34,    35,    -1,
      -1,    -1,    39,    40,    41,    42,    -1,    -1,    45,    -1,
      20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,
      -1,    58,    -1,    60,    34,    35,    -1,    -1,    -1,    39,
      40,    41,    42,    -1,    -1,    45,    -1,    20,    21,    -1,
      -1,    -1,    -1,    -1,    27,    28,    29,    -1,    58,    -1,
      60,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      -1,    -1,    45,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    -1,    58,    -1,    60,    94,    95,
      96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     3,    20,    21,    27,    28,    29,    34,    35,    39,
      40,    41,    42,    43,    45,    47,    48,    50,    51,    52,
      53,    54,    55,    58,    59,    60,    66,    69,    70,    75,
      76,    77,    83,    86,    89,    92,    96,    98,    99,   102,
     103,   107,   108,   109,   110,   111,   115,   116,   117,    34,
      58,    60,   107,   109,   114,    42,    90,    91,    45,    45,
      69,    69,   107,    16,    34,    35,    42,    58,    83,    16,
      34,    36,    42,     0,    46,   103,    93,    35,    42,    84,
      88,    69,    77,    34,    34,    70,     5,    69,   116,     6,
      61,    62,    63,    64,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    72,    28,    29,    31,    32,    33,    34,
      35,    36,   107,    37,     5,    38,    34,     5,    69,    69,
      58,   107,    58,     6,    34,    42,    58,    42,    60,    73,
      80,   110,   112,   113,    42,    34,    35,   105,   106,    42,
      81,    42,    58,    60,    79,    80,     6,    35,     5,    71,
      97,    83,   100,   107,   107,    77,   107,   109,   109,   109,
     109,   109,   111,   111,   111,   111,   111,   111,   111,   111,
     111,   111,   111,   111,   111,   111,   111,   111,   111,   111,
     111,   110,   115,   115,   113,    38,   112,    42,    37,   109,
      42,    60,    82,   110,    90,    17,    37,     5,    38,   107,
      35,    60,    78,    80,    83,    34,    17,     6,    73,     5,
      37,    42,   114,   103,     6,    34,     5,    69,    38,     5,
      38,   110,   110,    35,    84,   101,    69,    69,    37,    71,
      73,    37,     5,    38,   103,     6,    42,     5,    37,    58,
      69,    78,    12,    80,     5,    37,     6,    12,    80,    78,
     110,   110,    60,    80,   112,     6,    38,   109,    82,    81,
      35,    42,    85,    87,    80,     6,    35,    38,    81,   107,
     107,   103,    56,   104,   110,   112,     4,   110,    60,    80,
     110,    38,    38,    80,    35,    60,    80,    83,    95,   109,
      80,    37,    73,     6,   107,    69,    37,    79,    34,     5,
     110,   110,    69,    69,    34,    38,   103,     6,    78,    12,
      80,     6,    12,    80,     6,    94,   110,   110,     5,    38,
     110,    85,    38,   107,   107,    36,   110,    38,    80,   109,
      80,   109,    70,   107,    34,    37,    37,    37,    36,    69,
      77,    37,   110,   103,   103,    36,    71,    37,    37,   103
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
#line 231 "lg.y"
    {
	            
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

                        if(verbosity)  cout << "times: compile "<< CPUcompile-CPUcompileInit <<"s, execution " <<  CPUtime()-CPUcompile << "s\n";
                        deleteStack(stack);
                        //debugstack.clear() 
                        } 
                        fingraphique();
                        NbPtr = ShowAlloc("end execution -- ",lg1) - NbPtr;
                        
                        if (NbPtr) { cout << " ######## We forget of deleting   " << NbPtr << " Nb pointer,   " <<  lg1-lg0 << "Bytes\n" ;}
  return 0;;}
    break;

  case 4:
#line 270 "lg.y"
    {(yyval.cinst)=(yyvsp[0].cexp);;;;}
    break;

  case 5:
#line 271 "lg.y"
    { (yyval.cinst)= ((yyvsp[-1].cinst)+=(yyvsp[0].cexp)) ;}
    break;

  case 6:
#line 274 "lg.y"
    { (yyval.clist_id)=new ListOfId();;}
    break;

  case 7:
#line 275 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str)));}
    break;

  case 8:
#line 276 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str),(yyvsp[0].cexp))) ;}
    break;

  case 9:
#line 277 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double> **>()));}
    break;

  case 10:
#line 278 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double> **>(),true));}
    break;

  case 11:
#line 279 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-1].type)->right())) ;}
    break;

  case 12:
#line 280 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-2].type),true)) ;}
    break;

  case 13:
#line 281 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id))) ;}
    break;

  case 14:
#line 282 "lg.y"
    { (yyval.clist_id) = (yyvsp[-2].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str))) ;}
    break;

  case 15:
#line 283 "lg.y"
    { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-1].clist_id))) ;}
    break;

  case 16:
#line 284 "lg.y"
    { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[-2].str),(yyvsp[0].cexp))) ;}
    break;

  case 17:
#line 285 "lg.y"
    { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-1].str)),atype<FE<double> **>())) ;}
    break;

  case 18:
#line 286 "lg.y"
    { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),Find((yyvsp[-2].str)),atype<FE<double> **>(),true)) ;}
    break;

  case 19:
#line 287 "lg.y"
    { (yyval.clist_id) = (yyvsp[-3].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-1].type)->right())) ;}
    break;

  case 20:
#line 288 "lg.y"
    { (yyval.clist_id) = (yyvsp[-4].clist_id); (yyval.clist_id)->push_back(UnId((yyvsp[0].str),C_F0(),(yyvsp[-2].type),true)) ;}
    break;

  case 21:
#line 291 "lg.y"
    { (yyval.clist_id) = new ListOfId(); (yyval.clist_id)->push_back(UnId((yyvsp[0].str))); ;}
    break;

  case 22:
#line 292 "lg.y"
    { (yyval.clist_id)=(yyvsp[-2].clist_id)  ; (yyval.clist_id)->push_back(UnId((yyvsp[0].str))); ;}
    break;

  case 25:
#line 297 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[0].str),dcltype);}
    break;

  case 26:
#line 298 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[-2].str),dcltype,(yyvsp[0].cexp));}
    break;

  case 27:
#line 299 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariable>((yyvsp[-3].str),dcltype,(yyvsp[-1].args));
                                              (yyvsp[-1].args).destroy();}
    break;

  case 28:
#line 301 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 29:
#line 308 "lg.y"
    {(yyval.args)=(yyvsp[0].cexp);}
    break;

  case 30:
#line 309 "lg.y"
    {(yyval.args)=Find((yyvsp[-1].str));}
    break;

  case 31:
#line 310 "lg.y"
    { (yyval.args)=make_pair<const char *,const C_F0>((yyvsp[-2].str),(yyvsp[0].cexp));}
    break;

  case 32:
#line 311 "lg.y"
    { (yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp)) ;}
    break;

  case 33:
#line 312 "lg.y"
    { (yyval.args)= ((yyvsp[-4].args)+= make_pair<const char *,const C_F0>((yyvsp[-2].str),(yyvsp[0].cexp)));}
    break;

  case 35:
#line 316 "lg.y"
    {(yyval.type)=TypeArray((yyvsp[-3].type),(yyvsp[-1].type));}
    break;

  case 36:
#line 317 "lg.y"
    {(yyval.type)=TypeArray((yyvsp[-5].type),(yyvsp[-3].type),(yyvsp[-1].type));}
    break;

  case 37:
#line 318 "lg.y"
    {(yyval.type)=TypeTemplate((yyvsp[-3].type),(yyvsp[-1].type));}
    break;

  case 38:
#line 325 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[0].str),currentblock,fespacetype,fespacecomplex); ;}
    break;

  case 39:
#line 326 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[-3].str),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex); ;}
    break;

  case 40:
#line 327 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[-2].str),currentblock,fespacetype,(yyvsp[0].cexp),fespacecomplex) ;}
    break;

  case 41:
#line 328 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[-1].clist_id),currentblock,fespacetype,fespacecomplex) ;}
    break;

  case 42:
#line 329 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[-4].clist_id),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex) ;}
    break;

  case 43:
#line 330 "lg.y"
    { (yyval.cexp) =  NewFEvariable((yyvsp[-3].clist_id),currentblock,fespacetype,(yyvsp[0].cexp),fespacecomplex) ;}
    break;

  case 44:
#line 333 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[-3].str),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex); ;}
    break;

  case 45:
#line 334 "lg.y"
    { (yyval.cexp) =  NewFEarray((yyvsp[-4].clist_id),currentblock,fespacetype,(yyvsp[-1].cexp),fespacecomplex) ;}
    break;

  case 46:
#line 338 "lg.y"
    {fespacecomplex=false;  fespacetype = Find((yyvsp[0].str));;}
    break;

  case 47:
#line 339 "lg.y"
    {
             if ((yyvsp[-1].type) != typevarreal && (yyvsp[-1].type) != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=((yyvsp[-1].type)==typevarcomplex);
             fespacetype = Find((yyvsp[-3].str));;}
    break;

  case 48:
#line 344 "lg.y"
    {  (yyval.cexp) = (yyvsp[0].cexp)  ;}
    break;

  case 49:
#line 345 "lg.y"
    { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));;}
    break;

  case 50:
#line 347 "lg.y"
    {  (yyval.cexp) = (yyvsp[0].cexp)  ;}
    break;

  case 51:
#line 348 "lg.y"
    { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));;}
    break;

  case 52:
#line 350 "lg.y"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[0].cexp);}
    break;

  case 53:
#line 351 "lg.y"
    { (yyval.cexp)=0;  (yyval.cexp) = (yyvsp[0].cexp);}
    break;

  case 54:
#line 356 "lg.y"
    {(yyval.cexp)=currentblock->NewVar<LocalVariableFES,size_t>((yyvsp[-3].str),atype<pfes*>(),(yyvsp[-1].args),dimFESpaceImage((yyvsp[-1].args)));
     (yyvsp[-1].args).destroy(); ;}
    break;

  case 56:
#line 360 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 57:
#line 363 "lg.y"
    {dcltype=(yyvsp[0].type);}
    break;

  case 58:
#line 363 "lg.y"
    {(yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 59:
#line 364 "lg.y"
    {(yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 60:
#line 365 "lg.y"
    { (yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 61:
#line 366 "lg.y"
    {(yyval.cexp)=currentblock->NewID((yyvsp[-4].type),(yyvsp[-3].str),(yyvsp[-1].cexp));;}
    break;

  case 62:
#line 368 "lg.y"
    {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = (yyvsp[-4].type)->right();
                      routineinblock[kkembtype] = currentblock;
                      (yyvsp[-1].routine)=new Routine((yyvsp[-5].type),(yyvsp[-4].type)->right(),(yyvsp[-3].str),(yyvsp[-1].clist_id),currentblock);
                     // cout << " \n after new routine \n " << endl;                      
                      ;}
    break;

  case 63:
#line 376 "lg.y"
    { currentblock=(yyvsp[-5].routine)->Set((yyvsp[-1].cinst));
                       currentblock->Add((yyvsp[-7].str),"(",(yyvsp[-5].routine));
                       kkembtype--;
                       (yyval.cexp)=0;
                    
                        ;}
    break;

  case 64:
#line 383 "lg.y"
    {Block::open(currentblock); (yyvsp[-4].type)->SetArgs((yyvsp[-1].clist_id));;}
    break;

  case 65:
#line 385 "lg.y"
    {  (yyval.cinst)=currentblock->close(currentblock);
                         (yyval.cexp)=currentblock->NewID((yyvsp[-8].type),(yyvsp[-7].str),(yyvsp[-1].cexp),*(yyvsp[-5].clist_id));
                         delete (yyvsp[-5].clist_id); //  FH 23032005
                         ;}
    break;

  case 66:
#line 391 "lg.y"
    {  Block::open(currentblock);}
    break;

  case 67:
#line 392 "lg.y"
    {  (yyval.cexp)=currentblock->close(currentblock);}
    break;

  case 68:
#line 394 "lg.y"
    {ffassert(inloopcount<sizeStackOfLoop);  // modif FH july 2005
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 69:
#line 396 "lg.y"
    {ffassert(inloopcount<sizeStackOfLoop);
                StackOfLoop[inloopcount++]=currentblock;;}
    break;

  case 70:
#line 401 "lg.y"
    {dcltype=(yyvsp[0].type); Block::open(currentblock);  ;}
    break;

  case 71:
#line 402 "lg.y"
    {(yyval.cexp)=(yyvsp[0].cexp);}
    break;

  case 72:
#line 404 "lg.y"
    { Block::open(currentblock) ;}
    break;

  case 73:
#line 406 "lg.y"
    {(yyval.cexp)=0;;}
    break;

  case 74:
#line 407 "lg.y"
    {zzzfff->input((yyvsp[0].str));(yyval.cexp)= 0; ;}
    break;

  case 75:
#line 408 "lg.y"
    {load((yyvsp[0].str));(yyval.cexp)= 0; ;}
    break;

  case 76:
#line 409 "lg.y"
    {(yyval.cexp)=Try((yyvsp[-2].cinst),(yyvsp[0].cexp),currentblock->close(currentblock));;}
    break;

  case 77:
#line 410 "lg.y"
    {(yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 78:
#line 411 "lg.y"
    {(yyval.cexp)=(yyvsp[0].cexp);}
    break;

  case 79:
#line 412 "lg.y"
    {inloopcount--; (yyval.cexp)=For((yyvsp[-6].cexp),(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 80:
#line 414 "lg.y"
    {inloopcount--; 
                (yyval.cexp)=C_F0(For((yyvsp[-6].cexp),(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp)),currentblock->close(currentblock));}
    break;

  case 81:
#line 417 "lg.y"
    {inloopcount--;(yyval.cexp)=While((yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 82:
#line 418 "lg.y"
    {(yyval.cexp)=FIf((yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 83:
#line 419 "lg.y"
    {(yyval.cexp)=FIf((yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 84:
#line 420 "lg.y"
    { 
                      (yyval.cexp)=C_F0(new E_block((yyvsp[-1].cinst),(yyvsp[0].cexp)),atype<void>()) ;}
    break;

  case 85:
#line 422 "lg.y"
    {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[-1].str),C_F0(TheOperators,"[border]",(yyvsp[0].args)));}
    break;

  case 86:
#line 424 "lg.y"
    {
                      (yyval.cexp)=0;currentblock->NewID(atype<const E_Border *>(),(yyvsp[-4].str),C_F0(TheOperators,"[border]",(yyvsp[-2].args)));}
    break;

  case 87:
#line 427 "lg.y"
    {
                    if(inloopcount) 
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;}
    break;

  case 88:
#line 431 "lg.y"
    { 
                    if(inloopcount)
                        (yyval.cexp)= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");}
    break;

  case 89:
#line 435 "lg.y"
    { 
                    if (kkembtype>=0)
                      (yyval.cexp)= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo((yyvsp[-1].cexp))) ,atype<void>());
                     else lgerror(" return not in routine ") ;}
    break;

  case 90:
#line 442 "lg.y"
    {(yyval.cexp) =  (yyvsp[0].cexp); ;}
    break;

  case 91:
#line 445 "lg.y"
    { 
   Block::open(currentblock);
   (yyval.args) = currentblock->NewVar<LocalVariable>((yyvsp[-5].str),atype<double*>());
   (yyval.args)+= (yyvsp[-3].cexp);
   (yyval.args)+= (yyvsp[-1].cexp) ;}
    break;

  case 92:
#line 451 "lg.y"
    {   
   (yyval.args) = ((yyvsp[-1].args) += (yyvsp[0].cexp));
   currentblock->close(currentblock);}
    break;

  case 94:
#line 458 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));;}
    break;

  case 101:
#line 472 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 102:
#line 473 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"+=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 103:
#line 474 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"-=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 104:
#line 475 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"*=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 105:
#line 476 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"/=",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 107:
#line 482 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,"?:",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 109:
#line 486 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 110:
#line 487 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 111:
#line 488 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 112:
#line 489 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 113:
#line 490 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 114:
#line 491 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 115:
#line 492 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 116:
#line 493 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 117:
#line 494 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 118:
#line 495 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 119:
#line 496 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 120:
#line 497 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 121:
#line 498 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 122:
#line 499 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 123:
#line 500 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 124:
#line 501 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 125:
#line 502 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 126:
#line 503 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 127:
#line 504 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 128:
#line 509 "lg.y"
    {(yyval.cexp)=(yyvsp[0].cexp);}
    break;

  case 129:
#line 510 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,":");}
    break;

  case 130:
#line 511 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 131:
#line 512 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,":",(yyvsp[-4].cexp),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 132:
#line 515 "lg.y"
    {(yyval.args)=0;}
    break;

  case 133:
#line 516 "lg.y"
    {(yyval.args)=Find((yyvsp[0].str));}
    break;

  case 134:
#line 517 "lg.y"
    { (yyval.args)=make_pair<const char *,const C_F0>((yyvsp[-2].str),(yyvsp[0].cexp));}
    break;

  case 135:
#line 518 "lg.y"
    {(yyval.args)=(yyvsp[0].cexp);}
    break;

  case 136:
#line 519 "lg.y"
    { (yyval.args) = ((yyvsp[-2].args) += Find((yyvsp[0].str))) ;}
    break;

  case 137:
#line 520 "lg.y"
    { (yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp)) ;}
    break;

  case 138:
#line 521 "lg.y"
    { (yyval.args)= ((yyvsp[-4].args)+= make_pair<const char *,const C_F0>((yyvsp[-2].str),(yyvsp[0].cexp))) ;}
    break;

  case 139:
#line 524 "lg.y"
    {(yyval.args)=(yyvsp[0].cexp);}
    break;

  case 140:
#line 525 "lg.y"
    {(yyval.args) = ((yyvsp[-2].args) += (yyvsp[0].cexp)) ;}
    break;

  case 142:
#line 530 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[0].cexp));}
    break;

  case 144:
#line 534 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 145:
#line 535 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[-1].oper),(yyvsp[-2].cexp),(yyvsp[0].cexp));}
    break;

  case 146:
#line 536 "lg.y"
    {(yyval.cexp)=C_F0(TheOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
    break;

  case 147:
#line 540 "lg.y"
    {(yyval.cexp)=Find((yyvsp[0].str));;}
    break;

  case 148:
#line 541 "lg.y"
    {(yyval.cexp)= CConstant((yyvsp[0].lnum));}
    break;

  case 149:
#line 542 "lg.y"
    {(yyval.cexp)= CConstant((yyvsp[0].dnum));}
    break;

  case 150:
#line 543 "lg.y"
    {(yyval.cexp)= CConstant(complex<double>(0,(yyvsp[0].dnum)));}
    break;

  case 151:
#line 544 "lg.y"
    {(yyval.cexp)= CConstant<const char *>((yyvsp[0].str));}
    break;

  case 152:
#line 545 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-3].cexp),(yyvsp[-2].oper),(yyvsp[-1].args));;}
    break;

  case 153:
#line 546 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-3].cexp),(yyvsp[-2].oper),(yyvsp[-1].cexp));}
    break;

  case 154:
#line 547 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-5].cexp),(yyvsp[-4].oper),(yyvsp[-3].cexp),(yyvsp[-1].cexp));}
    break;

  case 155:
#line 548 "lg.y"
    {(yyval.cexp)=C_F0((yyvsp[-2].cexp),"[]");}
    break;

  case 156:
#line 549 "lg.y"
    { (yyval.cexp)=C_F0((yyvsp[-2].cexp),(yyvsp[0].str)) ;;}
    break;

  case 157:
#line 550 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[-2].str)),(yyvsp[0].str)) ;;}
    break;

  case 158:
#line 551 "lg.y"
    { (yyval.cexp)=C_F0(Find((yyvsp[-3].str)),(yyvsp[-2].oper),(yyvsp[-1].args)) ;;}
    break;

  case 159:
#line 552 "lg.y"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
    break;

  case 160:
#line 553 "lg.y"
    {(yyval.cexp)=C_F0(TheRightOperators,(yyvsp[0].oper),(yyvsp[-1].cexp));}
    break;

  case 161:
#line 554 "lg.y"
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

  case 162:
#line 563 "lg.y"
    {(yyval.cexp)=(yyvsp[-1].cexp);}
    break;

  case 163:
#line 564 "lg.y"
    { (yyval.cexp)=C_F0(TheOperators,"[]",(yyvsp[-1].args));}
    break;


    }

/* Line 1037 of yacc.c.  */
#line 2446 "lg.tab.cpp"

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


#line 569 "lg.y"
 


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
	retvalue=1,cerr <<  "Error:a block is not close" << endl;   
      else {
	  if( verbosity) {
	      cerr << " CodeAlloc : nb ptr  "<< CodeAlloc::nb << ",  size :"  <<  CodeAlloc::lg << endl;
	      cerr <<  "Bien: On a fini Normalement" << endl; }
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

int mainff (int  argc, char **argv)
{
GetEnvironment();     
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
#ifdef HAVE_LIBUMFPACK   
  if(verbosity)  cout << " UMFPACK ";  
#endif
 // callInitsFunct(); Pb opimisation 
  if(verbosity)  cout << endl;
  zzzfff->input(cc);
  EnvironmentLoad(); // just before compile     
  retvalue= Compile(); 
      
#ifdef PARALLELE
  end_parallele();
#endif
  //  currentblock->close(currentblock).eval(thestack);
  fingraphique();
  Destroylex( zzzfff);
  
   // ClearMem();
  return retvalue;
}


 

