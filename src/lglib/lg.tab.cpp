/* A Bison parser, made by GNU Bison 1.875a.  */

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

int ShowAlloc(char *s);
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
#define YYFINAL  66
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   948

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  70
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  41
/* YYNRULES -- Number of rules. */
#define YYNRULES  155
/* YYNRULES -- Number of states. */
#define YYNSTATES  332

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
     103,   106,   110,   114,   120,   122,   127,   134,   136,   141,
     145,   149,   156,   162,   167,   174,   176,   181,   183,   187,
     189,   193,   196,   202,   207,   209,   213,   214,   219,   223,
     226,   232,   233,   244,   245,   255,   257,   259,   261,   263,
     264,   268,   270,   273,   276,   279,   281,   291,   301,   307,
     313,   321,   325,   329,   336,   339,   342,   346,   354,   357,
     359,   363,   365,   367,   369,   371,   373,   375,   379,   383,
     387,   391,   395,   397,   401,   405,   409,   413,   417,   421,
     425,   429,   433,   437,   441,   445,   449,   453,   457,   461,
     465,   469,   473,   475,   477,   481,   487,   488,   490,   494,
     496,   500,   504,   510,   512,   516,   518,   521,   523,   527,
     531,   534,   536,   538,   540,   542,   544,   549,   554,   561,
     565,   569,   572,   575,   580,   584
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
       5,    55,    38,    -1,    42,    -1,    42,    35,   104,    38,
      -1,    42,     6,   104,    -1,    35,    75,    38,    -1,    35,
      75,    38,    35,   104,    38,    -1,    35,    75,    38,     6,
     104,    -1,    42,    34,   104,    37,    -1,    35,    75,    38,
      34,   104,    37,    -1,    57,    -1,    57,    16,    55,    17,
      -1,    81,    -1,    83,     5,    81,    -1,    80,    -1,    84,
       5,    80,    -1,    82,    84,    -1,    82,    35,    55,    38,
      83,    -1,    42,    34,    78,    37,    -1,    86,    -1,    87,
       5,    86,    -1,    -1,    79,    89,    77,    66,    -1,    43,
      87,    66,    -1,    85,    66,    -1,    56,    42,     6,   101,
      66,    -1,    -1,    56,    79,    42,    34,    74,    37,    90,
      67,    73,    68,    -1,    -1,    56,    42,    34,    74,    37,
      91,     6,   103,    66,    -1,    67,    -1,    68,    -1,    50,
      -1,    51,    -1,    -1,    79,    97,    77,    -1,    66,    -1,
      47,    45,    -1,    48,    45,    -1,   101,    66,    -1,    88,
      -1,    94,    34,   101,    66,   101,    66,   101,    37,    98,
      -1,    94,    34,    96,    66,   101,    66,   101,    37,    98,
      -1,    95,    34,   101,    37,    98,    -1,     3,    34,   101,
      37,    98,    -1,     3,    34,   101,    37,    98,     4,    98,
      -1,    92,    73,    93,    -1,    63,    42,   100,    -1,    63,
      42,    35,   107,    38,    66,    -1,    52,    66,    -1,    53,
      66,    -1,    54,   101,    66,    -1,    34,    42,     6,   101,
       5,   101,    37,    -1,    99,    98,    -1,   103,    -1,   101,
       5,   101,    -1,    21,    -1,    20,    -1,    27,    -1,    29,
      -1,    28,    -1,   104,    -1,   104,     6,   103,    -1,   104,
      58,   103,    -1,   104,    59,   103,    -1,   104,    60,   103,
      -1,   104,    61,   103,    -1,   108,    -1,   104,    22,   104,
      -1,   104,    26,   104,    -1,   104,    25,   104,    -1,   104,
      23,   104,    -1,   104,    24,   104,    -1,   104,    20,   104,
      -1,   104,    21,   104,    -1,   104,     9,   104,    -1,   104,
       8,   104,    -1,   104,    12,   104,    -1,   104,    13,   104,
      -1,   104,    10,   104,    -1,   104,    11,   104,    -1,   104,
      16,   104,    -1,   104,    19,   104,    -1,   104,    17,   104,
      -1,   104,    18,   104,    -1,   104,    15,   104,    -1,   104,
      14,   104,    -1,   104,    -1,    69,    -1,   104,    69,   104,
      -1,   104,    69,   104,    69,   104,    -1,    -1,    57,    -1,
      76,     6,   104,    -1,   105,    -1,   106,     5,    57,    -1,
     106,     5,   105,    -1,   106,     5,    76,     6,   104,    -1,
     103,    -1,   107,     5,   103,    -1,   109,    -1,   102,   109,
      -1,   110,    -1,   110,    31,   108,    -1,   110,    33,   108,
      -1,   110,    32,    -1,    42,    -1,    39,    -1,    40,    -1,
      41,    -1,    45,    -1,   110,    34,   106,    37,    -1,   110,
      35,   105,    38,    -1,   110,    35,   105,     5,   105,    38,
      -1,   110,    35,    38,    -1,   110,    36,    42,    -1,   110,
      29,    -1,   110,    28,    -1,    55,    34,   101,    37,    -1,
      34,   101,    37,    -1,    35,   107,    38,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   198,   198,   233,   236,   237,   240,   241,   242,   243,
     244,   245,   246,   247,   248,   249,   250,   251,   252,   253,
     254,   257,   258,   261,   261,   263,   264,   265,   266,   273,
     274,   275,   276,   277,   280,   281,   282,   289,   290,   291,
     292,   293,   294,   297,   298,   302,   303,   308,   309,   311,
     312,   314,   315,   319,   322,   323,   326,   326,   327,   328,
     329,   331,   330,   341,   340,   347,   348,   350,   351,   354,
     354,   357,   358,   359,   360,   361,   362,   363,   367,   368,
     369,   370,   372,   374,   377,   381,   385,   393,   399,   405,
     406,   411,   412,   413,   414,   415,   419,   420,   421,   422,
     423,   424,   428,   429,   430,   431,   432,   433,   434,   435,
     436,   437,   438,   439,   440,   441,   442,   443,   444,   445,
     446,   447,   452,   453,   454,   455,   458,   459,   460,   461,
     462,   463,   464,   467,   468,   472,   473,   476,   477,   478,
     479,   483,   484,   485,   486,   487,   488,   489,   490,   491,
     492,   493,   494,   495,   504,   505
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
      78,    78,    78,    78,    79,    79,    79,    80,    80,    80,
      80,    80,    80,    81,    81,    82,    82,    83,    83,    84,
      84,    85,    85,    86,    87,    87,    89,    88,    88,    88,
      88,    90,    88,    91,    88,    92,    93,    94,    95,    97,
      96,    98,    98,    98,    98,    98,    98,    98,    98,    98,
      98,    98,    98,    98,    98,    98,    98,    99,   100,   101,
     101,   102,   102,   102,   102,   102,   103,   103,   103,   103,
     103,   103,   104,   104,   104,   104,   104,   104,   104,   104,
     104,   104,   104,   104,   104,   104,   104,   104,   104,   104,
     104,   104,   105,   105,   105,   105,   106,   106,   106,   106,
     106,   106,   106,   107,   107,   108,   108,   109,   109,   109,
     109,   110,   110,   110,   110,   110,   110,   110,   110,   110,
     110,   110,   110,   110,   110,   110
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     2,     1,     1,     2,     0,     1,     3,     2,
       3,     2,     3,     3,     3,     5,     5,     4,     5,     4,
       5,     1,     3,     1,     1,     1,     3,     4,     3,     1,
       2,     3,     3,     5,     1,     4,     6,     1,     4,     3,
       3,     6,     5,     4,     6,     1,     4,     1,     3,     1,
       3,     2,     5,     4,     1,     3,     0,     4,     3,     2,
       5,     0,    10,     0,     9,     1,     1,     1,     1,     0,
       3,     1,     2,     2,     2,     1,     9,     9,     5,     5,
       7,     3,     3,     6,     2,     2,     3,     7,     2,     1,
       3,     1,     1,     1,     1,     1,     1,     3,     3,     3,
       3,     3,     1,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     1,     3,     5,     0,     1,     3,     1,
       3,     3,     5,     1,     3,     1,     2,     1,     3,     3,
       2,     1,     1,     1,     1,     1,     4,     4,     6,     3,
       3,     2,     2,     4,     3,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       0,     0,    92,    91,    93,    95,    94,     0,     0,   142,
     143,   144,   141,     0,   145,     0,     0,    67,    68,     0,
       0,     0,    34,     0,    45,     0,    71,    65,     0,     0,
       3,    56,     0,     0,    75,     0,     0,     0,     4,     0,
       0,    89,    96,   102,   135,   137,     0,     0,     0,   133,
       0,     0,    54,     0,    72,    73,    84,    85,     0,     0,
       0,     0,    34,     0,     0,     0,     1,     2,     5,     0,
       0,    37,    49,    51,    59,     0,     0,     0,     0,    74,
     136,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   152,   151,     0,   140,     0,
     126,     0,     0,     0,   154,     0,   155,     0,     0,    58,
      86,     0,     0,     0,     6,     0,     0,     0,     0,     0,
      82,    25,     0,    23,     0,    24,     0,    21,     0,     0,
       0,    66,    81,    69,     0,     0,     0,    90,    97,   111,
     110,   114,   115,   112,   113,   121,   120,   116,   118,   119,
     117,   108,   109,   103,   106,   107,   105,   104,    98,    99,
     100,   101,   138,   139,   141,   127,   123,     0,   122,   129,
       0,   149,     0,   150,     0,   134,   141,     0,     0,    29,
      55,   153,     0,    35,     0,     6,    24,     0,     7,     0,
       6,    46,     0,     0,    88,     0,     0,     0,    57,     0,
       0,    40,    39,     0,     0,    50,     0,     0,     0,     0,
       0,     0,     0,   146,     0,   147,    79,     0,    30,     0,
      53,     0,    60,     0,     0,     9,     0,    63,     0,     0,
      11,     0,     0,     0,    26,     0,    28,     0,     0,    47,
      52,    22,     0,     0,    38,    70,     0,     0,    78,   128,
     124,   130,     0,   131,     0,     0,    31,     0,    32,    36,
      13,    10,     6,    24,    14,     0,     0,     8,    12,    61,
       0,    83,    27,     0,     0,     0,    42,     0,     0,     0,
       0,     0,   148,    80,     0,     0,     0,    17,     0,     0,
      19,     0,     0,     0,     0,     0,    48,    41,     0,     0,
     125,   132,    33,    15,    18,    16,    20,     0,     0,    90,
       0,    43,     0,     0,    64,     0,    87,     0,    77,    76,
      62,    44
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,    28,    29,    30,   197,   136,   198,   132,   188,    31,
      72,   249,    32,   250,    73,    33,    52,    53,    34,    69,
     302,   276,    35,   142,    36,    37,   144,   216,    38,   129,
     130,    39,    40,    41,    42,   179,   180,    50,    43,    44,
      45
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -192
static const short yypact[] =
{
     451,   -22,  -192,  -192,  -192,  -192,  -192,   697,   697,  -192,
    -192,  -192,  -192,    21,  -192,    -3,     5,  -192,  -192,   -21,
      -1,   697,   142,    12,    52,    69,  -192,  -192,    88,    92,
     451,  -192,   128,    68,  -192,   451,   106,   112,  -192,     0,
     154,  -192,   513,  -192,  -192,   169,   697,   115,    55,  -192,
       3,   126,  -192,     6,  -192,  -192,  -192,  -192,     8,   697,
      99,   125,   130,   132,   124,   147,  -192,  -192,  -192,   143,
      82,   109,  -192,   182,  -192,   288,   720,   697,   697,  -192,
    -192,   697,   697,   697,   697,   697,   697,   697,   697,   697,
     697,   697,   697,   697,   697,   697,   697,   697,   697,   697,
     697,   697,   697,   697,   697,  -192,  -192,   697,  -192,   697,
     520,   558,   148,    75,  -192,   697,  -192,   635,    21,  -192,
    -192,    77,    35,   697,    93,   158,   208,   184,   697,   451,
    -192,   127,    29,  -192,   189,  -192,    39,  -192,   697,   697,
     131,  -192,  -192,  -192,   162,    30,    85,  -192,  -192,   907,
     907,   922,   922,   227,   227,   366,   366,   192,   192,   192,
     192,   198,   198,  -192,  -192,  -192,  -192,  -192,  -192,  -192,
    -192,  -192,  -192,  -192,   224,   225,  -192,   228,     7,  -192,
     104,  -192,    46,  -192,   451,  -192,   230,   191,   105,   890,
    -192,  -192,   183,  -192,    31,    93,    36,   108,   233,    37,
      93,  -192,   248,    48,  -192,   697,   635,   143,  -192,   133,
       4,   117,   890,   768,     4,  -192,   143,   697,   697,   451,
     697,   697,   581,  -192,   612,  -192,   252,   697,  -192,   666,
    -192,   238,  -192,    51,     4,  -192,   129,  -192,   697,     4,
    -192,   114,   697,   211,  -192,   116,  -192,     4,   244,  -192,
     274,  -192,   697,   697,  -192,   275,    32,    33,  -192,   890,
     351,   225,   277,  -192,   249,   451,   890,   282,   890,  -192,
    -192,  -192,    93,    43,   283,    45,   284,  -192,  -192,  -192,
     293,  -192,  -192,    53,   697,   133,   890,   799,   697,   697,
     697,   697,  -192,  -192,   697,    70,     4,  -192,   697,     4,
    -192,   697,   226,   697,   265,   830,  -192,  -192,   120,   121,
     890,   890,   890,  -192,  -192,  -192,  -192,   235,   451,   267,
     697,  -192,   451,   451,  -192,   390,  -192,   860,  -192,  -192,
    -192,  -192
};

/* YYPGOTO[NTERM-NUM].  */
static const short yypgoto[] =
{
    -192,  -192,  -192,   -34,  -191,    58,   -67,   -80,    96,   -17,
     166,    22,  -192,  -192,  -192,  -192,   193,  -192,  -192,  -192,
    -192,  -192,  -192,  -192,  -192,  -192,  -192,  -192,   -28,  -192,
    -192,    -7,  -192,     2,   175,  -104,  -192,   190,   -45,   270,
    -192
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -25
static const short yytable[] =
{
      48,    75,    68,   137,   233,    78,    63,   182,   115,   241,
      49,   118,    46,    78,    58,    82,    83,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
      97,    98,    99,   100,   207,    78,    78,    78,    78,   113,
     192,   116,    54,   177,   210,    56,   133,    68,   234,   239,
      55,   224,   121,   115,    61,   296,   236,   299,   210,   143,
      78,   135,   172,    51,   173,    57,    79,    62,    64,   145,
     146,   147,   119,   193,   120,   236,   221,   211,   133,   133,
      78,   295,    78,   148,   225,   133,   243,   133,    66,   270,
      78,   304,   114,   135,   135,   208,   218,   232,   288,   289,
     135,   204,   135,   168,   169,   170,   171,   199,   313,   222,
     229,    65,   184,   236,   191,   138,   194,   185,   263,   236,
     264,   229,   219,   252,   133,    78,    78,   246,   195,   235,
      49,   123,   240,   205,    74,   133,   255,   134,    67,   135,
      76,   223,   230,   251,   139,   237,    77,   137,    62,    59,
     196,   279,   253,   282,   122,   262,   226,   322,   323,   124,
     117,   206,   267,    70,   272,    60,   214,   271,   247,   274,
      71,   133,   278,    71,   125,   248,    59,    60,   199,   126,
     137,   127,   128,   199,    62,   131,   273,   140,     7,     8,
     183,   258,   200,     9,    10,    11,    12,   105,   106,    14,
     107,   108,   109,   110,   111,   112,   297,   244,   300,    47,
     256,   257,    94,    95,    96,    97,    98,    99,   100,   275,
      96,    97,    98,    99,   100,   201,   202,   209,   217,   314,
     -23,   -24,   316,   228,   220,   280,   227,   293,   231,   238,
     277,    88,    89,    90,    91,    92,    93,    94,    95,    96,
      97,    98,    99,   100,   242,   199,   265,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   161,
     162,   163,   164,   165,   166,   167,   269,   281,   284,   285,
     207,   308,   309,   291,   325,   178,   178,   292,   294,   298,
     301,     1,   189,   318,   328,   329,   319,    68,   303,   320,
     315,   324,   245,   317,   326,   283,   215,   306,     2,     3,
      80,   190,     0,   212,   213,     4,     5,     6,   203,     0,
       0,     0,     7,     8,     0,     0,     0,     9,    10,    11,
      12,    13,     0,    14,     0,    15,    16,     0,    17,    18,
      19,    20,    21,    22,    23,    24,     0,     0,     0,     0,
       0,    25,     0,     0,    26,    27,   141,     0,     0,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    92,
      93,    94,    95,    96,    97,    98,    99,   100,     0,     0,
       0,   189,    90,    91,    92,    93,    94,    95,    96,    97,
      98,    99,   100,     1,     0,   259,   260,   178,     0,   178,
       0,     0,   266,     0,   268,     0,     0,     0,     0,     0,
       2,     3,     0,     0,     0,     0,     0,     4,     5,     6,
     290,     0,     0,     0,     7,     8,     0,   286,   287,     9,
      10,    11,    12,    13,     0,    14,     0,    15,    16,     0,
      17,    18,    19,    20,    21,    22,    23,    24,     0,     0,
       0,     0,     0,    25,     1,     0,    26,    27,   330,   305,
       0,     0,     0,     0,     0,   310,   311,     0,     0,   312,
       0,     2,     3,     0,     0,     0,     0,     0,     4,     5,
       6,     0,     0,     0,     0,     7,     8,     0,     0,     0,
       9,    10,    11,    12,    13,   327,    14,     0,    15,    16,
       0,    17,    18,    19,    20,    21,    22,    23,    24,     0,
       0,     0,     0,     0,    25,     0,     0,    26,    27,    81,
       0,    82,    83,    84,    85,    86,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,    97,    98,    99,   100,
       2,     3,     0,     0,     0,     0,     0,     4,     5,     6,
       0,     0,     0,     0,     7,     8,     0,     0,     0,     9,
      10,    11,   174,     0,     0,    14,     0,     0,     0,     0,
       0,   101,   102,   103,   104,    47,     0,   175,     2,     3,
       0,     0,     0,     0,     0,     4,     5,     6,     0,   176,
       0,     0,     7,     8,     0,     0,   181,     9,    10,    11,
      12,     2,     3,    14,     0,     0,     0,     0,     4,     5,
       6,     0,     0,    47,     0,     7,     8,     0,     0,     0,
       9,    10,    11,   174,     0,     0,    14,   176,     0,     0,
       0,     0,     2,     3,     0,     0,    47,     0,   261,     4,
       5,     6,     0,     0,     0,     0,     7,     8,     0,     0,
     176,     9,    10,    11,    12,     2,     3,    14,     0,     0,
       0,     0,     4,     5,     6,     0,     0,    47,     0,     7,
       8,     0,     0,     0,     9,    10,    11,   186,     0,     0,
      14,   176,     0,     0,     0,     0,     2,     3,     0,     0,
      47,     0,   187,     4,     5,     6,     0,     0,     0,     0,
       7,     8,     0,     0,     0,     9,    10,    11,   174,     0,
       0,    14,     0,     0,     0,     0,     0,     2,     3,     0,
       0,    47,     0,   135,     4,     5,     6,     0,     0,     0,
       0,     7,     8,     0,     0,     0,     9,    10,    11,    12,
       2,     3,    14,     0,     0,     0,     0,     4,     5,     6,
       0,     0,    47,     0,     7,     8,     0,     0,     0,     9,
      10,    11,    12,     0,     0,    14,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    22,    82,    83,    84,    85,
      86,    87,    88,    89,    90,    91,    92,    93,    94,    95,
      96,    97,    98,    99,   100,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   254,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   307,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
      94,    95,    96,    97,    98,    99,   100,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   321,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
      94,    95,    96,    97,    98,    99,   100,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   331,    82,    83,
      84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
      94,    95,    96,    97,    98,    99,   100,    84,    85,    86,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
      97,    98,    99,   100,    86,    87,    88,    89,    90,    91,
      92,    93,    94,    95,    96,    97,    98,    99,   100
};

static const short yycheck[] =
{
       7,    35,    30,    70,   195,     5,    23,   111,     5,   200,
       8,     5,    34,     5,    21,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,     5,     5,     5,     5,     5,    46,
       5,    38,    45,   110,     5,    66,    42,    75,    12,    12,
      45,     5,    59,     5,    42,    12,     5,    12,     5,    76,
       5,    57,   107,    42,   109,    66,    66,    55,    16,    76,
      77,    78,    66,    38,    66,     5,    69,    38,    42,    42,
       5,   272,     5,    81,    38,    42,    38,    42,     0,    38,
       5,    38,    37,    57,    57,    66,    66,    66,    66,    66,
      57,   129,    57,   101,   102,   103,   104,   124,    38,     5,
       5,    42,    37,     5,    37,     6,   123,   115,   222,     5,
     224,     5,    37,     6,    42,     5,     5,   207,    35,   196,
     128,     6,   199,     6,    66,    42,   216,    55,    46,    57,
      34,    37,    37,   210,    35,    37,    34,   214,    55,    34,
      57,    37,    35,    37,    55,   222,   184,    37,    37,    34,
      34,    34,   229,    35,    35,    35,    35,   234,    35,   236,
      42,    42,   239,    42,    42,    42,    34,    35,   195,    55,
     247,    34,    35,   200,    55,    42,    57,     5,    34,    35,
      42,   219,    34,    39,    40,    41,    42,    28,    29,    45,
      31,    32,    33,    34,    35,    36,   273,   205,   275,    55,
     217,   218,    20,    21,    22,    23,    24,    25,    26,   236,
      22,    23,    24,    25,    26,    17,    42,    38,    66,   296,
       6,     6,   299,    42,     6,   242,     6,   265,    55,     6,
     238,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,     6,   272,     4,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,    38,    66,    34,     5,
       5,   288,   289,     6,   318,   110,   111,    38,     6,     6,
       6,     3,   117,    67,   322,   323,   303,   325,     5,    34,
     298,    66,   206,   301,    37,   247,   140,   285,    20,    21,
      40,   118,    -1,   138,   139,    27,    28,    29,   128,    -1,
      -1,    -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,
      42,    43,    -1,    45,    -1,    47,    48,    -1,    50,    51,
      52,    53,    54,    55,    56,    57,    -1,    -1,    -1,    -1,
      -1,    63,    -1,    -1,    66,    67,    68,    -1,    -1,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    -1,    -1,
      -1,   206,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,     3,    -1,   220,   221,   222,    -1,   224,
      -1,    -1,   227,    -1,   229,    -1,    -1,    -1,    -1,    -1,
      20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,
      69,    -1,    -1,    -1,    34,    35,    -1,   252,   253,    39,
      40,    41,    42,    43,    -1,    45,    -1,    47,    48,    -1,
      50,    51,    52,    53,    54,    55,    56,    57,    -1,    -1,
      -1,    -1,    -1,    63,     3,    -1,    66,    67,    68,   284,
      -1,    -1,    -1,    -1,    -1,   290,   291,    -1,    -1,   294,
      -1,    20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,
      29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    43,   320,    45,    -1,    47,    48,
      -1,    50,    51,    52,    53,    54,    55,    56,    57,    -1,
      -1,    -1,    -1,    -1,    63,    -1,    -1,    66,    67,     6,
      -1,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      20,    21,    -1,    -1,    -1,    -1,    -1,    27,    28,    29,
      -1,    -1,    -1,    -1,    34,    35,    -1,    -1,    -1,    39,
      40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,    -1,
      -1,    58,    59,    60,    61,    55,    -1,    57,    20,    21,
      -1,    -1,    -1,    -1,    -1,    27,    28,    29,    -1,    69,
      -1,    -1,    34,    35,    -1,    -1,    38,    39,    40,    41,
      42,    20,    21,    45,    -1,    -1,    -1,    -1,    27,    28,
      29,    -1,    -1,    55,    -1,    34,    35,    -1,    -1,    -1,
      39,    40,    41,    42,    -1,    -1,    45,    69,    -1,    -1,
      -1,    -1,    20,    21,    -1,    -1,    55,    -1,    57,    27,
      28,    29,    -1,    -1,    -1,    -1,    34,    35,    -1,    -1,
      69,    39,    40,    41,    42,    20,    21,    45,    -1,    -1,
      -1,    -1,    27,    28,    29,    -1,    -1,    55,    -1,    34,
      35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,    -1,
      45,    69,    -1,    -1,    -1,    -1,    20,    21,    -1,    -1,
      55,    -1,    57,    27,    28,    29,    -1,    -1,    -1,    -1,
      34,    35,    -1,    -1,    -1,    39,    40,    41,    42,    -1,
      -1,    45,    -1,    -1,    -1,    -1,    -1,    20,    21,    -1,
      -1,    55,    -1,    57,    27,    28,    29,    -1,    -1,    -1,
      -1,    34,    35,    -1,    -1,    -1,    39,    40,    41,    42,
      20,    21,    45,    -1,    -1,    -1,    -1,    27,    28,    29,
      -1,    -1,    55,    -1,    34,    35,    -1,    -1,    -1,    39,
      40,    41,    42,    -1,    -1,    45,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    55,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    38,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    38,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    37,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    37,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     3,    20,    21,    27,    28,    29,    34,    35,    39,
      40,    41,    42,    43,    45,    47,    48,    50,    51,    52,
      53,    54,    55,    56,    57,    63,    66,    67,    71,    72,
      73,    79,    82,    85,    88,    92,    94,    95,    98,   101,
     102,   103,   104,   108,   109,   110,    34,    55,   101,   103,
     107,    42,    86,    87,    45,    45,    66,    66,   101,    34,
      35,    42,    55,    79,    16,    42,     0,    46,    98,    89,
      35,    42,    80,    84,    66,    73,    34,    34,     5,    66,
     109,     6,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    58,    59,    60,    61,    28,    29,    31,    32,    33,
      34,    35,    36,   101,    37,     5,    38,    34,     5,    66,
      66,   101,    55,     6,    34,    42,    55,    34,    35,    99,
     100,    42,    77,    42,    55,    57,    75,    76,     6,    35,
       5,    68,    93,    79,    96,   101,   101,   101,   103,   104,
     104,   104,   104,   104,   104,   104,   104,   104,   104,   104,
     104,   104,   104,   104,   104,   104,   104,   104,   103,   103,
     103,   103,   108,   108,    42,    57,    69,    76,   104,   105,
     106,    38,   105,    42,    37,   103,    42,    57,    78,   104,
      86,    37,     5,    38,   101,    35,    57,    74,    76,    79,
      34,    17,    42,   107,    98,     6,    34,     5,    66,    38,
       5,    38,   104,   104,    35,    80,    97,    66,    66,    37,
       6,    69,     5,    37,     5,    38,    98,     6,    42,     5,
      37,    55,    66,    74,    12,    76,     5,    37,     6,    12,
      76,    74,     6,    38,   103,    78,    77,    35,    42,    81,
      83,    76,     6,    35,    38,    77,   101,   101,    98,   104,
     104,    57,    76,   105,   105,     4,   104,    76,   104,    38,
      38,    76,    35,    57,    76,    79,    91,   103,    76,    37,
     101,    66,    37,    75,    34,     5,   104,   104,    66,    66,
      69,     6,    38,    98,     6,    74,    12,    76,     6,    12,
      76,     6,    90,     5,    38,   104,    81,    38,   101,   101,
     104,   104,   104,    38,    76,   103,    76,   103,    67,   101,
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
  unsigned int yylineno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylineno);
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
                        int NbPtr = ShowAlloc("init execution "); // number of un delele ptr
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
                        NbPtr = ShowAlloc("end execution -- ") - NbPtr;
                        
                        if (NbPtr) { cout << " ######## We forget of deleting   " << NbPtr << " Nb pointer  " << endl;}
  return 0;;}
    break;

  case 4:
#line 236 "lg.y"
    {yyval.cinst=yyvsp[0].cexp;;;;}
    break;

  case 5:
#line 237 "lg.y"
    { yyval.cinst= (yyvsp[-1].cinst+=yyvsp[0].cexp) ;}
    break;

  case 6:
#line 240 "lg.y"
    { yyval.clist_id=new ListOfId();;}
    break;

  case 7:
#line 241 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str));}
    break;

  case 8:
#line 242 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;}
    break;

  case 9:
#line 243 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>()));}
    break;

  case 10:
#line 244 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true));}
    break;

  case 11:
#line 245 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;}
    break;

  case 12:
#line 246 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;}
    break;

  case 13:
#line 247 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;}
    break;

  case 14:
#line 248 "lg.y"
    { yyval.clist_id = yyvsp[-2].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str)) ;}
    break;

  case 15:
#line 249 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-1].clist_id)) ;}
    break;

  case 16:
#line 250 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[-2].str,yyvsp[0].cexp)) ;}
    break;

  case 17:
#line 251 "lg.y"
    { yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-1].str),atype<FE<double> **>())) ;}
    break;

  case 18:
#line 252 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,Find(yyvsp[-2].str),atype<FE<double> **>(),true)) ;}
    break;

  case 19:
#line 253 "lg.y"
    { yyval.clist_id = yyvsp[-3].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-1].type->right())) ;}
    break;

  case 20:
#line 254 "lg.y"
    { yyval.clist_id = yyvsp[-4].clist_id; yyval.clist_id->push_back(UnId(yyvsp[0].str,C_F0(),yyvsp[-2].type,true)) ;}
    break;

  case 21:
#line 257 "lg.y"
    { yyval.clist_id = new ListOfId(); yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;}
    break;

  case 22:
#line 258 "lg.y"
    { yyval.clist_id=yyvsp[-2].clist_id  ; yyval.clist_id->push_back(UnId(yyvsp[0].str)); ;}
    break;

  case 25:
#line 263 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[0].str,dcltype);}
    break;

  case 26:
#line 264 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-2].str,dcltype,yyvsp[0].cexp);}
    break;

  case 27:
#line 265 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariable>(yyvsp[-3].str,dcltype,yyvsp[-1].args);}
    break;

  case 28:
#line 266 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 29:
#line 273 "lg.y"
    {yyval.args=yyvsp[0].cexp;}
    break;

  case 30:
#line 274 "lg.y"
    {yyval.args=Find(yyvsp[-1].str);}
    break;

  case 31:
#line 275 "lg.y"
    { yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);}
    break;

  case 32:
#line 276 "lg.y"
    { yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;}
    break;

  case 33:
#line 277 "lg.y"
    { yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp));}
    break;

  case 35:
#line 281 "lg.y"
    {yyval.type=TypeArray(yyvsp[-3].type,yyvsp[-1].type);}
    break;

  case 36:
#line 282 "lg.y"
    {yyval.type=TypeArray(yyvsp[-5].type,yyvsp[-3].type,yyvsp[-1].type);}
    break;

  case 37:
#line 289 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[0].str,currentblock,fespacetype,fespacecomplex); ;}
    break;

  case 38:
#line 290 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;}
    break;

  case 39:
#line 291 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[-2].str,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;}
    break;

  case 40:
#line 292 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[-1].clist_id,currentblock,fespacetype,fespacecomplex) ;}
    break;

  case 41:
#line 293 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;}
    break;

  case 42:
#line 294 "lg.y"
    { yyval.cexp =  NewFEvariable(yyvsp[-3].clist_id,currentblock,fespacetype,yyvsp[0].cexp,fespacecomplex) ;}
    break;

  case 43:
#line 297 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-3].str,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex); ;}
    break;

  case 44:
#line 298 "lg.y"
    { yyval.cexp =  NewFEarray(yyvsp[-4].clist_id,currentblock,fespacetype,yyvsp[-1].cexp,fespacecomplex) ;}
    break;

  case 45:
#line 302 "lg.y"
    {fespacecomplex=false;  fespacetype = Find(yyvsp[0].str);;}
    break;

  case 46:
#line 303 "lg.y"
    {
             if (yyvsp[-1].type != typevarreal && yyvsp[-1].type != typevarcomplex) yyerror(" type of finite element <real> or <complex>");
             fespacecomplex=(yyvsp[-1].type==typevarcomplex);
             fespacetype = Find(yyvsp[-3].str);;}
    break;

  case 47:
#line 308 "lg.y"
    {  yyval.cexp = yyvsp[0].cexp  ;}
    break;

  case 48:
#line 309 "lg.y"
    { yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;}
    break;

  case 49:
#line 311 "lg.y"
    {  yyval.cexp = yyvsp[0].cexp  ;}
    break;

  case 50:
#line 312 "lg.y"
    { yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);;}
    break;

  case 51:
#line 314 "lg.y"
    { yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;}
    break;

  case 52:
#line 315 "lg.y"
    { yyval.cexp=0;  yyval.cexp = yyvsp[0].cexp;}
    break;

  case 53:
#line 320 "lg.y"
    {yyval.cexp=currentblock->NewVar<LocalVariableFES,size_t>(yyvsp[-3].str,atype<pfes*>(),yyvsp[-1].args,dimFESpaceImage(yyvsp[-1].args));}
    break;

  case 55:
#line 323 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 56:
#line 326 "lg.y"
    {dcltype=yyvsp[0].type;}
    break;

  case 57:
#line 326 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 58:
#line 327 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 59:
#line 328 "lg.y"
    { yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 60:
#line 329 "lg.y"
    {yyval.cexp=currentblock->NewID(yyvsp[-4].type,yyvsp[-3].str,yyvsp[-1].cexp);;}
    break;

  case 61:
#line 331 "lg.y"
    {   /* use the stack to store the prev return type*/
                      assert(kkembtype+1<nbembtype);
                      rettype[++kkembtype] = yyvsp[-4].type->right();
                      yyvsp[-1].routine=new Routine(yyvsp[-5].type,yyvsp[-4].type->right(),yyvsp[-3].str,yyvsp[-1].clist_id,currentblock);;}
    break;

  case 62:
#line 336 "lg.y"
    { currentblock=yyvsp[-5].routine->Set(yyvsp[-1].cinst);
                       currentblock->Add(yyvsp[-7].str,"(",yyvsp[-5].routine);
                       kkembtype--;
                       yyval.cexp=0 ;}
    break;

  case 63:
#line 341 "lg.y"
    {currentblock = new Block(currentblock); yyvsp[-4].type->SetArgs(yyvsp[-1].clist_id);;}
    break;

  case 64:
#line 343 "lg.y"
    {  yyval.cinst=currentblock->close(currentblock);
                         yyval.cexp=currentblock->NewID(yyvsp[-8].type,yyvsp[-7].str,yyvsp[-1].cexp,*yyvsp[-5].clist_id);}
    break;

  case 65:
#line 347 "lg.y"
    {  currentblock = new Block(currentblock);}
    break;

  case 66:
#line 348 "lg.y"
    {  yyval.cexp=currentblock->close(currentblock);}
    break;

  case 67:
#line 350 "lg.y"
    {inloopcount++;;}
    break;

  case 68:
#line 351 "lg.y"
    {inloopcount++;}
    break;

  case 69:
#line 354 "lg.y"
    {dcltype=yyvsp[0].type;currentblock = new Block(currentblock);}
    break;

  case 70:
#line 355 "lg.y"
    {yyval.cexp=yyvsp[0].cexp;}
    break;

  case 71:
#line 357 "lg.y"
    {yyval.cexp=0;;}
    break;

  case 72:
#line 358 "lg.y"
    {zzzfff->input(yyvsp[0].str);yyval.cexp= 0; ;}
    break;

  case 73:
#line 359 "lg.y"
    {load(yyvsp[0].str);yyval.cexp= 0; ;}
    break;

  case 74:
#line 360 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 75:
#line 361 "lg.y"
    {yyval.cexp=yyvsp[0].cexp;}
    break;

  case 76:
#line 362 "lg.y"
    {inloopcount--; yyval.cexp=For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 77:
#line 364 "lg.y"
    {inloopcount--; 
                yyval.cexp=C_F0(For(yyvsp[-6].cexp,yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp),currentblock->close(currentblock));}
    break;

  case 78:
#line 367 "lg.y"
    {inloopcount--;yyval.cexp=While(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 79:
#line 368 "lg.y"
    {yyval.cexp=FIf(yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 80:
#line 369 "lg.y"
    {yyval.cexp=FIf(yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 81:
#line 370 "lg.y"
    { 
                      yyval.cexp=C_F0(new E_block(yyvsp[-1].cinst,yyvsp[0].cexp),atype<void>()) ;}
    break;

  case 82:
#line 372 "lg.y"
    {
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-1].str,C_F0(TheOperators,"[border]",yyvsp[0].args));}
    break;

  case 83:
#line 374 "lg.y"
    {
                      yyval.cexp=0;currentblock->NewID(atype<const E_Border *>(),yyvsp[-4].str,C_F0(TheOperators,"[border]",yyvsp[-2].args));}
    break;

  case 84:
#line 377 "lg.y"
    {
                    if(inloopcount) 
                      yyval.cexp= C_F0(new E_throw(E_exception::e_break),atype<void>()); 
                    else lgerror("break not in loop") ;}
    break;

  case 85:
#line 381 "lg.y"
    { 
                    if(inloopcount)
                        yyval.cexp= C_F0(new E_throw(E_exception::e_continue),atype<void>()) ;
                    else lgerror("continue not in loop");}
    break;

  case 86:
#line 385 "lg.y"
    { 
                    if (kkembtype>=0)
                      yyval.cexp= C_F0(new E_throw(E_exception::e_return,rettype[kkembtype]->CastTo(yyvsp[-1].cexp)) ,atype<void>());
                     else lgerror(" return not in routine ") ;}
    break;

  case 87:
#line 393 "lg.y"
    { 
   currentblock = new Block(currentblock);
   yyval.args = currentblock->NewVar<LocalVariable>(yyvsp[-5].str,atype<double*>());
   yyval.args+= yyvsp[-3].cexp;
   yyval.args+= yyvsp[-1].cexp ;}
    break;

  case 88:
#line 399 "lg.y"
    {   
   yyval.args = (yyvsp[-1].args += yyvsp[0].cexp);
   currentblock->close(currentblock);}
    break;

  case 90:
#line 406 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);;}
    break;

  case 97:
#line 420 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 98:
#line 421 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"+=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 99:
#line 422 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"-=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 100:
#line 423 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"*=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 101:
#line 424 "lg.y"
    {yyval.cexp=C_F0(TheOperators,"/=",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 103:
#line 429 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 104:
#line 430 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 105:
#line 431 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 106:
#line 432 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 107:
#line 433 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 108:
#line 434 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 109:
#line 435 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 110:
#line 436 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 111:
#line 437 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 112:
#line 438 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 113:
#line 439 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 114:
#line 440 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 115:
#line 441 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 116:
#line 442 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 117:
#line 443 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 118:
#line 444 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 119:
#line 445 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 120:
#line 446 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 121:
#line 447 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 122:
#line 452 "lg.y"
    {yyval.cexp=yyvsp[0].cexp;}
    break;

  case 123:
#line 453 "lg.y"
    {yyval.cexp=C_F0(TheOperators,":");}
    break;

  case 124:
#line 454 "lg.y"
    {yyval.cexp=C_F0(TheOperators,":",yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 125:
#line 455 "lg.y"
    {yyval.cexp=C_F0(TheOperators,":",yyvsp[-4].cexp,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 126:
#line 458 "lg.y"
    {yyval.args=0;}
    break;

  case 127:
#line 459 "lg.y"
    {yyval.args=Find(yyvsp[0].str);}
    break;

  case 128:
#line 460 "lg.y"
    { yyval.args=make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp);}
    break;

  case 129:
#line 461 "lg.y"
    {yyval.args=yyvsp[0].cexp;}
    break;

  case 130:
#line 462 "lg.y"
    { yyval.args = (yyvsp[-2].args += Find(yyvsp[0].str)) ;}
    break;

  case 131:
#line 463 "lg.y"
    { yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;}
    break;

  case 132:
#line 464 "lg.y"
    { yyval.args= (yyvsp[-4].args+= make_pair<const char *,const C_F0>(yyvsp[-2].str,yyvsp[0].cexp)) ;}
    break;

  case 133:
#line 467 "lg.y"
    {yyval.args=yyvsp[0].cexp;}
    break;

  case 134:
#line 468 "lg.y"
    {yyval.args = (yyvsp[-2].args += yyvsp[0].cexp) ;}
    break;

  case 136:
#line 473 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[0].cexp);}
    break;

  case 138:
#line 477 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 139:
#line 478 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[-1].oper,yyvsp[-2].cexp,yyvsp[0].cexp);}
    break;

  case 140:
#line 479 "lg.y"
    {yyval.cexp=C_F0(TheOperators,yyvsp[0].oper,yyvsp[-1].cexp);}
    break;

  case 141:
#line 483 "lg.y"
    {yyval.cexp=Find(yyvsp[0].str);;}
    break;

  case 142:
#line 484 "lg.y"
    {yyval.cexp= CConstant(yyvsp[0].lnum);}
    break;

  case 143:
#line 485 "lg.y"
    {yyval.cexp= CConstant(yyvsp[0].dnum);}
    break;

  case 144:
#line 486 "lg.y"
    {yyval.cexp= CConstant(complex<double>(0,yyvsp[0].dnum));}
    break;

  case 145:
#line 487 "lg.y"
    {yyval.cexp= CConstant<const char *>(yyvsp[0].str);}
    break;

  case 146:
#line 488 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].args);;}
    break;

  case 147:
#line 489 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-3].cexp,yyvsp[-2].oper,yyvsp[-1].cexp);}
    break;

  case 148:
#line 490 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-5].cexp,yyvsp[-4].oper,yyvsp[-3].cexp,yyvsp[-1].cexp);}
    break;

  case 149:
#line 491 "lg.y"
    {yyval.cexp=C_F0(yyvsp[-2].cexp,"[]");}
    break;

  case 150:
#line 492 "lg.y"
    { yyval.cexp=C_F0(yyvsp[-2].cexp,yyvsp[0].str) ;;}
    break;

  case 151:
#line 493 "lg.y"
    {yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);}
    break;

  case 152:
#line 494 "lg.y"
    {yyval.cexp=C_F0(TheRightOperators,yyvsp[0].oper,yyvsp[-1].cexp);}
    break;

  case 153:
#line 495 "lg.y"
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

  case 154:
#line 504 "lg.y"
    {yyval.cexp=yyvsp[-1].cexp;}
    break;

  case 155:
#line 505 "lg.y"
    { yyval.cexp=C_F0(TheOperators,"[]",yyvsp[-1].args);}
    break;


    }

/* Line 999 of yacc.c.  */
#line 2344 "lg.tab.cpp"

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
	  char *yymsg;
	  int yyx, yycount;

	  yycount = 0;
	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  for (yyx = yyn < 0 ? -yyn : 0;
	       yyx < (int) (sizeof (yytname) / sizeof (char *)); yyx++)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      yysize += yystrlen (yytname[yyx]) + 15, yycount++;
	  yysize += yystrlen ("syntax error, unexpected ") + 1;
	  yysize += yystrlen (yytname[yytype]);
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yycount = 0;
		  for (yyx = yyn < 0 ? -yyn : 0;
		       yyx < (int) (sizeof (yytname) / sizeof (char *));
		       yyx++)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			const char *yyq = ! yycount ? ", expecting " : " or ";
			yyp = yystpcpy (yyp, yyq);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yycount++;
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


#line 510 "lg.y"
 


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
extern void ShowAlloc(const char *s);
extern void ShowNbAlloc(const char *s);
void init_lgfem() ;
void init_lgmesh() ;
void init_algo();
bool withrgraphique = false;
const char * StrVersionNumber();

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


 

