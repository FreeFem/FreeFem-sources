/* A Bison parser, made by GNU Bison 3.5.1.  */

/* Bison interface for Yacc-like parsers in C

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

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

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

#line 148 "lg.tab.hpp"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE lglval;

int lgparse (void);

#endif /* !YY_LG_LG_TAB_HPP_INCLUDED  */
