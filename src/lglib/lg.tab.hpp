typedef union{ 
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
#define	IF	257
#define	ELSE	258
#define	SET	259
#define	LTLT	260
#define	GTGT	261
#define	OR	262
#define	AND	263
#define	EQ	264
#define	NE	265
#define	LE	266
#define	GE	267
#define	DOTSTAR	268
#define	DOTSLASH	269
#define	UNARY	270
#define	PLUSPLUS	271
#define	MOINSMOINS	272
#define	LNUM	273
#define	DNUM	274
#define	CNUM	275
#define	ID	276
#define	FESPACEID	277
#define	IDPARAM	278
#define	STRING	279
#define	ENDOFFILE	280
#define	INCLUDE	281
#define	LOAD	282
#define	BIDON	283
#define	FOR	284
#define	WHILE	285
#define	BREAK	286
#define	CONTINUE	287
#define	RETURN	288
#define	TRY	289
#define	CATCH	290
#define	THROW	291
#define	TYPE	292
#define	FUNCTION	293
#define	FESPACE	294
#define	PLUSEQ	295
#define	MOINSEQ	296
#define	MULEQ	297
#define	DIVEQ	298
#define	ARROW	299
#define	BORDER	300
#define	CURVE	301
#define	SOLVE	302


extern YYSTYPE lglval;
