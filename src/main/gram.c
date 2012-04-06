extern char *malloc(), *realloc();

# line 2 "gram.y"
/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "Defn.h"

extern SEXP listAppend(SEXP,SEXP);
void pushCmt();
void popCmt();

static int eatln;
extern SEXP R_CommentSxp;

#define YYSTYPE		SEXP
#ifdef YYBYACC
#define YYRETURN(x)	{ return(x); }
#else
#define YYRETURN(x)	{ free((void*)yys); free((void*)yyv); return(x); }
#endif

# define STR_CONST 257
# define NUM_CONST 258
# define NULL_CONST 259
# define SYMBOL 260
# define FUNCTION 261
# define LEX_ERROR 262
# define LBB 263
# define ERROR 264
# define LEFT_ASSIGN 265
# define RIGHT_ASSIGN 266
# define FOR 267
# define IN 268
# define IF 269
# define ELSE 270
# define WHILE 271
# define NEXT 272
# define BREAK 273
# define REPEAT 274
# define GT 275
# define GE 276
# define LT 277
# define LE 278
# define EQ 279
# define NE 280
# define AND 281
# define OR 282
# define LOW 283
# define TILDE 284
# define UNOT 285
# define NOT 286
# define SPECIAL 287
# define UMINUS 288
# define UPLUS 289
#define yyclearin yychar = -1
#define yyerrok yyerrflag = 0
extern int yychar;
extern int yyerrflag;
#ifndef YYMAXDEPTH
#define YYMAXDEPTH 150
#endif
#ifndef YYSTYPE
#define YYSTYPE int
#endif
YYSTYPE yylval, yyval;
# define YYERRCODE 256

# line 166 "gram.y"


SEXP tagarg(SEXP arg, SEXP tag)
{
	switch (TYPEOF(tag)) {
	case NILSXP:
	case SYMSXP:
	case STRSXP:
		return lang2(arg, tag);
	default:
		error("incorrect tag type\n");
	}
}

/* Lists are created and grown using a special dotted pair. */
/* The CAR of the list points to the last cons-cell in the list */
/* and the CDR points to the first.  The list can be extracted */
/* from the pair by taking its CDR, while the CAR gives fast access */
/* to the end of the list. */

/* Create a stretchy-list dotted pair */

SEXP newlist(void)
{
	SEXP s = CONS(R_NilValue, R_NilValue);
	CAR(s) = s;
	return s;
}

/* Add a new element at the end of a stretchy list */

SEXP growlist(SEXP l, SEXP s)
{
	SEXP tmp;
	PROTECT(l);
	tmp = CONS(s, R_NilValue);
	UNPROTECT(1);
	SETCDR(CAR(l), tmp);
	CAR(l) = tmp;
	return l;
}


	/* Comment Handling */

	/* R_CommentSxp is of the same form as an expression list, */
	/* each time a new { is encountered a new element is placed */
	/* in the R_CommentSxp and when a } is encountered it is */
	/* removed. */

extern void ResetComment()
{
	R_CommentSxp = CONS(R_NilValue, R_NilValue);
}

void pushCmt()
{
	R_CommentSxp = CONS(R_NilValue, R_CommentSxp);
}

void popCmt()
{
	R_CommentSxp = CDR(R_CommentSxp);
}

int isComment(SEXP l)
{
	if (isList(l) && isString(CAR(l)) && !strncmp(CHAR(STRING(CAR(l))[0]), "#", 1))
		return 1;
	else
		return 0;
}

void addcomment(SEXP l)
{
	SEXP tcmt, cmt;
	int i, ncmt;

	tcmt = CAR(R_CommentSxp);

		/* Return if there are no comments */

	if (tcmt == R_NilValue || l == R_NilValue)
		return;

		/* Attach the comments as a comment attribute */

	ncmt = length(tcmt);
	cmt = allocVector(STRSXP, ncmt);
	for(i=0 ; i<ncmt ; i++) {
		STRING(cmt)[i] = CAR(tcmt);
		tcmt = CDR(tcmt);
	}
	PROTECT(cmt);
	setAttrib(l, R_CommentSymbol, cmt);
	UNPROTECT(1);
	
		/* Reset the comment accumulator */

	CAR(R_CommentSxp) = R_NilValue;
}

SEXP firstarg(SEXP s, SEXP tag)
{
	SEXP tmp;
	PROTECT(s);
	PROTECT(tag);
	tmp = newlist();
	tmp = growlist(tmp, s);
	TAG(CAR(tmp)) = tag;
	UNPROTECT(2);
	return tmp;
}

SEXP nextarg(SEXP l, SEXP s, SEXP tag)
{
	PROTECT(tag);
	l = growlist(l, s);
	TAG(CAR(l)) = tag;
	UNPROTECT(1);
	return l;
}


SEXP mkString(char *);
SEXP mkInteger(char *);
SEXP mkFloat(char *);
SEXP mkComplex(char *);
SEXP mkNA(void);
SEXP mkTrue(void);
SEXP mkFalse(void);


/*
//	Basic File IO:
//
//	This code is here because at this particular instant it
//	seems closely related to cget(), which appears below.
*/


int R_fgetc(FILE *fp)
{
	int c = fgetc(fp);
	return feof(fp) ? R_EOF : c;
}


static char *buf;		/* The input stream buffer */
static char *bufp;		/* Pointer within current buffer */
static int cnt = 0;		/* Pointer to character count */

static char buf1[MAXELTSIZE];	/* File or text buffer */

static char buf0[MAXELTSIZE];	/* Console buffer */
static char *bufp0;		/* Pointer within the console buffer */
static int cnt0 = 0;		/* Characters in the console buffer */

static current_input;

/*
//	Set the input stream for the parser
//	    input = 0	initialize
//	    input = 1	console
//	    input = 2	file
//	    input = 3	text
*/

void R_SetInput(int input)
{
	switch (input) {

	case 0:			/* Initialization / Reset */
		cnt = cnt0 = 0;
		buf = buf0;
		bufp = buf0;
		break;

	case 1:			/* Restore console values */
		cnt = cnt0;
		buf = buf0;
		bufp = bufp0;
		R_Console = 1;
		break;

	case 2:			/* Text or file input */
	case 3:
		if(R_Console == 1) {
			cnt0 = cnt;
			bufp0 = bufp;
		}
		cnt = 0;
		buf = buf1;
		bufp = buf1;
		R_Console = 0;
		break;
	}
	current_input = input;
}


/*
//	Fetch a single character from the current input stream.
//	The stream has been set by a call to R_SetInput().
*/

int cget()
{
	if (--cnt < 0) {
		switch(current_input) {

		case 1:
			if (ReadKBD(buf, MAXELTSIZE) == 0) {
				ClearerrConsole();
				return R_EOF;
			}
			break;

		case 2:
			if (fgets(buf, MAXELTSIZE, R_Inputfile) == NULL) {
				ResetConsole();
				return R_EOF;
			}
			break;

		case 3:
			if (R_ParseCnt < LENGTH(R_ParseText)) {
				strcpy(buf, CHAR(STRING(R_ParseText)[(R_ParseCnt)]));
				strcat(buf, "\n");
			}
			else return R_EOF;
			break;

		}
		R_ParseCnt++;
		bufp = buf;
		cnt = strlen(buf);
		cnt--;
	}
	return *bufp++;
}


/*
//	Push n characters back onto the input stream.
//	This is only called when the characters are
//	currently in the input buffer so pushing back
//	beyond the start of the buffer is impossible.
*/

void uncget(int n)
{
	cnt += n;
	bufp -= n;
}


/*
//	Lexical Analyzer:
//
//	Basic lexical analysis is performed by the following
//	routines.  Input is read a line at a time, and, if the
//	program is in batch mode, each input line is echoed to
//	standard output after it is read.
//
//	The function yylex() scans the input, breaking it into
//	tokens which are then passed to the parser.  The lexical
//	analyser maintains a symbol table (in a very messy fashion).
//
//	The fact that if statements need to parse differently
//	depending on whether the statement is being interpreted or
//	part of the body of a function causes the need for ifpop
//	and ifpush. When an if statement is encountered an 'i' is
//	pushed on a stack (provided there are parentheses active).
//	At later points this 'i' needs to be popped off of the if
//	stack.
*/

static int newline = 0;
static int reset = 1;
#ifndef DEBUG_LEX
static
#endif
char *parenp, parenstack[50];

static void ifpush(void)
{
	if (*parenp == '{' || *parenp == '[' || *parenp == '(' || *parenp == 'i')
		*++parenp = 'i';
}

static void ifpop(void)
{
	if (*parenp == 'i')
		parenp--;
}

static int typeofnext(void)
{
	int k, c;

	c = cget();
	if (isdigit(c)) k = 1;
	else if (isalpha(c) || c == '.')
		k = 2;
	else
		k = 3;
	uncget(1);
	return k;
}

static int nextchar(int expect)
{
	int c = cget();

	if (c == expect)
		return 1;
	else
		uncget(1);
	return 0;
}

		/* Special Symbols */
		/* Syntactic Keywords + Symbolic Constants */

struct {
	char *name;
	int token;
} keywords[] = {
	{ "NULL",	NULL_CONST	},
	{ "NA",		NUM_CONST	},
	{ "TRUE",	NUM_CONST	},
	{ "FALSE",	NUM_CONST	},
	{ "GLOBAL.ENV",	NUM_CONST	},
	{ "function",	FUNCTION	},
	{ "while",	WHILE		},
	{ "repeat",	REPEAT		},
	{ "for",	FOR		},
	{ "if",		IF		},
	{ "in",		IN		},
	{ "else",	ELSE		},
	{ "next",	NEXT		},
	{ "break",	BREAK		},
	{ "...",	SYMBOL		},
	{ 0,		0		}
};


	/* klookup has side effects, it sets yylval */

int klookup(s)
char *s;
{
	int i;

	for (i = 0; keywords[i].name; i++) {
		if (strcmp(keywords[i].name, s) == 0) {
			switch (keywords[i].token) {
			case NULL_CONST:
				PROTECT(yylval = R_NilValue);
				eatln = 0;
				break;
			case NUM_CONST:
				switch(i) {
				case 1:
					PROTECT(yylval = mkNA());
					break;
				case 2:
					PROTECT(yylval = mkTrue());
					break;
				case 3:
					PROTECT(yylval = mkFalse());
					break;
				case 4:
					PROTECT(yylval = R_GlobalEnv);
				}
				eatln = 0;
				break;
			case FUNCTION:
			case WHILE:
			case REPEAT:
			case FOR:
			case IF:
				eatln = 1;
				yylval = install(s);
				break;
			case IN:
				eatln = 1;
				break;
			case ELSE:
				ifpop();
				eatln = 1;
				break;
			case NEXT:
			case BREAK:
				eatln = 0;
				yylval = install(s);
				break;
			case SYMBOL:
				PROTECT(yylval = install(s));
				eatln = 0;
				break;
			}
			return keywords[i].token;
		}
	}
	return 0;
}

static void prompt()
{
	if (R_ParseText == R_NilValue && R_Console == 1)
		yyprompt(CHAR(STRING(GetOption(install("continue"), R_NilValue))[0]));
}

SEXP mkString(char *s)
{
	SEXP t;

	PROTECT(t = allocVector(STRSXP, 1));
	STRING(t)[0] = mkChar(s);
	UNPROTECT(1);
	return t;
}

SEXP mkFloat(char *s)
{
	SEXP t = allocVector(REALSXP, 1);
	REAL(t)[0] = atof(s);
	return t;
}

SEXP mkComplex(char *s)
{
	SEXP t = allocVector(CPLXSXP, 1);
	COMPLEX(t)[0].r = 0;
	COMPLEX(t)[0].i = atof(s);
	return t;
}

SEXP mkNA(void)
{
	SEXP t = allocVector(LGLSXP, 1);
	LOGICAL(t)[0] = NA_LOGICAL;
	return t;
}

SEXP mkTrue(void)
{
	SEXP s = allocVector(LGLSXP, 1);
	LOGICAL(s)[0] = 1;
	return s;
}

SEXP mkFalse(void)
{
	SEXP s = allocVector(LGLSXP, 1);
	LOGICAL(s)[0] = 0;
	return s;
}

void yyinit(void)
{
	newline = 0;
	reset = 1;
}

int yywrap()
{
	return feof(R_Inputfile);
}

#ifdef HAVE_LIBREADLINE
extern char R_prompt_buf[512];
#endif


void yyprompt(char *format, ...)
{
	va_list(ap);
	va_start(ap, format);
#ifdef HAVE_LIBREADLINE
	vsprintf(R_prompt_buf, format, ap);
#else
	REvprintf(format, ap);
#endif
	va_end(ap);
	fflush(stdout);
	RBusy(0);
}

void yyerror(char *s)
{
	int i;

	R_CommentSxp = R_NilValue;
	REprintf("%s", buf);
	for (i = 1; i < bufp - buf; i++) {
		REprintf(" ");
	}
	REprintf("^\n");
	if (R_Console == 0) {
		fclose(R_Inputfile);
		ResetConsole();
	}
	else {
		FlushConsole();
		REprintf("Error: %s\n", s);
	}
	newline = 0;
	reset = 1;
	cnt = 0;
}

void check_formals(SEXP formlist, SEXP new)
{
	int i;

	while( formlist != R_NilValue ) {
		if(TAG(formlist) == new ) {
			REprintf("%s", buf);
			for (i = 2; i < bufp - buf; i++) 
				REprintf(" ");
			REprintf("^\n");
			newline = 0;
			reset = 1;
			cnt = 0;
			error("Repeated formal argument.\n");
		}
		formlist=CDR(formlist);
	}
}

int yylex()
{
	SEXP f;
	int c, quote, kw;
	char *p, yytext[MAXELTSIZE];

	if (newline) {
		newline = 0;
		prompt();
	}

    again:
	if (reset) {
		parenp = parenstack;
		*parenp = ' ';
		reset = 0;
		eatln = 0;
		ResetComment();
	}

	while ((c = cget()) == ' ' || c == '\t' || c == '');

	if (c == '#') {
		p = yytext;
		*p++ = c;
		while ((c = cget()) != '\n' && c != R_EOF)
			*p++ = c;
		*p = '\0';
		if(R_CommentSxp != R_NilValue) {
			f = mkChar(yytext);
			f = CONS(f, R_NilValue);
			CAR(R_CommentSxp) = listAppend(CAR(R_CommentSxp), f);
		}
	}


	if (c == R_EOF) {
		return EOF;
	}

		/* This code deals with context sensitivity to      */
		/* newlines.  The main problem is finding out       */
		/* whether a newline is followed by an ELSE clause. */
		/* This is only of importance if we are inside one  */
		/* of "(", "[", or "{".			     */

	if (c == '\n') {
		if (eatln || *parenp == '[' || *parenp == '(') {
			prompt();
			goto again;
		}
		if (*parenp == 'i') {
			prompt();
			while ((c = cget()) == ' ' || c == '\t');
			if (c == R_EOF) {
				error("unexpected end-of-file in parse\n");
			}
			if (c == '#') {
				p = yytext;
				*p++ = c;
				while ((c = cget()) != '\n' && c != R_EOF)
					*p++ = c;
				*p = '\0';
				if(R_CommentSxp != R_NilValue) {
					f = mkChar(yytext);
					f = CONS(f, R_NilValue);
					CAR(R_CommentSxp) = listAppend(CAR(R_CommentSxp), f);
				}
			}
			if (c == '\n') {
				prompt();
				uncget(1);
				goto again;
			}
			if (c == '}') {
				while (*parenp == 'i')
					ifpop();
				parenp--;
				return c;
			}
			if (c == ',') {
				ifpop();
				return c;
			}
			uncget(1);
			if (!strncmp(bufp, "else", 4) && !isalnum(bufp[4]) && bufp[4] != '.') {
				eatln = 1;
				bufp += 4;
				cnt -= 4;
				ifpop();
				return ELSE;
			}
			ifpop();
		}
		else newline = 1;
		return '\n';
	}

		/* These are needed because both ";" and "," can */
		/* end an "if" clause without a newline.  Ifpop  */
		/* only does its thing in the right context.     */

	if (c == ';' || c == ',') {
		ifpop();
		return c;
	}

		/* Either digits or symbols can start with a "." */
		/* so we need to decide which it is and jump to  */
		/* the correct spot. */

	if (c == '.') {
		kw = typeofnext();
		if (kw >= 2) goto symbol;
	}

		/* literal numbers */

	if (c == '.' || isdigit(c)) {
		int seendot = (c == '.');
		int seenexp = 0;
		p = yytext;
		*p++ = c;
		while (isdigit(c = cget()) || c == '.' || c == 'e' || c == 'E') {
			if (c == 'E' || c == 'e') {
				if (seenexp)
					break;
				seenexp = 1;
				seendot = 1;
				*p++ = c;
				c = cget();
				if (!isdigit(c) && c != '+' && c != '-')
					break;
			}
			if (c == '.') {
				if (seendot)
					break;
				seendot = 1;
			}
			*p++ = c;
		}
		*p = '\0';
		if(c == 'i') {
			PROTECT(yylval = mkComplex(yytext));
		}
		else {
			PROTECT(yylval = mkFloat(yytext));
			uncget(1);
		}
		eatln = 0;
		return NUM_CONST;
	}

	/* literal strings */

	if (c == '\"' || c == '\'') {
		quote = c;
		p = yytext;
		while ((c = cget()) != R_EOF && c != quote) {
			if (c == '\n') {
				uncget(1);
				return ERROR;
			}
			if (c == '\\') {
				c = cget();
				switch (c) {
				case 'a':
					c = '\a';
					break;
				case 'b':
					c = '\b';
					break;
				case 'f':
					c = '\f';
					break;
				case 'n':
					c = '\n';
					break;
				case 'r':
					c = '\r';
					break;
				case 't':
					c = '\t';
					break;
				case 'v':
					c = '\v';
					break;
				case '\\':
					c = '\\';
					break;
				}
			}
			*p++ = c;
		}
		*p = '\0';
		PROTECT(yylval = mkString(yytext));
		eatln = 0;
		return STR_CONST;
	}

	/* special functions */
	if (c == '%') {
		p = yytext;
		*p++ = c;
		while ((c = cget()) != R_EOF && c != '%') {
			if (c == '\n') {
				uncget(1);
				return ERROR;
			}
			*p++ = c;
		}
		if (c == '%')
			*p++ = c;
		*p++ = '\0';
		PROTECT(yylval = install(yytext));
		eatln=1;
		return SPECIAL;
	}


	/* functions, constants and variables */

	/* gag, barf, but the punters want it */
	if (c == '_') {
		eatln = 1;
		yylval = install("<-");
		return LEFT_ASSIGN;
	}

    symbol:
	if (c == '.' || isalpha(c)) {
		p = yytext;
		do {
			*p++ = c;
		} while ((c = cget()) != R_EOF && (isalnum(c) || c == '.'));
		uncget(1);
		*p = '\0';

		if ((kw = klookup(yytext))) {
			if(kw == FUNCTION) pushCmt();
			return kw;
		}

		PROTECT(yylval = install(yytext));
		eatln = 0;
		return SYMBOL;
	}

	/* compound tokens */

	switch (c) {
	case '<':
		eatln = 1;
		if (nextchar('=')) {
			yylval = install("<=");
			return LE;
		}
		if (nextchar('-')) {
			yylval = install("<-");
			return LEFT_ASSIGN;
		}
		if (nextchar('<'))
			if (nextchar('-')) {
				yylval = install("<<-");
				return LEFT_ASSIGN;
			}
			else
				return ERROR;
		yylval = install("<");
		return LT;
	case '-':
		eatln = 1;
		if (nextchar('>'))
			if (nextchar('>')) {
				yylval = install("<<-");
				return RIGHT_ASSIGN;
			}
			else {
				yylval = install("<-");
				return RIGHT_ASSIGN;
			}
		yylval = install("-");
		return '-';
	case '>':
		eatln = 1;
		if (nextchar('=')) {
			yylval = install(">=");
			return GE;
		}
		yylval = install(">");
		return GT;
	case '!':
		eatln = 1;
		if (nextchar('=')) {
			yylval = install("!=");
			return NE;
		}
		yylval = install("!");
		return '!';
	case '=':
		eatln = 1;
		if (nextchar('=')) {
			yylval = install("==");
			return EQ;
		}
		return '=';
	case ':':
		eatln = 1;
		if (nextchar('=')) {
			yylval = install(":=");
			return LEFT_ASSIGN;
		}
		yylval = install(":");
		return ':';
	case '&':
		eatln = 1;
		if (nextchar('&')) {
			yylval = install("&&");
			return AND;
		}
		yylval = install("&");
		return AND;
	case '|':
		eatln = 1;
		if (nextchar('|')) {
			yylval = install("||");
			return OR;
		}
		yylval = install("|");
		return OR;
	case '{':
		*++parenp = c;
		yylval = install("{");
		pushCmt();
		return c;
	case '}':
		ifpop();
		if(*parenp == '{')
			popCmt();
		parenp--;
		return c;
	case '(':
		*++parenp = c;
		yylval = install("(");
		return c;
	case ')':
		eatln = 0;
		ifpop();
		parenp--;
		return c;
	case '[':
		*++parenp = c;
		if (nextchar('[')) {
			*++parenp = c;
			yylval = install("[[");
			return LBB;
		}
		yylval = install("[");
		return c;
	case ']':
		ifpop();
		eatln = 0;
		parenp--;
		return c;
	case '?':
		eatln = 1;
		strcpy(yytext, "help");
		yylval = install(yytext);
		return c;
	case '*':
		eatln=1;
		if (nextchar('*'))
			c='^';
		yytext[0] = c;
		yytext[1] = '\0';
		yylval = install(yytext);
		return c;
	case '+':
	case '/':
	case '^':
	case '~':
	case '$':
		eatln = 1;
		yytext[0] = c;
		yytext[1] = '\0';
		yylval = install(yytext);
		return c;
	default:
		return c;
	}
}
int yyexca[] ={
-1, 1,
	0, -1,
	-2, 0,
-1, 54,
	126, 0,
	-2, 15,
-1, 72,
	126, 0,
	-2, 25,
-1, 73,
	275, 0,
	276, 0,
	277, 0,
	278, 0,
	279, 0,
	280, 0,
	-2, 26,
-1, 74,
	275, 0,
	276, 0,
	277, 0,
	278, 0,
	279, 0,
	280, 0,
	-2, 27,
-1, 75,
	275, 0,
	276, 0,
	277, 0,
	278, 0,
	279, 0,
	280, 0,
	-2, 28,
-1, 76,
	275, 0,
	276, 0,
	277, 0,
	278, 0,
	279, 0,
	280, 0,
	-2, 29,
-1, 77,
	275, 0,
	276, 0,
	277, 0,
	278, 0,
	279, 0,
	280, 0,
	-2, 30,
-1, 78,
	275, 0,
	276, 0,
	277, 0,
	278, 0,
	279, 0,
	280, 0,
	-2, 31,
	};
# define YYNPROD 74
# define YYLAST 882
int yyact[]={

    23,    47,    32,   119,   127,    44,   134,    28,    26,   102,
    27,    92,    29,    98,    91,   125,    95,    84,   111,   106,
   110,    47,   133,    25,   116,    44,    47,    32,   109,   108,
    44,   107,    28,    26,   105,    27,    83,    29,   114,   121,
    62,   115,    60,    25,    58,    56,    47,    32,    25,    24,
    44,   120,    28,    26,    61,    27,    46,    29,    59,    30,
    57,    97,    48,     1,     0,    94,    47,    32,    25,     0,
    44,   118,    28,    26,     0,    27,    46,    29,     0,    30,
     0,    46,    89,    90,    30,     0,    47,    32,    25,     0,
    44,    33,    28,    26,     0,    27,    47,    29,     0,     0,
    44,    46,     0,     0,    30,     0,    47,    32,    25,     0,
    44,    96,    28,    26,     0,    27,    33,    29,     0,     0,
     0,    46,     0,     0,    30,     0,    47,    32,    25,     0,
    44,    93,    28,    26,   126,    27,    33,    29,     0,   131,
     0,    46,     0,     0,    30,     0,    47,    32,    25,     0,
    44,    46,    28,    26,    30,    27,    33,    29,     0,     0,
     0,    46,     0,     0,    30,     0,    47,    32,    25,     0,
    44,     0,    28,    26,     0,    27,    33,    29,     0,     0,
     0,    46,     0,     0,    30,     0,     0,     0,    25,     0,
     0,     0,     0,     0,     0,     0,    33,     0,     0,     0,
     0,    46,     0,    47,    30,     0,     0,    44,     0,     0,
     0,     0,     0,     0,     0,     0,    33,     0,     0,     0,
     0,    46,     0,     0,    30,    25,     0,     0,    45,     0,
    42,    43,     0,     0,     0,     0,    33,     0,     0,     0,
    39,    38,    34,    35,    36,    37,    40,    41,    45,     0,
     0,     0,    31,    45,     0,    42,    43,     0,    46,     0,
     0,    30,     0,     0,     0,    39,    38,    34,    35,    36,
    37,    40,    41,    45,     0,    42,    43,    31,     0,     0,
     0,     0,     0,     0,     0,    39,    38,    34,    35,    36,
    37,    40,    41,    45,     0,    42,    43,    31,     0,     0,
     0,     0,     0,     0,     0,    39,    38,    34,    35,    36,
    37,    40,    41,    45,     0,    42,    43,    31,     0,     0,
   117,     0,     0,    45,     0,    39,    38,    34,    35,    36,
    37,    40,    41,    45,     0,    42,    43,    31,     0,     0,
     0,     0,     0,     0,     0,    39,    38,    34,    35,    36,
    37,    40,    41,    45,     0,    42,    43,    31,     0,     0,
     0,     0,     0,     0,     0,    39,    38,    34,    35,    36,
    37,    40,    41,    45,     0,     0,     0,    31,     0,     0,
     0,     0,     0,     0,     0,    39,    38,    34,    35,    36,
    37,    40,    41,    45,     0,     0,     0,    31,     0,     0,
     0,     0,     0,     0,     0,    39,    38,    34,    35,    36,
    37,    40,    41,     0,    47,    32,     0,    31,    44,     0,
    28,    26,     0,    27,     0,    29,     0,     0,     0,     0,
    45,     0,     0,    47,    32,     0,    25,    44,     0,    28,
    26,     0,    27,     0,    29,     0,     0,     0,     2,     0,
     0,     0,    47,    32,    31,    25,    44,     0,    28,    26,
     0,    27,     0,    29,    47,    32,     0,     0,    44,    46,
    28,    13,    30,     0,    25,    29,     0,     0,    10,     0,
     0,    12,     0,    11,     0,     0,    25,     0,    46,     0,
     0,    30,     0,    13,     0,     0,     0,     0,     0,     0,
    10,    15,     0,    12,     0,    11,     0,    46,     0,     0,
    30,     0,     0,     0,    13,     0,     0,     0,     0,    46,
     0,    10,    30,    15,    12,     0,    11,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,    15,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     9,     0,     0,    14,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     9,     0,     0,    14,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     9,     0,     0,    14,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,    45,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,    39,    38,    34,    35,    36,    37,    40,
    45,     0,     0,     0,     0,    31,     0,     0,     0,     0,
     0,     0,    39,    38,    34,    35,    36,    37,     0,    45,
     0,     0,     0,     0,    31,     0,     0,     0,     0,     0,
     0,    45,     0,     0,     4,     6,     5,     7,     8,    16,
     0,     0,     0,    31,     0,    18,     0,    17,     0,    19,
    21,    22,    20,     0,     0,    31,     0,     6,     5,     7,
     8,    16,     0,     0,     0,     0,     0,    18,     0,    17,
     0,    19,    21,    22,    20,     0,     0,     0,    87,     5,
    88,    86,    16,     0,     0,     0,     0,    85,    18,     3,
    17,     0,    19,    21,    22,    20,     0,    49,    50,    51,
    52,    53,    54,    55,     0,     0,     0,     0,    63,     0,
     0,     0,     0,    64,    65,    66,    67,    68,    69,    70,
    71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
    81,    82,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,    99,   100,   101,     0,   103,
   104,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,   112,   113,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,   122,   123,   124,     0,     0,
     0,     0,     0,     0,   128,   129,     0,   130,     0,     0,
     0,     0,     0,     0,   132,     0,     0,     0,     0,     0,
     0,   135 };
int yypact[]={

 -1000,   438, -1000,   -10, -1000, -1000, -1000, -1000, -1000,   460,
   460,   460,   460,   460,   460,   460,     5,     4,     2,     0,
   460, -1000, -1000, -1000, -1000,   460,   460,   460,   460,   460,
   460,   460,   460,   460,   460,   460,   460,   460,   460,   460,
   460,   460,   460,   460,   481,   481,   481,  -246,     6,    90,
    70,    60,    60,   397,   130,    90,  -247,   460,   460,   460,
  -251,   460,   460,    90,    60,   428,   428,   167,   167,    60,
   -15,   167,   130,   416,   416,   416,   416,   416,   416,   397,
   378,    90,   110,    -7, -1000,    90,   -30,   -32,   -33,   -73,
   -75, -1000, -1000, -1000,   460,   460, -1000,    -3,   -37,    50,
    30,    90,  -265,    90,    10, -1000,    -5,   460,   460,   460,
   -78, -1000,    90,    90, -1000,  -256,   460,   460, -1000,   460,
 -1000,   481,    90,    90,    90, -1000,   460,   -39,    90,    90,
   -35, -1000,    90,   460, -1000,    90 };
int yypgo[]={

     0,    63,   747,    62,    61,    19,    36,    60,    58,    54,
    17 };
int yyr1[]={

     0,     1,     1,     1,     1,     1,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     9,
     7,     8,     3,     3,     3,     3,     3,     3,     6,     6,
    10,    10,    10,    10,    10,    10,    10,    10,     4,     4,
     4,     4,     4,     5 };
int yyr2[]={

     0,     1,     5,     7,     7,     5,     3,     3,     3,     3,
     7,     7,     5,     5,     5,     5,     5,     7,     7,     7,
     7,     7,     7,     7,     7,     7,     7,     7,     7,     7,
     7,     7,     7,     7,     7,     7,    13,     9,     7,    11,
     7,     7,     5,    11,     9,     7,     7,     3,     3,     7,
     7,    11,     1,     3,     7,     5,     7,     5,     3,     9,
     1,     3,     7,     5,     7,     5,     7,     5,     1,     3,
     7,     7,    11,     1 };
int yychk[]={

 -1000,    -1,    10,    -2,   256,   258,   257,   259,   260,   123,
    40,    45,    43,    33,   126,    63,   261,   269,   267,   271,
   274,   272,   273,    10,    59,    58,    43,    45,    42,    47,
    94,   287,    37,   126,   277,   278,   279,   280,   276,   275,
   281,   282,   265,   266,    40,   263,    91,    36,    -3,    -2,
    -2,    -2,    -2,    -2,    -2,    -2,    40,    -7,    40,    -8,
    40,    -9,    40,    -2,    -2,    -2,    -2,    -2,    -2,    -2,
    -2,    -2,    -2,    -2,    -2,    -2,    -2,    -2,    -2,    -2,
    -2,    -2,    -2,    -6,   -10,    -2,   260,   257,   259,    -6,
    -6,   260,   257,   125,    59,    10,    41,    -4,   260,    -2,
    -2,    -2,   260,    -2,    -2,    41,    -5,    61,    61,    61,
    93,    93,    -2,    -2,    41,    44,    61,   270,    41,   268,
    41,    44,    -2,    -2,    -2,    93,    -5,   260,    -2,    -2,
    -2,   -10,    -2,    61,    41,    -2 };
int yydef[]={

     1,    -2,     2,     0,     5,     6,     7,     8,     9,    52,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,    47,    48,     3,     4,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,    60,    60,    60,     0,     0,    53,
     0,    12,    13,    14,    -2,    16,    68,     0,     0,     0,
     0,     0,     0,    42,    17,    18,    19,    20,    21,    22,
    23,    24,    -2,    -2,    -2,    -2,    -2,    -2,    -2,    32,
    33,    34,    35,    73,    58,    61,     9,     7,     8,    73,
    73,    45,    46,    10,    55,    57,    11,     0,    69,    38,
     0,    40,     0,    41,     0,    37,     0,    63,    65,    67,
     0,    44,    54,    56,    73,     0,     0,     0,    50,     0,
    49,    60,    62,    64,    66,    43,     0,    71,    70,    39,
     0,    59,    36,     0,    51,    72 };
typedef struct { char *t_name; int t_val; } yytoktype;
#ifndef YYDEBUG
#	define YYDEBUG	0	/* don't allow debugging */
#endif

#if YYDEBUG

yytoktype yytoks[] =
{
	"STR_CONST",	257,
	"NUM_CONST",	258,
	"NULL_CONST",	259,
	"SYMBOL",	260,
	"FUNCTION",	261,
	"LEX_ERROR",	262,
	"LBB",	263,
	"ERROR",	264,
	"LEFT_ASSIGN",	265,
	"RIGHT_ASSIGN",	266,
	"FOR",	267,
	"IN",	268,
	"IF",	269,
	"ELSE",	270,
	"WHILE",	271,
	"NEXT",	272,
	"BREAK",	273,
	"REPEAT",	274,
	"GT",	275,
	"GE",	276,
	"LT",	277,
	"LE",	278,
	"EQ",	279,
	"NE",	280,
	"AND",	281,
	"OR",	282,
	"?",	63,
	"LOW",	283,
	"~",	126,
	"TILDE",	284,
	"UNOT",	285,
	"NOT",	286,
	"+",	43,
	"-",	45,
	"*",	42,
	"/",	47,
	"%",	37,
	"SPECIAL",	287,
	":",	58,
	"UMINUS",	288,
	"UPLUS",	289,
	"^",	94,
	"$",	36,
	"(",	40,
	"[",	91,
	"-unknown-",	-1	/* ends search */
};

char * yyreds[] =
{
	"-no such reduction-",
	"prog : /* empty */",
	"prog : prog '\n'",
	"prog : prog expr '\n'",
	"prog : prog expr ';'",
	"prog : prog error",
	"expr : NUM_CONST",
	"expr : STR_CONST",
	"expr : NULL_CONST",
	"expr : SYMBOL",
	"expr : '{' exprlist '}'",
	"expr : '(' expr ')'",
	"expr : '-' expr",
	"expr : '+' expr",
	"expr : '!' expr",
	"expr : '~' expr",
	"expr : '?' expr",
	"expr : expr ':' expr",
	"expr : expr '+' expr",
	"expr : expr '-' expr",
	"expr : expr '*' expr",
	"expr : expr '/' expr",
	"expr : expr '^' expr",
	"expr : expr SPECIAL expr",
	"expr : expr '%' expr",
	"expr : expr '~' expr",
	"expr : expr LT expr",
	"expr : expr LE expr",
	"expr : expr EQ expr",
	"expr : expr NE expr",
	"expr : expr GE expr",
	"expr : expr GT expr",
	"expr : expr AND expr",
	"expr : expr OR expr",
	"expr : expr LEFT_ASSIGN expr",
	"expr : expr RIGHT_ASSIGN expr",
	"expr : FUNCTION '(' formlist ')' gobble expr",
	"expr : expr '(' sublist ')'",
	"expr : IF ifcond expr",
	"expr : IF ifcond expr ELSE expr",
	"expr : FOR forcond expr",
	"expr : WHILE cond expr",
	"expr : REPEAT expr",
	"expr : expr LBB sublist ']' ']'",
	"expr : expr '[' sublist ']'",
	"expr : expr '$' SYMBOL",
	"expr : expr '$' STR_CONST",
	"expr : NEXT",
	"expr : BREAK",
	"cond : '(' expr ')'",
	"ifcond : '(' expr ')'",
	"forcond : '(' SYMBOL IN expr ')'",
	"exprlist : /* empty */",
	"exprlist : expr",
	"exprlist : exprlist ';' expr",
	"exprlist : exprlist ';'",
	"exprlist : exprlist '\n' expr",
	"exprlist : exprlist '\n'",
	"sublist : sub",
	"sublist : sublist gobble ',' sub",
	"sub : /* empty */",
	"sub : expr",
	"sub : SYMBOL '=' expr",
	"sub : SYMBOL '='",
	"sub : STR_CONST '=' expr",
	"sub : STR_CONST '='",
	"sub : NULL_CONST '=' expr",
	"sub : NULL_CONST '='",
	"formlist : /* empty */",
	"formlist : SYMBOL",
	"formlist : SYMBOL '=' expr",
	"formlist : formlist ',' SYMBOL",
	"formlist : formlist ',' SYMBOL '=' expr",
	"gobble : /* empty */",
};
#endif /* YYDEBUG */
#line 1 "/usr/lib/yaccpar"
/*	@(#)yaccpar 1.10 89/04/04 SMI; from S5R3 1.10	*/

/*
** Skeleton parser driver for yacc output
*/

/*
** yacc user known macros and defines
*/
#define YYERROR		goto yyerrlab
#define YYACCEPT	{ free(yys); free(yyv); return(0); }
#define YYABORT		{ free(yys); free(yyv); return(1); }
#define YYBACKUP( newtoken, newvalue )\
{\
	if ( yychar >= 0 || ( yyr2[ yytmp ] >> 1 ) != 1 )\
	{\
		yyerror( "syntax error - cannot backup" );\
		goto yyerrlab;\
	}\
	yychar = newtoken;\
	yystate = *yyps;\
	yylval = newvalue;\
	goto yynewstate;\
}
#define YYRECOVERING()	(!!yyerrflag)
#ifndef YYDEBUG
#	define YYDEBUG	1	/* make debugging available */
#endif

/*
** user known globals
*/
int yydebug;			/* set to 1 to get debugging */

/*
** driver internal defines
*/
#define YYFLAG		(-1000)

/*
** static variables used by the parser
*/
static YYSTYPE *yyv;			/* value stack */
static int *yys;			/* state stack */

static YYSTYPE *yypv;			/* top of value stack */
static int *yyps;			/* top of state stack */

static int yystate;			/* current state */
static int yytmp;			/* extra var (lasts between blocks) */

int yynerrs;			/* number of errors */

int yyerrflag;			/* error recovery flag */
int yychar;			/* current input token number */


/*
** yyparse - return 0 if worked, 1 if syntax error not recovered from
*/
int
yyparse()
{
	register YYSTYPE *yypvt;	/* top of value stack for $vars */
	unsigned yymaxdepth = YYMAXDEPTH;

	/*
	** Initialize externals - yyparse may be called more than once
	*/
	yyv = (YYSTYPE*)malloc(yymaxdepth*sizeof(YYSTYPE));
	yys = (int*)malloc(yymaxdepth*sizeof(int));
	if (!yyv || !yys)
	{
		yyerror( "out of memory" );
		return(1);
	}
	yypv = &yyv[-1];
	yyps = &yys[-1];
	yystate = 0;
	yytmp = 0;
	yynerrs = 0;
	yyerrflag = 0;
	yychar = -1;

	goto yystack;
	{
		register YYSTYPE *yy_pv;	/* top of value stack */
		register int *yy_ps;		/* top of state stack */
		register int yy_state;		/* current state */
		register int  yy_n;		/* internal state number info */

		/*
		** get globals into registers.
		** branch to here only if YYBACKUP was called.
		*/
	yynewstate:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;
		goto yy_newstate;

		/*
		** get globals into registers.
		** either we just started, or we just finished a reduction
		*/
	yystack:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;

		/*
		** top of for (;;) loop while no reductions done
		*/
	yy_stack:
		/*
		** put a state and value onto the stacks
		*/
#if YYDEBUG
		/*
		** if debugging, look up token value in list of value vs.
		** name pairs.  0 and negative (-1) are special values.
		** Note: linear search is used since time is not a real
		** consideration while debugging.
		*/
		if ( yydebug )
		{
			register int yy_i;

			(void)printf( "State %d, token ", yy_state );
			if ( yychar == 0 )
				(void)printf( "end-of-file\n" );
			else if ( yychar < 0 )
				(void)printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				(void)printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ++yy_ps >= &yys[ yymaxdepth ] )	/* room on stack? */
		{
			/*
			** reallocate and recover.  Note that pointers
			** have to be reset, or bad things will happen
			*/
			int yyps_index = (yy_ps - yys);
			int yypv_index = (yy_pv - yyv);
			int yypvt_index = (yypvt - yyv);
			yymaxdepth += YYMAXDEPTH;
			yyv = (YYSTYPE*)realloc((char*)yyv,
				yymaxdepth * sizeof(YYSTYPE));
			yys = (int*)realloc((char*)yys,
				yymaxdepth * sizeof(int));
			if (!yyv || !yys)
			{
				yyerror( "yacc stack overflow" );
				return(1);
			}
			yy_ps = yys + yyps_index;
			yy_pv = yyv + yypv_index;
			yypvt = yyv + yypvt_index;
		}
		*yy_ps = yy_state;
		*++yy_pv = yyval;

		/*
		** we have a new state - find out what to do
		*/
	yy_newstate:
		if ( ( yy_n = yypact[ yy_state ] ) <= YYFLAG )
			goto yydefault;		/* simple state */
#if YYDEBUG
		/*
		** if debugging, need to mark whether new token grabbed
		*/
		yytmp = yychar < 0;
#endif
		if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
			yychar = 0;		/* reached EOF */
#if YYDEBUG
		if ( yydebug && yytmp )
		{
			register int yy_i;

			(void)printf( "Received token " );
			if ( yychar == 0 )
				(void)printf( "end-of-file\n" );
			else if ( yychar < 0 )
				(void)printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				(void)printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ( ( yy_n += yychar ) < 0 ) || ( yy_n >= YYLAST ) )
			goto yydefault;
		if ( yychk[ yy_n = yyact[ yy_n ] ] == yychar )	/*valid shift*/
		{
			yychar = -1;
			yyval = yylval;
			yy_state = yy_n;
			if ( yyerrflag > 0 )
				yyerrflag--;
			goto yy_stack;
		}

	yydefault:
		if ( ( yy_n = yydef[ yy_state ] ) == -2 )
		{
#if YYDEBUG
			yytmp = yychar < 0;
#endif
			if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
				yychar = 0;		/* reached EOF */
#if YYDEBUG
			if ( yydebug && yytmp )
			{
				register int yy_i;

				(void)printf( "Received token " );
				if ( yychar == 0 )
					(void)printf( "end-of-file\n" );
				else if ( yychar < 0 )
					(void)printf( "-none-\n" );
				else
				{
					for ( yy_i = 0;
						yytoks[yy_i].t_val >= 0;
						yy_i++ )
					{
						if ( yytoks[yy_i].t_val
							== yychar )
						{
							break;
						}
					}
					(void)printf( "%s\n", yytoks[yy_i].t_name );
				}
			}
#endif /* YYDEBUG */
			/*
			** look through exception table
			*/
			{
				register int *yyxi = yyexca;

				while ( ( *yyxi != -1 ) ||
					( yyxi[1] != yy_state ) )
				{
					yyxi += 2;
				}
				while ( ( *(yyxi += 2) >= 0 ) &&
					( *yyxi != yychar ) )
					;
				if ( ( yy_n = yyxi[1] ) < 0 )
					YYACCEPT;
			}
		}

		/*
		** check for syntax error
		*/
		if ( yy_n == 0 )	/* have an error */
		{
			/* no worry about speed here! */
			switch ( yyerrflag )
			{
			case 0:		/* new error */
				yyerror( "syntax error" );
				goto skip_init;
			yyerrlab:
				/*
				** get globals into registers.
				** we have a user generated syntax type error
				*/
				yy_pv = yypv;
				yy_ps = yyps;
				yy_state = yystate;
				yynerrs++;
			skip_init:
			case 1:
			case 2:		/* incompletely recovered error */
					/* try again... */
				yyerrflag = 3;
				/*
				** find state where "error" is a legal
				** shift action
				*/
				while ( yy_ps >= yys )
				{
					yy_n = yypact[ *yy_ps ] + YYERRCODE;
					if ( yy_n >= 0 && yy_n < YYLAST &&
						yychk[yyact[yy_n]] == YYERRCODE)					{
						/*
						** simulate shift of "error"
						*/
						yy_state = yyact[ yy_n ];
						goto yy_stack;
					}
					/*
					** current state has no shift on
					** "error", pop stack
					*/
#if YYDEBUG
#	define _POP_ "Error recovery pops state %d, uncovers state %d\n"
					if ( yydebug )
						(void)printf( _POP_, *yy_ps,
							yy_ps[-1] );
#	undef _POP_
#endif
					yy_ps--;
					yy_pv--;
				}
				/*
				** there is no state on stack with "error" as
				** a valid shift.  give up.
				*/
				YYABORT;
			case 3:		/* no shift yet; eat a token */
#if YYDEBUG
				/*
				** if debugging, look up token in list of
				** pairs.  0 and negative shouldn't occur,
				** but since timing doesn't matter when
				** debugging, it doesn't hurt to leave the
				** tests here.
				*/
				if ( yydebug )
				{
					register int yy_i;

					(void)printf( "Error recovery discards " );
					if ( yychar == 0 )
						(void)printf( "token end-of-file\n" );
					else if ( yychar < 0 )
						(void)printf( "token -none-\n" );
					else
					{
						for ( yy_i = 0;
							yytoks[yy_i].t_val >= 0;
							yy_i++ )
						{
							if ( yytoks[yy_i].t_val
								== yychar )
							{
								break;
							}
						}
						(void)printf( "token %s\n",
							yytoks[yy_i].t_name );
					}
				}
#endif /* YYDEBUG */
				if ( yychar == 0 )	/* reached EOF. quit */
					YYABORT;
				yychar = -1;
				goto yy_newstate;
			}
		}/* end if ( yy_n == 0 ) */
		/*
		** reduction by production yy_n
		** put stack tops, etc. so things right after switch
		*/
#if YYDEBUG
		/*
		** if debugging, print the string that is the user's
		** specification of the reduction which is just about
		** to be done.
		*/
		if ( yydebug )
			(void)printf( "Reduce by (%d) \"%s\"\n",
				yy_n, yyreds[ yy_n ] );
#endif
		yytmp = yy_n;			/* value to switch over */
		yypvt = yy_pv;			/* $vars top of value stack */
		/*
		** Look in goto table for next state
		** Sorry about using yy_state here as temporary
		** register variable, but why not, if it works...
		** If yyr2[ yy_n ] doesn't have the low order bit
		** set, then there is no action to be done for
		** this reduction.  So, no saving & unsaving of
		** registers done.  The only difference between the
		** code just after the if and the body of the if is
		** the goto yy_stack in the body.  This way the test
		** can be made before the choice of what to do is needed.
		*/
		{
			/* length of production doubled with extra bit */
			register int yy_len = yyr2[ yy_n ];

			if ( !( yy_len & 01 ) )
			{
				yy_len >>= 1;
				yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
				yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
					*( yy_ps -= yy_len ) + 1;
				if ( yy_state >= YYLAST ||
					yychk[ yy_state =
					yyact[ yy_state ] ] != -yy_n )
				{
					yy_state = yyact[ yypgo[ yy_n ] ];
				}
				goto yy_stack;
			}
			yy_len >>= 1;
			yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
			yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
				*( yy_ps -= yy_len ) + 1;
			if ( yy_state >= YYLAST ||
				yychk[ yy_state = yyact[ yy_state ] ] != -yy_n )
			{
				yy_state = yyact[ yypgo[ yy_n ] ];
			}
		}
					/* save until reenter driver code */
		yystate = yy_state;
		yyps = yy_ps;
		yypv = yy_pv;
	}
	/*
	** code supplied by user is placed in this switch
	*/
	switch( yytmp )
	{
		
case 1:
# line 67 "gram.y"
{ newline = 0; } break;
case 2:
# line 68 "gram.y"
{ R_CurrentExpr = NULL; return 2; } break;
case 3:
# line 69 "gram.y"
{ R_CurrentExpr = yypvt[-1]; UNPROTECT(1); YYRETURN(3); } break;
case 4:
# line 70 "gram.y"
{ R_CurrentExpr = yypvt[-1]; UNPROTECT(1); YYRETURN(4); } break;
case 5:
# line 71 "gram.y"
{ YYABORT; } break;
case 6:
# line 74 "gram.y"
{ yyval = yypvt[-0]; } break;
case 7:
# line 75 "gram.y"
{ yyval = yypvt[-0]; } break;
case 8:
# line 76 "gram.y"
{ yyval = yypvt[-0]; } break;
case 9:
# line 77 "gram.y"
{ yyval = yypvt[-0]; } break;
case 10:
# line 78 "gram.y"
{ UNPROTECT(1); TYPEOF(yypvt[-1]) = LANGSXP; CAR(yypvt[-1]) = yypvt[-2]; yyval = yypvt[-1]; PROTECT(yyval); eatln = 0; } break;
case 11:
# line 79 "gram.y"
{ UNPROTECT(1); yyval = lang2(yypvt[-2], yypvt[-1]); PROTECT(yyval); } break;
case 12:
# line 80 "gram.y"
{ UNPROTECT(1); yyval = lang2(yypvt[-1], yypvt[-0]); PROTECT(yyval); } break;
case 13:
# line 81 "gram.y"
{ UNPROTECT(1); yyval = lang2(yypvt[-1], yypvt[-0]); PROTECT(yyval); } break;
case 14:
# line 82 "gram.y"
{ UNPROTECT(1); yyval = lang2(yypvt[-1], yypvt[-0]); PROTECT(yyval); } break;
case 15:
# line 83 "gram.y"
{ UNPROTECT(1); yyval = lang2(yypvt[-1], yypvt[-0]); PROTECT(yyval); } break;
case 16:
# line 84 "gram.y"
{ UNPROTECT(1); yyval = lang2(yypvt[-1], yypvt[-0]); PROTECT(yyval); } break;
case 17:
# line 85 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 18:
# line 86 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 19:
# line 87 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 20:
# line 88 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 21:
# line 89 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 22:
# line 90 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 23:
# line 91 "gram.y"
{ UNPROTECT(3); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 24:
# line 92 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 25:
# line 93 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 26:
# line 94 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 27:
# line 95 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 28:
# line 96 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 29:
# line 97 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 30:
# line 98 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 31:
# line 99 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 32:
# line 100 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 33:
# line 101 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 34:
# line 102 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 35:
# line 103 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-1], yypvt[-0], yypvt[-2]); PROTECT(yyval); } break;
case 36:
# line 105 "gram.y"
{ addcomment(yypvt[-0]); UNPROTECT(2); yyval = lang3(yypvt[-5], CDR(yypvt[-3]), yypvt[-0]); PROTECT(yyval); popCmt();} break;
case 37:
# line 106 "gram.y"
{ if(isString(yypvt[-3])) yypvt[-3]=install(CHAR(STRING(yypvt[-3])[0])); UNPROTECT(2); if(length(CDR(yypvt[-1])) == 1 && CADR(yypvt[-1]) == R_MissingArg )
										yyval = lang1(yypvt[-3]);
									else
										yyval = LCONS(yypvt[-3], CDR(yypvt[-1]));
									PROTECT(yyval); } break;
case 38:
# line 111 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-2], yypvt[-1], yypvt[-0]); PROTECT(yyval);  } break;
case 39:
# line 112 "gram.y"
{ UNPROTECT(3); yyval = lang4(yypvt[-4], yypvt[-3], yypvt[-2], yypvt[-0]); PROTECT(yyval); } break;
case 40:
# line 113 "gram.y"
{ UNPROTECT(2); yyval = lang4(yypvt[-2], CAR(yypvt[-1]), CDR(yypvt[-1]), yypvt[-0]); PROTECT(yyval); } break;
case 41:
# line 114 "gram.y"
{ UNPROTECT(2); yyval = lang3(yypvt[-2], yypvt[-1], yypvt[-0]); PROTECT(yyval); } break;
case 42:
# line 115 "gram.y"
{ UNPROTECT(1); yyval = lang2(yypvt[-1], yypvt[-0]); PROTECT(yyval); } break;
case 43:
# line 116 "gram.y"
{ UNPROTECT(2); yyval = LCONS(yypvt[-3], LCONS(yypvt[-4], CDR(yypvt[-2]))); PROTECT(yyval); } break;
case 44:
# line 117 "gram.y"
{ UNPROTECT(2); yyval = LCONS(yypvt[-2], LCONS(yypvt[-3], CDR(yypvt[-1]))); PROTECT(yyval); } break;
case 45:
# line 118 "gram.y"
{ yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); UNPROTECT(2); PROTECT(yyval); } break;
case 46:
# line 119 "gram.y"
{ yyval = lang3(yypvt[-1], yypvt[-2], yypvt[-0]); UNPROTECT(2); PROTECT(yyval); } break;
case 47:
# line 120 "gram.y"
{ yyval = lang1(yypvt[-0]); PROTECT(yyval); } break;
case 48:
# line 121 "gram.y"
{ yyval = lang1(yypvt[-0]); PROTECT(yyval); } break;
case 49:
# line 125 "gram.y"
{ yyval = yypvt[-1];  eatln = 1; } break;
case 50:
# line 128 "gram.y"
{ yyval = yypvt[-1]; ifpush(); eatln = 1; } break;
case 51:
# line 131 "gram.y"
{ UNPROTECT(2); yyval = LCONS(yypvt[-3],yypvt[-1]); PROTECT(yyval); eatln=1;} break;
case 52:
# line 135 "gram.y"
{ yyval = newlist(); PROTECT(yyval); } break;
case 53:
# line 136 "gram.y"
{ addcomment(yypvt[-0]); UNPROTECT(1); yyval = growlist(newlist(), yypvt[-0]); PROTECT(yyval);} break;
case 54:
# line 137 "gram.y"
{ addcomment(yypvt[-0]); UNPROTECT(2); yyval = growlist(yypvt[-2], yypvt[-0]); PROTECT(yyval);} break;
case 55:
# line 138 "gram.y"
{ yyval = yypvt[-1]; addcomment(CAR(yyval));} break;
case 56:
# line 139 "gram.y"
{ addcomment(yypvt[-0]); UNPROTECT(2); yyval = growlist(yypvt[-2], yypvt[-0]); PROTECT(yyval);} break;
case 57:
# line 140 "gram.y"
{ yyval = yypvt[-1];} break;
case 58:
# line 143 "gram.y"
{ UNPROTECT(1); yyval = firstarg(CAR(yypvt[-0]),CADR(yypvt[-0])); PROTECT(yyval); } break;
case 59:
# line 144 "gram.y"
{ UNPROTECT(2); yyval = nextarg(yypvt[-3], CAR(yypvt[-0]), CADR(yypvt[-0])); PROTECT(yyval); } break;
case 60:
# line 147 "gram.y"
{ yyval = lang2(R_MissingArg,R_NilValue); PROTECT(yyval); } break;
case 61:
# line 148 "gram.y"
{ UNPROTECT(1); yyval = tagarg(yypvt[-0], R_NilValue); PROTECT(yyval); } break;
case 62:
# line 149 "gram.y"
{ UNPROTECT(2); yyval = tagarg(yypvt[-0], yypvt[-2]); PROTECT(yyval); } break;
case 63:
# line 150 "gram.y"
{ UNPROTECT(1); yyval = tagarg(R_MissingArg, yypvt[-1]); PROTECT(yyval); } break;
case 64:
# line 151 "gram.y"
{ UNPROTECT(2); yyval = tagarg(yypvt[-0], yypvt[-2]); PROTECT(yyval); } break;
case 65:
# line 152 "gram.y"
{ UNPROTECT(1); yyval = tagarg(R_MissingArg, yypvt[-1]); PROTECT(yyval); } break;
case 66:
# line 153 "gram.y"
{ UNPROTECT(2); yyval = tagarg(yypvt[-0], install("NULL")); PROTECT(yyval); } break;
case 67:
# line 154 "gram.y"
{ UNPROTECT(1); yyval = tagarg(R_MissingArg, install("NULL")); PROTECT(yyval); } break;
case 68:
# line 157 "gram.y"
{ yyval = R_NilValue; PROTECT(yyval); } break;
case 69:
# line 158 "gram.y"
{ UNPROTECT(1); yyval = firstarg(R_MissingArg, yypvt[-0]); PROTECT(yyval); } break;
case 70:
# line 159 "gram.y"
{ UNPROTECT(2); yyval = firstarg(yypvt[-0], yypvt[-2]); PROTECT(yyval); } break;
case 71:
# line 160 "gram.y"
{ UNPROTECT(2); check_formals(yypvt[-2],yypvt[-0]); yyval = nextarg(yypvt[-2], R_MissingArg, yypvt[-0]); PROTECT(yyval); } break;
case 72:
# line 161 "gram.y"
{ UNPROTECT(3); check_formals(yypvt[-4],yypvt[-2]); yyval = nextarg(yypvt[-4], yypvt[-0], yypvt[-2]); PROTECT(yyval); } break;
case 73:
# line 164 "gram.y"
{eatln = 1;} break;
	}
	goto yystack;		/* reset registers in driver code */
}
