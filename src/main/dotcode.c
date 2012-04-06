/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995  Robert Gentleman and Ross Ihaka
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
#include "Mathlib.h"
#include <string.h>
#include <stdlib.h>

typedef int (*FUNC)();

	/* The following code loads in a compatibility module */
	/* written by Luke Tierney to support S version 4 on */
	/* Hewlett-Packard machines.  The relevant defines are */
	/* set up by autoconfigure */


#ifdef HAVE_DLFCN_H

#include <dlfcn.h>

#else

#ifdef HAVE_DL_H
#include "compat/hp_dlfcn.c"
#define HAVE_DLFCN_H
#endif

#endif

#ifdef HAVE_DLFCN_H

	/* This code adds dynamic loading of code to R.  It */
	/* replaces the original dotcode.c and contains code */
	/* from that file. This code should work on a system */
	/* that implements the Solaris dlopen interface */
	/* standard. Shared libraries are loaded by dlopen */
	/* and symbols are looked for by dlsym. Unloading of */
	/* code is not implemented. This code was developed */
	/* on ELF-Linux. */
	/* Modifications Copyright (C) 1996 Heiner Schwarte */
	/* Minor cleanup and bug fix (C) 1996 Ross Ihaka */


	/* The dlopen function returns a handle after successfully */
	/* loading a library. These handles are collected in a list. */

struct voidptrlist{
	void *ptr;
	struct voidptrlist *next;
};

static struct voidptrlist *list=NULL;

static void add_ptr(void *p)
{
	struct voidptrlist *tmp;
	if (p==NULL)
		return;
	tmp=(struct voidptrlist *)malloc(sizeof(struct voidptrlist));
	tmp->ptr=p;
	tmp->next=list;
	list=tmp;
}


	/* findDynProc checks whether one of the libraries */
	/* that have been loaded contains the symbol name and */
	/* returns a pointer to that symbol upon success. */


FUNC findDynProc(char const *name) {
	struct voidptrlist *tmp;
	FUNC fcnptr;
	char buf[MAXIDSIZE+1];

#ifdef HAVE_NO_SYMBOL_UNDERSCORE
	sprintf(buf, "%s", name);
#else
	sprintf(buf, "_%s", name);
#endif

	for (tmp = list; tmp != NULL; tmp = tmp->next) {
	    /* The following line is not legal ANSI C. */
	    /* It is only meant to be used in systems supporting */
	    /* the dlopen() interface, in which systems data and  */
	    /* function pointers _are_ the same size and _can_   */
	    /* be cast without loss of information.              */
	    fcnptr = (FUNC)dlsym(tmp->ptr, buf);
	    if (fcnptr != (FUNC)0) return fcnptr;
	}
	return (FUNC)0;
}


	/* do_dynload implements the R-Interface for the */
	/* loading of libraries */

#ifndef RTLD_LAZY
#define RTLD_LAZY 1
#endif

SEXP do_dynload(SEXP call, SEXP op, SEXP args, SEXP env)
{
	void *handle;
	checkArity(op,args);
	if (!isString(CAR(args)))
		errorcall(call, "character argument expected\n");
	handle = dlopen(CHAR(STRING(CAR(args))[0]), RTLD_LAZY);
	if(!handle) {
		errorcall(call, (char*)dlerror());
	}
	add_ptr(handle);
	return R_NilValue;
}

#else

SEXP do_dynload(SEXP call, SEXP op, SEXP args, SEXP env)
{
	error("no dyn.load support in this R version\n");
}

#endif




#define NIL -1

typedef struct {
	char *name;
	FUNC func;
} CFunTabEntry;

#include "ForeignDecl.h"

static CFunTabEntry CFunTab[] =
{
#include "ForeignTab.h"
	{NULL, NULL}
};
#undef F77_SUBROUTINE
#undef C_FUNCTION


	/* HashTable stores name - pointer pairs. Open addressing */
	/* with linear probing is used. Sometimes the hashtable */
	/* will be expanded and reorganized. The implementation */
	/* is entirely elementary.  Possible sizes of the */
	/* table are 2^p+1 where p is a positive integer. */

static CFunTabEntry *HashTable;
static int HASHSIZE;
static int NumberElem;


static int HashCode(char *symbol)
{
	unsigned int  code=0;
	char *p = symbol;

	while (*p)
		code = 8 * code + *p++;
	return code % HASHSIZE;
}


static void HashInstall(char *name, FUNC func)
{
	int key;
	NumberElem++;
	key = HashCode(name);
	while (HashTable[key].name != NULL)
		key = (key + 1) % HASHSIZE;
	HashTable[key].name = (char*)malloc(strlen(name)+1);
	strcpy(HashTable[key].name,name);
	HashTable[key].func = func;
}

static void HashExpand()
{
	int oldsize;
	int i;
	CFunTabEntry *OldTable;
	oldsize=HASHSIZE;
	OldTable=HashTable;
	HASHSIZE=2*HASHSIZE-1;
	NumberElem=0;
	HashTable = (CFunTabEntry *) malloc(HASHSIZE * sizeof(CFunTabEntry));
	for (i = 0; i < HASHSIZE; i++)
		HashTable[i].name = NULL;
	for(i=0;i<oldsize;i++) {
		if(OldTable[i].name!=NULL)
			HashInstall(OldTable[i].name,OldTable[i].func);
	}	
	for(i=0;i<oldsize;i++)
		free(OldTable[i].name);
	free(OldTable);
}

static FUNC HashLookup(char *symbol)
{
	int key;
	key = HashCode(symbol);
	while (HashTable[key].name != NULL) {
		if (strcmp(symbol, HashTable[key].name) == 0)
			return HashTable[key].func;
		else
			key = (key + 1) % HASHSIZE;
	}
	return (FUNC)0;
}


	/* Initialization of the hashed load table */

static SEXP NaokSymbol;
static SEXP DupSymbol;

void InitFunctionHashing()
{
	int n;
	int i, size = 3;
	NaokSymbol = install("NAOK");
	DupSymbol = install("DUP");
	n = sizeof(CFunTab)/sizeof(CFunTabEntry);
	while(size < n)
		size = 2*size-1;
	HASHSIZE = size;	
	NumberElem = 0;
	HashTable = (CFunTabEntry*) malloc(HASHSIZE * sizeof(CFunTabEntry));
	for (i = 0; i < HASHSIZE; i++)
		HashTable[i].name = NULL;
	for (i = 0; CFunTab[i].name; i++)
		HashInstall(CFunTab[i].name, CFunTab[i].func);
	HashExpand();
}

	/* Convert an R object to a non-moveable C object */
	/* and return a pointer to it.  This leaves pointers */
	/* for anything other than vectors and lists unaltered. */

static void *RObjToCPtr(SEXP s, int naok, int dup, int narg)
{
	int *iptr;
	double *rptr;
	char **cptr;
	complex *zptr;
	int i, l, n;

	switch(TYPEOF(s)) {
	case LGLSXP:
	case FACTSXP:
	case ORDSXP:
	case INTSXP:
		n = LENGTH(s);
		iptr = INTEGER(s);
		for(i=0 ; i<n ; i++) {
			if(!naok && iptr[i] == NA_INTEGER)
				error("NAs in foreign function call (arg %d)\n", narg);
		}
		if(dup) {
			iptr = (int*)R_alloc(n, sizeof(int));
			for(i=0 ; i<n ; i++)
				iptr[i] = INTEGER(s)[i];
		}
		return (void*)iptr;
		break;
	case REALSXP:
		n = LENGTH(s);
		rptr = REAL(s);
		for(i=0 ; i<n ; i++) {
			if(!naok && !FINITE(rptr[i]))
				error("NAs in foreign function call (arg %d)\n", narg);
		}
		if(dup) {
			rptr = (double*)R_alloc(n, sizeof(double));
			for(i=0 ; i<n ; i++)
				rptr[i] = REAL(s)[i];
		}
		return (void*)rptr;
		break;
#ifdef COMPLEX_DATA
	case CPLXSXP:
		n = LENGTH(s);
		zptr = COMPLEX(s);
		for(i=0 ; i<n ; i++) {
			if(!naok && (!FINITE(zptr[i].r || !FINITE(zptr[i].i))))
				error("NAs in foreign function call (arg %d)\n", narg);
		}
		if(dup) {
			zptr = (complex*)R_alloc(n, sizeof(complex));
			for(i=0 ; i<n ; i++)
				zptr[i] = COMPLEX(s)[i];
		}
		return (void*)zptr;
		break;
#endif
	case STRSXP:
		if(!dup)
			error("character variables must be duplicated in .C/.Fortran\n");
		n = LENGTH(s);
		cptr = (char**)R_alloc(n, sizeof(char*));
		for(i=0 ; i<n ; i++) {
			l = strlen(CHAR(STRING(s)[i]));
			cptr[i] = (char*)R_alloc(l+1, sizeof(char));
			strcpy(cptr[i], CHAR(STRING(s)[i]));
		}
		return (void*)cptr;
		break;
	case LISTSXP:
		if(!dup) return (void*)s;
		n = length(s);
		cptr = (char**)R_alloc(n, sizeof(char*));
		for(i=0 ; i<n ; i++) {
			cptr[i] = (char*)s;
			s = CDR(s);
		}
		return (void*)cptr;
	default:
		return (void*)s;
	}
}

static SEXP CPtrToRObj(void *p, int n, SEXPTYPE type)
{
	int *iptr;
	double *rptr;
	char **cptr;
	complex *zptr;
	SEXP *lptr;
	int i;
	SEXP s, t;

	switch(type) {
		case LGLSXP:
		case FACTSXP:
		case ORDSXP:
		case INTSXP:
			s = allocVector(type, n);
			iptr = (int*)p;
			for(i=0 ; i<n ; i++) {
				INTEGER(s)[i] = iptr[i];
			}
			break;
		case REALSXP:
			s = allocVector(type, n);
			rptr = (double*)p;
			for(i=0 ; i<n ; i++) {
				REAL(s)[i] = rptr[i];
			}
			break;
#ifdef COMPLEX_DATA
		case CPLXSXP:
			s = allocVector(type, n);
			zptr = (complex*)p;
			for(i=0 ; i<n ; i++) {
				COMPLEX(s)[i] = zptr[i];
			}
			break;
#endif
		case STRSXP:
			PROTECT(s = allocVector(type, n));
			cptr = (char**)p;
			for(i=0 ; i<n ; i++) {
				STRING(s)[i] = mkChar(cptr[i]);
			}
			UNPROTECT(1);
			break;
		case LISTSXP:
			PROTECT(t = s = allocList(n));
			lptr = (SEXP*)p;
			for(i=0 ; i<n ; i++) {
				CAR(t) = lptr[i];
				t = CDR(t);
			}
			UNPROTECT(1);
		default:
			s = (SEXP)p;
	}
	return s;
}

	/* Foreign Function Interface.  This code allows a */
	/* user to call C or Fortran code which is either */
	/* statically or dynamically linked into R. */


static SEXP naoktrim(SEXP s, int * len, int *naok, int *dup)
{
	SEXP value;

	if(s == R_NilValue) {
		value = R_NilValue;
		*naok = 0;
		*len = 0;
	}
	else if(TAG(s) == NaokSymbol) {
		value = naoktrim(CDR(s), len, naok, dup);
		*naok = asLogical(CAR(s));
	}
	else if(TAG(s) == DupSymbol) {
		value = naoktrim(CDR(s), len, naok, dup);
		*dup = asLogical(CAR(s));
	}
	else {
		CDR(s) = naoktrim(CDR(s), len, naok, dup);
		*len = *len + 1;
	}
	return s;
}

#define MAX_ARGS 65

SEXP do_dotCode(SEXP call, SEXP op, SEXP args, SEXP env)
{
	void **cargs;
	int dup, naok, nargs, which;
	FUNC fun;
	SEXP pargs, s;
	char buf[128], *p, *q, *vmax;
	
	vmax = vmaxget();
	which = PRIMVAL(op);

	op = CAR(args);
	if (!isString(op))
		errorcall(call, "function name must be a string\n");

		/* The following code modifies the argument list */
		/* We know this is ok because do_dotcode is entered */
		/* with its arguments evaluated. */
	
	args = naoktrim(CDR(args), &nargs, &naok, &dup);
	if(naok == NA_LOGICAL)
		errorcall(call, "invalid naok value\n");

	if(nargs > MAX_ARGS)
		errorcall(call, "too many arguments in foreign function call\n");
	cargs = (void**)R_alloc(nargs, sizeof(void*));
	
		/* Convert the arguments for use in foreign */
		/* function calls.  Note that we copy twice */
		/* once here, on the way into the call, and */
		/* once below on the way out. */

	nargs = 0;
	for(pargs = args ; pargs != R_NilValue; pargs = CDR(pargs)) {
		cargs[nargs] = RObjToCPtr(CAR(pargs), naok, dup, nargs+1);
		nargs++;
	}
	
	/* make up load symbol & look it up */
	
	p = CHAR(STRING(op)[0]);
	q = buf;
	while ((*q = *p) != '\0') {
		p++;
		q++;
	}
#ifdef HAVE_F77_UNDERSCORE
	if (which)
		*q++ = '_';
	*q = '\0';
#endif
	if (!(fun = HashLookup(buf))){
#ifdef HAVE_DLFCN_H
		if(!(fun = findDynProc(buf))) {
			errorcall(call, "C/Fortran function not in load table\n");
		}
		else {
			if((1.0*NumberElem)/HASHSIZE > 0.5)
				HashExpand();
			HashInstall(buf,fun);
		}
#else
		errorcall(call, "C/Fortran function not in load table\n");
#endif
	}
	
	switch (nargs) {
	case 0:
		/* Silicon graphics C chokes if there is */
		/* no argument to fun */
		fun(0);
		break;
	case 1:
		fun(cargs[0]);
		break;
	case 2:
		fun(cargs[0], cargs[1]);
		break;
	case 3:
		fun(cargs[0], cargs[1], cargs[2]);
		break;
	case 4:
		fun(cargs[0], cargs[1], cargs[2], cargs[3]);
		break;
	case 5:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4]);
		break;
	case 6:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5]);
		break;
	case 7:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6]);
		break;
	case 8:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7]);
		break;
	case 9:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8]);
		break;
	case 10:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9]);
		break;
	case 11:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10]);
		break;
	case 12:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11]);
		break;
	case 13:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12]);
		break;
	case 14:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13]);
		break;
	case 15:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14]);
		break;
	case 16:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15]);
		break;
	case 17:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16]);
		break;
	case 18:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17]);
		break;
	case 19:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18]);
		break;
	case 20:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19]);
		break;
	case 21:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20]);
		break;
	case 22:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21]);
		break;
	case 23:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22]);
		break;
	case 24:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23]);
		break;
	case 25:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24]);
		break;
	case 26:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25]);
		break;
	case 27:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26]);
		break;
	case 28:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27]);
		break;
	case 29:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28]);
		break;
	case 30:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29]);
		break;
	case 31:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30]);
		break;
	case 32:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31]);
		break;
	case 33:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32]);
		break;
	case 34:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33]);
		break;
	case 35:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34]);
		break;
	case 36:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35]);
		break;
	case 37:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36]);
		break;
	case 38:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37]);
		break;
	case 39:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38]);
		break;
	case 40:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39]);
		break;
	case 41:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40]);
		break;
	case 42:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41]);
		break;
	case 43:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42]);
		break;
	case 44:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43]);
		break;
	case 45:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44]);
		break;
	case 46:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45]);
		break;
	case 47:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46]);
		break;
	case 48:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47]);
		break;
	case 49:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48]);
		break;
	case 50:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49]);
		break;
	case 51:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50]);
		break;
	case 52:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51]);
		break;
	case 53:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52]);
		break;
	case 54:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53]);
		break;
	case 55:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54]);
		break;
	case 56:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55]);
		break;
	case 57:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55], cargs[56]);
		break;
	case 58:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55], cargs[56], cargs[57]);
		break;
	case 59:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55], cargs[56], cargs[57], cargs[58]);
		break;
	case 60:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55], cargs[56], cargs[57], cargs[58], cargs[59]);
		break;
	case 61:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55], cargs[56], cargs[57], cargs[58], cargs[59],
		    cargs[60]);
		break;
	case 62:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55], cargs[56], cargs[57], cargs[58], cargs[59],
		    cargs[60], cargs[61]);
		break;
	case 63:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55], cargs[56], cargs[57], cargs[58], cargs[59],
		    cargs[60], cargs[61], cargs[62]);
		break;
	case 64:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55], cargs[56], cargs[57], cargs[58], cargs[59],
		    cargs[60], cargs[61], cargs[62], cargs[63]);
		break;
	case 65:
		fun(cargs[0],  cargs[1],  cargs[2],  cargs[3],  cargs[4],
		    cargs[5],  cargs[6],  cargs[7],  cargs[8],  cargs[9],
		    cargs[10], cargs[11], cargs[12], cargs[13], cargs[14],
		    cargs[15], cargs[16], cargs[17], cargs[18], cargs[19],
		    cargs[20], cargs[21], cargs[22], cargs[23], cargs[24],
		    cargs[25], cargs[26], cargs[27], cargs[28], cargs[29],
		    cargs[30], cargs[31], cargs[32], cargs[33], cargs[34],
		    cargs[35], cargs[36], cargs[37], cargs[38], cargs[39],
		    cargs[40], cargs[41], cargs[42], cargs[43], cargs[44],
		    cargs[45], cargs[46], cargs[47], cargs[48], cargs[49],
		    cargs[50], cargs[51], cargs[52], cargs[53], cargs[54],
		    cargs[55], cargs[56], cargs[57], cargs[58], cargs[59],
		    cargs[60], cargs[61], cargs[62], cargs[63], cargs[64]);
		break;
	default:
		errorcall(call, "too many arguments, sorry\n");
	}

	if(dup) {
		nargs = 0;
		for(pargs=args ; pargs != R_NilValue ; pargs=CDR(pargs)) {
			PROTECT(s = CPtrToRObj(cargs[nargs],
					LENGTH(CAR(pargs)),
					TYPEOF(CAR(pargs))));
			ATTRIB(s) = duplicate(ATTRIB(CAR(pargs)));
			CAR(pargs) = s;
			nargs++;
			UNPROTECT(1);
		}
	}
	vmaxset(vmax);
	return (args);
}

static struct {
	char *name;
	SEXPTYPE type;
} typeinfo[] = {
	{"logical",	LGLSXP},
	{"integer",	INTSXP},
	{"double",	REALSXP},
#ifdef COMPLEX_DATA
	{"complex",	CPLXSXP},
#endif
	{"character",	STRSXP},
	{"list",	LISTSXP},
	{NULL,		0}
};

static int string2type(char *s)
{
	int i;
	for(i=0 ; typeinfo[i].name ; i++) {
		if(!strcmp(typeinfo[i].name, s)) {
			return typeinfo[i].type;
		}
	}
	error("type \"%s\" not supported in interlanguage calls\n", s);
}

void call_R(char *func, long nargs, void **arguments, char **modes,
	long *lengths, char **names, long nres, char **results)
{
	SEXP call, pcall, s;
	SEXPTYPE type;
	int i, j, n;

	if(!isFunction((SEXP)func))
		error("invalid function in call_R\n");
	if(nargs < 0)
		error("invalid argument count in call_R\n");
	if(nres < 0)
		error("invalid return value count in call_R\n");
	PROTECT(pcall = call = allocList(nargs+1));
	TYPEOF(call) = LANGSXP;
	CAR(pcall) = (SEXP)func;

	for(i=0 ; i<nargs ; i++) {
		pcall = CDR(pcall);
		type = string2type(modes[i]);
		switch(type) {
		case LGLSXP:
		case INTSXP:
			CAR(pcall) = allocSExp(type);
			INTEGER(CAR(pcall)) = (int*)(arguments[i]);
			LENGTH(CAR(pcall)) = lengths[i];
			break;
		case REALSXP:
			CAR(pcall) = allocSExp(REALSXP);
			REAL(CAR(pcall)) = (double*)(arguments[i]);
			LENGTH(CAR(pcall)) = lengths[i];
			break;
#ifdef COMPLEX_DATA
		case CPLXSXP:
			CAR(pcall) = allocSExp(CPLXSXP);
			COMPLEX(CAR(pcall)) = (complex*)(arguments[i]);
			LENGTH(CAR(pcall)) = lengths[i];
			break;
#endif
		case STRSXP:
			n = lengths[i];
			CAR(pcall) = allocVector(STRSXP, n);
			for(j=0 ; j<n ; j++) {
				s = allocSExp(CHARSXP);
				CHAR(s) = (char*)(arguments[i]);
				LENGTH(s) = strlen(CHAR(s));
				STRING(CAR(pcall))[i] = s;
			}
			break;
		case LISTSXP:
			n = lengths[i];
			CAR(pcall) = allocList(n);
			s = CAR(pcall);
			for(j=0 ; j<n ; j++) {
				CAR(s) = (SEXP)(arguments[i]);
				s = CDR(s);
			}
			break;
		}
		if(names && names[i])
			TAG(pcall) = install(names[i]);
		NAMED(CAR(pcall)) = 2;
	}

	PROTECT(s = eval(call, R_GlobalEnv));

	switch(TYPEOF(s)) {
	case LGLSXP:
	case INTSXP:
	case REALSXP:
#ifdef COMPLEX_DATA
	case CPLXSXP:
#endif
	case STRSXP:
		if(nres > 0)
			results[0] = RObjToCPtr(s, 1, 1, 0);
		break;
	case LISTSXP:
		n = length(s);
		if(nres < n) n = nres;
		for(i=0 ; i<n ; i++) {
			results[i] = RObjToCPtr(s, 1, 1, 0);
			s = CDR(s);
		}
	}
	UNPROTECT(2);
	return;
}

void call_S(char *func, long nargs, void **arguments, char **modes,
	long *lengths, char **names, long nres, char **results)
{
	call_R(func, nargs, arguments, modes,
		lengths, names, nres, results);
}

