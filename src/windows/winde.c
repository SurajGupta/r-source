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

#include "winde.h"

static int PUPcol;
static POINT DEMousePos;
static BOOL  ShiftState, ControlState;
static BOOL  NotSP;
BOOL CALLBACK DEVar(HWND, UINT, WPARAM, LPARAM);
RECT DERect;
HWND RDEWnd;
HMENU RMenuDE, RMenuDEWin;
COLORREF RDEPCol;            /* the current color of the pen */

/* 
   The spreadsheet function returns a list of vectors. The types of these
   vectors can be specified by the user as can their names. It the names
   are specified they are set during initialization. The user can change
   these via a menu interface, they can also change the type.

   The vectors are created too long and if they need to be increased this
   is done by using the next higher power of 2. They start 100 long. To cut
   them to the correct length for return you need to know the largest row number
   that was assigned to. LEVELS (sxpinfo.gp) is used to keep track of this, 
   separately for each vector. Vectors are initialized to NA when they are
   created so that NA is returned for any cell that was not set by the user.
   So that coercion back and forth maintains values of ssNA_REAL and ssNA_STRING
   I have set ssNA_STRING to be coerceVector(ssNA_REAL), very weird but easy.

        At some point the text drawing, line drawing etc should simply call
        the appropriate graphics functions rather than have all this duplicate
        code.

 */


/* 
   ssNewVector is just an interface to allocVector but it lets us
   set the fields to NA. We need to have a special NA for reals and
   strings so that we can differentiate between uninitialized elements
   in the vectors and user supplied NA's; hence ssNA_REAL and ssNA_STRING
 */

SEXP ssNewVector(SEXPTYPE type, int vlen)
{
        SEXP tvec;
        int j;

        tvec = allocVector(type, vlen);
        for (j = 0; j < vlen; j++)
                if (type == REALSXP)
                        REAL(tvec)[j] = ssNA_REAL;
                else if (type == STRSXP)
                        STRING(tvec)[j] = STRING(ssNA_STRING)[0];
        LEVELS(tvec) = 0;
        return (tvec);
}

SEXP do_dataentry(SEXP call, SEXP op, SEXP args, SEXP rho)
{
        SEXP  tvec2, tvec, colmodes, indata;
        SEXPTYPE type;
        int i, j, len;
        RCNTXT cntxt;

        PROTECT(indata = CAR(args));
        PROTECT(colmodes = CADR(args));

        if (!isList(indata) || !isList(colmodes))
                errorcall(call, "invalid argument\n");

        /* initialize the constants */

        bufp = buf;
        ne = 0;
        currentexp = 0;
        nneg = 0;
        ndecimal = 0;
        clength = 0;
        ccol = 1;
        crow = 1;
        colmin = 1;
        rowmin = 1;
        ssNA_REAL = -NA_REAL;
        tvec = allocVector(REALSXP, 1);
        REAL(tvec)[0] = ssNA_REAL;
        PROTECT(ssNA_STRING = coerceVector(tvec, STRSXP));
        ShiftState=0;
        ControlState=0;

        NotSP=0;
        RDEPCol=RRGB(0,0,0);
        doneSpread=1;
/* **** */
/* for windows bwidth and hwidth should be zero */
        bwidth = 0;
        hwidth = 0;



        /* setup inputlist  */

        if (indata != R_NilValue) {
                PROTECT(inputlist = duplicate(indata));
                tvec2 = colmodes;
                for (tvec = inputlist; tvec != R_NilValue; tvec = CDR(tvec)) {
                        type = TYPEOF(CAR(tvec));
                        if (CAR(tvec2) != R_NilValue)
                                type = str2type(CHAR(STRING(CAR(colmodes))[0]));
                        if (type != STRSXP)
                                type = REALSXP;
                        if (CAR(tvec) == R_NilValue) {
                                if (type == NILSXP)
                                        type = REALSXP;
                                CAR(tvec) = ssNewVector(type, 100);
                                TAG(tvec) = install("var1");
                                LEVELS(CAR(tvec)) = 0;
                        }
                        else if (!isVector(CAR(tvec)))
                                errorcall(call, "invalid type for value \n");
                        else {
                                if (TYPEOF(CAR(tvec)) != type)
                                        CAR(tvec) = coerceVector(CAR(tvec), type);
                                LEVELS(CAR(tvec)) = LENGTH(CAR(tvec));
                        }
                        tvec2 = CDR(tvec2);
                }
        }
        else {
                errorcall(call, "invalid parameter \n");
        }


        /* start up the window, more initializing in here */
        if (!DEInit())
                errorcall(call, "wrong device\n");

        /* set up a context which will close the window if there is an error */
        begincontext(&cntxt, 3, R_NilValue, R_NilValue, R_NilValue, R_NilValue);
        cntxt.cend = &CloseDE;

        highlightrect();

        while (doneSpread)
                EventLoop();

        endcontext(&cntxt);
        CloseDE();

        /* drop out unused columns */
        i = 0;
        for (tvec = inputlist; tvec != R_NilValue; tvec = CDR(tvec))
                if (CAR(tvec) == R_NilValue) {
                        if (i == 0)
                                inputlist = CDR(inputlist);
                        else {
                                tvec2 = nthcdr(inputlist, (i - 1));
                                SETCDR(tvec2, CDR(tvec));
                        }
                }
                else
                        i++;

        for (tvec = inputlist; tvec != R_NilValue; tvec = CDR(tvec)) {
                len = LEVELS(CAR(tvec));
                if (LENGTH(CAR(tvec)) != len) {
                        tvec2 = ssNewVector(TYPEOF(CAR(tvec)), len);
                        PROTECT(tvec);
                        for (j = 0; j < len; j++)
                                if (TYPEOF(CAR(tvec)) == REALSXP)
                                        if (REAL(CAR(tvec))[j] != ssNA_REAL)
                                                REAL(tvec2)[j] = REAL(CAR(tvec))[j];
                                        else
                                                REAL(tvec2)[j] = NA_REAL;
                                else if (TYPEOF(CAR(tvec)) == STRSXP)
                                        if (!streql(CHAR(STRING(CAR(tvec))[j]), CHAR(STRING(ssNA_STRING)[0])))
                                                STRING(tvec2)[j] = STRING(CAR(tvec))[j];
                                        else
                                                STRING(tvec2)[j] = NA_STRING;
                                else
                                        error("spreadsheet: internal memory problem");
                        CAR(tvec) = tvec2;
                        UNPROTECT(1);
                }
        }

        UNPROTECT(4);
        return inputlist;
}


static void doDEKey(WPARAM wParam)
{
        switch(wParam) {
                case VK_RETURN:
                case VK_DOWN:
                        advancerect(DOWN);
                        break;
                case VK_LEFT:
                        advancerect(LEFT);
                        break;
                case VK_RIGHT: 
                case VK_TAB:
                        advancerect(RIGHT);
                        break;
                case VK_UP:
                        advancerect(UP);
                        break;
                case VK_BACK: 
                case VK_DELETE: 
                        if (clength > 0) {
                                buf[clength - 1] = ' ';
                                printstring(buf, clength, crow, ccol);
                                clength--;
                                bufp--;
                        }
                        else
                                SysBeep();
                        break;  
                case VK_HOME:
                        ccol=1;
                        crow=1;
                        jumpwin(1, 1);
                        break;
                case VK_SHIFT:
                        ShiftState=1;
                        break;
                case VK_CONTROL:
                        ControlState=1;
                        break;
        }
}


/* Window Drawing Routines */

void drawwindow()
{
        int i;

        /* if there is an active cell enter the data in it */
        closerect();

        /* now set up the window with the new dimensions */
        GetClientRect(RDEWnd, &DERect);
        clearwindow(); 
        setattribsfromwindow();

        nwide = (windowWidth - 2 * bwidth) / box_w;

        setlineattribs(1);

        for (i = 1; i <= nwide; i++)
                drawline(i * box_w, hwidth, i * box_w, windowHeight);
        nhigh = (windowHeight - 2 * bwidth - hwidth) / box_h;
        for (i = 1; i <= nhigh; i++)
                drawline(0, hwidth + i * box_h, windowWidth, hwidth + i * box_h);
        colmax = colmin + (nwide - 2);  /* so row 0 and col 0 are reserved for labels */
        rowmax = rowmin + (nhigh - 2);
        printlabs();
        if (inputlist != R_NilValue)
                for (i = colmin; i <= colmax; i++)
                        drawcol(i);
        highlightrect();
}

/* find_coords finds the coordinates of the upper left corner of the given square on the screen */

void find_coords(int row, int col, int *xcoord, int *ycoord)
{
        *xcoord = bwidth + box_w * col;
        *ycoord = bwidth + hwidth + box_h * row;
}

/* 
   draw the window with the top left box at column wcol and 
   row wrow 
 */

void jumpwin(int wcol, int wrow)
{
        if (wcol < 0 || wrow < 0) {
                SysBeep();
                return;
        }
        closerect();
        colmin = wcol;
        rowmin = wrow;
        drawwindow();
}



void advancerect(int which)
{

        /* if we are in the header, changing a name then only down is allowed */
        if (crow < 1 && which != DOWN) {
                SysBeep();
                return;
        }

        closerect();

        switch (which) {
        case UP:
                if (crow == 1)
                        if (rowmin == 1)
                                SysBeep();
                        else
                                jumppage(UP);
                else
                        crow--;
                break;
        case DOWN:
                if (crow >= (nhigh - 1))
                        jumppage(DOWN);
                else
                        crow++;
                break;
        case RIGHT:
                if (ccol >= (nwide - 1))
                        jumppage(RIGHT);
                else
                        ccol++;
                break;
        case LEFT:
                if (ccol == 1)
                        if (colmin == 1)
                                SysBeep();
                        else
                                jumppage(LEFT);
                else
                        ccol--;
                break;
        default:
                abort();
        }

        highlightrect();
}

void drawrow(int whichrow)
{
        int i, src_x, src_y, lenip;
        char rlab[15];
        SEXP tvec;

        find_coords(whichrow, 0, &src_x, &src_y);
        cleararea(src_x, src_y, windowWidth, box_h);
        setlineattribs(1);
        for (i = 0; i <= nwide; i++)
                drawrectangle(i * box_w, src_y, box_w, box_h);

        sprintf(rlab, "R %d", rowmin + whichrow - 1);
        printstring(rlab, strlen(rlab), whichrow, 0);

        lenip = length(inputlist);
        for (i = colmin; i <= colmax; i++) {
                if (i > lenip)
                        break;
                tvec = CAR(nthcdr(inputlist, i - 1));
                if (tvec != R_NilValue)
                        if (whichrow + rowmin - 1 <= LEVELS(tvec))
                                printelt(tvec, whichrow + rowmin - 2, whichrow, i - colmin + 1);
        }

}

/* 
   printelt: print the correct value from vector[vrow] into the
   spread sheet in row ssrow and col sscol
 */
void printelt(SEXP invec, int vrow, int ssrow, int sscol)
{
        char *pp;

        if (TYPEOF(invec) == REALSXP) {
                if (REAL(invec)[vrow] != ssNA_REAL) {
                        pp=EncodeElement(invec, vrow, 0);
                        printstring(pp, strlen(pp), ssrow, sscol);
                }
        }
        else if (TYPEOF(invec) == STRSXP) {
                if (!streql(CHAR(STRING(invec)[vrow]), CHAR(STRING(ssNA_STRING)[0]))) {
                        pp=EncodeElement(invec, vrow, 0);
                        printstring(pp, strlen(pp), ssrow, sscol);
                }
        }
        else
                error("spreadsheet: internal memory error\n");
}

void drawcol(int whichcol)
{
        int i, src_x, src_y, len;
        char clab[15];
        SEXP tmp;

        find_coords(0, whichcol, &src_x, &src_y);
        cleararea(src_x, src_y, box_w, windowHeight);
        setlineattribs(1);
        for (i = 0; i <= nhigh; i++)
                drawrectangle(src_x, i * box_h, box_w, box_h);

        /* now fill it in if it is active */

        if (length(inputlist) >= whichcol + colmin - 1) {
                tmp = nthcdr(inputlist, whichcol + colmin - 2);
                if (TAG(tmp) != R_NilValue)
                        printstring(CHAR(PRINTNAME(TAG(tmp))),
                         strlen(CHAR(PRINTNAME(TAG(tmp)))), 0, whichcol);
                else {
                        sprintf(clab, "var%d", whichcol + colmin - 1);
                        printstring(clab, strlen(clab), 0, whichcol);
                }
                if (CAR(tmp) != R_NilValue) {
                        len = (LEVELS(CAR(tmp)) > rowmax) ? rowmax : LEVELS(CAR(tmp));
                        for (i = (rowmin - 1); i < len; i++)
                                printelt(CAR(tmp), i, i - rowmin + 2, whichcol);
                }
        }
        else {
                sprintf(clab, "var%d", whichcol + colmin - 1);
                printstring(clab, strlen(clab), 0, whichcol);
        }
}

void jumppage(int dir)
{
        switch (dir) {
        case UP:
                rowmin--;
                rowmax--;
                copyarea(0, hwidth + box_h, 0, hwidth + 2 * box_h);
                drawrow(1);
                break;
        case DOWN:
                rowmin++;
                rowmax++;
                copyarea(0, hwidth + 2 * box_h, 0, hwidth + box_h);
                drawrow((nhigh - 1));
                if (2 * bwidth + box_h * nhigh + hwidth != windowHeight)
                        drawrow(nhigh);
                break;
        case LEFT:
                colmin--;
                colmax--;
                copyarea(box_w, hwidth, 2 * box_w, hwidth);
                drawcol(1);
                break;
        case RIGHT:
                colmin++;
                colmax++;
                copyarea(2 * box_w, hwidth, box_w, hwidth);
                drawcol((nwide - 1));
                if (2 * bwidth + nwide * box_w != windowWidth)
                        drawcol(nwide);
                break;
        }
}

/* draw a rectangle, used to highlight/downlight the current box */
/* windows draws the line outside the rectangle so we need to subtract lwd */
void printrect(int lwd)
{
        setlineattribs(lwd);
        lwd--;
        drawrectangle(ccol * box_w+lwd, crow * box_h+lwd, box_w-lwd, box_h-lwd);
}

void downlightrect()
{
        RDEPCol=RRGB(255,255,255);
        printrect(2);
        RDEPCol=RRGB(0,0,0);
        printrect(1);
}

void highlightrect()
{
        RDEPCol=RRGB(0,0,0);
        printrect(2);
}

/* 
        when a buttonpress event happens find the square that is being pointed to

 */

int findsquare()
{

        int xw, yw, xr, yr, wcol, wrow;

        closerect();
        querypointer(&xr, &yr, &xw, &yw);

        /* translate to box coordinates */

        wcol = (xw - bwidth) / box_w;
        wrow = (yw - bwidth - hwidth) / box_h;

        /* see if it is in the row labels */
        if (wcol == 0) {
                SysBeep();
                highlightrect();
                return 0;
        }

        /* next check to see if it is in the column labels */

        if (yw < hwidth + bwidth + box_h)
                if (xw > bwidth + box_w) {
                        if( ccol != wcol ) {   /* if it is a new col set the focus to it */
                                closerect();
                                ccol=wcol;
                                crow=1;
                                highlightrect();
                        }
                        popupmenu(wcol);
                }
                else {
                        highlightrect();
                        SysBeep();
                }
        else if (wcol != ccol || wrow != crow) {
                ccol = wcol;
                crow = wrow;
        }
        highlightrect();
        return 0;
}

static SEXP getccol()
{
        SEXP tmp, tmp2;
        int i, len, wcol, wrow;
        SEXPTYPE type;
        char cname[10];

        wcol = ccol + colmin - 1;
        wrow = crow + rowmin - 1;
        if (length(inputlist) < wcol)
                inputlist = listAppend(inputlist, allocList(wcol - length(inputlist)));
        tmp = nthcdr(inputlist, wcol - 1);
        if (CAR(tmp) == R_NilValue) {
                len = (wrow < 100) ? 100 : wrow;
                CAR(tmp) = ssNewVector(REALSXP, len);
                if (TAG(tmp) == R_NilValue) {
                        sprintf(cname, "var%d", wcol);
                        TAG(tmp) = install(cname);
                }
        }
        if (!isVector(CAR(tmp)))
                error("internal type error in spreadsheet\n");
        len = LENGTH(CAR(tmp));
        type = TYPEOF(CAR(tmp));
        if (len < wrow) {
                tmp2 = ssNewVector(type, 2 * len);
                for (i = 0; i < len; i++)
                        if (type == REALSXP)
                                REAL(tmp2)[i] = REAL(CAR(tmp))[i];
                        else if (type == STRSXP)
                                STRING(tmp2)[i] = STRING(CAR(tmp))[i];
                        else
                                error("internal type error in spreadsheet\n");
                LEVELS(tmp2) = LEVELS(CAR(tmp));
                CAR(tmp) = tmp2;
        }
        return (CAR(tmp));
}

/* 
   close up the entry to a square, put the value that has been entered
   into  the correct place and as the correct type
 */

void closerect()
{
        SEXP cvec, tvec;
        
        *bufp = '\0';

        /* first check to see if anything has been entered */
        if (clength != 0) {
                if (crow == 0) {        /* then we are entering a new column name */
                        if (length(inputlist) < ccol + colmin - 1)
                                inputlist = listAppend(inputlist, allocList((ccol - colmin - 1 + length(inputlist))));
                        tvec = nthcdr(inputlist, ccol + colmin - 2);
                        TAG(tvec) = install(buf);
                }
                else {
                        cvec = getccol();
                        if ((crow + rowmin - 1) > LEVELS(cvec))
                                LEVELS(cvec) = (crow + rowmin - 1);
                        if (TYPEOF(cvec) == STRSXP) {
                                tvec = allocString(strlen(buf));
                                strcpy(CHAR(tvec), buf);
                                STRING(cvec)[(rowmin + crow - 2)] = tvec;
                        }
                        else
                                REAL(cvec)[(rowmin + crow - 2)] = atof(buf);
                }
        }
        else if (crow == 0) {
                sprintf(buf, "var%d", ccol);
                printstring(buf, strlen(buf), 0, ccol - colmin + 1);
        }

        downlightrect();

        ndecimal = 0;
        nneg = 0;
        ne = 0;
        currentexp = 0;
        clength = 0;
        bufp = buf;
}

/*
   print a null terminated string, check to see if it is longer than the print area and print
        it, left adjusted if necessary; clear the area of previous text;
        Windoze
 */
void printstring(char *ibuf, int buflen, int row, int col)
{
        int len, x_pos, y_pos;

        find_coords(row, col, &x_pos, &y_pos);
        cleararea(col * box_w + text_offset, hwidth + row * box_h + text_offset,
                  box_w - 2 * text_offset, box_h - 2 * text_offset);
        len = nchars(ibuf, buflen);
        drawtext(x_pos + text_offset, y_pos + text_offset, ibuf, len);
}

int nchars(char *ibuf, int len)
{
        int i;

        for (i = len; i > 1; i--)
                if (textwidth(ibuf, i) < (box_w - text_offset))
                        break;
        return i;
}

void clearrect()
{
        cleararea(ccol * box_w, hwidth + crow * box_h, box_w, box_h);
}

/* 
   handlechar has to be able to parse decimal numbers and strings,
   depending on the current column type, only printing characters should get this far
 */

void handlechar(WPARAM c)
{
        SEXP tvec;


        if (clength == 0) {
                if (length(inputlist) >= ccol + colmin - 1)
                        tvec = nthcdr(inputlist, ccol + colmin - 2);
                else
                        tvec = R_NilValue;
                if (crow == 0)  /* variable name */
                        currentexp = 3;
                else if (TYPEOF(CAR(tvec)) == STRSXP)   /* character data */
                        currentexp = 2;
                else
                        currentexp = 1;         /* numeric data */
                clearrect();
                highlightrect();
        }

        if (currentexp == 1)    /* we are parsing a number */
                switch (c) {
                case '-':
                        if (nneg == 0)
                                nneg++;
                        else
                                goto donehc;
                        break;
                case '.':
                        if (ndecimal == 0)
                                ndecimal++;
                        else
                                goto donehc;
                        break;
                case 'e':
                case 'E':
                        if (ne == 0) {
                                nneg = ndecimal = 0;    /* might have decimal in exponent */
                                ne++;
                        }
                        else
                                goto donehc;
                        break;
                default:
                        if (!isdigit(c))
                                goto donehc;
                        break;
                }
        if (currentexp == 3) {
                if (isspace(c))
                        goto donehc;
                if (clength == 0)
                        if (c != '.' && !isalpha(c))
                                goto donehc;
                        else if (c != '.' && !isalnum(c))
                                goto donehc;
        }

        if (clength++ > 29) {
                warning("spreadsheet: expression too long");
                clength--;
                goto donehc;
        }

        *bufp++ = c;
        printstring(buf, clength, crow, ccol);
        return;

      donehc:SysBeep();
}

void printlabs()
{
        char clab[10];
        int i;
        SEXP tppoint;

        if (length(inputlist) > colmin)
                tppoint = nthcdr(inputlist, colmin - 1);
        else
                tppoint = R_NilValue;

        for (i = colmin; i <= colmax; i++)
                if (TAG(tppoint) != R_NilValue) {
                        printstring(CHAR(PRINTNAME(TAG(tppoint))),
                                    strlen(CHAR(PRINTNAME(TAG(tppoint)))), 0, i - colmin + 1);
                        tppoint = CDR(tppoint);
                }
                else {
                        sprintf(clab, "var%d", i);
                        printstring(clab, strlen(clab), 0, i - colmin + 1);
                }
        for (i = rowmin; i <= rowmax; i++) {
                sprintf(clab, "R %d", i);
                printstring(clab, strlen(clab), i - rowmin + 1, 0);
        }
}

/* Menu Functions */

static void doDEMenu(WPARAM wParam, LPARAM lParam)
{
        switch (GET_WM_COMMAND_ID(wParam,lParam)) {
                case RRR_QUIT:
                        doneSpread=0;
        }
}

static int validName(char *text)
{
        char tmp;


        tmp=*text++;
        if (tmp != '.' && !isalpha(tmp))
                return 0;

        while (*text != '\0') {
                tmp=*text++;
                if ( tmp != '.' && !isalnum(tmp) )
                        return 0;
        }
        return 1;
}


/* 
   copyarea is a lot more complicated than you would expect 
   so you need to be sure the two regions are the same size
 */

static void copyarea(int src_x, int src_y, int dest_x, int dest_y)
{
        int destw, desth;
        HDC Ihdc;

        destw = (src_x < dest_x) ? windowWidth - dest_x : windowWidth - src_x;
        desth = (src_y < dest_y) ? windowHeight - dest_y : windowHeight - src_y;
        Ihdc=GetDC(RDEWnd);
        BitBlt(Ihdc, dest_x, dest_y, destw, desth, Ihdc, src_x, src_y, SRCCOPY);
        ReleaseDC(RDEWnd, Ihdc);
}

static void drawline(int fromx, int fromy, int tox, int toy)
{
        HDC devHdc;
        POINT lp;

        devHdc=GetDC(RDEWnd);
        if (xlast != fromx || ylast != fromy)
                MoveToEx(devHdc, fromx, fromy, &lp);
        LineTo(devHdc,tox, toy);
        xlast = tox;
        ylast = toy;
        ReleaseDC(RDEWnd, devHdc);
}


static void drawrectangle(int xpos, int ypos, int width, int height)
{
        HDC devHdc;

        devHdc=GetDC(RDEWnd);
        SelectObject(devHdc, GetStockObject(NULL_BRUSH));
        Rectangle(devHdc, xpos, ypos, xpos + width+1, ypos + height+1);
        SelectObject(devHdc, GetStockObject(WHITE_BRUSH));
        ReleaseDC(RDEWnd, devHdc);
}

static void setattribsfromwindow()
{
        windowWidth = DERect.right - DERect.left;
        windowHeight = DERect.bottom - DERect.top;
        box_w = textwidth(digits, 10);
        box_h = 25;
        bwidth = 0;
}

/* set the line width; use the current pen color */
static void setlineattribs(int width)
{
        HDC devHdc;
        HPEN hPen;

        devHdc=GetDC(RDEWnd);
        hPen = SelectObject(devHdc, CreatePen(PS_SOLID, width, RDEPCol));
        if (NotSP)
                DeleteObject(hPen);
        ReleaseDC(RDEWnd, devHdc);
        NotSP=1;
}

/* Text Drawing; */
static void drawtext(int xpos, int ypos, char *text, int len)
{
        HDC devHdc;

        devHdc=GetDC(RDEWnd);
        TextOut(devHdc,xpos, ypos, text, len);
        ReleaseDC(RDEWnd, devHdc);
}

/* cant query the pointer in Windoze */
static void querypointer(int *xglobal, int *yglobal, int *xlocal, int *ylocal)
{
        POINT tp1;

        *xlocal = DEMousePos.x;
        *ylocal = DEMousePos.y;
        tp1=DEMousePos;
        ClientToScreen(RDEWnd, &tp1);
        *xglobal = tp1.x;
        *yglobal = tp1.y;
}

/* find the width of a text string */
static int textwidth(char *text, int nchar)
{
        SIZE lpSize;
        HDC devHdc;

        devHdc=GetDC(RDEWnd);
        GetTextExtentPoint(devHdc,text, nchar, &lpSize);
        ReleaseDC(RDEWnd,devHdc);
        return lpSize.cx;
}


/* Open/Close Windows */
/* to make it the "right" size first set it up and then twiddle the width
   once we can get the text metrics
*/
static int DEInit(void)
{
        MDICREATESTRUCT mdicreate;
        HDC             Ihdc;
        TEXTMETRIC      Itm;
        RECT r;

        GetClientRect(RClient, (LPRECT) &r);
        mdicreate.szClass = RDEClass;
        mdicreate.szTitle = "R DataEntry";
        mdicreate.hOwner =(HINSTANCE) RInst;
                  mdicreate.x = r.left+20;
                  mdicreate.y = 20;
        mdicreate.cx = CW_USEDEFAULT;
        mdicreate.cy = CW_USEDEFAULT;
        mdicreate.style = 0;
        mdicreate.lParam=NULL;
        RDEWnd = (HWND) (UINT) SendMessage(RClient, WM_MDICREATE,0,
                (LONG) (LPMDICREATESTRUCT) &mdicreate);
        if( RDEWnd == NULL )
                return 0;



        Ihdc=GetDC(RDEWnd);
        GetTextMetrics(Ihdc, &Itm);
        SetBkMode(Ihdc, TRANSPARENT);
        ReleaseDC(RDEWnd, Ihdc);

        fh = Itm.tmHeight;
        fw = Itm.tmMaxCharWidth;

        /* find out how wide the input boxes should be and set up the window 
                size defaults */

        box_w = textwidth(digits, strlen(digits)) + 4;
        box_h = Itm.tmAscent + Itm.tmDescent + 4;
        nwide = 6; /* box_w/i;*/
        nhigh = 26; /*box_h/i;*/
        text_offset = 2 + Itm.tmDescent;
        windowWidth = 6 * box_w;
        windowHeight = 26 * box_h;
        SendMessage(RDEWnd, WM_SIZE, SIZE_RESTORED, 
                                MAKELONG(windowWidth, windowHeight));

        ShowWindow(RDEWnd, SW_SHOW);

        return 1;
}

LRESULT FAR PASCAL RDEWndProc(HWND hWnd, UINT message, WPARAM wParam,
        LPARAM lParam)
{
        HDC hdc;
        PAINTSTRUCT ps;
        int i;

        switch(message) {
                case WM_CREATE:
                        return 0;
                case WM_LBUTTONDOWN:
                case WM_RBUTTONDOWN:
                case WM_MBUTTONDOWN:
                        DEMousePos.x= LOWORD(lParam);
                        DEMousePos.y= HIWORD(lParam);
                        findsquare();
                        break;
                case WM_COMMAND:    
                        doDEMenu(wParam, lParam);
                        return 0;
                case WM_CHAR:
                if( ControlState ) {
                        if( wParam == 'f' ) 
                                jumpwin(colmin, rowmax);
                        if( wParam == 'b' ) {
                                i = (1 > rowmin - nhigh) ? 1 : rowmin - nhigh;
                                jumpwin(colmin, i);
                        }
                }
                else {
                        if( wParam == '\r' || wParam == '\t' || wParam == '\n' )
                         return 0;
                        else
                                handlechar(wParam);
                }
                return 0;
                case WM_KEYUP:
                        if( wParam == VK_SHIFT )
                                ShiftState=0;
                        if( wParam == VK_CONTROL )
                                ControlState=0;
                        break;
                case WM_KEYDOWN:
                        doDEKey( wParam);
                        break;
                case WM_PAINT:
                        hdc=BeginPaint(hWnd, &ps);
                        if( IsWindow(RDEWnd) )
                                drawwindow();
                        EndPaint(hWnd, &ps);
                        return 0;
                case WM_SETFOCUS:
                        SetFocus(hWnd);
                        if (IsWindow(RDEWnd))
                                SendMessage(RDEWnd, WM_MDIACTIVATE, (WPARAM) NULL,
                                        (LPARAM) RDEWnd);
                        break;
                case WM_MDIACTIVATE:
                        if((HWND) lParam == hWnd ) {
                                SendMessage(RClient,WM_MDISETMENU, (WPARAM) RMenuDE, (LPARAM) RMenuDEWin);
                                DrawMenuBar(RFrame);
                                return 0;
                        }
                        else if( doneSpread ) {
                            MessageBox(hWnd, "You must quit the data entry window before doing anything else",
                                "R Data Entry", MB_OK | MB_ICONEXCLAMATION);
                            PostMessage(RClient, WM_MDIACTIVATE, (UINT) RDEWnd,0);
                            return 0;
                        }
                        break;
                case WM_CLOSE:
                case WM_DESTROY:
                        doneSpread=0;
                        return 0;
        }
        return(DefMDIChildProc(hWnd, message, wParam, lParam));
}



static void CloseDE()
{
        closerect();
        SendMessage(RClient, WM_MDIDESTROY, (WPARAM) (HWND) RDEWnd, 0);
}

/* clear the rectangle with left top corner at xpos,ypos 
   and with width and height given */
static void cleararea(int xpos, int ypos, int width, int height)
{
        HDC Nhdc;
        RECT tre;

        tre.left=xpos;
        tre.top=ypos;
        tre.right=xpos+width;
        tre.bottom=ypos+height;
        Nhdc=GetDC(RDEWnd);
        FillRect(Nhdc, &tre, GetStockObject(WHITE_BRUSH));
        ReleaseDC(RDEWnd, Nhdc);
}

static void clearwindow()
{
        HDC Nhdc;

        Nhdc=GetDC(RDEWnd);
        FillRect(Nhdc, &DERect, GetStockObject(WHITE_BRUSH));
        ReleaseDC(RDEWnd, Nhdc);
}

/* Menus  */

/* comment this out for now */
void popupmenu(int col)
{
        PUPcol=col;
        DialogBox(RInst, "DEVarBox", RDEWnd, (DLGPROC) DEVar);
        SetFocus(RDEWnd);
}

extern BOOL CALLBACK DEVar(HWND hDlg, UINT message,
        WPARAM wParam,LPARAM lParam)
{
        char name[20];
        static SEXP tvec;
        int type, levs;


        switch (message) {
                case    WM_INITDIALOG:
                        if (length(inputlist) < PUPcol + colmin - 1)
                                inputlist = listAppend(inputlist, allocList(PUPcol + colmin - 1 - length(inputlist)));
                        tvec = nthcdr(inputlist, PUPcol + colmin - 2);
                        if (TAG(tvec) != R_NilValue)
                                sprintf(name,"%s", CHAR(PRINTNAME(TAG(tvec))));
                        else
                                sprintf(name, "Var%d", PUPcol + colmin - 1);
                        SetDlgItemText(hDlg, RDD_NAME, name);
                        SetWindowText(hDlg, name);
                        if( TYPEOF(CAR(tvec)) == STRSXP )
                                CheckRadioButton( hDlg, RDD_NUM, RDD_CHAR, RDD_CHAR);
                        else
                                CheckRadioButton( hDlg, RDD_NUM, RDD_CHAR, RDD_NUM);
                        return TRUE;
                case WM_COMMAND:
                        switch (wParam) {
                                case IDOK:
                                        GetDlgItemText(hDlg, RDD_NAME, name, 20);
                                        if (!validName(name) ) {
                                                MessageBox(RDEWnd,  "The variable name is not valid",
                                                        "R Application", MB_ICONEXCLAMATION | MB_OK);
                                                name[0]='\0';
                                                SetDlgItemText(hDlg, RDD_NAME, name);
                                                return 1;
                                        }
                                TAG(tvec) = install(name);
                                printstring(name, strlen(name), 0, PUPcol);
                                type = SendDlgItemMessage(hDlg, RDD_NUM, BM_GETCHECK,0,0);
                                if (type == 1) {
                                        if (CAR(tvec) == R_NilValue)
                                                CAR(tvec) = ssNewVector(REALSXP, 100);
                                        else {
                                                levs = LEVELS(CAR(tvec));
                                                CAR(tvec) = coerceVector(CAR(tvec), REALSXP);
                                                LEVELS(CAR(tvec)) = levs;
                                        }
                                }
                                else {
                                        if (CAR(tvec) == R_NilValue)
                                                CAR(tvec) = ssNewVector(STRSXP, 100);
                                        else {
                                                levs = LEVELS(CAR(tvec));
                                                CAR(tvec) = coerceVector(CAR(tvec), STRSXP);
                                                LEVELS(CAR(tvec)) = levs;
                                        }
                                }
                                        EndDialog(hDlg, TRUE);
                                        return TRUE;
                                case IDCANCEL:
                                        EndDialog(hDlg, FALSE);
                                        return TRUE;
                                case RDD_NUM:
                                case RDD_CHAR:
                                        CheckRadioButton(hDlg, RDD_NUM, RDD_CHAR, wParam);
                                        return TRUE;   
                                }

                }
                return FALSE;
}
