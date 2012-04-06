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

#include "Graphics.h"
#include <stdio.h>
#include <math.h>
#include "generic.h"
#include <windows.h>
/* no printing yet 
#include <Printing.h>
*/

/* OUTSTANDING ISSUES		*/
/* how to handle the locator under Windows? */
/* Only True Type Fonts can be rotated so we need to choose one
    -GetRasterizerCaps indicates whether TrueType is available on
 		the machine
    -EnumFontFamilies can be called to enumerate available fonts
	and you can check these to see if they are TrueType
    -CHOOSEFONT (the interactive dialog) has a flag that allows
	only TrueType fonts to be enumerated (should we allow
	only TrueType??)

   MetaFiles: a means of enabling print (the easy way is bitmaps
	but I doubt all printers will support this).
	-perhaps this will allow us to redraw the screen when it is
		resized

*/

#define G_PI          3.141592653589793238462643383276
int FontSize=9;
static LOGFONT RGraphLF;
static char szRectClass[] = "MdiRectClass";

extern void SysBeep(void)

/* Basic Win Objects */
RECT			graphicsRect = { 0, 0, 300, 400 };
int			fontw, fonth;

static void CopyGraphToScrap(void);

/* need a function to change fonts because things like fontw and fonth
   need to change at the same time. Because there is some approximating
   of fonts you need to first select the new font into the Display Context
   and then check to see how big it is.
*/

static void Win_RGSetFont(int fontsize)
{
	HDC		Ihdc;
	TEXTMETRIC	Itm;

	Ihdc=GetDC(RGraphWnd);
	RGraphLF.lfHeight= -fontsize*GetDeviceCaps(Ihdc,LOGPIXELSY)/72;
	SelectObject(Ihdc, CreateFontIndirect(&RGraphLF));
	GetTextMetrics(Ihdc, &Itm);
	fonth = Itm.tmHeight;
	fontw = Itm.tmMaxCharWidth;
	ReleaseDC(RGraphWnd);
}

LRESULT FAR PASCAL GraphWndProc(HWND hWnd, UINT message, WPARAM wParam,
        LPARAM lParam)
{
        HDC hdc;
        PAINTSTRUCT ps;

        switch(message) {
                case WM_CREATE:
			/* initialize the Graphics Logfont */
                        lstrcpy(RGraphLF.lfFaceName, "Ariel");
			Win_RGSetFont(FontSize);
                        break;
                case WM_SIZE:
			/* reset the graphicsRect; can we redraw the picture  */
			GetClientRect(RGraphWnd, &graphicsRect);
			Win_NewPlot();
                        break;
                case WM_LBUTTONDOWN:
                        break;
                case WM_RBUTTONDOWN:
                        break;
                case WM_COMMAND:
			doGraphicsMenu(hWnd, wParam, lParam);
                case WM_PAINT:
                        hdc=BeginPaint(hWnd, &ps);
                        EndPaint(hWnd, &ps);
                        return 0;
                case WM_SETFOCUS:
                        SetFocus(RGraphWnd);
                        return 0;
                case WM_MDIACTIVATE:
                        if((HWND) lParam == hWnd ) {
                            SendMessage(hWndClient, WM_MDISETMENU, 
				(WPARAM) RMenuGraph, (LPARAM) RMenuGraphWin);
                            DrawMenuBar(hWndFrame);
                        }
                        return(0);
                case WM_QUERYENDSESSION:
                        break;
                case WM_CLOSE:
                        ShowWindow(hWnd, SW_MINIMIZE);
                        return(0);
                case WM_DESTROY:
                        return(0);
         }
  return(DefMDIChildProc(hWnd, message, wParam, lParam));
}


static void doGraphicsMenu(HWND GWnd, WPARAM wParam, LPARAM lParam)
{
	switch (wParam) {
		case IDM_SETUP:
			SysBeep(2);
			/*
			SetupPrinter(); */
			break;
		case IDM_PRINT:
			SysBeep(2);
			/* DoPlotPrint(graphicsPicture); */
			break;
		case IDM_COPY:		
			CopyGraphToScrap();
			break;
		}
}

void InitGraphicsContext(void)
{
	GetClientRect(RGraphWnd, &graphicsRect);
	Win_RGSetFont(FontSize);
	
	/* graphContext.active = 0;	1 only when locator in action */
}

/* Current Clipping Rectangle */
static double Clipxl, Clipxr, Clipyb, Clipyt;

/* Last Point Coordinates */
static int xlast = 0;
static int ylast = 0;

/* Open a Graphics Window and set Graphics State */
static int Win_Open(void)
{
	MDICREATESTRUCT mdicreate;
	HINSTANCE hInst;

	hInst=GetWindowLong(hWndFrame, GWL_HINSTANCE);
	mdicreate.szClass = szRectClass;
        mdicreate.szTitle = "R Graphics";
        mdicreate.hOwner = hInst;
        mdicreate.x = CW_USEDEFAULT;
        mdicreate.y = CW_USEDEFAULT;
        mdicreate.cx = min(r.right-r.left,r.bottom-r.top);
        mdicreate.cy = min(r.right-r.left,r.bottom-r.top);
        mdicreate.style = 0;
        mdicreate.lParam=NULL;
        RGraphWnd = (HWND) (UINT) SendMessage(hWndClient, WM_MDICREATE,0,
                (LONG) (LPMDICREATESTRUCT) &mdicreate);
	if( RGraphWnd == NULL )
		return 0;

	ShowWindow(RGraphWnd, SW_SHOW);

	InitGraphicsContext();
	DevInit = 1;
	return 1;
}

/* Set the Clipping Rectangle */
static void Win_Clip(int x0, int x1, int y0, int y1)
{
	if(x0 < x1) {
		Clipxl = x0;
		Clipxr = x1;
	}
	else {
		Clipxl = x1;
		Clipxr = x0;
	}
	if(y0 < y1) {
		Clipyb = y0;
		Clipyt = y1;
	}
	else {
		Clipyb = y1;
		Clipyt = y0;
	}
}

/* Actions on Window Resize */
static void Win_Resize()
{
	DP->right = graphicsRect.right;
	DP->bottom = graphicsRect.bottom;
}

/* Begin a New Plot; under Win32 just clear the graphics rectangle. */

static void Win_NewPlot()
{
	HDC Nhdc;
	
	Nhdc=GetDC(RGraphWnd);
	FillRect(Nhdc, &graphicsRect, GetStockObject(WHITE_BRUSH));
	ReleaseDC(RGraphWnd, Nhdc);
}

/* Close the Graphics Window */
static void Win_Close()
{
}

/* MoveTo */
static void Win_MoveTo(int x, int y)
{
	xlast = x;
	ylast = y;
	devHdc=GetDC(RGraphWnd);
	MoveTo(devHdc,x,y);
	ReleaseDC(RGraphWnd,devHdc);
}

/* Dot */
static void Win_Dot(void)
{
}



/* LineTo */
static void Win_LineTo(int x, int y)
{
	devHdc=GetDC(RGraphWnd);
	LineTo(devHdc,x, y);
	ReleaseDC(RGraphWnd,devHdc);
	xlast = x;
	ylast = y;
}

/* Draw a Filled Rectangle */
static void Win_Rect(int x0, int y0, int x1, int y1, int col, int fill)
{
}

static double cex;
static int fontsize;

/* Horizontal Text Drawing */
static void Win_Text(char *str, double xc, double yc)
{
	int x, y;
	LPSIZE lpSize;

	if(cex != GP->cex) {
		cex = GP->cex;
		fontsize = FontSize * cex;
		Win_RGSetFont(fontsize);
	}

	devHdc=GetDC(RGraphWnd);
	GetTextExtentPoint(devHdc,str,strlen(str),lpSize);
	x = -xc*lpSize.cx;
	y = yc*lpSize.cy;
	MoveTo(xlast+x, ylast+y);
	TextOut(devHdc,x,y,str,strlen(str));
	ReleaseDC(RGraphWnd,devHdc);
}

/* Rotated Text   
   xlast and ylast are the the lower left corner of the box
   xc and yc indicate text justification
   0 means left justified,
   1 means right justified
   0.5 means centered

  Under Win32 only TrueType Fonts can be rotated. This means
  that the user will have to select one of these for us to have
  rotated text.
*/
static void Win_RText(char *str, double xc, double yc, int rot)
{
	int fontsize,x,y,nstr;
	double rotrad, xoff, yoff;
	HDC devHdc;
	LPSIZE lpSize;

	if( RGraphLF.lfEscapement != rot*10 || cex != GP->cex) {
		RGraphLF.lfEscapement=rot*10;
		cex = GP->cex;
		fontsize = FontSize * cex;
		Win_RGSetFont(fontsize);
	}

	devHdc=GetDC(RGraphWnd);
	nstr = strlen(str);
	GetTextExtentPoint(devHdc,str,nstr,lpSize);
	x = xc * lpSize.cx;
	y = yc * lpSize.cy;
	rotrad=rot*(2*G_PI)/360;
	xoff = x*cos(rotrad)+y*sin(rotrad);
	yoff = x*sin(rotrad)+y*cos(rotrad);
	x= (int) xoff;
	y= (int) yoff;
	TextOut(devHdc, xlast- x, ylast+y, str, nstr);
	ReleaseDC(RGraphWnd,devHdc);
	return;
}

/* Return the Pointer Location */
static int Win_Locator(int *x, int *y)
{
	return 0;
}

static WindowPtr savePort;

/* Set the Graphics Mode */
static void Win_Mode(int mode)
{
}

/* Keep the Graphics Window in Front */
static void Win_Hold(void)
{
}

static void CopyGraphToScrap(void)
{
	HDC hDC, hDCMem;
	HBITMAP hBM, hOldBM;
	int left, top, width, height;

	hDC = GetDC(RGraphWnd);
	hDCMem = CreateCompatibleDC(hDC);
	left=graphicsRect.left;
	top=graphicsRect.top;
	width=graphicsRect.right-graphicsRect.left;
	height=graphicsRect.top-graphicsRect.bottom;
	hBM = CreateCompatibleBitmap(hDC, width, height);
	if (hBM) {
		SelectObject(hDCMem, hBM);
		BitBlt(hDCMem, left, top, width, height, hDC, left, 
			top, SRCCOPY);
		OpenClipboard(RGraphWnd);
		EmptyClipboard();
		SetClipboardData(CF_BITMAP, hBM);
		CloseClipboard();
	}
	else SysBeep(10);

	DeleteDC(hDCMem);
	ReleaseDC(RGraphWnd, hDC);
  }
}

/* Device Driver */
WinDeviceDriver()
{
	HDC DDhdc;

	DevInit = 0;
	if( ! Win_Open() ) return 0;

	DevOpen = Win_Open;
	DevClose = Win_Close;
	DevResize = Win_Resize;
	DevNewPlot = Win_NewPlot;
	DevClip = Win_Clip;
	DevMoveTo = Win_MoveTo;
	DevLineTo = Win_LineTo;
	DevText = Win_Text;
	DevRText = Win_RText;
	DevDot = Win_Dot;
	DevRect = Win_Rect;
	DevLocator = Win_Locator;
	DevMode = Win_Mode;
	DevHold = Win_Hold;

	GP->left = 0;
	GP->right = graphicsRect.right - graphicsRect.left;
	GP->bottom = graphicsRect.bottom - graphicsRect.top;
	GP->top = 0;

	DDhdc=GetDC(RGraphWnd);

	/* character size in raster */

	GP->cra[0] = fontw;
	GP->cra[1] = FontSize;

	GP->xCharOffset = 0.0;
	GP->yCharOffset = 0.0;

	/* inches per raster x and then y */
	GP->ipr[0] = 1/GetDeviceCaps(DDhdc, LOGPIXELSX);
	GP->ipr[1] = 1/GetDeviceCaps(DDhdc, LOGPIXELSY);

	GP->canResizePlot = 1;
	GP->canChangeFont = 0;
	GP->canRotateText = 1;
	GP->canResizeText = 0;
	GP->canClip = 0;

	DevInit = 1;
	cex = 1.0;
	
	ReleaseDC(RGraphWnd, DDhdc);

	return 1;
}
