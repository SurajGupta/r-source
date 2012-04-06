 /*  R : A Computer Langage for Statistical Data Analysis
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


#include "wincons.h"

/*  TO DO

  Need to adjust RPrintText somewhat so that it only tries to print
  lines of text that are as wide as the printer. Perhaps some sort of
  line breaking procedure is needed.
  The user should be able to set the font somehow (perhaps we should use
  the same font as is current in the display window).

  RPrintGraph: should only need the Metafile and the size of the PlotWin
  to do its thing.
*/

BOOL bPrint;
HWND hdlgCancel;

void RPrintGraph(HWND hwnd, HANDLE RMhmf, HWND PlotWin)
{
        PRINTDLG pd;
        DOCINFO di;
        SIZE p1;
        RECT rect;
        int nError;


        pd.lStructSize = sizeof(PRINTDLG);
        pd.hDevMode = (HANDLE) NULL;
        pd.hDevNames = (HANDLE) NULL;
        pd.Flags = PD_RETURNDC;
        pd.hwndOwner = hwnd;
        pd.hDC = (HDC) NULL;
        pd.nFromPage = 1;
        pd.nToPage = 1;
        pd.nMinPage = 0;
        pd.nMaxPage = 0;
        pd.nCopies = 1;
        pd.hInstance = (HANDLE) NULL;
        pd.lpfnPrintHook = (LPPRINTHOOKPROC) NULL;
        pd.lpfnSetupHook = (LPPRINTHOOKPROC) NULL;
        pd.lpPrintTemplateName = (LPSTR) NULL;
        pd.lpSetupTemplateName = (LPSTR) NULL;
        pd.hPrintTemplate = (HANDLE) NULL;
        pd.hSetupTemplate = (HANDLE) NULL;

        PrintDlg(&pd);

        bPrint = TRUE;

        SetAbortProc(pd.hDC, AbortProc);

        hdlgCancel = CreateDialog(RInst, (LPSTR) "PAbortDlg", hwnd,
                        (DLGPROC) AbortPrintJob);

        EnableWindow(hwnd, FALSE);

        di.cbSize = sizeof(DOCINFO);
        di.lpszDocName = "Test 1";
        di.lpszOutput = (LPSTR) NULL;

        nError = StartDoc(pd.hDC, &di);
        if( nError == SP_ERROR )
                goto Error;

        nError = StartPage(pd.hDC);
        if( nError <= 0 )
                goto Error;

        SetMapMode(pd.hDC, MM_ISOTROPIC);
        SetBkMode(pd.hDC, TRANSPARENT);
        GetClientRect(PlotWin, &rect);
        SetWindowExtEx(pd.hDC,rect.right,rect.bottom,&p1);
        SetViewportExtEx(pd.hDC, GetDeviceCaps(pd.hDC, HORZRES), GetDeviceCaps(pd.hDC, VERTRES),NULL);
        PlayMetaFile(pd.hDC, RMhmf);

        nError = EndPage(pd.hDC);
        if(nError <= 0 )
                goto Error;

        EndDoc(pd.hDC);

Error:
        EnableWindow(hwnd, TRUE);
        DestroyWindow(hdlgCancel);
        DeleteDC(pd.hDC);
}

void RPrintText(HWND hwnd, HWND TextWin)
{
        PRINTDLG pd;
        DOCINFO di;
        TEXTMETRIC tm;
        int nError, nTotalLines, ychar, nLinesPerPage, nTotalPages, nPage;
        int nNonColCopy, nColCopy, nLine, nLineNum, nCharsPerLine, nchars;
        char rbuf[256];

        pd.lStructSize = sizeof(PRINTDLG);
        pd.hDevMode = (HANDLE) NULL;
        pd.hDevNames = (HANDLE) NULL;
        pd.Flags = PD_ALLPAGES | PD_COLLATE | PD_RETURNDC;
        pd.hwndOwner = hwnd;
        pd.hDC = (HDC) NULL;
        pd.nFromPage = 0;
        pd.nToPage = 0;
        pd.nMinPage = 0;
        pd.nMaxPage = 0;
        pd.nCopies = 1;
        pd.hInstance = (HANDLE) NULL;
        pd.lpfnPrintHook = (LPPRINTHOOKPROC) NULL;
        pd.lpfnSetupHook = (LPPRINTHOOKPROC) NULL;
        pd.lpPrintTemplateName = (LPSTR) NULL;
        pd.lpSetupTemplateName = (LPSTR) NULL;
        pd.hPrintTemplate = (HANDLE) NULL;
        pd.hSetupTemplate = (HANDLE) NULL;



        if( !PrintDlg(&pd))
                goto Error2;

   SetAbortProc(pd.hDC, AbortProc);

        hdlgCancel = CreateDialog(RInst, (LPSTR) "PAbortDlg", hwnd,
                        (DLGPROC) AbortPrintJob);

        EnableWindow(hwnd, FALSE);
        nTotalLines=Edit_GetLineCount(TextWin);
        if (nTotalLines==0)
                goto Error;

        GetTextMetrics( pd.hDC, &tm);
        ychar=tm.tmHeight+tm.tmExternalLeading;
        nCharsPerLine=GetDeviceCaps(pd.hDC, HORZRES)/tm.tmAveCharWidth;
        nLinesPerPage=GetDeviceCaps(pd.hDC, VERTRES)/ychar;
        nTotalPages=(nTotalLines+nLinesPerPage -1)/nLinesPerPage;

        bPrint = TRUE;


        di.cbSize = sizeof(DOCINFO);
        di.lpszDocName = "Test 1";
        di.lpszOutput = (LPSTR) NULL;

        nError = StartDoc(pd.hDC, &di);
        if( nError == SP_ERROR )
                goto Error;


        for( nColCopy=0; nColCopy< (pd.Flags & PD_COLLATE ? pd.nCopies : 1);
                        nColCopy++) {
                                for( nPage=0 ; nPage < nTotalPages ; nPage++ ) {
                                        for (nNonColCopy=0; nNonColCopy < (pd.Flags & PD_COLLATE ? 1 : pd.nCopies);
                                                nNonColCopy++ )
                                                {
                                                        nError = StartPage(pd.hDC);
                                                        if( nError <= 0 )
                                                                goto Error;
                                                        for(nLine=0;nLine<nLinesPerPage; nLine++) {
                                                                nLineNum=nLinesPerPage*nPage+nLine;
                                                                if( nLineNum >= nTotalLines )
                                                                        break;
                                                                nchars=Edit_GetLine(TextWin, nLineNum, rbuf,256);
                                                                TextOut(pd.hDC, 0, ychar*nLine, rbuf, nchars);
                                                        }
                                                        nError = EndPage(pd.hDC);
                                                        if(nError <= 0 )
                                                                goto Error;
                                        }
                                }
        }
        EndDoc(pd.hDC);

Error:
        EnableWindow(hwnd, TRUE);
        DestroyWindow(hdlgCancel);
        DeleteDC(pd.hDC);
Error2:
        return;
}

BOOL CALLBACK AbortProc(HDC hdcPrn, int nCode)
{
        MSG msg;

        while (bPrint && PeekMessage((LPMSG) &msg, (HWND) NULL, 0,0, PM_REMOVE))
        {
                if(!IsDialogMessage(hdlgCancel, (LPMSG) &msg)) {
                        TranslateMessage((LPMSG) &msg);
                        DispatchMessage((LPMSG) &msg);
                }
        }
        return bPrint;
}

LRESULT CALLBACK AbortPrintJob(HWND hDlg, UINT message, WPARAM wParam,
                        LPARAM lParam)
{
        switch(message) {
                case WM_INITDIALOG:
                        /* SetDlgItemText(hDlg, RDD_FILE, ofn.lpstrFile); */
                        return TRUE;
                case WM_COMMAND:
                        bPrint = FALSE;
                        return TRUE;
                default:
                        return FALSE;
        }
}
