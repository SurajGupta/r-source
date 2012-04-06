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

#include "wincons.h"
#include "Graphics.h"

static char szDirName[RBuffLen];


        /*--- I n i t i a l i z a t i o n -- C o d e ---*/

static BOOL CheckSystem(void)
{
  DWORD vinfo = GetVersion();

  if (! ((LOBYTE(LOWORD(vinfo)) > 3) ||
         (LOBYTE(LOWORD(vinfo)) == 3 && HIBYTE(LOWORD(vinfo)) >= 10))) {
    MessageBox((HWND) NULL,
               "R requires Windows 3.1 or higher",
               NULL,
               MB_ICONHAND);
    return(FALSE);
  }
  return(TRUE);
}


extern HMENU RMenuEdit;

int WINAPI WinMain(HANDLE hinstCurrent, HANDLE hinstPrevious, 
LPSTR lpszCmdParam, int nCmdShow)
{
        int i;
        char *exe;

        if (! CheckSystem()) return(FALSE);

        /* Create the Windows */
        if( !hinstPrevious )
                if( !InitApplication(hinstCurrent) )
                        return FALSE;

        if( !InitInstance(hinstCurrent, nCmdShow) )
                return FALSE;

        /* need to get the working directory */
        i = GetModuleFileName(NULL, szDirName, RBuffLen);
        if (i > RBuffLen || i==0)
                return FALSE;
        exe = strrchr(szDirName,'\\');
        *exe = '\0';
        setenv("RHOME",szDirName,1);
        setenv("HOME",szDirName,1);
        strcpy(R_ImageName,szDirName);
        strcat(R_ImageName,"\\.RData");

        mainloop();

        DestroyMenu(RMenuDE);
        DestroyMenu(RMenuGraph);
        DestroyMenu(RMenuInit);
        DestroyMenu(RMenuConsole);
        DestroyMenu(RMenuEdit);
}

void R_StartUp(void)
{
        R_Init = 1;
}

void R_Busy(int yes)
{
        /* Placeholder */
}

        /*--- I / O --S u p p o r t -- C o d e ---*/


        /*--- P l a t f o r m -- D e p e n d e n t -- F u n c t i o n s ---*/



SEXP do_machine(SEXP call, SEXP op, SEXP args, SEXP env)
{
        return mkString("Win32");
}

SEXP do_system(SEXP call, SEXP op, SEXP args, SEXP rho)
{
        errorcall(call, "\"system\" is only available on Unix");
}

/* this function should set up the file pointer and open the file
   on some systems it should check to make sure that there is room to
   store the image
*/


void dump_image(char* fname, int jump)
{
                FILE *fp;
                long Vsize;
                DWORD Clust, FreeClust, SectPerClust, BytesPerSect;
                char    tstr[2];

        fp = fopen(fname, "wb");
        if( !fp )
                error("can't save data -- unable to open file\n");

        /* check to see if another drive was specified; 
           must be a right way to do this one;
           currently (1997) this size stuff is being ignored
        */
        if(fname[1] == ':' )
                strncpy(tstr,fname,2);
        else
                tstr[0]='\0';

        if( strlen(tstr) > 0 )
                GetDiskFreeSpace(tstr,&BytesPerSect, &SectPerClust, &FreeClust, &Clust);
        else
                GetDiskFreeSpace(NULL,&BytesPerSect, &SectPerClust, &FreeClust, &Clust);
        Vsize=BytesPerSect*SectPerClust*FreeClust;


        R_WriteMagic(fp, R_MAGIC_BINARY);
        BinarySave(FRAME(R_GlobalEnv), fp);
        fclose(fp);
        if(jump)
                jump_to_toplevel();
}

void R_InitialData(void)
{
        R_RestoreGlobalEnv();
}

void RBusy(int which)
{
}

void R_SaveGlobalEnv(void)
{
        FILE *fp = fopen(R_ImageName, "w");
        if (!fp)
                error("can't save data -- unable to open %s\n",R_ImageName);
        R_WriteMagic(fp, R_MAGIC_BINARY);
        BinarySave(FRAME(R_GlobalEnv), fp);
        fclose(fp);
}               

void R_RestoreGlobalEnv(void)
{                       
        FILE *fp = fopen(R_ImageName,"r");
        if (!fp) {      
                /* warning here perhaps */
                return;
        }               
        if(!R_Quiet) 
                Rprintf("[Previously saved workspace restored]\n\n");
                        
        switch(R_ReadMagic(fp)) {
        case R_MAGIC_BINARY:
                FRAME(R_GlobalEnv) = BinaryLoad(fp);
                break;
        case R_MAGIC_ASCII:
                FRAME(R_GlobalEnv) = AsciiLoad(fp);
                break;
        default:
                fclose(fp);
                error("workspace file corrupted -- no data loaded\n");
        }
}

