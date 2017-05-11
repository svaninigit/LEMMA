/****************************************************************************
 *
 *    main.c -- A Tool for Muon Tomography 
 *
 *    This program includes the ROOT libraries 
 *
 *    This program was developed for the Muon Tomography Experiment at Padua
 *    University, by Sara Vanini.   
 *
 *    S. Vanini January 2010
 ****************************************************************************
 */

/*****************************************************************************
 *       FUNCTIONALITY                                  VERSION 1.10
 *****************************************************************************
 *
 *    This program, can accept ... command line arguments. The 1st is
 *      the name....  The second argument, if given, instructs the 
 *      program to..
 *
 *    This program generates 2 windows:  a large canvas drawing area and a 
 *      command window.  By selecting options in the command window and clicking
 *      one of the draw buttons, a group of graphs and/or histograms will be 
 *      displayed to the canvas window.
 *
 *    THE CANVAS WINDOW:
 *      This window is constructed as a standard ROOT canvas.  As such, it has
 *      the ROOT canvas toolbar.  From this toolbar it is possible to save
 *      canvas as a .gif, .ps or .eps file, view the colors and markers available,
 *      and several other things.  Since this program is built as a 'Stand-Alone'
 *      type and does not operate out of ROOT's cint, you may find that some 
 *      of the toolbar options do not work.
 *
 *    THE COMMAND WINDOW:
 *      This window consists of five tabs, a menu bar, and a bar of selection
 *      boxes and buttons (the 'Button Bar').
 *
 *      TAB 1 -- ....
 *         From here, you can 
 *
 *      TAB 2 -- ....
 *         From here, ...  
 *
 *
 *      MENU BAR
 *         This bar currently only has one pop-up menu.  From this menu, you can .....
 *         You can also close the menu and exit the program from here.
 *         
 *   OBJECTS IN MEMORY
 *     When a Histogram or graph is drawn, the object is kept in memory purposefully.  This
 *     allows the user to move the contents of the canvas around with the mouse.  For Lego
 *     and Surface plots, the plot can be rotated to provide a better view.  If the objects
 *     are deleted after they are drawn, the picture will remain, but clicking anywhere
 *     in the canvas will cause the plot to disappear.
 *
 *   NETSCAPE NOTE
 *     If Netscape is active when this program is run, this program will not be able to select
 *     and use many of the normal colors.  Netscape can use a command line option to prevent
 *     it from eating the colormaps completely, but I don't recall exactly what it is.  I have
 *     selected as selectable colors those colors which seem to be available regardless of
 *     whether Netscape is running.  If you run this program and see much more white than usual,
 *     close Netscape!
 *
 *
 *                                                     Sara Vanini
 *                                                     15 January 2010
 *   
 *   LATEST IMPROVEMENT -- 100115
 */

/* system headers */
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h> 
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>

/*  Headers for ROOT */
#include <TROOT.h>
#include <TApplication.h>

#include "TOMMainFrame.h"

//////////////////////////////////////////////////////////////
//Statements that are included in almost every ROOT GUI app //
//////////////////////////////////////////////////////////////
extern void InitGui();
VoidFuncPtr_t initfuncs[] = {InitGui, 0 };
TROOT root("GUI", "Muon TOMography tool", initfuncs);

int main(int argc, char **argv)
{
  for (int i=0; i<argc; i++)
    if ( strstr(argv[i], "-help")!=NULL ) {
      cout << "Usage:    runTOM -f [file name] " << endl;
      cout << "Example:  runTOM -f r395.i1 opens the file r395.i1 " << endl;
      return 0;
    }

  //Standard main function
  TApplication theApp("App", &argc, argv);
  TOMMainFrame mainWindow(gClient->GetRoot(), 500, 400, argc, argv);
  theApp.Run();
  return 0;
}


