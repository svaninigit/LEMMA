#ifndef TOMMAINFRAME_h
#define TOMMAINFRAME_h

// ROOT classes
#include <TGFrame.h>
#include <TGButton.h>
#include <TGMenu.h>
#include <TGTab.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGButtonGroup.h>
#include <TGToolBar.h>
#include <TGCanvas.h>

// TOM classes
#include "Monitor/RawHistos.h"
#include "PattRec/RawAnalyzer.h"

///////////////////////////////////////////////////////////////////////
// Muon TOMography Tool:                                             //   
// Class that defines the main Window                               //
// Sara Vanini 100110                                                //
///////////////////////////////////////////////////////////////////////
class TOMMainFrame : public TGMainFrame {

public:
  TOMMainFrame(const TGWindow *p, UInt_t w, UInt_t h, int argc, char** argv);
  virtual ~TOMMainFrame();
  void CloseWindow();
  Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

private:
  //MEMBER FUNCTIONS
  // raw
  int   OpenRawFile(int, char*);
  int   CloseRawFile();
  int   ReadFile();
  int   FillRawHistos(bool flagDown);
  //calib
  int   ComputeTTrig();
  int   DumpCalibHistoFile();
  int   DumpTTrigs(); 
  int   DumpRawHistoFile();

  //general
  int   HandleRadio(Long_t parm1);
  int   HandleCheck(Long_t parm1, bool flag);
  int   FillHistogram(Int_t padindex, const char *chname);
  int   Fill2DHistogram(Int_t padindex, const char *chname);
  int   FillGraph(Int_t padindex, const char *chname);
  int   GetStatus();

private:  
   //GUI OBJECTS --  Main Window
  TGButton              *fCloseButton;
  TGCompositeFrame      *fButtonFrame;
  TGMenuBar             *fMenuBar;
  TGPopupMenu           *fMenuFile;
  TGTab 		*fTab;
  TGToolBar 		*tb;

    //GUI OBJECTS -- Tabs
  TGCompositeFrame 	*fCompositeFrameTab1, *fCompositeFrameTab2, 
			*fCompositeFrameTab3, *fCompositeFrameTab4, *fCompositeFrameTab5;
  TGTextButton 		*fTextButtonTab1, *fTextButtonTab2n1, *fTextButtonTab2n2, *fTextButtonTab2n3, 
			*fTextButtonTab3, *fTextButtonTab4, *fTextButtonTab5;
  TGCheckButton		*fCheckButtonTab1, *fCheckButtonTab2n1, 
			*fCheckButtonTab3, *fCheckButtonTab4, *fCheckButtonTab5;
  TGCheckButton		*fCheckButtonFileTab1, *fCheckButtonFileTab2, *fCheckButtonFileTab3, 
			*fCheckButtonFileTab5;
  TGIcon 		*fIconTab1, *fIconTab2;
  TGLabel 		*fLabelTab1, *fLabelTab2;
  TGButtonGroup         *fRadioFrameTab1, *fRadioFrameTab2, *fRadioFrameTab3;
  TGRadioButton         *fRadiosTab1[4], *fRadiosTab2[2], *fRadiosTab3[3];


    //OBJECTS USED TO ORGANIZE AND PRESENT DATA
  TCanvas               *fC1;

    // other TOM classes
  RawHistos 	     * _rawHistos; 
  TTrigCalibration   * _ttrigCalib; 
  RawAnalyzer  	     * _rawAnalyzer;

   //NON-OBJECT VARIABLES -- uses ROOT defined data types for portability
  char  rawFileName[128];


};

#endif /*TOMMAINFRAME_h*/
