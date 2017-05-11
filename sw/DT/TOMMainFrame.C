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
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGMsgBox.h>
#include <TGDockableFrame.h>
#include <TG3DLine.h>
#include <TGButtonGroup.h>
#include <TGToolBar.h>
#include <TGFileDialog.h>

#include "TOMMainFrame.h"

bool DEBUG_TOM = true;

enum EMessageID {
  M_FILE_OPEN, M_FILE_CLOSE, M_QUIT, B_QUIT
};

const char *filetypes[] = { "TOMography raw files", "*.i*",
			    "All files",     "*",
			    0,               0 };

TOMMainFrame::TOMMainFrame(const TGWindow *p, UInt_t w, UInt_t h, int argc, char** argv)
: TGMainFrame(p,w,h) {  
 
  //Init variables
  _rawHistos 	= 0;
  _ttrigCalib 	= 0;
  fC1 		= 0;
  _rawAnalyzer = new RawAnalyzer();
  strcpy(rawFileName, "");

  // menu file
  fMenuFile = new TGPopupMenu(gClient->GetRoot());
  fMenuFile->AddEntry("&Open Raw File...", M_FILE_OPEN);
  fMenuFile->AddEntry("&Close Raw File...", M_FILE_CLOSE);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("&Exit", M_QUIT);
  fMenuFile->Associate(this);

  // menu bar and add popup menus
  fMenuBar = new TGMenuBar(this,100,20,kHorizontalFrame);
  fMenuBar->AddPopup("&File",fMenuFile, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
  AddFrame(fMenuBar, new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1));

  // line
  TGHorizontal3DLine *fHLine = new TGHorizontal3DLine(this,500,2); // p,w,h
  AddFrame(fHLine, new TGLayoutHints(kLHintsLeft | kLHintsBottom | kLHintsExpandX,2,2,2,2));
  fHLine->MoveResize(0,22);

  // tool bar
  tb = new TGToolBar(this,500,30);
  tb->MoveResize(0,30);
  AddFrame(tb, new TGLayoutHints(kLHintsTop | kLHintsExpandX ));

  // quit button
  fButtonFrame = new TGHorizontalFrame(this, 500, 30,  kFixedWidth);
  fCloseButton = new TGTextButton(fButtonFrame, "&Quit", B_QUIT);
  fCloseButton->SetToolTipText("Closes and exits program");
  fCloseButton->Associate(this);
  fButtonFrame->AddFrame(fCloseButton, new TGLayoutHints(kLHintsBottom | kLHintsRight, 2, 2, 5, 1));
  AddFrame(fButtonFrame, new TGLayoutHints(kLHintsBottom, 0, 0, 0, 0));
  fButtonFrame->MoveResize(0,370);
  
  // main frame already done in main...
  SetLayoutBroken(kTRUE);

  // tab widget
  fTab = new TGTab(this,500,320);

  // TAB1
  	// container of "Tab1"
  fCompositeFrameTab1 = fTab->AddTab("Raw Data");
  fCompositeFrameTab1->SetLayoutManager(new TGVerticalLayout(fCompositeFrameTab1));
  fCompositeFrameTab1->SetLayoutBroken(kTRUE);

	//button
  fTextButtonTab1 = new TGTextButton(fCompositeFrameTab1,"Read Raw Data", 100);
  fTextButtonTab1->SetTextJustify(36);
  fTextButtonTab1->SetMargins(0,0,0,0);
  fTextButtonTab1->SetWrapLength(-1);
  fCompositeFrameTab1->AddFrame(fTextButtonTab1, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fTextButtonTab1->MoveResize(80,30,300,50);
  fTextButtonTab1->Associate(this);

	//check button
  fCheckButtonTab1 = new TGCheckButton(fCompositeFrameTab1, "Fill Histograms", 101);
  fCheckButtonTab1->MoveResize(80,100,150,20);
  fCheckButtonTab1->Associate(this);

  	//radio buttons
  fRadioFrameTab1 = new TGButtonGroup(fCompositeFrameTab1, "Histograms: occupancy & time boxes", kVerticalFrame);
  //fRadioFrameTab1->SetLayoutBroken(kTRUE);
  fRadioFrameTab1->MoveResize(80,130,300,110);

  fRadiosTab1[0]  = new TGRadioButton(fRadioFrameTab1, "TDC channels", 110);
  fRadiosTab1[1]  = new TGRadioButton(fRadioFrameTab1, "CHAMBER 1", 111);
  fRadiosTab1[2]  = new TGRadioButton(fRadioFrameTab1, "CHAMBER 2", 112);
  fRadiosTab1[3]  = new TGRadioButton(fRadioFrameTab1, "SUPERLAYERS", 113);

  for(int i=0; i<4; i++) {
	fRadiosTab1[i]->SetTextJustify(36);
	fRadiosTab1[i]->SetMargins(0,0,0,0);
	fRadiosTab1[i]->SetWrapLength(-1);
        fRadiosTab1[i]->Associate(this);
  }
  fCompositeFrameTab1->AddFrame(fRadioFrameTab1, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fRadiosTab1[0]->SetState(kButtonDown);  //Default
  fRadioFrameTab1->Show();

	//check button
  fCheckButtonFileTab1 = new TGCheckButton(fCompositeFrameTab1, "Dump Histograms to File", 102);
  fCheckButtonFileTab1->MoveResize(80,250,180,20);
  fCheckButtonFileTab1->Associate(this);


  // TAB2: Calibration
  	// container of "Tab2"
  fCompositeFrameTab2 = fTab->AddTab("Calibration");
  fCompositeFrameTab2->SetLayoutManager(new TGVerticalLayout(fCompositeFrameTab2));
  fCompositeFrameTab2->SetLayoutBroken(kTRUE);

	//button 1
  fTextButtonTab2n1 = new TGTextButton(fCompositeFrameTab2,"Load t0s", 200);
  fTextButtonTab2n1->SetTextJustify(36);
  fTextButtonTab2n1->SetMargins(0,0,0,0);
  fTextButtonTab2n1->SetWrapLength(-1);
  fCompositeFrameTab2->AddFrame(fTextButtonTab2n1, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fTextButtonTab2n1->MoveResize(80,30,90,50);
  fTextButtonTab2n1->Associate(this);

	//button 2
  fTextButtonTab2n2 = new TGTextButton(fCompositeFrameTab2,"Load tTrigs", 201);
  fTextButtonTab2n2->SetTextJustify(36);
  fTextButtonTab2n2->SetMargins(0,0,0,0);
  fTextButtonTab2n2->SetWrapLength(-1);
  fCompositeFrameTab2->AddFrame(fTextButtonTab2n2, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fTextButtonTab2n2->MoveResize(190,30,90,50);
  fTextButtonTab2n2->Associate(this);

	//button 3
  fTextButtonTab2n3 = new TGTextButton(fCompositeFrameTab2,"Compute tTrigs", 202);
  fTextButtonTab2n3->SetTextJustify(36);
  fTextButtonTab2n3->SetMargins(0,0,0,0);
  fTextButtonTab2n3->SetWrapLength(-1);
  fCompositeFrameTab2->AddFrame(fTextButtonTab2n3, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fTextButtonTab2n3->MoveResize(290,30,90,50);
  fTextButtonTab2n3->Associate(this);

  	//radio buttons
  fRadioFrameTab2 = new TGButtonGroup(fCompositeFrameTab2, "Histograms: time boxes", kVerticalFrame);
  //fRadioFrameTab1->SetLayoutBroken(kTRUE);
  fRadioFrameTab2->MoveResize(80,130,300,110);

  fRadiosTab2[0]  = new TGRadioButton(fRadioFrameTab2, "per ROBs", 210);
  fRadiosTab2[1]  = new TGRadioButton(fRadioFrameTab2, "per SLs", 211);

  for(int i=0; i<2; i++) {
	fRadiosTab2[i]->SetTextJustify(36);
	fRadiosTab2[i]->SetMargins(0,0,0,0);
	fRadiosTab2[i]->SetWrapLength(-1);
	fRadiosTab2[i]->Associate(this);
  }
  fRadioFrameTab2->Show();
  fCompositeFrameTab2->AddFrame(fRadioFrameTab2, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fRadiosTab2[0]->SetState(kButtonDown);  //Default

	//check button 1
  fCheckButtonTab2n1 = new TGCheckButton(fCompositeFrameTab2, "Dump tTrigs value to File", 221);
  fCheckButtonTab2n1->MoveResize(80,220,200,20);
  fCheckButtonTab2n1->Associate(this);

	//check button 2
  fCheckButtonFileTab2 = new TGCheckButton(fCompositeFrameTab2, "Dump Histograms to File", 222);
  fCheckButtonFileTab2->MoveResize(80,250,200,20);
  fCheckButtonFileTab2->Associate(this);


  // TAB3
  	// container of "Tab3"
  fCompositeFrameTab3 = fTab->AddTab("Muon Tracks");
  fCompositeFrameTab3->SetLayoutManager(new TGVerticalLayout(fCompositeFrameTab3));
  fCompositeFrameTab3->SetLayoutBroken(kTRUE);

	//button
  fTextButtonTab3 = new TGTextButton(fCompositeFrameTab3,"Start Pattern Recognition", 300);
  fTextButtonTab3->SetTextJustify(36);
  fTextButtonTab3->SetMargins(0,0,0,0);
  fTextButtonTab3->SetWrapLength(-1);
  fCompositeFrameTab3->AddFrame(fTextButtonTab3, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fTextButtonTab3->MoveResize(80,30,300,50);
  fTextButtonTab3->Associate(this);

	//check button
  fCheckButtonTab3 = new TGCheckButton(fCompositeFrameTab3, "Fill Histograms", 301);
  fCheckButtonTab3->MoveResize(80,100,150,20);
  fCheckButtonTab3->Associate(this);

  	//radio buttons
  fRadioFrameTab3 = new TGButtonGroup(fCompositeFrameTab3, "Histograms: slope & position", kVerticalFrame);
  //fRadioFrameTab1->SetLayoutBroken(kTRUE);
  fRadioFrameTab3->MoveResize(80,130,300,110);

  fRadiosTab3[0]  = new TGRadioButton(fRadioFrameTab3, "CHAMBER 1", 310);
  fRadiosTab3[1]  = new TGRadioButton(fRadioFrameTab3, "CHAMBER 2", 311);
  fRadiosTab3[2]  = new TGRadioButton(fRadioFrameTab3, "SUPERLAYERS", 312);

  for(int i=0; i<3; i++) {
	fRadiosTab3[i]->SetTextJustify(36);
	fRadiosTab3[i]->SetMargins(0,0,0,0);
	fRadiosTab3[i]->SetWrapLength(-1);
	fRadiosTab3[i]->Associate(this);
  }
  fRadioFrameTab3->Show();
  fCompositeFrameTab3->AddFrame(fRadioFrameTab3, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fRadiosTab3[0]->SetState(kButtonDown);  //Default

	//check button
  fCheckButtonFileTab3 = new TGCheckButton(fCompositeFrameTab3, "Dump Histograms to File", 302);
  fCheckButtonFileTab3->MoveResize(80,250,180,20);
  fCheckButtonFileTab3->Associate(this);

  // TAB4
  	// container of "Tab4"
  fCompositeFrameTab4 = fTab->AddTab("Event Display");
  fCompositeFrameTab4->SetLayoutManager(new TGVerticalLayout(fCompositeFrameTab4));
  fCompositeFrameTab4->SetLayoutBroken(kTRUE);

	//button
  fTextButtonTab4 = new TGTextButton(fCompositeFrameTab4,"Start Event Display", 40000);
  fTextButtonTab4->SetTextJustify(36);
  fTextButtonTab4->SetMargins(0,0,0,0);
  fTextButtonTab4->SetWrapLength(-1);
  fCompositeFrameTab4->AddFrame(fTextButtonTab4, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fTextButtonTab4->MoveResize(80,30,300,50);
  fTextButtonTab4->Associate(this);

  // TAB5
  	// container of "Tab5"
  fCompositeFrameTab5 = fTab->AddTab("Image");
  fCompositeFrameTab5->SetLayoutManager(new TGVerticalLayout(fCompositeFrameTab5));
  fCompositeFrameTab5->SetLayoutBroken(kTRUE);

	//button
  fTextButtonTab5 = new TGTextButton(fCompositeFrameTab5,"Start Image Reconstruction", 500);
  fTextButtonTab5->SetTextJustify(36);
  fTextButtonTab5->SetMargins(0,0,0,0);
  fTextButtonTab5->SetWrapLength(-1);
  fCompositeFrameTab5->AddFrame(fTextButtonTab5, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fTextButtonTab5->MoveResize(80,30,300,50);
  fTextButtonTab5->Associate(this);

	//check button
  fCheckButtonTab5 = new TGCheckButton(fCompositeFrameTab5, "Fill Histograms", 501);
  fCheckButtonTab5->MoveResize(80,100,150,20);
  fCheckButtonTab5->Associate(this);

	//check button
  fCheckButtonFileTab5 = new TGCheckButton(fCompositeFrameTab5, "Dump Image File", 502);
  fCheckButtonFileTab5->MoveResize(80,250,180,20);
  fCheckButtonFileTab5->Associate(this);

  // Standard Operations
  fTab->SetTab(0);
  fTab->Move(0,50);
  fTab->Resize(fTab->GetDefaultSize());

  AddFrame(fTab, new TGLayoutHints(kLHintsBottom | kLHintsExpandX | kLHintsExpandY, 2, 2, 5, 1));

  //check status
  int retval;
  int status = GetStatus();
  if (status==-2) {
     new  TGMsgBox(fClient->GetRoot(), this,
         "Error in file!", "Access Error -- Reading Stopped.",
         kMBIconStop, kMBOk, &retval);
  }

   // sets window name and show the main frame
   // for style
   SetWindowName("Muon TOMography tool");
   SetMWMHints(kMWMDecorAll,kMWMFuncAll,kMWMInputModeless);
   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();

};

TOMMainFrame::~TOMMainFrame()
{
  //Destructor gets rid of widgets
  delete fCloseButton;
  delete fButtonFrame;
  delete fMenuBar;
  delete fMenuFile;
  delete fTab;

  delete fCompositeFrameTab1;
  delete fTextButtonTab1;
  delete fIconTab1;
  delete fLabelTab1;
  delete fRadioFrameTab1;
  delete [] fRadiosTab1;
  delete fCheckButtonTab1;
  delete fCheckButtonFileTab1;

  delete fCompositeFrameTab2;  
  delete fTextButtonTab2n1;
  delete fTextButtonTab2n2;
  delete fTextButtonTab2n3;
  delete fIconTab2;
  delete fLabelTab2;
  delete fRadioFrameTab2;
  delete [] fRadiosTab2;
  delete fCheckButtonTab2n1;
  delete fCheckButtonFileTab2;

  delete fCompositeFrameTab3;
  delete fTextButtonTab3;
  delete fRadioFrameTab3;
  delete [] fRadiosTab3;
  delete fCheckButtonTab3;
  delete fCheckButtonFileTab3;

  delete fCompositeFrameTab4;
  delete fTextButtonTab4;
  delete fCheckButtonTab4;

  delete fCompositeFrameTab5;
  delete fTextButtonTab5;
  delete fCheckButtonTab5;
  delete fCheckButtonFileTab5;

  delete fC1;
  delete _rawAnalyzer;
}

void TOMMainFrame::CloseWindow()
{
  //Standard close window routine
  TGMainFrame::CloseWindow();
  if(fCheckButtonFileTab1->GetState() == kButtonDown)
    DumpRawHistoFile();
  gApplication->Terminate(0);
}

int TOMMainFrame::GetStatus()
{
  return 1;
}

int TOMMainFrame::OpenRawFile(int mode, char *filename)
{
  if (mode==0) {
    TGFileInfo fi;
    fi.fFileTypes = (const char **)filetypes;
    fi.fIniDir    = StrDup("/data/radmu/flat/"); 
    new TGFileDialog(fClient->GetRoot(), this, kFDOpen,&fi);

    if(DEBUG_TOM)
      cout << "You selected the file " << fi.fFilename << endl;

    if (fi.fFilename == 0) return -1;  //User hit 'Cancel'
    if (strstr(fi.fFilename, ".i")==NULL) {
      // NOT A CONFIGURATION FILE
      return -2;
    }
    strcpy(rawFileName, fi.fFilename);
  }
  else strcpy(rawFileName, filename);


  ifstream inpFile( rawFileName, ios::in );
  if (!inpFile) { //i.e., file does not exist
    return -3;
  }
  return 1;
}

int TOMMainFrame::CloseRawFile()
{
  strcpy(rawFileName, "");
  if(DEBUG_TOM)
	cout << "Raw File : " << rawFileName << endl;
  return 10;
} 

int TOMMainFrame::ReadFile()
{
  int retval;

  if(strcmp(rawFileName,"")==0){
     new  TGMsgBox(fClient->GetRoot(), this,
         "Error in Reading!", "Read Error -- No File Specified!\n Select from File Menu - Open Raw File",
         kMBIconStop, kMBOk, &retval);
  }
  else
     _rawAnalyzer->goAnalysis(rawFileName,40000);

  if(DEBUG_TOM)
	cout << "END Reading Raw File : " << rawFileName << endl;
  
  return 1;
}

int TOMMainFrame::FillRawHistos(bool flagDown){
  
  if(DEBUG_TOM) 
    cout << "Fill raw histos" << endl;
  
  if(flagDown){
    //Set up the Drawing Canvas & Pad Configuration
    fC1 = new TCanvas("ReadbackDisplayProgram", "Tomography Display", 5, 5, 800, 900);
    
    if(!_rawHistos)
      _rawHistos = new RawHistos();
    _rawHistos->setFlagFillHistos(true);
    // pass RawHistos pointer to RawAnalyzer
    _rawAnalyzer->setRawHistosPtr(_rawHistos);
  }
  else{
    if(fC1)
      delete fC1;
    _rawHistos->setFlagFillHistos(false);
    delete _rawHistos;
    _rawAnalyzer->setRawHistosPtr(0);
  }
  return 1;
}

int TOMMainFrame::ComputeTTrig(){
  
  if(DEBUG_TOM) 
    cout << "Fill histos to fit" << endl;
  
  //Set up the Drawing Canvas & Pad Configuration
  if(!fC1)
    fC1 = new TCanvas("ReadbackDisplayProgram", "Tomography Display", 5, 5, 800, 900);
  if(!_ttrigCalib)
    _ttrigCalib = new TTrigCalibration();
  _ttrigCalib->setFlagFillHistos(true);
  // pass CalibHistos pointer to CalibAnalyzer
  _rawAnalyzer->setTTrigCalibPtr(_ttrigCalib);
  
  return 1;
}

int TOMMainFrame::DumpTTrigs(){
  char  calibFileName[128];
  strcpy(calibFileName,"ttrigs.txt"); 

/*  if(strcmp(rawFileHistoName,"")==0){
     new  TGMsgBox(fClient->GetRoot(), this,
         "Error in Reading!", "Read Error -- No File Specified!\n Select from File Menu - Open Raw File",
         kMBIconStop, kMBOk, &retval);
  }
*/
  _ttrigCalib->dumpTTrigs(calibFileName);

  if(DEBUG_TOM)
	cout << "Dumping Calibration TTrig File : " << calibFileName << endl;
  
  return 1;
}

int TOMMainFrame::DumpCalibHistoFile()
{
  char  calibHistoFileName[128];
  strcpy(calibHistoFileName,"calibHistos.root"); 

/*  if(strcmp(rawFileHistoName,"")==0){
     new  TGMsgBox(fClient->GetRoot(), this,
         "Error in Reading!", "Read Error -- No File Specified!\n Select from File Menu - Open Raw File",
         kMBIconStop, kMBOk, &retval);
  }
*/
  _ttrigCalib->dumpHistos(calibHistoFileName);

  if(DEBUG_TOM)
	cout << "Dumping Calibration Histo File : " << calibHistoFileName << endl;
  
  return 1;
}

int TOMMainFrame::DumpRawHistoFile()
{
  char  rawHistoFileName[128];
  strcpy(rawHistoFileName,"rawHistos.root"); 

/*  if(strcmp(rawFileHistoName,"")==0){
     new  TGMsgBox(fClient->GetRoot(), this,
         "Error in Reading!", "Read Error -- No File Specified!\n Select from File Menu - Open Raw File",
         kMBIconStop, kMBOk, &retval);
  }
*/
  _rawHistos->dumpHistos(rawHistoFileName);

  if(DEBUG_TOM)
	cout << "Dumping Raw Histo File : " << rawHistoFileName << endl;
  
  return 1;
}

Bool_t TOMMainFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
  //Function that defines what to do when window widgets are activated
  int status = 0;

  //   BELOW IS A USEFUL DEBUGGING STATMENT WHEN WIDGETS DON'T WORK CORRECTLY
  //   MOST COMMON ERROR -- failing to call 'Associate(this)' for that widget

  if(DEBUG_TOM){
    cout << "Msg = " << msg << "    GET_MSG = " << GET_MSG(msg) 
       << "    SUB_MSG = " << GET_SUBMSG(msg) << "  parm1 = "
       << parm1 << "  parm2 = " << parm2 << endl;
  } 

  switch (GET_MSG(msg)) {
    case kC_COMMAND:  //command type event

      switch (GET_SUBMSG(msg)) {

        case kCM_MENU:
	  switch (parm1) {
  	    case M_FILE_OPEN:    status = OpenRawFile(0, "");  break;
            case M_FILE_CLOSE:   status = CloseRawFile();  break;
	    case M_QUIT:   CloseWindow();  break;
	  }
	break;

        case kCM_BUTTON:  	// button event
	  switch (parm1) {
	    case B_QUIT:   CloseWindow(); break;
            case 100: ReadFile(); break; 
	    //case 200: 
	    //case 201:	
	    case 202: ComputeTTrig(); break;
            default: break;
          }
        break;

        case kCM_CHECKBUTTON:	// check button
	  if(parm1==101 && fCheckButtonTab1->GetState() == kButtonDown) 
		status = HandleCheck(parm1,true);
          if(parm1==101 && fCheckButtonTab1->GetState() == kButtonUp) 
		status = HandleCheck(parm1,false);
	  if(parm1==221 && fCheckButtonTab2n1->GetState() == kButtonDown) 
		status = HandleCheck(parm1,true);
	  if(parm1==222 && fCheckButtonFileTab2->GetState() == kButtonDown) 
		status = HandleCheck(parm1,true);

	  break;
        case kCM_RADIOBUTTON: //radio button event
	  status = HandleRadio(parm1);
	  break;
        default:      break;
      }
    break;
    default:      break;
  }

  int retval;

  if (status==-3) {
    new  TGMsgBox(fClient->GetRoot(), this,
         "Error in Opening!", "Opening Error -- File Does Not Exist!",
         kMBIconStop, kMBOk, &retval);
  }
  else if (status==-2) {
    new  TGMsgBox(fClient->GetRoot(), this,
         "Error in Opening!", "Opening Error -- Not a flat file!",
         kMBIconStop, kMBOk, &retval);
  }
  else if (status==-1) {
     new  TGMsgBox(fClient->GetRoot(), this,
         "Error in Opening!", "Opening Error -- No File Specified!",
         kMBIconStop, kMBOk, &retval);
  }
  else if (status==1) {
     new  TGMsgBox(fClient->GetRoot(), this,
         "Raw File Selected!", rawFileName,
         kMBIconAsterisk, kMBOk, &retval);
  }

  return kTRUE;

}

int TOMMainFrame::HandleRadio(Long_t parm1)
{
  switch (parm1) {
    case 110: 	// TDC histos
        if(DEBUG_TOM)
	  cout << "HandleRadio 110: ..." << endl;
	fC1->Clear();
	fC1->cd();
	_rawHistos->buildCanvasTDC(fC1);
	_rawHistos->setCanvas(0);
        _rawHistos->updateCanvas(0);
        fC1->Update();
        break;
    case 111: 	// CH1
	fC1->Clear();
	fC1->cd();
	_rawHistos->buildCanvasCH1(fC1);
        _rawHistos->setCanvas(1);
        _rawHistos->updateCanvas(1);
        fC1->Update();
	break;
    case 112: 	// CH2
	fC1->Clear();
	fC1->cd();
	_rawHistos->buildCanvasCH2(fC1);
        _rawHistos->setCanvas(2);
        _rawHistos->updateCanvas(2);
        fC1->Update();
	break;
    case 113: 	// SLs
	fC1->Clear();
	fC1->cd();
	_rawHistos->buildCanvasSLS(fC1);
        _rawHistos->setCanvas(3);
        _rawHistos->updateCanvas(3);
        fC1->Update();
   	break;
    case 210: 	// ROS
	fC1->Clear();
	fC1->cd();
	_ttrigCalib->buildCanvasTROB(fC1);
        _ttrigCalib->setCanvas(0);
        _ttrigCalib->computeTTrig(0);
        fC1->Update();
	break;
    case 211: 	// SLs
	fC1->Clear();
	fC1->cd();
	_ttrigCalib->buildCanvasTSL(fC1);
        _ttrigCalib->setCanvas(1);
        _ttrigCalib->computeTTrig(1);
        fC1->Update();
   	break;

  default: break;
  }
  return 110;
}

int TOMMainFrame::HandleCheck(Long_t parm1, bool flagDown)
{
  switch (parm1) {
  case 101: 
	FillRawHistos(flagDown);	
	break;
  case 221:
        if(flagDown)
	  DumpTTrigs();
  case 222:
        if(flagDown)
          DumpCalibHistoFile();
	
  default:   break;
  }


  return 100;
}


int TOMMainFrame::FillGraph(Int_t padindex, const char *chname )
{
 return 1;
}

int TOMMainFrame::Fill2DHistogram(Int_t padindex, const char *chname)
{
  return 1;
}

int TOMMainFrame::FillHistogram(Int_t padindex, const char *chname)
{

  return 1;
}


