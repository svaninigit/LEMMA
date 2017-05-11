// 2008/01/28
// Sara Vanini & Silvia Pesente: 
// histograms for muon radiography data acquisition monitor
// calibration tools 

#ifndef CALIBHISTOS_h
#define CALIBHISTOS_h

#include <iostream.h>
#include <stdio.h>
#include <iomanip.h>
#include <map.h>
#include <fstream.h>
#include <algorithm>

// ROOT utilities
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGCanvas.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TFrame.h"
#include "TString.h"
#include "TFile.h"

// TOM classes
#include "../PattRec/HITCollection.h"


class CalibHistos {

public:
  CalibHistos();
  ~CalibHistos();

  // member functions
  void initHistos();
  void fillHistos(HITCollection * hits);
  void dumpTTrigs(char * name);
  void dumpHistos(char * name);
  void computeTTrigs();

  // canvas
  void buildCanvasROS(TCanvas * c=0);
  void buildCanvasSLS(TCanvas * c=0);
  void updateCanvas(int flag);

  // flag for histos filling and CalibHisto pointer
  inline void setFlagFillHistos(bool flag) { _flagFillHistos = flag; return; };
  inline bool flagFillHistos() { return _flagFillHistos; }
  inline bool flagCanvasROS() { return _flagCanvasROS; }
  inline bool flagCanvasSLS() { return _flagCanvasSLS; }

  inline void setCanvas(int flag) { _nCanvas = flag; return; };
  inline int getCanvas() { return _nCanvas; }
 
private:
  //for debugging
  bool DEBUG_FLAG;

  //utilities
  bool _flagFillHistos;
  bool _flagCanvasROS, _flagCanvasSLS;
  int _nCanvas;

  //occupancy histo cuts
  int minTime;
  int maxTime;

  //graphic tools
  TCanvas 	*c1, *c2;
  TPad 		*pad1, *pad2;

  //histos
};
#endif /*CALIBHISTOS_h*/

