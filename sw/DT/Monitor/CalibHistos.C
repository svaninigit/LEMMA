#include "CalibHistos.h"
#include "../PattRec/HITCollection.h"
#include "../PattRec/HIT.h"


CalibHistos::CalibHistos() {

  // init useful variables
  //for debugging
  DEBUG_FLAG = 1;

  //occupancy histo cuts --> FIX tdc units or nseconds?
  minTime = 1200;
  maxTime = 2200;

  // init histograms
  initHistos();

  // reset flags
  _flagCanvasROS=0;
  _flagCanvasSLS=0;

  return;
};

CalibHistos::~CalibHistos() {

   if(_flagCanvasROS){
    delete pad1;
  }

  if(_flagCanvasSLS){
    delete pad2;
  }

  if(c1)
    delete c1;

  if(c2)
    delete c2;

  return;
};

void CalibHistos::initHistos(){
  
  return;
};

void CalibHistos::buildCanvasROS(TCanvas * c1) {

  if(DEBUG_FLAG)
    cout << "CalibHistos::buildCanvasROS" << endl; 
  //build the MINICRATE canvas with pads
  if(c1==0)
    c1 = new TCanvas("c1","RADMU monitor: ROS",1,10,800,800);
  c1->SetFillColor(10);
  //gStyle->SetFrameFillColor(18); 
  c1->SetBorderSize(2);
  TPaveText * title = new TPaveText(.2,0.96,.8,.995);
  title->SetFillColor(33);
  title->AddText("RADMU MONITOR: ROS time boxes");
  title->Draw();
  pad1 = new TPad("pad1","The pad with ROS time boxes",0.05,0.50,0.95,0.92,19,1,1);
  pad1->Draw();

  // flag
  _flagCanvasROS = true;

  return;
}

void CalibHistos::buildCanvasSLS(TCanvas * c2) {

  //SLS canvas
  if(c2==0)
    c2 = new TCanvas("c2","RADMU monitor: SLS time boxes",160,10,800,800);
  c2->SetFillColor(10);
  c2->SetBorderSize(1);
  
  TPaveText * title_c2 = new TPaveText(.2,0.96,.8,.995);
  title_c2->SetFillColor(33);
  title_c2->AddText("RADMU MONITOR: SLS time boxes");
  title_c2->Draw();

  pad2 = new TPad("pad2","The pad with SLS time boxes",0.01,0.01,0.5,0.95,19,1,1);
  pad2->Divide(1,12);  
  pad2->Draw();

  // flag
  _flagCanvasSLS = true;

  return;
}


void CalibHistos::fillHistos(HITCollection * hits){

  if(DEBUG_FLAG)
    cout << "CalibHistos::fillHistos : number of hits= " << hits->Get_NumberHITS() << endl;

  // loop on hits
  for(int i=0; i<hits->Get_NumberHITS(); i++) {

	// get hit and fill histos
   	HIT *hit = hits->hit(i);
	int se = hit->CH_ID();
        int sl = hit->SL_ID();
	int lay = hit->L_ID();
	int wire = hit->wire_ID();
	int TDC_Id = hit->TDC_ID();
	int ROB_Id = hit->ROB_ID();
	int channel = hit->channel_ID();

	float time = hit->rtime();
	//float dtime = hit->dtime();


	//fill tdc histograms
	// ROB: 0..24;  TDC:0..3;  ch:0..31
	int tdcflag = TDC_Id + 4*ROB_Id;
	int chflag = channel + 32*tdcflag;
				
	if(DEBUG_FLAG)
		cout << "Filling : lay " << lay << "  sl " << sl << 
			" tube " << wire << endl;				

	if(DEBUG_FLAG)
		cout 	<< " Leading measurment " << " TDC " << setbase(10) << TDC_Id 
			<< " ch " << channel << " time " << time << endl;

  } //end hit loop

  return;

}

void CalibHistos::updateCanvas(int flag) {
  // NB flag: 0=ROS, 1=SLS

  if(DEBUG_FLAG)
    cout << "CalibHistos::updateCanvas() " << ", _flagCanvasROS " << _flagCanvasROS
		<< ", _flagCanvasSLS " << _flagCanvasSLS << endl; 

  //ROS
  if(flag==0){
    pad1->cd();
 
    pad1->Update();
  }

  //SLS
  if(flag==1){
    pad2->cd();

  }

  if(DEBUG_FLAG)
    cout << "CalibHistos::updateCanvas() : END" << endl;

  return;
}; 

void CalibHistos::dumpHistos(char* fileName){

  if(DEBUG_FLAG)
    cout << "CalibHistos::dumpHistos() " << endl;  

  TFile * file = new TFile(fileName,"RECREATE");
  file->cd();


  return;
};

void CalibHistos::dumpTTrigs(char * name){

  if(DEBUG_FLAG)
    cout << "CalibHistos::dumpTTrigs() " << endl;

  return;
}

void CalibHistos::computeTTrigs(){

  if(DEBUG_FLAG)
    cout << "CalibHistos::computeTTrigs() " << endl;

  return;
}



