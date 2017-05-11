// 2008/01/28
// Sara Vanini : simple root monitor for muon radiography data acquisition


#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <map>
#include <fstream>
#include <string>

#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGCanvas.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TFrame.h"
#include "TFrame.h"
#include "TROOT.h" 

//occupancy histo cuts
//int minTime = 1200;
//int maxTime = 2200;
int minTime = 0;
int maxTime = 22000;

//define map array and function to fill it
int nros =  1;
int nrob =  19;
int ntdc =  4;
int ncha = 32;
map<int, int> chmap;
void fillMap();
int getTDCid( int ros, int rob, int tdc, int cha );
int getTubeId( int se, int sl, int lay, int tube );
int getTube(int idTube);
int getLay(int idTube);
int getSL(int idTube);
int getSe(int idTube);  

//for debugging
bool DEBUG_FLAG = 0;

//graphic tools
TCanvas * c1;
TCanvas * ch1;
TCanvas * ch2;
TCanvas * sls;
TPad * pad1;
TPad * pad2;
TPad * pad3;
TPad * pad_hits_ch1;
TPad * pad_time_ch1;
TPad * pad_hits_ch2;
TPad * pad_time_ch2;
TPad * pad_hits_sls;
TPad * pad_time_sls;

//histos
TH1F * hocc_lay[32]; 

void radmuMonitor(char *filename)
//void radmuMonitor()
{
  //reset 
  gROOT->Reset();

  //fill channel map
  fillMap();

  /*
  //map test
  int rawdata = getTDCid(1,12,0,4);
  int chdata = getTubeId(10,3,2,2);
  map<int,int>::iterator iter = chmap.find(rawdata);
  if( iter != chmap.end() ) {
    	cout << "Key : " << iter->first << " entry : " << iter->second << endl;
	cout << "se :" << getSe(iter->second) <<
		  ", SL :" << getSL(iter->second) <<
		  ", Lay :" << getLay(iter->second) <<
		  ", Tube :" << getTube(iter->second) << endl;
  }
*/

  //build the MINICRATE canvas with pads
  c1 = new TCanvas("c1","RADMU monitor: TDC",100,10,800,800);
  c1->SetFillColor(10);
  //gStyle->SetFrameFillColor(18); 
  c1->SetBorderSize(2);
  TPaveText * title = new TPaveText(.2,0.96,.8,.995);
  title->SetFillColor(33);
  title->AddText("RADMU MONITOR: TDC hits");
  title->Draw();
  pad1 = new TPad("pad1","The pad with hits",0.05,0.50,0.95,0.92,19,1,1);
  pad1->Draw();
  pad2 = new TPad("pad2","The pad with times",0.05,0.05,0.50,0.45,19,1,1);
  pad2->Draw();
  pad3 = new TPad("pad3","The pad with times",0.50,0.05,0.95,0.45,19,1,1);
  pad3->Draw();
 // pad1->cd();

  //Chamber 2 canvas
  ch2 = new TCanvas("ch2","RADMU monitor: CHAMBER 2",100,10,800,800);
  ch2->SetFillColor(10);
  ch2->SetBorderSize(1);
  
  TPaveText * title_ch2 = new TPaveText(.2,0.96,.8,.995);
  title_ch2->SetFillColor(33);
  title_ch2->AddText("RADMU MONITOR: chamber 2");
  title_ch2->Draw();

  pad_hits_ch2 = new TPad("pad_hits_ch2","The pad with hits",0.01,0.01,0.5,0.95,19,1,1);
  pad_hits_ch2->Divide(1,12);  
  pad_hits_ch2->Draw();

  pad_time_ch2 = new TPad("pad_time_ch2","The pad with times",0.5,0.01,0.99,0.95,19,1,1);
  pad_time_ch2->Divide(1,3);
  pad_time_ch2->Draw();

  //Chamber 1 canvas
  ch1 = new TCanvas("ch1","RADMU monitor: CHAMBER 1",100,10,800,800);
  ch1->SetFillColor(10);
  ch1->SetBorderSize(1);
 
  TPaveText * title_ch1 = new TPaveText(.2,0.96,.8,.995);
  title_ch1->SetFillColor(33);
  title_ch1->AddText("RADMU MONITOR: chamber 1");
  title_ch1->Draw();
  
  pad_hits_ch1 = new TPad("pad_hits_ch1","The pad with hits",0.01,0.01,0.5,0.95,19,1,1);
  pad_hits_ch1->Divide(1,12);  
  pad_hits_ch1->Draw();

  pad_time_ch1 = new TPad("pad_time_ch1","The pad with times",0.5,0.01,0.99,0.95,19,1,1);
  pad_time_ch1->Divide(1,3);
  pad_time_ch1->Draw();


  //SLs canvas
  sls = new TCanvas("sls","RADMU monitor:SLs",100,10,800,800);
  sls->SetFillColor(10);
  sls->SetBorderSize(1);
 
  TPaveText * title_sls = new TPaveText(.2,0.96,.8,.995);
  title_sls->SetFillColor(33);
  title_sls->AddText("RADMU MONITOR: SL2");
  title_sls->Draw();
  
  pad_hits_sls = new TPad("pad_hits_sls","The pad with hits",0.01,0.01,0.5,0.95,19,1,1);
  pad_hits_sls->Divide(1,12);  
  pad_hits_sls->Draw();

  pad_time_sls = new TPad("pad_time_sls","The pad with times",0.5,0.01,0.99,0.95,19,1,1);
  pad_time_sls->Divide(1,3);
  pad_time_sls->Draw();

  //occupancy histograms: all layers from bottom to top 4+4+12+12
  int nCell = 72;
  for(int i=0; i<32; i++)
  {
  	string nameLayer( "hocc" );
	int layNum = i+1;	
	if(i<4){
		nameLayer += "_sl1_lay_";
		nameLayer += layNum + '0';
	}
	if(i>=4 && i<8){
		nameLayer += "_sl2_lay_";
		nameLayer += layNum - 4  + '0';
	}
	if(i>=8 && i<20){
		nameLayer += "_ch1_lay_";
		if(layNum<18)
			nameLayer += layNum - 8 + '0';
		else{
			nameLayer += 1 + '0';
			nameLayer += layNum - 18  + '0';
		}		
	}
	if(i>=20){
		nameLayer += "_ch2_lay_";
		if(layNum<30)
			nameLayer += layNum - 20 + '0';
		else{
			nameLayer += 1 + '0';
			nameLayer += layNum - 30  + '0';
		}				
	}
	
	hocc_lay[   i] = new TH1F( nameLayer.c_str(),
                             nameLayer.c_str(), nCell, 0.5, nCell + 0.5 );
  }

   //TDC histograms
  TH2F * hocc  = new TH2F("hocc","Channel vs TDC",200,0,99,66,-1,32);
  TH1F * htime = new TH1F("htime","TimeBox",2000,1000,3000); 
  TH2F * htimech = new TH2F("htimech","Channel vs time",2000,1000,3000,8000,0,4000); 

  //time boxes in SLs histograms from bottom to top 1+1+3+3
  int histoWidth = 2000;
  int histoEdgeL = 1000;
  int histoEdgeH = 3000;

  TH1F * htime_sl[8];
  for(int j; j<8; j++)
  {
	string name( "htime" );
  	name += "_sl_";
  	name += (j+1) + '0';

    	htime_sl[  j] = new TH1F( name.c_str(),
                                name.c_str(),
                                histoWidth, histoEdgeL, histoEdgeH );
  }


  //FILE is a C-type file pointer included in stdio.h
  FILE *infile;
  infile = fopen(filename,"r");
  if(infile==0)
    cout << "ERROR opening file..." << endl;

  //infile = fopen("r395.i1","r");
  //long int is 32-bit=4-byte word
  long evbu; 
  //NB size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream); 
  //Reads data from the given stream into the array pointed to by ptr. 
  //It reads nmemb number of elements of size size (in byte).
  //The total number of bytes read is (size*nmemb). 
  //On success the number of elements read is returned. 
  //On error or end-of-file the total number of elements successfully 
  //read (which may be zero) is returned. 

  //NB ftell(infile) returns the current file position. For binary stream, 
  //then the value is the number of bytes from the beginning of the file.  

  for(int words=1; words<1000000; words++)
  {	
	// get the 32-bit word
    	fread(&evbu,4,1,infile);

      // DMA: swap 16 lsb with 16 msb
      long word = (evbu>>16) | (evbu<<16);

      if(DEBUG_FLAG)
        cout << "Reading word " << word << endl;

	// type of packet and TTCcounts
	long type = word & 0xFF000000;
	type = type >> 24;			// type of data packet: bit 24-31
	long TTCcount = word & 0x00FFFFFF;	// TTC count: bit 0-23

	if(DEBUG_FLAG)
	{
	if(type == 0x1F)
		cout << " ---> Event header ! " << endl;
              
	if(type == 0x3F)
		cout << " ---> Event trailer ! " << endl;
              
	if(type == 0xDF)
		cout << " ---> Error Flag ! " << endl;

	if(type == 0xFF)
		cout << " ---> Debugging data ! " << endl;
	}

	//variable declaration: the value is filled at the first occurance
	long Event_Id, Bunch_Id, TDC_Id, channel, time, ROB_Id;
	map<int,int>::iterator iter;

	if(type != 0x1F && type != 0x3F && type != 0xDF && type != 0xFF )
	{
		long group = type & 0xE0;
		group = group >> 5;
		
		//variable declaration
		long Event_Id, Bunch_Id, TDC_Id, channel, time, ROB_Id;

		switch(group)
		{
		case 0:
			Event_Id = TTCcount & 0xFFF000; 
			Event_Id  = Event_Id >> 12;
			Bunch_Id = TTCcount & 0xFFF;
			ROB_Id = type & 0x1F;

			if(DEBUG_FLAG)
			{
				cout << "Group header: ";
				cout << " Event " << setbase(10) << Event_Id << 
					" bunch " << Bunch_Id << 
					" rob " << ROB_Id << endl;	
			}
			break;
		case 1:
			Event_Id = TTCcount & 0xFFF000; 
			Event_Id  = Event_Id >> 12;
			Bunch_Id = TTCcount & 0xFFF;
			if(DEBUG_FLAG)
			{
				cout << "Group trailer: ";
				cout << " Event " << setbase(10) << Event_Id << 
					" word count " << Bunch_Id << endl;	
			}
			break;
		case 2:
			Event_Id = TTCcount & 0xFFF000; 
			Event_Id  = Event_Id >> 12;
			Bunch_Id = TTCcount & 0xFFF;
			TDC_Id = ROB_Id & 0x3;	
			if(DEBUG_FLAG)
			{
				cout << "TDC Id " << setbase(10) << TDC_Id << " header: ";
				cout << " Event " << setbase(10) << Event_Id << 
					" bunch " << Bunch_Id << endl;	
			}
			break;
		case 3:
			Event_Id = TTCcount & 0xFFF000; 
			Event_Id  = Event_Id >> 12;
			Bunch_Id = TTCcount & 0xFFF;
			TDC_Id = ROB_Id & 0x3;		
			if(DEBUG_FLAG)
			{
				cout << "TDC Id " << setbase(10) << TDC_Id << " trailer: ";
				cout << " Event " << setbase(10) << Event_Id << 
					" word count " << Bunch_Id << endl;	
			}
			break;
		case 4:				
			TDC_Id = type & 0x3;	
			channel = TTCcount & 0xF80000;
			channel = channel >> 19;
			time = TTCcount & 0x7FFFF;

			htime->Fill(time/4.);
                  {
			// ROB: 0..24;  TDC:0..3;  ch:0..31
			int tdcflag = TDC_Id + 4*ROB_Id;
			int chflag = channel + 32*tdcflag;

			if(DEBUG_FLAG)
				cout <<	"tdc " << TDC_Id << " rob " << ROB_Id << 
					" ch " << channel << " tdcflag " << tdcflag << 
					" chflag " << chflag << endl;
		
			//fill tdc histograms
            //if(time/4. > minTime && time/4. < maxTime)
				hocc->Fill(tdcflag,channel);
			htimech->Fill(time/4.,chflag);
                  
                  //fill ch histograms
   			int rawdata = getTDCid(1,ROB_Id,TDC_Id,channel);
			//cout << "TDCid " << rawdata << endl;

			int lay;
			int sl;
		        iter = chmap.find(rawdata);
  			if( iter != chmap.end() )
			{
				if(DEBUG_FLAG)
				{
					cout << "SE " << getSe(iter->second) <<  
					" sl " << getSL(iter->second) << 
					" lay " << getLay(iter->second) << 
					" tube " << getTube(iter->second) << endl;
				}
 			  	if( getSe(iter->second) == 8 || getSe(iter->second) == 9
				||getSe(iter->second) == 10 || getSe(iter->second) == 11 ){
				if(getSe(iter->second) == 8){
					lay = getLay(iter->second) - 1;	//0,1,2,3
					sl = 0;
				}
				if(getSe(iter->second) == 9){
					lay = getLay(iter->second) + 3;	//4,5,6,7
					sl = 1;
				}
				if(getSe(iter->second) == 10){
					lay = getLay(iter->second) + 7 + (getSL(iter->second)-1)*4;	//8,9,..19
					sl = getSL(iter->second) + 1;
				}
				if(getSe(iter->second) == 11){
					lay = getLay(iter->second) + 19 + (getSL(iter->second)-1)*4;	//20,21,..
					sl = getSL(iter->second) + 4;
				}
				
				
				if(DEBUG_FLAG)
					cout << "Filling : lay " << lay << "  sl " << sl << 
						" tube " << getTube(iter->second) << endl;				
                //if(time/4. > minTime && time/4. < maxTime)
					hocc_lay[lay]->Fill(getTube(iter->second));
				htime_sl[sl]->Fill(time/4.);
			  }
			}
                  }
 
			if(DEBUG_FLAG)
				cout 	<< " Leading measurment " << " TDC " << setbase(10) << TDC_Id 
					<< " ch " << channel << " time " << time << endl;
			break;
		case 5:
			if(DEBUG_FLAG)
				cout << " Trailing measurment " << endl;
			break;	
		case 6:
			if(DEBUG_FLAG)
				cout << " Errors " << endl;
			break;
		case 7:	
			if(DEBUG_FLAG)
				cout << " Debugging data " << endl;	
			break;
		}//end switch


	}//end no ROS type

    	//update TDC plots
	if(words/100000. == int(words/100000.))
	{
    		pad1->cd();
    		hocc->Draw();	
    		pad1->Update();
    		pad2->cd();
    		htimech->Draw();
    		pad2->Update(); 
    		pad3->cd();
    		htime->Draw();
    		pad3->Update();

		//chamber 2 plots
		ch2->cd();
		pad_hits_ch2->cd(1);
		hocc_lay[31]->Draw();
		pad_hits_ch2->cd(2);
		hocc_lay[30]->Draw();
		pad_hits_ch2->cd(3);
		hocc_lay[29]->Draw();
		pad_hits_ch2->cd(4);
		hocc_lay[28]->Draw();
		pad_hits_ch2->cd(5);
		hocc_lay[27]->Draw();
		pad_hits_ch2->cd(6);
		hocc_lay[26]->Draw();
		pad_hits_ch2->cd(7);
		hocc_lay[25]->Draw();
		pad_hits_ch2->cd(8);
		hocc_lay[24]->Draw();
		pad_hits_ch2->cd(9);
		hocc_lay[23]->Draw();
		pad_hits_ch2->cd(10);
		hocc_lay[22]->Draw();
		pad_hits_ch2->cd(11);
		hocc_lay[21]->Draw();
		pad_hits_ch2->cd(12);
		hocc_lay[20]->Draw();
		pad_hits_ch2->Update();

		pad_time_ch2->cd(1);
		htime_sl[7]->Draw();
		pad_time_ch2->cd(2);
		htime_sl[6]->Draw();
		pad_time_ch2->cd(3);
		htime_sl[5]->Draw();

		//chamber 1 plots
		ch1->cd();
		pad_hits_ch1->cd(1);
		hocc_lay[19]->Draw();
		pad_hits_ch1->cd(2);
		hocc_lay[18]->Draw();
		pad_hits_ch1->cd(3);
		hocc_lay[17]->Draw();
		pad_hits_ch1->cd(4);
		hocc_lay[16]->Draw();
		pad_hits_ch1->cd(5);
		hocc_lay[15]->Draw();
		pad_hits_ch1->cd(6);
		hocc_lay[14]->Draw();
		pad_hits_ch1->cd(7);
		hocc_lay[13]->Draw();
		pad_hits_ch1->cd(8);
		hocc_lay[12]->Draw();
		pad_hits_ch1->cd(9);
		hocc_lay[11]->Draw();
		pad_hits_ch1->cd(10);
		hocc_lay[10]->Draw();
		pad_hits_ch1->cd(11);
		hocc_lay[9]->Draw();
		pad_hits_ch1->cd(12);
		hocc_lay[8]->Draw();
		pad_hits_ch1->Update();

		pad_time_ch1->cd(1);
		htime_sl[4]->Draw();
		pad_time_ch1->cd(2);
		htime_sl[3]->Draw();
		pad_time_ch1->cd(3);
		htime_sl[2]->Draw();


		//SLs plots
		sls->cd();
		pad_hits_sls->cd(1);
		hocc_lay[7]->Draw();
		pad_hits_sls->cd(2);
		hocc_lay[6]->Draw();
		pad_hits_sls->cd(3);
		hocc_lay[5]->Draw();
		pad_hits_sls->cd(4);
		hocc_lay[4]->Draw();
		pad_hits_sls->cd(9);
		hocc_lay[3]->Draw();
		pad_hits_sls->cd(10);
		hocc_lay[2]->Draw();
		pad_hits_sls->cd(11);
		hocc_lay[1]->Draw();
		pad_hits_sls->cd(12);
		hocc_lay[0]->Draw();
		pad_hits_sls->Update();

		pad_time_sls->cd(1);
		htime_sl[1]->Draw();
		pad_time_sls->cd(3);
		htime_sl[0]->Draw();

	}	 

    }//end loop on words
  return;
}

void fillMap()
{
  if(DEBUG_FLAG) 
	cout << "reading map" << endl;

  ifstream file("legnaro2ROS25.txt");

  int ros;
  int rob;
  int tdc;
  int cha;

  int se;
  int sl;
  int lay;
  int tube;

  int ddu;
  int wh;
  int st;

  if(DEBUG_FLAG) 
	cout << "map loop ... " << endl;
  while( file >> ddu >> ros >> rob >> tdc >> cha >> wh >> st >> se >> sl >> lay >> tube ) 
  {
  //	if(DEBUG_FLAG) cout << "fillMap "
  //               << ddu << " " << ros << " " << rob << " " << tdc << " " << cha  << " => "
  //               << wh << " " << st << " " << se << " " << " " <<  sl << " " << lay << " " << tube << endl;

	int idTDC = getTDCid(ros,rob,tdc,cha);
        int idTube = getTubeId(se,sl,lay,tube);

	chmap.insert( make_pair( idTDC, idTube ) ); 

  }
  if(DEBUG_FLAG) 
	cout << "TDC map read" << endl;
  return;
}


int getTDCid( int ros, int rob, int tdc, int cha ) {
  return ( ( ( ( ( ros * nrob ) + rob ) * ntdc ) + tdc ) * ncha ) + cha;
}

int getTubeId( int se, int sl, int lay, int tube ) {
  return (se * 10000 + sl * 1000 + lay * 100 + tube);
}

int getTube(int idTube) {
  return (idTube - int(idTube/100.)*100.);
}

int getLay(int idTube) {
  return (int(idTube/100.) - int(idTube/1000.)*10);
}

int getSL(int idTube) {
  return (int(idTube/1000.) - int(idTube/10000.)*10);
}

int getSe(int idTube) {
  return int(idTube/10000.);
}


