#include <algorithm>

#include "EventBuilder.h"
#include "TSystem.h"

EventBuilder::EventBuilder()
{
    m_debug = false;

    Int_t iLay=0;
    map_detID[10]=iLay++;
    map_detID[20]=iLay++;
    map_detID[30]=iLay++;
    map_detID[40]=iLay++;
    map_detID[50]=iLay++;
    map_detID[51]=iLay++;
    map_detID[55]=iLay++;
    map_detID[56]=iLay++;
    map_detID[70]=iLay++;

    return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::setDebug(bool debug){
    m_debug = debug;

    std::cout << "Debug flag  " << m_debug << std::endl;
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::setAlignment(bool doAli){
    doAlignment = doAli;

    std::cout << "Set Local->Global alignment flag  " << doAlignment << std::endl;
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::openDataFiles(std::string inputDTFile, std::string inputSiFile, std::string outputFile){
//     std::cout << "Debug flag  " << m_debug << std::endl;
    /// DT root file
  m_dtFile = new TFile((inputDTFile+".root").c_str());
  std::cout << "Opening mu data file " << (inputDTFile+".root").c_str() << std::endl;
  m_tree = (TTree*)m_dtFile->Get("RADMU");
  
//   m_tree->SetBranchAddress( "EVENT"  ,  &dtEvent);
//   m_tree->SetBranchAddress( "SEG_ns" ,  &nseg   );
//   m_tree->SetBranchAddress( "SEG_sn" ,  &segN   );
//   m_tree->SetBranchAddress( "SEG_xh1",  &s1p    );
//   m_tree->SetBranchAddress( "SEG_xh2",  &s2p    );
//   m_tree->SetBranchAddress( "SEG_xh3",  &s3p    );
//   m_tree->SetBranchAddress( "SEG_xh4",  &s4p    );
//   m_tree->SetBranchAddress( "SEG_xh5",  &s5p    );
//   m_tree->SetBranchAddress( "SEG_xh6",  &s6p    );
//   m_tree->SetBranchAddress( "SEG_xh7",  &s7p    );
//   m_tree->SetBranchAddress( "SEG_xh8",  &s8p    );
  m_tree->SetBranchAddress("EVENT"  , &dtEvent, &b_EVENT  );
  m_tree->SetBranchAddress("SEG_ns" , &nseg   , &b_ntes   );
  m_tree->SetBranchAddress("SEG_sn" ,  segN   , &b_SEG_sn );
  m_tree->SetBranchAddress("SEG_xh1",  s1p    , &b_SEG_xh1);
  m_tree->SetBranchAddress("SEG_xh2",  s2p    , &b_SEG_xh2);
  m_tree->SetBranchAddress("SEG_xh3",  s3p    , &b_SEG_xh3);
  m_tree->SetBranchAddress("SEG_xh4",  s4p    , &b_SEG_xh4);
  m_tree->SetBranchAddress("SEG_xh5",  s5p    , &b_SEG_xh5);
  m_tree->SetBranchAddress("SEG_xh6",  s6p    , &b_SEG_xh6);
  m_tree->SetBranchAddress("SEG_xh7",  s7p    , &b_SEG_xh7);
  m_tree->SetBranchAddress("SEG_xh8",  s8p    , &b_SEG_xh8);
  m_idt = 0;
  
  /// Si root file;
  //   TFile f2("treeSiliconAligned.root","recreate");
  m_siFile = new TFile((inputSiFile).c_str());
  std::cout << "Opening Si data file " << (inputSiFile).c_str() << std::endl;
  t2       = (TTree*)m_siFile->Get("LEMMA");

  t2->SetBranchAddress("iev"   , &si_iev,    &b_iev);
  t2->SetBranchAddress("iev_mu", &si_ievmu,  &b_ievmu);
  t2->SetBranchAddress("nhits" , &si_nhits,  &b_nhits);
  t2->SetBranchAddress("subdet",  si_subdet, &b_subdet);
  t2->SetBranchAddress("xh"    ,  si_xh,     &b_xh);
  t2->SetBranchAddress("yh"    ,  si_yh,     &b_yh);
  t2->SetBranchAddress("zh"    ,  si_zh,     &b_zh);
  t2->SetBranchAddress("itrack",  si_itrack, &b_itrack);
  m_isi = 0;

    /// output file
    m_outFile = new TFile(outputFile.c_str(),"recreate");
    m_outTree = new TTree("LEMMA","Event Builder LEMMA data output");
    m_outTree->Branch("siiev" , &siiev , "siiev/I");
    m_outTree->Branch("iev"   , &iev   , "iev/I");
    m_outTree->Branch("nhits" , &nhits , "nhits/I");
    m_outTree->Branch("subdet", &subdet, "subdet[nhits]/I");
    m_outTree->Branch("xh"    , &xh    , "xh[nhits]/F");
    m_outTree->Branch("yh"    , &yh    , "yh[nhits]/F");
    m_outTree->Branch("zh"    , &zh    , "zh[nhits]/F");
    m_outTree->Branch("itrack", &itrack, "itrack[nhits]/F");

    
    openSiTxtFiles();


    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::openSiTxtFiles(){
  std::string line;
  myfile.open("../../data/run200308_multi.dat");
  isSiFileOpen = false;
  if (myfile.is_open()) {
    isSiFileOpen = true; 
    std::cout << "Si data file loaded." << std::endl;
  }
  else std::cout << "\033[1;31mERROR: could not open Si data file.\033[0m" << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::matchEvents(int nEvents){
  int nsi=0;
  int ndt=0;
   getDtEvent();
   getDtEvent();
  while(m_idt<std::min(nEvents,int(m_tree->GetEntries())) && int(m_tree->GetEntries())&&
	m_isi<std::min(nEvents,int(t2    ->GetEntries())) && int(t2    ->GetEntries())){

        std::cout << "Debug flag  " << m_debug << std::endl;
        std::cout << "\n --- Matching event " << m_idt << std::endl;

        int dtN = getDtEvent(); if (dtN<0) { std::cout << "\033[1;35m CORBEZZOLI \033[0m" << std::endl; dtN = getDtEvent();}// PAPOCCHIO
	int siN = getSiEvent();
	//int siN = getSiTxtEvent();

	std::cout << "DtN " << dtN << "\t SiN " << siN << std::endl;

        while(siN < dtN) {
	  siN=getSiEvent();
	  nsi++;
	  std::cout << "siN++" << std::endl;
	if (nsi> int(t2    ->GetEntries()) ||
	    ndt> int(m_tree->GetEntries())) break;
	}

        while(siN > dtN) {
	  ndt++;
	  dtN=getDtEvent(); if (dtN<0) { std::cout << "\033[1;35m CORBEZZOLI \033[0m" << std::endl; dtN = getDtEvent();}// PAPOCCHIO
	  std::cout << "dtN++" << std::endl;
	if (nsi> int(t2    ->GetEntries()) ||
	    ndt> int(m_tree->GetEntries())) break;
	}

        if   (siN == dtN){
            std::cout << "\033[1;31mEvent " << dtN << " matched! \033[0m" << std::endl;
	    dumpGlobalEvent();
	}
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int EventBuilder::getDtEvent(){

    m_tree->GetEntry(m_idt);
    std::cout << "DT event " << dtEvent << std::endl;
    m_idt++;
    return dtEvent;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int EventBuilder::getSiEvent(){

  t2->GetEntry(m_isi);
  std::cout << "Si event " << si_ievmu << std::endl;
  ++m_isi;
  return si_ievmu;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
int EventBuilder::getSiTxtEvent(){

  if (isSiFileOpen) {
    std::string line;
    Int_t nLines=0;
    std::vector<float> vRecord;
    while ( getline (myfile,line) ) {
      TString sline(line);
      TString token;
      
      TObjArray* slineArray = sline.Tokenize(" ");
      Int_t ntokens = slineArray->GetEntries();
      for (Int_t it=0; it<ntokens; ++it) { // Store all tokes of each line
	TObjString* objString = (TObjString*) (slineArray->At(it));
	TString stoken = objString->String();
	// 	cout << stoken.Data() << endl;
	if (!stoken.IsFloat()) {
	  std::cout << stoken.Data() << std::endl;
	  std::cout << "\033[1;31mERROR in reading. Non float entry in record\033[0m" << std::endl;
	  return 2;
	}
	vRecord.push_back(stoken.Atof()+1E-12); 
      }
      std::cout << line << std::endl;
      if (++nLines==13) break;
    }
    for (Int_t ihit=0; ihit<nMaxHits; ++ihit) {zh[ihit] = 0; zh[ihit] = 0; zh[ihit] = 0;}
    si_nhits=0;
    Ssiz_t it=0;
    si_iev = vRecord.at(it++);
    si_xh[0]     = vRecord.at(it++);
    si_yh[1]     = vRecord.at(it++);
    si_xh[2]     = vRecord.at(it++);
    si_yh[3]     = vRecord.at(it++);
    si_subdet[0] = detIDs[0];
    si_subdet[1] = detIDs[1];
    si_subdet[2] = detIDs[2];
    si_subdet[3] = detIDs[3];
    si_nhits+=4;
    it++; //triggerInfo
    si_ievmu = vRecord.at(it++);
    it+=48; //digis
    for (Int_t iLayer=2; iLayer<nSiLayers; iLayer++) {
      si_nhits += vRecord.at(it++);
      for (Int_t iCluster=0; iCluster<nMaxClusters; iCluster++) {
	si_subdet[it] = detIDs[iLayer];
	si_xh    [it] = vRecord.at(it++);
      }
      si_nhits += vRecord.at(it++);
      for (Int_t iCluster=0; iCluster<nMaxClusters; iCluster++) {
	si_subdet[it] = detIDs[iLayer];
	si_xh    [it] = vRecord.at(it++);
      }
    }
    std::cout << it << " " << vRecord.size() << std::endl;
    vRecord.clear();
  }
  return si_iev;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::calibrate(){

    /// TO BE COMPLETED!

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::openAlignments(std::string alignment_file){
  std::cout << "Opening alignment file: " << alignment_file << std::endl;

//   float x0SiAndMu[nSiLayers+1][3];
//   float phi0Theta0Mu        [2];
  
  std::string line;
  ifstream myfile (alignment_file);
  std::vector<float> vAli;
  if (myfile.is_open())
    {
      while( getline(myfile, line) )
	{
	  std::stringstream ss(line); 
	  std::string token;
	  while (std::getline(ss, token, ' ')){
	    if (!TString(token).IsFloat()) {
	      std::cout << "\033[1;31mERROR. Non float entry in alignment file\033[0m" << std::endl;
	      return;
	    }
	    vAli.push_back(TString(token).Atof());
	  }   
	}
      if ((int)vAli.size() != 3*(nSiLayers+1)+2) {
	std::cout << "\033[1;31mERROR in parsing alignment file\033[0m" << std::endl;
	return;
      }
      for (int iAli=0; iAli<3*(nSiLayers+1); ++iAli) {
	int iLay=iAli/3;
	int xORy=iAli%3;
	x0SiAndMu[iLay][xORy] = vAli.at(iAli);
      }
      phi0Theta0Mu[0] = vAli.at(3*(nSiLayers+1)  );
      phi0Theta0Mu[1] = vAli.at(3*(nSiLayers+1)+1);

    }

  std::cout << "xh \t yh \t zh" << std::endl;
  for (Int_t iLayer=0; iLayer<nSiLayers+1; iLayer++) {
    std::cout << x0SiAndMu[iLayer][0] << "\t"<< x0SiAndMu[iLayer][1] << "\t"<< x0SiAndMu[iLayer][2] << std::endl;
  }
  std::cout << "phi0mu \t theta0mu" << std::endl;
  std::cout <<  phi0Theta0Mu[0] << "\t" << phi0Theta0Mu[1] << std::endl;

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::dumpGlobalEvent(){

  /// Si data
  // n hits from Si
  int nSiHits = si_nhits;

  if (nSiHits+12>nMaxHits) {
    std::cout << "\033[1;31m Warning: too many hits. Skipping event. \033[0m" << std::endl;
    return;
  }

  siiev = si_iev;
  for (Int_t iSiHit=0; iSiHit<nSiHits; ++iSiHit) {
   subdet [iSiHit] = si_subdet[iSiHit];
   xh     [iSiHit] = si_xh    [iSiHit];
   yh     [iSiHit] = si_yh    [iSiHit];
   zh     [iSiHit] = si_zh    [iSiHit];
   itrack [iSiHit] = si_itrack[iSiHit];
  }
  /// DT data
    iev    = dtEvent;
    nhits = nSiHits + 12;
    if(nseg!=2)
        std::cout << "\033[1;34m SOMETHING NASTY..... more than 2 segments fitted? \033[0m" << std::endl;

    // set default for variables x,y, subdet
    for(int ih=0 + nSiHits; ih<12 + nSiHits; ih++){
      xh[ih] = -999.;
      yh[ih] = -999.;
      zh[ih] = -999.;
      subdet[ih] = 70;
    }

    // set z of the layers
    float z_layers[12] = {-10.75, -9.45, -8.15, -6.85,  7.45, 8.75, 10.05, 11.35,12.85, 14.15, 15.45, 16.75};
    
    // fill HITS - Z
    for(int ih=0 ; ih<12 ; ih++) {
//       xh[ih+nSiHits] = dt_xh[iSiHit];
//       yh[ih+nSiHits] = dt_yh[iSiHit];
      zh[ih+nSiHits] = z_layers[ih];
    }

    for ( int is = 0; is < nseg; is++ ){
        if(segN[is]<0){
            // SL THETA hits
            yh[4 + nSiHits] = s1p[is];
            yh[5 + nSiHits] = s2p[is];
            yh[6 + nSiHits] = s3p[is];
            yh[7 + nSiHits] = s4p[is];
        } else {
            // SL PHI hits
            xh[0 + nSiHits] = s1p[is];
            xh[1 + nSiHits] = s2p[is];
            xh[2 + nSiHits] = s3p[is];
            xh[3 + nSiHits] = s4p[is];
            xh[8 + nSiHits] = s5p[is];
            xh[9 + nSiHits] = s6p[is];
            xh[10 + nSiHits] = s7p[is];
            xh[11 + nSiHits] = s8p[is];
        }
    }

    if (doAlignment) // Apply alignment to get global coordinates
    for (Int_t ih = 0; ih < nhits; ++ih) {
      std::map<Int_t,Int_t>::iterator it_map_detID;
      it_map_detID = map_detID.find(subdet[ih]);
      Int_t iLayer=-1;
      if (it_map_detID != map_detID.end()) iLayer = it_map_detID->second;
//       std::cout << iLayer << " " << x0SiAndMu[iLayer][0] << " " << x0SiAndMu[iLayer][1] << " " << x0SiAndMu[iLayer][2] << " " << zh[ih] << std::endl;
      xh[ih] += x0SiAndMu[iLayer][0];
      yh[ih] += x0SiAndMu[iLayer][1];
      zh[ih] += x0SiAndMu[iLayer][2];
    }


    if(m_debug==true){
        std::cout << "---- EventBuilder::dumpGlobalEvent " << std::endl;
        std::cout << "Filling DT data..." << std::endl;
        std::cout << " ---> event ID " << dtEvent << std::endl;
        std::cout << " ---> muon track nhits " << nhits << std::endl;
    }

    /// fill tree
    m_outTree->Fill();

    return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::dumpOutput(){

    m_outTree->Write();
    m_outFile->Close();

    return;
}
