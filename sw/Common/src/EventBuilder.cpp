#include <algorithm>

#include "EventBuilder.h"
#include "TSystem.h"

EventBuilder::EventBuilder()
{
    m_debug = false;
    return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::setDebug(bool debug){
    m_debug = debug;

    std::cout << "Debug flag  " << m_debug << std::endl;
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::openDataFiles(std::string inputDTFile, std::string inputSiFile, std::string outputFile){
    std::cout << "Debug flag  " << m_debug << std::endl;
    /// DT root file
    m_dtFile = new TFile(inputDTFile.c_str());
    m_tree = (TTree*)m_dtFile->Get("RADMU");

    m_tree->SetBranchAddress( "EVENT",  &dtEvent);
    m_tree->SetBranchAddress( "SEG_ns", &nseg);
    m_tree->SetBranchAddress( "SEG_sn",  &segN);
    m_tree->SetBranchAddress( "SEG_1r",  &s1p);
    m_tree->SetBranchAddress( "SEG_2r",  &s2p);
    m_tree->SetBranchAddress( "SEG_3r",  &s3p);
    m_tree->SetBranchAddress( "SEG_4r",  &s4p);
    m_tree->SetBranchAddress( "SEG_5r",  &s5p);
    m_tree->SetBranchAddress( "SEG_6r",  &s6p);
    m_tree->SetBranchAddress( "SEG_7r",  &s7p);
    m_tree->SetBranchAddress( "SEG_8r",  &s8p);
    m_idt = 0;

    /// Si txt file
    m_siFile.open(inputSiFile.c_str());

    /// output file
    m_outFile = new TFile(outputFile.c_str(),"recreate");
    m_outTree = new TTree("LEMMA","Event Builder LEMMA data output");
    m_outTree->Branch("iev", &iev, "iev/I");
    m_outTree->Branch("nhits", &nhits, "nhits/I");
    m_outTree->Branch("subdet", &subdet, "subdet[100]/I");
    m_outTree->Branch("xh", &xh, "xh[100]/F");
    m_outTree->Branch("yh", &yh, "yh[100]/F");
    m_outTree->Branch("zh", &zh, "zh[100]/F");
    m_outTree->Branch("itrack", &itrack, "itrack[100]/F");

    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::matchEvents(int nEvents){
    while(m_idt<std::min(nEvents,int(m_tree->GetEntries())) && !m_siFile.eof()){

        std::cout << "Debug flag  " << m_debug << std::endl;
        std::cout << "\n --- Matching event " << m_idt << std::endl;

        int dtN = getDtEvent();
        int siN = getSiEvent();

        while(siN < dtN)
            siN=getSiEvent();

//        if(dtN==siN){
            std::cout << "Event " << dtN << " matched! " << std::endl;
            dumpGlobalEvent();
//        }
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

    float Sievent;
    float Sidata;

    for(int i=0; i<93; i++)
        m_siFile >> Sidata;

    m_siFile >> Sievent;

    // just for testing....
    Sievent -= 20670;

    std::cout << "Si event " << Sievent << std::endl;

    return Sievent;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::calibrate(){

    /// TO BE COMPLETED!

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::dumpGlobalEvent(){

    /// Si data -----> TO BE COMPLETED !!!
    // n hits from Si
    int nSiHits = 0;

    /// DT data
    iev    = dtEvent;
    nhits = nSiHits + 12;
    if(nseg!=2)
        std::cout << "SOMETHING NASTY..... more than 2 segments fitted?" << std::endl;

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
    for(int ih=0 + nSiHits; ih<12 + nSiHits; ih++)
        zh[ih] = z_layers[ih];

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

    if(m_debug==true){
        std::cout << "---- EventBuilder::dumpGlobalEvent " << std::endl;
        std::cout << "Flling DT data..." << std::endl;
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
