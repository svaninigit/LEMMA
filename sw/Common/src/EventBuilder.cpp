#include <algorithm>

#include "EventBuilder.h"
#include "TSystem.h"

EventBuilder::EventBuilder()
{
    //gDebug=2;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void EventBuilder::openDataFiles(std::string inputDTFile, std::string inputSiFile, std::string outputFile){

    /// DT root file
    m_dtFile = new TFile(inputDTFile.c_str());
    m_tree = (TTree*)m_dtFile->Get("RADMU");

    m_tree->SetBranchAddress( "EVENT",  &dtEvent);
    m_tree->SetBranchAddress( "SEG_ns", &nseg);
    m_tree->SetBranchAddress( "SEG_sx",  &segX);
    m_tree->SetBranchAddress( "SEG_ss",  &segS);
    m_tree->SetBranchAddress( "SEG_sk",  &segK);
    m_tree->SetBranchAddress( "SEG_sn",  &segNpoints);
    m_tree->SetBranchAddress( "SEG_t0",  &segT0);
    m_tree->SetBranchAddress( "SEG_1r",  &s1r);
    m_tree->SetBranchAddress( "SEG_2r",  &s2r);
    m_tree->SetBranchAddress( "SEG_3r",  &s3r);
    m_tree->SetBranchAddress( "SEG_4r",  &s4r);
    m_tree->SetBranchAddress( "SEG_5r",  &s5r);
    m_tree->SetBranchAddress( "SEG_6r",  &s6r);
    m_tree->SetBranchAddress( "SEG_7r",  &s7r);
    m_tree->SetBranchAddress( "SEG_8r",  &s8r);
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
        std::cout << "\n --- Matching event " << m_idt << std::endl;

        int dtN = getDtEvent();
        int siN = getSiEvent();

        while(siN < dtN)
            siN=getSiEvent();

        if(dtN==siN){
            std::cout << "Events matched! " << std::endl;
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

    /// fill DT data
    iev    = dtEvent;
    nhits = abs(segNpoints) - 1300;

    for(int i=0; i<nhits; i++)
        subdet[i] = 70;

    // to be fixed with correct coordinates instead of residuals
    if(segNpoints > 0){
        xh[0]     = s1r;
        xh[1]     = s2r;
        xh[2]     = s3r;
        xh[3]     = s4r;
        xh[4]     = -999;
        xh[5]     = -999;
        xh[6]     = -999;
        xh[7]     = -999;
        xh[8]     = s5r;
        xh[9]     = s6r;
        xh[10]     = s7r;
        xh[11]     = s8r;
    }
    else{
        xh[0]     = -999;
        xh[1]     = -999;
        xh[2]     = -999;
        xh[3]     = -999;
        xh[4]     = s1r;
        xh[5]     = s2r;
        xh[6]     = s3r;
        xh[7]     = s4r;
        xh[8]     = -999;
        xh[9]     = -999;
        xh[10]     = -999;
        xh[11]     = -999;
    }

    // to be fixed with correct coordinates
    float yh_init = 20.;
    for(int i=0; i<12;i++){
        yh[i] = yh_init + 2.;
        zh[i] = 25.60;
    }

    /// SV FIX fill itrack, track ID 

    /// fill Si data -----> TO BE COMPLETED !!!

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
