#ifndef EVENTBUILDER_H
#define EVENTBUILDER_H

//---- STL
#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TTree.h"


class EventBuilder
{
public:
    EventBuilder();
    ~EventBuilder(){}

    void setDebug(bool debug);
    void openDataFiles(std::string inputDTFile, std::string inputSiFile, std::string outputFile);
    void matchEvents(int nEvents);
    int getDtEvent();
    int getSiEvent();
    void dumpGlobalEvent();
    void dumpOutput();
    void calibrate();

private:
    /// si data file
    std::ifstream m_siFile;

    /// output file
    TFile *       m_outFile;
    TTree *     m_outTree;
    Int_t           iev;
    Int_t           nhits;
    Int_t           subdet[100];
    Double_t        xh[100];
    Double_t        yh[100];
    Double_t        zh[100];
    Int_t           itrack[100];
    
    /// dt data file
    TFile * m_dtFile;
    TTree * m_tree;
    int m_idt;
    int dtEvent, nseg;
    Int_t segN[2];
    Float_t s1p[50],s2p[50], s3p[50], s4p[50], s5p[50], s6p[50], s7p[50], s8p[50];

    bool m_debug;
};

#endif // EVENTBUILDER_H
