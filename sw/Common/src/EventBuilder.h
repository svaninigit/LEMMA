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
    TFile * m_outFile;
    TTree * m_outTree;
    Int_t     iev;
    Int_t     nhits;
    Int_t           subdet[100];
    Double_t        xh[100];
    Double_t        yh[100];
    Double_t        zh[100];
    Int_t           itrack[100];
    
    /// dt data file
    TFile * m_dtFile;
    TTree * m_tree;
    int m_idt;
    int dtEvent, nseg, segNpoints;
    float segX, segS, segK, segT0, s1r, s2r, s3r, s4r, s5r, s6r, s7r, s8r;
};

#endif // EVENTBUILDER_H
