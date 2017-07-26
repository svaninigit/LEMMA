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
    Float_t        xh[100];
    Float_t        yh[100];
    Float_t        zh[100];
    Int_t           itrack[100];
    
    /// dt data file
    TFile * m_dtFile;
    TTree * m_tree;
    int m_idt;
    int dtEvent, nseg;
    Int_t segN[2];
    Float_t s1p[2],s2p[2], s3p[2], s4p[2], s5p[2], s6p[2], s7p[2], s8p[2];

    bool m_debug;
};

#endif // EVENTBUILDER_H
