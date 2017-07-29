#ifndef EVENTBUILDER_H
#define EVENTBUILDER_H

//---- STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "TFile.h"
#include "TTree.h"

using namespace std;

class EventBuilder
{
public:
    EventBuilder();
    ~EventBuilder(){}

    void setDebug(bool debug);
    void openAlignments(std::string alignment_file);
    void openDataFiles(std::string inputDTFile, std::string inputSiFile, std::string outputFile);
    void matchEvents(int nEvents);
    int getDtEvent();
    int getSiEvent();
    void dumpGlobalEvent();
    void dumpOutput();
    void calibrate();
    void setAlignment(bool doAli);
    void openSiTxtFiles();
    int getSiTxtEvent();
private:
    static const Int_t nMaxHits=100;
    static const Int_t nSiLayers = 8; // Si layers
    static const Int_t nMaxClusters = 5;

    
  //// si data file
  TFile * m_siFile;
  TTree * t2;
  // Declaration of leaf types
  Int_t           si_iev;
  Int_t           si_nhits;
  Int_t           si_ievmu;
  Int_t           si_subdet[nMaxHits];
  Float_t         si_xh[nMaxHits];
  Float_t         si_yh[nMaxHits];
  Float_t         si_zh[nMaxHits];
  Int_t           si_itrack[nMaxHits];

  Int_t           m_isi;

  // List of branches
  TBranch        *b_iev;   //!
  TBranch        *b_ievmu;   //!
  TBranch        *b_nhits;   //!
  TBranch        *b_subdet;   //!
  TBranch        *b_xh;   //!
  TBranch        *b_yh;   //!
  TBranch        *b_zh;   //!
  TBranch        *b_itrack;   //!
  

  /// output file
    TFile *       m_outFile;
    TTree *     m_outTree;
    Int_t           iev;
    Int_t           siiev;
    Int_t           nhits;
    Int_t           subdet[nMaxHits];
    Float_t         xh[nMaxHits];
    Float_t         yh[nMaxHits];
    Float_t         zh[nMaxHits];
    Int_t           itrack[nMaxHits];
    
    /// dt data file
    TFile * m_dtFile;
    TTree * m_tree;
    int m_idt;

   Int_t           dtEvent;
   Int_t           nseg;
   Int_t           segN[2];
   Float_t         s1p[2];   //[ntes]
   Float_t         s2p[2];   //[ntes]
   Float_t         s3p[2];   //[ntes]
   Float_t         s4p[2];   //[ntes]
   Float_t         s5p[2];   //[ntes]
   Float_t         s6p[2];   //[ntes]
   Float_t         s7p[2];   //[ntes]
   Float_t         s8p[2];   //[ntes]
   // List of branches
   TBranch        *b_EVENT;   //!
   TBranch        *b_ntes;   //!
   TBranch        *b_SEG_sn;   //!
   TBranch        *b_SEG_xh1;   //!
   TBranch        *b_SEG_xh2;   //!
   TBranch        *b_SEG_xh3;   //!
   TBranch        *b_SEG_xh4;   //!
   TBranch        *b_SEG_xh5;   //!
   TBranch        *b_SEG_xh6;   //!
   TBranch        *b_SEG_xh7;   //!
   TBranch        *b_SEG_xh8;   //!

    bool m_debug;
    bool doAlignment;

    bool isSiFileOpen;

    // Alignment variables
    float x0SiAndMu[nSiLayers+1][3];
    float phi0Theta0Mu          [2];

    std::map<Int_t,Int_t> map_detID;
    ifstream myfile;
    const Int_t detIDs[nSiLayers+1]={10,20,30,40,50,51,55,56,70};
};

#endif // EVENTBUILDER_H
