#include "Track_IO.h"

using namespace std;


Track_IO::Track_IO(){
  
  infile=NULL;
  txtfile=NULL;
  fo_tree=NULL;
  fo_histo==NULL;
  HBfile==NULL;

  return;
}

Track_IO::~Track_IO(){
  delete fo_tree;
  delete fo_histo;
  delete HBfile;
  fo_tree=NULL;
  fo_histo==NULL;
  HBfile==NULL;
  return;
}


FILE * Track_IO::openINFile(char *fin, int runID, int ID){
  
  const char *fileName;
  TString fname(fin);
  TString fname_fin="";
  if(runID==-1){
    char * pch;
    pch=strrchr(fin,'.');
    //printf ("Last occurence of '.' found at %d \n",pch-fin+1);
    if(pch==NULL){
      printf("You have to pass file with the right extension (***.i*)!!!\n");
      return NULL;
    }
    int num=pch-fin+2;
    fname.Remove(num,3);
    fname += ID;
    fileName=fname.Data();

    cout<<fileName<<endl;
  }
  else {
    TString path;
    
    path = "/"; ///mudata/LNL/muonegrafia/"; //path over network
    fname_fin += path;
    fname_fin += fname;
    fname_fin += ".i";
    fname_fin += ID;
    
    fileName=fname_fin.Data();
    
    cout<<fileName<<endl;
  }
  
  infile = fopen(fileName,"rb");
  printf("Opening file: %s\n",fileName);
  if(infile==NULL){
    printf("file %s doesn't exist\n",fileName);
    return infile;
  }
  cout << "File open: " << fileName << endl;
  
  return infile;
}



FILE * Track_IO::openTXTFile(int runN, int maxEvent, int runID){

  char fileName[300];

//  if(lxradiomu)
//    sprintf(fileName,"/data/radmu/Patt_Rec/INFO/Statistics_r%d_%d_%dev.txt",runN,runID,maxEvent);
//  else
//  sprintf(fileName,"/data/tom_data/PattRec/INFO/Statistics_r%d_%d_%dev.txt",runN,runID,maxEvent);
  sprintf(fileName,"./OUTPUT/Statistics_r%d_%d_%dev.txt",runN,runID,maxEvent);

  txtfile = fopen(fileName,"w");
  printf("Opening file: %s\n",fileName);

  cout << "File open: " << fileName << endl;

  return txtfile;
}

TFile * Track_IO::openOUTRootFile(int runN, int maxEvent, int runID){
  
  char fileNameT[100];

//  if(lxradiomu)
//    sprintf(fileNameT,"/data/radmu/Patt_Rec/Radmufit_r%d_%d_%dev_PR.root",runN,runID,maxEvent);
//  else
//      sprintf(fileNameT,"/data/tom_data/PattRec/RootFile/Radmufit_r%d_%d_%dev_PR.root",runN,runID,maxEvent);
  sprintf(fileNameT,"./OUTPUT/Radmufit_r%d_%d_%dev_PR.root",runN,runID,maxEvent);

  
  cout<<"Opening file to store Tree: "<<fileNameT<<endl;
  // create a new ROOT file with a tree to store variables for analysis
  fo_tree = new TFile(fileNameT,"RECREATE");
  
  return fo_tree; 
}

TFile * Track_IO::openOUTHistoFile(int runN, int maxEvent, int runID){
  
  char nome[100];

//  if(lxradiomu)
//    sprintf(nome,"/data/radmu/Patt_Rec/INFO/prova_HIT_r%d_%d_%dev.root",runN,runID,maxEvent);
//  else
//      sprintf(nome,"/data/tom_data/PattRec/INFO/prova_HIT_r%d_%d_%dev.root",runN,runID,maxEvent);
  sprintf(nome,"./OUTPUT/prova_HIT_r%d_%d_%dev.root",runN,runID,maxEvent);

  cout<<"Opening file to store Histos: "<<nome<<endl;
  fo_histo = new TFile(nome,"RECREATE");
  
  return fo_histo; 
}

ofstream * Track_IO::openOUTHBFile(int runN, int maxEvent, int runID){
  
  char nome[100];

//  if(lxradiomu)
//    sprintf(nome,"/data_03/HitBank_files/HitBank_r%d_%d_%dev.dat",runN,runID,maxEvent);
//  else
//      sprintf(nome,"/data/tom_data/PattRec/HitBank/HitBank_r%d_%d_%dev.dat",runN,runID,maxEvent);
  sprintf(nome,"./HitBank/HitBank_r%d_%d_%dev.dat",runN,runID,maxEvent);

  cout<<"Opening file to store HitBank: "<<nome<<endl;
  HBfile = new ofstream();
  HBfile->open(nome, ios::out | ios::binary);
  
  
  return HBfile; 
}
