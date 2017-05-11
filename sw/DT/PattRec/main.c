// 16th January 2010
// Silvia Pesente : simple program to read raw data and do Patter Recognition
// SP usage example: ./runPR file_name_without_extension runN RunTrig MAXEvent IF_compute_TTrig
//                   ./runPR trg_camera_sup_bti_opened_thr25mv_r966 966 966 500000 0
//

#include "RawAnalyzer.h"

RawAnalyzer *analyze=NULL;

int main(int argc, char *argv[])
{
  
  char *fin=NULL;
  int runN=0;      // Run number to be read
  int runTrig=0;   // Number of run to be used as TRIG
  int maxEvent=0;  // Max. number of events to read
  int runID_min=0; // First run to read
  bool ttrig=0;    // if 0, go with Pattern Recognition, if 1 compute ttrig
  if (argc <7) {
    printf("Too few parameters!\nYou must write:\n   ./runPR nome_file (without extension) runN runTrig MAXEvent runID_min IF_compute_ttrig\n");
    return(1);
  } else if (argc == 7) {
    fin = argv[1];
    runN = atoi(argv[2]);
    runTrig = atoi(argv[3]);
    maxEvent = atoi(argv[4]);
    runID_min = atoi(argv[5]);
    ttrig = atoi(argv[6]);
  } else {
    printf("Too many parameters!\n");
    return(1);
  }  
  
  analyze=new RawAnalyzer();
  analyze->goAnalysis(fin, maxEvent, runN, runID_min, runTrig, ttrig);  
  delete analyze;
  
  return 0;
}

