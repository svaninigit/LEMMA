////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 20170515 Sara Vanini DT software for LEMMA testbeam
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//---- STL
#include <iostream>
#include <fstream>
#include <map>

//---- Core
#include "src/Options.h"

//---- DT classes
#include "src/RawAnalyzer.h"
RawAnalyzer *analyze=NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RECIPE PARAMETERS //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
    static struct Parameters : Options {
        std::string inFileName;
        std::string outFileName;
        int runNum;
        int ttrigRunNum;
        int nEvents;
        std::string execute;

        struct PattRec {
            bool debug;
        } pattrec;

        Parameters(const char *hello = "Program options") : Options(hello) {
        add_options()
        ("help", "printout help")

        // GENERAL //
        ("inFileName",           &inFileName,           std::string("test"),              "Input file name")
        ("outFileName",         &outFileName,        std::string("test"),              "Output root file name")
        ("runNum",                 &runNum,                (int)0.,                                 "Run number ID")
        ("ttrigRunNum",         &ttrigRunNum,        (int)0.,                                 "Run number ID to compute time-trig calibration")
        ("nEvents",                 &nEvents,                (int)10.,                                "number of events to read")
        ("execute",                 &execute,                 std::string("pattrec"),         "Execution options: ttrig, pattrec, analysis, dump")

        // PATTREC //
        ("pattrec.debug",       &pattrec.debug,     false,                                      "enable debugging dumps in pattrec code")
        ;
        }
    } p;   // <-- INSTANCE //
} // namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

    //--- Avoid annoying ROOT reader warnings
    gErrorIgnoreLevel=kError;

    //--- Parameters
    std::string config_file("parameters.cfg");
    p.add_options() ("config",&config_file,"set config file");
    p.parse_command_line(argc,argv);
    p.parse_config_file(config_file);
    p.parse_command_line(argc,argv);

    //--- Recipe
    analyze=new RawAnalyzer();
    analyze->goAnalysis(fin, maxEvent, runN, runID_min, runTrig, ttrig);
    delete analyze;

    return 0;
}





