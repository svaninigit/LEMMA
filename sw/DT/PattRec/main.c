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
#include "src/ReaderROS8.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RECIPE PARAMETERS //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
    static struct Parameters : Options {
        std::string inFileName;
        std::string outFileName;
        int runNum;
        int nEvents;
        int ttrigRunNum;
        bool ttrigFlag;
        std::string execute;

        struct Unpack {
            bool debug;
        } unpack;

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
        ("nEvents",                 &nEvents,                (int)10.,                                "number of events to read")
        ("ttrigRunNum",         &ttrigRunNum,        (int)0.,                                 "Run number ID to compute time-trig calibration")
        ("ttrigFlag",                &ttrigFlag,                 (bool)0.,                              "Flag for activate the ttrig computation")
        ("execute",                 &execute,                 std::string("pattrec"),         "Execution options: ttrig, pattrec, analysis, dump - TO BE IMPLEMENTED")

        // UNPACK //
        ("unpack.debug",       &unpack.debug,     false,                                     "enable debugging dumps in unpacking code")

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
    std::string config_file("utils/parameters.cfg");
    p.add_options() ("config",&config_file,"set config file");
    p.parse_command_line(argc,argv);
    p.parse_config_file(config_file);
    p.parse_command_line(argc,argv);

//    //--- Recipe for ROS25 data
//   RawAnalyzer * analyze=new RawAnalyzer();
//   analyze->goAnalysis(p.inFileName, p.nEvents, p.runNum, p.ttrigRunNum, p.ttrigFlag);

    //--- Recipe for ROS8 data
    ReaderROS8 *analyze=new ReaderROS8();
    analyze->setDebug(p.unpack.debug);
    analyze->goAnalysis(p.inFileName, p.nEvents, p.runNum, p.ttrigRunNum, p.ttrigFlag);

   delete analyze;

    return 0;
}





