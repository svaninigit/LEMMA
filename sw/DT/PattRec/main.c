////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 20170515 Sara Vanini DT software for LEMMA testbeam
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//---- STL
#include <iostream>
#include <fstream>
#include <map>

//---- Core
#include "../../Common/src/Options.h"

//---- DT classes
#include "src/RawAnalyzer.h"
#include "src/ReaderROS8.h"
#include "../../Common/src/EventBuilder.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RECIPE PARAMETERS //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
    static struct Parameters : Options {
        std::string inputDTFile;
        std::string rootDTFile;
        std::string inputSiFile;
        std::string outputFile;
        int runNum;
        int nEvents;
        int ttrigRunNum;
        bool ttrigFlag;
        bool n2chambers;
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
        ("inputDTFile",           &inputDTFile,          std::string("test"),              "DT  chamber raw data file")
        ("rootDTFile",           &rootDTFile,            std::string("test"),              "DT  chamber reconstructed data file")
        ("inputSiFile",            &inputSiFile,             std::string("test"),              "Si detectors raw data file")
        ("outputFile",             &outputFile,              std::string("test"),              "Output file name")
        ("runNum",                 &runNum,                (int)0.,                                  "Run number ID")
        ("nEvents",                 &nEvents,                (int)10.,                                "number of events to read")
        ("ttrigRunNum",         &ttrigRunNum,        (int)0.,                                 "Run number ID to compute time-trig calibration")
        ("n2chambers",         &n2chambers,         (bool)0,                               "Flag for activate 2 chambers analysis")
        ("execute",                 &execute,                 std::string("pattrec"),         "Execution options: ttrig, pattrec, eventbuild")

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
    ReaderROS8 *reader=new ReaderROS8();
    reader->setDebug(p.unpack.debug);

    /// ttrig computation
    if(p.execute=="ttrig"){
        std::cout << "\n\n *** RUNNING time-trigger computation *** \n\n";
        reader->goAnalysis(p.inputDTFile, p.nEvents, p.runNum, p.ttrigRunNum, 1, p.n2chambers);
        delete reader;
    }

    /// patter recognition
    else if (p.execute=="pattrec"){
        std::cout << "\n\n *** RUNNING pattern recognition *** \n\n";
        reader->goAnalysis(p.inputDTFile, p.nEvents, p.runNum, p.ttrigRunNum, 0, p.n2chambers);
        delete reader;
    }

    /// event builder
    else if (p.execute == "eventbuilder"){
        std::cout << "\n\n *** RUNNING event builder *** \n\n";

        EventBuilder *builder = new EventBuilder();
        builder->openDataFiles(p.rootDTFile,p.inputSiFile,p.outputFile);
        builder->matchEvents(p.nEvents);
        builder->dumpOutput();

        delete builder;
    }
    else if (p.execute == "geotest"){
        Geom * geo = new Geom();

        std::cout << "X WIRE POSITIONS " << std::endl;
        for(int isl=1; isl <= 3; isl++)
            for(int iw=1; iw < 60; iw++)
                for(int il=1; il <= 1; il++)
                        std::cout << "SL " << isl << ", LAY " << il << ", WIRE " << iw << " ---> X wire " << geo->get_x_wire(11,isl,il,iw) << std::endl;
        
        std::cout << "Y WIRE POSITIONS " << std::endl;
         for(int isl=1; isl <= 3; isl++)
            for(int il=1; il <= 4; il++)
                 std::cout << "SL " << isl << ", Lay " << il << " ---> " << geo->get_y_wire(11,isl,il,1) << std::endl;
        }

    return 0;
}





