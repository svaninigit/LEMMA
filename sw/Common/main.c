////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 20170720 Sara Vanini Common software for LEMMA testbeam
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//---- STL
#include <iostream>
#include <fstream>
#include <map>

//---- Core
#include "src/Options.h"

//---- ROOT
#include <TError.h>

//---- LEMMA classes
#include "src/EventBuilder.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RECIPE PARAMETERS //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
    static struct Parameters : Options {
        std::string dirData;
        std::string rootDTFile;
        std::string inputSiFile;
        std::string outputFile;
        std::string runNum;
        int nEvents;
        std::string execute;
        bool doAlign;
        bool debug;

        Parameters(const char *hello = "Program options") : Options(hello) {
        add_options()
        ("help", "printout help")

        // GENERAL //
        ("dirData",                &dirData,               std::string("test"),              "dir data files")
/*         ("rootDTFile",             &rootDTFile,            std::string("test"),              "DT  chamber reconstructed data file") */
/*         ("inputSiFile",            &inputSiFile,           std::string("test"),              "Si detectors raw data file") */
        ("outputFile",             &outputFile,            std::string("test"),              "Output file name")
        ("runNum",                 &runNum,                std::string("4390"),              "Run number ID")
        ("nEvents",                &nEvents,               (int)10.,                                "number of events to read")
        ("execute",                &execute,                std::string("eventbuilder"), "Execution options: eventbuilder, PUT MORE IF NEEDED")
        ("doAlign",                &doAlign,               (bool)1,                                "Convert Local->Global")
        ("debug",                  &debug,                 (bool)0,                                 "Debugging flag")
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
    std::string alignment_file("utils/alignments.dat");
    p.add_options() ("config",&config_file,"set config file");
    p.parse_command_line(argc,argv);
    p.parse_config_file(config_file);
    p.parse_command_line(argc,argv);

    /// event builder
    if (p.execute == "eventbuilder"){
        std::cout << "\n\n *** RUNNING event builder *** \n\n";

        EventBuilder *builder = new EventBuilder();
        std::cout << "debug " << p.debug << std::endl;

        builder->setDebug(p.debug);
	builder->openAlignments(alignment_file);
	builder->setAlignment(p.doAlign);
        builder->openDataFiles(p.dirData+"/Run_"+p.runNum+"_DT_pos",p.dirData+"/Run_"+p.runNum+"_Si.root",p.outputFile);
        builder->matchEvents(p.nEvents);
        builder->dumpOutput();

        delete builder;
    }

    return 0;
}





