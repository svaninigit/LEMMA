# Makefile for the ROOT RADMU programs.
# Expectation Maximization - EM
# Author: Sara Vanini 10/6/2008

ObjSuf        = o
SrcSuf        = cxx
ExeSuf        = x
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

# ROOT
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

# Linux
CXX           = g++
#flags for debugging (e.g. with Valgrind)
#CXXFLAGS      = -g -O0 -O -Wno-deprecated -Wall -fPIC -I$(ROOTSYS)/include
#flags for maximum optimization (but no debugging....)
CXXFLAGS      = -O3 -Wno-deprecated -Wall -fPIC -I$(ROOTSYS)/include
CXXFLAGS     += $(ROOTCFLAGS)
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -Wl,-soname,libEvent.so -shared
LIBS          = $(ROOTLIBS) -lgcc -lm -ldl -rdynamic
#GLIBS         = $(ROOTLIBS) $(ROOTGLIBS) -L/usr/X11R6/lib \
#                -lXpm -lX11 -lgcc -lm -ldl -rdynamic
GLIBS         = $(ROOTGLIBS)

#------------------------------------------------------------------------------

OBJS = main.o TOMMainFrame.o  Track_IO.o TTrigCalibration.o Save_HistosAndTree.o RawHistos.o HITCollection.o HIT.o HITColl_Seg.o HITColl_Layer.o Segment.o Track.o FIT.o TimeCorr.o Geom.o RawAnalyzer.o


runTOM: ${OBJS}
	${CXX} -o runTOM ${CXXFLAGS} ${OBJS} $(LIBS) $(GLIBS)

# ---- EM classes -----
Voxel.o: Voxel.C Voxel.h
	${CXX} ${CXXFLAGS} -c Voxel.C

Muon.o: Muon.C Muon.h
	${CXX} ${CXXFLAGS} -c Muon.C

VoxCollection.o: VoxCollection.C VoxCollection.h
	${CXX} ${CXXFLAGS} -c VoxCollection.C

MuonCollection.o: MuonCollection.C MuonCollection.h
	${CXX} ${CXXFLAGS} -c MuonCollection.C

# ---- PattRec classes -----
HITCollection.o: PattRec/HITCollection.C PattRec/HITCollection.h
	${CXX} ${CXXFLAGS} -c PattRec/HITCollection.C

HIT.o: PattRec/HIT.C PattRec/HIT.h
	${CXX} ${CXXFLAGS} -c PattRec/HIT.C

HITColl_Seg.o: PattRec/HITColl_Seg.C PattRec/HITColl_Seg.h
	${CXX} ${CXXFLAGS} -c PattRec/HITColl_Seg.C

Segment.o: PattRec/Segment.C PattRec/Segment.h
	${CXX} ${CXXFLAGS} -c PattRec/Segment.C

Track.o: PattRec/Track.C PattRec/Track.h
	${CXX} ${CXXFLAGS} -c PattRec/Track.C

HITColl_Layer.o: PattRec/HITColl_Layer.C PattRec/HITColl_Layer.h
	${CXX} ${CXXFLAGS} -c PattRec/HITColl_Layer.C

FIT.o: PattRec/FIT.C PattRec/FIT.h
	${CXX} ${CXXFLAGS} -c PattRec/FIT.C

TimeCorr.o: PattRec/TimeCorr.C PattRec/TimeCorr.h
	${CXX} ${CXXFLAGS} -c PattRec/TimeCorr.C

Geom.o: PattRec/Geom.C PattRec/Geom.h
	${CXX} ${CXXFLAGS} -c PattRec/Geom.C

RawAnalyzer.o: PattRec/RawAnalyzer.C PattRec/RawAnalyzer.h
	${CXX} ${CXXFLAGS} -c PattRec/RawAnalyzer.C

Save_HistosAndTree.o: PattRec/Save_HistosAndTree.C PattRec/Save_HistosAndTree.h
	${CXX} ${CXXFLAGS} -c PattRec/Save_HistosAndTree.C

#RawHistos.o: RawHistos.C RawHistos.h
#        ${CXX} ${CXXFLAGS} -c RawHistos.C

TTrigCalibration.o: Monitor/TTrigCalibration.C Monitor/TTrigCalibration.h
	${CXX} ${CXXFLAGS} -c Monitor/TTrigCalibration.C

Track_IO.o: PattRec/Track_IO.C PattRec/Track_IO.h
	${CXX} ${CXXFLAGS} -c PattRec/Track_IO.C



# ---- TOM classes -----
TOMMainFrame.o: TOMMainFrame.C TOMMainFrame.h
	${CXX} ${CXXFLAGS} -c TOMMainFrame.C

# ---- Monitor classes -----
RawHistos.o: Monitor/RawHistos.C Monitor/RawHistos.h
	${CXX} ${CXXFLAGS} -c Monitor/RawHistos.C

#CalibHistos.o: Monitor/CalibHistos.C Monitor/CalibHistos.h
#	${CXX} ${CXXFLAGS} -c Monitor/CalibHistos.C

# ---- THE MAIN -----
main.o: main.c
	${CXX} ${CXXFLAGS} -c main.c

clean:
	rm -f runTOM ${OBJS}
	@echo "all cleaned up!"
