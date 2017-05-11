# Makefile for the ROOT RADMU programs.
# Pattern Recognition - PR
# Author: Silvia Pesente march 2009

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
CXX           = g++ -m64
#flags for debugging with Valgrind
#CXXFLAGS      = -g -O -O -Wall -fPIC -I$(ROOTSYS)/include
CXXFLAGS      = -O2 -Wall -Wno-deprecated -fPIC -I$(ROOTSYS)/include
CXXFLAGS     += $(ROOTCFLAGS)
LD            = g++
#LDFLAGS       = -g
SOFLAGS       = -Wl,-soname,libEvent.so -shared
LIBS          = $(ROOTLIBS) -lgcc -lm -ldl -rdynamic
#GLIBS         = $(ROOTLIBS) $(ROOTGLIBS) -L/usr/X11R6/lib \
#                -lXpm -lX11 -lgcc -lm -ldl -rdynamic
GLIBS         = $(ROOTGLIBS)

#------------------------------------------------------------------------------

OBJS = HITCollection.o HITColl_Seg.o HITColl_Layer.o HIT.o Segment.o Track.o FIT.o TimeCorr.o Geom.o RawAnalyzer.o Save_HistosAndTree.o RawHistos.o TTrigCalibration.o Track_IO.o main.o

#OBJS = HITCollection.o HITColl_Seg.o HITColl_Layer.o HIT.o Segment.o Track.o FIT.o TimeCorr.o Geom.o RawAnalyzer.o Save_HistosAndTree.o RawHistos.o ImgAnalyzer.o Muon.o VoxCollection.o MuonCollection.o TTrigCalibration.o Track_IO.o main.o

#runPR: ${OBJS} Constants.h Debugs.h ../EM/header.h
runPR: ${OBJS} Constants.h Debugs.h
	${CXX} -o runPR ${CXXFLAGS} ${OBJS} $(LIBS) $(GLIBS)

HIT.o: HIT.C HIT.h
	${CXX} ${CXXFLAGS} -c HIT.C

HITCollection.o: HITCollection.C HITCollection.h
	${CXX} ${CXXFLAGS} -c HITCollection.C

HITColl_Seg.o: HITColl_Seg.C HITColl_Seg.h
	${CXX} ${CXXFLAGS} -c HITColl_Seg.C

Segment.o: Segment.C Segment.h
	${CXX} ${CXXFLAGS} -c Segment.C

Track.o: Track.C Track.h
	${CXX} ${CXXFLAGS} -c Track.C

HITColl_Layer.o: HITColl_Layer.C HITColl_Layer.h
	${CXX} ${CXXFLAGS} -c HITColl_Layer.C

FIT.o: FIT.C FIT.h
	${CXX} ${CXXFLAGS} -c FIT.C

TimeCorr.o: TimeCorr.C TimeCorr.h
	${CXX} ${CXXFLAGS} -c TimeCorr.C

Geom.o: Geom.C Geom.h
	${CXX} ${CXXFLAGS} -c Geom.C

RawAnalyzer.o: RawAnalyzer.C RawAnalyzer.h
	${CXX} ${CXXFLAGS} -c RawAnalyzer.C

Save_HistosAndTree.o: Save_HistosAndTree.C Save_HistosAndTree.h
	${CXX} ${CXXFLAGS} -c Save_HistosAndTree.C

#ImgAnalyzer.o: ../EM/ImgAnalyzer.C ../EM/ImgAnalyzer.h
#	${CXX} ${CXXFLAGS} -c ../EM/ImgAnalyzer.C

#Muon.o: ../EM/Muon.C ../EM/Muon.h
#	${CXX} ${CXXFLAGS} -c ../EM/Muon.C

#VoxCollection.o: ../EM/VoxCollection.C ../EM/VoxCollection.h
#	${CXX} ${CXXFLAGS} -c ../EM/VoxCollection.C

#MuonCollection.o: ../EM/MuonCollection.C ../EM/MuonCollection.h
#	${CXX} ${CXXFLAGS} -c ../EM/MuonCollection.C

RawHistos.o: ../Monitor/RawHistos.C ../Monitor/RawHistos.h
	${CXX} ${CXXFLAGS} -c ../Monitor/RawHistos.C

TTrigCalibration.o: ../Monitor/TTrigCalibration.C ../Monitor/TTrigCalibration.h
	${CXX} ${CXXFLAGS} -c ../Monitor/TTrigCalibration.C

Track_IO.o: Track_IO.C Track_IO.h
	${CXX} ${CXXFLAGS} -c Track_IO.C

main.o: main.c
	${CXX} ${CXXFLAGS} -c main.c

clean:
	rm -f myprogram ${OBJS}
	@echo "all cleaned up!"
