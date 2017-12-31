CC = gcc
CXX = g++
CCFLAGS = -g -O1 -W -Wall -Wunused-parameter #-pedantic -fPIC
ROOTFLAGS = `root-config --cflags --glibs`
ROOTVERSION = -D ROOT5
PHAST = /sps/compass/npierre/PHAST
PHAST_LIBS =
PHAST_INCL =
LHAPDF = /sps/compass/npierre/LHAPDF6
LHAPDF_INCL += -I$(LHAPDF)/include
LHAPDF_LIBS += -L$(LHAPDF)/lib

ifeq ($(norc),1)
CCFLAGS += -DNORC
else
PHAST_INCL += -I$(PHAST)/include -lGeom -lMathMore $(shell cernlib kernlib mathlib)
PHAST_LIBS += -L$(PHAST)/lib -lPhast
endif

ifeq ($(debug),1)
CCFLAGS += -DDEBUG
endif

all : analySIDIS acceptance comparison extractor
analySIDIS : analySIDIS_split analySIDIS_collect
acceptance : accsplit acccollect
comparison : compMCRD
extractor : FFExtractor

%.o: %.cc %.h
	$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -c -o $@ $<

analySIDIS_split: analySIDIS_split.cc analySIDIS_split.h
	$(CXX) $(CCFLAGS) -Wno-ignored-qualifiers $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(PHAST_LIBS) $(PHAST_INCL)

analySIDIS_collect: analySIDIS_collect.cc analySIDIS_collect.h
	$(CXX) $(CCFLAGS) -Wno-ignored-qualifiers $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

accsplit: acceptance_split.cc acceptance_split.h
	$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

acccollect: acceptance_collect.cc acceptance_collect.h
	$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

compMCRD: compMCRD.cc compMCRD.h
	$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

FFExtractor: FFExtractor.cc FFExtractor.h
	$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(LHAPDF_LIBS) $(LHAPDF_INCL)

clean :
	rm -rf *.o accsplit acccollect analySIDIS_split analySIDIS_collect compMCRD FFExtractor
