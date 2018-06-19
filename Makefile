CC = gcc
CXX = g++
CCFLAGS = -g -O1 -W -Wall -Wno-unused-parameter -Wno-ignored-qualifiers #-pedantic -fPIC
ROOTFLAGS = `root-config --cflags --glibs`
ROOTVERSION = -D ROOT5
PHAST = /sps/compass/npierre/PHAST
PHAST_LIBS =
PHAST_INCL =
LHAPDF = /sps/compass/npierre/LHAPDF6
LHAPDF_INCL += -I$(LHAPDF)/include
LHAPDF_LIBS += -L$(LHAPDF)/lib -lLHAPDF

ifeq ($(DEBUG),1)
CCFLAGS += -DDEBUG
endif

all : analySIDIS acceptance comparison extractor plotter
analySIDIS : analySIDIS_split analySIDIS_collect
acceptance : accsplit acccollect
comparison : compMCRD compMCMC
extractor : FFExtractor
plotter : FFPlotter


%.o: %.cc %.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -c -o $@ $<

analySIDIS_split: analySIDIS_split.cc analySIDIS_split.h
	@echo 'Building SIDIS analysis package..'
	@$(CXX) $(CCFLAGS) -Wno-ignored-qualifiers $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(PHAST_LIBS) $(PHAST_INCL)

analySIDIS_collect: analySIDIS_collect.cc analySIDIS_collect.h
	@$(CXX) $(CCFLAGS) -Wno-ignored-qualifiers $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	@echo 'SIDIS analysis package built !'

accsplit: acceptance_split.cc acceptance_split.h
	@echo 'Building acceptance package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

acccollect: acceptance_collect.cc acceptance_collect.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	@echo 'Acceptance package built !'

compMCRD: compMCRD.cc compMCRD.h
	@echo 'Building RD/MC.MC/MC package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

compMCMC: compMCMC.cc compMCMC.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	@echo 'RD/MC.MC/MC package built !'

FFExtractor: FFExtractor.cc FFExtractor.h
	@echo 'Building FF extraction package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(LHAPDF_LIBS) $(LHAPDF_INCL)

FFPlotter: FFPlotter.cc FFPlotter.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(LHAPDF_LIBS) $(LHAPDF_INCL)
	@echo 'FF extraction package built !'

clean :
	@rm -rf *.o accsplit acccollect analySIDIS_split analySIDIS_collect compMCRD compMCMC FFExtractor FFPlotter
