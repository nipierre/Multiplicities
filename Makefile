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
acceptance : accsplit accfuse acccollect
comparison : compRDRD compMCRD compMCMC compMult
extractor : FFExtractor FFPlotter
plotter : plotMult


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

accfuse: acceptance_fuse.cc acceptance_fuse.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

acccollect: acceptance_collect.cc acceptance_collect.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	@echo 'Acceptance package built !'

compRDRD: compRDRD.cc compRDRD.h
	@echo 'Building RD/RD.RD/MC.MC/MC.Mult/Mult package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

compMCRD: compMCRD.cc compMCRD.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

compMCMC: compMCMC.cc compMCMC.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

compMult: compMult.cc compMult.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	@echo 'RD/RD.RD/MC.MC/MC.Mult/Mult package built !'

FFExtractor: FFExtractor.cc FFExtractor.h
	@echo 'Building FF extraction package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(LHAPDF_LIBS) $(LHAPDF_INCL)

FFPlotter: FFPlotter.cc FFPlotter.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(LHAPDF_LIBS) $(LHAPDF_INCL)
	@echo 'FF extraction package built !'

plotMult: plotMult.cc plotMult.h
	@echo 'Building plotting device package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(LHAPDF_LIBS) $(LHAPDF_INCL)
	@echo 'Plotting device package built !'

clean :
	@rm -rf *.o accsplit accfuse acccollect analySIDIS_split analySIDIS_collect compRDRD compMCRD compMCMC compMult FFExtractor FFPlotter plotMult
