# Multiplicities Analysis

## Building

 - MAKEFILE PHONY TARGETS:
  - `make`: builds all targets
  - `make analySIDIS`: builds the SIDIS analysis framework.
  - `make acceptance`: builds the acceptance analysis framework
  - `make comparison`: builds several comparison applets
  - `make extractor`: builds the Fragmentation Functions extraction applet
  - `make plotter`: builds the Fragmentation Functions plotter applet
  - `make dvm`: builds the Diffractive Vector Mesons applet
  - `make clean`

## Flow

**Acceptance Calculation** _[acceptance_split --> acceptance_fuse --> acceptance_collect]_ **-->** **Multiplicities Calculation** _[analySIDIS_split --> analySIDIS_collect]_

## Usage

### [C++] analySIDIS_split

**Description:**
Takes the TTree and does the cut of the analysis. Outputs DIS Event and Hadron counts.

**Requires:**
 - **ROOT TTree from Phast User Event 20**
 - **RICH matrices in data/rich_mat_2016.txt and data/rich_mat_error.txt**
 - **Target description in data/target-274508-274901.dat**

**User Dependence:**
 - **data_path to Real Data @ line 25**

**In File Flags:**
 - **Momentum Boundaries**
 - **X Boundaries**
 - **Y Boundaries**
 - **Muon charge separation [YES/NO | 1/0]**

**Call:**
```Bash
./analySIDIS_split [PeriodFile] [OPTIONAL FLAG]
```
where `[PeriodFile]` is a file containing the periods to treat. It has the following structure:

```
P01 0
P02 0
P03 0
P04 0
P05 0
P06 0
P07 1
P08 1
P09 1
P10 1
P11 1
```

`[OPTIONAL FLAG]`
 - `-k` : draw kinematic plots

**Outputs:**
 - **Count files in $PWD/rawmult**

### [C++] analySIDIS_collect

**Description:**
Takes the output of `analySIDIS_split`, computes the Multiplicities and stores/plots them.

**Requires:**
 - **Output of `analySIDIS_split`**
 - **Output from Acceptance framework**
 - **Inclusive and Semi-inclusive Radiative Corrections**
 - **Diffractive Vector Meson Correction**

**In File Flags**
 - **Diffraction Vector Meson [YES/NO | 1/0]**
 - **Radiative Corrections [YES/NO | 1/0]**
 - **No Acceptance [YES/NO | 1/0]**
 - **Y Integration or mean [Mean/Weighted Mean/Integration | 1/2/3]**
 - **Staggered Multiplicities in plot [YES/NO | 1/0]**

**Call:**
 ```Bash
 ./analySIDIS_collect [PeriodFile]
 ```
 where `[PeriodFile]` is a file containing the periods to treat.

**Outputs:**
 - **Multiplicity text files in $data_path:**
  - **multiplicities_{}.txt for Hadron, Pion, Kaon and Proton**
  - **multiplicities_{}_yavg.txt for Hadron, Pion, Kaon and Proton**
  - **multiplicities_raw.txt**
  - **multiplicities_h{}.txt for + and -**
  - **multiplicities_hadron_{}.txt for pt and theta**
  - **reldiff.txt**
 - **Multiplicity plots in data_path:**
  - **{}_multiplicity_file.pdf for Hadron, Pion, Kaon and Proton**
  - **{}_multiplicity_zvtx_file.pdf for Hadron, Pion, Kaon and Proton**
  - **{}_multiplicity_yavg_file.pdf for Hadron, Pion, Kaon and Proton**
  - **{}_multiplicity_sum_file.pdf for Hadron, Pion, Kaon and Proton**
  - **{}_multiplicity_ratio_file.pdf for Hadron, Pion, Kaon and Proton**

## [C++] acceptance_split

**Description:**
Takes the Monte Carlo files and does the cut of the analysis. Outputs DIS event and Hadron counts.

**Requires:**
 - **Target description in data/target-274508-274901.dat**

**User Dependence**
 - **data_path to MC data @ line 14**

**In File Flags:**
 - **Momentum Boundaries**
 - **X Boundaries**
 - **Y Boundaries**
 - **W Boundaries**
 - **XX0 Limit**

**Call:**
```Bash
./acceptance_split [PeriodFile] [OPTIONAL FLAG]
```
where `[PeriodFile]` is a file containing the periods to treat.

`[OPTIONAL FLAG]`
 - `-k` : draw kinematic plots

**Outputs:**
 - **Count files in $PWD/acceptance**

**NB: Separate mu+/- charge when treating**

## [C++] acceptance_fuse

**Description:**
Fuses mu+/- counts data per period

**Requires:**
 - **Outputs from `acceptance_split`**
 - **Target description in data/target-274508-274901.dat**

**User Dependence**
 - **data_path to MC data @ line 14**

**In File Flags:**
 - **Momentum Boundaries**
 - **X Boundaries**
 - **Y Boundaries**
 - **W Boundaries**
 - **XX0 Limit**

**Call:**
```Bash
./acceptance_fuse [PeriodName] [MU+ FILELIST] [MU- FILELIST]
```
where `[PeriodName]` is the name of the period eg. `P07`.

**Outputs:**
 - **Fused count files in $PWD/acceptance**

## [C++] acceptance_collect

**Description:**
Takes the output of `analySIDIS_split`, computes the Multiplicities and stores/plots them.

**Requires:**
 - **Output of `analySIDIS_split`**
 - **Output from Acceptance framework**
 - **Inclusive and Semi-inclusive Radiative Corrections**
 - **Diffractive Vector Meson Correction**

**User Dependence:**
- **dirroot to acceptance counts @ line 43**

**In File Flags**
 - **Y STAGGERING (SPREAD) [YES/NO | 1/0]**

**Call:**
 ```Bash
 ./acceptance_collect [PeriodFile]
 ```
 where `[PeriodFile]` is a file containing the periods to treat.

**Outputs:**
 - **Acceptance text files in `$dirroot/acceptance/$YEAR`:**
  - **acceptance_{}.txt per Period**
  - **acceptance_yavg_{}.txt per Period**
  - **acceptance_vtx_{}.txt per Period**
  - **acceptance_theta_{}.txt per Period**
  - **reldiff_vtx_{}.txt per Period**
  - **reldiff.txt**
 - **Acceptance plots in `$dirroot/acceptance/$YEAR`:**
  - **{}_acceptance_{}.pdf for Hadron, Pion, Kaon, Proton and per Period**
  - **{}_acceptance_corr_{}.pdf for Hadron, Pion, Kaon, Proton and per Period**

## [C++] compMult

## [C++] compRDRD

## [C++] compMCRD

## [C++] compMCMC

## [C++] DVM

## [C++] FFExtractor

## [C++] FFPlotter

## Multiple Quick Studies Scripts
 - [Julia] RichStudy
 - [Julia] TargetStudy
 - [Julia] RadiativeCorrectionFactors
 - [Julia] MultPredictions
 - [Julia] MultStudies
 - [Julia] MultVertexed
 - [Julia] MultXCheck
 - [Julia] Import2006HEP
 - [Julia] HadronCount
 - [Julia] DVMXC
 - [Julia] DVMComparison

## Multiple plotting devices
 - [C] PlotAccComp
 - [C] PlotDVM
 - [C++] plotMult
 - [C] PlotMultComp
 - [C] PlotRC
 - [C] FitElectron
 - [C] ElectronMacro

## More Infos
 - [Julia](https://julialang.org)
 - [ROOT](http://ROOT.cern.ch)
