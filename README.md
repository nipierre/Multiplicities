# Multiplicities Analysis

## Summary
 1. [Building](#Building)
 2. [Flow](#Flow)
 3. [Usage](#Usage)
 4. [analySIDIS_split](#analySIDIS-split)
 5. [analySIDIS_collect](#analySIDIS-collect)
 6. [acceptance_split](#acceptance-split)
 7. [acceptance_fuse](#acceptance-fuse)
 8. [acceptance_collect](#acceptance-collect)
 9. [compMult](#compMult)
 10. [compRDRD](#compRDRD)
 11. [compMCRD](#compMCRD)
 12. [compMCMC](#compMCMC)
 13. [DVM](#DVM)
 14. [FFExtractor](#FFExtractor)
 15. [FFPlotter](#FFPlotter)
 16. [Multiple Quick Studies Scripts](#quick-scripts)
 17. [Multiple plotting devices](#plotting-device)
 18. [More Infos](#more-infos)

## Building

 MAKEFILE PHONY TARGETS:
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


## [C++] analySIDIS_split<a name="analySIDIS-split" />

**Description:**
Takes the TTree and does the cut of the analysis. Outputs DIS Event and Hadron counts.

**Requires:**
 - **ROOT TTree from Phast User Event 20**
 - **RICH matrices in `data/rich_mat_2016.txt` and `data/rich_mat_error.txt`**
 - **Target description in `data/target-274508-274901.dat`**

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
 - **Count files in `$PWD/rawmult`**

## [C++] analySIDIS_collect<a name="analySIDIS-collect" />

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
 - **Multiplicity text files in `$data_path`:**
  - **`multiplicities_{}.txt` for Hadron, Pion, Kaon and Proton**
  - **`multiplicities_{}_yavg.txt` for Hadron, Pion, Kaon and Proton**
  - **`multiplicities_raw.txt`**
  - **`multiplicities_h{}.txt` for + and -**
  - **`multiplicities_hadron_{}.txt` for pt and theta**
  - **`reldiff.txt`**
 - **Multiplicity plots in `data_path`:**
  - **`{}_multiplicity_file.pdf` for Hadron, Pion, Kaon and Proton**
  - **`{}_multiplicity_zvtx_file.pdf` for Hadron, Pion, Kaon and Proton**
  - **`{}_multiplicity_yavg_file.pdf` for Hadron, Pion, Kaon and Proton**
  - **`{}_multiplicity_sum_file.pdf` for Hadron, Pion, Kaon and Proton**
  - **`{}_multiplicity_ratio_file.pdf` for Hadron, Pion, Kaon and Proton**

## [C++] acceptance_split<a name="acceptance-split" />


**Description:**
Takes the Monte Carlo files and does the cut of the analysis. Outputs DIS event and Hadron counts.

**Requires:**
 - **Target description in `data/target-274508-274901.dat`**

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
 - **Count files in `$PWD/acceptance`**

**NB: Separate mu+/- charge when treating**

## [C++] acceptance_fuse<a name="acceptance-fuse" />


**Description:**
Fuses mu+/- counts data per period

**Requires:**
 - **Outputs from `acceptance_split`**
 - **Target description in `data/target-274508-274901.dat`**

**User Dependence**
 - **`data_path` to MC data @ line 14**

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
 - **Fused count files in `$PWD/acceptance`**

## [C++] acceptance_collect<a name="acceptance-collect" />


**Description:**
Takes the output of `analySIDIS_split`, computes the Multiplicities and stores/plots them.

**Requires:**
 - **Output of `analySIDIS_split`**
 - **Output from Acceptance framework**
 - **Inclusive and Semi-inclusive Radiative Corrections**
 - **Diffractive Vector Meson Correction**

**User Dependence:**
- **`dirroot` to acceptance counts @ line 43**

**In File Flags**
 - **Y STAGGERING (SPREAD) [YES/NO | 1/0]**

**Call:**
 ```Bash
 ./acceptance_collect [PeriodFile]
 ```
 where `[PeriodFile]` is a file containing the periods to treat.

**Outputs:**
 - **Acceptance text files in `$dirroot/acceptance/$YEAR`:**
  - **`acceptance_{}.txt` per Period**
  - **`acceptance_yavg_{}.txt` per Period**
  - **`acceptance_vtx_{}.txt` per Period**
  - **`acceptance_theta_{}.txt` per Period**
  - **`reldiff_vtx_{}.txt` per Period**
 - **Acceptance plots in `$dirroot/acceptance/$YEAR`:**
  - **`{}_acceptance_{}.pdf` for Hadron, Pion, Kaon, Proton and per Period**
  - **`{}_acceptance_corr_{}.pdf` for Hadron, Pion, Kaon, Proton and per Period**

## [C++] compMult<a name="compMult" />

**Call:**
```Bash
./compMult [MULT_FILE_1] [MULT_FILE_2] [CUTFILE]
```

where `[CUTFILE]` is a file with kinematical cuts. The structure is the following:

```
Xmin  0.004
Xmax  0.4
Ymin  0.1
Ymax  0.9
Wmin  5
Wmax  17
Pmin  3
Pmax  40
```

## [C++] compRDRD<a name="compRDRD" />

**Call:**
```Bash
./compRDRD [RD_FILELIST_1] [RD_FILELIST_2] [CUTFILE]
```

where `[CUTFILE]` is a file with kinematical cuts.

## [C++] compMCRD<a name="compMCRD" />

**Call:**
```Bash
./compMCRD [RD_FILELIST] [MC_FILELIST] [CUTFILE]
```

where `[CUTFILE]` is a file with kinematical cuts.

## [C++] compMCMC<a name="compMCMC" />

**Call:**
```Bash
./compMCMC [MC_FILELIST_1] [MC_FILELIST_2] [CUTFILE]
```

where `[CUTFILE]` is a file with kinematical cuts.

## [C++] DVM<a name="DVM" />

**Call:**
```Bash
./DVM [SIDIS_FILELIST] [RHO_FILELIST] [PHI_FILELIST] [CUTFILE]
```

where `[CUTFILE]` is a file with kinematical cuts.

## [C++] FFExtractor<a name="FFExtractor" />

**Call:**
```Bash
./FFExtractor [OPTIONS]
```
with `[OPTIONS]` being:
 - ```-pion-deut [PI+_FILE] [PI-_FILE]```
 - ```-pion-prot [PI+_FILE] [PI-_FILE]```
 - ```-kaon-3 [K+_PROT] [K-_PROT] [K+_DEUT] [K-_DEUT]```
 - ```-kaon-4 [K+_PROT] [K-_PROT] [K+_DEUT] [K-_DEUT]```
 - ```-dummy-data [MULT_BASE_FILE]```

## [C++] FFPlotter<a name="FFPlotter" />

**Call:**
```Bash
./FFPlotter [OPTIONS]
```
with `[OPTIONS]` being:
 - ```-pion [PI_FILE]```
 - ```-pion-next [PI_FILE]```
 - ```-kaon [K_FILE]```
 - ```-kaon-next [K_FILE]```

## Multiple Quick Studies Scripts<a name="quick-scripts" />
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

## Multiple plotting devices<a name="plotting_device" />
 - [C] PlotAccComp
 - [C] PlotDVM
 - [C++] plotMult
 - [C] PlotMultComp
 - [C] PlotRC
 - [C] FitElectron
 - [C] ElectronMacro

## More Infos<a name="more-infos" />
 - [Julia](https://julialang.org)
 - [ROOT](http://ROOT.cern.ch)
