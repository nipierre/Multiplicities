#ifndef __PLOTMULT_H_INCLUDED__
#define __PLOTMULT_H_INCLUDED__

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <TLatex.h>

using namespace std;

struct Multiplicities { Double_t tab[2][3][4]; };

Multiplicities fMultiplicities[9][6][12];
Multiplicities fMultiplicities_yavg[9][12];

double z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};

//Graphic Style

int fMarkerColor[6] = {2,95,209,226,4,221};
int fMarkerStyle[6][2] = {{24,20},{26,22},{25,21},{27,33},{28,34},{30,29}};

#endif
