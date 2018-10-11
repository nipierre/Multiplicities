#ifndef __COMPMULT_H_INCLUDED__
#define __COMPMULT_H_INCLUDED__

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <TLatex.h>

using namespace std;

struct Multiplicities { Double_t tab[2][3][4]; };

Multiplicities fMultiplicities1[9][6][12];
Multiplicities fMultiplicities2[9][6][12];
Multiplicities fMultiplicities1_yavg[9][12];
Multiplicities fMultiplicities2_yavg[9][12];
Multiplicities fMultiplicities1_zavg[9];
Multiplicities fMultiplicities2_zavg[9];


#endif
