#pragma once

#include <vector>
#include "peakparameter.h"
#include "settings.h"
#include "util.h"

//#include "mainflow.h"

using namespace std;

int detectPeaks(const vector<double> &bindingEnergy, const vector<double> &smoothedSpectrum, vector<double> &thirdDerivative, const Settings &settings, vector<PeakParameter*> &peakList);
void adjustThirdDeriv(vector<double> &thirdDerivative);
void findExtremalIndexList(const vector<double> &smoothedSpectrum, const vector<double> &thirdDerivative, vector<int> &extremalIndexList);
void pickupPeakParameter(const vector<double> &xx, const vector<double> &smoothedSpectrum, const vector<int> &extremalIndexList, vector<PeakParameter*> &peakList, const Settings &settings);

