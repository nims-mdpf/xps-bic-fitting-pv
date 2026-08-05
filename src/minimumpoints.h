#pragma once
//#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <numeric>
#include "minimumpointsparameter.h"
#include "util.h"
#include "settings.h"
#include "output.h"

using namespace std;

int minimumPoints(const int ndata, const int m_SG, const vector<double> &bindingEnergy, const vector<double> &spectrum, const Settings &settings, vector<double> &smoothedSpectrum, vector<double> &firstDerivative, vector<double> &secondDerivative, vector<MaximumPoints*> &maximumPointsList, vector<MinimumPoints*> &minimumPointsList, vector<double> &sumSpectrum);
//int minimumPointsMain(const int ndata, const vector<double> &bindingEnergy, const vector<double> &spectrum, vector<MinimumPoints*> &minimumPointsList, const Settings &settings);
int savitzkyGolayFirst(const int ndata, const vector<double> &diffSpectrum, const int numSmooth, vector<double> &smoothedSpectrum, vector<double> &firstDerivative);
void adjustFirstDeriv(vector<double> &firstDerivative);
int detectMinimumPoints(const int ndata, const vector<double> &xx, const vector<double> smoothedSpectrum, vector<double> &firstDerivative, const Settings &settings, vector<MaximumPoints*> &maximumPointsList, vector<MinimumPoints*> &minimumPointsList);
void pickupMaximumPointsList(const vector<double> &xx, const vector<double> smoothedSpectrum, const vector<double> &firstDerivative, vector<MaximumPoints*> &maximumPointsList, const vector<MinimumPoints*> &minimumPointsList);
void pickupMinimumPointsList(const vector<double> &xx, const vector<double> smoothedSpectrum, const vector<double> &firstDerivative, vector<MinimumPoints*> &minimumPointsList);
int removeMinimumPointsByEnergyIntervalLimit(const vector<double> &bindingEnergy, vector<MinimumPoints*> &minimumPointsList, const vector<double> &smoothedSpectrum);
int removeMinimumPointsByHeightLimit(const vector<double> &bindingEnergy, vector<MaximumPoints*> &maximumPointsList, vector<MinimumPoints*> &minimumPointsList, const vector<double> &smoothedSpectrum);
int removeMinimumPointsByHeightLimit2(const vector<double> &bindingEnergy, vector<MaximumPoints*> &maximumPointsList, vector<MinimumPoints*> &minimumPointsList, const vector<double> &smoothedSpectrum);
