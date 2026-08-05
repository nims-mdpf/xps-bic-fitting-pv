#pragma once
#define _USE_MATH_DEFINES
//#include <vector>

//#include "util.h"
#include "settings.h"
#include "minimumpoints.h"

int output_result(const string &kwd, vector<double> bindingEnergy, vector<double> spectrum,
	vector<double> background, vector<PeakParameter*> peakList,
	vector<double> originalBindingEnergy, vector<double> originalSpectrum,
	const int idHigher, const int idLower, vector<double> smoothedSpectrum,
	const Settings &settings);
int outputBackgroundList(const vector<double> &bindingEnergy, const vector<double> &spectrum, const vector<MinimumPoints*> &minimumPointsList, const vector<vector<double>> &backgroundList);
void outputBackground(const vector<double> &bindingEnergy, const vector<double> &spectrum, const vector<double> &background, const vector<double> &diffSpectrum);
void outputSavitzkyGolay(const vector<double> &bindingEnergy, const vector<double> &diffSpectrum,
	const vector<double> &smoothedSpectrum, const vector<double> &thirdDerivative);
void outputMaximumPoints(const string &kwd, const vector<double> &bindingEnergy, const vector<double> &spectrum, vector<MaximumPoints*> &maximumPointsList);
void outputMinimumPoints(const string &kwd, const vector<double> &bindingEnergy, const vector<double> &originalSpectrum, vector<MinimumPoints*> &minimumPointsList);

