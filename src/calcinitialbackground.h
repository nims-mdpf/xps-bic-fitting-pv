#pragma once
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "settings.h"
#include "minimumpointsparameter.h"
#include "optimizeparameters.h" // To use make_background_Linear()

using namespace std;

//int calcInitialBackground(const vector<MinimumPoints*> &minimumPointsList, const int ndata,
//	const vector<double> &spectrum, const vector<double> &bindingEnergy,
//	vector<vector<double>> &backgroundList, vector<double> &background,
//	vector<double> &diffSpectrum, const Settings &settings, const vector<double> &smoothedSpectrum);
//
//int calcInitialBackgroundShirley(const vector<MinimumPoints*> &minimumPointsList, const int ndata,
//	const vector<double> &spectrum, const vector<double> &bindingEnergy,
//	vector<vector<double>> &backgroundList, vector<double> &background,
//	vector<double> &diffSpectrum, const Settings &settings, const double initBGLower);

int calcInitialBackground(const int ndata, const vector<double> &bindingEnergy,
	const vector<double> &spectrum, const vector<double> &smoothedSpectrum,
	const vector<MinimumPoints*> &minimumPointsList, const Settings &settings,
	vector<double> &background, vector<double> &diffSpectrum);

int calcInitialBackgroundShirley(const int ndata, const vector<double> &bindingEnergy,
	const vector<double> &spectrum, const vector<MinimumPoints*> &minimumPointsList,
	const Settings &settings, vector<double> &background,
	vector<double> &diffSpectrum, const double initBGLower);

int calcInitialBackgroundLinear(const int ndata,
	const vector<double> &spectrum, const vector<double> &smoothedSpectrum,
	const Settings &settings, vector<double> &background, vector<double> &diffSpectrum);

int calcInitialBackgroundLinear2(const int ndata, const vector<double> &bindingEnergy,
	const vector<double> &spectrum, const vector<double> &smoothedSpectrum,
	const vector<MinimumPoints*> &minimumPointsList, const Settings &settings,
	vector<double> &background, vector<double> &diffSpectrum);
