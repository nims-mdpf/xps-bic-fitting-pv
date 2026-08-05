#pragma once
#define _USE_MATH_DEFINES
#include <string.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
// #include <vector>

using namespace std;

class Settings
{
private:
	int readSettings();
	void setDefaultSettings();

public:
	Settings();
	~Settings() {};
	void setOptionalArgs(int argc, char *argv[]);
	void outputSettings();

	// int nSavitzkyGolaySmoothingPoint;
	double ratioPeakHeightThreshold;
	double bindingEnergyCutoffs[2];
	int maxLoopIterativeShirley;
	bool allowNegativeIterativeShirley;
	// bool allowHigherEdgeActive;
	// bool allowLowerEdgeActive;
	double lorentzRatioCutoffs[2];
	double minimumPeakFWHM;		 // shinotsuka20180228
	double minimumPeakAreaRatio; // shinotsuka20180530
	double allowedNegativeRatio; // shinotsuka20180601

	double spectrum_min;
	double spectrum_max;
	double spectrum_range;
	double bindingEnergy_min;
	double bindingEnergy_max;
	double bindingEnergy_range;
	string backgroundType; // "Shirley" (TODO: "Linear")
};
