#pragma once
#include <vector>
#include "minimumpointsparameter.h"
#include "peakparameter.h"

using namespace std;

class OptimalModel
{
public:
	int mSG;
	int sg_window;
	vector<double> bindingEnergy;
	vector<double> spectrum;
	vector<double> smoothedSpectrum;
	vector<double> firstDerivative;
	vector<double> thirdDerivative;
	vector<double> background;
	vector<double> diffSpectrum;
	vector<MinimumPoints *> minimumPointsList;
	vector<MaximumPoints *> maximumPointsList;
	int idHigher;
	int idLower;
	vector<PeakParameter *> peakList;
	double nrmseSmoothedSpectrum;
	const vector<double> *originalBindingEnergy;
	const vector<double> *originalSpectrum;
	int idHigherGlobal;
	int idLowerGlobal;
	Settings settings;
};
