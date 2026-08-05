#pragma once
// #include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <numeric>
#include "peakparameter.h"
#include "minimumpointsparameter.h"
#include "optimizeparameters.h"
#include "util.h"
#include "settings.h"

using namespace std;

int optimizeFlow(const int ndata, const int m_SG, vector<double> &bindingEnergy,
				 vector<double> &spectrum, vector<double> &background, const Settings &settings,
				 vector<PeakParameter *> &peakList,
				 vector<double> &sumSpectrum, double &chi2, double &aic, double &bic, double &reducedChi2,
				 double &bgHigher, double &bgLower,
				 double &scalingfactor, double &aic_gauss, double &bic_gauss, double &caic_poisson, double &caic_gauss,
				 bool &isConverged, int &numTotalOptimizeLoop, int &outerOptimizeLoop);

void normalize_spectrum_and_parameters(const int ndata, vector<double> &bindingEnergy, vector<double> &spectrum,
									   vector<double> &background, vector<PeakParameter *> &peakList,
									   vector<MinimumPoints *> &minimumPointsList, double &bgHigher, double &bgLower,
									   double &spectrum_min, double &spectrum_range, double &bindingEnergy_min, double &bindingEnergy_range);
void denormalize_spectrum_and_parameters(const int ndata, vector<double> &bindingEnergy, vector<double> &spectrum,
										 vector<double> &background, vector<PeakParameter *> &peakList,
										 vector<MinimumPoints *> &minimumPointsList, double &bgHigher, double &bgLower,
										 double &spectrum_min, double &spectrum_range, double &bindingEnergy_min, double &bindingEnergy_range);

int savitzkyGolay(const int ndata, const vector<double> &diffSpectrum, const int numSmooth, vector<double> &smoothedSpectrum, vector<double> &thirdDerivative);
void normalizeHeight(const int ndata, const vector<double> &diffSpectrum, const vector<double> &smoothedSpectrum, vector<PeakParameter *> &peakList);

int adaptiveSmoothing(const int numSmooth, const vector<double> &bindingEnergy, const vector<double> &diffSpectrum, vector<double> &smoothedSpectrum);
// double calc_negative_log_Likelyhood(const vector<double> &spectrum, const vector<double> &sumSpectrum,
//	const vector<double> &background, const int idHigher, const int idLower);
double calc_scalingfactor(const vector<double> &spectrum, const vector<double> &sumSpectrum,
						  const vector<double> &background);
// double calc_deviance_poisson(const vector<double> &spectrum, const vector<double> &sumSpectrum,
//	const vector<double> &background, const int idHigher, const int idLower, double &scalingfactor);
double calc_deviance_poisson(const vector<double> &spectrum, const vector<double> &sumSpectrum,
							 const vector<double> &background, double &scalingfactor, const bool do_scaling);
double calc_deviance_gauss(const vector<double> &spectrum, const vector<double> &sumSpectrum,
						   const vector<double> &background);
// double calc_deviance_poissongauss(const vector<double> &spectrum, const vector<double> &sumSpectrum,
//	const vector<double> &background, double &scalingfactor);
