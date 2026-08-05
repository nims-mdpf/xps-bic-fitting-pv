#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include "peakparameter.h"
#include "minimumpointsparameter.h"
#include "util.h"
#include <iostream>
#include <algorithm>
#include "settings.h"

#define MYLIBDIR "C:\\mylibs\\clapack_vs2015_x64"
#pragma comment(lib, MYLIBDIR "\\lapack.lib")
#pragma comment(lib, MYLIBDIR "\\blas.lib")
#pragma comment(lib, MYLIBDIR "\\libf2c.lib")
extern "C"
{
	void dgesv_(const int &N, const int &NRHS,
				double **A, const int &LDA, int *IPIV,
				double **B, const int &LDB, int &INFO);
};

using namespace std;

void calculateMatrixAtA(vector<vector<double>> &a, vector<double> &s, vector<vector<double>> &ata, vector<double> &ats);
void calculateMatrixAtAD(const double lambda, vector<vector<double>> &ata);
void solveMatrixEquation(vector<vector<double>> &ata, vector<double> &ats, vector<double> &e);
void restrictPeaks(vector<PeakParameter *> &peakList, const Settings &settings);
int removeInfinitesimalPeaks(vector<PeakParameter *> &peakList);
int removeSmallPeaks(vector<PeakParameter *> &peakList);
int removeSmallPeaksByWholeArea(vector<PeakParameter *> &peakList, const double threshold);
int removeSmallPeaksByPeakWidth(vector<PeakParameter *> &peakList, const double energy_step, const double threshold_fwhm);
int removePeaksByEnergyRange(const vector<double> &bindingEnergy, vector<PeakParameter *> &peakList);
int removePeaksByHWHM(vector<PeakParameter *> &peakList, const double threshold_width);
// int make_background(const double &bgHigher, const double &bgLower,
//	const vector<double> &ffp, vector<double> &qpq, vector<double> &background);
int make_background(const double &bgHigher, const double &bgLower,
					const vector<double> &ffp, const Settings &settings,
					vector<double> &qpq, vector<double> &background);
int make_background_Shirley(const double &bgHigher, const double &bgLower,
							const vector<double> &ffp, vector<double> &qpq, vector<double> &background);
int make_background_Linear(const double &yH, const double &yL, vector<double> &background);
int calcChiSquare(const vector<double> &x, const vector<double> &spectrum,
				  const vector<PeakParameter *> &peakList,
				  const double &bgHigher, const double &bgLower, const Settings &settings,
				  vector<double> &background, vector<double> &qpq, vector<double> &s,
				  double &newChiSquare, double &newReducedChiSquare);
// void outputChiSquare(const vector<double> &x, const vector<double> &spectrum, const vector<double> &s);
void outputChiSquare(const vector<double> &x, const vector<double> &spectrum, const vector<double> &ffp, const vector<double> &bg, const vector<double> &s);

int output_history(const int m_SG, vector<double> x, vector<double> y, vector<double> background, vector<PeakParameter *> peakList, const Settings &settings);

void calculateMatrixA_for_peak(const int ndata, const vector<double> &x, const vector<PeakParameter *> &peakList, vector<vector<double>> &a);
void calculateMatrixA_for_background_Shirley(const int ndata, const vector<double> &qpq, vector<vector<double>> &a);
void calculateMatrixA_for_background_Linear(vector<vector<double>> &a);

void copy_peakparameter(const vector<PeakParameter *> &peakList, const double bgHigher, const double bgLower,
						vector<PeakParameter *> &peakList_saved, double &bgHigher_saved, double &bgLower_saved);
int optimizeParameters(const int ndata, const int m_SG, const vector<double> &x, const vector<double> &spectrum,
					   const bool isPeakOptimize, const bool isBgOptimize, const bool isFixedLorentzRatio, const double fixedLorentzRatio,
					   const Settings &settings, vector<double> &background, vector<PeakParameter *> &peakList,
					   bool isNumberPeakFixed, int loopMaxOptimisze,
					   double &optimizedChiSquare, double &optimizedReducedChiSquare,
					   double &bgHigher, double &bgLower, int &loop);

class Spectrum
{
public:
	int ndata;
	vector<double> x;
	vector<double> y;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	void setdata(const vector<double> &xi, const vector<double> &yi);
	// PeakParameter(double _h, double _p, double _s, double _g) :
	// 	height(_h), position(_p), sigma(_s), gamma(_g) {}
	// ~PeakParameter() {}
};

extern Spectrum spect;
