#pragma once
#include <vector>

#include "settings.h"
#include "minimumpoints.h"
#include "calcinitialbackground.h"
#include "output.h"
#include "optimizeflow.h"
#include "optimalmodel.h"
#include "detectpeaks.h"
#include "util.h"
#include "mytimer.h"

#include <iostream>
#include <iomanip>
#include <thread>

using namespace std;

void mainFlow(const vector<double> &originalBindingEnergy,
			  const vector<double> &originalSpectrum, const Settings &settings, MyTimer &mytimer);
void optimizeFlow_wrapper(OptimalModel *om);
int get_numThreadLimit();
int get_omp_num_threads();
vector<int> make_smoothingList();
vector<int> make_pickupList(const int n);
vector<double> cutoffEnergyRange(const vector<double> &x, const int idS, const int idE);
void getGlobalEnergyRange(const vector<OptimalModel *> &optModelList, int &idHigherGlobal, int &idLowerGlobal);
void summarize_peakList(const vector<PeakParameter *> &peakList_Temp,
						double &min_width, double &max_width);
void getBoundary(vector<MinimumPoints *> &minimumPointsList, int &idHigher, int &idLower,
				 double &initBGHigher, double &initBGLower);
int hasNegative(const vector<double> &v);
void pickupModel(vector<OptimalModel *> &optModelList, const double nrmseThreshold = 0.08);
