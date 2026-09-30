#include "mainflow.h"
#include <algorithm>
#include <omp.h>

/*
 * ---------------------------------------------------------------------
 * Main flow of this program
 * There are two major loops for the round-robin smoothing.
 * 1st loop:
 *   - Savitzky-Golay smoothing
 *   - Detect minimum points
 *   - Remove unnecessary minimum points
 *   - Get a boundary of energy range
 *   - Make initial background spectrum by the iterative Shirley method
 *   - Detect initial peaks
 *
 * After 1st loop, get a global boundary of energy range
 *
 * 2nd loop:
 *   - Optimize the peaks and background
 * ---------------------------------------------------------------------
 * INPUT
 *    originalBindingEnergy: Original binding energies (scaled)
 *    originalSpectrum:      Original spectrum (scaled)
 *    settings:              Settings
 * ---------------------------------------------------------------------
 * H. Shinotsuka 2019/08/08
 * ---------------------------------------------------------------------
 */
void mainFlow(const vector<double> &originalBindingEnergy,
			  const vector<double> &originalSpectrum, const Settings &settings, MyTimer &mytimer)
{
	cout << "mainFlow()" << endl;
	cout << endl;
	cout << "#=========================================================#" << endl;
	cout << "#  During round-robin smoothing, make initial background, #" << endl;
	cout << "#  initial peaks, then get global energy range            #" << endl;
	cout << "#=========================================================#" << endl;
	cout << endl;

	// Number of data points
	const int originalNdata = (int)originalBindingEnergy.size();

	// shinotsuka20181204, Tune the strength of smoothing
	vector<int> mSGList = make_smoothingList();
	int num_mSGList = (int)mSGList.size();

	vector<int> pickupList = make_pickupList(num_mSGList);

	stringstream ss;

	ofstream ofs2;
	ofs2.open("out_auto_smoothing_minpoints.csv");
	ofs2 << "m_SG,min_points_position(eV),min_points_height" << endl;

	vector<OptimalModel *> optModelList;

	vector<double> smoothedSpectrum;
	smoothedSpectrum.resize(originalNdata);
	smoothedSpectrum = originalSpectrum;

	for (int m_SG = 1; m_SG <= num_mSGList; m_SG++)
	{
		int sg_window = mSGList[m_SG - 1];
		cout << endl
			 << m_SG << " / " << num_mSGList << " Round-robin for the smoothing number. "
			 << "sg_window= " << sg_window << endl;

		vector<double> smoothedSpectrum_copy = smoothedSpectrum;
		savitzkyGolayGeneral(smoothedSpectrum_copy, sg_window, 0, smoothedSpectrum);

		if (pickupList[m_SG - 1] == 0)
		{
			continue;
		}

		vector<double> firstDerivative;
		vector<double> thirdDerivative;
		firstDerivative.resize(originalNdata);
		thirdDerivative.resize(originalNdata);

		savitzkyGolayGeneral(smoothedSpectrum_copy, sg_window, 1, firstDerivative);
		savitzkyGolayGeneral(smoothedSpectrum_copy, sg_window, 3, thirdDerivative);

		OptimalModel *optModel = new OptimalModel();
		optModel->mSG = m_SG;
		optModel->sg_window = sg_window;
		optModel->bindingEnergy = originalBindingEnergy;
		optModel->spectrum = originalSpectrum;
		optModel->smoothedSpectrum = smoothedSpectrum;
		optModel->firstDerivative = firstDerivative;
		optModel->thirdDerivative = thirdDerivative;

		vector<MinimumPoints *> minimumPointsList;
		vector<MaximumPoints *> maximumPointsList;

		// Detect minimum points from a smoothed spectrum, in order to achieve the "wise" optimization.
		detectMinimumPoints(originalNdata, originalBindingEnergy, smoothedSpectrum, firstDerivative, settings, maximumPointsList, minimumPointsList);

		copyVectorT(maximumPointsList, optModel->maximumPointsList);
		copyVectorT(minimumPointsList, optModel->minimumPointsList);

		if (minimumPointsList.size() != maximumPointsList.size() + 1)
		{
			cout << "Warning! Inconsitent between minimum and maximum points" << endl;
			outputMinimumPoints("minimumPoints_search_debug", originalBindingEnergy, originalSpectrum, minimumPointsList);
			outputMaximumPoints("maximumPoints_search_debug", originalBindingEnergy, originalSpectrum, maximumPointsList);
			continue;
		}

		// Remove unnecessary minimum points
		removeMinimumPointsByHeightLimit2(originalBindingEnergy, maximumPointsList, minimumPointsList, smoothedSpectrum);

		cout << "m_SG, nMin, nMax= " << m_SG << ", " << minimumPointsList.size() << ", " << maximumPointsList.size() << endl;

		cout << "minimumPointList:" << endl;
		cout << "ID,position,height" << endl;
		for (int i = 0; i < minimumPointsList.size(); i++)
		{
			cout << i << "," << minimumPointsList[i]->position << "," << minimumPointsList[i]->height << endl;
			ofs2 << m_SG << "," << minimumPointsList[i]->position << "," << minimumPointsList[i]->height << endl;
		}

		// shinotsuka20181106
		if (minimumPointsList.size() < 2)
		{
			cout << "Error: minimumPointList does not have two or more elements. Skip this loop" << endl;
			continue;
		}

		// Get a minimum point at the most higher B.E. side,
		// and the weekest minimum point
		int idHigher, idLower;
		double bgHigher, bgLower;
		if (settings.backgroundType == "Shirley")
		{
			getBoundary(minimumPointsList, idHigher, idLower, bgHigher, bgLower);
		}
		else if (settings.backgroundType == "Linear")
		{
			idHigher = minimumPointsList.front()->id;
			idLower = minimumPointsList.back()->id;
			bgHigher = minimumPointsList.front()->height;
			bgLower = minimumPointsList.back()->height;
		}
		else
		{
			cout << "Error at mainFlow(). Unknown settings.backgroundType = " << settings.backgroundType << endl;
			exit(1);
		}
		cout << "Background boundary" << endl;
		cout << idHigher << " " << originalBindingEnergy[idHigher] << " " << bgHigher << endl;
		cout << idLower << " " << originalBindingEnergy[idLower] << " " << bgLower << endl;

		// shinotsuka20181106
		if (idHigher == idLower)
		{
			cout << "Error: minimumPointList does not have a range. Skip this loop" << endl;
			continue;
		}

		int ndata = idLower - idHigher + 1;

		cout << "originalNdata=" << originalNdata << endl
			 << "ndata=" << ndata << endl;

		vector<double> spectrumSlice = cutoffEnergyRange(originalSpectrum, idHigher, idLower);
		vector<double> bindingEnergySlice = cutoffEnergyRange(originalBindingEnergy, idHigher, idLower);
		vector<double> smoothedSpectrumSlice = cutoffEnergyRange(smoothedSpectrum, idHigher, idLower);
		vector<double> thirdDerivativeSlice = cutoffEnergyRange(thirdDerivative, idHigher, idLower);

		// Estimate an initial background spectrum by the iterative Shirley method
		vector<double> background(ndata), diffSpectrum(ndata);

		// Make a background[] by the iterative Shirley method.
		// Also make a diffSpectrum[] <-- spectrum[] - background[].
		if (!calcInitialBackground(ndata, bindingEnergySlice, spectrumSlice, smoothedSpectrumSlice, minimumPointsList, settings, background, diffSpectrum))
		{
			cout << "Error: Failed to get initial background. Skip this loop" << endl;
			continue;
		}

		// Extend the background[] and diffSpectrum to whole energy range
		vector<double> background_extend(originalBindingEnergy.size());
		vector<double> diffSpectrum_extend(originalBindingEnergy.size());
		for (int i = 0; i < originalBindingEnergy.size(); i++)
		{
			if (i <= idHigher)
			{
				background_extend[i] = background.front();
			}
			else if (i <= idLower)
			{
				background_extend[i] = background[i - idHigher];
			}
			else
			{
				background_extend[i] = background.back();
			}
			diffSpectrum_extend[i] = originalSpectrum[i] - background_extend[i];
		}

		bgHigher = background.front();
		bgLower = background.back();
		cout << "Debug. bgHigher = " << bgHigher << endl;
		cout << "Debug. bgLower = " << bgLower << endl;

		optModel->idHigher = idHigher;
		optModel->idLower = idLower;
		optModel->background = background_extend;
		optModel->diffSpectrum = diffSpectrum_extend;

		vector<PeakParameter *> peakList;
		if (detectPeaks(bindingEnergySlice, diffSpectrum, thirdDerivativeSlice, settings, peakList))
		{
			cout << "Warning at detectPeaks()." << endl;
		}

		if (peakList.size() == 0)
		{
			continue;
		}

		// if (1) {
		//	//Output all initial models
		//	ostringstream oss;
		//	oss << setw(4) << setfill('0') << m_SG;
		//	output_result("peak_initial_mSG" + oss.str(), bindingEnergySlice, spectrumSlice, background, peakList, originalBindingEnergy, originalSpectrum, idHigher, idLower, smoothedSpectrum, settings);
		// }

		optModel->peakList = peakList;

		// NRMSE (Normalized root mean square error) between original spectrum and smoothed spectrum
		optModel->nrmseSmoothedSpectrum = calcRMSE(optModel->spectrum, optModel->smoothedSpectrum, 1);

		optModelList.push_back(optModel);
	}

	ofs2.close();

	int numModel = optModelList.size();
	printf("Number of models = %d\n", numModel);
	if (numModel == 0)
	{
		printf("Error. No initial model was prepared.\n");
		exit(1);
	}

	int idHigherGlobal, idLowerGlobal;
	getGlobalEnergyRange(optModelList, idHigherGlobal, idLowerGlobal);
	cout << "idHigherGlobal = " << idHigherGlobal << endl;
	cout << "idLowerGlobal  = " << idLowerGlobal << endl;

	for (OptimalModel *om : optModelList)
	{
		om->idHigherGlobal = idHigherGlobal;
		om->idLowerGlobal = idLowerGlobal;
		om->originalBindingEnergy = &originalBindingEnergy;
		om->originalSpectrum = &originalSpectrum;
		om->settings = settings;
	}

	//// Remove optModels whose NRMSE are greater than a threshold
	// const double nrmseThreshold = 0.08;
	// pickupModel(optModelList, nrmseThreshold);
	// cout << "\nNumber of optModelList after pickup = " << optModelList.size() << endl;

	ofs2.open("out_rmspeSmoothedSpectrum.csv");
	ofs2 << "mSG,RMSPE" << endl;
	for (OptimalModel *om : optModelList)
	{
		ofs2 << om->mSG << "," << om->nrmseSmoothedSpectrum << endl;
	}
	ofs2.close();

	mytimer.show("before optimization");

	cout << endl;
	cout << "#=========================================================#" << endl;
	cout << "#                  Start optimization                     #" << endl;
	cout << "#=========================================================#" << endl;
	cout << endl;

	ofstream ofs;
	ofs.open("out_auto_smoothing2.csv");
	ofs << "m_SG,num_peaks,chi2_old,chi2_new,reduced_chi2_new,aic_poisson,bic_poisson,aic_gauss,bic_gauss,min_width,max_width,scaling_factor,cAIC_poisson,cAIC_gauss,smoothing_window,isConverged,numTotalOptimizeLoop,outerOptimizeLoop" << endl;

	// Slice between [idHigherGlobal, idLowerGlobal]
	// int ndata = idLowerGlobal - idHigherGlobal + 1;
	vector<double> bindingEnergySlice = cutoffEnergyRange(originalBindingEnergy, idHigherGlobal, idLowerGlobal);
	vector<double> spectrumSlice = cutoffEnergyRange(originalSpectrum, idHigherGlobal, idLowerGlobal);

	// Set global spect
	spect.setdata(bindingEnergySlice, spectrumSlice);

	// Start optimization using thread
	int numThreadLimit = get_numThreadLimit();
	omp_set_num_threads(numThreadLimit);
// printf("omp_get_num_threads() = %d\n", omp_get_num_threads( ));
// Parallelization by openMP
// 序盤にピーク本数の多いモデルが集まっているのでチャンクサイズ＝１がベスト
// #pragma omp parallel for schedule(static, 1)
#pragma omp parallel for
	for (int ll = 0; ll < numModel; ll++)
	{
		optimizeFlow_wrapper(optModelList[ll]);
	}

	mytimer.show("after optimization");
	cout << "mainFlow() finished" << endl;
}

void optimizeFlow_wrapper(OptimalModel *om)
{
	cout << endl
		 << "Optimization for m_SG = " << om->mSG << endl;

	// if (om->mSG == 2400) {
	//	cout << "Press any key";
	//	getchar();
	// }

	// Slice between [idHigherGlobal, idLowerGlobal]
	int ndata = om->idLowerGlobal - om->idHigherGlobal + 1;
	int m_SG = om->mSG;
	// vector<double> bindingEnergySlice = cutoffEnergyRange(originalBindingEnergy, idHigherGlobal, idLowerGlobal);
	// vector<double> spectrumSlice = cutoffEnergyRange(originalSpectrum, idHigherGlobal, idLowerGlobal);
	vector<double> backgroundSlice = cutoffEnergyRange(om->background, om->idHigherGlobal, om->idLowerGlobal);
	// vector<double> smoothedSpectrumSlice = cutoffEnergyRange(om->smoothedSpectrum, idHigherGlobal, idLowerGlobal);
	// vector<double> thirdDerivativeSlice = cutoffEnergyRange(om->thirdDerivative, idHigherGlobal, idLowerGlobal);
	// vector<double> diffSpectrumSlice = cutoffEnergyRange(om->diffSpectrum, idHigherGlobal, idLowerGlobal);

	// vector<PeakParameter *> peakList_Temp = om->peakList;
	vector<PeakParameter *> peakList_Temp;

	for (auto p : om->peakList)
	{
		peakList_Temp.push_back(new PeakParameter(*p));
	}

	vector<double> sumSpectrum;
	double chi2_new, aic_poisson_new, bic_poisson_new, reduced_chi2_new, scalingfactor, aic_gauss_new, bic_gauss_new, caic_poisson, caic_gauss;
	// int idHigher = om->idHigher;
	// int idLower = om->idLower;
	double bgHigher = om->background.front();
	double bgLower = om->background.back();

	// make sure to get positive chi-square, BG must be positive.
	bgHigher = max(bgHigher, 1e-15);
	bgLower = max(bgLower, 1e-15);

	bool isConverged = false;
	int numTotalOptimizeLoop = 0;
	int outerOptimizeLoop = 0;
	// Optimize paramters between [idHigherGlobal, idLowerGlobal]
	if (optimizeFlow(ndata, m_SG, spect.x, spect.y, backgroundSlice,
					 om->settings, peakList_Temp,
					 sumSpectrum, chi2_new, aic_poisson_new, bic_poisson_new, reduced_chi2_new,
					 bgHigher, bgLower, scalingfactor, aic_gauss_new, bic_gauss_new, caic_poisson, caic_gauss,
					 isConverged, numTotalOptimizeLoop, outerOptimizeLoop))
	{
		cout << "Warning at peakSeparationMain()." << endl;
		exit(1);
	}

	// Output all figures for m_SG
	if (1)
	{
		ostringstream oss;
		oss << setw(4) << setfill('0') << m_SG;
		output_result("peak_search_mSG" + oss.str(), spect.x, spect.y, backgroundSlice, peakList_Temp, *(om->originalBindingEnergy), *(om->originalSpectrum), om->idHigherGlobal, om->idLowerGlobal, om->smoothedSpectrum, om->settings);
	}

	return;
}

int get_numThreadLimit()
{
	// Start optimization using thread
	int NumOfCores = (int)thread::hardware_concurrency();
	printf("Total number of cores = %d\n", NumOfCores);
	int omp_num_threads = get_omp_num_threads();
	printf("OMP_NUM_THREADS = %d\n", omp_num_threads);

	int numThreadLimit = min(max(1, omp_num_threads), max(1, NumOfCores - 1));
	printf("Upper limit of number of cores = %d\n", numThreadLimit);
	return numThreadLimit;
}

int get_omp_num_threads()
{
	int num_threads = 0;
	const char *tmp = getenv("OMP_NUM_THREADS");
	string env_var(tmp ? tmp : "");
	if (env_var.empty())
	{
		cerr << "[ERROR] No such variable found!" << endl;
		// exit(EXIT_FAILURE);
	}
	else
	{
		num_threads = stod(env_var);
	}
	return num_threads;
}

/*
 * ---------------------------------------------------------------------
 * Make a list for round-robin smoothing.
 * Each element means a windows size of a Savitzky-Golay smoothing.
 * ---------------------------------------------------------------------
 */
vector<int> make_smoothingList()
{
	// Make a list of {5, 5, 5, ..., 7, 7, 7, ..., ..., 15, 15, 15}
	vector<int> smList;
	vector<vector<int>> smoothingWindowList = {{5, 480}, {7, 480}, {9, 480}, {11, 480}, {13, 480}, {15, 480}};
	for (int i = 0; i < smoothingWindowList.size(); i++)
	{
		smoothingWindowList[i][0], smoothingWindowList[i][1];
		for (int j = 0; j < smoothingWindowList[i][1]; j++)
		{
			smList.push_back(smoothingWindowList[i][0]);
		}
	}

	cout << "Number of smoothing = " << smList.size() << endl;
	for (auto &sg_window : smList)
	{
		cout << sg_window << " ";
	}
	cout << endl;

	return smList;
}

/*
 * ---------------------------------------------------------------------
 * Make a pick up list, which contains {0, 1}.
 * The size of the list is the same as the number of round-robin smoothing.
 * When "0", smoothing is done, but optimization is skipped.
 * When "1", both smoothing and optimization is done.
 * ---------------------------------------------------------------------
 * e.g. {0, 1, 0, 1, ..., 0, 0, 0, 1, 0, 0, 0, 1, ...}
 * ---------------------------------------------------------------------
 */
vector<int> make_pickupList(const int n)
{
	vector<int> pickupList;
	for (int i = 1; i <= n; i++)
	{
		int pickup = 1;
		if (i <= 40)
		{
			if (i % 2 != 0)
			{
				pickup = 0;
			}
		}
		else if (i <= 80)
		{
			if (i % 4 != 0)
			{
				pickup = 0;
			}
		}
		else if (i <= 240)
		{
			if (i % 8 != 0)
			{
				pickup = 0;
			}
		}
		else if (i <= 480)
		{
			if (i % 16 != 0)
			{
				pickup = 0;
			}
		}
		else if (i <= 960)
		{
			if (i % 16 != 0)
			{
				pickup = 0;
			}
		}
		else
		{
			if (i % 32 != 0)
			{
				pickup = 0;
			}
		}
		pickupList.push_back(pickup);
	}

	int numPickupList = (int)pickupList.size();
	cout << "Number of pickupList = " << numPickupList << endl;
	int numDoOptimize = accumulate(pickupList.begin(), pickupList.end(), 0);
	int numSkipOptimize = n - numDoOptimize;
	cout << "Number of DoOptimize = " << numDoOptimize << endl;
	cout << "Number of SkipOptimize = " << numSkipOptimize << endl;

	if (n != numPickupList)
	{
		cout << "Error: Number of pickupList is inconsitent with number of smoothing." << endl;
		exit(1);
	}

	return pickupList;
}

/*
 * ---------------------------------------------------------------------
 * Cut off energy range
 * INPUT
 *    x: vector data
 *    idS: ID of a starting point
 *    idE: ID of a end point
 * RETURN
 *    A sliced data
 * ---------------------------------------------------------------------
 */
vector<double> cutoffEnergyRange(const vector<double> &x, const int idS, const int idE)
{
	return vector<double>(x.begin() + idS, x.begin() + idE + 1);
}

/*
 * Get two IDs of global energy range which will be used at the parameter optimization
 */
void getGlobalEnergyRange(const vector<OptimalModel *> &optModelList, int &idHigherGlobal, int &idLowerGlobal)
{
	idHigherGlobal = optModelList[0]->idHigher;
	idLowerGlobal = optModelList[0]->idLower;
	// cout << "--------------------------------------------" << endl;
	for (auto om : optModelList)
	{
		// cout << om->idHigher << " " << om->idLower << endl;
		idHigherGlobal = min(idHigherGlobal, om->idHigher);
		idLowerGlobal = max(idLowerGlobal, om->idLower);
	}
	// cout << "--------------------------------------------" << endl;
	// cout << idHigherGlobal << " " << idLowerGlobal << endl;
	// cout << "--------------------------------------------" << endl;
}

// Get a minimum and maximum peak width from a peakList_Temp
void summarize_peakList(const vector<PeakParameter *> &peakList_Temp,
						double &min_width, double &max_width)
{
	min_width = 0.0;
	max_width = 0.0;

	if (peakList_Temp.size() == 0)
	{
		return;
	}

	min_width = peakList_Temp[0]->width;
	max_width = min_width;

	for (int j = 0; j < peakList_Temp.size(); j++)
	{
		min_width = min(min_width, peakList_Temp[j]->width);
		max_width = max(max_width, peakList_Temp[j]->width);
	}
}

void getBoundary(vector<MinimumPoints *> &minimumPointsList, int &idHigher, int &idLower,
				 double &bgHigher, double &bgLower)
{
	if (minimumPointsList.size() < 2)
	{
		cout << "Error: minimumPointList does not have two or more elements." << endl;
		exit(1);
	}

	double mostMinimumPoint_height = minimumPointsList.front()->height;
	idHigher = minimumPointsList.front()->id;
	bgHigher = minimumPointsList.front()->height;
	idLower = idHigher;
	for (int i = 1; i < minimumPointsList.size(); i++)
	{
		if (minimumPointsList[i]->height < mostMinimumPoint_height)
		{
			mostMinimumPoint_height = minimumPointsList[i]->height;
			idLower = minimumPointsList[i]->id;
			bgLower = minimumPointsList[i]->height;
		}
	}
	// if (idHigher == idLower) {
	//	cout << "Warning: minimumPointList does not have a range." << endl;
	// }
}

int hasNegative(const vector<double> &v)
{
	for (int i = 0; i < v.size(); i++)
	{
		if (v[i] < 0.0)
		{
			return 1;
		}
	}
	return 0;
}

void pickupModel(vector<OptimalModel *> &optModelList, const double nrmseThreshold)
{
	// Remove optModels whose NRMSE are greater than a threshold
	int cnt = 0;
	const int numModel = (int)optModelList.size();
	for (int i = numModel - 1; i >= 0; i--)
	{
		if (optModelList[i]->nrmseSmoothedSpectrum > nrmseThreshold)
		{
			optModelList.erase(optModelList.begin() + i);
			cnt++;
		}
	}
	cout << cnt << " models are removed because of a large NRMSE values" << endl;
	return;
}