#include "calcinitialbackground.h"

/*
* ---------------------------------------------------------------------
* Calculate initial background
* INPUT
*    ndata:              Number of data points
*    spectrum:           Raw spectrum
*    settings:           Settings
* OUTPUT
*    backgroundList:     An iteratively determined Shirley background
*    diffSpectrum:       The difference spectrum <-- (raw spectrum) - (background)
* ---------------------------------------------------------------------
*/
int calcInitialBackground(const int ndata, const vector<double> &bindingEnergy, 
	const vector<double> &spectrum, const vector<double> &smoothedSpectrum, 
	const vector<MinimumPoints*> &minimumPointsList, const Settings &settings, 
	vector<double> &background,	vector<double> &diffSpectrum)
{
	const double initBGLower = smoothedSpectrum.back();
	const double initBGHigher = smoothedSpectrum.front();

	if (settings.backgroundType == "Shirley") {
		if (!calcInitialBackgroundShirley(ndata, bindingEnergy, spectrum, minimumPointsList, settings, background, diffSpectrum, initBGLower)) {
			cout << "Warning! calcInitialBackgroundShirley() not converged." << endl;
			return false;
		}
	}
	else if (settings.backgroundType == "Linear") {
		// Linear background using two points;
		// (bindingEnergy.front(), smoothedSpectrum.front()) and (bindingEnergy.back(), smoothedSpectrum.back())
		if (!calcInitialBackgroundLinear(ndata, spectrum, smoothedSpectrum, settings, background, diffSpectrum)) {
			cout << "Warning! calcInitialBackgroundLinear() failed." << endl;
			return false;
		}

		//if (!calcInitialBackgroundLinear2(ndata, bindingEnergy, spectrum, smoothedSpectrum, minimumPointsList, settings, background, diffSpectrum)) {
		//	cout << "Warning! calcInitialBackgroundLinear2() failed." << endl;
		//	return false;
		//}

		//ofstream ofs;
		//ofs.open("linearBG.csv");
		//ofs << "bindingEnergy,spectrum,smoothedSpectrum,background,diffSpectrum" << endl;
		//for (int i = 0; i < ndata; i++) {
		//	ofs << bindingEnergy[i] << "," << spectrum[i] << "," << smoothedSpectrum[i] << "," << background[i] << "," << diffSpectrum[i] << endl;
		//}
		//ofs.close();

		//exit(1);

	}
	else {
		cout << "Error at calcInitialBackground(). Unknown settings.backgroundType = " << settings.backgroundType << endl;
		exit(1);
	}

	return true;
}


/*
* ---------------------------------------------------------------------
* The iterative shirley method.
* INPUT
*    ndata:              Number of data points
*    spectrum:           Raw spectrum
*    settings:           Settings
* OUTPUT
*    backgroundList:     An iteratively determined Shirley background
*    diffSpectrum:       The difference spectrum <-- (raw spectrum) - (background)
* ---------------------------------------------------------------------
*/
int calcInitialBackgroundShirley(const int ndata, const vector<double> &bindingEnergy, 
	const vector<double> &spectrum, const vector<MinimumPoints*> &minimumPointsList,
	const Settings &settings, vector<double> &background,
	vector<double> &diffSpectrum, const double initBGLower)
{
	cout << "Iterative Shirley method" << endl;
	cout << "count,conversion,sumNew" << endl;

	const double deltaConversion = 0.002;
	bool isConverged = false;

	double intensityLeft;
	double intensityRight;
	
	int ndataNarrow;	// number of points to narrow
	int ndataExtend;	// Number of points to extend
	int nloop = 0;
	for (int i = 0; i < minimumPointsList.size(); i++) {
		if (minimumPointsList[i]->position > bindingEnergy.back()) {
			nloop++;
		}
	}
	cout << endl << " nloop = " << nloop << endl;
	vector<double> pnRatio(nloop);
	vector<double> backgroundArea(nloop);
	vector<int> goodPointsList;
	goodPointsList.clear();

	int numConvergedBG = 0;

	vector<vector<double>> backgroundList(minimumPointsList.size(), vector<double>(ndata));

	for (int k = 0; k < nloop; k++) {
		bool isSkipCalculation = false;
			intensityLeft = minimumPointsList[k]->height;
			//intensityRight = spectrum.back();
			//intensityRight = minimumPointsList[idLower]->height;
			intensityRight = initBGLower;
			double minimumPointsXaxis = (double)minimumPointsList[k]->position;
			ndataNarrow = 0;
			for (int i = 0; i < ndata; i++) {
				if (minimumPointsXaxis > bindingEnergy[i]) {
					ndataNarrow++;
				}
			}
			ndataExtend = ndata - ndataNarrow;

		// Initialize diffSpectrum[]
		for (int i = 0; i < ndata; i++) {
			diffSpectrum[i] = spectrum[i];
		}
		double sumOld = 1e300; // The integral of diffSpectrum[], an old value.
		double sumNew = 1e300; // The integral of diffSpectrum[],  a new value.
		vector<double> areaP(ndataNarrow); // The cummulative integration of diffSpectrum[].

		for (int count = 0; count < settings.maxLoopIterativeShirley * 2; count++) {
			// Calculate the areas P_i step by step from the left side.
			// The steps of abscissa are all the same, assumed to be 1.0. No meaning.

			// Calculate the background of a narrow range
			areaP[0] = 0.0;
			for (int j = 1; j < ndataNarrow; j++) {
				areaP[j] = areaP[j - 1] + (diffSpectrum[j - 1 + ndataExtend] + diffSpectrum[j + ndataExtend]) * 0.5;
			}
			sumNew = areaP[ndataNarrow - 1];

			for (int j = 0; j < ndataNarrow; j++) {
				backgroundList[k][j + ndataExtend] = (intensityLeft - intensityRight) * (sumNew - areaP[j]) / sumNew + intensityRight;
				diffSpectrum[j + ndataExtend] = spectrum[j + ndataExtend] - backgroundList[k][j + ndataExtend];
			}

			if (!settings.allowNegativeIterativeShirley) {
				// Negative values of diffSpectrum[] are not allowed.
				for (int j = 0; j < ndataNarrow; j++) {
					diffSpectrum[j + ndataExtend] = max(diffSpectrum[j + ndataExtend], 0.0);
				}
			}

			double conversion = abs((sumNew - sumOld) / sumOld);
			cout << count << "," << conversion << "," << sumNew << endl;

			isConverged = (conversion < deltaConversion);
			if (isConverged) break;
			sumOld = sumNew;
			/**************************************************************/
			if (count >= settings.maxLoopIterativeShirley && sumNew < 0) {
				//cout << "Remove a minimumPointList[" << k << "] because it does not converge and the area is a negative value"
				//	<< minimumPointsList[k]->height << " " << minimumPointsList[k]->position << endl;
				//minimumPointsList.erase(minimumPointsList.begin() + k);
				isSkipCalculation = true;
				break;
			}
			/***************************************************************/
		}

		//if (isSkipCalculation) {
		//	// Initialize diffSpectrum[]
		//	for (int i = 0; i < ndata; i++) {
		//		diffSpectrum[i] = spectrum[i];
		//	}
		//	sumOld = 1e300; // The integral of diffSpectrum[], an old value.
		//	sumNew = 1e300; // The integral of diffSpectrum[],  a new value.
		//	for (int count = 0; count < settings.maxLoopIterativeShirley; count++) {
		//		areaP[0] = 0.0;
		//		for (int j = 1; j < ndataNarrow; j++) {
		//			areaP[j] = areaP[j - 1] + (diffSpectrum[j - 1 + ndataExtend] + diffSpectrum[j + ndataExtend]) * 0.5;
		//		}
		//		sumNew = areaP[ndataNarrow - 1];

		//		for (int j = 0; j < ndataNarrow; j++) {
		//			backgroundList[k][j + ndataExtend] = (intensityLeft - intensityRight) * (sumNew - areaP[j]) / sumNew + intensityRight;
		//			diffSpectrum[j + ndataExtend] = spectrum[j + ndataExtend] - backgroundList[k][j + ndataExtend];
		//		}

		//		// Negative values of diffSpectrum[] are not allowed.
		//		for (int j = 0; j < ndataNarrow; j++) {
		//			diffSpectrum[j + ndataExtend] = max(diffSpectrum[j + ndataExtend], 0.0);
		//		}

		//		double conversion = abs((sumNew - sumOld) / sumOld);
		//		cout << count << "," << conversion << "," << sumNew << endl;

		//		isConverged = (conversion < deltaConversion);
		//		if (isConverged) break;
		//		sumOld = sumNew;
		//	}
		//}

		if (isSkipCalculation) {
			// Give up to get Shirley background
			cout << "Error: Iterative Shirley method does not converged." << endl;
			//return false;
			continue;
		}

		numConvergedBG += 1;

		// Extend background
		// Calculate coefficients
		double c = (intensityLeft - intensityRight) / sumNew;

		vector<double> areaR(ndataExtend + 1); // The cummulative integration of diffSpectrum[].
		areaR[0] = 0;
		for (int n = 1; n <= ndataExtend; n++) {
			areaR[n] = areaR[n - 1] + (diffSpectrum[ndataExtend - n + 1] + diffSpectrum[ndataExtend - n + 2]) * 0.5;
			backgroundList[k][ndataExtend - n] = areaR[n] * c + intensityLeft;
			diffSpectrum[ndataExtend - n] = spectrum[ndataExtend - n] - backgroundList[k][ndataExtend - n];
		}

		// Calculate positive and negative areas and their ratios
		double positiveAreaSum = 0.0;
		double negativeAreaSum = 0.0;
		for (int i = 0; i < ndata; i++) {
			double spectrum_from_BG = spectrum[i] - backgroundList[k][i];
			if (spectrum_from_BG > 0) {
				positiveAreaSum += abs(spectrum_from_BG);
			}
			else {
				negativeAreaSum += abs(spectrum_from_BG);
			}
		}
		pnRatio[k] = negativeAreaSum / positiveAreaSum;
		cout << "k=" << k << "," << "pnRatio=" << pnRatio[k] << "," << "positiveArea=" << positiveAreaSum << "," << "negativeArea=" << negativeAreaSum << endl;

		for (int i = 0; i < ndata; i++) {
			backgroundArea[k] += backgroundList[k][i];
		}
		cout << "k=" << k << "," << "backgroundArea=" << backgroundArea[k] << endl;

		//shinotsuka20180601
		//if (pnRatio[k] < 0.01) {
		//	goodPointsList.push_back(k);
		//}
		if (pnRatio[k] < settings.allowedNegativeRatio) {
			goodPointsList.push_back(k);
		}
	}

	if (numConvergedBG == 0) {
		return false;
	}

	// Determine the optimal ID of the minimum point by the following condition.
	// If there are "good" points whose negative areas are small enough,
	// the optimal ID is determined so that the background area becomes maximum.
	// If there is NO "good" points, the optimal ID is determined so that
	// the negative area becomes minimum.
	int optimalID = 0;
	double maxBackgroundArea = 0;
	double minPnRatio = 1;

	if (goodPointsList.size() > 0) {
		for (int i = 0; i < goodPointsList.size(); i++) {
			if (maxBackgroundArea < backgroundArea[goodPointsList[i]]) {
				maxBackgroundArea = backgroundArea[goodPointsList[i]];
				optimalID = goodPointsList[i];
			}
		}
	}

	else {
		for (int k = 0; k < nloop; k++) {
			if (minPnRatio > pnRatio[k]) {
				minPnRatio = pnRatio[k];
				optimalID = k;
			}
		}
	}

	// Get the optimal background result
	for (int i = 0; i < ndata; i++) {
		background[i] = backgroundList[optimalID][i];
		diffSpectrum[i] = spectrum[i] - backgroundList[optimalID][i];
	}
	cout << " ***** optimalID = " << optimalID << " ******" << endl;

	//return isConverged;
	return true;
}

/*
* ---------------------------------------------------------------------
* Linear background
* INPUT
*    ndata:              Number of data points
*    spectrum:           Raw spectrum
*    smoothedSpectrum:   Smoothed spectrum
* OUTPUT
*    background:         Linear background
*    diffSpectrum:       Difference spectrum <-- (raw spectrum) - (background)
* ---------------------------------------------------------------------
*/
int calcInitialBackgroundLinear(const int ndata, 
	const vector<double> &spectrum, const vector<double> &smoothedSpectrum, 
	const Settings &settings, vector<double> &background,	vector<double> &diffSpectrum)
{
	make_background_Linear(smoothedSpectrum.front(), smoothedSpectrum.back(), background);
	for (int i = 0; i < ndata; i++) {
		diffSpectrum[i] = spectrum[i] - background[i];
	}

	// When the ratio (negative area / positive area) < settings.allowedNegativeRatio,
	// the background is appropriate.
	double positiveAreaSum = 0.0;
	double negativeAreaSum = 0.0;
	for (int i = 0; i < ndata; i++) {
		double spectrum_from_BG = diffSpectrum[i];
		if (spectrum_from_BG > 0) {
			positiveAreaSum += abs(spectrum_from_BG);
		}
		else {
			negativeAreaSum += abs(spectrum_from_BG);
		}
	}
	double pnRatio = negativeAreaSum / positiveAreaSum;
	cout << "pnRatio=" << pnRatio << "," << "positiveArea=" << positiveAreaSum << "," << "negativeArea=" << negativeAreaSum << endl;

	return (pnRatio < settings.allowedNegativeRatio);
}

int calcInitialBackgroundLinear2(const int ndata, const vector<double> &bindingEnergy,
	const vector<double> &spectrum, const vector<double> &smoothedSpectrum, 
	const vector<MinimumPoints*> &minimumPointsList, const Settings &settings, 
	vector<double> &background, vector<double> &diffSpectrum)
{
	vector<double> bgtemp(ndata);
	double maxBackgroundArea = 0;

	for (int k = 0; k < minimumPointsList.size() - 1; k++) {
		double x1 = minimumPointsList[k]->position;
		double y1 = minimumPointsList[k]->height;
		for (int j = k + 1; j < minimumPointsList.size(); j++) {
			double x2 = minimumPointsList[j]->position;
			double y2 = minimumPointsList[j]->height;
			double a = (y2 - y1)/(x2- x1);
			for (int i = 0; i < ndata; i++) {
				bgtemp[i] = a * (bindingEnergy[i] - x1) + y1;
				diffSpectrum[i] = spectrum[i] - bgtemp[i];
			}
			// When the ratio (negative area / positive area) < settings.allowedNegativeRatio,
			// the background is appropriate.
			double positiveAreaSum = 0.0;
			double negativeAreaSum = 0.0;
			for (int i = 0; i < ndata; i++) {
				double spectrum_from_BG = diffSpectrum[i];
				if (spectrum_from_BG > 0) {
					positiveAreaSum += abs(spectrum_from_BG);
				}
				else {
					negativeAreaSum += abs(spectrum_from_BG);
				}
			}
			double pnRatio = negativeAreaSum / positiveAreaSum;
			cout << "pnRatio=" << pnRatio << "," << "positiveArea=" << positiveAreaSum << "," << "negativeArea=" << negativeAreaSum << endl;

			double backgroundArea = (bgtemp.front() + bgtemp.back()) * 0.5;

			if (pnRatio < settings.allowedNegativeRatio && backgroundArea > maxBackgroundArea) {
				maxBackgroundArea = backgroundArea;
				background = bgtemp;
			}

		}



	}

	return (maxBackgroundArea > 0.0);
}
