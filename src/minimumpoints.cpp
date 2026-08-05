#include "minimumpoints.h"

/*
* ---------------------------------------------------------------------
* Detect minimum points from a raw spectrum, in order to achieve the "wise" optimization.
* INPUT
*    ndata:              Number of data points
*    bindingEnergy:      Binding energies (eV)
*    spectrum:           Raw spectrum
*    settings:           Settings
* OUTPUT
*    minimumPointsList:  List of detected minimum points
* ---------------------------------------------------------------------
*/
//int minimumPointsMain(const int ndata, const vector<double> &bindingEnergy, const vector<double> &spectrum, vector<MinimumPoints*> &minimumPointsList, const Settings &settings)
//{
//	cout << "minimumPointsMain() starts." << endl;
//	vector<double> smoothedSpectrum;
//	smoothedSpectrum.resize(ndata);
//	vector<double> firstDerivative;
//	firstDerivative.resize(ndata);
//	vector<double> secondDerivative;
//	secondDerivative.resize(ndata);
//	vector<MaximumPoints*> maximumPointsList;
//
//	vector<double> sumSpectrum;
//	vector<double> smoothedDiffSpectrum(ndata);
//	int m_SG = 0;
//	int nMin = max((int)(bindingEnergy.front() - bindingEnergy.back()) / 10 + 1, 2);
//
//	ofstream ofs1;
//
//	for (m_SG = 5; m_SG < ndata / 3; m_SG += 2) {
//	//	if (minimumPoints(ndata, m_SG, bindingEnergy, spectrum, settings, smoothedSpectrum, firstDerivative, minimumPointsList, sumSpectrum)) {
//	//		cout << "Error at detectMinimumPoints()." << endl;
//	//	}
//	//	nMin[m_SG] = (int)minimumPointsList.size();
//	//for (int i = 0; i < ndata; i++) {
//		//	for (int j = 0; j < minimumPointsList.size(); j++) {
//		//		if (i == minimumPointsList[j]->id) {
//		//			minimumPointsHistogram[i] += m_SG;
//		//		}
//		//	}
//		//}
//		if (minimumPoints(ndata, m_SG, bindingEnergy, spectrum, settings, smoothedSpectrum, firstDerivative, secondDerivative, maximumPointsList, minimumPointsList, sumSpectrum)) {
//			cout << "Error at detectMinimumPoints()." << endl;
//		}
//
//		//int numRemovedMinimumPointsByHeightLimit
//		//shinotsuka20171205
//		//removeMinimumPointsByHeightLimit(bindingEnergy, maximumPointsList, minimumPointsList, smoothedSpectrum);
//		removeMinimumPointsByHeightLimit2(bindingEnergy, maximumPointsList, minimumPointsList, smoothedSpectrum);
//		//cout << "The number of removed minimum points by height limit = " << numRemovedMinimumPointsByHeightLimit << endl;
//
//		// numRemovedMinimumPointsByEnergyIntervalLimit
//		removeMinimumPointsByEnergyIntervalLimit(bindingEnergy, minimumPointsList, smoothedSpectrum);
//		//cout << "The number of removed minimum points by binding energy interval limit = " << numRemovedMinimumPointsByIntervalLimit << endl;
//
//		cout << "m_SG, minimumPointsList.size(), nMin= " << m_SG << "," << minimumPointsList.size() << "," << nMin << endl;
//
//		if (minimumPointsList.size() <= nMin) {
//		//if (minimumPointsList.size() < 6) {
//			break;
//		}
//		//ofs4 << m_SG;
//			//for (int i = 0; i < minimumPointsList.size(); i++) {
//			//	ofs4 << "," << bindingEnergy[minimumPointsList[i]->id];
//			//}
//			//ofs4 << endl;
//	}
//	//ofs4.close();
//	cout << "Smoothing points when finding minimum points = " << m_SG << endl;
//	//if (minimumPoints(ndata, m_SG, bindingEnergy, spectrum, settings, smoothedSpectrum, firstDerivative, secondDerivative, minimumPointsList, sumSpectrum)) {
//	//	cout << "Error at detectMinimumPoints()." << endl;
//	//}
//	ofs1.open("out_Savitzky_Golay_first.csv");
//	ofs1 << "Binding_Energy(eV),spectrum,firstDerivative" << endl;
//	for (int i = 0; i < ndata; i++) {
//		ofs1 << bindingEnergy[i] << "," << spectrum[i] << "," << firstDerivative[i] << endl;
//	}
//	//ofs2.open("out_minimumPointsHistogram.csv");
//	//ofs2 << "id,spectrum,minimumPointsHistogram" << endl;
//	//for (int i = 0; i < ndata; i++) {
//	//	ofs2 << i << "," << spectrum[i] << "," << minimumPointsHistogram[i] << endl;
//	//}
//	//ofs3.open("m_SG__nmdata.csv");
//	//ofs3 << "m_SG,minimumPointsList.size()" << endl;
//	//for (m_SG = 5; m_SG < 92; m_SG += 2) {
//	//	ofs3 << m_SG << "," << nMin[m_SG] << endl;
//	//}
//	
//	//smoothedDiffSpectrum
//	//ofs5.open("out_smoothedDiffSpectrum.csv");
//	//ofs5 << "Binding_Energy(eV),spectrum,smoothedSpectrum,diff" << endl;
//	//for (int i = 0; i < ndata; i++) {
//	//	smoothedDiffSpectrum[i] = spectrum[i] - smoothedSpectrum[i];
//	//	ofs5 << bindingEnergy[i] << "," << spectrum[i] << "," << smoothedSpectrum[i] << "," << smoothedDiffSpectrum[i] << endl;
//	//}
//	outputMaximumPoints("maximumPoints_search", bindingEnergy, spectrum, maximumPointsList);
//
//	return 0;
//}

/*
* ---------------------------------------------------------------------
* Detect minimum points from a spectrum by the first derivative
* INPUT
*    ndata:              Number of data points
*    spectrum:           Raw spectrum
*    numSmooth:          Smoothing data points around a target point.
*                        It must be 3, 5, 7, 9, ... .
* OUTPUT
*    smoothedSpectrum:   The smoothed spectrum.
*    firstDerivative:    The first derivative of the smoothed spectrum.
* ---------------------------------------------------------------------
*/
int minimumPoints(const int ndata, const int m_SG, const vector<double> &bindingEnergy, const vector<double> &spectrum, const Settings &settings, vector<double> &smoothedSpectrum, vector<double> &firstDerivative, vector<double> &secondDerivative, vector<MaximumPoints*> &maximumPointsList, vector<MinimumPoints*> &minimumPointsList, vector<double> &sumSpectrum)
{
	if (smoothedSpectrum.size() != ndata) { smoothedSpectrum.resize(ndata); }
	if (firstDerivative.size() != ndata) { firstDerivative.resize(ndata); }
	if (secondDerivative.size() != ndata) { secondDerivative.resize(ndata); }

	//cout << "minimumPointsMain() starts." << endl;
	//if (savitzkyGolayFirst(ndata, spectrum, m_SG, smoothedSpectrum, firstDerivative)) {
	//	cout << "Error at savitzkyGolayFirst()." << endl;
	//	exit(1);
	//}
	if (savitzkyGolayGeneral(spectrum, m_SG, 0, smoothedSpectrum)
		|| savitzkyGolayGeneral(spectrum, m_SG, 1, firstDerivative)) {
		cout << "Error after savitzkyGolayGeneral() at minimumPoints()." << endl;
		exit(1);
	}
	//if (savitzkyGolayGeneral(spectrum, m_SG + 10, 2, secondDerivative)) {
	//	cout << "Error after savitzkyGolayGeneral() at minimumPoints()." << endl;
	//	exit(1);
	//}


	if(detectMinimumPoints(ndata, bindingEnergy, smoothedSpectrum, firstDerivative, settings, maximumPointsList, minimumPointsList)){
		cout << "Error at detectMinimumPoints()." << endl;
	}
	return 0;
}

/*
* ---------------------------------------------------------------------
* The Savitzky-Golay method to make a smoothed spectrum and its first derivative.
* INPUT
*    ndata:              Number of data points
*    spectrum:           Raw spectrum
*    numSmooth:          Smoothing data points around a target point.
*                        It must be 3, 5, 7, 9, ... .
* OUTPUT
*    smoothedSpectrum:   The smoothed spectrum.
*    firstDerivative:    The first derivative of the smoothed spectrum.
* ---------------------------------------------------------------------
*/
int savitzkyGolayFirst(const int ndata, const vector<double> &spectrum, const int numSmooth, vector<double> &smoothedSpectrum, vector<double> &firstDerivative)
{
	if (numSmooth % 2 == 0 || numSmooth < 3) {
		cout << "Smoothing data points must be positive integers 3, 5, 7, ...." << endl;
		return 1;
	}

	const int m = (numSmooth - 1) / 2;
	//cout << "savitzkyGolay(), m=" << m << endl;

	// Smoothing
	double normFactor = (4.0 * m * m - 1.0) * (2.0 * m + 3.0) / 3.0;
	double wtFactor = 3.0 * m * (m + 1.0);
	vector<double> wt0(m+1); // Weight factor for smoothing
	vector<double> wt1(m+1); // Weight factor for smoothed 1st derivative
	for (int j = 0; j <= m; j++) {
		double temp = (wtFactor - 1.0 - 5.0*j*j);
		wt0[j] = temp / normFactor;
		wt1[j] = (double)j;
		//cout << "wt0[" << j << "]=" << wt0[j] << "  wt3[" << j << "]="  << endl;
	}
	int j, i_plus_j, i_minus_j;
	double sum0, sum1;
//#pragma omp parallel for private(j,i_plus_j, i_minus_j,sum0,sum1,sum2)
	for (int i = 0; i < ndata; i++) {
		sum0 = wt0[0] * spectrum[i];
		sum1 = 0.0;

//#pragma omp parallel for reduction (+:sum1,sum2)
		for (j = 1; j <= m; j++) {
			i_plus_j = min(ndata - 1, i + j);
			i_minus_j = max(0, i - j);
			sum0 += wt0[j] * (spectrum[i_plus_j] + spectrum[i_minus_j]);
			sum1 += wt1[j] * (spectrum[i_plus_j] - spectrum[i_minus_j]);
		}
		smoothedSpectrum[i] = sum0;
		firstDerivative[i] = sum1;
	}

	return 0;
}

/*
* ---------------------------------------------------------------------
* Peak detection from a smoothed spectrum and its third derivative.
* INPUT
*    ndata:              Number of data points
*    xx:                 x-axis
*    smoothedSpectrum:   The smoothed spectrum.
*    thirdDerivative:    The third derivative of the smoothed spectrum.
*    settings:           Settings
* OUTPUT
*    peakList:           The peak parameters
* ---------------------------------------------------------------------
*/
int detectMinimumPoints(const int ndata, const vector<double> &xx, const vector<double> smoothedSpectrum, vector<double> &firstDerivative, const Settings &settings, vector<MaximumPoints*> &maximumPointsList, vector<MinimumPoints*> &minimumPointsList)
{
	// Adjust the smoothed first derivative spectrum.
	adjustFirstDeriv(firstDerivative);
	// Find minimum points of the smoothed spectrum.
	pickupMinimumPointsList(xx, smoothedSpectrum, firstDerivative, minimumPointsList);
	// Find maximum points of the smoothed spectrum.
	pickupMaximumPointsList(xx, smoothedSpectrum, firstDerivative, maximumPointsList, minimumPointsList);

	//shinotsuka20180115, rare case
	// ... -> min[99] -> max[99] -> min[100] -> min[101], then remove min[101]
	if (minimumPointsList.size() > 2 &&
		minimumPointsList[minimumPointsList.size()-2]->id > maximumPointsList[maximumPointsList.size()-1]->id) {
		// Remove minimumPointsList.end()
		minimumPointsList.pop_back();
	}

	return 0;
}

/*
* ---------------------------------------------------------------------
* Adjust the smoothed first derivative spectrum.
* When a value is zero, then modify it to a "very small positive value".
* INPUT/OUTPUT
*    firstDerivative:    The first derivative of the smoothed spectrum.
* ---------------------------------------------------------------------
*/
void adjustFirstDeriv(vector<double> &firstDerivative)
{
	const int ndata = (int)firstDerivative.size();

	bool hasZero_AbsFirstDeriv = false;
	for (int i = 0; i < ndata; i++) {
		if (firstDerivative[i] == 0.0) {
			hasZero_AbsFirstDeriv = true;
			break;
		}
	}

	if (hasZero_AbsFirstDeriv) {
		double minAbsFirstDeriv = 1e300;
		for (int i = 0; i < ndata; i++) {
			if (firstDerivative[i] != 0.0) {
				minAbsFirstDeriv = min(minAbsFirstDeriv, abs(firstDerivative[i]));
			}
		}

		for (int i = 0; i < ndata; i++) {
			if (firstDerivative[i] == 0.0) {
				//cout << "First derivative is zero at " << i << ". Change to " << abs(minAbsFirstDeriv * 1e-3 << endl;
				firstDerivative[i] = abs(minAbsFirstDeriv * 1e-3);
			}
		}
	}
}

/*
* ---------------------------------------------------------------------
* Pick up local minimum points.
* In case of a local minimum, the first derivative crosses zero
* from negative to positive. 
*
* INPUT
*    xx:                 x-axis
*    smoothedSpectrum:   The smoothed spectrum.
*    firstDerivative:    The first derivative of the smoothed spectrum.
*    settings:           Settings
*
* OUTPUT
*    minimumPointsList:  The local minimum points list.
* ---------------------------------------------------------------------
*/
void pickupMinimumPointsList(const vector<double> &xx, const vector<double> smoothedSpectrum, const vector<double> &firstDerivative, vector<MinimumPoints*> &minimumPointsList)
{
	//shinotsuka20190807
	//minimumPointsList.clear();
	releaseMemoryT(minimumPointsList);
	const int ndata = (int)firstDerivative.size();

	//add pseudo minimum point left
	if (firstDerivative[0] > 0) {
		double height = smoothedSpectrum[0];
		double position = xx[0];
		minimumPointsList.push_back(new MinimumPoints(height, position, 0));
	}

	for (int i = 1; i < ndata - 1; i++) {
		if (firstDerivative[i - 1] < 0 && firstDerivative[i] > 0) {
			int j = 0;
			if (abs(smoothedSpectrum[i - 1]) > abs(smoothedSpectrum[i])) {
				j = i;
			}
			else {
				j = i - 1;
			}
			double height = smoothedSpectrum[j];
			double position = xx[j];
			minimumPointsList.push_back(new MinimumPoints(height, position, j));
		}
	}

	//add pseudo minimum point right
	if (firstDerivative[ndata - 1] < 0) {
		double height = smoothedSpectrum[ndata - 1];
		double position = xx[ndata - 1];
		minimumPointsList.push_back(new MinimumPoints(height, position, ndata - 1));
	}
}
void pickupMaximumPointsList(const vector<double> &xx, const vector<double> smoothedSpectrum, const vector<double> &firstDerivative, vector<MaximumPoints*> &maximumPointsList, const vector<MinimumPoints*> &minimumPointsList)
{
	//shinotsuka20190807
	//maximumPointsList.clear();
	releaseMemoryT(maximumPointsList);
	if(minimumPointsList.size() < 2){
		return;
	}
	const int ndata = (int)firstDerivative.size();

	for (int i = minimumPointsList.front()->id + 1; i < ndata - 1; i++) {
		if (firstDerivative[i] < 0 && firstDerivative[i - 1] > 0) {
			int j = 0;
			if (abs(smoothedSpectrum[i - 1]) < abs(smoothedSpectrum[i])) {
				j = i;
			}
			else {
				j = i - 1;
			}
			double height = smoothedSpectrum[j];
			double position = xx[j];
			maximumPointsList.push_back(new MaximumPoints(height, position, j));
		}
	}
}


int removeMinimumPointsByEnergyIntervalLimit(const vector<double> &bindingEnergy, vector<MinimumPoints*> &minimumPointsList, const vector<double> &smoothedSpectrum)
{
	//shinotsuka20181106, Br3d of KBr�ɑΉ������邽�ߐ��l��5.0����4.0�ɏC������
	//const double bindingEnergyThreshold = 5.0;
	const double bindingEnergyThreshold = 4.0;
	int j;
	int numRemovedMinimumPoints = 0;
	const int ndata = (int)bindingEnergy.size();
	for (int i = (int)minimumPointsList.size() - 2; i >= 0; i--) {

		if (minimumPointsList[i+1]->id == ndata - 1 || minimumPointsList[i]->id == 0) {
			continue;
		}
		else {
			if (abs(minimumPointsList[i + 1]->position - minimumPointsList[i]->position) < bindingEnergyThreshold) {
				cout << minimumPointsList[i + 1]->position << " " << minimumPointsList[i]->position << endl;
				j = 0;
				//if (secondDerivative[minimumPointsList[i + 1]->id] < secondDerivative[minimumPointsList[i]->id]) {
				//if (abs(secondDerivative[minimumPointsList[i + 1]->id]) < abs(secondDerivative[minimumPointsList[i]->id])) {
				if (smoothedSpectrum[minimumPointsList[i + 1]->id] > smoothedSpectrum[minimumPointsList[i]->id]) {
					j = i + 1;
				}
				else {
					j = i;
				}
				cout << "Remove a minimumPointList[" << j << "] by binding energy interval limit "
					<< minimumPointsList[j]->height << " " << minimumPointsList[j]->position << endl
					<< " because minimum points are too dense " << endl;
				minimumPointsList.erase(minimumPointsList.begin() + j);
				numRemovedMinimumPoints++;

			}
		}
	}
	//getchar();
	//cout << minimumPointsList.size() << endl;
	return numRemovedMinimumPoints;
}

int removeMinimumPointsByHeightLimit(const vector<double> &bindingEnergy, vector<MaximumPoints*> &maximumPointsList, vector<MinimumPoints*> &minimumPointsList, const vector<double> &smoothedSpectrum)
{
	//Height threshold
	const double heightThreshold = 0.02;

	const int ndata = (int)bindingEnergy.size();
	double hMax = 0.0;
	double hMin = 1e300;
	for (int i = 0; i < ndata; i++) {
		if (smoothedSpectrum[i] > hMax) {
			hMax = smoothedSpectrum[i];
		}
		else {
			if (smoothedSpectrum[i] < hMin) {
				hMin = smoothedSpectrum[i];
			}
		}
	}
	double hStandard = hMax - hMin;
	int numRemovedMinimumPoints = 0;

	for (int i = (int)minimumPointsList.size() - 1; i >= 0; i--) {
		/****************************************************/
		if (minimumPointsList.size() - maximumPointsList.size() > 1) {
			//cout << "i=" << i <<  "minimumPointsList.size()=" << minimumPointsList.size() << ",maximumPointsList.size()" << maximumPointsList.size() << endl;
			cout << "Remove a minimumPointList[" << i + 1 << "] (Exception handling)"
				<< minimumPointsList[i + 1]->height << " " << minimumPointsList[i + 1]->position << endl;
			minimumPointsList.erase(minimumPointsList.begin() + i + 1);
		}
		else if (minimumPointsList.size() - maximumPointsList.size() < 1) {
			cout << "i=" << i << "minimumPointsList.size()=" << minimumPointsList.size() << ",maximumPointsList.size()" << maximumPointsList.size() << endl;
			cout << "Remove a maximumPointList[" << i << "] (Exception handling)"
				<< maximumPointsList[i]->height << " " << maximumPointsList[i]->position << endl;
			maximumPointsList.erase(maximumPointsList.begin() + i);
		}
		/****************************************************/

		double hRatio = 0.0;
		if (i == 0) {
			hRatio = (max(smoothedSpectrum[0], maximumPointsList.front()->height) - minimumPointsList.front()->height) / hStandard;
		}
		else if (i == minimumPointsList.size() - 1) {
			hRatio = (max(smoothedSpectrum[ndata - 1], maximumPointsList.back()->height) - minimumPointsList.back()->height) / hStandard;
		}
		else {
			maximumPointsList[i]->height;
			//cout << minimumPointsList.size() << endl;
			minimumPointsList[i]->height;
			hRatio = (max(maximumPointsList[i]->height, maximumPointsList[i - 1]->height) - minimumPointsList[i]->height) / hStandard;
		}

		//cout << i << " hRatio=" << hRatio << endl;
		if (hRatio < 0.0) {
			//cout << "check2 " << maximumPointsList.back()->id << " " << bindingEnergy[maximumPointsList.back()->id] << endl;
			//cout << "check0 " << smoothedSpectrum[ndata - 1] << " " << maximumPointsList.back()->height << " " << minimumPointsList.back()->height << endl;
			//cout << "check1 " << max(maximumPointsList[i]->height, maximumPointsList[i - 1]->height) 
			//	<< " " << (max(maximumPointsList[i]->height, maximumPointsList[i - 1]->height) - minimumPointsList[i]->height) / hStandard
			//	<< " "  << maximumPointsList[i]->height << " " << maximumPointsList[i - 1]->height << " " << minimumPointsList[i]->height << endl;
			//cout << endl;
			//cout << maximumPointsList[i - 1]->id << " " << minimumPointsList[i]->id << " " << maximumPointsList[i]->id << endl;

			//int ndata = (int)bindingEnergy.size();
			//ofstream ofs("debug.csv");
			//ofs << "id,Binding_Energy(eV),smoothedSpectrum" << endl;

			//for (int j = 0; j < ndata; j++) {
			//	ofs << j << "," << bindingEnergy[j] << "," << smoothedSpectrum[j] << endl;
			//}
			//ofs.close();
			//outputMaximumPoints("maximum", bindingEnergy, smoothedSpectrum, maximumPointsList);
			//outputMinimumPoints("minimum", bindingEnergy, smoothedSpectrum, minimumPointsList);
			//exit(1);
		}
		if (hRatio < heightThreshold) {
			cout << "Remove a minimumPointList[" << i << "] by height limit"
				<< minimumPointsList[i]->height << " " << minimumPointsList[i]->position << endl;
			minimumPointsList.erase(minimumPointsList.begin() + i);
			if (i > 0) {
				maximumPointsList.erase(maximumPointsList.begin() + i - 1);
			}
			numRemovedMinimumPoints++;

		}
	}
	
	//cout << minimumPointsList.size() << endl;
	return numRemovedMinimumPoints;
}

int removeMinimumPointsByHeightLimit2(const vector<double> &bindingEnergy, vector<MaximumPoints*> &maximumPointsList, vector<MinimumPoints*> &minimumPointsList, const vector<double> &smoothedSpectrum)
{
	const double hCriteria = 0.03;  // 0.03 default
	//const double hCriteria = 0.01;  // 0.03 default //shinotsuka20180221, test
	cout << "Height criteria for minimum points deletion= " << hCriteria << endl;

	const int ndata = (int)bindingEnergy.size();
	double hMax = 0.0;
	double hMin = 1e300;
	for (int i = 0; i < ndata; i++) {
		hMax = max(hMax, smoothedSpectrum[i]);
		hMin = min(hMin, smoothedSpectrum[i]);
	}
	double hStandard = hMax - hMin;
	int numRemovedMinimumPoints = 0;

	cout << "Before deletion, minimumPointsList.size()= " << minimumPointsList.size() << endl;

	for (int k = (int)maximumPointsList.size() - 1; k >= 0; k--) {
		//
		// Assuming that there are always two minimum points next to each maximum point
		//
		double hLeft = maximumPointsList[k]->height - minimumPointsList[k]->height;
		double hRight = maximumPointsList[k]->height - minimumPointsList[k+1]->height;

		bool isLeftSmall = (hLeft < hStandard * hCriteria);
		bool isRightSmall = (hRight < hStandard * hCriteria);

		// search a deletion candidate among minimum points
		int deletionIDMin = -1;
		if (isLeftSmall && isRightSmall) {
			//if (hLeft < hRight) { deletionIDMin = k; } else {deletionIDMin = k+1;}
			deletionIDMin = (hLeft < hRight) ? k : k + 1;
		}
		else if (isLeftSmall) {
			deletionIDMin = k;
		}
		else if (isRightSmall) {
			deletionIDMin = k + 1;
		}
		else {
			// Nothing to do
			continue;
		}

		// which maximum point to be removed? next to the minimum point.
		int deletionIDMax = -1;
		if (deletionIDMin == 0) {
			deletionIDMax = 0;
		}
		else if (deletionIDMin == minimumPointsList.size() - 1) {
			deletionIDMax = (int)maximumPointsList.size() - 1;
		}
		else {
			hLeft = maximumPointsList[deletionIDMin - 1]->height - minimumPointsList[deletionIDMin]->height;
			hRight = maximumPointsList[deletionIDMin]->height - minimumPointsList[deletionIDMin]->height;
			deletionIDMax = (hLeft < hRight) ? deletionIDMin - 1 : deletionIDMin;
		}

		// delete a pair of {a minimum point and a maximum point}
		cout << "Delete miminumPointsList[" << deletionIDMin << "]"
			<< ", position= " << minimumPointsList[deletionIDMin]->position 
			<< ", height= " << minimumPointsList[deletionIDMin]->height << endl;

		minimumPointsList.erase(minimumPointsList.begin() + deletionIDMin);
		maximumPointsList.erase(maximumPointsList.begin() + deletionIDMax);
		numRemovedMinimumPoints++;

	}

	cout << "After deletion, minimumPointsList.size()= " << minimumPointsList.size() << endl;

	return numRemovedMinimumPoints;
}

