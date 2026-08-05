#include "detectpeaks.h"

/*
* ---------------------------------------------------------------------
* Peak detection from a smoothed spectrum and its third derivative.
* INPUT
*    ndata:              Number of data points
*    xx:                 x-axis
*    smoothedSpectrum:   The smoothed spectrum without background
*    thirdDerivative:    The third derivative of the smoothed spectrum.
*    settings:           Settings
* OUTPUT
*    peakList:           The peak parameters
* ---------------------------------------------------------------------
*/
int detectPeaks(const vector<double> &xx, const vector<double> &smoothedSpectrum, vector<double> &thirdDerivative, const Settings &settings, vector<PeakParameter*> &peakList)
{
	// Adjust the smoothed third derivative spectrum.
	adjustThirdDeriv(thirdDerivative);

	// Find extremal indecies of the smoothed spectrum.
	vector<int> extremalIndexList;
	findExtremalIndexList(smoothedSpectrum, thirdDerivative, extremalIndexList);

	// Find peaks
	pickupPeakParameter(xx, smoothedSpectrum, extremalIndexList, peakList, settings);

	return 0;
}

/*
* ---------------------------------------------------------------------
* Adjust the smoothed third derivative spectrum.
* When a value is zero, then modify it to a "very small positive value".
* INPUT/OUTPUT
*    thirdDerivative:    The third derivative of the smoothed spectrum.
* ---------------------------------------------------------------------
*/
void adjustThirdDeriv(vector<double> &thirdDerivative)
{
	const int ndata = (int)thirdDerivative.size();

	bool hasZero_AbsThirdDeriv = false;
	for (int i = 0; i < ndata; i++) {
		if (thirdDerivative[i] == 0.0) {
			hasZero_AbsThirdDeriv = true;
			break;
		}
	}

	if (hasZero_AbsThirdDeriv) {
		double minAbsThirdDeriv = 1e300;
		for (int i = 0; i < ndata; i++) {
			if (thirdDerivative[i] != 0.0) {
				minAbsThirdDeriv = min(minAbsThirdDeriv, abs(thirdDerivative[i]));
			}
		}

		for (int i = 0; i < ndata; i++) {
			if (thirdDerivative[i] == 0.0) {
				//cout << "Third derivative is zero at " << i << ". Change to " << abs(minAbsThirdDeriv * 1e-3 << endl;
				thirdDerivative[i] = abs(minAbsThirdDeriv * 1e-3);
			}
		}
	}
}

/*
* ---------------------------------------------------------------------
* Find extremal indecies of the smoothed spectrum.
* In case of a local maximum, the third derivative crosses zero
* from negative to positive. Then add a negative sign to the index.
* In case of a local minimum, the third derivative crosses zero
* from positive to negative. Then add a positive sign to the index.
*
* INPUT
*    smoothedSpectrum:   The smoothed spectrum.
*    thirdDerivative:    The third derivative of the smoothed spectrum.
* OUTPUT
*    extremalIndexList:  Extremal indecies
* ---------------------------------------------------------------------
*/
void findExtremalIndexList(const vector<double> &smoothedSpectrum, const vector<double> &thirdDerivative, vector<int> &extremalIndexList)
{
	extremalIndexList.clear();
	const int ndata = (int)thirdDerivative.size();

	extremalIndexList.push_back(0); //dummy

	for (int i = 1; i < ndata - 1; i++) {
		if (thirdDerivative[i - 1] * thirdDerivative[i] < 0.0) {
			int iSign = (thirdDerivative[i - 1] < 0.0) ? -1 : 1;
			int iTemp;
			// Take a better choice,  xx[i] or xx[i-1].
			if (abs(smoothedSpectrum[i]) > abs(smoothedSpectrum[i - 1])) {
				iTemp = (i - 1) * iSign;
			}
			else {
				iTemp = i * iSign;
			}
			extremalIndexList.push_back(iTemp);
		}
	}
	extremalIndexList.push_back(ndata - 1); //dummy
}


/*
* ---------------------------------------------------------------------
* Pick up peak parameters
* INPUT
*    xx:                 x-axis
*    smoothedSpectrum:   The smoothed spectrum.
*    extremalIndexList:  Extremal indecies
*    settings:           Settings
* OUTPUT
*    peakList:           The peak parameters
* ---------------------------------------------------------------------
*/
void pickupPeakParameter(const vector<double> &xx, const vector<double> &smoothedSpectrum, const vector<int> &extremalIndexList, vector<PeakParameter*> &peakList, const Settings &settings)
{
	releaseMemoryT(peakList);
	// Maximum height of a smoothed spectrum
	double heightMax = *max_element(smoothedSpectrum.begin(), smoothedSpectrum.end());
	double heightThreshold = heightMax * settings.ratioPeakHeightThreshold;

	for (int k = 1; k < extremalIndexList.size() - 1; k++) {
		int id = abs(extremalIndexList[k]);
		int idm = abs(extremalIndexList[k - 1]);
		int idp = abs(extremalIndexList[k + 1]);
		if (extremalIndexList[k] < 0 && smoothedSpectrum[id] > heightThreshold) {
			double height = smoothedSpectrum[id];
			double position = xx[id];
			double width = 0.5 * max(abs(xx[idm] - xx[id]), abs(xx[id] - xx[idp]));
			double lorentzRatio = 0.0;
			peakList.push_back(new PeakParameter(height, position, width, lorentzRatio));
		}
	}
}




