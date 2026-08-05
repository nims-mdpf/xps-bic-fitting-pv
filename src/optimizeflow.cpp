#define _USE_MATH_DEFINES
#include "optimizeflow.h"
#include "output.h"


/*
* ---------------------------------------------------------------------
* The peak separation of an input diffSpectrum[].
* INPUT
*    ndata:              Number of data points
*    m_SG:               Number of smoothing data points
*    bindingEnergy:      Binding energies (eV)
*    diffSpectrum:       The difference spectrum <-- (raw spectrum[]) - (background[])
*    settings:           Settings
* OUTPUT
*    smoothedSpectrum:   The smoothed spectrum.
*    peakList:           The peak parameters
*    sumSpectrum:        The spectrum superposed by peakList[]
*    chi2:               The chi square between input spectrum and synthesized spectrum
* ---------------------------------------------------------------------
*/
int optimizeFlow(const int ndata, const int m_SG, vector<double> &bindingEnergy,
	vector<double> &spectrum, vector<double> &background, const Settings &settings, 
	vector<PeakParameter*> &peakList, 
	vector<double> &sumSpectrum, double &chi2, double &aic, double &bic, double &reducedChi2,
	double &bgHigher, double &bgLower,
	double &scalingfactor, double &aic_gauss, double &bic_gauss, double &caic_poisson, double & caic_gauss,
	bool &isConverged, int &numTotalOptimizeLoop, int &outerOptimizeLoop)
{
	//Remove peaks outside the energy range.
	int numRemovedPeaks3 = removePeaksByEnergyRange(bindingEnergy, peakList);
	cout << "The number of removed peaks by energy range = " << numRemovedPeaks3 << endl;

	//shinotsuka20180216, output all figures for m_SG
	//if (1) {
	//	cout << "Debug_20180215. Output peak search result for all smoothing points." << endl;
	//	output_result("peak_get_mSG" + to_string(m_SG), bindingEnergy, spectrum, background, peakList, bindingEnergy, spectrum, idHigher, idLower);
	//}

	cout << "Number of initially detected peaks = " << peakList.size() << endl;

	bool isPeakOptimize = true;
	bool isBgOptimize = true;
	bool isFixedLorentzRatio = false;  //Originally, true. But false is better.
	double fixedLorentzRatio = 0.0;
	bool isNumberPeakFixed = true;
	int loopMaxOptimisze = 100;

	const double convergence = 0.001;
	//const double convergence = 0.0001;
	double chi2old = 1e300;
	isConverged = false;

	//shinotsuka20180601
	int numRemovedPeaks1 = 0;
	int numRemovedPeaks2 = 0;
	const double energy_step = abs(bindingEnergy[1] - bindingEnergy[0]);
	numTotalOptimizeLoop = 0;
	outerOptimizeLoop = 1;
	for (outerOptimizeLoop = 1; outerOptimizeLoop <= 20; outerOptimizeLoop++)
	{
		cout << "\nOuter optimization loop = " << outerOptimizeLoop << endl;

		// Optimize peak parameters, while background parameters are fixed
		int numLoopPeak = 0;
		isPeakOptimize = true;
		isBgOptimize = false;
		optimizeParameters(ndata, m_SG, bindingEnergy, spectrum,
			isPeakOptimize, isBgOptimize, isFixedLorentzRatio, fixedLorentzRatio,
			settings, background, peakList, 
			isNumberPeakFixed, loopMaxOptimisze, chi2, reducedChi2,
			bgHigher, bgLower, numLoopPeak);
		numTotalOptimizeLoop += numLoopPeak;

		// Optimize background parameters, while peak parameters are fixed
		int numLoopBG = 0;
		isPeakOptimize = false;
		isBgOptimize = true;
		optimizeParameters(ndata, m_SG, bindingEnergy, spectrum,
			isPeakOptimize, isBgOptimize, isFixedLorentzRatio, fixedLorentzRatio,
			settings, background, peakList, 
			isNumberPeakFixed, loopMaxOptimisze, chi2, reducedChi2,
			bgHigher, bgLower, numLoopBG);
		numTotalOptimizeLoop += numLoopBG;

		//if (outerOptimizeLoop > 10) {
		//	int numLoopBoth = 0;
		//	// Optimize both peak and background parameters
		//	isPeakOptimize = true;
		//	isBgOptimize = true;
		//	optimizeParameters(ndata, m_SG, bindingEnergy, spectrum,
		//		isPeakOptimize, isBgOptimize, isFixedLorentzRatio, fixedLorentzRatio,
		//		settings, background, peakList,
		//		isNumberPeakFixed, loopMaxOptimisze, chi2, reducedChi2,
		//		bgHigher, bgLower, numLoopBoth);
		//	numTotalOptimizeLoop += numLoopBoth;
		//}

		// Remove very small peaks and very narrow peaks
		numRemovedPeaks1 = removeSmallPeaksByWholeArea(peakList, settings.minimumPeakAreaRatio);
		numRemovedPeaks2 = removeSmallPeaksByPeakWidth(peakList, energy_step, settings.minimumPeakFWHM / settings.bindingEnergy_range);

		// When some peak removed, the fitting result gets worse.
		// Then optimize peaks again based on the rest prominent peaks.
		
		// It is converged 
		//     when the number of peaks is the same 
		//     and chi-square is almost the same as the previous model
		isConverged = (numRemovedPeaks1 + numRemovedPeaks2 == 0) && (abs((chi2 - chi2old) / chi2) < convergence);

		if (isConverged) { 
			cout << "\nOuter optimization loop converged after " << outerOptimizeLoop << " times" << endl;
			cout << "chi2old, chi2 = " << chi2old << ", " << chi2 << endl;
			break;
		}

		chi2old = chi2;
	}

	if (!isConverged) {
		cout << "\nOuter optimization loop not converged" << endl;
	}

	cout << "m_SG, numTotalOptimizeLoop, outerOptimizeLoop, isConverged = " << m_SG << ", " << numTotalOptimizeLoop << ", " << outerOptimizeLoop << ", " << isConverged << endl;
	

	make_spectrum(bindingEnergy, peakList, sumSpectrum);

	int num_parameter;
	bool do_scaling = false;
	if (do_scaling) {
		num_parameter = 4 * (int)peakList.size() + 3;
	}
	else {
		num_parameter = 4 * (int)peakList.size() + 2;
	}

	double deviance = calc_deviance_poisson(spectrum, sumSpectrum, background, scalingfactor, do_scaling);
	aic = deviance + 2.0 * num_parameter;
	bic = deviance + (double)num_parameter * log((double)ndata);

	caic_poisson = deviance + 2.0 * num_parameter * ndata / (ndata - num_parameter - 1.0);

	num_parameter = 4 * (int)peakList.size() + 2;
	double deviance_gauss = calc_deviance_gauss(spectrum, sumSpectrum, background); //shinotsuka20180607
	aic_gauss = deviance_gauss + 2.0 * num_parameter;
	bic_gauss = deviance_gauss + (double)num_parameter * log((double)ndata);

	caic_gauss = deviance_gauss + 2.0 * num_parameter * ndata / (ndata - num_parameter - 1.0);

	return 0;
}

void normalize_by_min_and_range(double &v, const double vmin, const double vrange)
{
	v = (v - vmin) / vrange;
}

void denormalize_by_min_and_range(double &v, const double vmin, const double vrange)
{
	v = vmin + vrange * v;
}

void normalize_spectrum_and_parameters(const int ndata, vector<double> &bindingEnergy, vector<double> &spectrum,
	vector<double> &background, vector<PeakParameter*> &peakList,
	vector<MinimumPoints*> &minimumPointsList, double &bgHigher, double &bgLower, 
	double &spectrum_min, double &spectrum_range, double &bindingEnergy_min, double &bindingEnergy_range)
{
	double spectrum_max = *max_element(spectrum.begin(), spectrum.end());
	spectrum_min = *min_element(spectrum.begin(), spectrum.end());
	spectrum_range = spectrum_max - spectrum_min;
	//cout << "spectrum_min, spectrum_max = " << spectrum_min << "," << spectrum_max << endl;
	//cout << *spectrum.begin() << endl;
	//cout << peakList[0]->width << endl;
	for (int i = 0; i < ndata; i++) {
		normalize_by_min_and_range(spectrum[i], spectrum_min, spectrum_range);
		normalize_by_min_and_range(background[i], spectrum_min, spectrum_range);
	}
	for (int j = 0; j < peakList.size(); j++) {
		//normalize_by_min_and_range(peakList[j]->height, spectrum_min, spectrum_range);
		peakList[j]->height /= spectrum_range;
	}
	for (int k = 0; k < minimumPointsList.size(); k++) {
		normalize_by_min_and_range(minimumPointsList[k]->height, spectrum_min, spectrum_range);
	}
	normalize_by_min_and_range(bgHigher, spectrum_min, spectrum_range);
	normalize_by_min_and_range(bgLower, spectrum_min, spectrum_range);

	double bindingEnergy_max = *max_element(bindingEnergy.begin(), bindingEnergy.end());
	bindingEnergy_min = *min_element(bindingEnergy.begin(), bindingEnergy.end());
	bindingEnergy_range = bindingEnergy_max - bindingEnergy_min;
	cout << "bindingEnergy_min, bindingEnergy_max = " << bindingEnergy_min << "," << bindingEnergy_max << endl;
	for (int i = 0; i < ndata; i++) {
		normalize_by_min_and_range(bindingEnergy[i], bindingEnergy_min, bindingEnergy_range);
	}
	for (int j = 0; j < peakList.size(); j++) {
		normalize_by_min_and_range(peakList[j]->position, bindingEnergy_min, bindingEnergy_range);
		peakList[j]->width /= bindingEnergy_range;
	}
	for (int k = 0; k < minimumPointsList.size(); k++) {
		normalize_by_min_and_range(minimumPointsList[k]->position, bindingEnergy_min, bindingEnergy_range);
	}
}

void denormalize_spectrum_and_parameters(const int ndata, vector<double> &bindingEnergy, vector<double> &spectrum,
	vector<double> &background, vector<PeakParameter*> &peakList,
	vector<MinimumPoints*> &minimumPointsList, double &bgHigher, double &bgLower,
	double &spectrum_min, double &spectrum_range, double &bindingEnergy_min, double &bindingEnergy_range)
{
	//cout << "spectrum_min, spectrum_max = " << *min_element(spectrum.begin(), spectrum.end()) << "," << *max_element(spectrum.begin(), spectrum.end()) << endl;
	//cout << *spectrum.begin() << endl;
	//cout << peakList[0]->width << endl;
	for (int i = 0; i < ndata; i++) {
		denormalize_by_min_and_range(spectrum[i], spectrum_min, spectrum_range);
		denormalize_by_min_and_range(background[i], spectrum_min, spectrum_range);
	}
	for (int j = 0; j < peakList.size(); j++) {
		//denormalize_by_min_and_range(peakList[j]->height, spectrum_min, spectrum_range);
		peakList[j]->height *= spectrum_range;
	}
	for (int k = 0; k < minimumPointsList.size(); k++) {
		denormalize_by_min_and_range(minimumPointsList[k]->height, spectrum_min, spectrum_range);
	}
	denormalize_by_min_and_range(bgHigher, spectrum_min, spectrum_range);
	denormalize_by_min_and_range(bgLower, spectrum_min, spectrum_range);

	for (int i = 0; i < ndata; i++) {
		denormalize_by_min_and_range(bindingEnergy[i], bindingEnergy_min, bindingEnergy_range);
	}
	for (int j = 0; j < peakList.size(); j++) {
		denormalize_by_min_and_range(peakList[j]->position, bindingEnergy_min, bindingEnergy_range);
		peakList[j]->width *= bindingEnergy_range;
	}
	for (int k = 0; k < minimumPointsList.size(); k++) {
		denormalize_by_min_and_range(minimumPointsList[k]->position, bindingEnergy_min, bindingEnergy_range);
	}
}


/*
* ---------------------------------------------------------------------
* The Savitzky-Golay method to make a smoothed spectrum and its third derivative.
* INPUT
*    ndata:              Number of data points
*    diffSpectrum:       The difference spectrum <-- (raw spectrum) - (background)
*    numSmooth:          Smoothing data points around a target point.
*                        It must be 3, 5, 7, 9, ... .
* OUTPUT
*    smoothedSpectrum:   The smoothed spectrum.
*    thirdDerivative:    The third derivative of the smoothed spectrum.
* ---------------------------------------------------------------------
*/
int savitzkyGolay(const int ndata, const vector<double> &diffSpectrum, const int numSmooth, vector<double> &smoothedSpectrum, vector<double> &thirdDerivative)
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
	vector<double> wt3(m+1); // Weight factor for smoothed 3rd derivative
	for (int j = 0; j <= m; j++) {
		double temp = (wtFactor - 1.0 - 5.0*j*j);
		wt0[j] = temp / normFactor;
		wt3[j] = - (double)j * temp;
		//cout << "wt0[" << j << "]=" << wt0[j] << "  wt3[" << j << "]=" << wt3[j] << endl;
	}
	int j, i_plus_j, i_minus_j;
	double sum1, sum2;
//#pragma omp parallel for private(j,i_plus_j, i_minus_j,sum1,sum2)
	for (int i = 0; i < ndata; i++) {
		sum1 = wt0[0] * diffSpectrum[i];
		sum2 = 0.0;
//#pragma omp parallel for reduction (+:sum1,sum2)
		for (j = 1; j <= m; j++) {
			i_plus_j = min(ndata - 1, i + j);
			i_minus_j = max(0, i - j);
			sum1 += wt0[j] * (diffSpectrum[i_plus_j] + diffSpectrum[i_minus_j]);
			sum2 += wt3[j] * (diffSpectrum[i_plus_j] - diffSpectrum[i_minus_j]);
		}
		smoothedSpectrum[i] = sum1;
		thirdDerivative[i] = sum2;
	}

	return 0;
}


/*
* ---------------------------------------------------------------------
* Normalize the height of the peaks so that the areas of the original
* spectrum and the smoothed spectrum become equal.
* INPUT
*    ndata:              Number of data points
*    diffSpectrum:       The original different spectrum
*    smoothedSpectrum:   The smoothed spectrum
*    peakList:           The peak parameters
* OUTPUT
*    peakList:           The peak parameters, whose heights are modified
* ---------------------------------------------------------------------
*/
void normalizeHeight(const int ndata, const vector<double> &diffSpectrum, const vector<double> &smoothedSpectrum, vector<PeakParameter*> &peakList)
{
	double areaOrg = (diffSpectrum[0] + diffSpectrum[ndata - 1])*0.5;
	double areaSmoothed = (smoothedSpectrum[0] + smoothedSpectrum[ndata - 1])*0.5;
	for (int i = 1; i < ndata - 1; i++) {
		areaOrg += diffSpectrum[i];
		areaSmoothed += smoothedSpectrum[i];
	}
	double areaRatio = areaSmoothed / areaOrg;
	cout << "areaOrg, areaSmoothed, areaRatio= " << areaOrg << " " << areaSmoothed << " " << areaRatio << endl;
	for (int k = 0; k < peakList.size(); k++) {
		peakList[k]->height *= areaRatio;
	}
}

// Do not use
int adaptiveSmoothing(const int numSmooth, const vector<double> &bindingEnergy, const vector<double> &diffSpectrum, vector<double> &smoothedSpectrum)
{
	double sigman2 = 300.0;

	if (numSmooth % 2 == 0 || numSmooth < 3)
		return 1;

	const int m = (numSmooth - 1) / 2;

	const int ndata = (int)bindingEnergy.size();
	vector<double> sigmax2(ndata);
	vector<double> meanx(ndata);


	for (int i = 0; i < ndata; i++) {
		meanx[i] = 0.0;
		for (int j = -m; j <= m; j++) {
			int i_plus_j = max(min(ndata - 1, i + j), 0);
			meanx[i] += diffSpectrum[i_plus_j];
		}
		meanx[i] /= (double)numSmooth;

		sigmax2[i] = 0.0;
		for (int j = -m; j <= m; j++) {
			int i_plus_j = max(min(ndata - 1, i + j), 0);
			sigmax2[i] += pow(diffSpectrum[i_plus_j] - meanx[i], 2);
		}
		sigmax2[i] /= (double)numSmooth;

		smoothedSpectrum[i] = (sigmax2[i] - sigman2) / sigmax2[i] * (diffSpectrum[i] - meanx[i]) + meanx[i];
	}
	ofstream ofs;
	ofs.open("out_adaptiveSmoothing.csv");
	ofs << "Binding_Energy(eV),diff,smoothedSpectrum,meanx,sigmax2" << endl;
	for (int i = 0; i < ndata; i++) {
		ofs << bindingEnergy[i] << "," << diffSpectrum[i] << "," << smoothedSpectrum[i] << "," << meanx[i] << "," << sigmax2[i] << endl;
	}
	return 0;
}


double calc_scalingfactor(const vector<double> &spectrum, const vector<double> &sumSpectrum,
	const vector<double> &background)
{
	//shinotsuka20180423, ���U�̍Ŗސ���
	double b = 0.0;
	const int ndata = (int)spectrum.size();

	for (int i = 0; i < ndata; i++) {
		double fitted = sumSpectrum[i] + background[i];
		double tmp = spectrum[i] - fitted;
		b += tmp * tmp / fitted;
	}
	b /= ndata;
	double scalingfactor = 1.0 / b;
	cout << "b=" << b << endl;
	cout << "c=" << scalingfactor << endl;
	return scalingfactor;
}

double calc_deviance_poisson(const vector<double> &spectrum, const vector<double> &sumSpectrum,
	const vector<double> &background, double &scalingfactor, const bool do_scaling)
{
	const int ndata = (int)spectrum.size();
	scalingfactor = 1.0;
	if (do_scaling) {
		scalingfactor = calc_scalingfactor(spectrum, sumSpectrum, background);
	}

	//�X�y�N�g���f�[�^�ɃX�P�[�����O���q���|���āA��E�x�̌v�Z
	double deviance = 0.0;
	for (int i = 0; i < ndata; i++) {
		double c_fitted = scalingfactor * (sumSpectrum[i] + background[i]);
		double c_spectrum = scalingfactor * spectrum[i];
		//double c_spectrum = (double)((int)(scalingfactor * spectrum[i]));
		deviance += c_fitted - c_spectrum * log(c_fitted) + lgamma(c_spectrum + 1.0);
		//cout << scalingfactor << "\t" << c_fitted << "\t" << c_spectrum << endl;
	}
	deviance *= 2.0;
	cout << "deviance=" << deviance << endl;
	return deviance;
}

double calc_deviance_gauss(const vector<double> &spectrum, const vector<double> &sumSpectrum,
	const vector<double> &background)
{
	const int ndata = (int)spectrum.size();

	double variance = 0.0;
	for (int i = 0; i < ndata; i++) {
		double tmp = spectrum[i] - sumSpectrum[i] - background[i];
		variance += tmp * tmp;
	}
	variance /= (double)ndata;
	cout << "variance(Gauss)=" << variance << endl;
	double deviance = (double)ndata * (1.0 + log(2.0 * M_PI * variance));
	cout << "deviance(Gauss)=" << deviance << endl;
	return deviance;
}
