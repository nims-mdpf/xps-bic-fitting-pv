#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include "optimizeparameters.h"

Spectrum spect;

void Spectrum::setdata(const vector<double> &xi, const vector<double> &yi)
{
	ndata = xi.size();
	x = xi;
	y = yi;
	xmax = *max_element(x.begin(), x.end());
	xmin = *min_element(x.begin(), x.end());
	ymax = *max_element(y.begin(), y.end());
	ymin = *min_element(y.begin(), y.end());
}

/*
 * ---------------------------------------------------------------------
 * Matrix {A^(trans) A} and vector {A^(trans) S}
 * INPUT
 *    a:                  Matrix A
 *	 s:                  Vector s
 *    lambda:             Stabilization factor
 * OUTPUT
 *    ata:                A^(trans) A
 *    ats:                A^(trans) S
 * ---------------------------------------------------------------------
 */
void calculateMatrixAtA(vector<vector<double>> &a, vector<double> &s, vector<vector<double>> &ata, vector<double> &ats)
{
	const int ndata = (int)a.size();
	const int num_parameter = (int)a[0].size();
	// error check for the size of a[][], s[], ata[][], ats[]?

	// Calculate A^(trans) * A
	for (int j = 0; j < num_parameter; j++)
	{
		for (int k = j; k < num_parameter; k++)
		{
			ata[j][k] = 0.0;
			for (int i = 0; i < ndata; i++)
			{
				ata[j][k] += a[i][j] * a[i][k];
			}
			ata[k][j] = ata[j][k];
		}
	}

	// Calculate A^(trans) * S
	for (int j = 0; j < num_parameter; j++)
	{
		ats[j] = 0.0;
		for (int i = 0; i < ndata; i++)
		{
			ats[j] += a[i][j] * s[i];
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * Matrix {A^(trans) A + lambda * D}, where D is a diagonal matrix of {A^(trans) A}
 * INPUT
 *    lambda:             Stabilization factor
 *    ata:                A^(trans) A
 * OUTPUT
 *    ata:                A^(trans) A + lambda * D
 * ---------------------------------------------------------------------
 */
void calculateMatrixAtAD(const double lambda, vector<vector<double>> &ata)
{
	// for (int j = 0; j < (int)ata[0].size(); j++) {
	//	ata[j][j] *= (1.0 + lambda);
	// }

	// shinotsuka20171206
	// Marquardt�@�̊��� r, ���艻Gauss-Newton�@�̊��� (1-r)
	// const double r = 1.0;
	// for (int j = 0; j < (int)ata[0].size(); j++) {
	//	//ata[j][j] = r * (1.0 + lambda) * ata[j][j] + (1.0 - r) * (ata[j][j] + lambda);
	//	ata[j][j] += lambda * (r * ata[j][j] + 1.0 - r);
	// }

	// shinotsuka20180209, test
	for (int j = 0; j < (int)ata[0].size(); j++)
	{
		ata[j][j] += lambda * (ata[j][j] + 1.0);
	}
}

/*
 * ---------------------------------------------------------------------
 * Solve a matrix equation [A^(trans)A + lambda * D][e] = [-A^(trans) * S]
 * INPUT
 *    ata:                A^(trans)A + lambda * D, where D is a diagonal matrix of A^(trans) A
 *    ats:                A^(trans) * S
 * OUTPUT
 *    e:                  e
 * ---------------------------------------------------------------------
 */
void solveMatrixEquation(vector<vector<double>> &ata, vector<double> &ats, vector<double> &e)
{
	// Prepare to use dgesv_()
	const int n = (int)ata.size();
	int *piv = new int[n];
	double *aa = new double[n * n];
	double *b = new double[n];
	int info;
	for (int j = 0; j < n; j++)
	{
		for (int k = 0; k < n; k++)
		{
			aa[n * k + j] = ata[j][k];
		}
		b[j] = -ats[j];
	}
	dgesv_(n, 1, (double **)aa, n, piv, (double **)b, n, info);
	// shinotsuka20171206
	// if (info != 0) {
	//	cout << "Error at solveMatrixEquation(), dgesv_() info= " << info << endl;
	// }
	for (int j = 0; j < n; j++)
	{
		e[j] = b[j];
	}
	// shinotsuka20190730
	delete[] piv;
	delete[] aa;
	delete[] b;
}

void restrictPeaks(vector<PeakParameter *> &peakList, const Settings &settings)
{
	for (int j = 0; j < peakList.size(); j++)
	{
		// shinotsuka20180209
		// peakList[j]->height = max(peakList[j]->height, 0.0);
		// peakList[j]->width = max(peakList[j]->width, 0.0);
		peakList[j]->height = max(peakList[j]->height, 1e-15);
		peakList[j]->width = max(peakList[j]->width, 1e-15);
		// peakList[j]->width = min(peakList[j]->width, 5.0);   //shinotsuka20180219, test. not used
		peakList[j]->lorentzRatio = max(peakList[j]->lorentzRatio, settings.lorentzRatioCutoffs[0]);
		peakList[j]->lorentzRatio = min(peakList[j]->lorentzRatio, settings.lorentzRatioCutoffs[1]);
	}
}

int removeInfinitesimalPeaks(vector<PeakParameter *> &peakList)
{
	int numRemovedPeaks = 0;
	const double threshold = 1e-14; // Threshold of the peak area to be removed

	if (peakList.size() == 0)
	{
		return 0;
	}

	for (int j = (int)peakList.size() - 1; j >= 0; j--)
	{
		if (peakList[j]->width <= threshold || peakList[j]->height <= threshold)
		{
			cout << "Remove a peakList[" << j << "] "
				 << peakList[j]->height << " " << peakList[j]->position
				 << " " << peakList[j]->width << " " << peakList[j]->lorentzRatio
				 << " because its height or width is infinitesimal small." << endl;
			peakList.erase(peakList.begin() + j);
			numRemovedPeaks++;
		}
		if (peakList.size() == 1)
			break;
	}
	return numRemovedPeaks;
}

// shinotsuka20180530
//
//  Remove very small peaks relative to a whole peak area
//  threshold : Threshold of the peak area to be removed
//
int removeSmallPeaksByWholeArea(vector<PeakParameter *> &peakList, const double threshold = 0.01)
{
	if (peakList.size() == 1)
	{
		return 0;
	}
	int numRemovedPeaks = 0;

	// Remove a peak whose area is much less than the maximum peak area.
	vector<double> areas(peakList.size());
	double stdArea = 0.0;
	for (int j = 0; j < peakList.size(); j++)
	{
		areas[j] = peakList[j]->height * peakList[j]->width * (peakList[j]->lorentzRatio * M_PI + (1.0 - peakList[j]->lorentzRatio) * sqrt(M_PI / log(2.0)));
		// stdArea = max(stdArea, areas[j]);
		stdArea += areas[j];
	}
	for (int j = (int)peakList.size() - 1; j >= 0; j--)
	{
		if (areas[j] / stdArea < threshold)
		{
			cout << "Remove a peakList[" << j << "] "
				 << peakList[j]->height << " " << peakList[j]->position
				 << " " << peakList[j]->width << " " << peakList[j]->lorentzRatio
				 << " because its area " << areas[j] << " / " << stdArea << " (whole area) < " << threshold << " (threshold)" << endl;
			peakList.erase(peakList.begin() + j);
			numRemovedPeaks++;
		}
		if (peakList.size() == 1)
			break;
	}
	return numRemovedPeaks;
}

// shinotsuka20180601
//
//  Remove very narrow peaks
//
int removeSmallPeaksByPeakWidth(vector<PeakParameter *> &peakList, const double energy_step, const double threshold_fwhm = 0.02)
{
	if (peakList.size() == 1)
	{
		return 0;
	}
	int numRemovedPeaks2 = 0;
	for (int i = (int)peakList.size() - 1; i >= 0; i--)
	{
		if (peakList[i]->width < 0.5 * max(energy_step, threshold_fwhm))
		{
			cout << "Remove a peakList[" << i << "] by HWHM "
				 << peakList[i]->height << " " << peakList[i]->position
				 << " " << peakList[i]->width << " " << peakList[i]->lorentzRatio
				 << " because HWHM < threshold width " << energy_step << " " << threshold_fwhm << endl;
			peakList.erase(peakList.begin() + i);
			numRemovedPeaks2++;
		}
		if (peakList.size() == 1)
			break;
	}
	return numRemovedPeaks2;
}

int removeSmallPeaks(vector<PeakParameter *> &peakList)
{
	if (peakList.size() == 1)
	{
		return 0;
	}
	int numRemovedPeaks = 0;
	double threshold = 0.01; // Threshold of the peak area to be removed

	// Remove a peak whose area is much less than the maximum peak area.
	vector<double> areas(peakList.size());
	double stdArea = 0.0;
	for (int j = 0; j < peakList.size(); j++)
	{
		areas[j] = peakList[j]->height * peakList[j]->width * (peakList[j]->lorentzRatio * M_PI + (1.0 - peakList[j]->lorentzRatio) * sqrt(M_PI / log(2.0)));
		stdArea = max(stdArea, areas[j]);
	}
	for (int j = (int)peakList.size() - 1; j >= 0; j--)
	{
		if (areas[j] / stdArea < threshold)
		{
			cout << "Remove a peakList[" << j << "] "
				 << peakList[j]->height << " " << peakList[j]->position
				 << " " << peakList[j]->width << " " << peakList[j]->lorentzRatio
				 << " because its area " << areas[j] << " / " << stdArea << " (std. area) < " << threshold << " (threshold)" << endl;
			peakList.erase(peakList.begin() + j);
			numRemovedPeaks++;
		}
		if (peakList.size() == 1)
			break;
	}
	return numRemovedPeaks;
}

int removePeaksByHWHM(vector<PeakParameter *> &peakList, const double threshold_width)
{
	if (peakList.size() == 1)
	{
		return 0;
	}
	int numRemovedPeaks2 = 0;
	for (int i = (int)peakList.size() - 1; i >= 0; i--)
	{
		// cout << i << " " << areas[j] << endl;
		if (peakList[i]->width < threshold_width)
		{
			cout << "Remove a peakList[" << i << "] by HWHM"
				 << peakList[i]->height << " " << peakList[i]->position
				 << " " << peakList[i]->width << " " << peakList[i]->lorentzRatio << endl
				 << " because the HWHM is smaller than threshold width " << endl;
			peakList.erase(peakList.begin() + i);
			numRemovedPeaks2++;
		}
		if (peakList.size() == 1)
			break;
	}
	return numRemovedPeaks2;
}

int removePeaksByEnergyRange(const vector<double> &bindingEnergy, vector<PeakParameter *> &peakList)
{
	if (peakList.size() == 1)
	{
		return 0;
	}
	int numRemovedPeaks3 = 0;
	for (int i = (int)peakList.size() - 1; i >= 0; i--)
	{
		if (peakList[i]->position > bindingEnergy.front() || peakList[i]->position < bindingEnergy.back())
		{
			cout << "Remove a peakList[" << i << "] by energy range "
				 << peakList[i]->height << " " << peakList[i]->position
				 << " " << peakList[i]->width << " " << peakList[i]->lorentzRatio << endl
				 << " because its position " << peakList[i]->position << " is out of chiSquare calculation range " << endl;
			peakList.erase(peakList.begin() + i);
			numRemovedPeaks3++;
		}
		if (peakList.size() == 1)
			break;
	}
	return numRemovedPeaks3;
}

/*
 * ---------------------------------------------------------------------
 * Jacobian Matrix A with bgHigher and bgLower
 * INPUT
 *    ndata:              Number of data points
 *    x:                  Binding energies (eV)
 *    peakList:           The peak parameters
 *    bgHigher:           background height at the left
 *    bgLower:            background height at the right
 * OUTPUT
 *    a:                  Jacobian Matrix A
 * ---------------------------------------------------------------------
 */
void calculateMatrixA2(const int ndata, const vector<double> &x, vector<PeakParameter *> &peakList, vector<double> &qpq, vector<vector<double>> &a)
{
	int cnt = 0;

	for (int i = 0; i < ndata; i++)
	{
		cnt = -1;
		for (int k = 0; k < peakList.size(); k++)
		{
			PeakParameter pp = *peakList[k];

			double mu = pp.lorentzRatio;
			double gauss = gaussFunc(x[i], pp.position, pp.width);
			double lorentz = lorentzFunc(x[i], pp.position, pp.width);

			double xcw = (x[i] - pp.position) / pp.width;
			double xcw2 = xcw / pp.width;
			double c_dLdc = 2.0 * xcw2;
			double c_dGdc = 2.0 * log(2.0) * xcw2;
			// double c_dGdw = c_dGdc * xcw2;

			double dLdc = c_dLdc * pow(lorentz, 2);
			double dLdw = xcw * dLdc;

			double dGdc = c_dGdc * gauss;
			double dGdw = xcw * dGdc;

			double dVdA = mu * lorentz + (1.0 - mu) * gauss;
			double dVdc = pp.height * (mu * dLdc + (1.0 - mu) * dGdc);
			double dVdw = pp.height * (mu * dLdw + (1.0 - mu) * dGdw);
			double dVdm = pp.height * (lorentz - gauss);

			a[i][++cnt] = dVdA;
			a[i][++cnt] = dVdc;
			a[i][++cnt] = dVdw;
			a[i][++cnt] = dVdm;
		}

		a[i][++cnt] = qpq[i];
		a[i][++cnt] = 1.0 - qpq[i];
	}
}

/*
 * ---------------------------------------------------------------------
 * Jacobian Matrix A for peak parameters
 * INPUT
 *    ndata:              Number of data points
 *    x:                  Binding energies (eV)
 *    peakList:           The peak parameters
 * OUTPUT
 *    a:                  Jacobian Matrix A
 * ---------------------------------------------------------------------
 */
void calculateMatrixA_for_peak(const int ndata, const vector<double> &x, const vector<PeakParameter *> &peakList, vector<vector<double>> &a)
{
	int cnt = 0;

	for (int i = 0; i < ndata; i++)
	{
		cnt = -1;
		for (int k = 0; k < peakList.size(); k++)
		{
			PeakParameter pp = *peakList[k];

			double mu = pp.lorentzRatio;
			double gauss = gaussFunc(x[i], pp.position, pp.width);
			double lorentz = lorentzFunc(x[i], pp.position, pp.width);

			double xcw = (x[i] - pp.position) / pp.width;
			double xcw2 = xcw / pp.width;
			double c_dLdc = 2.0 * xcw2;
			double c_dGdc = 2.0 * log(2.0) * xcw2;
			// double c_dGdw = c_dGdc * xcw2;

			double dLdc = c_dLdc * pow(lorentz, 2);
			double dLdw = xcw * dLdc;

			double dGdc = c_dGdc * gauss;
			double dGdw = xcw * dGdc;

			double dVdA = mu * lorentz + (1.0 - mu) * gauss;
			double dVdc = pp.height * (mu * dLdc + (1.0 - mu) * dGdc);
			double dVdw = pp.height * (mu * dLdw + (1.0 - mu) * dGdw);
			double dVdm = pp.height * (lorentz - gauss);

			a[i][++cnt] = dVdA;
			a[i][++cnt] = dVdc;
			a[i][++cnt] = dVdw;
			a[i][++cnt] = dVdm;
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * Jacobian Matrix A for Shirley background parameters
 * INPUT
 *    ndata:              Number of data points
 *    x:                  Binding energies (eV)
 *    qpq:                Cummulative integration Q/(P+Q) for Shirley background
 * OUTPUT
 *    a:                  Jacobian Matrix A
 * ---------------------------------------------------------------------
 */
void calculateMatrixA_for_background_Shirley(const int ndata, const vector<double> &qpq, vector<vector<double>> &a)
{
	for (int i = 0; i < ndata; i++)
	{
		a[i][0] = qpq[i];
		a[i][1] = 1.0 - qpq[i];
	}
}

/*
 * ---------------------------------------------------------------------
 * Jacobian Matrix A for linear background parameters
 * INPUT
 *    ndata:              Number of data points
 * OUTPUT
 *    a:                  Jacobian Matrix A
 * ---------------------------------------------------------------------
 */
void calculateMatrixA_for_background_Linear(vector<vector<double>> &a)
{
	const int ndata = (int)a.size();
	for (int i = 0; i < ndata; i++)
	{
		a[i][0] = ((double)(i - ndata + 1)) / ((double)(0 - ndata + 1));
		a[i][1] = 1.0 - a[i][0];
	}

	// for (int i = 0; i < ndata; i++) {
	//	printf("%f,%f\n", a[i][0], a[i][1]);
	// }
	// printf("dbg\n"); getchar();
}

/*
 * ---------------------------------------------------------------------
 * Make a background with bgHigher, bgLower and ffp[]
 * INPUT
 *    bgHigher:           background height at the higher binding energy side
 *    bgLower:            background height at the lower binding energy side
 *    ffp:                superposition of the peaks free from any background
 * OUTPUT
 *    qpq:                Q(x)/(P(x)+Q(x)) in the Shirley equation
 *    background:         The Shirley background
 * ---------------------------------------------------------------------
 */
int make_background(const double &bgHigher, const double &bgLower,
					const vector<double> &ffp, const Settings &settings,
					vector<double> &qpq, vector<double> &background)
{
	int ret = -1;
	if (settings.backgroundType == "Shirley")
	{
		ret = make_background_Shirley(bgHigher, bgLower, ffp, qpq, background);
	}
	else if (settings.backgroundType == "Linear")
	{
		// Linear background using two points {bgHigher, bgLower}
		ret = make_background_Linear(bgHigher, bgLower, background);
	}
	else
	{
		cout << "Error at make_background(). Unknown settings.backgroundType = " << settings.backgroundType << endl;
		exit(1);
	}
	return ret;
}

int make_background_Shirley(const double &bgHigher, const double &bgLower,
							const vector<double> &ffp, vector<double> &qpq, vector<double> &background)
{
	const int ndata = (int)ffp.size();
	if (ndata == 0)
	{
		return 1;
	}
	qpq.resize(ndata);
	background.resize(ndata);
	vector<double> areaQ(ndata);
	areaQ[ndata - 1] = 0.0;
	for (int i = ndata - 2; i >= 0; i--)
	{
		areaQ[i] = areaQ[i + 1] + (ffp[i] + ffp[i + 1]) * 0.5;
	}

	double sumNew = areaQ[0];
	// Background as functions of bgHigher and bgLower
	for (int i = 0; i < ndata; i++)
	{
		qpq[i] = areaQ[i] / sumNew;
		background[i] = (bgHigher - bgLower) * qpq[i] + bgLower;
	}
	return 0;
}

int make_background_Linear(const double &yH, const double &yL, vector<double> &background)
{
	const int ndata = (int)background.size();
	for (int i = 0; i < ndata; i++)
	{
		background[i] = ((double)(-i) * yL + (double(i - ndata + 1)) * yH) / (double)(-ndata + 1);
	}
	return 0;
}

/*
 * ---------------------------------------------------------------------
 * Calculate a chi square value.
 * INPUT
 *    x[]:                Binding energies (eV)
 *    spectrum[]:         Raw spectrum
 *    peakList[]:         The peak parameters
 *    bgHigher:           background height at the higher binding energy side
 *    bgLower:            background height at the lower binding energy side
 * OUTPUT
 *    qpq[]:              Q(x)/(P(x)+Q(x)) in the Shirley equation
 *    s[]:                A residual spectrum
 *    newChiSquare:       The chi square value
 * ---------------------------------------------------------------------
 */
int calcChiSquare(const vector<double> &x, const vector<double> &spectrum,
				  const vector<PeakParameter *> &peakList,
				  const double &bgHigher, const double &bgLower, const Settings &settings,
				  vector<double> &background, vector<double> &qpq, vector<double> &s,
				  double &newChiSquare, double &newReducedChiSquare)
{
	const int ndata = (int)x.size();

	if (ndata != spectrum.size())
	{
		cout << "Error. Invalid array size, x[] and spectrum[]. " << ndata << " " << spectrum.size() << endl;
		return 1;
	}

	// Update a superposition of the peaks free from any background
	vector<double> ffp;
	if (make_spectrum(x, peakList, ffp))
	{
		cout << "Error. Failed to make a superposed spectrum." << endl;
		return 1;
	}
	if (make_background(bgHigher, bgLower, ffp, settings, qpq, background))
	{
		cout << "Error. Failed to make a background." << endl;
		return 1;
	}

	// Vector s[]
	s.resize(ndata);

	for (int i = 0; i < ndata; i++)
	{
		s[i] = (ffp[i] + background[i]) - spectrum[i];
	}

	newChiSquare = 0.0;
	// for (int i = idHigher; i <= idLower; i++) {
	for (int i = 0; i < ndata; i++)
	{
		newChiSquare += pow(s[i], 2.0) / (ffp[i] + background[i]);
	}

	// int ndata_for_chi2 = (idLower - idHigher + 1);
	int ndata_for_chi2 = ndata;
	int ndf_for_chi2 = ndata_for_chi2 - (4 * (int)peakList.size() + 2);
	newReducedChiSquare = newChiSquare / ndf_for_chi2;

	// This is time consuming.
	// outputChiSquare(x, spectrum, ffp, background, s);

	// if (newChiSquare < 0.0) {
	//	cout << "dbg20190809b, negative chi-square" << endl;
	//	cout << "bgHigher, bgLower = " << bgHigher << ", " << bgLower << endl;
	//	for (int i = 0; i < ndata; i++) {
	//		cout << s[i] << " " << ffp[i] << " " << background[i] << endl;
	//	}
	// }

	return 0;
}

void outputChiSquare(const vector<double> &x, const vector<double> &spectrum, const vector<double> &ffp, const vector<double> &bg, const vector<double> &s)
{
	static int call_count = 0;
	const int ndata = (int)x.size();
	string ofname;
	ofstream ofs;
	stringstream ss;

	ss << "x,spectrum,background,fitted,chiSquare" << endl;
	string kwd = "output_chiSquare";

	ofname = kwd + ".csv";

	// shinotsuka20181122
	char ofnamec[1024];
	sprintf(ofnamec, "output_chiSquare_%05d.csv", ++call_count);
	ofname = string(ofnamec);
	// cout << ofname << endl;
	// return;

	ofs.open(ofname);
	for (int i = 0; i < ndata; i++)
	{
		ss << x[i]
		   << "," << spectrum[i]
		   << "," << bg[i]
		   << "," << ffp[i] + bg[i]
		   << "," << pow(s[i] / spectrum[i], 2.0) << endl;
	}
	ss << endl;

	if (ofs)
	{
		cout << ofname << ", call_count=" << call_count << endl;
		ofs << ss.str();
		ofs.close();
	}
	else
	{
		cout << "\n!!! Failed to create a file: " << ofname << " !!!\n"
			 << endl;
	}
}

int output_history(const int m_SG, vector<double> x, vector<double> y, vector<double> background, vector<PeakParameter *> peakList, const Settings &settings)
{
	string ofname;
	ofstream ofs;
	static int call_count = 0;
	const int ndata = (int)x.size();
	const int originalNdata = (int)x.size();
	cout << "<< m_SG = " << m_SG << "  history No. " << call_count << ">>" << endl;

	// Output a list of peak parameters and their areas
	cout << "Number of peaks=" << peakList.size() << endl;

	if (peakList.size() == 0)
	{
		cout << "Skip to output file because there is no peaks." << endl;
		return 0;
	}

	// shinotsuka20181205a
	//  Before output the result, denormalize spectra and parameters
	// denormalize_spectra_and_parameters(x, y, background, peakList,
	//	x, y, smoothedSpectrum, settings);

	stringstream ss;
	ss.str("");

	// Output fitted spectra
	vector<double> sumSpectrum;
	vector<vector<double>> peaks;
	make_spectrum_and_peak(x, peakList, sumSpectrum, peaks);

	ss << "Binding_Energy(eV),spectrum,background,fitting";
	for (int j = 0; j < peakList.size(); j++)
	{
		ss << ",peak[" << j << "]";
	}
	ss << endl;

	for (int i = 0; i < x.size(); i++)
	{
		ss << x[i]
		   << "," << y[i]
		   << "," << background[i]
		   << "," << background[i] + sumSpectrum[i];
		for (int j = 0; j < peakList.size(); j++)
		{
			ss << "," << background[i] + peaks[j][i];
		}
		ss << endl;
	}

	// ofname = "result_" + kwd + ".csv";
	char ofnamec[1024];
	// sprintf(ofnamec, "optimize_history_No_%05d.csv", ++call_count);
	sprintf(ofnamec, "optimize_history_mSG_%04d_No_%05d.csv", m_SG, ++call_count);
	ofname = string(ofnamec);

	ofs.open(ofname);
	if (ofs)
	{
		cout << "Output to " + ofname << endl;
		ofs << ss.str();
		ofs.close();
	}
	else
	{
		cout << "\n!!! Failed to create a file: " << ofname << " !!!\n"
			 << endl;
	}
	ss.str("");

	return 0;
}

/*
 * ---------------------------------------------------------------------
 * Copy from {peakList_in, bgHigher_in, bgLower_in}
 *        to {peakList_out, bgHigher_out, bgLower_out}
 * INPUT
 *    peakList_in[]:         Peak parameters
 *    bgHigher_in:           Background height at the higher binding energy side
 *    bgLower_in:            Background height at the lower binding energy side
 * OUTPUT
 *    peakList_out[]:        Peak parameters
 *    bgHigher_out:          Background height at the higher binding energy side
 *    bgLower_out:           Background height at the lower binding energy side
 * ---------------------------------------------------------------------
 */
void copy_peakparameter(const vector<PeakParameter *> &peakList_in, const double bgHigher_in, const double bgLower_in,
						vector<PeakParameter *> &peakList_out, double &bgHigher_out, double &bgLower_out)
{
	// shinotsuka20190807
	copyVectorT(peakList_in, peakList_out);
	bgHigher_out = bgHigher_in;
	bgLower_out = bgLower_in;
}

/*
 * ---------------------------------------------------------------------
 * Peak optimization by the Levenberg-Marquard method (damped Gauss-Newton method)
 * INPUT
 *    ndata:                Number of data points
 *    x[]:                  Binding energies (eV)
 *    spectrum[]:           Raw spectrum
 *    isFixedBgHigher:      True if the background height at the higher binding energy side is fixed
 *    isFixedBgLower:       True if the background height at the lower binding energy side is fixed
 *    isFixedLorentzRatio:  True if the lorentzRatio of each peak is fixed
 *    fixedLorentzRatio:    Fixed lorentzRatio value when isFixedLorentzRatio = true
 *    background[]:         Initial background determined by iterative Shirley method
 *    peakList[]:           The initial peak parameters
 *    isNumberPeakFixed:    True if the number of peaks keeps unchanged even though a small peak appears.
 *    loopMaxOptimisze:     The maximum number of loops for the optimization
 * INPUT/OUTPUT
 *    peakList[]:           Optimized peak parameters
 *    background[]:         Optimized background determined by active Shirley method
 *    optimizedChiSquare:   Optimized chi^2
 *    bgHigher:             Background height at the higher binding energy side
 *    bgHigher:             Background height at the lower binding energy side
 * ---------------------------------------------------------------------
 */
int optimizeParameters(const int ndata, const int m_SG, const vector<double> &x, const vector<double> &spectrum,
					   const bool isPeakOptimize, const bool isBgOptimize, const bool isFixedLorentzRatio, const double fixedLorentzRatio,
					   const Settings &settings, vector<double> &background, vector<PeakParameter *> &peakList,
					   bool isNumberPeakFixed, int loopMaxOptimisze,
					   double &optimizedChiSquare, double &optimizedReducedChiSquare,
					   double &bgHigher, double &bgLower, int &loop)
{
	cout << "optimizeParameters()" << endl;
	const double convergence = 0.001; // 0.01
	// const double convergence = 0.0001;
	// const double convergence = 0.0001;
	double oldChiSquare = 1e300;
	double newChiSquare = 1e300;
	double oldReducedChiSquare = 1e300; // shinotsuka20180220
	double newReducedChiSquare = 1e300; // shinotsuka20180220

	vector<int> outputHistoryList{2, 216, 608, 1600, 2880}; // Indicate a list of m_SG to output optimization history
	outputHistoryList.clear();								// Usually, not output optimization history

	// const int num_parameter = (isFixedBg) ? ((int)peakList.size() * 4) : 2;
	//  Number of parameters is
	//      4K when optimizing peak parameters where K is the number of peaks
	//          The order of a parameter vector is
	//              e[0-3]:           {height, position, width, lorentzRatio} of peakList[0]
	//              e[4-7]:           {height, position, width, lorentzRatio} of peakList[1]
	//              ......
	//              e[(4K-4)-(4K-1)]: {height, position, width, lorentzRatio} of peakList[K-1]
	//      2 when optimizing background parameters
	//          The order of a parameter vector is
	//              e[0]:           background height at the higher binding energy
	//              e[1]:           background height at the lower binding energy
	int num_parameter = 0;
	if (isPeakOptimize)
	{
		num_parameter += (int)peakList.size() * 4;
	}
	if (isBgOptimize)
	{
		num_parameter += 2;
	}

	if (num_parameter == 0)
	{
		cout << "Warning: There is no parameter to be optimized." << endl;
		return 0;
	}

	vector<vector<double>> a(ndata, vector<double>(num_parameter));
	// Rectangular matrix A[ndata][num_parameter]
	vector<vector<double>> ata(num_parameter, vector<double>(num_parameter));
	// Square matrix {A^(trans) A}[num_parameter][num_parameter]

	vector<double> s(ndata);		   // Vector s[ndata], residual spectrum
	vector<double> ats(num_parameter); // Vector {A^(trans) s}
	vector<double> e(num_parameter);   // The steepest descent direction
									   // which is solved by [A^(trans) A + lambda D] e = A^(trans) s
	vector<double> qpq(ndata);		   // Cummulative integration Q(x)/(P(x)+Q(x)) in the Shirley's equation

	if (ndata <= 0)
	{
		cout << "Error. ndata <= 0. " << ndata << endl;
		return 1;
	}
	if (ndata != x.size())
	{
		cout << "Error. x.size() is not ndata. " << x.size() << " " << ndata << endl;
		return 1;
	}
	if (ndata != spectrum.size())
	{
		cout << "Error. spectrum.size() is not ndata." << spectrum.size() << " " << ndata << endl;
		return 1;
	}

	cout << "isPeakOptimize, isBgOptimize, isFixedLorentzRatio, fixedLorentzRatio= "
		 << isPeakOptimize << " " << isBgOptimize << " " << isFixedLorentzRatio << " " << fixedLorentzRatio << endl;

	if (calcChiSquare(x, spectrum, peakList, bgHigher, bgLower, settings, background, qpq, s, oldChiSquare, oldReducedChiSquare))
	{
		cout << "Error. Failed to calculate newChiSquare." << endl;
		return 1;
	}
	vector<PeakParameter *> peakList_saved;
	double bgHigher_saved = 0.0;
	double bgLower_saved = 0.0;

	copy_peakparameter(peakList, bgHigher, bgLower, peakList_saved, bgHigher_saved, bgLower_saved);

	bool isPeakListUpdated = true;
	double lambda = 1.0;
	const double lambda_ascend = 5.0; // 3.0
	const double lambda_decend = 0.1;
	cout << "lambda,lambda_ascend,lambda_decend= " << lambda << "," << lambda_ascend << "," << lambda_decend << endl;
	int countSlowConvergence = 0;

	cout << ",loop,lambda,chi2_old,chi2_new,red.chi2" << endl;

	loop = 1;

	for (loop = 1; loop <= loopMaxOptimisze; loop++)
	{
		if (isPeakListUpdated)
		{
			// When a peakList was updated, calculate matrix A,
			// vector qpq, matrix {A^(trans) A} and vector {A^(trans) s}.
			// On the other hand, when a peakList was not updated,
			// we can use the previous values.

			// calculateMatrixA2(ndata, x, peakList, qpq, a);

			if (isPeakOptimize && !isBgOptimize)
			{
				calculateMatrixA_for_peak(ndata, x, peakList, a);
			}
			else if (isBgOptimize && !isPeakOptimize)
			{
				if (settings.backgroundType == "Shirley")
				{
					calculateMatrixA_for_background_Shirley(ndata, qpq, a);
				}
				else if (settings.backgroundType == "Linear")
				{
					// Jacobian Matrix A for linear background parameters
					calculateMatrixA_for_background_Linear(a);
				}
				else
				{
					cout << "Error at optimizeParameters(). Unknown settings.backgroundType = " << settings.backgroundType << endl;
					exit(1);
				}
			}
			else if (isPeakOptimize && isBgOptimize)
			{
				// test shinotsuka20190822
				if (settings.backgroundType == "Shirley")
				{
					calculateMatrixA2(ndata, x, peakList, qpq, a);
				}
			}
			else
			{
			}
		}

		calculateMatrixAtA(a, s, ata, ats);

		// Calculate matrix {A^(trans) A + lambda D}
		calculateMatrixAtAD(lambda, ata);

		// Solve a matrix equation [A^(trans) A + lambda D] e = A^(trans) s
		solveMatrixEquation(ata, ats, e);

		// Update parameters;  p[] <- p[] + e[]
		if (isPeakOptimize && !isBgOptimize)
		{
			int cnt = -1;
			for (int j = 0; j < peakList.size(); j++)
			{
				peakList[j]->height += e[++cnt];
				peakList[j]->position += e[++cnt];
				peakList[j]->width += e[++cnt];
				if (isFixedLorentzRatio)
				{
					peakList[j]->lorentzRatio = fixedLorentzRatio;
				}
				else
				{
					peakList[j]->lorentzRatio += e[++cnt];
				}
			}
			restrictPeaks(peakList, settings);
		}
		else if (isBgOptimize && !isPeakOptimize)
		{
			bgHigher += e[0];
			bgHigher = max(bgHigher, 1e-15);
			bgLower += e[1];
			bgLower = max(bgLower, 1e-15);
		}
		else if (isPeakOptimize && isBgOptimize)
		{
			int cnt = -1;
			for (int j = 0; j < peakList.size(); j++)
			{
				peakList[j]->height += e[++cnt];
				peakList[j]->position += e[++cnt];
				peakList[j]->width += e[++cnt];
				if (isFixedLorentzRatio)
				{
					peakList[j]->lorentzRatio = fixedLorentzRatio;
				}
				else
				{
					peakList[j]->lorentzRatio += e[++cnt];
				}
			}
			bgHigher += e[++cnt];
			bgHigher = max(bgHigher, 1e-15);
			bgLower += e[++cnt];
			bgLower = max(bgLower, 1e-15);
			restrictPeaks(peakList, settings);
		}
		else
		{
		}

		// restrictPeaks(peakList, settings);

		// Caluculate a chi square for updated parameters
		if (calcChiSquare(x, spectrum, peakList, bgHigher, bgLower, settings, background, qpq, s, newChiSquare, newReducedChiSquare))
		{
			cout << "Error. Failed to calculate newChiSquare." << endl;
			return 1;
		}
		// cout << loop << "," << lambda << "," << oldChiSquare << "," << newChiSquare << "," << newReducedChiSquare << endl;
		cout << "[Opt]," << loop << "," << lambda << "," << oldChiSquare << "," << newChiSquare << "," << newReducedChiSquare << endl;

		// Update a damping factor lambda
		if (abs(newChiSquare) < abs(oldChiSquare))
		{
			countSlowConvergence = 0;
			// Check a convergence
			if (abs((oldChiSquare - newChiSquare) / newChiSquare) < convergence)
			{
				break;
			}
			// Update peakList, bgHigher, bgLower and chiSquare.
			// Let lambda small, to achieve a fast convergence.
			lambda *= lambda_decend;

			if (isNumberPeakFixed == false)
			{
				// Remove small peaks and update chiSquare
				if (removeSmallPeaks(peakList))
				{
					if (calcChiSquare(x, spectrum, peakList, bgHigher, bgLower, settings, background, qpq, s, newChiSquare, newReducedChiSquare))
					{
						cout << "Error. Failed to calculate newChiSquare." << endl;
						return 1;
					}
					cout << loop << "," << lambda << "," << oldChiSquare << "," << newChiSquare << "," << newReducedChiSquare << endl;
				}
			}

			oldChiSquare = newChiSquare;
			oldReducedChiSquare = newReducedChiSquare; // shinotsuka20180220
			copy_peakparameter(peakList, bgHigher, bgLower, peakList_saved, bgHigher_saved, bgLower_saved);
			isPeakListUpdated = true;

			// shinotsuka20190208, output update history for movie
			// output_history(m_SG, x, spectrum, background, peakList, settings);
			// if (m_SG == 1152 || m_SG == 1184 || m_SG == 688 || m_SG == 256 || m_SG == 76) { //C1s_CuO
			// if (m_SG == 528 || m_SG == 2 || m_SG == 2880) {
			//	output_history(m_SG, x, spectrum, background, peakList, settings);
			// }
			if (!outputHistoryList.empty())
			{
				if (exist_in(outputHistoryList, m_SG))
				{
					output_history(m_SG, x, spectrum, background, peakList, settings);
				}
			}
		}
		else
		{
			// Check a slow convergence
			if (abs((oldChiSquare - newChiSquare) / newChiSquare) < convergence)
			{
				if (++countSlowConvergence >= 10)
				{
					cout << "Found a slow convergence. break this loop." << endl;
					break;
				}
			}
			else
			{
				countSlowConvergence = 0;
			}
			// Restore peakList, bgHigher, bgLower and chiSquare.
			// Let lambda large, to achieve a slow convergence.
			lambda *= lambda_ascend;
			copy_peakparameter(peakList_saved, bgHigher_saved, bgLower_saved, peakList, bgHigher, bgLower);
			isPeakListUpdated = false;
		}
	} // Next loop

	// Optimization finished

	optimizedChiSquare = min(oldChiSquare, newChiSquare);
	optimizedReducedChiSquare = min(oldReducedChiSquare, newReducedChiSquare);

	restrictPeaks(peakList, settings);

	int numRemovedPeaks4 = removeInfinitesimalPeaks(peakList);

	// Remove peaks outside the minimum points region.
	int numRemovedPeaks3 = removePeaksByEnergyRange(x, peakList);

	if ((numRemovedPeaks3 > 0) || (numRemovedPeaks4 > 0))
	{
		if (calcChiSquare(x, spectrum, peakList, bgHigher, bgLower, settings, background, qpq, s, newChiSquare, newReducedChiSquare))
		{
			cout << "Error. Failed to calculate newChiSquare." << endl;
			return 1;
		}
		optimizedChiSquare = newChiSquare;
		optimizedReducedChiSquare = newReducedChiSquare;
	}
	cout << "chi2, red.chi2= " << optimizedChiSquare << " " << optimizedReducedChiSquare << endl;

	if (!outputHistoryList.empty())
	{
		if (exist_in(outputHistoryList, m_SG))
		{
			output_history(m_SG, x, spectrum, background, peakList, settings);
		}
	}

	return 0;
}
