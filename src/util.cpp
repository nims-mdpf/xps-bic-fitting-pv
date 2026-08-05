#define _CRT_SECURE_NO_WARNINGS
#include "util.h"

using namespace std;
/*
 * ---------------------------------------------------------------------
 * Lorentz function L(x; c, w), where the maximum height = 1.0 at x = c.
 * INPUT
 *    x:  Main variable
 *    c:  Peak position
 *    w:  Half width at half maximum (HWHM)
 * ---------------------------------------------------------------------
 */
double lorentzFunc(const double x, const double c, const double w)
{
	return 1.0 / (1.0 + pow((x - c) / w, 2));
}

/*
 * ---------------------------------------------------------------------
 * Gaussian function G(x; c, w), where the maximum height = 1.0 at x = c.
 * INPUT
 *    x:  Main variable
 *    c:  Peak position
 *    w:  Half width at half maximum (HWHM)
 * ---------------------------------------------------------------------
 */
double gaussFunc(const double x, const double c, const double w)
{
	return pow(2.0, -pow((x - c) / w, 2));
}

/*
 * ---------------------------------------------------------------------
 * Pseudo-Voigt function V(x), where the maximum height = pp->height
 * at x = pp->position.
 * INPUT
 *    x:  Main variable
 *    pp: Peak parameter
 * ---------------------------------------------------------------------
 */
double voigtFunc(const double x, const PeakParameter *pp)
{
	double voigt = 0.0;
	voigt = pp->height * (pp->lorentzRatio * lorentzFunc(x, pp->position, pp->width) + (1.0 - pp->lorentzRatio) * gaussFunc(x, pp->position, pp->width));
	return voigt;
}

int make_spectrum_and_peak(const vector<double> &xx, const vector<PeakParameter *> &peakList, vector<double> &sumSpectrum, vector<vector<double>> &peaks)
{
	const int ndata = (int)xx.size();
	if (ndata == 0)
	{
		return 1;
	}
	sumSpectrum.clear();
	sumSpectrum.resize(ndata);
	peaks.clear();
	peaks.resize(peakList.size());
	for (int j = 0; j < peakList.size(); j++)
	{
		peaks[j].resize(ndata);
	}
	for (int i = 0; i < ndata; i++)
	{
		double x = xx[i];
		sumSpectrum[i] = 0.0;
		for (int j = 0; j < peakList.size(); j++)
		{
			peaks[j][i] = voigtFunc(x, peakList[j]);
			sumSpectrum[i] += peaks[j][i];
		}
	}
	return 0;
}
int make_spectrum(const vector<double> &xx, const vector<PeakParameter *> &peakList, vector<double> &sumSpectrum)
{
	const int ndata = (int)xx.size();
	if (ndata == 0)
	{
		return 1;
	}
	sumSpectrum.resize(ndata);
	for (int i = 0; i < ndata; i++)
	{
		double x = xx[i];
		sumSpectrum[i] = 0.0;
		for (int j = 0; j < peakList.size(); j++)
		{
			sumSpectrum[i] += voigtFunc(x, peakList[j]);
		}
	}
	return 0;
}

/*
 * ---------------------------------------------------------------------
 * Uniform random number 0.0 < ran0() < 1.0
 * ---------------------------------------------------------------------
 */
double ran0()
{
	// return ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
	static mt19937 mt(12345);
	return ((double)mt() + 1.0) / ((double)mt.max() + 2.0);
}

/*
 * ---------------------------------------------------------------------
 * Noised intensity that conforms to the Poisson distribution, where
 * the mean and the covariance are both lambda.
 * The returned value is 0, 1, 2, 3, ... .
 * ---------------------------------------------------------------------
 */
int poissonNoise(double lambda)
{
	int k = 0;
	double ln_v = lambda + log(ran0());
	while (ln_v >= 0.0)
	{
		ln_v += log(ran0());
		k++;
	}
	return k;
}

/*
 * ---------------------------------------------------------------------
 * Noised intensity that conforms to the Gauss distribution, where
 * the mean is lambda and the standard deviation is sigma.
 * The returned value is 0, 1, 2, 3, ... .
 * ---------------------------------------------------------------------
 */
double gaussNoise(double sigma, double lambda)
{
	double r1 = ran0();
	double r2 = ran0();
	double v = sqrt(-2.0 * log(r1)) * cos(2.0 * M_PI * r2);
	// v = v * sqrt(sigma) + lambda;
	v = v * sigma + lambda;
	return v;
}

/*
 * ---------------------------------------------------------------------
 * The Savitzky-Golay method to make a smoothed spectrum or 1st, 2nd, 3rd derivative.
 * Only for order = 0, the absolute intensities have meaning. Others are relative intensities.
 * INPUT
 *    ndata:              Number of data points
 *    y:                  Input spectrum
 *    numSmooth:          Smoothing data points around a target point.
 *                        It must be {3, 5, 7, 9, ...} .
 *    order:              Order of the derivative {0, 1, 2, 3}
 *                        0 ==> Smoothing, polynomial degree 2,3
 *                        1 ==> 1st derivative, polynomial degree 1,2
 *                        2 ==> 2nd derivative, polynomial degree 2,3
 *                        3 ==> 3rd derivative, polynomial degree 3,4
 * OUTPUT
 *    s:                  Smoothed spectrum of the (order)-th derivative
 * ---------------------------------------------------------------------
 */
int savitzkyGolayGeneral(const vector<double> &y, const int numSmooth, const int &order, vector<double> &s)
{
	if (numSmooth % 2 == 0 || numSmooth < 3)
	{
		cout << "Error at savitzkyGolayGeneral(). Smoothing data points must be 3, 5, 7, ...." << endl;
		return 1;
	}
	const int ndata = (int)y.size();
	s.resize(ndata);
	if (ndata == 0)
	{
		cout << "Error at savitzkyGolayGeneral(). Number of data = 0." << endl;
		return 1;
	}

	const int m = (numSmooth - 1) / 2;
	// cout << "savitzkyGolay(), m=" << m << endl;

	// Smoothing coefficients
	vector<double> wt(m + 1); // Weight factor for smoothing
	if (order == 0)
	{
		double normFactor = (4.0 * m * m - 1.0) * (2.0 * m + 3.0) / 3.0;
		double wtConst = 3.0 * m * (m + 1.0);
		for (int j = 0; j <= m; j++)
		{
			double temp = (wtConst - 1.0 - 5.0 * j * j);
			wt[j] = temp / normFactor;
		}
	}
	else if (order == 1)
	{
		for (int j = 0; j <= m; j++)
		{
			wt[j] = (double)j;
		}
	}
	else if (order == 2)
	{
		double wtConst = -(double)m * (m + 1.0);
		for (int j = 0; j <= m; j++)
		{
			wt[j] = (1.0 + 2.0 * m) * (wtConst + 3.0 * j * j);
		}
	}
	else if (order == 3)
	{
		double wtConst = 3.0 * m * (m + 1.0);
		for (int j = 0; j <= m; j++)
		{
			double temp = (wtConst - 1.0 - 5.0 * j * j);
			wt[j] = -(double)j * temp;
		}
	}
	else
	{
		cout << "Error at savitzkyGolayGeneral(). Order of the derivative should be {0, 1, 2, 3}." << endl;
		return 1;
	}

	int j, i_plus_j, i_minus_j;
	double sum1;
	double orderSign = 1.0;
	if (order % 2 == 1)
	{
		orderSign = -1.0;
	}
	for (int i = 0; i < ndata; i++)
	{
		sum1 = wt[0] * y[i];
		for (j = 1; j <= m; j++)
		{
			i_plus_j = min(ndata - 1, i + j);
			i_minus_j = max(0, i - j);
			sum1 += wt[j] * (y[i_plus_j] + orderSign * y[i_minus_j]);
		}
		s[i] = sum1;
	}

	// double normFactor = (4.0 * m * m - 1.0) * (2.0 * m + 3.0) / 3.0;
	// double wtFactor = 3.0 * m * (m + 1.0);
	// vector<double> wt0(m + 1); // Weight factor for smoothing
	// vector<double> wt3(m + 1); // Weight factor for smoothed 3rd derivative
	// for (int j = 0; j <= m; j++) {
	//	double temp = (wtFactor - 1.0 - 5.0*j*j);
	//	wt0[j] = temp / normFactor;
	//	wt3[j] = -(double)j * temp;
	//	//cout << "wt0[" << j << "]=" << wt0[j] << "  wt3[" << j << "]=" << wt3[j] << endl;
	// }
	// int j, i_plus_j, i_minus_j;
	// double sum1, sum2;
	////#pragma omp parallel for private(j,i_plus_j, i_minus_j,sum1,sum2)
	// for (int i = 0; i < ndata; i++) {
	//	sum1 = wt0[0] * y[i];
	//	sum2 = 0.0;
	//	//#pragma omp parallel for reduction (+:sum1,sum2)
	//	for (j = 1; j <= m; j++) {
	//		i_plus_j = min(ndata - 1, i + j);
	//		i_minus_j = max(0, i - j);
	//		sum1 += wt0[j] * (y[i_plus_j] + y[i_minus_j]);
	//		sum2 += wt3[j] * (y[i_plus_j] - y[i_minus_j]);
	//	}
	//	s[i] = sum1;
	//	thirdDerivative[i] = sum2;
	// }

	return 0;
}

/*
 * ---------------------------------------------------------------------
 * Check if a file exists or not
 * INPUT
 *    path: File path
 * RETURN:
 *    1    if a file exists
 *    0    if a file does not exist
 * ---------------------------------------------------------------------
 */
int existFile(const char *path)
{
	FILE *fp = fopen(path, "r");
	if (fp == NULL)
	{
		return 0;
	}

	fclose(fp);
	return 1;
}

/*
 * ---------------------------------------------------------------------
 * RMSE or its analogy
 * RMSE (Root mean square error) between two vectors v1[] and v2[]
 * NRMSE (Normalized RMSE): RMSE normalized by the difference between maximum and minimum
 * RMSPEiRoot Mean Square Percentage Error): RMSE normalized by v1[]
 *
 * INPUT
 *    vector<double> v1: original spectrum data
 *    vector<double> v2: smoothed spectrum data
 * RETURN:
 *    RMSE value
 * ---------------------------------------------------------------------
 */
double calcRMSE(const vector<double> &v1, const vector<double> &v2, const int errType)
{
	const int n1 = (int)v1.size();
	const int n2 = (int)v2.size();
	if (n1 != n2)
	{
		cout << "Error at calcRMSE(). Length mismatch " << n1 << " " << n2 << endl;
		exit(1);
	}
	else if (n1 == 0)
	{
		cout << "Error at calcRMSE(). Length of vector is zero " << endl;
		exit(1);
	}

	double s = 0.0;
	if (errType == 0 || errType == 1)
	{
		// RMSE (Root mean square error)
		for (int i = 0; i < n1; i++)
		{
			s += pow(v1[i] - v2[i], 2);
		}
		s = sqrt(s / n1);

		if (errType == 1)
		{
			// NRMSE (Normalized RMSE)
			double minmaxrange = (*max_element(v1.begin(), v1.end())) - (*min_element(v1.begin(), v1.end()));
			s = s / minmaxrange;
		}
	}
	else
	{
		cout << "Error at calcRMSE(). Invalid errType " << errType << endl;
		exit(1);
	}

	return s;
}