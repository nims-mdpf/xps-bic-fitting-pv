#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include "peakparameter.h"

using namespace std;

double lorentzFunc(const double x, const double c, const double w);
double gaussFunc(const double x, const double c, const double w);
double voigtFunc(const double x, const PeakParameter *pp);
int make_spectrum_and_peak(const vector<double> &xx, const vector<PeakParameter *> &peakList, vector<double> &sumSpectrum, vector<vector<double>> &peaks);
int make_spectrum(const vector<double> &xx, const vector<PeakParameter *> &peakList, vector<double> &sumSpectrum);
double ran0();
int poissonNoise(double lambda);
double gaussNoise(double sigma, double lambda);
int savitzkyGolayGeneral(const vector<double> &y, const int numSmooth, const int &order, vector<double> &s);
int existFile(const char *path);
double calcRMSE(const vector<double> &v1, const vector<double> &v2, const int errType = 1);

/*
 * ---------------------------------------------------------------------
 * Release memory of a vector variable like vector<T*>
 * T is an arbitrary class, e.g. int, double, MyClass, etc.
 * ---------------------------------------------------------------------
 * // Example
 * vector<MyClass*> vec;
 * vec.push_back(new MyClass(...));
 * ...
 * // When the vector is no longer needed
 * releaseMemory(vec);
 * ---------------------------------------------------------------------
 */
template <typename T>
void releaseMemoryT(vector<T *> &vList)
{
	for (T *v : vList)
	{
		delete v;
	}
	vector<T *>().swap(vList);
}

/*
 * ---------------------------------------------------------------------
 * Make a copy of a vector variable like vector<T*>
 * T is an arbitrary class, e.g. int, double, MyClass, etc.
 * ---------------------------------------------------------------------
 */
template <typename T>
void copyVectorT(const vector<T *> &vin, vector<T *> &vout)
{
	releaseMemoryT(vout);
	for (T *v : vin)
	{
		vout.push_back(new T(*v));
	}
}

/*
 * "in" function like python
 */
template <typename Container>
bool exist_in(const Container &c, const typename Container::value_type &v)
{
	return (c.end() != std::find(c.begin(), c.end(), v));
}
