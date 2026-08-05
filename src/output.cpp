#define _USE_MATH_DEFINES
#include "output.h"

inline void denormalize_by_min_and_range2(double &v, const double vmin, const double vrange)
{
	v = vmin + vrange * v;
}

void denormalize_spectra_and_parameters(vector<double> &bindingEnergy, vector<double> &spectrum,
										vector<double> &background, vector<PeakParameter *> &peakList,
										vector<double> &originalBindingEnergy, vector<double> &originalSpectrum,
										vector<double> &smoothedSpectrum, const Settings &settings)
{
	for (int j = 0; j < peakList.size(); j++)
	{
		peakList[j]->height *= settings.spectrum_range;
		denormalize_by_min_and_range2(peakList[j]->position, settings.bindingEnergy_min, settings.bindingEnergy_range);
		peakList[j]->width *= settings.bindingEnergy_range;
	}
	for (int i = 0; i < originalBindingEnergy.size(); i++)
	{
		denormalize_by_min_and_range2(originalBindingEnergy[i], settings.bindingEnergy_min, settings.bindingEnergy_range);
		denormalize_by_min_and_range2(originalSpectrum[i], settings.spectrum_min, settings.spectrum_range);
		denormalize_by_min_and_range2(smoothedSpectrum[i], settings.spectrum_min, settings.spectrum_range);
	}
	for (int i = 0; i < bindingEnergy.size(); i++)
	{
		denormalize_by_min_and_range2(bindingEnergy[i], settings.bindingEnergy_min, settings.bindingEnergy_range);
	}
	for (int i = 0; i < background.size(); i++)
	{
		denormalize_by_min_and_range2(spectrum[i], settings.spectrum_min, settings.spectrum_range);
	}
	for (int i = 0; i < background.size(); i++)
	{
		denormalize_by_min_and_range2(background[i], settings.spectrum_min, settings.spectrum_range);
	}
}

/*
 * ---------------------------------------------------------------------
 * Output a result
 * INPUT
 *    kwd:                Keyword (Ex. 1st, 2nd or final)
 *    bindingEnergy:      Binding energies (eV)
 *    spectrum:           Raw spectrum
 *    background:         Background spectrum
 *    peakList:           Peak parameters
 * ---------------------------------------------------------------------
 */
int output_result(const string &kwd, vector<double> bindingEnergy, vector<double> spectrum,
				  vector<double> background, vector<PeakParameter *> peakList,
				  vector<double> originalBindingEnergy, vector<double> originalSpectrum,
				  const int idHigher, const int idLower, vector<double> smoothedSpectrum,
				  const Settings &settings)
{
	FILE *fp;
	string ofname;
	ofstream ofs;
	const int ndata = (int)bindingEnergy.size();
	const int originalNdata = (int)originalBindingEnergy.size();
	cout << "<<result " + kwd + ">>" << endl;

	// Output a list of peak parameters and their areas
	cout << "Number of peaks=" << peakList.size() << endl;

	if (peakList.size() == 0)
	{
		cout << "Skip to output file because there is no peaks." << endl;
		return 0;
	}

	// shinotsuka20181205a
	//  Before output the result, denormalize spectra and parameters
	denormalize_spectra_and_parameters(bindingEnergy, spectrum, background, peakList,
									   originalBindingEnergy, originalSpectrum, smoothedSpectrum, settings);

	stringstream ss;
	// ss << "peakID,height,position,width,lorentzRatio,area" << endl;

	// for (int j = 0; j < peakList.size(); j++) {
	// 	double area = peakList[j]->height * peakList[j]->width * (peakList[j]->lorentzRatio * M_PI + (1.0 - peakList[j]->lorentzRatio) * sqrt(M_PI / log(2.0)));
	// 	ss << j
	// 		<< "," << peakList[j]->height
	// 		<< "," << peakList[j]->position
	// 		<< "," << peakList[j]->width
	// 		<< "," << peakList[j]->lorentzRatio
	// 		<< "," << area
	// 		<< endl;
	// }

	// cout << ss.str();

	// ofname = "parameters_" + kwd + ".csv";
	// ofs.open(ofname);
	// if (ofs) {
	// 	cout << "Output to " + ofname << endl;
	// 	ofs << ss.str();
	// 	ofs.close();
	// }
	// else {
	// 	cout << "\n!!! Failed to create a file: " << ofname << " !!!\n" << endl;
	// }

	// ss.str("");

	ofname = "parameters_" + kwd + ".csv";
	if ((fp = fopen(ofname.c_str(), "w")) == NULL)
	{
		printf("Error. Cannot open %s\n", ofname.c_str());
		exit(1);
	}
	fprintf(fp, "peakID,height,position,width,lorentzRatio,area\n");
	for (int j = 0; j < peakList.size(); j++)
	{
		double area = peakList[j]->height * peakList[j]->width * (peakList[j]->lorentzRatio * M_PI + (1.0 - peakList[j]->lorentzRatio) * sqrt(M_PI / log(2.0)));
		fprintf(fp, "%d,%.5e,%.3f,%.3f,%.3f,%.5e\n", j, peakList[j]->height, peakList[j]->position, peakList[j]->width, peakList[j]->lorentzRatio, area);
	}
	fclose(fp);

	// Output fitted spectra
	vector<double> sumSpectrum;
	vector<vector<double>> peaks;
	// shinotsuka20180611
	// make_spectrum_and_peak(bindingEnergy, peakList, sumSpectrum, peaks);
	make_spectrum_and_peak(originalBindingEnergy, peakList, sumSpectrum, peaks);

	// shinotsuka20180611
	ss << "# idHigher = " << idHigher << endl;
	ss << "# idLower = " << idLower << endl;

	ss << "Binding_Energy(eV),spectrum,smoothing,background,sum_peaks,fitting";
	for (int j = 0; j < peakList.size(); j++)
	{
		ss << ",peak[" << j << "]";
	}
	ss << endl;

	// shinotsuka20180613
	vector<double> background_extend(originalBindingEnergy.size());
	// for (int i = 0; i < originalBindingEnergy.size(); i++) {
	//	if (i <= idLower) {
	//		background_extend[i] = background[i];
	//	}
	//	else {
	//		background_extend[i] = background_extend[idLower];
	//	}
	// }
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
	}
	for (int i = 0; i < originalBindingEnergy.size(); i++)
	{
		ss << originalBindingEnergy[i]
		   << "," << originalSpectrum[i]
		   << "," << smoothedSpectrum[i]
		   << "," << background_extend[i]
		   << "," << sumSpectrum[i]
		   << "," << background_extend[i] + sumSpectrum[i];
		for (int j = 0; j < peakList.size(); j++)
		{
			ss << "," << background_extend[i] + peaks[j][i];
		}
		ss << endl;
	}

	ofname = "result_" + kwd + ".csv";
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

int outputBackgroundList(const vector<double> &bindingEnergy, const vector<double> &spectrum, const vector<MinimumPoints *> &minimumPointsList, const vector<vector<double>> &backgroundList)
{
	string ofname;
	ofstream ofs;
	string kwd = "backgroundlist";
	const int ndata = (int)bindingEnergy.size();
	int nloop = 0;
	for (int i = 0; i < minimumPointsList.size(); i++)
	{
		if (minimumPointsList[i]->position > bindingEnergy.back())
		{
			nloop++;
		}
	}

	// test
	cout << "minimumPointsList.size=" << minimumPointsList.size() << endl;

	cout << "<<result " + kwd + ">>" << endl;
	stringstream ss;

	ofname = "result_" + kwd + ".csv";
	ofs.open(ofname);

	ss << "Binding_Energy(eV),spectrum";
	for (int j = 0; j < nloop; j++)
	{
		ss << ",BG[" << j << "]";
	}
	ss << endl;
	for (int i = 1; i < ndata; i++)
	{
		ss << bindingEnergy[i] << "," << spectrum[i];
		for (int k = 0; k < nloop; k++)
		{
			ss << "," << backgroundList[k][i];
		}
		ss << endl;
	}

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

void outputBackground(const vector<double> &bindingEnergy, const vector<double> &spectrum, const vector<double> &background, const vector<double> &diffSpectrum)
{
	int ndata = (int)bindingEnergy.size();
	ofstream ofs("out_background.csv");
	ofs << "Binding_Energy(eV),spectrum,background,diff" << endl;

	for (int j = 0; j < ndata; j++)
	{
		ofs << bindingEnergy[j] << "," << spectrum[j] << "," << background[j] << "," << diffSpectrum[j] << endl;
	}
	ofs.close();
}

void outputSavitzkyGolay(const vector<double> &bindingEnergy, const vector<double> &diffSpectrum,
						 const vector<double> &smoothedSpectrum, const vector<double> &thirdDerivative)
{
	ofstream ofs("out_Savitzky_Golay.csv");
	ofs << "Binding_Energy(eV),diff,smoothedSpectrum,thirdDerivative" << endl;
	for (int i = 0; i < bindingEnergy.size(); i++)
	{
		ofs << bindingEnergy[i] << "," << diffSpectrum[i] << "," << smoothedSpectrum[i] << "," << thirdDerivative[i] << endl;
	}
	ofs.close();
}

void outputMaximumPoints(const string &kwd, const vector<double> &bindingEnergy, const vector<double> &spectrum, vector<MaximumPoints *> &maximumPointsList)
{
	string ofname;
	ofstream ofs;
	const int ndata = (int)bindingEnergy.size();
	cout << "<<result " + kwd + ">>" << endl;

	stringstream ss;

	ss << "maximumPointsID,height,position,id" << endl;

	// Output a list of minimum points parameters and their areas
	cout << "Number of maximum points=" << maximumPointsList.size() << endl;

	for (int k = 0; k < maximumPointsList.size(); k++)
	{
		ss << k
		   << "," << maximumPointsList[k]->height
		   << "," << maximumPointsList[k]->position
		   << "," << maximumPointsList[k]->id
		   << endl;
	}

	cout << ss.str();

	ofname = "parameters_" + kwd + ".csv";
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
}

/*
 * ---------------------------------------------------------------------
 * Results of local minimum points
 * INPUT
 *    kwd:                Keyword (Ex. 1st, 2nd or final)
 *    bindingEnergy:      Binding energies (eV)
 *    spectrum:           Raw spectrum
 *    minimumPointsList:  Minimum points parameters
 * ---------------------------------------------------------------------
 */
void outputMinimumPoints(const string &kwd, const vector<double> &bindingEnergy, const vector<double> &spectrum, vector<MinimumPoints *> &minimumPointsList)
{
	string ofname;
	ofstream ofs;
	const int ndata = (int)bindingEnergy.size();
	cout << "<<result " + kwd + ">>" << endl;

	stringstream ss;

	ss << "minimumPointsID,height,position,id" << endl;

	// Output a list of minimum points parameters and their areas
	cout << "Number of minimum points=" << minimumPointsList.size() << endl;

	for (int k = 0; k < minimumPointsList.size(); k++)
	{
		ss << k
		   << "," << minimumPointsList[k]->height
		   << "," << minimumPointsList[k]->position
		   << "," << minimumPointsList[k]->id
		   << endl;
	}

	cout << ss.str();

	ofname = "parameters_" + kwd + ".csv";
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
}
