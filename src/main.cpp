/*
* ---------------------------------------------------------------------
* Automated peak separation program
*  - It produces at most 155 optimized models
*  - a model is composed by a superposition of pseudo-Voigt functions and a Shirley background
*  - initial models are created by changing the smoothing degree, 
*    candidates of peaks are detected using 3rd derivative of a smoothed spectrum
*  - parameters of a model is optimized by the Levenberg-Marquardt method
*  - background is also optimized by the active Shirley method
* 
* H. Shinotsuka 2019/08/19
* ---------------------------------------------------------------------
*/
#define _USE_MATH_DEFINES
#include <cmath>
#include "minimumpointsparameter.h"

#include "main.h"
//#include <omp.h>
#include "mytimer.h"

#include "minimumpoints.h"

using namespace std;

int main(int argc, char *argv[])
{
	show_version();
	MyTimer mytimer;
	//int isDebug = 1; // 0 = normal mode. 1 = debug mode.

	string fname;

	if (argc >= 2) {
		fname = string(argv[1]);
	}
	else {
		show_usage(argv[0]);
		exit(1);
	}

	// Load a setting file
	Settings settings;
	settings.setOptionalArgs(argc, argv);
	settings.outputSettings();

	// Load an input csv file as a raw spectrum.
	vector<double> originalBindingEnergy, originalSpectrum;
	if (readcsv(fname, originalBindingEnergy, originalSpectrum, settings)) {
		cout << "Failed to read " << fname << "." << endl;
		exit(1);
	}

	mainFlow(originalBindingEnergy, originalSpectrum, settings, mytimer);

	//mytimer.show();

	cout << "Main() Finished." << endl;
	//ofstream ofs("as_exe.done"); ofs.close();
	return 0;
}

void show_version()
{
	cout << "#***********************************************************" << endl;
	cout << "#    Active Shirley method" << endl;
	cout << "#    Compiled at " << __TIME__ << " on " << __DATE__ << endl;
	cout << "#***********************************************************" << endl;
}

vector<string> split(string& input, char delimiter)
{
	istringstream stream(input);
	string field;
	vector<string> result;
	while (getline(stream, field, delimiter)) {
		result.push_back(field);
	}
	return result;
}

void show_usage(char *argv0)
{
	cout << "Usage: " << endl << "    a.exe <csv_file>" << endl << endl;
	//printf("Usage: \n    %s <csv_file>\n\n", argv0);
}

void normalize_by_min_and_range2(double &v, const double vmin, const double vrange)
{
	v = (v - vmin) / vrange;
}


void normalize_spectrum(vector<double> &x, vector<double> &y, Settings &settings)
{
	const int ndata = (int)x.size();
	double y_max = *max_element(y.begin(), y.end());
	//double y_min = *min_element(y.begin(), y.end());
	double y_min = 0.0;
	double y_range = y_max - y_min;

	double x_max = *max_element(x.begin(), x.end());
	double x_min = *min_element(x.begin(), x.end());
	double x_range = x_max - x_min;

	cout << "x_min, x_max, x_range = " << x_min << "," << x_max << "," << x_range << endl;
	cout << "y_min, y_max, y_range = " << y_min << "," << y_max << "," << y_range << endl;

	for (int i = 0; i < ndata; i++) {
		//normalize_by_min_and_range2(y[i], y_min, y_range);
		y[i] /= y_max;
		normalize_by_min_and_range2(x[i], x_min, x_range);
	}

	settings.bindingEnergy_min = x_min;
	settings.bindingEnergy_max = x_max;
	settings.bindingEnergy_range = x_range;
	settings.spectrum_min = y_min;
	settings.spectrum_max = y_max;
	settings.spectrum_range = y_range;

}

/*
* ---------------------------------------------------------------------
* Read a csv file. Lines of strings are neglected.
* INPUT
*    fname:              File name
*    settings:           Settings
* OUTPUT
*    x:                  The 1st column
*    y:                  The 2nd column
* ---------------------------------------------------------------------
*/
//int readcsv(const string &fname, vector<double> &x, vector<double> &y, const Settings &settings)
int readcsv(const string &fname, vector<double> &x, vector<double> &y, Settings &settings)
{
	cout << "Reading " << fname << endl;
	cout << "bindingEnergyCutoffs= " << settings.bindingEnergyCutoffs[0] << " " << settings.bindingEnergyCutoffs[1] << endl;
	ifstream ifs;
	ifs.open(fname.c_str());
	if (ifs.fail()) {
		return 1;
	}
	x.clear();
	y.clear();
	string line;

	while (getline(ifs, line)) {
		replace(line.begin(), line.end(), ',', ' ');
		if (line.empty() || line[0] == '#') continue;
		stringstream buf(line);
		double xt = 0.0, yt = 0.0;
		if (buf >> xt >> yt) {
			//cout << "xt,yt= " << xt << " " << yt << endl;
			if (xt >= settings.bindingEnergyCutoffs[0] && xt <= settings.bindingEnergyCutoffs[1]) {
				x.push_back(xt);
				y.push_back(yt);
			}
		}
		else {
			cout << "Error to read a line= " << buf.str() << endl;
		}
	}

	const int originalNdata = (int)y.size(); // Number of data poinsts
	cout << "Number of original data= " << originalNdata << endl;
	if (originalNdata == 0) {
		cout << "Error. Number of original data points is zero." << endl;
		return 1;
	}

	//shinotsuka20181204, normalize spectrum
	normalize_spectrum(x, y, settings);

	return 0;
}

//shinotsuka20180228
// Read the result from out_auto_smoothing2.csv,
// Ignore a row which has a very narrow peak width,
// Sort by the reduced chi square.
void sort_result(const Settings &settings)
{
	cout << endl << "sort_result() start" << endl;
	const int ncol = 9;  // Number of columns

	struct Outauto {
		double data[ncol];
	};

	const string fname = "out_auto_smoothing2.csv";
	vector<Outauto> oaList;

	cout << "Reading " << fname << endl;
	ifstream ifs;
	ifs.open(fname.c_str());
	if (ifs.fail()) {
		//return 1;
		exit(1);
	}

	string line;
	string header;
	getline(ifs, header); // Skip the header

	while (getline(ifs, line)) {
		replace(line.begin(), line.end(), ',', ' ');
		if (line.empty() || line[0] == '#') continue;
		bool ok = true;
		Outauto oa;
		stringstream buf(line);
		for (int j = 0; j<ncol; j++) {
			if (!(buf >> oa.data[j])) {
				cout << "Error to read a line= " << buf.str() << endl;
				//exit(1);
				ok = false;
				break;
			}
		}
		if (!ok) continue;
		oaList.push_back(oa);
	}
	ifs.close();

	//cout << "Number of data before= " << oaList.size() << endl;

	// Erase a row where a minimum HWHM <= settings.minimumPeakFWHM / 2
	vector<Outauto>::iterator it = oaList.begin();
	while (it != oaList.end()) {
		if ((it->data[7]) * 2.0 <= settings.minimumPeakFWHM) {
			it = oaList.erase(it);
		}
		else {
			++it;
		}
	}
	cout << "Number of data= " << oaList.size() << endl;

	// data[4]すなわちreduced chi squareをキーにして昇順ソート
	sort(oaList.begin(), oaList.end(),
		[](const Outauto& x, const Outauto& y) { return x.data[4] < y.data[4]; });

	// Output
	const string fname2 = "out_auto_smoothing3.csv";
	cout << "Output to " << fname2 << endl;
	ofstream ofs;
	ofs.open(fname2.c_str());
	if (ofs.fail()) {
		exit(1);
	}

	ofs << header << endl;
	for (int i = 0; i < oaList.size(); i++) {
		ofs << oaList[i].data[0];
		for (int j = 1; j<ncol; j++) {
			ofs << "," << oaList[i].data[j];
		}
		ofs << endl;
	}
	ofs.close();

	cout << "sort_result() end" << endl;
}


ofstream fs("as_exe.done");
//fs.close();
