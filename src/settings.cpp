#include "settings.h"

Settings::Settings()
{
	//if (readSettings()) {
	//	setDefaultSettings();
	//}
	setDefaultSettings();
	readSettings();
	//outputSettings();
}

int Settings::readSettings()
{
	const string settingfile = "settings.inp";
	ifstream ifs;
	ifs.open(settingfile.c_str());
	if (ifs.fail()) {
		return 1;
	}

	cout << "Read settings from " << settingfile << endl;

	string line;
	string str1, str2;
	while (getline(ifs, line)) {
		replace(line.begin(), line.end(), ',', ' ');
		if (line.empty() || line[0] == '#') continue;
		stringstream buf(line);
		if (!(buf >> str1)) {
			cout << "  Skip. " << buf.str() << endl;
			continue;
		}

		bool ok = false;
		//if (str1 == "SavitzkyGolaySmoothingPoint") {
		//	ok = (bool)(buf >> nSavitzkyGolaySmoothingPoint);
		//}
		if (str1 == "ratioPeakHeightThreshold") {
			ok = (bool)(buf >> ratioPeakHeightThreshold);
		}
		else if (str1 == "bindingEnergyCutoffs") {
			ok = (bool)(buf >> bindingEnergyCutoffs[0] >> bindingEnergyCutoffs[1]);
		}
		else if (str1 == "maxLoopIterativeShirley") {
			ok = (bool)(buf >> maxLoopIterativeShirley);
		}
		else if (str1 == "allowNegativeIterativeShirley") {
			ok = (bool)(buf >> allowNegativeIterativeShirley);
		}
		//else if (str1 == "allowHigherEdgeActive") {
		//	ok = (bool)(buf >> allowHigherEdgeActive);
		//}
		//else if (str1 == "allowLowerEdgeActive") {
		//	ok = (bool)(buf >> allowLowerEdgeActive);
		//}
		else if (str1 == "lorentzRatioCutoffs") {
			ok = (bool)(buf >> lorentzRatioCutoffs[0] >> lorentzRatioCutoffs[1]);
		}
		else if (str1 == "minimumPeakFWHM") {
			ok = (bool)(buf >> minimumPeakFWHM);
		}
		else if (str1 == "minimumPeakAreaRatio") {
			ok = (bool)(buf >> minimumPeakAreaRatio);
		}
		else if (str1 == "allowedNegativeRatio") {
			ok = (bool)(buf >> allowedNegativeRatio);
		}
		else if (str1 == "backgroundType") {
			ok = (bool)(buf >> backgroundType);
		}
		else {
			cout << "  Unknown key: " << str1 << endl;
			continue;
		}

		if (!ok) {
			cout << "  Skip. " << buf.str() << endl;
			continue;
		}
	} // Next line

	return 0;
}

void Settings::setDefaultSettings()
{
	//nSavitzkyGolaySmoothingPoint = 0;
	ratioPeakHeightThreshold = 0.01;
	bindingEnergyCutoffs[0] = -1e300;
	bindingEnergyCutoffs[1] = 1e300;
	maxLoopIterativeShirley = 15;
	allowNegativeIterativeShirley = 1;
	//allowHigherEdgeActive = 1;
	//allowLowerEdgeActive = 1;
	lorentzRatioCutoffs[0] = 0.0;
	lorentzRatioCutoffs[1] = 1.0;
	minimumPeakFWHM = 0.2;
	minimumPeakAreaRatio = 0.01;  // or 0.005
	allowedNegativeRatio = 0.005;  // or 0.01
	backgroundType = "Shirley";
}


/*
* ---------------------------------------------------------------------
* Analyze optional command line arguments
* argv[1] is reserved as a file name of input spectrum data
* argv[2], argv[3], ... are analyzed here
*
*	-s    Use Shirley background (Default)
*   -l    Use Linear background
*   
* ---------------------------------------------------------------------
*/
#define ARG_FLAG(s) (strcmp(s,argv[i])==0)
#define ARG_SPARAM(s) ((strcmp(s,argv[i])==0)&&(i+1)<argc&&((r=argv[++i])!=NULL))
#define ARG_LPARAM(s) ((strncmp(s,argv[i],strlen(s))==0)&&((l=strlen(s))<strlen(argv[i]))&&((r=&argv[i][l])!=NULL))
void Settings::setOptionalArgs(int argc, char* argv[])
{
	cout << "Settings::setOptionalArgs" << endl;
	cout << "Number of command line arguments = " << argc << endl;

	if (argc <= 2) {
		return;
	}

	for (int i = 2; i<argc; i++) {
		//cout << "Check " << argv[i] << endl;
		char*r = NULL;
		int l = 0;
		if (argv[i][0] != '-') {
			cout << "Error. Unrecognized option " << argv[i] << endl;
			exit(1);
		}
		else if (ARG_FLAG("-l")) {
			printf("Option l: backgroundType = 'Linear'\n");
			backgroundType = "Linear";
		}
		else if (ARG_FLAG("-s")) {
			printf("Option s: backgroundType = 'Shirley'\n");
			backgroundType = "Shirley";
		}
		//else if (ARG_SPARAM("-s")) {
		//	printf("option s with value %s\n", r);
		//}
		else if (ARG_LPARAM("--allowedNegativeRatio=")) {
			printf("Option allowedNegativeRatio with value %s\n", r);
			//sscanf(r, "%lf", &allowedNegativeRatio);
			stringstream buf(r);
			if (!(bool)(buf >> allowedNegativeRatio)) {
				cout << "Error. Unrecognized allowedNegativeRatio value " << r << endl;
				exit(1);
			}
			cout << "allowedNegativeRatio = " << allowedNegativeRatio << endl;
		}
		else {
			cout << "Error. Unrecognized option " << argv[i] << endl;
			exit(1);
		}
	}

	return;
}

void Settings::outputSettings()
{
	cout << "Settings:" << endl;
	//cout << "  nSavitzkyGolaySmoothingPoint= " << nSavitzkyGolaySmoothingPoint << endl;
	cout << "  ratioPeakHeightThreshold= " << ratioPeakHeightThreshold << endl;
	cout << "  bindingEnergyCutoffs= " << bindingEnergyCutoffs[0] << " " << bindingEnergyCutoffs[1] << endl;
	cout << "  maxLoopIterativeShirley= " << maxLoopIterativeShirley << endl;
	cout << "  allowNegativeIterativeShirley= " << allowNegativeIterativeShirley << endl;
	//cout << "  allowHigherEdgeActive= " << allowHigherEdgeActive << endl;
	//cout << "  allowLowerEdgeActive= " << allowLowerEdgeActive << endl;
	cout << "  lorentzRatioCutoffs= " << lorentzRatioCutoffs[0] << " " << lorentzRatioCutoffs[1] << endl;
	cout << "  minimumPeakFWHM= " << minimumPeakFWHM << endl;
	cout << "  minimumPeakAreaRatio= " << minimumPeakAreaRatio << endl;
	cout << "  allowedNegativeRatio= " << allowedNegativeRatio << endl;
	cout << "  backgroundType= " << backgroundType << endl;
	cout << endl;
}
