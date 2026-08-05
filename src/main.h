#pragma once

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "settings.h"
#include "mainflow.h"

using namespace std;

int main(int argc, char *argv[]);
void show_version();
vector<string> split(string& input, char delimiter);
void show_usage(char *argv0);
void normalize_by_min_and_range2(double &v, const double vmin, const double vrange);
void normalize_spectrum(vector<double> &x, vector<double> &y, Settings &settings);
int readcsv(const string &fname, vector<double> &x, vector<double> &y, Settings &settings);
void sort_result(const Settings &settings);
