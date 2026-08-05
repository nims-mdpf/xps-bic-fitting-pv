#pragma once
#include <cstdio>
#include <string>
#include <chrono>
using namespace std::chrono;

class MyTimer
{
private:
	time_point<system_clock> start_time, end_time;
	long long elapsed_time, split_time;
	void calc_elapsed_time()
	{
		time_point<system_clock> previous_time = end_time;
		end_time = system_clock::now();
		elapsed_time = duration_cast<microseconds>(end_time - start_time).count();
		split_time = duration_cast<microseconds>(end_time - previous_time).count();
	}

public:
	MyTimer()
	{
		reset();
	}
	~MyTimer() {}
	void reset()
	{
		start_time = system_clock::now();
		end_time = start_time;
	}
	long long get_elapsed_time()
	{
		// calc_elapsed_time();
		return elapsed_time;
	}
	string time_to_str(long long residue)
	{
		// calc_elapsed_time();
		// long long residue = elapsed_time;
		long long cnst[5] = {3600000000, 60000000, 1000000, 1000, 1};
		int hmmss[5] = {0};
		for (int i = 0; i < 5; i++)
		{
			hmmss[i] = (int)(residue / cnst[i]);
			residue = residue % cnst[i];
		}
		char chr[100];
#ifdef _MSC_VER
		sprintf_s(chr, 100, "%d:%02d:%02d.%03d%03d", hmmss[0], hmmss[1], hmmss[2], hmmss[3], hmmss[4]);
#else
		sprintf(chr, "%d:%02d:%02d.%03d%03d", hmmss[0], hmmss[1], hmmss[2], hmmss[3], hmmss[4]);
#endif
		return string(chr);
	}
	//	string get_elapsed_time_str()
	//	{
	//		//calc_elapsed_time();
	//		long long residue = elapsed_time;
	//		long long cnst[5] = { 3600000000, 60000000, 1000000, 1000, 1 };
	//		int hmmss[5] = { 0 };
	//		for (int i = 0; i < 5; i++) {
	//			hmmss[i] = (int)(residue / cnst[i]);
	//			residue = residue % cnst[i];
	//		}
	//		char chr[100];
	// #ifdef _MSC_VER
	//		sprintf_s(chr, 100, "%d:%02d:%02d.%03d%03d", hmmss[0], hmmss[1], hmmss[2], hmmss[3], hmmss[4]);
	// #else
	//		sprintf(chr, "%d:%02d:%02d.%03d%03d", hmmss[0], hmmss[1], hmmss[2], hmmss[3], hmmss[4]);
	// #endif
	//		return string(chr);
	//	}
	void show(string kwd = "")
	{
		calc_elapsed_time();
		cout << "Elapsed time(string) (" << kwd << ") = " << time_to_str(elapsed_time) << ", Split time = " << time_to_str(split_time) << endl;
		cout << "Elapsed time(microseconds) (" << kwd << ") = " << elapsed_time << ", Split time = " << split_time << endl;
		// cout << "Elapsed time(microseconds)= " << get_elapsed_time() << endl;
		// cout << "Elapsed time= " << get_elapsed_time_str() << endl;
	}
};
