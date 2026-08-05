#pragma once

class MinimumPoints
{
public:
	double height;	 // Height
	double position; // position (eV)
	int id;			 // ID that corresponds to the position
	MinimumPoints(double _h, double _p, int _id) : height(_h), position(_p), id(_id) {}
	~MinimumPoints() {}
};

class MaximumPoints
{
public:
	double height;	 // Height
	double position; // position (eV)
	int id;			 // ID that corresponds to the position
			// MaximumPoints(double _h, double _p, double _id) :
			//	height(_h), position(_p), id(_id) {}
	MaximumPoints(double _h, double _p, int _id) : height(_h), position(_p), id(_id) {}
	~MaximumPoints() {}
};