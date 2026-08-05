#pragma once

class PeakParameter
{
public:
	double height;		 // Height
	double position;	 // Position (eV)
	double width;		 // Half width at half maximum (HWHM) (eV)
	double lorentzRatio; // mu, ratio of Lorentz function in the pseudo-Voigt function
	PeakParameter(double _h, double _p, double _w, double _mu) : height(_h), position(_p), width(_w), lorentzRatio(_mu) {}
	~PeakParameter() {}
};
