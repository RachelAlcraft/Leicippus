#include "Xray.h"

#define M_PI 3.14159265358979323846  /* pi */

Xray::Xray(double amplitude, double phase, double lamda, VectorThree direction)
{
	_amplitude = amplitude;
	_phase = phase;
	_direction = direction;
	Lamda = lamda;
}

double Xray::displacement(double x)
{
	double y = _amplitude * cos((2 * M_PI * x) - _phase);
	return y;
}

Xray Xray::scatter(VectorThree position, double theta, double electron_density)
{
	double newAmplitude = _amplitude * electron_density; //or some other function
	//also dependent on theta and lamda

	return *this;
}

/*
* The imaginary number is:
* Since cos(phi) + i*sin(phi) = e^(i*phi)
Z = _amplitude*e^(i*_phase)      = _amplitude(cos(_phase) + i*sin(_phase))
https://onlinelibrary.wiley.com/doi/pdf/10.1002/9783527619139.app1
See A3, A is the complex amplitude, I am not happy
*/