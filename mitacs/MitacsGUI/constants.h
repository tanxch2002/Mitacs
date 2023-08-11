// physical and other constants - global variables
//#ifndef CONSTANTS_H
//#define CONSTANTS_H

//#include <cmath>

extern const double Boltzmann_constant = 0.69503; 
// this is actually the conversion factor , 1 K = 0.69 cm-1, NOT a Boltzmann constant
//New 2021: So energy(cm-1) with Boltzman_constant(cm-1)
extern const double pi = 3.141592;
extern const double round_off = 0.0000001;
// as simulation time and step are defined as double, the round off issue makes problem and the last step doesn't work.
extern const double round_off_limit = 1e-8; 
// I may use this limit to set some very small numbers to zero.
extern const double accuracy_measure = 1e-1; 
// too large
extern const double extreme_small_number = 1e-70;

extern const double energy_conversion_factor = 1e3;
// to convert energy (frequency) from GHz to MHz (by multipilication).
// All equation have been drived by assumption that energy is in GHz and time in nano second. 
//However, later Mehdi decided to change the scales to MHz and microseconds.
extern const double time_conversion_factor = 1e-3;
// to convert time from nano second to micro second (by multipilication).
extern const double light_speed = 2.99792458 * 1e8;
// m/s, the speed of light in vacuum
extern const double h_planck = 6.626068 * 1e-34;
// J.S.  planck's constant
extern const double k_b = 1.38e-23;
extern const int max_ensemble_for_writing = 100; 
// Some files were designed to check the validity
//const double w_0 = 1e3 * energy_conversion_factor; //MHz. i.e. 10^12 Hz.

 

//#endif // CONSTANTS_H
