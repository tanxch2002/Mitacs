// physical and other constants
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

// 2021: Added some constants.

const double Boltzmann_constant = 0.69503; // this is actually the conversion factor , 1 K = 0.69 cm-1, NOT a Boltzmann constant
//New 2021: So energy(cm-1) with Boltzman_constant(cm-1)
const double pi = 3.141592;
const double round_off = 0.0000001;// as simulation time and step are defined as double, the round off issue makes problem and the last step doesn't work.
const double round_off_limit = 1e-8; // I may use this limit to set some very small numbers to zero.
const double accuracy_measure = 1e-1; // too large
const double extreme_small_number = 1e-70;
int balancing;// balancing in solving eigensystem is on (1) or off (0).
const double energy_conversion_factor = 1e3;
// to convert energy (frequency) from GHz to MHz (by multipilication).
// All equation have been drived by assumption that energy is in GHz and time in nano second. 
//However, later Mehdi decided to change the scales to MHz and microseconds.
const double time_conversion_factor = 1e-3;
// to convert time from nano second to micro second (by multipilication).
const double light_speed = 2.99792458 * 1e8;// m/s, the speed of light in vacuum
const double h_planck = 6.626068 * 1e-34;// J.S.  planck's constant
const double k_b = 20.8366455 * energy_conversion_factor;
//in MHz over Kelvin. k_b(J/K) / h (J.s) = 1.3806503*1e-23 / 6.626068*1e-34 => 
//k_b = 2.08366455 * 1e10  in Hz/K or 2.08366455 * 1e1 GHz/K

//const double w_0 = 1e3 * energy_conversion_factor; //MHz. i.e. 10^12 Hz.

double mass; //in kg
// decreased from 100 amu to 10 amu (look at my log book volume II, Feb 8, 2012) and one amu is 1.66* 10^-27 kg

double d_distance_of_two_wells; //In this version this is the "primary" quantity (one accessible to user via interface) and other quantities depend on this one
//one angstroem or 0.1 nm.  Thickness of the barrier and width of the well in the model with rectangular barriers.
double width_well; //the width of well or thickness of the barrier

int ensemble_number; // number of ensembles used. More ensembles are tested until we find 5000 with resonant wells.
const int max_ensemble_for_writting=100; // Some files were designed to check the validity
int Hamilton_dimension; // Should be 512/1024/2048 for optimal matrix solving - the size of the high resolution energy landscape used in Miller-Colbert algorithm mentioned in Garashchuk paper
										   //

int well_numbers;// let's just make it straightforward...
int point_numbers;
/* the points which are needed for plotting the main curve (parabola). The distance between two adjacent points is 2*interval/(point_numbers-1)
                               (or (max-min)/(point_numbers-1) if the interval is not symmetrical with recpect to zero(. For example we
                               want to have 5 points for drawing parabola so each step should be one unit ([2*2]/[5-1]) to span -2 to 2 (our interval in this case) is five steps.*/
double interval; // as of July 22, 2021 this is determined by the well thickness
//from center of the landscape to edge, in meters. Essentially it is a pseudo-shift of the parabola, allowing for loop 

double step; // in general, (max-min)/(point_numbers);


//coupling between the ground and excited state.
int bottom_decoupling; 
//if this value is 1 then a new set of random variable is used for the bottom of wells in the excited state; otherwise, the same value from the ground state is used (by multipling it by stretch)
int barrier_decoupling; 
//if this value is 1 then a new set of random variable is used for the the barrier heights in the excited state; otherwise, the same value from the ground state is used (by multipling it by stretch)
//
// double burn_wavelength = 686.8131; //nm  or 14560 cm-1 or 4.368 x10^14 Hz
//
//
//parameters for distribution of energy gap between the ground and excited states baselines
double inhomogeneous_bandwidth = 180;// cm^-1.  need conversion to MHz
double standard_deviation_of_translation ; // in MHz
// the relationship between bandwidth in cm^-1 and standard deviation is FWHM = 2.355 * sigma. 
//It should be noted that 1 cm^-1 = 29.9792458 GHz then FWHM=180*29.97 = 5394.6 GHz.
double mean_value_of_translation; // in cm-1, need conversion to MHz
//based on the last paper burn frequency was 14560 wavenumbers or 686.8131 nm 
//438794.499 and 5995.84916 GHz are corresponding to 14641 and 180 wavenumbers from the SDF peak and width of CP43, 
//A-state, respectively, in the last paper. 180 comes from measuring inhomogeneous line broadening experiments.
double burn_frequency;//frequency in cm-1, needs conversion to MHz.

double burn_power;//mW. 1microW is 1e-3, as an example. 5, 10, and 15 W are taken from Kohler's paper.
double beam_diameter ;//1/(2/sqrt(pi)); //cm. NB: In Kohler's paper the power is given in W/cm^2, so I have to conver the diameter in such a way that it gets rid of other constants in my equation (done on on September 12, 2012, in the second log book). In that calculation, area = pi r^2 = pi d^2/4. by assuming d = 2/sqrt(pi) the area would be 1cm^2. BTW sqrt(pi)/2 = 0.886 cm.
double absorption_cross_section_at_peak_RT = 2.258e-16;//cm^2// September 12, 2012 (second log book)
double homogeneous_width_at_RT;//cm^-1 for Cp43
double homogeneous_width_at_BT;//cm^-1 for Cp43
double phonon_coupling ;// this is S or Huang–Rhys factor for Cp43.
double absorption_cross_section_at_peak_BT;//=1.054e-12;//cm^2, see the notes on September 12, 2012 (in the second log book).
double photon_flux ; //burn frequency in MHz by now, Planck constant in SI units
double time_interval_between_two_consecutive_burns ;
//1*(beam_diameter*beam_diameter/(absorption_cross_section_at_peak_BT*burn_power*burn_wavelength))*1.56015e-13;
//in seconds, see the notes on September 12, 2012 (in the second log book).
double scan_time ;
//second. For SMS, 2.6 second comes from the scan speed of 44cm^-1/s (Kohler's paper). 
//Our range here is approximately between 434733-438215GHz (or 14505-14621 wavenumbers) so delta_t = 116/44 = 2.6 the time to scan the whole range. I assumed that laser stays there for 2.6/well_numbers as an approximate value.
double whole_burning_time ;
//scan_time/well_numbers ; // in second, for burning experiment it was 19min = 19 *60 seconds in one of our experiments.

//scan parameters:
int number_of_scan; //one if we do not want to scan, only want the holeburn?
double minimum_frequency;
// 380000 * energy_conversion_factor; // MHz or 14505.5 cm^-1 or 689.4nm. 
//This the minimum frequency that I found it postburn_spectrum.txt in a typlical hole burning calculation for flat landscape with 3000 systems.
double maximum_frequency;
//460000 * energy_conversion_factor; // MHz or 14621.7 cm^-1 or 683.9nm. 
//This the maximum frequency that I found it postburn_spectrum.txt in a typlical hole burning calculation for flat landscape with 3000 systems.
double frequency_step;
// in MHz. I played with this value such that all wells are chosen once only.
int number_of_frequency_steps ;
double each_scan_time;// micro second. 
//It means that 2.6 seconds takes to scan the whole band from minimum_frequency to maximum_frequency 
//(see the note for whole_burning_time, or my notes in the second log book on April 9th, 2013).
double time_interval_between_each_frequency_step;
//
double scaling_factor ;// this factor scales the energy landscape (barrier heights and bottom of wells) in the ground state 
//which in turn affects the barrier heights in the excited state.
double scaling_parameters ;// this is the same as scaling_factor but I want to use this parameter 
//to change mu, sigma, mean, and deviation at the same time. 
//The former factor is reserved for scaling an energy landscape when it is imported in the program.
//
//
// Ground state barrier distribution:
double mu; // this is currently 37,100,000 MHz
// or 37100 GHz or 1233 cm-1
// mu should be around 1000 wavenumbers (look at my log book volume II). 
//This corresponds to the mean value of lambda close to 22 in the ground state.
double sigma;
//NOTICE: Following value is used in most of simulatios: 7990.0 * scaling_parameters * energy_conversion_factor;//mu/4.5;
// I changed this parameter to find lambda distribution in the ground state close to 2 based on the last paper. 
//When denominator is less than 5, 
//random generator makes negative values for potential which is consequently become zero in lambda. 
//This makes non gaussian distribution for lambda.
double stretch; // ratio of excited and ground state BARRIERS
//The stretch should never exceed 1. This resuls in the mean value of lambda close to 11 (10.94) and sigma almost 1.1 in the excited state. 
//Changing this value affects the average energy difference between the bottom of two wells in the ground state 
// i.e. "asymmetry". In the ground state it should be ~10 wavenumbers. see also Curvature Calculation.nb
//

//
// Bottoms of the wells in the ground state
double mean;// the mean value for the bottom of wells; 600 GHz is 20 cm-1
double deviation;
// the standard deviation value for the bottom of wells
//
int generate_energy_landscape ;
//1 means that energy landscape is generated and 0 means that energy landscape (and Translation in case that bottom of wells are independent in two states) is/are imported.
double room_temperature ;//Kelvin.// it could be the same as initial_temperature.

//parameters related to the cooling for finding the distribution before the start of burning.
double extra_temperature; 
// I would like to add extra temperature to measure HGK at different temperatures while not changing the temperature step. 
//It should be noted that the distribution at higher temperatures is close to each other.
double initial_temperature; //Kelvin.
double final_temperature; // Kelvin.
double temperature_step; // the step of reducing of temperature. NB. Changing of this parameter should not be done unless cooling time changes properly.
double cooling_time;  // 2 minutes, time of one cooling step
// In case of burning at different temperatures, the whole cooling time cannot be the same. 
//Thus I assumed that the system cools down by 5 degrees every 2min.//  whole_cooling_time /(int ((initial_temperature - final_temperature)/temperature_step) + 1);// nano second. Taken time to cool down the sample in each step (temperature_step). It should be in such a way that the total amount of time to cool the sample down from initial T to burn T be almost an hour.
double whole_cooling_time;
//2*60*60*1e9 * time_conversion_factor;
// 1 hour in nano second is 60*60*1e9 and in micro second is 60*60*1e6.  30*24*60*60*1e9 * time_conversion_factor;//
//const double cooling_time_step = cooling_time / 1.0; // in this way instead of using an aribitrary time step (1e-3) 
//I am using such a time step which results in having fixed row number (e.g. 5001).
//NOTICE: There is a problem with the choice of time interval and step due to the round-off. 
//For instance, for burn_time = 3.5 and burn_time_step = 0.1, A0_excited for last calculation is 0 because the loop time does not continue till the end, 
//meaning simulation_time (in this example, burn_time). It is better to choose the time step as a multiplier of the total time.
double burn_time; 
// micro second. As rates are in GHz, the simulation time should be in ns and the maximum of it should be as equal as the life time of excited state.
//const double burn_time_step = burn_time/1.0;
int number_of_burn;
//int(whole_burning_time/time_interval_between_two_consecutive_burns);//20000;//corresponding to a 37% hole. 
//multiplication of this number by burn_time gives the total desired burn time, e.g. 1e6 * 3.5 ns = 3.5 ms.
int number_of_HGK_measurement;//0.5*number_of_burn;
//NB: This number should be chosen in shuch a way that number_of_burn is divisible by it. 
//The number of measurement to create the HGK curve
int HGK_measurement_frequency;
//NB:There is a round off error when I use int(whole_burning_time/time_interval_between_two_consecutive_burns) which I did not spend time to catch it. 
//since measuring the hole depth after each act of burn consumes a lot of memory, this parameter is defined.
//It means, how frequently the hole depth should be measured (e.g. every 100 acts of burn)
//
int analytical_burn_calculation;
// 0 means calculating through multi-steps burning and 1 means calculating the result of post burn through analytical model.
int burn_in_one_step ;
// in case that I don't need HGK, I can use the analytical function to calculate the result of burning in one step.
int recovery_between_two_acts_of_burn;
//recovery within two consecutive absorption is activated if recovery_between_two_acts_of_burn = 1 otherwise it is not considered.
int recovery_correction_during_frequency_scan;
// 1 means considering the recovery for the time that system is not in resonance with the laser during frequency scan.

double whole_recovery_time;
//2*167*60*1e9 * time_conversion_factor;//1000*1e9 * time_conversion_factor;//167*60*1e9 * time_conversion_factor;//13*365*24*60*60*1e9 * time_conversion_factor;//time_interval_between_two_consecutive_burns*1e6;// an hour in nano second is 1 * 60*60*1e9.
int number_of_recovery_measurement;
// it means that during the whole_recovery_time, we will measure the spectum number_of_recovery_measurement (e.g. 60) times.
double recovery_time;
// it means that calculation for each step will take recovery_time.
//const double recovery_time_step = recovery_time / 1.0;

double number_of_longtime_measurement; // this is for some simulations of ultra-long time evolution...
double whole_long_time ;
// this time (13*365*24*60*60*1e9 = 13 years) is conidered as a longtime to calculate stationary state and then the recovery depth is compared with this distribution.
double long_time;
//const double long_time_step = long_time / 1.0;

double burn_temperature; //Kelvin.
int starting_from_equilibrium;
//It could be either 0 (non-equilibrium) or 1 (equilibrium). 
//This parameter determins whether the pre burn distribution is the result of cooling cycle (starting_from_equilibrium = 0) or is 
//Boltzman distribution at burn temperature (or equilibrium, starting_from_equilibrium = 1).
int starting_from_uniform_distribution;
//It could be either 0 (non uniform initial distribution, it could be either non-equilibrium or the results of cooling) or 1 (uniform distribution, i.e. 1/22).
int writing_details_of_burning_in_reson_well; 
// Writing of prob_in_reson_well_after_each_act_burn needs a great deal of memory. 0 means avoid doing that.

double highest_experimental_cycling_temperature;// This is the highest temperature during thermocycling experiment.
double highest_cycling_temperature;
// Since the first temperature in thermocycling loop should be greater than burn temperature which results in having final temperature higher than what we have 
//in practice (experimental_cycling_temperature).


double whole_thermocycling_time;   
// The measurment time for warming the sample up is about 250 minutes 
//(this is half of the whole thermocycling experiment that we usually do in the lab). This estimation is not completely correct since it does not take the same time for warming up and cooling down. In the first approximation, I assumed that temperature profile is linear and the same for both warming and cooling. Also, the sample does not stay at each cycling maximum temperature for the same time, however, I belive that it is not a big issue since the rate at low temperature (which sample reaches to them faster) is low enough that there is no significant difference between 5 or 20 minutes.
int number_of_thermocycling_measurement ;
//const double cycle_temperature_step = cycle_temperature-burn_temperature / number_of_thermocycling_measurement;
double thermocycling_time ;// it means that calculation for each step will take thermocycling_time.
//const double thermocycling_time_step = thermocycling_time / 1.0;
double recovery_temperature ;//Kelvin.
int number_of_row_burning ;
//burn_time/burn_time_step;// by omitting the loop in the SolveRateMatrix, there is no need to have this line. 
//the number of row in the rate_file (or rates_of_ensemble.txt) is related to the integer answer of dividing burn_time by time_step
int number_of_row_cooling ;
//cooling_time/cooling_time_step;
// the number of row in the rates_file.txt is related to the integer answer of dividing cooling_time by cooling_time_step 
//(this is the step numbers which takes the sample to be cooled down)
int number_of_row_recovery;//recovery_time/recovery_time_step;
int number_of_row_longtime;//long_time/long_time_step;
int number_of_row_thermocycling ;//thermocycling_time/thermocycling_time_step;


int number_of_column;// the first column is for simulation time so the number of column should be increased by 1

double curvature_of_parabola ; // 3 * 1e22 * energy_conversion_factor; //in MHz / m2
// MHz, the curvature of parabola. This value results in having smallest asymmetry delta of about 10 cm^-1 or ~300 GHz. 
//kd2=300 GHz, d=10-10, so k=3x1031 Hz/m2, or 3x10^22 Ghz/m2 or 1000 more MHz
int dimension ; 
double lambda_coefficient ; //3861418 in unknown units
int all_frequency ;// number of all possible frequencies which come out of the whole ensemble.
int total_bins ; // a fixed number of bins for creating spectra. sqrt (all_frequency) is very big for number of bins, 
//I changed it to ensemble number. Square-root choice is used in Excel and many others
int averaged_probability_SMS_bin_number ;

#endif // CONSTANTS_H
