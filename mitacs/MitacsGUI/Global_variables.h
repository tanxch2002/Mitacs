
// this file is attached to the beginning of various cpp files to indicate that these variables are global
// and were defined somewhere and one has to look for that place.
// //
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>// for producing gaussian distribution in gsl
#include <gsl/gsl_math.h>// for eigensystem
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>// for producing gaussian distribution
#include <gsl/gsl_linalg.h>// for solving UB=A0 via linear algebra.
#include <string>
using namespace std;
#pragma once


#ifdef mainfile             // This portion applies to the MitacsGUIDlg.cpp only

//variables that contain energy landscapes modulated by randomness
	//New 2021: gsl-type variables that are shared between cool / burn / recovery subroutines in the mainfile
gsl_matrix* excited_parabola_for_rate; // zig-zag landscapes with randomness  we are not usiung them in this cycle...
gsl_matrix* ground_parabola_for_rate;
//
gsl_matrix* eigenve; //excited state eigenvalues - as of November 26, 2021 they are matrixes
gsl_matrix* eigenvg; //ground state eigenvalues
// eigenvectors are not passed to other programs, they are used internally in generate_energy_landscape.
//
// 3D arrays converted to 1D arrays.
gsl_vector* coeff_matrixe;
//coefficients c_n for expanding localized states into eigenstates
gsl_vector* coeff_matrixg ;
//
gsl_vector* attempt_freq_matrix_e;
gsl_vector* attempt_freq_matrix_g;
//
//T-independent one-act tunneling probabilities T(E):
gsl_vector* TE_matrix_ground;
gsl_vector* TE_matrix_excited;
//NB: no need for separate left to right and right to left matrices, one just needs to shift them by one, as these probabilities
// actually correspond to barriers, not wells
gsl_vector* A0_ground;
gsl_vector* A0_excited;
gsl_matrix* pre_burn;
gsl_matrix* post_burn;
gsl_matrix* post_burn_for_recovery;
gsl_vector* starting_well;
gsl_vector* v_barrier_height_excited;
gsl_vector* v_barrier_height_ground;
gsl_matrix* energy_difference; // all transition energies, defined as difference between energy levels with largest rhos
gsl_vector* bin; //supposedly contains probabilities binned according to frequencies, for making spectra
//
double const Boltzmann_constant = 0.69503;
// this is actually the conversion factor , 1 K = 0.69 cm-1, NOT a Boltzmann constant
//New 2021: So energy(cm-1) with Boltzman_constant(cm-1)
double pi = 3.141592;
double round_off = 0.0000001;
// as simulation time and step are defined as double, the round off issue makes problem and the last step doesn't work.
double round_off_limit = 1e-8;
// I may use this limit to set some very small numbers to zero.
double accuracy_measure = 0.01;
// too large
double extreme_small_number = 1e-70;

double energy_conversion_factor = 1e3;
// to convert energy (frequency) from GHz to MHz (by multipilication).
// All equation have been drived by assumption that energy is in GHz and time in nano second. 
//However, later Mehdi decided to change the scales to MHz and microseconds.
double time_conversion_factor = 1e-3;
// to convert time from nano second to micro second (by multipilication).
double light_speed = 2.99792458 * 1e8;
// m/s, the speed of light in vacuum
double h_planck = 6.626068 * 1e-34;
// J.S.  planck's constant
double k_b = 1.38e-23;
int max_ensemble_for_writing = 100;
// Some files were designed to check the validity
//const double w_0 = 1e3 * energy_conversion_factor; //MHz. i.e. 10^12 Hz.
CString situation;
double bin_step; // parameters of frequency binning.
double smallest_frequency;

double mass = 0;
double d_distance_of_two_wells = 0; //in nm
double width_well = 0; //the width of a well or thickness of a barrier
double inhomogeneous_bandwidth = 0;// cm^-1.  need conversion to MHz
double standard_deviation_of_translation = 0; // in MHz
double mean_value_of_translation = 0; // in cm-1, need conversion to MHz
double burn_frequency = 0;//frequency in cm-1, needs conversion to MHz.
double burn_power = 0;//in mW. 1microW is 1e-3, as an example. 5, 10, and 15 W are taken from Kohler's paper.
double beam_diameter = 0; //1/(2/sqrt(pi)); //cm. NB: In Kohler's paper the power is given in W/cm^2, so I have to conver the diameter in such a way that it gets rid of other constants in my equation (done on on September 12, 2012, in the second log book). In that calculation, area = pi r^2 = pi d^2/4. by assuming d = 2/sqrt(pi) the area would be 1cm^2. BTW sqrt(pi)/2 = 0.886 cm.
double absorption_cross_section_at_peak_RT = 0;//cm^2// September 12, 2012 (second log book)
double homogeneous_width_at_RT = 0;//cm^-1 for Cp43
double homogeneous_width_at_BT = 0;//cm^-1 for Cp43
double phonon_coupling = 0;// this is S or Huang–Rhys factor for Cp43.
double absorption_cross_section_at_peak_BT = 0;//=1.054e-12;//cm^2, see the notes on September 12, 2012 (in the second log book).
double photon_flux = 0; //burn frequency in MHz by now, Planck constant in SI units
double time_interval_between_two_consecutive_burns = 0; // in seconds
//time_interval_between_two_consecutive_burns = time_interval_between_two_consecutive_burns * 1e6; // in mks
double whole_burning_time = 0; //in seconds
double burn_time = 0; //excited state lifetime in microseconds, time for one act of burn
int number_of_burn = 0; // number of acts of burn in HGK, both time intervals are in seconds
int HGK_measurement_frequency = 0;//It means, how frequently the hole depth should be measured (e.g. every 100 acts of burn) 
int number_of_HGK_measurement = 0;
//The number of measurements to create the HGK curve

//cooling
double room_temperature = 0;//Kelvin.// it could be the same as initial_temperature.
	//parameters related to the cooling for finding the distribution before the start of burning.
double final_temperature = 0; // Kelvin.
double initial_temperature = 0; // Kelvin.
double temperature_step = 0;// the step of reducing of temperature. NB. Changing of this parameter should not be done unless cooling time changes properly.
double whole_cooling_time = 0; //in minutes
double cooling_time = 0; //cooling time at each cooling step
double burn_temperature = 0; // Kelvin.

//Single Molecule Experiment
int number_of_scan = 0; //one if we do not want to scan, only want the holeburn? Should we change it to zero?
double scan_time = 0;
double each_scan_time = 0;// micro seconds. 
double minimum_frequency = 0;

double maximum_frequency = 0;

double frequency_step = 0;

// in cm-1. Mehdi played with this value such that all wells are chosen once only.
int number_of_frequency_steps = 0;
double time_interval_between_each_frequency_step = 0; // in microseconds...

//Landscapes

int ensemble_number = 0; // number of ensembles used. More ensembles are tested until we find 5000 with resonant wells.
int Hamiltonian_well_numbers = 0; // Should be 512/1024/2048 for optimal matrix solving - the size of the high resolution energy landscape used in Miller-Colbert algorithm mentioned in Garashchuk paper
int well_numbers = 0;
int point_numbers = 0;
double interval = 0; // as of July 22, 2021 this is determined by the well thickness
//from center of the landscape to edge, in meters. Essentially it is a pseudo-shift of the parabola, allowing for loop 
double step = 0;
double mu = 0;
double sigma = 0;
double stretch = 0; //0.19; // ratio of excited and ground state BARRIERS
// Bottoms of the wells in the ground state
double mean = 0;
double deviation = 0;
double curvature_of_parabola = 0;
double lambda_coefficient = 0;
double well_STD_ratio = 6;
double cutoff_energy_difference;
// lambda_coefficient is the constant part of lambda, which is calculated in the createRateMatrix.cpp; 

//Recovery

double recovery_temperature = 0;//Kelvin.
double whole_recovery_time = 0; //in minutes
int number_of_recovery_measurement = 0;
// it means that during the whole_recovery_time, we will stop the evolution of probabilities and calculate the spectum number_of_recovery_measurement (e.g. 60) times.
double recovery_time = 0;
// it means that calculation for each step will take recovery_time.

int number_of_longtime_measurement = 1; // this is for simulations of ultra-long time evolution...
double whole_long_time = 0;
// this time (13*365*24*60*60*1e6 = 13 years in microseconds) 
//is conidered as a longtime to calculate stationary state and then the recovery depth is compared with this distribution.
double long_time = 0;
//const double long_time_step = long_time / 1.0;

double highest_experimental_cycling_temperature = 0;// This is the highest temperature during thermocycling experiment.
double highest_cycling_temperature = 0;
// Since the first temperature in thermocycling loop should be greater than burn temperature which results in having final temperature higher than what we have 
//in practice (experimental_cycling_temperature).

double whole_thermocycling_time = 0;
int number_of_thermocycling_measurement = 0;
//const double cycle_temperature_step = cycle_temperature-burn_temperature / number_of_thermocycling_measurement;
double thermocycling_time = 0;// it means that calculation for each step will take thermocycling_time.

int all_frequency = 0;// number of all possible frequencies which come out of the whole ensemble.
int total_bins; // a fixed number of bins for creating spectra. sqrt (all_frequency) is very big for number of bins, 
//I changed it to ensemble number. Square-root choice is used in Excel and many others
int averaged_probability_SMS_bin_number = 500;

//Checkboxes:

int generate_energy_landscape = 1;
//1 means that energy landscape is generated and 0 means that energy landscape (and Translation in case that bottom of wells are independent in two states) is/are imported.
int analytical_burn_calculation = 0;
// 0 means calculating through multi-steps burning and 1 means calculating the result of post burn through analytical model.
int burn_in_one_step = 0;
// in case that I don't need HGK, I can use the analytical function to calculate the result of burning in one step.
int recovery_between_two_acts_of_burn = 0;
//recovery within two consecutive absorption is activated if recovery_between_two_acts_of_burn = 1 otherwise it is not considered.
int recovery_correction_during_frequency_scan = 0;
// 1 means considering the recovery for the time that system is not in resonance with the laser during frequency scan.
int writing_details_of_burning_in_reson_well = 0;
// Writing of prob_in_reson_well_after_each_act_burn needs a great deal of memory. 0 means avoid doing that.
int starting_from_equilibrium = 0;
int starting_from_uniform_distribution = 0;
int barrier_decoupling = 1;
int bottom_decoupling = 1;
int balancing = 1;
int binning = 0;

double scaling_factor = 1;
double scaling_parameters = 1;
double extra_temperature = 0; 

//Something related to writing to files, did not change so far
int number_of_row_burning = 1;
//burn_time/burn_time_step;// by omitting the loop in the SolveRateMatrix, there is no need to have this line. 
//the number of row in the rate_file (or rates_of_ensemble.txt) is related to the integer answer of dividing burn_time by time_step
int number_of_row_cooling = 1;
//cooling_time/cooling_time_step;
// the number of row in the rates_file.txt is related to the integer answer of dividing cooling_time by cooling_time_step 
//(this is the step numbers which takes the sample to be cooled down)
int number_of_row_recovery = 1;//recovery_time/recovery_time_step;
int number_of_row_longtime = 1;//long_time/long_time_step;
int number_of_row_thermocycling = 1;//thermocycling_time/thermocycling_time_step;
int number_of_column = 0;// the first column is for simulation time so the number of column should be increased by 1

double sum_resonant_probabilities_preburn;
double sum_resonant_probabilities_postburn;
double sum_resonant_probabilities_recovery;


//the following functions are declared as extern here since they are external to the main file

extern void GenerateEnergyLandscape(gsl_vector* coeff_matrixe, gsl_vector* coeff_matrixg, // coefficients c_n
    gsl_matrix* eigenve, gsl_matrix* eigenvg, //eigenvalues
    gsl_vector* TE_matrix_excited, gsl_vector* TE_matrix_ground, // transmission probabilities
    gsl_vector* classical_attempt_freq_matrix_e, gsl_vector* classical_attempt_freq_matrix_g, //attempt frequencies
    gsl_vector* parabola_ground, gsl_matrix* energy_difference, int ensemble, gsl_vector* starting_well, 
    int& binning, gsl_vector* v_barrier_height_excited, gsl_vector* v_barrier_height_ground,
    gsl_matrix* excited_parabola_for_rate, gsl_matrix* ground_parabola_for_rate, //low-res energy landscapes
    gsl_rng* r, int& produced_ensemble, gsl_vector* trans, double cutoff_energy_difference);

extern void CalculateBoltzmanDist(gsl_matrix* parabola_for_rate, gsl_matrix* boltzmn_distribution, int column_number, double temperature);

extern void CreateRateMatrix(gsl_vector* coeff_matrix, gsl_matrix* eigenv, gsl_vector* attempt_frequency, gsl_vector* TE,
    gsl_matrix* parabola_for_rate, gsl_matrix* rate_matrix, gsl_matrix* temp_matrix, double temperature, int ensemble, gsl_vector* v_barrier_height);

extern void SolveRateMatrix(gsl_matrix* rate_matrix, int ensemble, gsl_matrix* rate_file, const double simulation_time, /*const double time_step,*/
    gsl_vector* A0, int number_of_row, double temperature/*, double initial_time, int &ODE_counter*/);

extern void SetOutputAtDifferentTemp(gsl_matrix* output_at_different_T, gsl_matrix* rate_file, int& temperature_counter, int ensemble,
    const int number_of_row /*either cooling or burning */);

extern void SolveEignesystem(gsl_matrix* a, gsl_vector_complex* eval, gsl_matrix_complex* evec, int& error);

extern void CheckSituation(double temperature, CString &situation);

extern int to1D(int x, int y, int z);

//extern void LoadExperimentalThermocyclingProfile(gsl_matrix* experimental_thermocycling_time_temperature);
extern void CheckFrequencyRange(gsl_matrix* en_difference);
extern void FindResonantWellandTimeDelay(gsl_matrix* en_difference, gsl_matrix* time_delay_in_frequency_scan, int ensemble);
extern void BinFrequency(gsl_matrix* total_frequency, gsl_vector* bin, double& bin_step, double& smallest_frequency);
extern void SetTotalFrequencies(gsl_matrix* wells/*either hole_wells or longtime_distribution*/, gsl_matrix* total_frequency_hole);
extern void OrderFrequenciesVsWell(gsl_matrix* ene_difference, gsl_matrix* sorted_frequency_vs_well, int ensemble);
extern void ConvertWellToFrequency(gsl_matrix* energy_difference, gsl_matrix* rate_file, gsl_matrix* last_rates, gsl_matrix* total_frequency, int number_of_row);
extern void SetDirectory();

extern void WriteParabola(gsl_vector* parabola, string filename);
extern void WriteEnergyLandscape(gsl_matrix* landscape_file);
extern void WriteRateMatrix(gsl_matrix* rate_matrix);
extern void WriteA0(gsl_vector* A0);
extern void WriteB(gsl_vector* B);//unnecessary
extern void WriteU(gsl_matrix* U);//unnecessary
extern void WriteConcentrationAfterSimulationTime(gsl_matrix* rate_file, double number_of_row, double time_step, string filename);
extern void WriteLastRowOfRateEquation(gsl_matrix* last_rates, string filename);
extern void WriteWholeFrequencyOfEnsemble(gsl_matrix* total_frequency, string filename);
extern void WriteDetailOfBinning(gsl_matrix* total_frequency, double bin_step, double smallest_frequency, string filename);
extern void WriteBinnedFrequency(gsl_vector* bin, double bin_step, double smallest_frequency, string filename);
extern void WriteDistribution(gsl_vector* random_numbers_binning, gsl_vector* potential_binning, gsl_vector* lambda_binning, string filename);
extern void WriteResonantWell(gsl_matrix* time_delay_in_frequency_scan);
extern void WriteProbabilitiesVsWellAfterEachScan(gsl_matrix* probabilities_in_each_scan, string filename);
extern void WriteProbabilitiesVsFrequencyAfterEachScan(gsl_vector* sort, gsl_matrix* probabilities_in_each_scan, string filename);
extern void WriteProbabilitiesInResonantWell(gsl_matrix* prob_in_reson_well_after_each_act_burn);
extern void WriteAvgProbabilitiesInResonantWell(gsl_matrix* averaged_probability_in_resonant_well);
extern void WriteAvgProbabilitiesInResonantWellVsFrequency(gsl_matrix* averaged_probability_in_resonant_well, gsl_vector* sorted_frequency);
extern void WriteBinnedAvgProProbabilitiesInResonantWellVsFrequency(gsl_matrix* averaged_probability_in_resonant_well, gsl_vector* sorted_frequency);
extern void WriteProbabilitiesAfterFrequencyScan(gsl_matrix* probabilities_in_each_frequency_step);
extern void WriteInitialWellNumber(gsl_vector* initial_well);
extern void WriteGeneralInformation(double calculation_time, double bin_step, double hole_burning_yield);
extern void WriteRateAtDifferentTemp(gsl_matrix* output_at_different_T);
extern void WriteRateElements(gsl_matrix* rate, double temperature);// this is only for dubuging purpose
extern void WriteLongtimeDist(gsl_matrix* longtime_distribution, string filename);
extern void WriteInitialwellsAndMostProbableWells(gsl_matrix* boltzman_distribution_BT, gsl_vector* initial_well, gsl_matrix* last_rates, gsl_vector* deepest_well);
extern void WriteDiscrepancyOfEquilibrium(gsl_matrix* boltzman_distribution, gsl_matrix* output, string filename);
extern void WriteProbabilityHeight(gsl_vector* probability_height, double post_burn_probability, string filename);
extern void WriteRecoveredSpectrum(gsl_matrix* recovered_probability, double smallest_frequency, double bin_step, string filename);
extern void WriteEquilibriumDifference(gsl_vector* equilibrium_ensemble, string filename);
extern void WriteBurningYieldvsStartingWell(gsl_vector* initial_well, gsl_vector* hole_burning_yield_starting_well, gsl_vector* pre_burn_at_starting_well,
    gsl_matrix* energy_difference);
extern void WriteHoleGrowthKinetics(gsl_matrix* hole_growth_kinetics, gsl_vector* pre_burn_at_starting_well);
extern void WriteRateAtDifferentTime(gsl_matrix* output_at_different_recovery_time, const double number_of_measurement, string filename);
extern void WriteRatioProbability(gsl_matrix* longtime_distribution_BT, gsl_matrix* output_at_different_recovery_time, double post_burn_probability);
extern void WriteProbabilityRecoveredWell(gsl_matrix* output_at_different_recovery_time, gsl_vector* initial_well, gsl_matrix* last_rates_preburn,
    gsl_matrix* last_rates_postburn, const int number_of_measurement, string filename);
extern void WriteProbabilityInBurnedWell(gsl_matrix* output_at_different_recovery_time, gsl_vector* initial_well, gsl_matrix* last_rates_preburn, gsl_matrix* last_rates_postburn);
extern void WriteProbabilityInOnlyBurnedWell(gsl_matrix* output_at_different_recovery_time, gsl_vector* initial_well, gsl_matrix* last_rates_preburn,
    gsl_matrix* last_rates_postburn, const int number_of_measurement, const double time/*either recovery_time or thermocycling_time*/, string filename);
extern void WriteRelativeHoleDepthRecovery(gsl_matrix* output_at_different_recovery_time, gsl_vector* initial_well,
    gsl_matrix* last_rates_preburn, gsl_matrix* last_rates_postburn);
extern void WriteAveragedProbabilityInOnlyBurnedWell(gsl_matrix* output_at_different_recovery_time, gsl_vector* initial_well,
    gsl_matrix* last_rates_preburn, gsl_matrix* last_rates_postburn, const int number_of_measurement,
    const double time/*either recovery_time or thermocycling_time*/, string filename);
extern void WriteAveragedProbabilityInDeepestWell(gsl_matrix* output_at_different_recovery_time, gsl_matrix* last_rates_preburn, gsl_matrix* last_rates_postburn,
    const int number_of_measurement, const double time, gsl_vector* deepest_well, gsl_vector* initial_well, string filename);
extern void WriteProbabilitiesAtDifferentT(gsl_matrix* output_at_different_T, gsl_vector* initial_well);
extern void WriteAveragedOutputAtdifferentT(gsl_matrix* output_at_different_T);
extern void WriteTranslationEnergy(double translation);
extern void WriteAllConstants();
extern void WriteCutOffEnergy(double cutoff_energy_difference);
extern void WriteBarrierHeightandParabolaforRate(gsl_vector* parabola_for_rate, string filename1, gsl_vector* barrier_height, string filename2, int ensemble, int recovery_counter);
//
#else                        // Applies for all other files (that are not mainfile / MitacsGUIDlg.cpp)

//New 2021: gsl-type variables that are shared between cool / burn / recovery subroutines in the mainfile
extern gsl_matrix* excited_parabola_for_rate; // zig-zag landscapes with randomness  we are not usiung them in this cycle...
extern gsl_matrix* ground_parabola_for_rate;
//
extern gsl_matrix* eigenve; //excited state eigenvalues - as of November 26, 2021 they are matrixes
extern gsl_matrix* eigenvg; //ground state eigenvalues
// eigenvectors are not passed to other programs, they are used internally in generate_energy_landscape.
//
// 3D arrays converted to 1D arrays.
extern gsl_vector* coeff_matrixe;
//coefficients c_n for expanding localized states into eigenstates
extern gsl_vector* coeff_matrixg;
//
extern gsl_vector* attempt_freq_matrix_e;
extern gsl_vector* attempt_freq_matrix_g;
//
//T-independent one-act tunneling probabilities T(E):
extern gsl_vector* TE_matrix_ground;
extern gsl_vector* TE_matrix_excited;
//NB: no need for separate left to right and right to left matrices, one just needs to shift them by one, as these probabilities
// actually correspond to barriers, not wells
extern gsl_vector* A0_ground;
extern gsl_vector* A0_excited;
extern gsl_matrix* pre_burn;
extern gsl_matrix* post_burn;
extern gsl_matrix* post_burn_for_recovery;
extern gsl_vector* v_barrier_height_excited;
extern gsl_vector* v_barrier_height_ground;
extern gsl_matrix* energy_difference;

extern CString situation;
extern double const Boltzmann_constant;
// this is actually the conversion factor , 1 K = 0.69 cm-1, NOT a Boltzmann constant
//New 2021: So energy(cm-1) with Boltzman_constant(cm-1)
extern double pi;
extern double  round_off;
// as simulation time and step are defined as double, the round off issue makes problem and the last step doesn't work.
extern double  round_off_limit;
// I may use this limit to set some very small numbers to zero.
extern double  accuracy_measure; // for solving Rate Matrix
// too large
extern double  extreme_small_number;

extern double energy_conversion_factor;
// to convert energy (frequency) from GHz to MHz (by multipilication).
// All equation have been drived by assumption that energy is in GHz and time in nano second. 
//However, later Mehdi decided to change the scales to MHz and microseconds.
extern double  time_conversion_factor;
// to convert time from nano second to micro second (by multipilication).
extern double  light_speed;
// m/s, the speed of light in vacuum
extern double  h_planck;
// J.S.  planck's constant
extern double   k_b;
extern int  max_ensemble_for_writing;
// Some files were designed to check the validity
//const double w_0 = 1e3 * energy_conversion_factor; //MHz. i.e. 10^12 Hz.
extern double mass; //in kg
// decreased from 100 amu to 10 amu (look at my log book volume II, Feb 8, 2012) and one amu is 1.66* 10^-27 kg

extern double d_distance_of_two_wells;
//In this version this is the "primary" quantity (one accessible to user via interface) and other quantities depend on this one
//Thickness of the barrier and width of the well in the model with rectangular barriers.
extern double width_well;
//the width of well or thickness of the barrier

extern int ensemble_number; // number of ensembles used. More ensembles are tested until we find 5000 with resonant wells.

extern int Hamiltonian_well_numbers; // Should be 512/1024/2048 for optimal matrix solving - the size of the high resolution energy landscape used in Miller-Colbert algorithm mentioned in Garashchuk paper
                                          //
extern int well_numbers;
extern int point_numbers;
/* the points which are needed for plotting the main curve (parabola). The distance between two adjacent points is 2*interval/(point_numbers-1)
                               (or (max-min)/(point_numbers-1) if the interval is not symmetrical with recpect to zero(. For example we
                               want to have 5 points for drawing parabola so each step should be one unit ([2*2]/[5-1]) to span -2 to 2 (our interval in this case) is five steps.*/
extern double interval;
// as of July 22, 2021 this is determined by the well thickness
//from center of the landscape to edge, in meters. Essentially it is a pseudo-shift of the parabola, allowing for loop 

extern double step; // in general, (max-min)/(point_numbers);


//coupling between the ground and excited state.
extern int bottom_decoupling;
//if this value is 1 then a new set of random variable is used for the bottom of wells in the excited state; otherwise, the same value from the ground state is used (by multipling it by stretch)
extern int barrier_decoupling;
//if this value is 1 then a new set of random variable is used for the the barrier heights in the excited state; otherwise, the same value from the ground state is used (by multipling it by stretch)
//
// double burn_wavelength = 686.8131; //nm  or 14560 cm-1 or 4.368 x10^14 Hz
//
//
//parameters for distribution of energy gap between the ground and excited states baselines
extern double inhomogeneous_bandwidth;
// cm^-1.  need conversion to MHz
extern double standard_deviation_of_translation;
// in MHz
// the relationship between bandwidth in cm^-1 and standard deviation is FWHM = 2.355 * sigma. 
//It should be noted that 1 cm^-1 = 29.9792458 GHz then FWHM=180*29.97 = 5394.6 GHz.
extern double mean_value_of_translation;
// in cm-1, need conversion to MHz
//based on the last paper burn frequency was 14560 wavenumbers or 686.8131 nm 
//438794.499 and 5995.84916 GHz are corresponding to 14641 and 180 wavenumbers from the SDF peak and width of CP43, 
//A-state, respectively, in the last paper. 180 comes from measuring inhomogeneous line broadening experiments.
extern double burn_frequency;
//frequency in cm-1, needs conversion to MHz.

extern double burn_power;
//mW. 1microW is 1e-3, as an example. 5, 10, and 15 W are taken from Kohler's paper.
extern double beam_diameter;
//1/(2/sqrt(pi)); //cm. NB: In Kohler's paper the power is given in W/cm^2, so I have to conver the diameter in such a way that it gets rid of other constants in my equation (done on on September 12, 2012, in the second log book). In that calculation, area = pi r^2 = pi d^2/4. by assuming d = 2/sqrt(pi) the area would be 1cm^2. BTW sqrt(pi)/2 = 0.886 cm.
extern double absorption_cross_section_at_peak_RT;
//cm^2// September 12, 2012 (second log book)
extern double homogeneous_width_at_RT;
//cm^-1 for Cp43
extern double homogeneous_width_at_BT;
//cm^-1 for Cp43
extern double phonon_coupling;
// this is S or Huang–Rhys factor for Cp43.
extern double absorption_cross_section_at_peak_BT;
//=1.054e-12;//cm^2, see the notes on September 12, 2012 (in the second log book).
extern double photon_flux;
//burn frequency in MHz by now, Planck constant in SI units
extern  double time_interval_between_two_consecutive_burns;
//1*(beam_diameter*beam_diameter/(absorption_cross_section_at_peak_BT*burn_power*burn_wavelength))*1.56015e-13;
//in seconds, see the notes on September 12, 2012 (in the second log book).
extern double whole_burning_time;
//scan_time/well_numbers ; // in second, for burning experiment it was 19min = 19 *60 seconds in one of our experiments.

//scan parameters:
extern double scan_time;
//second. For SMS, 2.6 second comes from the scan speed of 44cm^-1/s (Kohler's paper). 
//Our range here is approximately between 434733-438215GHz (or 14505-14621 wavenumbers) so delta_t = 116/44 = 2.6 the time to scan the whole range. I assumed that laser stays there for 2.6/well_numbers as an approximate value.

extern int number_of_scan;
//one if we do not want to scan, only want the holeburn?
extern double minimum_frequency;
// 380000 * energy_conversion_factor; // MHz or 14505.5 cm^-1 or 689.4nm. 
//This the minimum frequency that I found it postburn_spectrum.txt in a typlical hole burning calculation for flat landscape with 3000 systems.
extern double maximum_frequency;
//460000 * energy_conversion_factor; // MHz or 14621.7 cm^-1 or 683.9nm. 
//This the maximum frequency that I found it postburn_spectrum.txt in a typlical hole burning calculation for flat landscape with 3000 systems.
extern double frequency_step;
// in MHz. I played with this value such that all wells are chosen once only.
extern int number_of_frequency_steps;
extern double each_scan_time;
// micro second. 
//It means that 2.6 seconds takes to scan the whole band from minimum_frequency to maximum_frequency 
//(see the note for whole_burning_time, or my notes in the second log book on April 9th, 2013).
extern double time_interval_between_each_frequency_step;
//
extern double scaling_factor;
// this factor scales the energy landscape (barrier heights and bottom of wells) in the ground state 
//which in turn affects the barrier heights in the excited state.
extern double scaling_parameters;
// this is the same as scaling_factor but I want to use this parameter 
//to change mu, sigma, mean, and deviation at the same time. 
//The former factor is reserved for scaling an energy landscape when it is imported in the program.
//
//
// Ground state barrier distribution:
extern double mu;
// this is currently 37,100,000 MHz
// or 37100 GHz or 1233 cm-1
// mu should be around 1000 wavenumbers (look at my log book volume II). 
//This corresponds to the mean value of lambda close to 22 in the ground state.
extern double sigma;
//NOTICE: Following value is used in most of simulatios: 7990.0 * scaling_parameters * energy_conversion_factor;//mu/4.5;
// I changed this parameter to find lambda distribution in the ground state close to 2 based on the last paper. 
//When denominator is less than 5, 
//random generator makes negative values for potential which is consequently become zero in lambda. 
//This makes non gaussian distribution for lambda.
extern double stretch;
// ratio of excited and ground state BARRIERS
//The stretch should never exceed 1. This resuls in the mean value of lambda close to 11 (10.94) and sigma almost 1.1 in the excited state. 
//Changing this value affects the average energy difference between the bottom of two wells in the ground state 
// i.e. "asymmetry". In the ground state it should be ~10 wavenumbers. see also Curvature Calculation.nb
//
extern double well_STD_ratio;
//
// Bottoms of the wells in the ground state
extern double mean;
// the mean value for the bottom of wells; 600 GHz is 20 cm-1
extern double deviation;
// the standard deviation value for the bottom of wells
//
extern int generate_energy_landscape;
//1 means that energy landscape is generated and 0 means that energy landscape (and Translation in case that bottom of wells are independent in two states) is/are imported.
extern double room_temperature;
//Kelvin.// it could be the same as initial_temperature.
extern double cutoff_energy_difference;

//parameters related to the cooling for finding the distribution before the start of burning.
extern double extra_temperature;
// I would like to add extra temperature to measure HGK at different temperatures while not changing the temperature step. 
//It should be noted that the distribution at higher temperatures is close to each other.
extern double initial_temperature; //Kelvin.
extern double final_temperature; // Kelvin.
extern double temperature_step; // the step of reducing of temperature. NB. Changing of this parameter should not be done unless cooling time changes properly.
extern double cooling_time;  // 2 minutes, time of one cooling step
// In case of burning at different temperatures, the whole cooling time cannot be the same. 
//Thus I assumed that the system cools down by 5 degrees every 2min.//  whole_cooling_time /(int ((initial_temperature - final_temperature)/temperature_step) + 1);// nano second. Taken time to cool down the sample in each step (temperature_step). It should be in such a way that the total amount of time to cool the sample down from initial T to burn T be almost an hour.
extern double whole_cooling_time;
//2*60*60*1e9 * time_conversion_factor;
// 1 hour in nano second is 60*60*1e9 and in micro second is 60*60*1e6.  30*24*60*60*1e9 * time_conversion_factor;//
//const double cooling_time_step = cooling_time / 1.0; // in this way instead of using an aribitrary time step (1e-3) 
//I am using such a time step which results in having fixed row number (e.g. 5001).
//NOTICE: There is a problem with the choice of time interval and step due to the round-off. 
//For instance, for burn_time = 3.5 and burn_time_step = 0.1, A0_excited for last calculation is 0 because the loop time does not continue till the end, 
//meaning simulation_time (in this example, burn_time). It is better to choose the time step as a multiplier of the total time.
extern double burn_time;
// micro second. As rates are in GHz, the simulation time should be in ns and the maximum of it should be as equal as the life time of excited state.
//const double burn_time_step = burn_time/1.0;
extern int number_of_burn;
//int(whole_burning_time/time_interval_between_two_consecutive_burns);//20000;//corresponding to a 37% hole. 
//multiplication of this number by burn_time gives the total desired burn time, e.g. 1e6 * 3.5 ns = 3.5 ms.
extern int number_of_HGK_measurement;
//0.5*number_of_burn;
//NB: This number should be chosen in shuch a way that number_of_burn is divisible by it. 
//The number of measurement to create the HGK curve
extern int HGK_measurement_frequency;
//NB:There is a round off error when I use int(whole_burning_time/time_interval_between_two_consecutive_burns) which I did not spend time to catch it. 
//since measuring the hole depth after each act of burn consumes a lot of memory, this parameter is defined.
//It means, how frequently the hole depth should be measured (e.g. every 100 acts of burn)
//
extern int analytical_burn_calculation;
// 0 means calculating through multi-steps burning and 1 means calculating the result of post burn through analytical model.
extern int burn_in_one_step;
// in case that I don't need HGK, I can use the analytical function to calculate the result of burning in one step.
extern int recovery_between_two_acts_of_burn;
//recovery within two consecutive absorption is activated if recovery_between_two_acts_of_burn = 1 otherwise it is not considered.
extern  int recovery_correction_during_frequency_scan;
// 1 means considering the recovery for the time that system is not in resonance with the laser during frequency scan.

extern  double whole_recovery_time;
//2*167*60*1e9 * time_conversion_factor;//1000*1e9 * time_conversion_factor;//167*60*1e9 * time_conversion_factor;//13*365*24*60*60*1e9 * time_conversion_factor;//time_interval_between_two_consecutive_burns*1e6;// an hour in nano second is 1 * 60*60*1e9.
extern int number_of_recovery_measurement;
// it means that during the whole_recovery_time, we will measure the spectum number_of_recovery_measurement (e.g. 60) times.
extern double recovery_time;
// it means that calculation for each step will take recovery_time.
//const double recovery_time_step = recovery_time / 1.0;

extern int number_of_longtime_measurement; // this is for some simulations of ultra-long time evolution...
extern double whole_long_time;
// this time (13*365*24*60*60*1e9 = 13 years) is conidered as a longtime to calculate stationary state and then the recovery depth is compared with this distribution.
extern double long_time;
//const double long_time_step = long_time / 1.0;

extern  double burn_temperature;
//Kelvin.
extern int starting_from_equilibrium;
//It could be either 0 (non-equilibrium) or 1 (equilibrium). 
//This parameter determins whether the pre burn distribution is the result of cooling cycle (starting_from_equilibrium = 0) or is 
//Boltzman distribution at burn temperature (or equilibrium, starting_from_equilibrium = 1).
extern int starting_from_uniform_distribution;
//It could be either 0 (non uniform initial distribution, it could be either non-equilibrium or the results of cooling) or 1 (uniform distribution, i.e. 1/22).
extern int writing_details_of_burning_in_reson_well;
// Writing of prob_in_reson_well_after_each_act_burn needs a great deal of memory. 0 means avoid doing that.

extern double highest_experimental_cycling_temperature;
// This is the highest temperature during thermocycling experiment.
extern double highest_cycling_temperature;
// Since the first temperature in thermocycling loop should be greater than burn temperature which results in having final temperature higher than what we have 
//in practice (experimental_cycling_temperature).


extern double whole_thermocycling_time;
// The measurment time for warming the sample up is about 250 minutes 
//(this is half of the whole thermocycling experiment that we usually do in the lab). This estimation is not completely correct since it does not take the same time for warming up and cooling down. In the first approximation, I assumed that temperature profile is linear and the same for both warming and cooling. Also, the sample does not stay at each cycling maximum temperature for the same time, however, I belive that it is not a big issue since the rate at low temperature (which sample reaches to them faster) is low enough that there is no significant difference between 5 or 20 minutes.
extern int number_of_thermocycling_measurement;
//const double cycle_temperature_step = cycle_temperature-burn_temperature / number_of_thermocycling_measurement;
extern  double thermocycling_time;
// it means that calculation for each step will take thermocycling_time.
//const double thermocycling_time_step = thermocycling_time / 1.0;
extern double recovery_temperature;
//Kelvin.
extern  int number_of_row_burning;
//burn_time/burn_time_step;// by omitting the loop in the SolveRateMatrix, there is no need to have this line. 
//the number of row in the rate_file (or rates_of_ensemble.txt) is related to the integer answer of dividing burn_time by time_step
extern int number_of_row_cooling;
//cooling_time/cooling_time_step;
// the number of row in the rates_file.txt is related to the integer answer of dividing cooling_time by cooling_time_step 
//(this is the step numbers which takes the sample to be cooled down)
extern int number_of_row_recovery;
//recovery_time/recovery_time_step;
extern int number_of_row_longtime;
//long_time/long_time_step;
extern int number_of_row_thermocycling;
//thermocycling_time/thermocycling_time_step;
extern double sum_resonant_probabilities_preburn;
extern double sum_resonant_probabilities_postburn;
extern double sum_resonant_probabilities_recovery;

extern int number_of_column;
// the first column is for simulation time so the number of column should be increased by 1

extern double curvature_of_parabola;
// 3 * 1e22 * energy_conversion_factor; //in MHz / m2
// MHz, the curvature of parabola. This value results in having smallest asymmetry delta of about 10 cm^-1 or ~300 GHz. 
//kd2=300 GHz, d=10-10, so k=3x1031 Hz/m2, or 3x10^22 Ghz/m2 or 1000 more MHz
extern int well_numbers;
extern double lambda_coefficient;
//3861418 in unknown units
extern int all_frequency;
// number of all possible frequencies which come out of the whole ensemble.
extern  int total_bins;
// a fixed number of bins for creating spectra. sqrt (all_frequency) is very big for number of bins, 
//I changed it to ensemble number. Square-root choice is used in Excel and many others
extern int averaged_probability_SMS_bin_number;
extern int barrier_decoupling;
extern int bottom_decoupling;
extern int balancing;
extern int binning;

extern int to1D(int x, int y, int z);

#endif

//#endif
//#pragma once