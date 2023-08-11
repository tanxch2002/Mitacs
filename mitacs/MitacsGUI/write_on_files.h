#ifndef WRITE_ON_FILES_H
#define WRITE_ON_FILES_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition
#include <string>
#include <fstream>
#include <iostream>
#include "MitacsGUIDlg.h"
using namespace std;

void SetDirectory();
void WriteParabola (gsl_vector * parabola, string filename);
void WriteEnergyLandscape(gsl_matrix * landscape_file);
void WriteRateMatrix (gsl_matrix * rate_matrix);
void WriteA0 (gsl_vector * A0);
void WriteB (gsl_vector * B);//unnecessary
void WriteU (gsl_matrix * U);//unnecessary
void WriteConcentrationAfterSimulationTime (gsl_matrix * rate_file, double number_of_row, double time_step, string filename);
void WriteLastRowOfRateEquation (gsl_matrix * last_rates, string filename);
void WriteWholeFrequencyOfEnsemble(gsl_matrix * total_frequency, string filename);
void WriteDetailOfBinning(gsl_matrix * total_frequency, double bin_step, double smallest_frequency, string filename);
void WriteBinnedFrequency (gsl_vector * bin, double bin_step, double smallest_frequency, string filename);
void WriteDistribution (gsl_vector * random_numbers_binning, gsl_vector * potential_binning, gsl_vector * lambda_binning, string filename);
void WriteResonantWell (gsl_matrix * time_delay_in_frequency_scan);
void WriteProbabilitiesVsWellAfterEachScan (gsl_matrix * probabilities_in_each_scan, string filename);
void WriteProbabilitiesVsFrequencyAfterEachScan (gsl_vector * sort, gsl_matrix * probabilities_in_each_scan, string filename);
void WriteProbabilitiesInResonantWell (gsl_matrix * prob_in_reson_well_after_each_act_burn);
void WriteAvgProbabilitiesInResonantWell (gsl_matrix * averaged_probability_in_resonant_well);
void WriteAvgProbabilitiesInResonantWellVsFrequency (gsl_matrix * averaged_probability_in_resonant_well, gsl_vector * sorted_frequency);
void WriteBinnedAvgProProbabilitiesInResonantWellVsFrequency (gsl_matrix * averaged_probability_in_resonant_well, gsl_vector * sorted_frequency);
void WriteProbabilitiesAfterFrequencyScan (gsl_matrix * probabilities_in_each_frequency_step);
void WriteInitialWellNumber (gsl_vector * initial_well);
void WriteGeneralInformation (double calculation_time, double bin_step, double hole_burning_yield);
void WriteRateAtDifferentTemp (gsl_matrix* output_at_different_T);
void WriteRateElements (gsl_matrix* rate, double temperature);// this is only for dubuging purpose
void WriteLongtimeDist ( gsl_matrix * longtime_distribution, string filename);
void WriteInitialwellsAndMostProbableWells (gsl_matrix * boltzman_distribution_BT, gsl_vector * initial_well, gsl_matrix * last_rates, gsl_vector * deepest_well);
void WriteDiscrepancyOfEquilibrium (gsl_matrix * boltzman_distribution, gsl_matrix * output, string filename);
 void WriteProbabilityHeight (gsl_vector * probability_height, double post_burn_probability, string filename);
void WriteRecoveredSpectrum (gsl_matrix * recovered_probability, double smallest_frequency, double bin_step, string filename);
void WriteEquilibriumDifference (gsl_vector * equilibrium_ensemble, string filename);
void WriteBurningYieldvsStartingWell (gsl_vector * initial_well, gsl_vector * hole_burning_yield_starting_well, gsl_vector * pre_burn_at_starting_well,
                                          gsl_matrix * energy_difference);
 void WriteHoleGrowthKinetics (gsl_matrix * hole_growth_kinetics, gsl_vector * pre_burn_at_starting_well);
void WriteRateAtDifferentTime (gsl_matrix* output_at_different_recovery_time, const double number_of_measurement, string filename);
void WriteRatioProbability (gsl_matrix * longtime_distribution_BT, gsl_matrix * output_at_different_recovery_time, double post_burn_probability);
void WriteProbabilityRecoveredWell(gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well, gsl_matrix * last_rates_preburn,
                                   gsl_matrix *last_rates_postburn, const int number_of_measurement, string filename);
void WriteProbabilityInBurnedWell (gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well, gsl_matrix * last_rates_preburn, gsl_matrix *last_rates_postburn);
void WriteProbabilityInOnlyBurnedWell(gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well, gsl_matrix * last_rates_preburn,
                                      gsl_matrix *last_rates_postburn,  const int number_of_measurement, const double time/*either recovery_time or thermocycling_time*/, string filename);
void WriteRelativeHoleDepthRecovery (gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well,
                                     gsl_matrix * last_rates_preburn, gsl_matrix *last_rates_postburn);
void WriteAveragedProbabilityInOnlyBurnedWell(gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well,
                                              gsl_matrix * last_rates_preburn, gsl_matrix *last_rates_postburn, const int number_of_measurement,
                                              const double time/*either recovery_time or thermocycling_time*/, string filename);
void WriteAveragedProbabilityInDeepestWell(gsl_matrix * output_at_different_recovery_time, gsl_matrix * last_rates_preburn, gsl_matrix *last_rates_postburn,
                                      const int number_of_measurement, const double time, gsl_vector * deepest_well, gsl_vector * initial_well,  string filename);
void WriteProbabilitiesAtDifferentT (gsl_matrix * output_at_different_T, gsl_vector * initial_well);
void WriteAveragedOutputAtdifferentT (gsl_matrix * output_at_different_T);
void WriteTranslationEnergy (double translation);
void WriteAllConstants ();
void WriteCutOffEnergy (double cutoff_energy_difference);
void WriteBarrierHeightandParabolaforRate (gsl_vector * parabola_for_rate, string filename1, gsl_vector * barrier_height, string filename2, int ensemble, int recovery_counter);
#endif // WRITE_ON_FILES_H
