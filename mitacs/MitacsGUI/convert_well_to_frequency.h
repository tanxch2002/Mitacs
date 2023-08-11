#ifndef CONVERT_WELL_TO_FREQUENCY_H
#define CONVERT_WELL_TO_FREQUENCY_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition

void ConvertWellToFrequency (gsl_matrix * energy_difference, gsl_matrix * rate_file, gsl_matrix * last_rates, gsl_matrix * total_frequency,  int number_of_row);


#endif // CONVERT_WELL_TO_FREQUENCY_H
