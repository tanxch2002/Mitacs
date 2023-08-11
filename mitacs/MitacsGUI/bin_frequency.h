#ifndef BIN_FREQUENCY_H
#define BIN_FREQUENCY_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition

void BinFrequency (gsl_matrix * total_frequency, gsl_vector * bin, double& bin_step, double& smallest_frequency);

#endif // BIN_FREQUENCY_H
