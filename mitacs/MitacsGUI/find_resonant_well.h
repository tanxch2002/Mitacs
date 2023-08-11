#ifndef FIND_RESONANT_WELL_H
#define FIND_RESONANT_WELL_H
#include <gsl/gsl_matrix.h> // it is necessary for the matrix and vector definition in gsl.
#include "Global_variables.h"
void FindResonantWellandTimeDelay (gsl_matrix * energy_difference,  gsl_matrix * time_delay_in_frequency_scan, int ensemble);


#endif // FIND_RESONANT_WELL_H
