#ifndef ORDER_WELL_VS_FREQUENCY_H
#define ORDER_WELL_VS_FREQUENCY_H

#include "Global_variables.h"
#include <gsl/gsl_matrix.h>

void OrderFrequenciesVsWell (gsl_matrix * energy_difference, gsl_matrix * sorted_frequency_vs_well, int ensemble);


#endif // ORDER_WELL_VS_FREQUENCY_H
