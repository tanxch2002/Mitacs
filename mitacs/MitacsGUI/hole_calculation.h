#ifndef HOLE_CALCULATION_H
#define HOLE_CALCULATION_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition
void SetHoleProbability (gsl_vector * pre_burn, gsl_vector * post_burn, gsl_matrix * hole_wells, int ensemble);
void SetTotalFrequencies (gsl_matrix * hole_wells, gsl_matrix * total_frequency_hole);
#endif // HOLE_CALCULATION_H
