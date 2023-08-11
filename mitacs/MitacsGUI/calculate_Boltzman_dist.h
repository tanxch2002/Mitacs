#ifndef CALCULATE_BOLTZMAN_DIST_H
#define CALCULATE_BOLTZMAN_DIST_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition
#include "MitacsGUIDlg.h"

void CalculateBoltzmanDist(gsl_matrix * parabola_for_rate, gsl_matrix * boltzmn_distribution, int column_number, double temperature);
#endif // CALCULATE_BOLTZMAN_DIST_H
