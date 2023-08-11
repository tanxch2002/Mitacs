#ifndef CREATE_RATE_MATRIX_H
#define CREATE_RATE_MATRIX_H
#include <gsl/gsl_matrix.h>
void CreateRateMatrix(gsl_vector* coeff_matrix, gsl_matrix* eigenv, gsl_vector * attempt_frequency, gsl_vector* TE,
                 gsl_matrix * parabola_for_rate, gsl_matrix * rate_matrix, gsl_matrix * temp_matrix,
                      double temperature, int ensemble, gsl_vector * v_barrier_height);

#endif // CREATE_RATE_MATRIX_H
