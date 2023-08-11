#ifndef SOLVE_RATE_MATRIX_H
#define SOLVE_RATE_MATRIX_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition


void SolveRateMatrix(gsl_matrix * rate_matrix, int ensemble, gsl_matrix* rate_file, const double simulation_time, /*const double time_step,*/
                     gsl_vector * A0, int number_of_row, double temperature/*, double initial_time, int &ODE_counter*/);


#endif // SOLVE_RATE_MATRIX_H
