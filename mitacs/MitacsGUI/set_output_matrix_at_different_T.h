#ifndef SET_OUTPUT_MATRIX_AT_DIFFERENT_T_H
#define SET_OUTPUT_MATRIX_AT_DIFFERENT_T_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition

void SetOutputAtDifferentTemp (gsl_matrix* output_at_different_T, gsl_matrix* rate_file, int& temperature_counter, int ensemble,
                             const int number_of_row /*either cooling or burning */);
#endif // SET_RATE_MATRIX_AT_DIFFERENT_T_H
