#ifndef SOLVE_EIGEN_SYSTEM_H
#define SOLVE_EIGEN_SYSTEM_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition



void SolveEignesystem (gsl_matrix* a, gsl_vector_complex *eval, gsl_matrix_complex *evec, int &error);

#endif // SOLVE_EIGEN_SYSTEM_H
