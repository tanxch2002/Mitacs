#ifndef SOLVE_ODE_H
#define SOLVE_ODE_H
#include <gsl/gsl_matrix.h>
void SolveODE (gsl_vector * A0, double simulation_time);
#endif // SOLVE_ODE_H
