#ifndef GENERATE_ENERGY_LANDSCAPE_H
#define GENERATE_ENERGY_LANDSCAPE_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition
#include <gsl/gsl_rng.h>// for producing gaussian distribution

void GenerateEnergyLandscape(gsl_vector* coeff_matrixe, gsl_vector* coeff_matrixg, // coefficients c_n
    gsl_matrix* eigenve, gsl_matrix* eigenvg, //eigenvalues
    gsl_vector* TE_matrix_excited, gsl_vector* TE_matrix_ground, // transmission probabilities
    gsl_vector* classical_attempt_freq_matrix_e, gsl_vector* classical_attempt_freq_matrix_g, //attempt frequencies
    gsl_vector* parabola_ground, gsl_matrix* energy_difference, int ensemble, gsl_vector* starting_well, 
    int& binning, gsl_vector* v_barrier_height_excited, gsl_vector* v_barrier_height_ground,
    gsl_matrix* excited_parabola_for_rate, gsl_matrix* ground_parabola_for_rate, //low-res energy landscapes
    gsl_rng* r, int& produced_ensemble, gsl_vector* trans, double cutoff_energy_difference);

#endif // GENERATE_ENERGY_LANDSCAPE_H
