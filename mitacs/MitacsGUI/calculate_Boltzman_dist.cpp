#include "pch.h"
#include "calculate_Boltzman_dist.h"
#include "Global_variables.h"
//using namespace std;

void CalculateBoltzmanDist(gsl_matrix* parabola_for_rate, gsl_matrix* boltzmn_distribution, int row_number/*ensemble number*/, double temperature) {
    // modified by VZ in 2023 to ensure that units are ok
    // now boltzmann constant is in SI, so energies have to be converted from MHz to Joules

    double delta, test;
    double* nn;
    double boltzmann_dist;
    double normalization_coefficient = 0;
    gsl_vector* normalization = gsl_vector_calloc(boltzmn_distribution->size1);
    //
    gsl_vector* temp = gsl_vector_calloc(point_numbers);
    gsl_matrix_get_row(temp, parabola_for_rate, row_number); // in MHz
    // as many rows as there are molecules, as many columns as there are rough points
    // so tAKING THE ROW was a correct thiong to do, just the name used to be bad - column_number
    nn = (double*)calloc(point_numbers, sizeof(double));
    for (int i = 0; i < point_numbers; i++)
    {
        test = gsl_vector_get(temp, i);  //  in MHz
        nn[i] = test;
    }
    
    double bottom_of_parabola = gsl_vector_min(temp); // in MHz
   
    for (int row_counter = 0; row_counter < well_numbers; row_counter++)
    {
        delta = gsl_vector_get(temp, row_counter*2+1) - bottom_of_parabola;
        delta = delta * 1e6 * h_planck; // 2023: converting MHz to Joules
        boltzmann_dist = exp(-delta / (k_b * temperature));
        gsl_matrix_set(boltzmn_distribution, row_counter, row_number, boltzmann_dist); //well, ensemble
        normalization_coefficient = normalization_coefficient + boltzmann_dist;
    }
    free(nn);
    gsl_matrix_get_col(normalization, boltzmn_distribution, row_number); // now normalization contains the desired vector
    gsl_vector_scale(normalization, 1.0 / normalization_coefficient);
    gsl_matrix_set_col(boltzmn_distribution, row_number, normalization);
    nn = (double*)calloc(well_numbers, sizeof(double));
    for (int i = 0; i < well_numbers; i++)
    {
        test = gsl_matrix_get(boltzmn_distribution, i, row_number); // reporting probabilities
        nn[i] = test;
    }
    free(nn);
    gsl_vector_free(temp);
    gsl_vector_free(normalization);
}
