#include "pch.h"
#include "check_frequency_range.h"
#include <iostream>
#include "Global_variables.h"
using namespace std;
void CheckFrequencyRange (gsl_matrix * en_difference) // modified in 2022 to operate with new energy_difference matrix
{
    //gsl_matrix* energy_difference = gsl_matrix_calloc(well_numbers, ensemble_number);
    gsl_vector * optimum = gsl_vector_calloc (well_numbers*ensemble_number);
	//optimum is a vector of the size well_numbers; looks like it contains frequencies
    int a;
    for (int i = 0; i < ensemble_number; i++) // 5000 --> 0....4999
    {
        for (int j = 0; j < well_numbers; j++) // 20 --> 0...19
        {
            a = gsl_matrix_get(en_difference, j, i) /(h_planck*1e6); // converting energy_differences to MHz
            gsl_vector_set(optimum, i * well_numbers +j, a); //4999*20+19=99999
            // i=0: 0...19
            // i=1: 20..39, 
        }
    }

    double min, max;
    int test=0;
    gsl_vector_minmax(optimum, &min, &max);

    if (minimum_frequency > min){
        cout << "Please decrease the minimum frequency in frequency scan part since the minimum energy in the absorption band is smaller than chosen value" << endl;
        test = -1;
    }
    if (maximum_frequency < max){
        cout << "Please increase the maximum frequency in frequency scan part since the maximum energy in the absorption band is bigger than chosen value" << endl;
        test = -1;
    }

    if (test == -1)
        throw;

    gsl_vector_free (optimum);
}
