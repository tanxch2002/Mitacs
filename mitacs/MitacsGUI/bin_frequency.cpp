#include "pch.h"
#include "bin_frequency.h"
#include "Global_variables.h"
#include <gsl/gsl_sort_vector.h>// for finding the minimum and maximum frequency in binning process.
#include <iostream>
using namespace std;

void BinFrequency (gsl_matrix * total_frequency, gsl_vector * bin, double& bin_step, double& smallest_frequency)
// this thing uses as input the total_frequency file that contains frequencies and probabilities.
// it calculates minimal frequency, maximal frequency and the bin step
// AND actually bins the probabilites and puts the sums of probabilities into the vector bin. 
// the frequency information apparently gets recreated elsewhere from smallest_frequency and bin step
//
{// as we need the value of bin_step and smallest_frequency in the rest of main program, M passed them by reference
    //finding the minimum and maximum frequency for binnin purposes.
    gsl_vector * optimum_frequencies = gsl_vector_calloc (all_frequency);
    // for finding the maximum and minimum value of frequency among the whole of frequencies.
    gsl_matrix_get_col (optimum_frequencies, total_frequency,0);
    // filling the vector optimum_frequencies with the first column of total_frequency which is the all possible frequencies in the ensemble.
    gsl_sort_vector (optimum_frequencies); 
    //it sorts the optimum_friquencies vector then the first and second elements of this vector is used in bin_step (below)
    
    double largest_frequency = gsl_vector_max(optimum_frequencies); 
    // TO DO: as optimum_friquencies vector has been sorted, it may use another command to read the last element of this sorted vector as a largest one to make the program faster.
    smallest_frequency = gsl_vector_min(optimum_frequencies);
    //cout << "The maximum frequency is: " << largest_frequency << endl;
    //cout << "The minimum frequency is: " << smallest_frequency << endl;

    // binning frequency of all ensemble

    //double bin_step = abs((gsl_vector_get(optimum_frequencies,1)-gsl_vector_get(optimum_frequencies,0)));// the two smallest frequencies between all of frequencies in the ensemble.
    //int total_bins = (largest_frequency - smallest_frequency) / bin_step;// the denominator was chosen in this way since I belive most of the ensemble has more or less the same step between two first frequencies in the begining.

    //***** changing the binning way by rewriting two above lines as below*****//

     bin_step = (largest_frequency - smallest_frequency) / total_bins;
    //bin_step = 10547.849084372607;//WARNING: M changed this to have the same bin_step for both energy landscapes (parabolic and flat models)

    //cout << "bin step= "<< bin_step << "\ttotal bins: "<< total_bins<< endl;

    double summation_bin;

    for (int i = 0; i <= total_bins; i++)
    {// the less than or equal to sign applied here to account the maximum frequency.
        summation_bin = 0;
        for (int j = 0; j < all_frequency; j++)
        {
            if (gsl_matrix_get(total_frequency, j, 0) >= smallest_frequency + i * bin_step &&
                gsl_matrix_get(total_frequency, j, 0) < smallest_frequency + (i + 1) * bin_step)
            {
                summation_bin += gsl_matrix_get(total_frequency, j, 1);
                gsl_vector_set(bin, i, summation_bin);
            }
        }
    }


    gsl_vector_free (optimum_frequencies);
}
