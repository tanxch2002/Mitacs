#include "pch.h"
#include "convert_well_to_frequency.h"
#include "Global_variables.h"
//#include "write_on_files.h"--- not writing anything from inside this subroutine

using namespace std;
void ConvertWellToFrequency (gsl_matrix * energy_difference, gsl_matrix * rate_file, gsl_matrix * last_rates, gsl_matrix * total_frequency, int number_of_row)
//energy difference is a matrix that has well-number rows and ensemble_number columns
// we can send into it the pre-burn type of array instead of rate_file...
// pre_burn = gsl_matrix_calloc(ensemble_number, well_numbers);
{
    for (int row_counter = 0; row_counter < well_numbers; row_counter++)
    {
        for (int column_counter = 0; column_counter < ensemble_number; column_counter++)
        {
            gsl_matrix_set(last_rates, row_counter, column_counter, gsl_matrix_get(rate_file, number_of_row, (well_numbers * column_counter + row_counter + 1)));
            // looks like last_rates is just transposed pre-burn...
        }
    }
    //
    // the first column of rate_file is time.
    
    // The different frequency for each well is calculated and then the frequencies and corresponding probabilities should be sorted in ascending way.
    gsl_matrix * output = gsl_matrix_calloc (well_numbers, (2*ensemble_number));// one column is frequency (energy difference between excited and ground states for each well) and the other one is for probability of corresponiding ensemble.
    //NOTE I think output matrix is useless in this version.
    int total_frequency_counter = 0;
    for (int j = 0; j < ensemble_number; j++)
    {
        for (int i = 0; i < well_numbers; i++)
        {// X-axis of output matrix is frequency in GHz or MHZ - need to check VZ
            gsl_matrix_set (output, i, (2 * j), gsl_matrix_get(energy_difference, i, j));
            // filling the even columns of output by frequencies for each ensembles: 0 (frequency-axis for the first ensemble), 2 (frequency-axis for the second ensemble), 4 (frequency-axis for the third ensemble), ...
            gsl_matrix_set (output, i, (2*j + 1), gsl_matrix_get (last_rates, i, j));
            // filling the odd columns of output (1, 3, 5, ...) by probabilities of wells for each ensembles.
            gsl_matrix_set (total_frequency, total_frequency_counter, 0, gsl_matrix_get(energy_difference, i, j));
            // this matrix will be used to keep all frequencies of all molecules. The first column (0) is frequency (GHz) or X-axis. total_frequency_counter can be replaced by i + ensemle * dimension but I kept it for the sake of clarity.
            gsl_matrix_set (total_frequency, total_frequency_counter, 1, gsl_matrix_get(output, i, (2*j + 1)) / ensemble_number);
            // dividing by ensemble_number makes normalized output. The second column (1) is probability or y-axis
            ++total_frequency_counter;
        }

    }
    gsl_matrix_free (output);

}
