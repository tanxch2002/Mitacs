// in this function the energy difference (in terms of frequency) between well is ordered in ascending way (reminder: the order of well is not the same as order of frequency, meaning that well#0 could have greater energy than well#3 and so on)
#include "pch.h"
#include "order_well_vs_frequency.h"
#include <gsl/gsl_sort_vector.h>
#include "Global_variables.h"

using namespace std;

void OrderFrequenciesVsWell (gsl_matrix * ene_difference, gsl_matrix * sorted_frequency_vs_well, int ensemble)
{
    gsl_vector * sorted_frequency = gsl_vector_calloc(well_numbers);
    gsl_matrix_get_col (sorted_frequency, ene_difference, ensemble);
    // energy difference between each well in corresponding system is stored in sort.
    gsl_sort_vector (sorted_frequency);
    //"This function sorts the elements of the vector v into ascending numerical order."

    for (int i = 0 ; i < ensemble_number; i++)
	{
        double frequency;
        int j=0;// it scans the value in the frequency column from the first value till the last one.
        int k=0;// it scans sorted frequencies in vector "sort".
        while (k < well_numbers )
		{// This loops continues until all frequencies are sorted by comparing with vector "sort".

            if (gsl_matrix_get(ene_difference, j, ensemble) == gsl_vector_get(sorted_frequency,k))
			{
                frequency = gsl_matrix_get(ene_difference, j, ensemble);
                gsl_matrix_set(sorted_frequency_vs_well, k, 0, frequency);// this is the frequency corresponding to resonant well.
                gsl_matrix_set(sorted_frequency_vs_well, k, 1, j);// this is the resonant well number.
                j=-1;// All unsorted values should be compared with sorted one, so each time comparison will be reset. Since j should be 0 for the next k, It should be set to -1 here, since after satisfying of this if-condition one unit will be added to j (++j). NB:  It should be mentioned that this is not the optimum way to compare values because for ordered values there should not be any comparison again. This repetition is time consuming but I did not find easier/faster way to do that at this moment.
                ++k;
            } //end if
            ++j;
        } //end while
    } //end for

    gsl_vector_free (sorted_frequency);

} //end sub
