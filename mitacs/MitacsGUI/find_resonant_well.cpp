#include "pch.h"
#include "find_resonant_well.h"
#include "order_well_vs_frequency.h"
#include <gsl/gsl_sort_vector.h>
#include <iostream>
#include "Global_variables.h"

using namespace std;
//This function finds the resonant wells in the course of a SMS scan 
// and the time delay which is needed to reach them according to the scan speed
//this is used in single-molecule modeling
void FindResonantWellandTimeDelay (gsl_matrix * en_difference,  gsl_matrix * time_delay_in_frequency_scan, int ensemble)
{
    gsl_matrix * sorted_frequency_vs_well = gsl_matrix_calloc(well_numbers, 2); // rows, columns
	// one column for sorted frequency and another one for number of well for each system.
	//
    OrderFrequenciesVsWell (en_difference, sorted_frequency_vs_well, ensemble);
    // returns array where frequencies are sorted but wells are not
    double frequency = minimum_frequency; // starting frequency of scan
    int number_of_frequncy;
    //this is the number of steps to reach to resonant frequency while scanning.
    double deltaF, deltaTrecovery;
    int resonant_well_number;
    double closest_frequency_to_resonant_well;
	//
    for (int i = 0; i < well_numbers; i++) //dimension is number of wells+1
		//looks like this goes by index of sorted frequency column, which is not well number
	{
        deltaF = gsl_matrix_get(sorted_frequency_vs_well, i, 0) - frequency; 
        //wells sorted by frequency (0-th column), well numbers may be in random order...
	    //since we start at minimal frequency, first deltaF is going to be positive
	    // this is how far we 
        if (deltaF < 0)
	    {
            cout << "The frequency step is so big that it jumps over well " << i << ", please reduce it" << std::endl;
		//correct warning should be that we are past the resonant well and moving away from it, not that stepis too large...
            throw;
        }  //end if
	    resonant_well_number = (int)gsl_matrix_get (sorted_frequency_vs_well, i, 1); //resonant well number 
	    //deltaF is how far we are from a given frequency on a frequency scale
	    //resonant well number is the respective well number (not i)



	    //setting new matrix
        gsl_matrix_set (time_delay_in_frequency_scan, i, 0, resonant_well_number);
	    //time_delay_in_frequency_scan is a matrix containing the list of delays
	    //filling one column of the new matrix with well numbers, which will be out of order with respect to i

        number_of_frequncy = (int)(deltaF/frequency_step);
	    // the number of steps which it will take to reach the resonant well.
	
        deltaTrecovery = number_of_frequncy * time_interval_between_each_frequency_step;

        gsl_matrix_set (time_delay_in_frequency_scan, i, 1, deltaTrecovery);
	    // second coulum: this is how much time it will take to get there from current frequency


        closest_frequency_to_resonant_well =  (int)gsl_matrix_get (sorted_frequency_vs_well, i, 0);
	    
        gsl_matrix_set (time_delay_in_frequency_scan, i, 2, closest_frequency_to_resonant_well);

        frequency += number_of_frequncy * frequency_step; //frequency is reset to frequency of the current well
	    //frequencies are ordered, the well numbers are not

    } //end for

 //gsl_vector_free (resonant_frequency);
 gsl_matrix_free (sorted_frequency_vs_well);

}
