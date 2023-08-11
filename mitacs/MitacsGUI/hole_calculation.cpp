#include "pch.h"
#include "hole_calculation.h"
#include "Global_variables.h"
//#include<iostream>
// I believe we no longer call this one - 2022
using namespace std;
 /* void SetHoleProbability (gsl_vector * pre_burn, gsl_vector * post_burn, gsl_matrix * hole_wells, int ensemble)
{
    
    gsl_vector_sub (post_burn, pre_burn);// post_burn = post_burn - pre_burn
    gsl_matrix_set_col (hole_wells, ensemble, post_burn);//the column number of ensemble is equal to post_burn which is now the result of post_burn - pre_burn
    //upper lines are equivalent to the below lines.
    /*double difference;
    for (int i = 0; i < dimension; ++i){
        difference = gsl_vector_get (post_burn, i) - gsl_vector_get (pre_burn, i);// post_burn - pre_burn
        gsl_matrix_set (hole_wells, i, ensemble, difference);
        
    }*/
    



void SetTotalFrequencies (gsl_matrix * wells/*either hole_wells or longtime_distribution*/, gsl_matrix * total_frequency_hole)
{
    // this function convert a table of well_number vs probability of each well in an ensemble to total_frequency

    int total_frequency_counter = 0;
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        // SOMTHING IS WORNG here since sometimes ensemble is greater than ensemble_number but it enters into the loop. 
        //It can be resolved by changing ++ensemble to ensemble++ then one can change it again and it works !!!
        for (int i = 0; i < well_numbers; i++)
        {// X-axis of hole_wells matrix is frequency in GHz
            gsl_matrix_set (total_frequency_hole, total_frequency_counter, 1, gsl_matrix_get(wells, i, ensemble) / ensemble_number);
            // dividing by ensemble_number makes normalized output. The second column (1) is probability or y-axis
            //std::cout << ensemble;
            ++total_frequency_counter;
        }

    }
}
