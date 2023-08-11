//#include "stdafx.h"
#include "pch.h"
#include "set_output_matrix_at_different_T.h"
#include "Global_variables.h"
using namespace std;
void SetOutputAtDifferentTemp (gsl_matrix* output_at_different_T, gsl_matrix* rate_file, int& temperature_counter, int ensemble,
                               const int number_of_row /*either cooling or burning */)
{
//looks like this thing saves the matrix containing all rates 
    int well_numbers_counter = ensemble * well_numbers + 1; // ensemble is variable

    int A0_counter = 0;
    for (int i = 0; i < well_numbers; i++)
    {
        gsl_matrix_set (output_at_different_T, A0_counter, temperature_counter, gsl_matrix_get (rate_file, number_of_row, well_numbers_counter));
        // matrix to be modified, index, index, value to put there
        //so it gets something from rate_file, with rates at different T, and get the first element of some specific column... Not sure why it does it.
        
        // Mehdi used A0_counter to avoid defining another similar counter to count row numbers.

        ++A0_counter;
        ++well_numbers_counter;
    }

    ++temperature_counter;// it keeps the column number of output_at_different_T.




}
