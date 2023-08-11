#include "pch.h"
#include "MitacsGUI.h"
#include "MitacsGUIDlg.h"
#include "create_rate_matrix.h"
#include <cmath>
#include <iostream>
#include "write_on_files.h"// only for debugging
#include <conio.h>
#include <iomanip>
#include "Global_variables.h"



//using namespace std;
// 
// Technically, this 2023 incarnation of this subroutine produces not the matrix handled by solve_rate_matrix,
// but so called GENERATOR MATRIX where sum of each row is 0, and not sum of each column
// Also the output (matrix) contains rates in Hz, not MHz.
// Solve rate matrix takes care of the changes necessary to convert to mehdi's input
// it still uses time in microseconds for input
// 
//As you can see, we are passing pointers from one function to another, not variables themselves
//this allows for bi-directional flow of date in and out 

void CreateRateMatrix(gsl_vector* coeff_matrix, gsl_matrix* eigenvalues,
    gsl_vector * attempt_frequency, gsl_vector* TE, /*transition probabilities, left to right*/
    gsl_matrix * parabola_for_rate /*6; low-res energy landscape, for either in the ground or excited state, size is point numbers*/,
                      gsl_matrix * rate_matrix, gsl_matrix * temp_matrix/* 8*/,
                      double temperature  /*9*/, int ensemble, gsl_vector * v_barrier_height /*either in the ground or excited state*/)
    // the latter quantity is actually not used inside this subroutine in this version. 
{                                             
    double thermal_average_numerator, thermal_average_denominator, thermal_average_numerator_RtoL, thermal_average_denominator_RtoL;
    double inverse_temp = 1 / (k_b * temperature);
    double Energy, Joule_E, expo, prob_n, prob_n2, expo_phonon, n_phonon, delta_E; // I added the definition of prob_n2 here - Nov 2022
    //extern int well_numbers, Hamiltonian_well_numbers;
    int index, index2, index3, index4; // starting state, tunneling first, phonon first
    double current_phonon, current_T, current_freq;
    gsl_vector* Thermal_average_vector_LtoR=gsl_vector_calloc(well_numbers); //for tunneling
    gsl_vector* Thermal_average_vector_RtoL = gsl_vector_calloc(well_numbers);
    double* mm;
    double test, a=0;
    double T; //transition probability for barrier-hopping

    // since phonon population numbers depend only on the energy difference between two states, one can calculate them once for each pair of states
    // no need to recalculate them for each well

    //
    gsl_matrix* phonon_numbers = gsl_matrix_calloc(Hamiltonian_well_numbers, Hamiltonian_well_numbers); 
    // here we store all phonon-related factors, 
    //that are different for uphill and downhill processes
    // these are all temperature-dependent, so we have to recalculate them
    for (int j = 0; j < Hamiltonian_well_numbers; j++) //cycle over energy levels, initial states , phonons included
    {
        //thermal_average_numerator = 0;
        for (int k = 0; k < Hamiltonian_well_numbers; k++) //cycle over final states, phonons included 
        {
                if (k < j) // downhill tunneling
                {
                    delta_E = (gsl_matrix_get(eigenvalues, ensemble, j) - gsl_matrix_get(eigenvalues, ensemble, k)) ; 
                    //eignvals are in joules, positive asymmetry
                    expo_phonon = pow(2.718281828459, (inverse_temp * delta_E));
                    n_phonon = 1 / (expo_phonon - 1); //phonon population number
                    n_phonon = n_phonon + 1; //downhill tunneling involves giving energy to a phonon
                }
                if (k == j) n_phonon = 0;
                // forbidding non-phonon-assisted processes... if delta=0, then expo_phonon=1 and we would get division by zero
                // alternatively, if (k == j) thermal_average_numerator = 1; allowing no-phonon-assisted processes, need to make a checkbox on GUI...
                if (k > j) //uphill tunneling
                {
                    delta_E = (gsl_matrix_get(eigenvalues, ensemble, k) - gsl_matrix_get(eigenvalues, ensemble, j)); 
                    //eigenvals in joules, positive
                    expo_phonon = pow(2.718281828459, (inverse_temp * delta_E));
                    n_phonon = 1 / (expo_phonon - 1);
                    //thermal_average_numerator = n_phonon; //uphill tunneling involves taking energy from a phonon, so the phonon need to have that enery in the first place
                }
                gsl_matrix_set(phonon_numbers, j,k, n_phonon); //first index initial state, second index final state
         }// end of cycle over k (final states)
    } //end of the cycle over initial states j
    //
    // //
    // //
    // Calculating MODIFIED Garashchuk's stuff, thermally avearged tunneling rates (not probabilities), i.e. including the attempt frequencies
    // phonon scattering properly included as of November 2022. 
    // also included is the possiblity that for some final states in some wells particle with respective energy would be impossible

    //left to right
    for (int i = 0; i < well_numbers; i++) //if well_numbers=22, there are only 21 from which one can go left to right
    { //cycle over wells
        thermal_average_numerator = 0;
        thermal_average_denominator = 0;
        
        for (int j = 0; j < Hamiltonian_well_numbers; j++) //cycle over initial energy levels, with phonon scattering included
        {
            Energy = gsl_matrix_get(eigenvalues, ensemble, j); 
            // current / initial stationary state energy, in  MHz
            Joule_E = Energy ;
            expo = pow(2.718281828459, (-1 * inverse_temp * Joule_E)); // for higher levels everything is thermally forbidden even at 300K
            //
            index = to1D(ensemble, i, j); // initial well, initial state
            prob_n = gsl_vector_get(coeff_matrix, index); // rho for the initial state
            //probability of finding the particle in the initial state with particular energy in a particular well.

            for (int l = 0; l < Hamiltonian_well_numbers; l++) // final states
            {//left to right
                if (i == well_numbers - 1) 
                {
                    thermal_average_numerator = 0; // no movement further to the right from right-most well
                    thermal_average_denominator = 1; // to prevent division by zero
                }
                else {
                    index2 = to1D(ensemble, (i+1), l); // final state in product well, Molecules, wells, levels
                    prob_n2 = gsl_vector_get(coeff_matrix, index2); // rho for the final state in product well
                    index3 = to1D(ensemble, i, l); // for the case of phonon first
                    current_T = (gsl_vector_get(TE, index)+ gsl_vector_get(TE, index3))/2;
                    // average of tunneling first and phonon first
                    current_freq = gsl_vector_get(attempt_frequency, index);
                    current_phonon = gsl_matrix_get(phonon_numbers, j, l);
                //Calculating Thermally averaged tunneling RATES for each well (averaging the product of transmission PROBABILIY with attempt frequency and phonon factors). 
                     thermal_average_numerator += current_T * current_freq *current_phonon  * expo * prob_n * prob_n2;
                     thermal_average_denominator += prob_n * expo * prob_n2; // thermal and QM averaging for initial state only, but we still multiply by prob_n2 
                     //to prevent the denominator from getting unrealistically large.
                     // another possibility is to only sum terms for which the final well / state rho is large enough
                }
            
            } //end of the cycle over final states
        } //end of the cycle over initial states
        if (thermal_average_denominator > 0)
        {
            gsl_vector_set(Thermal_average_vector_LtoR, i, (thermal_average_numerator / thermal_average_denominator));
        }
        
    } //end of cycle over wells

    mm = (double*)calloc(well_numbers, sizeof(double));
    for (int i = 0; i < well_numbers; i++)
    {
        test = gsl_vector_get(Thermal_average_vector_LtoR, i); 
        // reporting rates in Hz, on the order of 4-8 x 10^9 at 300K and single Hz at 5 K
        mm[i] = test;
    }
       // and now right to left...
    for (int i = 0; i < well_numbers; i++) //if well_numbers=22, there are only 21 from which one can go left to right
    { //cycle over wells
        thermal_average_numerator_RtoL = 0;
        thermal_average_denominator_RtoL = 0;
        for (int j = 0; j < Hamiltonian_well_numbers; j++) //cycle over initial energy levels, with phonon scattering included
        {
            Energy = gsl_matrix_get(eigenvalues, ensemble, j);
            // current / initial stationary state energy, in  MHz
            Joule_E = Energy;
            expo = pow(2.718281828459, (-1 * inverse_temp * Joule_E)); // for higher levels everything is thermally forbidden even at 300K
            //
            index = to1D(ensemble, i, j); // initial well, initial state
            prob_n = gsl_vector_get(coeff_matrix, index); // rho for the initial state
            //probability of finding the particle in the initial state with particular energy in a particular well.

            for (int l = 0; l < Hamiltonian_well_numbers; l++) // final states
            {//left to right
                
                if (i == 0)
                {
                    thermal_average_numerator_RtoL = 0; // cannot go right to left from the left-most well
                    thermal_average_denominator_RtoL = 1; // to prevent division by zero
                }
                else
                {
                    index2 = to1D(ensemble, (i - 1), l); // final state in product well
                    prob_n2 = gsl_vector_get(coeff_matrix, index2); // rho for the final state in product well
                    //
                    index3 = to1D(ensemble, i-1, l); // FINAL STATE, barrier between i and i-1 - for phonon scattering first
                    // Taking into account that T(E)-s have zero for last well (cannot go further right
                    //therefore for last barrier one needs to take pre-last well
                    //  
                    
                    index4= to1D(ensemble, i - 1, j); // Initial state, barrier between i and i-1
                        current_T = (gsl_vector_get(TE, index4) + gsl_vector_get(TE, index3)) / 2;
                    // average of tunneling first and phonon first
                    current_freq = gsl_vector_get(attempt_frequency, index); // for starting well
                    current_phonon = gsl_matrix_get(phonon_numbers, j, l);
                    thermal_average_numerator_RtoL += current_T * current_freq * current_phonon * expo * prob_n * prob_n2;
                    thermal_average_denominator_RtoL += prob_n * expo * prob_n2;
                }
            } //end of the cycle over final states
           
        } //end of the cycle over initial states
        
        if (thermal_average_denominator_RtoL > 0)
        {
            gsl_vector_set(Thermal_average_vector_RtoL, i, (thermal_average_numerator_RtoL / thermal_average_denominator_RtoL));
            // 1 over 2pi removed from there on the grounds that T(E) are dimensionless and we know the units of the attempt frequency - Hz
        }
        a = a + 1;
    } //end of cycle over wells

    for (int i = 0; i < well_numbers; i++)
    {
        test = gsl_vector_get(Thermal_average_vector_RtoL, i); // reporting rates in Hz, on the order of 4-7 x 10^9 at 300K
        mm[i] = test;
    }
 //NB: there is a difference between tunneling first then phonon scattering and
// phonon scattering first, then tunneling. 
// Should not matter much, as either the phonon factor(for uphill tunneling)  or thermal population factor (for downhill tunneling) 
//will kill any terms for too large difference of energy levels, and for close levels (from the same group) order should not matter much as T(E) is not that different.
    // It does matter in a situation when current energy is impossible in the product well because bottom is higher.
    // then T(E) for initial energy is set top zero, but T(E) for energy after adding vibrational energy may not be zero
    // Taking average of T(E) for two paths
    // 
    // 
    // 2021：Calculating the barrier hopping rates. Most of the calculation is the same as above, but with rates determined by Bolzmann factors. No phonon populations or final state probabilities
    double energy_temp;
    double barrier_top_temp;
    //double barrier_hopping_temp;
    gsl_vector* Thermal_barrier_hopping_LtoR = gsl_vector_calloc(well_numbers); // 2021: Thermal average of barrier hopping
    gsl_vector* Thermal_barrier_hopping_RtoL = gsl_vector_calloc(well_numbers); // 2021: Thermal average of barrier hopping
    gsl_matrix* barrier_hopping = gsl_matrix_calloc((well_numbers - 1), Hamiltonian_well_numbers); 
    // 2021:Rates for classical thermally activated process (barrier - hopping rates)
    //
    //left to right
    for (int i = 0; i < well_numbers; i++) // cycle over wells
    {
        thermal_average_numerator = 0;
        thermal_average_denominator = 0;
        for (int j = 0; j < Hamiltonian_well_numbers; j++) 
            //cycle over initial eigenstates, final state is irrelevant, only reaching the top of the barrier is relevant
        {
            energy_temp = gsl_matrix_get(eigenvalues, ensemble, j) ; // in Joules
            expo = pow(2.718281828459, (-1 * inverse_temp * energy_temp)); // prob to be in initial state, thermodynamic
            index = to1D(ensemble, i, j);
            prob_n = gsl_vector_get(coeff_matrix, index); // prob to be in initial state, quantum
            //
            if (i == (well_numbers - 1))  //left to right
            {
                thermal_average_numerator = 0; //cannot go right from the right-most well
                thermal_average_denominator = 1;
            }
            else {
                barrier_top_temp = gsl_matrix_get(parabola_for_rate, ensemble, (2 * i + 2)); // in MHz
                //parabola for rate has point_numbers points
                delta_E = (barrier_top_temp *1e6*h_planck) - energy_temp; 
                // to top of the barrier from the current state, in joules
                if (delta_E > 0) T = pow(2.718281828459, ((-delta_E) * inverse_temp)); // hopping probability
                else T = 1; // if energy above the barrier
                current_freq = gsl_vector_get(attempt_frequency, index); // for starting well
                thermal_average_numerator += current_freq * T * prob_n * expo; //Gar, 12, numerator
                thermal_average_denominator += prob_n * expo; //Gar Eq. 12 denominatior
            }
        }
        if (thermal_average_denominator > 0) gsl_vector_set(Thermal_barrier_hopping_LtoR, i, (thermal_average_numerator / thermal_average_denominator));
        else gsl_vector_set(Thermal_barrier_hopping_LtoR, i, 0);
    } //end cycle over i (wells)
    //
    for (int i = 0; i < well_numbers; i++)
    {
        test = gsl_vector_get(Thermal_barrier_hopping_LtoR, i); // reporting rates in Hz 7-11x10^10 at 300 K
        mm[i] = test;
    }
    
    // 
    // right to left
    for (int i = 0; i < well_numbers; i++) // cycle over wells
    {
        thermal_average_numerator_RtoL = 0;
        thermal_average_denominator_RtoL = 0;

        for (int j = 0; j < Hamiltonian_well_numbers; j++)
            //cycle over initial eigenstates, final state is irrelevant, only reaching the top of the barrier is relevant
        {
            energy_temp = gsl_matrix_get(eigenvalues, ensemble, j); // in Joules
            expo = pow(2.718281828459, (-1 * inverse_temp * energy_temp));
            index = to1D(ensemble, i, j);
            prob_n = gsl_vector_get(coeff_matrix, index); // reactant well, initial state
            //
            if (i == 0) //first well
            {
                thermal_average_numerator_RtoL = 0; //cannot go left from the left-most well
                thermal_average_denominator_RtoL = 1;
            }
            else {
                barrier_top_temp = gsl_matrix_get(parabola_for_rate, ensemble, (2 * i)); //parabola for rate has point_numbers points
                // first barrier for second well, etc
                delta_E = (barrier_top_temp * 1e6 * h_planck) - energy_temp; // in joules
                if (delta_E > 0) T = pow(2.718281828459, ((-delta_E) * inverse_temp));
                else T = 1;
                current_freq = gsl_vector_get(attempt_frequency, index); // for starting well
                thermal_average_numerator_RtoL += current_freq * T * prob_n * expo; //Gar, 12, numerator
                thermal_average_denominator_RtoL += prob_n * expo;
                // probably do not need to modify this part, as it already contains prob_n that has already been modified to be based on rho's (projections)
            }
        }
        if (thermal_average_denominator_RtoL > 0) gsl_vector_set(Thermal_barrier_hopping_RtoL, i, (thermal_average_numerator_RtoL / thermal_average_denominator_RtoL));
        else gsl_vector_set(Thermal_barrier_hopping_RtoL, i, 0);
    } //end cycle over i (wells)
    //
    // now we have contributions to diagonals parallel to the main diagonal. 
    //Need to put them properly into a matrix.
    // first index=from, second index=to

    
    for (int i = 0; i < well_numbers; i++)
    {
        test = gsl_vector_get(Thermal_barrier_hopping_RtoL, i); // reporting rates in Hz 7-11x10^10 at 300 K
        mm[i] = test;
    }
    for (int i = 0; i < well_numbers-1; i++) //i is number of barrier, there is well_numbers-1 barriers
    {
        //off-diagonal elements
        gsl_matrix_set(rate_matrix, i, (i + 1), gsl_vector_get(Thermal_barrier_hopping_LtoR, i) + gsl_vector_get(Thermal_average_vector_LtoR, i));
        // 0:1... 1:2...2:3............18:19
        // zeroth rates are for the first inner barrier
            
        gsl_matrix_set(rate_matrix, (i + 1), i, gsl_vector_get(Thermal_barrier_hopping_RtoL, i+1) + gsl_vector_get(Thermal_average_vector_RtoL, i+1));
        // 1:0...2:1... 3:2... ... 19:18
        //
    }
    //
    for (int i = 0; i < well_numbers; i++)
    {
        test = gsl_matrix_get(rate_matrix, 1, i); // reporting rates in Hz, row, column
        mm[i] = test;
    }
    
    //
    // //
    // diagonal elements of the rate matrix, rates of staying in the same well
    gsl_matrix_set(rate_matrix, 0, 0, -gsl_matrix_get(rate_matrix, 0, 1)); //NB: changed on Feb 8, 2023
    gsl_matrix_set(rate_matrix, (well_numbers - 1), (well_numbers - 1), -gsl_matrix_get(rate_matrix, (well_numbers - 1), (well_numbers - 2)));

    for (int i = 1; i < well_numbers - 1; i++)
        gsl_matrix_set(rate_matrix, i, i, -(gsl_matrix_get(rate_matrix, i, (i - 1)) + gsl_matrix_get(rate_matrix, i, (i+1))));

    
    for (int i = 0; i < well_numbers; i++)
    {
        test = gsl_matrix_get(rate_matrix, 1, i); // reporting rates in Hz, row, column
        mm[i] = test;
    }
    
    for (int i = 0; i < well_numbers; i++)
    {
        test = gsl_matrix_get(rate_matrix, 2, i); // reporting rates in Hz, row, column
        mm[i] = test;
    }
    free(mm);

    gsl_matrix_memcpy(temp_matrix, rate_matrix);
    //
    gsl_matrix_free(phonon_numbers);
    gsl_vector_free(Thermal_average_vector_LtoR);
    gsl_vector_free(Thermal_average_vector_RtoL);
    gsl_vector_free(Thermal_barrier_hopping_LtoR);
    gsl_vector_free(Thermal_barrier_hopping_RtoL);
 }
 // Technically, this 2023 incarnation of this subroutine produces not the matrix handled by solve_rate_matrix,
// but so called GENERATOR MATRIX where the sum of each row is 0, and not the sum of each column
// The output contains rates in Hz, not MHz.
// Solve rate matrix takes care of the changes
// it still uses time in microseconds for input

