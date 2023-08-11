#include "pch.h"
#include "solve_eigen_system.h"
#include "solve_rate_matrix.h"
#include "Global_variables.h"
#include <iostream>
#include <cmath>
#include <stdio.h>
//#include <gsl/gsl_rng.h>// for producing gaussian distribution
//#include <gsl/gsl_randist.h>// for producing gaussian distribution
//#include <gsl/gsl_linalg.h>// for solving UB=A0 via linear algebra.
//#include <gsl/gsl_complex_math.h>// for complex arithmetic operations
//#include <gsl/gsl_errno.h>
#include "check_situation.h"
#include "solve_ODE.h"
#include "write_on_files.h"// only for debugging

// Modified in 2023 to allow for different format of input rate_matrix
// which now has rows summing up to zero and is in Hz.
// it is transposed and multiplead by 1e-6 below
// the time input still has to be in microseconds
using namespace std;

void SolveRateMatrix(gsl_matrix * rate_matrix, int ensemble, gsl_matrix* rate_file /*for cooling or burning*/,
                     const double evolution_time/*it could be  burn_time or cooling_time, */,
                     gsl_vector * A0, int number_of_row /*either cooling or burning */, double temperature) 
    // temperature seems to be used as some sort of a flag, not as real temperature
    // rate file does not contain rates, it contains probabilities
{
    int column_numbers_counter = ensemble * well_numbers + 1; 
    // for ensemble=ensemble_number it gives columns_number from outside

    int bad_condition = 0;// When the eigensystem method does not work correctly, this number will be something other than zero 
    //(according to different possible situations) and program chooses (slower) ODE solver method instead. 
    //For instance, when the probability in one well is more then one, the loop is broken immediately and ODE solver will be chosen. 
    //Also, the probabilities should sum up to 1, otherwise, ODE solver will be chosen, etc.
    double test, sum = 0;
    double* nn;

    gsl_matrix_transpose(rate_matrix); // introduced in 2023 to ensure compatibility with 
    // the rest of the original Mehdi's logic.
    // Below this point contents of each column should add to zero, see Jackson Biophysics Chapter 9 page 216
    gsl_matrix_scale(rate_matrix, 1e-6); // converting Hz to MHz
    
    nn = (double*)calloc(well_numbers, sizeof(double));
    for (int i = 0; i < well_numbers; i++)
    {
        test = gsl_matrix_get(rate_matrix, 1, i); // reporting rates in Hz, row, column
        nn[i] = test;
    }
    free(nn);

    for (int i = 0; i < well_numbers; i++)
	{
        if (abs(gsl_vector_get(A0,i)) < extreme_small_number)
            gsl_vector_set(A0,i,0);
        sum += gsl_vector_get(A0,i);
    }

    if (abs(sum-1)> accuracy_measure) //accuracy measure is 0.01
	{
        cout << "The sum of probabilities in system " << ensemble << " at temperature " << temperature << " is " << sum << " not 1 from the beginning." << endl;
        throw;
    }
    //
    gsl_vector * A02 = gsl_vector_calloc (A0->size);
	// Mehdi: It is used as a backup for A0 when we want to switch to ODE solver. A0 could be changed during the previous steps.
    gsl_vector_memcpy(A02, A0);
    gsl_matrix * rate_matrix2 = gsl_matrix_calloc (well_numbers, well_numbers);
    gsl_matrix_memcpy(rate_matrix2, rate_matrix);
    //
	//after solving the eigensystem "rate_matrix" will be changed to something completely different, so we need a backup.
    //
    gsl_vector *B = gsl_vector_calloc (well_numbers); 
	// an unkown vector which is equal to U^-1 * A0, however, it is found by solving following equation directly UB = A0.x
    gsl_vector_complex * eigenvalue = gsl_vector_complex_calloc (well_numbers);// eigen values
    gsl_matrix_complex * eigenvector = gsl_matrix_complex_calloc (well_numbers, well_numbers);// eigen vectors.


    double trace=0; // sum of the diagonal elements of the rate matrix
    for (int i = 0; i < well_numbers; i++)
        trace += gsl_matrix_get(rate_matrix, i, i);

    int error=0;

    SolveEignesystem(rate_matrix, eigenvalue, eigenvector, error);// rate_matrix has been changed from this point.

    gsl_matrix *U = gsl_matrix_calloc (well_numbers, well_numbers); // U matrix to keep eigen vectors
    gsl_matrix *UU = gsl_matrix_calloc (well_numbers, well_numbers); 
    // Mehdi: after decompostion for solving UB = A0, U will be changed, so I will make a copy of it temporarly. 
    //It should be solved in another way to make program more efficient.
    //
    gsl_vector * eigenval = gsl_vector_calloc (well_numbers);
    gsl_vector * multiplication = gsl_vector_calloc (well_numbers);

    double imaginary_part=0;
  
    if (error == -1)
	{// This condition should better be checked right after SolveEigensystem, however, 
        //after goto there should not be any variable declarations to avoid getting [-fpermisive] error.
        bad_condition = 1;
        goto ODE_solver;
    }

    for ( int i = 0; i < well_numbers; i++)
    {
        gsl_complex eigenvalue_i = gsl_vector_complex_get (eigenvalue, i); // i-th eigenvalue
        gsl_vector_complex_view eigenvector_i = gsl_matrix_complex_column (eigenvector, i); // i-th eigenvector

        gsl_vector_set(eigenval, i, GSL_REAL(eigenvalue_i)); // making eigenvalues real numbers by force
        
        for ( int j = 0; j < well_numbers; j++)
        {
            gsl_complex z = gsl_vector_complex_get(&eigenvector_i.vector, j);  //j-th element of i-th eigenvector
            imaginary_part = GSL_IMAG(z);
            //cout << GSL_REAL(z) << " + i " << GSL_IMAG(z) << "\n\t";

            if (to_string(abs(GSL_REAL(z))) == "nan"|| to_string(abs(GSL_IMAG (z))) == "nan")
			{// in case that eigenvectors are not calculated properly (any element equals to "nan"), ODE_solver should be used to solve the rate equations.
                bad_condition = 1;
                goto ODE_solver;
            }
            if(abs(imaginary_part)>1e-10) // if imaginary part is too large
			{
                bad_condition = 2;
                goto ODE_solver;
                //throw;
            }

            gsl_matrix_set (U, j, i, GSL_REAL(z)); 
            // when the rate matrix has only real positive elements and the summation of each column is 0, all eigen values are real, 
            //however the routine of gsl is functioning in such a manner that complex numbers are involved with solving a nonsymmetric matrix eigen value problem.
            gsl_matrix_set (UU, j, i, GSL_REAL(z));
        } // end of for j
    } // end of for i
    //
    double summation_over_eigenvalues, rel_trace;
    summation_over_eigenvalues = 0;
    rel_trace=abs(trace * 1e-6);
    for (int i = 0; i < well_numbers; i++)        summation_over_eigenvalues += gsl_vector_get(eigenval, i);
    if (abs(summation_over_eigenvalues - trace) > rel_trace) // made it relative to trace, not absolute 1e-6
    {
        bad_condition = 3;
        goto ODE_solver;
    }

      // to check that the magnetiude of each eigenvector is unit./////////////
    double summation;
    for (int i = 0; i < well_numbers; i++)
    {
        summation = 0;
        for (int k = 0; k < well_numbers; k++)  summation += gsl_matrix_get (UU, k, i) * gsl_matrix_get (UU, k, i);
        //if (abs(summation-1) > 1e-10 && temperature == -100) cout << "The magnetiude of eigenvector of corresponding eigenvalue in ensemble " << ensemble << " should be 1 however it is " << summation << " in well " << i << endl;
    }


    //are eigenvectors correct? I guess we are checking what happens if we substitute them back into thr equation we wanted to solve VZ
    double check_value;
    for (int k = 0; k < well_numbers; k++)
    {
        for (int i = 0; i < well_numbers; i++)
        {
            summation = 0;
            for (int j = 0; j < well_numbers; j++)            summation += gsl_matrix_get(rate_matrix2, i, j) * gsl_matrix_get(UU, j, k);
            //
            gsl_vector_set(multiplication, i, summation);
            check_value = abs(gsl_vector_get(multiplication, i) - gsl_vector_get(eigenval, k) * gsl_matrix_get(UU, i, k));
            if (check_value > 1e-6)
            {// if eigenvectors do not satisfy this condition, it means that eigenvectors are not calculated very well.
                bad_condition = 4;
                goto ODE_solver;
            } // end if
        } //end of for i
    }// end of for k
    gsl_vector_free (multiplication);

    //Is rate_matrix . U = U.Lambda ? (page 225 of Molecular and Cellular Biophysics, as a condition to make sure that U is correct)

    for (int i= 0; i < well_numbers; i++)
        for (int j = 0; j < well_numbers; j++)
		{
            summation = 0;
            for (int k = 0; k < well_numbers; k++)    summation += gsl_matrix_get(rate_matrix2, i, k) * gsl_matrix_get(UU, k, j);

            if (abs(summation-gsl_matrix_get(UU, i, j)*gsl_vector_get (eigenval, j)) > 1e-4)
		    {// in some systems the result of egensystem solver is not applicable which is mentioned as "nan" by gsl. in this case ODE_Solver should be used.
                //cout << "U in ensemble " <<  ensemble << " is not correct"<< " since U is " << gsl_matrix_get(UU, i, j) << " and QU is " << summation << endl;
                bad_condition = 5;
                goto ODE_solver;
            } //end if
        } //end for

    // At this point we are sure that eigenvalues and eigenvectors satisfay the general rule AX=lX 
    //(where l and X are eigenvalue and eigenvector respectively.)


    if (bad_condition==0) // if the above bad things did not happen, we are going to actually do some singular value decomposition
	{

        //singular value decomposition
        gsl_matrix * V = gsl_matrix_calloc (U->size2, U->size2);
        gsl_vector * S = gsl_vector_calloc (U->size2);
        gsl_vector * work = gsl_vector_calloc (U->size2);
        gsl_linalg_SV_decomp (U, V, S, work);
        gsl_linalg_SV_solve (U, V, S, A0, B);

        /* double condition_number;// wikipedia: "...the condition number associated with the linear equation Ax = b gives a bound on how inaccurate the solution x will be after approximation"
        condition_number = abs(gsl_vector_get(S,0)/gsl_vector_get(S,well_numbers-1));
        if (condition_number>1e12)cout << "The condition number for system " << ensemble << " at T = " << temperature << " is\t " << condition_number <<endl;
        */
        gsl_matrix_free (V);
        gsl_vector_free (S);
        gsl_vector_free (work);

        for (int i=0; i < well_numbers; i++)
	    {
            check_value = 0;
            for (int j = 0; j < well_numbers; j++)            check_value += gsl_matrix_get(UU, i, j) * gsl_vector_get(B, j);
            if (abs(check_value-gsl_vector_get(A0,i)) > 1e-3)
		    {
                bad_condition = 6;
                break; // break out of cycle
            } //end if
        } //end for i
    } //end if bad condition==0
    //
    //
    //int time_counter;// this counter keeps the number of time's step in terms of an integer for applying in the array for_write_file.
    //
    if(bad_condition == 0) // here we are actually using time and letting it evolve...
    {
        double A;// differnt rates A[0] is [A], A[1] is [B], and A[2] is [C]
        double Ai = 0;// A in each well at each time
        sum = 0;// it checks that at all the time the total probability is always (close to)1.
        double b,u,expon;
        int A0_counter = 0;
        for (int i = 0; i < well_numbers; i++)
        {
           // time_counter = 1;
            // if initial_time is equal to evolution_time then time_counter should start at 1 not 0.
            //for (double time = initial_time; time <= evolution_time + round_off;  time += time_step){
            //double time = evolution_time;// above loop could be replaced by this line.
            A = 0;

            for ( int k = 0; k < well_numbers; k++)
            {
                gsl_complex eigenvalue_k = gsl_vector_complex_get (eigenvalue,k);
                b = gsl_vector_get (B,k); 
                u = gsl_matrix_get (UU, i, k); 
                expon = exp (GSL_REAL(eigenvalue_k) * evolution_time); // using time
                Ai = b * u * expon;
                A += Ai;
            }

            A = abs(A);// cmath should be included otherwise this function returns zero.
            if (A>1.0001) 
            {// probabilty in each well must be less than 1. If some probability managed to evolve to be larger than one, we should stop and redo it
                bad_condition = 7;
                break; // get out of this loop, to line 280. Do not write anything to rate_file

            }
            // NB:  here we are writing probabilities to (badly named) rate file.
            // NB: it actually contains probabilities!
            
            // if bad_condition is still zero:
            gsl_matrix_set (rate_file, number_of_row, column_numbers_counter, A); // this clearly behaves as if there are many rows in this file and 
            // therefore we have intermediate steps stored here...
            // this makes sense, as for spectra calculations we are using last segment...
            // well_numbers_counter = ensemble*well_numbers+1
            //++time_counter; // time counter can reach number of wells
            //
            //}
            //time counter loop// unneccessary, I have to remove this loop since I am calulating the probability in specific moment.

            if (bad_condition == 0) // 
            {
                /*if (time_counter - 1 != number_of_row)
                {
                    CheckSituation (temperature, situation);
                    cout << "number of row in " << situation << " is " << number_of_row << " and time_counter is " << time_counter << endl;
                    cout << "There is a round-off problem in time loop, please chose another time interval or time step in "<< situation << " parameters" << endl;
                    throw;
                } */
                if (abs(gsl_matrix_get (rate_file, number_of_row, column_numbers_counter)) < extreme_small_number)
                    gsl_matrix_set (rate_file, number_of_row, column_numbers_counter,0);
                // in case that the probability in some well is so small the calculation error could be huge. 
                //Mehdi set the value smaller than some limite (which he found small enough by try-and-error) to zero.

                gsl_vector_set (A0, A0_counter, gsl_matrix_get (rate_file, number_of_row, column_numbers_counter));
                sum += A;
               ++A0_counter;
               ++column_numbers_counter;
           }else break;

        }// i loop (main loop over all wells)
    } // end of if bad condition=0

    if (abs(sum-1) > accuracy_measure)// in case that bad_condition is not 1, the summation over all probabilities should be one.
    // in some systems the result of egensystme solver is not applicable which is mentioned as "nan" by gsl. in this case ODE_Solver should be used.
        bad_condition = 8;




    ODE_solver:
    if (bad_condition != 0)
    {
        CheckSituation (temperature, situation);
        switch (bad_condition)
        {
        case 1:
            cout << "Eigensystem calculation cannot be completed through matrix calculation. The ODE solver is used." << endl;
            break;
        case 2:
            cout << "Eigenvalue is complex (imaginary part is " << imaginary_part << ") in system " << ensemble << " and " << situation << " process. so ODE solver is used." << endl;
            imaginary_part = 0;
            break;
        case 3:
            cout << "Trace is = " << trace << " and sum over eigenvalues is " << summation_over_eigenvalues << endl;
            cout << "Trace of rate matrix is not the same as sum over eigenvalues in system " << ensemble << " and " << situation << " process. Thus ODE solver is used." << endl;
            break;
        case 4:
            cout << "Eigenvectors are not correct in system " << ensemble << " and " << situation << " process. Thus ODE solver is used." << endl;
            break;
        case 5:
            cout << "Q.U is not equal to U.Lambda in system " << ensemble << " and " << situation << " process. Thus ODE solver is used." << endl;
            break;
        case 6:
            cout << "UB is not equal to A0 in system " << ensemble << " and " << situation << " process. Thus ODE solver is used." << endl;
            break;
        case 7:
            cout << "Probability in one well of system " << ensemble << " and " << situation << " process is greater than one, so ODE solver is used." << endl;
            break;
        case 8:
            cout << "The sum over probability in system " << ensemble << " and " << situation << " process is not one, so ODE solver is used." << endl;
            break;
        }

       column_numbers_counter = ensemble * well_numbers + 1; //reset
        // in any case that ODE_solver is called, this counter should be reset.
        gsl_vector_memcpy(A0, A02); // destination, source... copying original incoming set of probabilities from the backup
        
        //
        WriteRateMatrix (rate_matrix2); // writes rate matrix to file rate_matrix.dat
        // rate_matrix2 is the backup of the original matrix and is used in the ODE solver.
        
        SolveODE(A0, evolution_time); // reads file from rate_matrix.dat
        // returns modified A0
        // // I do not think we need to do it this way - writing to file and immediately reading from file. Waste of time.
        // on the other hand, solve_ODE is not very easy to comprehend or modify
        //
        //NB: sometimes for some reason (e.g. sqrt of a negative value,...) one of the elements of A0 
        // could be "nan", which unfortunatelly the program does not complain about and continue solving the problem. 
        // However, at this case, ODE solver cannot solve it. 
        // Try to check all elements of the matrix and/or vectors.
        // //
        //++ODE_counter;
        sum = 0;
        //time_counter = 1;
        for (int i = 0; i < well_numbers; ++i)
	    {
            if (to_string(abs(gsl_vector_get(A0,i))) == "nan")
                cout << "Even ODE solver cannot solve the rate equation in system " << ensemble << " and " << situation << " process." << endl;
            gsl_matrix_set (rate_file, number_of_row, column_numbers_counter, gsl_vector_get(A0,i));
            sum += gsl_vector_get(A0,i);
            ++column_numbers_counter;
	    }
        //if(temperature == -100) cout << "After solving the ODE directly, the probability in system " << ensemble << " at temperature " << temperature << " is " << sum << endl;

    } // end of if bad condition is not 0


    if (abs(sum-1) > accuracy_measure) // checking if sum of new probabilities is close to one
    {
        CheckSituation (temperature, situation);
        if (situation == "cooling" )
        {
            cout << "The probability in system " << ensemble << " and " <<  situation << " process and temperature " << temperature << " is " << sum << " , however, it should be 1" << endl;
            cout << "bad conditon is " << bad_condition << endl;
            throw;
        }    //end if
        else {
            cout << "The probability in system " << ensemble << " and " <<  situation << " process at time " << evolution_time << " is " << sum << " , however, it should be 1" << endl;
            cout << "bad conditon is " << bad_condition << endl;
             throw;
        } //end else

    } //end if


    gsl_vector_complex_free(eigenvalue);
    gsl_matrix_complex_free(eigenvector);
    gsl_matrix_free (U);
    gsl_matrix_free (UU);
    gsl_vector_free(B);
    gsl_matrix_free (rate_matrix2);
    gsl_vector_free (eigenval);
    gsl_vector_free (A02);
}

