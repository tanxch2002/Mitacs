#include "pch.h"
#include "solve_ODE.h"
#include "write_on_files.h"
//#include <gsl/gsl_odeiv2.h>
//#include <gsl/gsl_errno.h>
#include "Global_variables.h"
#include <malloc.h>
using namespace std;

void ImportRateMatrix(gsl_matrix * rate_matrix)
{
	FILE * f;
	fopen_s(&f, "rate_matrix.dat", "r");
    gsl_matrix_fscanf(f, rate_matrix); // file to matrix
    fclose (f);
}


int func (double , const double y[], double f[], void *)// it seems that there is no need to define t and params, so they are removed.
{
    gsl_matrix * rate_matrix = gsl_matrix_calloc (well_numbers, well_numbers);
    ImportRateMatrix(rate_matrix);
    for (int i = 0; i < well_numbers; i++)
	{
        f[i]=0;
        for (int j = 0; j < well_numbers; j++)
		{
            //cout << "rate_matrix : " << gsl_matrix_get(rate_matrix, i, j) << " Y[ " << j << " ] = " << y[j] << endl;
            f[i]+= gsl_matrix_get(rate_matrix, i, j) * y[j];
        }
    }
    gsl_matrix_free (rate_matrix);
    return GSL_SUCCESS;
}


int jac (double , const double* , double *dfdy, double dfdt[], void*)// Jacobian something...
{
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, well_numbers, well_numbers);// dfdy_mat is the Jacobian matrix
    gsl_matrix * m = &dfdy_mat.matrix;

    ImportRateMatrix(m);

    for (int i = 0; i < well_numbers; i++) dfdt[i]=0;

    return GSL_SUCCESS;
}
//#include <cmath>
//
//
void SolveODE(gsl_vector* A0, double simulation_time) 
{
    int* null_pointer;
    null_pointer = 0;
    double* y;
    y = (double*)calloc(well_numbers, sizeof(double));
    for (int i = 0; i < well_numbers; i++)
    {
        y[i] = gsl_vector_get(A0, i); // converts gsl array to regular array
    }
    gsl_odeiv2_system sys = { func, jac, well_numbers, &null_pointer }; // jac contains import rate matrix

    gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_bsimp, 1e-10, 1e-10, 0.0);
    //bsimp: Implicit Bulirsch-Stoer method or Bader and Deuflhard (slower than msbdf BUT it work better).
    //msbdf: A variable-coefficient linear multistep backward differentiation formula (BDF) method in Nordsieck form.
    
    double t = 0.0, t1 = simulation_time;
   
    //double y[well_numbers];
    
    double ti = t1 ;
    int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

    if (status != GSL_SUCCESS)
    {
        printf ("error, return value=%d\n", status);
        throw;
    }

    //printf ("%.5e \t%.5e \t %.5e \t %.5e \t%.5e\n", t, y[0], y[1], y[2], y[10]);
    double test;
    for (int i = 0; i < well_numbers; i++)
    {
        test = y[i];
        test=std::abs(test);

        if(std::abs(y[i]) < extreme_small_number) y[i] = 0;
            // in case that the probability in resonant well is so small the calculation error could be huge. 
            //Mehdi set the value smaller than some limit (which he found small enough by try-and-error) 
            //to zero.abs() returns zero even when the cmath is included.
         
        gsl_vector_set (A0, i, y[i]);
    }
    gsl_odeiv2_driver_free (d);

}
