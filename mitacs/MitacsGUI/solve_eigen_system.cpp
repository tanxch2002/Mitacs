#include "pch.h"
#include "solve_eigen_system.h"
//#include <gsl/gsl_math.h>// for eigensystem
//#include <gsl/gsl_eigen.h>// for eigensystem
//#include <gsl/gsl_complex.h>// for eigensystem
//#include <gsl/gsl_complex_math.h>// for eigensystem
//#include <gsl/gsl_errno.h>
#include "Global_variables.h"
#include "write_on_files.h"// only for debugging
#include <iostream>
//#include <string>
using namespace std;
void SolveEignesystem (gsl_matrix* a, gsl_vector_complex *eval, gsl_matrix_complex *evec, int &error)
{
    
    WriteRateMatrix(a); // writes rate matrix to file rate_matrix.dat
    //gsl_vector * eigenvec = gsl_vector_calloc (evec->size1);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (a->size1);
    // "This function allocates a workspace for computing eigenvalues and eigenvectors 
    //of n-by-n real nonsymmetric matrices. The size of the workspace is O(5n)".  [http://www.gnu.org/software/gsl/manual/html_node/Real-Nonsymmetric-Matrices.html]
    gsl_eigen_nonsymm_params (1, balancing, w->nonsymm_workspace_p);
    // it prevents finding negative concentraion due to having small rate.
    gsl_eigen_nonsymmv (a, eval, evec, w); 
    //  "This function computes eigenvalues and right eigenvectors 
    //of the n-by-n real nonsymmetric matrix A". [http://www.gnu.org/software/gsl/manual/html_node/Real-Nonsymmetric-Matrices.html]

    /*for (int i = 0; i < well_numbers; ++i){
        gsl_complex eigenvalue_i = gsl_vector_complex_get (eval,i);
        if (GSL_REAL(eigenvalue_i) > 0) {
            //cout << "eigen = " << GSL_REAL(eigenvalue_i) << endl;
            GSL_SET_COMPLEX(&eigenvalue_i, 0, 0);
            gsl_vector_complex_set(eval,i,eigenvalue_i);
        }
    }*/

    for (int i = 0; i < well_numbers; i++)
    {
        gsl_complex eigenvalue_i = gsl_vector_complex_get (eval,i);
        if (to_string(abs(GSL_REAL(eigenvalue_i))) == "nan"|| to_string(abs(GSL_IMAG (eigenvalue_i))) == "nan")
            error =-1;
        if (abs(GSL_REAL(eigenvalue_i)) < extreme_small_number) 
        {
            GSL_SET_COMPLEX(&eigenvalue_i, 0, 0);
            gsl_vector_complex_set(eval,i,eigenvalue_i);
        }
    }

    gsl_eigen_nonsymmv_free (w);

    gsl_eigen_nonsymmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_DESC); 
    // "This function simultaneously sorts the eigenvalues stored in the vector eval 
    //and the corresponding real eigenvectors stored in the columns of the matrix evec into ascending or descending order according to the value of the parameter sort_type". Here it sorts in descending order in magnitude. [http://www.gnu.org/software/gsl/manual/html_node/Sorting-Eigenvalues-and-Eigenvectors.html]

    int last_eigenvalue = well_numbers-1;
    gsl_complex eigenvalue =  gsl_vector_complex_get(eval, last_eigenvalue);
    // the first index in array is 0 so the last index is well_numbers-1.
    //cout << "last eigenvalue before setting to zero " << GSL_REAL(eigenvalue) << endl;
    GSL_SET_COMPLEX(&eigenvalue, 0, 0);
    gsl_vector_complex_set(eval,last_eigenvalue,eigenvalue);
    //cout << "last eigenvalue after setting to zero " << GSL_REAL(eigenvalue) << endl;

}
