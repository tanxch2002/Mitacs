

#include "pch.h"
//#include "MitacsGUI.h"
//#include "MitacsGUIDlg.h"
#include <iostream>
#include "generate_energy_landscape.h"
//#include "write_on_files.h"
#include <conio.h>
#include <iomanip>
#include "Global_variables.h"
using namespace std;

//This one calculates coefficients T(E), attempt frequencies,  stationary-state energies, etc
void GenerateEnergyLandscape(
	gsl_vector* coeff_matrixe, gsl_vector* coeff_matrixg, /* coefficients rho */
	gsl_matrix* eigenve, gsl_matrix* eigenvg, /* eigenvalues */
	gsl_vector* TE_matrix_excited, gsl_vector* TE_matrix_ground, /* transmission probabilities */
	gsl_vector* classical_attempt_freq_matrix_e, gsl_vector* classical_attempt_freq_matrix_g, /* attempt frequencies */
	gsl_vector* parabola_ground, gsl_matrix* energy_difference, int ensemble, gsl_vector* st_well,
	int& binning, gsl_vector* v_barrier_height_excited, gsl_vector* v_barrier_height_ground, /* barriers, for determination of barrier distributions, defined with respect to the lowest stationary state */
	gsl_matrix* excited_parabola_for_rate, gsl_matrix* ground_parabola_for_rate, /*low-res energy landscapes, one landscape */
	gsl_rng* r, int& produced_ensemble, gsl_vector* trans, double cutoff_energy_difference)
	//
	// most important variables for monitoring correct data flow: parabola_ground, ground_parabola_for_rate, excited_parabola_for_rate, 
	//ensemble is the number of the current system in the ensemble.
{
	// 
	// parabola_ground is unmodified baseline of the landscape (parabola or flat) coming from the outside
	gsl_vector* parabola_excited = gsl_vector_calloc(point_numbers);
	// these are the non-modified parabolas (including the possibility of flat baseline / parabola coefficient=0) to which one adds tops and bottoms:
	gsl_vector* temp_parabola_ground = gsl_vector_calloc(point_numbers);

	gsl_vector_memcpy(parabola_excited, parabola_ground); //destination, then source
	//temporary backup of baseline / parabola unaffected by random numbers
	gsl_vector_memcpy(temp_parabola_ground, parabola_ground); // backup without barriers

	double translation; // energy difference between baselines of the excited and ground state landscapes
	double energy_test; // the difference between the current transtion frequency and burn frequency
	double para_ground, para_excited, bottom_of_wells_ground, bottom_of_wells_excited, deltaE;

	double current_eigenval_gr, current_eigenval_ex, max_gr, max_ex;
	int max_index_gr, max_index_ex;
		
	double* nn;
	double delta = 0;
	double delta_square = 0;
	//asymmetry; the difference between the two neighboring well's levels. 
	double V_b, V_p, V_r, Energy, Length, h_bar, k_r, k_p, k_bar, k_r2, k_p2, k_b2, transmission_probability_LtoR;
	h_bar = 1.05 * pow(10, -34); // in SI
	
	double eigenfunction_integral, norm_eigenfunction_integral, eigenvector_norm_factor;
	// Now we have to take into account that as we add randomness to the bottoms of the wells, we may have a situation
	// 	  where no well has transition frequency close enough to the burn frequency.
	// 	   In this case we have to redo the random number generation for the well bottoms until we get some well in resonance with the laser
	
	gsl_matrix_complex* eigenmatrix = gsl_matrix_complex_alloc(Hamiltonian_well_numbers, Hamiltonian_well_numbers);
	//The corresponding eigenfunctions are stored in the COLUMNS of the matrix
		//creating the matrix for Miller-colbert algorithm
	gsl_matrix_complex* matrix1 = gsl_matrix_complex_calloc(Hamiltonian_well_numbers, Hamiltonian_well_numbers);
	gsl_vector* eigenv = gsl_vector_calloc(Hamiltonian_well_numbers);
	// Vector which contains eigenvalues / energies, in Joules
	gsl_eigen_hermv_workspace* eigenspace2 = gsl_eigen_hermv_alloc(Hamiltonian_well_numbers);

	// all kinds ofthings related to setting up the x-axis.
	// we need to convert zig-zag landscapes to 512/1024/2048-point landscapes suitable for Miller-Colbert

	gsl_vector* xcoordvector = gsl_vector_calloc(Hamiltonian_well_numbers);
	//Creating a length N vector for assigning x coordinates - long hi-res landscape
	double x_value;
	int number_of_full_flats_before_main = (point_numbers) / 2;
	int borders[2][30]; // Calculating borders of the well, needed to calculate projections rho. First index is 0 or 1 for beginning and end respectively, second index is number of the current well
	
	int pointsperflat = (int)(Hamiltonian_well_numbers / point_numbers); //12 from 12.488 THIS IS WHERE I CALUCLATE THE PONITS OF INTEGRAL
	//Calculates the number of points per flat line (barrier or well bottom); remainder is distributed between the outer top regions
	int leftoverpoints = Hamiltonian_well_numbers - (pointsperflat * point_numbers); //512-492=20 leftover points
	//Calculates number of extra points  that will be distributed to edge regions
	//(If there is an odd number of leftoverpoints, the extra point will go to the end.  
	int extrastartpoints = leftoverpoints / 2;
	int extraendpoints = leftoverpoints - extrastartpoints;
	
	double delta_x = width_well / pointsperflat;     //step of high resolution x-axis
	int center_flat_point;  //center flat point determines the number of the point 
	//within the Middle well/barrier to place the x=0
	int nth_point ;
	double;
	
	if (pointsperflat % 2 == 0)
	{ //divisible by 2
		center_flat_point = (pointsperflat / 2) + 1; //This is the number of the point within the center barrier/well where the x=0 will be centered;
		//If pointsperflat is even, the well will be slightly off center. Say there are 10 points, the 6th point will be chosen. Which will have 5 to its left and 4 points to its right (And therefore a proportional discrepency in length)
	}
	else {
		center_flat_point = (pointsperflat + 1) / 2;
		//If pointsperflat is odd, the well will be centered. Say there are 11 points, this will give the 6th point with even spacing on either side.
	}
	nth_point = extrastartpoints + (pointsperflat * ((point_numbers) / 2)) + center_flat_point;
	double xshift = delta_x * (nth_point - 1.0);
	x_value = (-xshift); //starting point, negative x
	for (int i = 0; i < Hamiltonian_well_numbers; i++) //creating the array of x-coordinates
	{
		gsl_vector_set(xcoordvector, i, x_value);
		x_value = x_value + delta_x;
	}

	for (int i = 0; i < point_numbers; i = i + 2) //had to add = here
	{
		borders[0][i / 2] = extrastartpoints + (i + 1) * pointsperflat; // beginnings
		borders[1][i / 2] = extrastartpoints + (i + 2) * pointsperflat; //ends
		// i=0 borders[0][0]=esp+1xpointsperflat
		//i=0 borders[1][0]esp+2xpoitsperflat
		// i=2 borders [0][1]=esp+3xpoitsperflat... OK
	}
	// supposedly we calculated all borders of the wells (defined as indexes, integers between zero and 512/1024/2048 etc...

	double scaling_constant = pow(h_planck, 2) / (pow(2 * pi * delta_x, 2) * 2 * mass);
	// looks like it is in SI units...
	gsl_complex vec;
	double test;

	gsl_vector* potentialvector = gsl_vector_calloc(Hamiltonian_well_numbers);
	//allocating space for long array containing potential V
	gsl_vector* Vb_vector = gsl_vector_calloc(well_numbers);
	gsl_vector* Vr_vector = gsl_vector_calloc(well_numbers);
	gsl_vector* Vp_vector = gsl_vector_calloc(well_numbers);
	gsl_vector* well_potential_U = gsl_vector_calloc(well_numbers);
	int counter = 0;
	double U, E, velocity;
	double cn_value, help_real; // NB: we did not rename the coefficient, but since the Fall of 2022 the meaning of cn is Garashchuk's projection rho.
	int index, t;
	gsl_complex help_complex;

	int kk = 0; // as of July 9 2021 this is a flag indicating if barriers are really barriers (are not lower than wells)
	double final_random;
	//random number for barrier height based on normal distribution. 
	//
	double barrier_height_ground, barrier_height_excited;
	int jj;

	// TODO: prediction of how much above the bottoms relevant energy levels will be for given md2, bottom distribution, etc.

	do //outside check for resonance, checks if after solving the matrix it is still in resonance
	{
		jj = 0;
		do //looking for ensembles that would have at least one well in resonance with the laser 
		{
			++produced_ensemble;
			// counting all ensembles that were ever tried, can be used to calculate absorption spectrum without any holes,
			//that was made from many systems without caring if they have any resonant wells or not...
			//
			gsl_vector_memcpy(parabola_ground, temp_parabola_ground); //destination then source
			gsl_vector_memcpy(parabola_excited, temp_parabola_ground);
			// in each new try for producing a landscape, baselines / parabolas  should be reset from backup to the "smooth" state without randomness.
			//
			translation = mean_value_of_translation + gsl_ran_gaussian(r, standard_deviation_of_translation);
			// in MHz, randomly generated
			// this one is energy difference between the baselines of excited and ground states if energy landscape baselines are flat
			//or between the bottoms of the parabolas if energy landscape baselines are parabolic
			// /this random number is generated separately using the parameters of inhomogeneous broadening
			//NB: this has to be inside the do-loop, as otherwise some transition energies are too far off 
			//and the program cannot find a match at all / gets into an infinite do loop

			energy_test = translation; // just some very large number to start with 
			//mean of the SDF, approximately

			for (int i = 1; i < point_numbers; i = i + 2) //1,3,7....43 - wells, 0,2,4,6.. 44 barriers
			{
				// going over bottoms, starting from 1, with step 2 //doing at as many times as there are wells
				para_ground = gsl_vector_get(parabola_ground, i);
				// bottom's value before modification with random number, ground state
				para_excited = gsl_vector_get(parabola_excited, i);
				// bottom's value before modification with random number, excited state

				bottom_of_wells_ground = mean + gsl_ran_gaussian(r, deviation);
				//NB: mean is positive, so the average of the bottoms is elevated by some positive number with respect to baseline / parabola
				//
				if (bottom_decoupling == 1) // if excited state is NOT a scaled copy of the ground state
					bottom_of_wells_excited = mean + gsl_ran_gaussian(r, deviation);
				// Nb: this is done with the same parameters as for the ground state. so it will be multipled by the stretch factor (<1) somewhere later / below

				else  bottom_of_wells_excited = bottom_of_wells_ground; // no rescale yet

				//resetting values of bottoms in the vectors that will eventually be used to add barriers:

				gsl_vector_set(parabola_ground, i, para_ground + bottom_of_wells_ground);
				gsl_vector_set(parabola_excited, i, para_excited + bottom_of_wells_excited);
				//parabola_ground and parabola_exited now contain random bottoms

				deltaE = stretch * (para_excited + bottom_of_wells_excited) + translation - (para_ground + bottom_of_wells_ground);
				// This is current transition energy, for this well of this system
				// NB: stretch is used to calculate transition energy, but the excited state landscape itself is not rescaled yet.
				// 
				// TODO: modify this deltaE based on corrections above

				//now we have to compare the current transition frequency with the burn freuency
				//to begin with, energy_test is set to full transition energy, i.e it is too large
				//
				if (abs(deltaE - burn_frequency) < energy_test)
					// the first time this condition is always satisfied.
					//100<15000
				{
					energy_test = abs(deltaE - burn_frequency); // energy_test is vastly reduced after the first step, to ~100

					gsl_vector_set(st_well, ensemble, jj); // remember the current position. If it is good one, we will go out of cycle and this will stay
					//otherwise it will be reset upon trying to remake the system
				} //end if
				++jj; //try next well of the same system

				gsl_matrix_set(ground_parabola_for_rate, ensemble, i, para_ground + bottom_of_wells_ground);
				gsl_matrix_set(excited_parabola_for_rate, ensemble, i, stretch * (para_excited + bottom_of_wells_excited));

			} //end for i
			//at the end of the for loop, we have the number of the well that is the BEST match to burn energy (frequency)
			// NB: potential problem if the first well we tried was best but not close enough
			//
			//Best match may not be good enough, and if it is not, we repeat the procedure

		} while (energy_test > cutoff_energy_difference);
		//
		// This was the part that ensures trying until we generate a system with one resonant-enough well
		//This loop keeps working while the difference between the transition frequency and burn frequency is bigger than cutoff_energy_difference (or bin_step).
		
		gsl_vector_set(trans, ensemble, translation);
		// vector containing difference between excited and ground state BASELINE 
		// so far so good, no strange numbers...
		// 

			//*****************BARRIER HEIGHTS *****************************************
			// we only generate barrier heights after it is ensured that the system actually likely has a resonant well
			// otherwise it is a waste of time
			// 

		for (int i = 0; i < point_numbers; i = i + 2) //0,2,4...40 (for point_numbers=41
		{
			do // the purpose of this extra loop is to remake the barrier if the random number generator
			// produces such a low barrier that it is below the adjacent two bottoms and therefore is not a barrier any more
				//Ground state
			{
				final_random = mu + gsl_ran_gaussian(r, sigma); // one barrier
				// producing a number for barrier top based on Gaussian distribution. 
				//
				if (i == 0) final_random = final_random * 20; // outside-left barrier
				if (i == point_numbers - 1) final_random = final_random * 20; // outside-right barrier
				if (i == 0 && gsl_vector_get(parabola_ground, 1) > final_random) kk = 1; // outer left wall too low
				if (i == 0 && gsl_vector_get(parabola_ground, 1) < final_random) kk = 0; // outer left wall OK
				if (i == point_numbers - 1 && gsl_vector_get(parabola_ground, (point_numbers - 2)) > final_random) kk = 1;
				// outer right wall too low
				if (i == point_numbers - 1 && gsl_vector_get(parabola_ground, (point_numbers - 2)) < final_random) kk = 0;
				// outer right wall OK
				if (i > 0 && i < point_numbers - 2 && (gsl_vector_get(parabola_ground, (i - 1)) > final_random || gsl_vector_get(parabola_ground, (i + 1)) > final_random)) kk = 1;
				// one of the middle barriers too low
				if (i > 0 && i < point_numbers - 2 && (gsl_vector_get(parabola_ground, (i - 1)) < final_random && gsl_vector_get(parabola_ground, (i + 1)) < final_random)) kk = 0;
				// one of the middle barriers OK
			} while (kk == 1);

			gsl_vector_set(parabola_ground, i, gsl_vector_get(parabola_ground, i) + final_random);
			//this vector contains all barrier tops for this particular system/ensemble. bottoms are undefined
			gsl_matrix_set(ground_parabola_for_rate, ensemble, i, gsl_vector_get(parabola_ground, i));

			//
			// EXCITED STATE
			//
			if (barrier_decoupling == 1) // if excited state is generated independently from the ground state... currently always true, need to work on checkboxes
			{
				do // the purpose of this extra loop is to remake the barrier if the random number generator
			// produces such a low barrier that it is below the adjacent well bottom(s) and therefore is not a barrier any more
				{
					final_random = mu + gsl_ran_gaussian(r, sigma);
					// producing a number for barrier top based on Gaussian distribution. 
					if (i == 0) final_random = final_random * 20; // outside-left barrier, deliberately made 20 times higher to ensure that states with energies above inner barriers are still bound states
					if (i == point_numbers - 1) final_random = final_random * 20; // outside-right barrier, deliberately made 20 times higher than it should be
					if (i == 0 && gsl_vector_get(parabola_excited, 1) > final_random) kk = 1; // outer left wall too low
					if (i == 0 && gsl_vector_get(parabola_excited, 1) < final_random) kk = 0; // outer left wall OK
					if (i == point_numbers - 1 && gsl_vector_get(parabola_excited, (point_numbers - 2)) > final_random) kk = 1; // outer right wall too low
					if (i == point_numbers - 1 && gsl_vector_get(parabola_excited, (point_numbers - 2)) < final_random) kk = 0; // outer right wall OK
					if (i > 0 && i < point_numbers - 2 && (gsl_vector_get(parabola_excited, (i - 1)) > final_random || gsl_vector_get(parabola_excited, (i + 1)) > final_random)) kk = 1; // one of the middle barriers too low
					if (i > 0 && i < point_numbers - 2 && (gsl_vector_get(parabola_excited, (i - 1)) < final_random || gsl_vector_get(parabola_excited, (i + 1)) < final_random)) kk = 0; // one of the middle barriers OK
				} while (kk == 1);

				gsl_vector_set(parabola_excited, i, gsl_vector_get(parabola_excited, i) + final_random);
				gsl_matrix_set(excited_parabola_for_rate, ensemble, i, stretch * (gsl_vector_get(parabola_excited, i)));
				// need those two outside to calculate initial populations based on Boltzmann    
				//
			}
			else // otherwise it stays the same as for the ground state 
				 //= excited state is compressed and shifted copy of the ground state
			{
				gsl_vector_set(parabola_excited, i, gsl_vector_get(parabola_ground, i));
				gsl_matrix_set(excited_parabola_for_rate, ensemble, i, stretch * (gsl_vector_get(parabola_ground, i)));
			}
		} // end of i loop


		gsl_vector_memcpy(parabola_ground, temp_parabola_ground);
		// since initial randonmess-free parabola part is calculated once only, outside, it should be reset here, 
		//otherwise for the next system the modified parabola will be used instead of smooth parabola.
		
		
		//anyway, now we have Mehdi-style zig-zag landscapes, for both ground and excited states

		
		// stuff specific to ground state starts here.

		t = 0;  ///big x vector point counter (512/1024/2048)
		// Filling the potential vector:
		for (int j = 0; j < point_numbers; j++)
		{
			if (j == 0)
			{ // beginning edge (barrier)
				while (t < (extrastartpoints + pointsperflat))
				{
					gsl_vector_set(potentialvector, t, gsl_matrix_get(ground_parabola_for_rate, ensemble, j));
					//must be still in MHz
					t = t + 1;
				}
			}

			if (j > 0 && j < (point_numbers - 1)) //segments other than end barriers
			{
				while (t <= ((extrastartpoints - 1) + (j + 1) * (pointsperflat)))
				{
					for (int k = 0; k < pointsperflat; k++)
					{
						gsl_vector_set(potentialvector, t, gsl_matrix_get(ground_parabola_for_rate, ensemble, j));
						t = t + 1;
					}
				}
			}
			if (j == (point_numbers - 1)) //end edge barrier
			{
				while (t < Hamiltonian_well_numbers)
				{
					gsl_vector_set(potentialvector, t, gsl_matrix_get(ground_parabola_for_rate, ensemble, j));
					t = t + 1;
				}
			}
		} // end of for

		nn = (double*)calloc(Hamiltonian_well_numbers, sizeof(double));
		for (int i = 0; i < Hamiltonian_well_numbers; i++)
		{
			test = gsl_vector_get(potentialvector, i) * 1e6 * h_planck / 1.98e-23;
			nn[i] = test;
		}

		//VZ: filling the Miller-colbert matrix for the ground state...
	for (int i = 0; i < Hamiltonian_well_numbers; i++)
		{
			for (int j = 0; j < Hamiltonian_well_numbers; j++)
			{
				if (j != i) //off-diagonal elements
				{
					//gsl_complex vec;
					GSL_SET_COMPLEX(&vec, scaling_constant * pow(-1, j - i) * 2 / pow(j - i, 2), 0); //scaling constant in SI
					//looks like this is in cm-1
					gsl_matrix_complex_set(matrix1, i, j, vec);
				}
				if (j == i) //diagonal
				{
					//gsl_complex vec;
					test = gsl_vector_get(potentialvector, j); // Mhz
					GSL_SET_COMPLEX(&vec, scaling_constant * pow(-1, j - i) * pi * pi / 3 + test * 1e6 * h_planck, 0);
					// in Joules
					gsl_matrix_complex_set(matrix1, i, j, vec);
				}
			}
		}
		gsl_eigen_hermv(matrix1, eigenv, eigenmatrix, eigenspace2); // diagonalizing the matrix
		gsl_eigen_hermv_sort(eigenv, eigenmatrix, GSL_EIGEN_SORT_VAL_ASC);
		gsl_matrix_set_row(eigenvg, ensemble, eigenv);
		// copying eigenvalues for this particular ensemble - 
		// //now eigenvalues are stored in a matrix. This matrix contains all eigenvalues for all molecules.
		// I presume we use it in create_rate_matrix to calculate asymmetries
		//

		//AG: Normalizing The Eigenfunctions
		//eigenfunctions should be normalized to integral of wavefunction squared =1 
		// // 
		// 
		for (int i = 0; i < Hamiltonian_well_numbers; i++)
		{
			test = gsl_vector_get(eigenv, i) / 1.98e-23;
			nn[i] = test;
		}
		for (int i = 0; i < Hamiltonian_well_numbers; i++)
		{
			test = GSL_REAL(gsl_matrix_complex_get(eigenmatrix, i, 0));
			nn[i] = test; 
		}
		for (int i = 0; i < Hamiltonian_well_numbers; i++)
		{
			test = GSL_REAL(gsl_matrix_complex_get(eigenmatrix, i, 1));
			nn[i] = test;
		}
		for (int i = 0; i < Hamiltonian_well_numbers; i++)
		{
			test = GSL_REAL(gsl_matrix_complex_get(eigenmatrix, i, 2));
			nn[i] = test;
		}
		for (int i = 0; i < Hamiltonian_well_numbers; i++)
		{
			test = GSL_REAL(gsl_matrix_complex_get(eigenmatrix, i, 3));
			nn[i] = test;
		}
		 
		// Set a breakpoint here and check if the wavefunction is now properly normalized
		// 
		// 2022: This version is not calculating C_n, it is calculating projections rho, as in Garashchuk's paper.
		// //so it is not calculating Gaussans either
		// 

		for (int i = 0; i < well_numbers; i++)
			//2022: if calculating projections rho, one needs to calculate limits for each well, and then use those limits in the inner cycle. These limits are calculated above.
		{
			for (int j = 0; j < Hamiltonian_well_numbers; j++)
			{ // for every eigenfunction
				cn_value = 0;
				// if calculating projections rho, this cycle involves integration in only one well.
				for (int k = borders[0][i]; k <= borders[1][i]; k++)
				{ //Integration within one well
					help_complex = gsl_matrix_complex_get(eigenmatrix, k, j); 
					// get value of the eigenfunction, that in general can be complex 
					help_real = GSL_REAL(gsl_complex_mul(gsl_complex_conjugate(help_complex), help_complex));
					// we are multiplying complex number with its own complex conjugate, so the result must definitely be real
						// we use GSL_REAL to simply change the type of the variable
					cn_value += help_real * delta_x;
					// 2022: we are re-using variable cn_value and coeff_matrixg and coeff_matrixe below, but its meaning has changed.
						// need to make sure there are no other parts of the program that are using these coefficients.
					//create_rate_matrix checked and changed
				}
				//gsl_matrix_set(coeff_matrixg, i, j, cn_value);
				index = to1D(ensemble, i, j); // molecule, well number, state number
				gsl_vector_set(coeff_matrixg, index, cn_value); // coeff_matrixg is 1D representation of the array that would better be a 3D matrix...
			}
		}
		// This is ground state... 
		// Classical Attempt Frequency, in Hz:
		//
		for (int i = 1; i < point_numbers; i = i + 2)
		{
			test = gsl_matrix_get(ground_parabola_for_rate, ensemble, i);
			gsl_vector_set(well_potential_U, counter, test);
			// potentialvector is in MHz; i is always odd, so i-1 is even, so divisible
			counter = counter + 1;
		}
		// E_kin=E-U (where U is the respective well bottom) and E_kin=mv^2/2.
		// E-U=mv^2/2 then v=sqrt(2(E-U)/m)
		for (int i = 0; i < well_numbers; i++)
		{
			U = gsl_vector_get(well_potential_U, i) * h_planck * 1e6; //in Joules
			for (int j = 0; j < Hamiltonian_well_numbers; j++)
			{
				E = gsl_vector_get(eigenv, j); //eigenvalues are in Joules
				if (E > U)
				{
					velocity = sqrt((2 / mass) * (E - U)); // ~136 m/s, not relativistic
					index = to1D(ensemble, i, j);
					gsl_vector_set(classical_attempt_freq_matrix_g, index, velocity / (2 * width_well));
				}
				else gsl_vector_set(classical_attempt_freq_matrix_g, index, 0);
				// attempt frequency array is 1D representation of the array that would better be 3D...
				// so there are attempt frequencies for every level, well and system in the ensemble
			}
		}
		
		nn = (double*)calloc(Hamiltonian_well_numbers * well_numbers, sizeof(double));
		for (int i = 0; i < Hamiltonian_well_numbers * well_numbers; i++)
		{
			test = gsl_vector_get(classical_attempt_freq_matrix_g, i) * h_planck / 1.98e-23;
			nn[i] = test;
		}
		// NB: here I made some changes on 16/6/2021 to make it easier to track what is initial well and what is final well

		counter = 0;
		// getting Vb and Vr, Vp values .
		for (int i = 2; i < point_numbers - 1; i = i + 2)
		{ // segment number 2 is the first barrier (0 is left wall, 1 is the first well
			gsl_vector_set(Vb_vector, counter, gsl_matrix_get(ground_parabola_for_rate, ensemble, i)); // barrier top
			//bottom of the well to the right of the barrier, product well for left to right tunneling:
			gsl_vector_set(Vp_vector, counter, gsl_matrix_get(ground_parabola_for_rate, ensemble, (i + 1))); //bottom of the well to the right of the barrier, product well for left to right tunneling
			gsl_vector_set(Vr_vector, counter, gsl_matrix_get(ground_parabola_for_rate, ensemble, (i - 1))); //bottom of the well to the left of the barrier, reactant well for left to right tunneling
			//cout << counter << endl;
			counter = counter + 1;
		}

		//calculating Transmission probablity for every eigenvector in each well (using all Vb Vr Vp values).
		//T(E) matrix is (q-1) rows times N columns, each row is one well and contains transmission probablity for each eigenvector energy (eigenvalue).
		
		for (int i = 0; i < well_numbers - 1; i++)
		{ //cycle over wells, right-most well cannot transfer further to the right
			V_b = gsl_vector_get(Vb_vector, i) * 1e6 * h_planck; //barrier height in joule
			V_p = gsl_vector_get(Vp_vector, i) * 1e6 * h_planck; // bottom of the well to the right of the barrier. 
			// for left to right tunneling, this would be "product well" using Garashchuk terminology
			V_r = gsl_vector_get(Vr_vector, i) * 1e6 * h_planck; // bottom of the well to the left of the barrier. 

			for (int j = 0; j < Hamiltonian_well_numbers; j++)
			{ //cycle over energy levels
				Energy = gsl_vector_get(eigenv, j); // current stationary state energy in Joules
				Length = width_well;
				//2021: k_r is reactant and k_p is product
				if (Energy > V_r)	k_r = sqrt(2 * mass * (Energy - V_r)) / h_bar; // reactant well, well to the left of the barrier for left to right transition
				else k_r = 0;
				k_bar = sqrt(2 * mass * abs(V_b - Energy)) / h_bar;
				if (Energy > V_p) k_p = sqrt(2 * mass * (Energy - V_p)) / h_bar; //product well
				else k_p = 0;
				//
				k_r2 = pow(k_r, 2);
				k_b2 = pow(k_bar, 2);
				k_p2 = pow(k_p, 2);

				if (Energy > V_r && Energy > V_p) // if sqrt are not negative and therefore transition meaningful
				{
					if (Energy > V_b)  //eq 13 of Garashchuk
						transmission_probability_LtoR = (4 * k_r * k_b2 * k_p) / (pow((k_r + k_p), 2) * k_b2 + (k_r2 * k_p2 + k_b2 * (k_b2 - k_r2 - k_p2)) * pow(sin(k_bar * Length), 2));
					//energy above the tip of the barrier
					else transmission_probability_LtoR = (4 * k_r * k_b2 * k_p) / (pow((k_r + k_p), 2) * k_b2 + (k_r2 * k_p2 + k_b2 * (k_b2 - k_r2 - k_p2)) * pow(sinh(k_bar * Length), 2));
					//energy below the top of the barrier
				}
				else transmission_probability_LtoR = 0;
				// on the order of 10^-25. effective lambda 28... Ok...
				index = to1D(ensemble, i, j);
				gsl_vector_set(TE_matrix_ground, index, transmission_probability_LtoR);
				// transmission probability array is 1D representation of the array that would better be 3D...

				//these are probabilities of one single act of tunneling and are not yet thermally averaged
				// probabilities are the same for the same barrier for right to left, but since the rest of our logic goes by the well, we have to take it into account
				//to save memory, we can do this simple shift in create_rate_matrix
			} //end of level loop
		//
		} // end of well loop
		//
		gsl_vector_set_zero(potentialvector); //this variable/vector of 512/1024/2048 points will be reused for the excited state
			//
			// End of ground state, beginning of the excited state 88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
			// stuff specific to excited state starts here.88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
			//
			// Mehdi-style excited state landscape is properly scaled by now.

		t = 0;  ///big x vector counter (512/1024/2048)
		for (int j = 0; j < point_numbers; j++)
		{
			if (j == 0)
			{ // beginning edge (barrier)
				while (t <= (extrastartpoints + pointsperflat - 1))
				{
					gsl_vector_set(potentialvector, t, gsl_matrix_get(excited_parabola_for_rate, ensemble, j));
					t = t + 1;
				}
			}

			if (j > 0 && j < (point_numbers - 1))
			{ //segments other than end barriers
				while (t <= ((extrastartpoints - 1) + (j + 1) * (pointsperflat)))
				{
					for (int k = 0; k < pointsperflat; k++)
					{
						gsl_vector_set(potentialvector, t, gsl_matrix_get(excited_parabola_for_rate, ensemble, j));
						t = t + 1;
					}
				}

			}
			if (j == (point_numbers - 1))
			{ //end edge barrier
				while (t <= Hamiltonian_well_numbers - 1)
				{
					gsl_vector_set(potentialvector, t, gsl_matrix_get(excited_parabola_for_rate, ensemble, j));
					t = t + 1;
				}
			}
		}
		//
		//creating the matrix for Miller-colbert algorithm
		gsl_matrix_complex_set_zero(matrix1);

		//VZ: filling the Miller-colbert matrix...
		// looks like it is in SI units...
		for (int i = 0; i < Hamiltonian_well_numbers; i++)
		{
			for (int j = 0; j < Hamiltonian_well_numbers; j++)
			{
				if (j != i)
				{
					GSL_SET_COMPLEX(&vec, scaling_constant * pow(-1, j - i) * 2 / pow(j - i, 2), 0); //scaling constant in SI
					gsl_matrix_complex_set(matrix1, i, j, vec);
				}
				if (j == i)
				{
					test = gsl_vector_get(potentialvector, j);
					GSL_SET_COMPLEX(&vec, scaling_constant * pow(-1, j - i) * pi * pi / 3 + test * 1e6 * h_planck, 0); // in Hz
					//note that we scale the excited state down only at this particular point
					gsl_matrix_complex_set(matrix1, i, j, vec);
				}
			}
		}

		gsl_matrix_complex_set_zero(eigenmatrix);
		gsl_vector_set_zero(eigenv);
		gsl_eigen_hermv(matrix1, eigenv, eigenmatrix, eigenspace2);
		gsl_eigen_hermv_sort(eigenv, eigenmatrix, GSL_EIGEN_SORT_VAL_ASC);
		gsl_matrix_set_row(eigenve, ensemble, eigenv); // copying eigenvalues for this particular ensemble - now eigenvalues are stored in a matrix
		//
		// //
		//AG: Normalizing The Eigenfunctions
		//eigenfunctions should be normalized to integral of wavefunction squared =1 
		// //NB: one needs to check what happens to the normalzation of the delocalized functions...
		//calculate integral
		for (int i = 0; i < Hamiltonian_well_numbers; i++)
		{
			eigenfunction_integral = 0;
			norm_eigenfunction_integral = 0;

			for (int j = 0; j < Hamiltonian_well_numbers; j++)
			{
				eigenfunction_integral += pow(GSL_REAL(gsl_matrix_complex_get(eigenmatrix, j, i)), 2) * delta_x;
			}
			//calculate normalization factor
			eigenvector_norm_factor = sqrt(1 / eigenfunction_integral);

			//normalize eigenmatrix
			for (int k = 0; k < Hamiltonian_well_numbers; k++)
			{
				gsl_matrix_complex_set(eigenmatrix, k, i, gsl_complex_mul_real(gsl_matrix_complex_get(eigenmatrix, k, i), eigenvector_norm_factor));
			}
			//printing normalized eigenfunction integral values.
			for (int u = 0; u < Hamiltonian_well_numbers; u++)
			{
				norm_eigenfunction_integral += pow(GSL_REAL(gsl_matrix_complex_get(eigenmatrix, u, i)), 2) * delta_x;
			}
			// just for verification purposes
		}
		//
		// 
		//VZ: creating matrix containing all coefficients rho (in this version)...
		//
		for (int i = 0; i < well_numbers; i++)
			// if calculating projections rho, one needs to calculate limits for each well, and then use those limits in the inner cycle. These limits are calculated above.
		{
			for (int j = 0; j < Hamiltonian_well_numbers; j++)
			{ // for every eigenfunction
				cn_value = 0;
				// if calculating projections rho, this cycle involves 
				   //integration in only one well.
				for (int k = borders[0][i]; k <= borders[1][i]; k++)
				{ //Integration within one well
					help_complex = gsl_matrix_complex_get(eigenmatrix, k, j); // get value of the eigenfunction, that in general can be complex 
					help_real = GSL_REAL(gsl_complex_mul(gsl_complex_conjugate(help_complex), help_complex)); // we are multiplying complex number with its own complex conjugate, so the result must definitely be real
						// we use GSL_REAL to simply change the type of the variable
					cn_value += help_real * delta_x;
					// we are re-using variable cn_value and coeff_matrixg and coeff_matrixe below, but its meaning has changed.
						// need to make sure there are no other parts of the program that are using these coefficients.
				}
				//gsl_matrix_set(coeff_matrixg, i, j, cn_value);
				index = to1D(ensemble, i, j);
				gsl_vector_set(coeff_matrixe, index, cn_value); // coeff_matrixg is 1D representation of the array that would better be 3D...
			}
		}
		//
		// Classical Attempt Frequency, in Hz:
		//
		gsl_vector_set_zero(well_potential_U);
		counter = 0;
		for (int i = 1; i < point_numbers; i = i + 2)
		{
			gsl_vector_set(well_potential_U, counter, gsl_matrix_get(excited_parabola_for_rate, ensemble, i)); // potentialvector is in MHz
			counter += 1;
		}
		// E_kin=E-U (where U is the respective well bottom) and E_kin=mv^2/2.
		// E-U=mv^2/2 then v=sqrt(2(E-U)/m)
		for (int i = 0; i < well_numbers; i++)
		{
			U = gsl_vector_get(well_potential_U, i) * h_planck * 1e6; //in Joule
			for (int j = 0; j < Hamiltonian_well_numbers; j++)
			{
				E = gsl_vector_get(eigenv, j);
				if (E > U)
				{
					velocity = sqrt((2 / mass) * (E - U));
					test = velocity / (2 * width_well);
					//gsl_matrix_set(classical_attempt_freq_matrix_e, i, j, test);
					index = to1D(ensemble, i, j);
					gsl_vector_set(classical_attempt_freq_matrix_e, index, test); // attempt frequency array is 1D representation of the array that would better be 3D...
				}
				else gsl_vector_set(classical_attempt_freq_matrix_e, index, 0);
			}
		}


		nn = (double*)calloc(Hamiltonian_well_numbers * well_numbers, sizeof(double));
		for (int i = 0; i < Hamiltonian_well_numbers * well_numbers; i++)
		{
			test = gsl_vector_get(classical_attempt_freq_matrix_e, i) * h_planck / 1.98e-23;
			nn[i] = test;
		}
		
		counter = 0;
		// getting Vb and Vr, Vp values .
		for (int i = 2; i < point_numbers - 1; i = i + 2)
		{ // segment number 2 is the first barrier (0 is left wall, 1 is the first well
			gsl_vector_set(Vb_vector, counter, gsl_matrix_get(excited_parabola_for_rate, ensemble, i)); // barrier top
			//bottom of the well to the right of the barrier, product well for left to right tunneling
			gsl_vector_set(Vp_vector, counter, gsl_matrix_get(excited_parabola_for_rate, ensemble, (i + 1))); //bottom of the well to the right of the barrier, product well for left to right tunneling
			gsl_vector_set(Vr_vector, counter, gsl_matrix_get(excited_parabola_for_rate, ensemble, (i - 1))); //bottom of the well to the left of the barrier, reactant well for left to right tunneling
			//cout << counter << endl;
			counter += 1;
		}

		//calculating Transmission probablity for every eigenvector in each well (using all Vb Vr Vp values).

		for (int i = 0; i < well_numbers - 1; i++)
		{ //cycle over wells
			V_b = gsl_vector_get(Vb_vector, i) * 1e6 * h_planck; //barrier height
		   //double V_0 = gsl_vector_get(V0_vector, i); // bottom of the well to the right of the barrier. 
			V_p = gsl_vector_get(Vp_vector, i) * 1e6 * h_planck; // bottom of the well to the right of the barrier. 
			// for left to right tunneling, this would be "product well" using Garashchuk terminology
			V_r = gsl_vector_get(Vr_vector, i) * 1e6 * h_planck; // bottom of the well to the left of the barrier. 

			for (int j = 0; j < Hamiltonian_well_numbers; j++)
			{ //cycle over energy levels
				Energy = gsl_vector_get(eigenv, j); // current stationary state energy in Hz
				Length = width_well;
				//double k_r = sqrt(2 * mass * Energy); // reactant well, well to the left of the barrier for left to right transition
				//2021: k_r is reactant and k_p is product
				// //
				if (Energy > V_r)	k_r = sqrt(2 * mass * (Energy - V_r)) / h_bar; // reactant well, well to the left of the barrier for left to right transition
				else k_r = 0;
				k_bar = sqrt(2 * mass * abs(V_b - Energy)) / h_bar;
				if (Energy > V_p) k_p = sqrt(2 * mass * (Energy - V_p)) / h_bar; //product well
				else k_r = 0;
				//
				k_r2 = pow(k_r, 2);
				k_b2 = pow(k_bar, 2);
				k_p2 = pow(k_p, 2);

				if (Energy > V_r && Energy > V_p) // if sqrt is not negative
				{
					if (Energy > V_b)  //eq 13 of Garashchuk
						transmission_probability_LtoR = (4 * k_r * k_b2 * k_p) / (pow((k_r + k_p), 2) * k_b2 + (k_r2 * k_p2 + k_b2 * (k_b2 - k_r2 - k_p2)) * pow(sin(k_bar * Length), 2));
					//energy above the tip of the barrier
					else transmission_probability_LtoR = (4 * k_r * k_b2 * k_p) / (pow((k_r + k_p), 2) * k_b2 + (k_r2 * k_p2 + k_b2 * (k_b2 - k_r2 - k_p2)) * pow(sinh(k_bar * Length), 2));
					//energy below the top of the barrier
				}
				else transmission_probability_LtoR = 0;
				//gsl_matrix_set(TE_matrix_excited, i, j, transmission_probability_LtoR);
				index = to1D(ensemble, i, j);
				gsl_vector_set(TE_matrix_excited, index, transmission_probability_LtoR);
				// transmission probability array is 1D representation of the array that would better be 3D...
				//these are probabilities of one single act of tunneling and are not yet thermally averaged
				// 1.3 e-12, so effective lambda ~13.7
			} // end of cycle over energies
		} //end of cycle over wells


		//2022: finally, calculating transition energies (defined as differences in stationary states with the largest rho
		// essentially, this is an analog of Franck-Condon principle - the most likely transition is the one with the largest wavefunction overlap...
		// lowest-energy stationary-state wavefunctions are fairly localized... At least for parameters so far used by Jing. Check needed for better md2.
		// largest rhos are 0.99
		
		nn = (double*)calloc(point_numbers, sizeof(double));
		for (int i = 0; i < point_numbers; i++)
		{
			test = gsl_matrix_get(ground_parabola_for_rate, ensemble, i) * 1e6 * h_planck / 1.98e-23;
			nn[i] = test;
		}

		for (int i = 0; i < point_numbers; i++)
		{
			test = gsl_matrix_get(excited_parabola_for_rate, ensemble, i) * 1e6 * h_planck / 1.98e-23;  // landscape shown in cm-1, but saved in MHz
			nn[i] = test;
		}
		//
		energy_test = translation; // resetting for future use in the test 
		jj = 0;
		for (int i = 0; i < well_numbers; i++) // for each well, find index of the state with the largest rho.
		{
			max_gr = 0;
			max_index_gr = 0;
			for (int j = 0; j < well_numbers; j++) // j can run over the number of states, but we are only looking within the lowest-energy group
			{
				index = to1D(ensemble, i, j); // molecule, well number, state number
				cn_value = gsl_vector_get(coeff_matrixg, index);
				if (cn_value > max_gr)
				{
					max_gr = cn_value;
					max_index_gr = j;
				}
			}
			current_eigenval_gr = gsl_matrix_get(eigenvg, ensemble, max_index_gr);
			// eigenvalue corresponding to most likely state

			//excited state
			max_ex = 0;
			max_index_ex = 0;
			for (int j = 0; j < well_numbers; j++) // j can run over the number of states, but we are only looking within the lowest-energy group
			{
				index = to1D(ensemble, i, j); // molecule, well number, state number
				cn_value = gsl_vector_get(coeff_matrixe, index);
				if (cn_value > max_ex)
				{
					max_ex = cn_value;
					max_index_ex = j;
				}
			}
			current_eigenval_ex = gsl_matrix_get(eigenve, ensemble, max_index_ex); // eigenvalue corresponding to most likely state

			deltaE = translation  + current_eigenval_ex/ (1e6 * h_planck) - current_eigenval_gr/ (1e6 * h_planck);
			// This is current transition energy, in MHz, for this well of this system
			
			if (abs(deltaE - burn_frequency) < energy_test)
				// the first time this condition is always satisfied.
				//100<15000
			{
				energy_test = abs(deltaE - burn_frequency); // energy_test is vastly reduced after the first step, to <100 cm-1

				gsl_vector_set(st_well, ensemble, jj); // remember the current position. If it is good one, we will go out of cycle and this will stay
				//otherwise it will be reset upon trying to remake the system
			} //end if
			++jj; //try next well of the same system
			//
			gsl_matrix_set(energy_difference, i, ensemble, (translation * 1e6 * h_planck + current_eigenval_ex - current_eigenval_gr)); // in joules
			// saving energy differences between ground and excited state, 
			
			//eigenvalues are calculated with stretch taken into account, so no further stretch is required
			gsl_matrix_set(ground_parabola_for_rate, ensemble, (2 * i + 1), current_eigenval_gr / (1e6 * h_planck));
			gsl_matrix_set(excited_parabola_for_rate, ensemble, (2 * i + 1), current_eigenval_ex / (1e6 * h_planck));
		} //end of for all wells
		// 
	} while (energy_test > cutoff_energy_difference);

	for (int i = 0; i < well_numbers; i++)
	{
		test = gsl_matrix_get(energy_difference, i, ensemble) / 1.98e-23; // transition energies shown in cm-1, 
		nn[i] = test;
	}
	free(nn); 

	binning = binning + well_numbers - 1; // this is likely the counter of items that have to be binned to create something
	// why -1 ? I guess this is binning specifically to create barrier distributions.
	// well_numbers-1 is the number of barriers between the wells
	gsl_vector_free(temp_parabola_ground);
	gsl_vector_free(parabola_excited);
	gsl_vector_free(Vb_vector);
	gsl_vector_free(Vp_vector);
	gsl_vector_free(Vr_vector);
	gsl_vector_free(eigenv);
	gsl_vector_free(potentialvector);
	gsl_vector_free(xcoordvector); // no exceptions up to here...
	gsl_matrix_complex_free(eigenmatrix);
	gsl_matrix_complex_free(matrix1);
	gsl_vector_free(well_potential_U);
	//gsl_matrix_free(normalized_gauss_matrix);
} // end sub
