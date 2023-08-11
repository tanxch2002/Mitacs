

#include "pch.h"
#include "MitacsGUI.h"
#include "MitacsGUIDlg.h"
#include <iostream>
#include <string>
#include "generate_energy_landscape.h"
#include <gsl/gsl_randist.h>// for producing gaussian distribution in gsl
#include "write_on_files.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <conio.h>
#include <iomanip>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "Global_variables.h"
using namespace std;

//This one calculates coefficients, T-S, attempt frequencies for 
void GenerateEnergyLandscape(
	gsl_matrix* coeff_matrixe, gsl_matrix* coeff_matrixg, /* coefficients c_n */
	gsl_vector* eigenve, gsl_vector* eigenvg, /* eigenvalues */
	gsl_matrix* TE_matrix_excited, gsl_matrix* TE_matrix_ground, /* transmission probabilities */
	gsl_matrix* classical_attempt_freq_matrix_e, gsl_matrix* classical_attempt_freq_matrix_g, /* attempt frequencies */
	gsl_vector* parabola_ground, gsl_matrix* energy_difference, int ensemble, int& starting_well, gsl_matrix* land_scape_file,	
	gsl_vector* random_numbers_binning_ground, gsl_vector* random_numbers_binning_excited, /* barrier heghts for stasitical purposes, all lendscapes in one vector */
	int& binning, gsl_vector* v_barrier_height_excited, gsl_vector* v_barrier_height_ground, /* barriers, for determination of distributions */
	gsl_vector* excited_parabola_for_rate, gsl_vector* ground_parabola_for_rate, /*low-res energy landscapes, one landscape */
	gsl_rng* r,
	int& produced_ensemble, gsl_vector* trans, double cutoff_energy_difference)
	//
	// most important variables for monitoring correct data flow: parabola_ground, ground_parabola_for_rate, excited_parabola_for_rate, 
	//ensemble is the number of the current system in the ensemble.
	//land_scape_file likely contains the same info as ..._parabola_for_rate arrays, so this may be a bit redundant...
{
	//
	//allocating memory...
	//
	gsl_vector* random_ground = gsl_vector_calloc(point_numbers); //barrier tops, 22x2+1 items
	gsl_vector* random_excited = gsl_vector_calloc(point_numbers);

	gsl_vector* bottom_ground = gsl_vector_calloc(point_numbers); //Bottoms of the wells only:
	gsl_vector* bottom_excited = gsl_vector_calloc(point_numbers);


	gsl_vector* some_parabola_ground = gsl_vector_calloc(point_numbers); //  baseline+random
	gsl_vector* some_parabola_excited = gsl_vector_calloc(point_numbers);
	// to make sure that this vector is reset for every ensemble. 
	//checking if input parabola is correct:

	// 
	// parabola_ground is unmodified baseline of the landscape (parabola or flat) coming from the outside
	gsl_vector* parabola_excited = gsl_vector_calloc(point_numbers);

	gsl_vector_memcpy(parabola_excited, parabola_ground); //destination, then source
	//the bottom of parabola in the excited state will be modifed but the places below the barriers will remain untouched.
	// 	   in the 2021 version of the program we should not care, as all rates are calculated from stationary state energies and barrier tops
	//
	// these are the non-modified parabolas (including the possibility of flat baseline / parabola coefficient=0) to which one adds tops and bottoms:
	gsl_vector* temp_parabola_ground = gsl_vector_calloc(point_numbers);
	gsl_vector* temp_parabola_excited = gsl_vector_calloc(point_numbers);
	//temporary backup of baseline / parabola unaffected by random numbers
	//
	gsl_vector_memcpy(temp_parabola_ground, parabola_ground);
	gsl_vector_memcpy(temp_parabola_excited, parabola_ground);
	//
	// to see the distribution of translation after satisfying of burn condition (energy_test > (e.g.)100)
	//
	//
	//
	// allocating memory is done, now opening the files where we remember some results.
	//
	fstream file;
	file.open("translation_distribution.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file translation_distribution.txt";
	}

	if (produced_ensemble == 0) file << "ensemble\ttranslation" << endl; //creating header for the file

	string file_name;
	file_name.append("translation_energy_");
	// the rest of the file name is generated in the main() based on the ground or excited state information.
	file_name.append(to_string(ensemble_number));
	if (curvature_of_parabola == 0)
		file_name.append("Flat");
	else
		file_name.append("Para");
	if (bottom_decoupling == 1)
		file_name.append("_decoupled bottoms");
	else file_name.append("_coupled bottoms");
	if (barrier_decoupling == 1)
		file_name.append("_decoupled barriers");
	else file_name.append("_coupled barriers");
	file_name.append(".dat");

	fstream file1;
	// Mehdi: This is a vector which contains the value for transition energies for each system 
	//(the first element is translation for the first system and so on). 
	//In case that bottom of wells in the excited state are independent of ones in the ground state, 
	//I can not find the translation by looking at the difference between two states. I have to keep this data in a separate file.
	file1.open(file_name, ios::out | ios::app);
	if (!file1) {
		cout << "Error opening the file " << file_name;
	}
	//
	//now we have two files open, not sure what we are writing into them. Looks like writing into files is commented out, so opening and closing files could be commented out too... VZ
	// 	   
	double translation; // energy difference between baselines of the excited and ground state landscapes
	double energy_test; // the difference between the current transtion frequency and burn frequency
	int i, j;
	double para_ground, para_excited, bottom_of_wells_ground, bottom_of_wells_excited, deltaE;


	// Now we have to take into account that as we add randomness to the bottoms of the wells, we may have a situation
	// 	  where no well has transition frequency close enough to the burn frequency.
	// 	   In this case we have to redo the random number generation for the well bottoms until we get some well in resonance with the laser
	// 
	//NB: here we are working with one system that has ground and excited states. 
	//. The loop over the number of systems / whole ensemble (currently 5000) is outside of this sub.


	do
	{
		//looking for ensembles that would have at least one well in resonance with the laser 
		++produced_ensemble; // counting all ensembles that were ever tried, originates from where we calculated absorption spectrum withiut any holes,
		//that was made from many systems without caring if they have any resonant wells or not...
		//
		gsl_vector_memcpy(parabola_ground, temp_parabola_ground); //destination then source
		gsl_vector_memcpy(parabola_excited, temp_parabola_excited);
		// in each new try for producing a landscape, baselines / parabolas  should be reset from backup to the "smooth" state without randomness.
		//

		gsl_vector_set_zero(bottom_ground);// these two vectors should also be reset for each try.
		gsl_vector_set_zero(bottom_excited);

		translation = mean_value_of_translation + gsl_ran_gaussian(r, standard_deviation_of_translation); // in MHz
		// this one is energy difference between the baselines of excited and ground states if energy landscape baselines are flat
		//or between the bottoms of the parabolas if energy landscape baselines are parabolic
		// /this random number is generated separately using the parameters of inhomogeneous broadening
		//NB: this has to be inside the loop, as otherwise some transition energie are too far off and the program cannot find a match at all / infinite do loop

		energy_test = translation; // to find the closest well (in terms of energy) to the burn frequency. 
		//mean of the SDF, approximately

		j = 0;
		for (i = 1; i < point_numbers ; i = i + 2) //1,3,7....43 - wells, 0,2,4,6.. 44 barriers
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
			{
				bottom_of_wells_excited = mean + gsl_ran_gaussian(r, deviation);
				// Nb: this is done with the same parameters as for the ground state. so it will be multipled by the stretch factor (<1) somewhere later / below
			}
			else  bottom_of_wells_excited = bottom_of_wells_ground; 

			gsl_vector_set(bottom_ground, i, bottom_of_wells_ground);
			gsl_vector_set(bottom_excited, i, bottom_of_wells_excited);
			//every second element of these arrays gets defined - randomness only

			// bottoms of the wells, randomness + original smooth baseline (parabola)
			gsl_vector_set(some_parabola_ground, i, para_ground + bottom_of_wells_ground);
			gsl_vector_set(some_parabola_excited, i, para_excited + bottom_of_wells_excited);
			// 

			//resetting values of bottoms in the vectors that will eventually be used to add barriers:
			gsl_vector_set(parabola_ground, i, gsl_vector_get(some_parabola_ground, i));
			gsl_vector_set(parabola_excited, i, gsl_vector_get(some_parabola_excited, i));
			//parabola_ground and parabola_exited now contain random bottoms

			deltaE = stretch * gsl_vector_get(some_parabola_excited, i) + translation - gsl_vector_get(some_parabola_ground, i);
			// This is current transition energy, for this well of this system
			// NB: stretch is used to calculate transition energy, but the excited state landscape itself is not modified.
			//
			gsl_matrix_set(energy_difference, j, ensemble, deltaE);
			//this looks like the matrix with energy differences i.e transition energies for all wells of all ensembles
			//this is later used for creating the pre-burn absorption spectra (binning)
			// NB: perhaps should use difference between lowest eigenstates instead...

			//now we have to compare the current transition frequency with the burn freuency
			//to begin with, energy_test is set to full transition energy, i.e it is too large
			if (abs(deltaE - burn_frequency) < energy_test) // if we are close enough, the first time this condition is always satisfied.
				//100<15000
			{
				energy_test = abs(deltaE - burn_frequency); // energy_test is vastly reduced after the first step, to ~100
				starting_well = j; // remember the current position. If it is good one, we will go out of cycle and this will stay
				//otherwise it will be reset upon trying to remake the system
			} //end if
			++j; //try next well of the same system
		} //end for
		//at the end of the for loop, we have the number of the well that is the best match to transition energy - j
		// NB: potential problem if the first well we tried was best but not close enough
		//
		//Best match may not be good enough, and if it is not, we repeat the procedure

	} while (energy_test > cutoff_energy_difference);
	// This is the part that ensures trying until we generate a system with one resonant-enough well
	//This loop keeps working while the difference between the transition frequency and burn frequency is bigger than cutoff_energy_difference (or bin_step).

	// Here we created just one system that satisfies the resonance condition...
	//the cycle over 5000 systems is outside, and this whole subroutine is called from there.

	file << ensemble << "\t" << translation << endl; // actually writing something to one of those files...
	//file << translation << endl;
	gsl_vector_set(trans, ensemble, translation);
	// vector containing difference between excited and ground state baselines 

	// so far so good, no strange numbers...
	// 
	//copying relevant numbers ot the final product:
	for (i = 1; i < point_numbers; i = i + 2)
	{
	gsl_vector_set(excited_parabola_for_rate, i, (stretch * (gsl_vector_get(some_parabola_excited, i)))); 
	// we do not not shift the excited state to higher energy...
	gsl_vector_set(ground_parabola_for_rate, i, gsl_vector_get(parabola_ground, i));
	}
	// 	   
		//*****************BARRIER HEIGHTS *****************************************
		// we only generate barrier heights after it is ensured that the system actually has a resonant well
		// otherwise it is a waste of time
		// 

    int k = 0; // as of July 9 2021 this is a flag indicating if barriers are really barriers (are not too low)
    gsl_vector_set_zero (random_ground);// to make sure that this vector is reset for every ensemble
    gsl_vector_set_zero (random_excited);// to make sure that this vector is reset for every ensemble
	// Even numbers are barriers, odd numbers are wells
    double final_random;
	//random number for barrier height based on normal distribution. 

    double average_ground, average_excited;
	// these are the variables for the bottom of each barrier. Since the bottoms of the wells are subject to distribution,
	//barrier height is not simply the random number generated below, it is random number minus the average 
	//value of the bottoms of two neighboring wells.
	// 	   This quantity is used for presentation purposes only and is not sent to rate generator
	//
    double barrier_height_ground, barrier_height_excited;
	// so far so good...
	//
    for ( i= 0; i< point_numbers; i=i+2) //0,2,4...40 (for point_numbers=41
	{	
		do // the purpose of this extra loop is to remake the barrier if the random number generator
		// produces such a low barrier that it is below the adjacent two bottoms and therefore is not a barrier any more
			//Ground state
		{
			final_random = mu + gsl_ran_gaussian(r, sigma);
			// producing a number for barrier top based on Gaussian distribution. 
			if (i == 0 && gsl_vector_get(parabola_ground, 1) > final_random) k = 1; // outer left wall too low
			if (i == 0 && gsl_vector_get(parabola_ground, 1) < final_random) k = 0; // outer left wall OK
			if (i == point_numbers - 1 && gsl_vector_get(parabola_ground,(point_numbers - 2)) > final_random) k = 1; 
			// outer right wall too low
			if (i == point_numbers - 1 && gsl_vector_get(parabola_ground, (point_numbers - 2)) < final_random) k = 0; 
			// outer right wall OK
			if (i > 0 && i < point_numbers - 2 && (gsl_vector_get(parabola_ground, (i - 1)) > final_random || gsl_vector_get(parabola_ground, (i + 1)) > final_random)) k = 1; 
			// one of the middle barriers too low
			if (i > 0 && i < point_numbers - 2 && (gsl_vector_get(parabola_ground, (i - 1)) < final_random || gsl_vector_get(parabola_ground, (i + 1)) < final_random)) k = 0; 
			// one of the middle barriers OK
		} while (k == 1);
		
		gsl_vector_set (random_numbers_binning_ground, binning, final_random);
		// creating distribution of ground-state barrier tops (not exactly the same as barrier height)
        gsl_vector_set (random_ground, i, final_random);
		//this vector contains all barrier tops for this particular system/ensemble. bottoms are undefined

		//parabola_ground contains modified bottoms... 
		//average of the bottoms of two wells surrounding the barrier
		if (i == 0) average_ground = gsl_vector_get(parabola_ground, 1); // height of left wall is counted from the bottom of left-most well
		if (i > 0 && i < point_numbers - 2) average_ground = ( gsl_vector_get (parabola_ground, (i-1) ) + gsl_vector_get (parabola_ground, (i+1))) / 2;
		if (i == point_numbers - 1) average_ground = gsl_vector_get(parabola_ground, (point_numbers - 2)); //height of the right wall is counted from the bottom of the right-most well

		barrier_height_ground = gsl_vector_get (random_ground, i) + gsl_vector_get (parabola_ground, i) - average_ground;
		//parabola_ground is not yet modified with barriers, only with bottoms, so even elements are just original baseline
		//
		//
		// EXCITED STATE
		//
		if (barrier_decoupling == 1) // if excited state is generated independently from the ground state
		{
			do // the purpose of this extra loop is to remake the barrier if the random number generator
		// produces such a low barrier that it is below the adjacent well bottom(s) and therefore is not a barrier any more
			{
				final_random = mu + gsl_ran_gaussian(r, sigma);
				// producing a number for barrier top based on Gaussian distribution. 
				if (i == 0 && gsl_vector_get(parabola_ground, 1) > final_random) k = 1; // outer left wall too low
				if (i == 0 && gsl_vector_get(parabola_ground, 1) < final_random) k = 0; // outer left wall OK
				if (i == point_numbers - 1 && gsl_vector_get(parabola_ground, (point_numbers - 2)) > final_random) k = 1; // outer right wall too low
				if (i == point_numbers - 1 && gsl_vector_get(parabola_ground, (point_numbers - 2)) < final_random) k = 0; // outer right wall OK
				if (i > 0 && i < point_numbers - 2 && (gsl_vector_get(parabola_ground, (i - 1)) > final_random || gsl_vector_get(parabola_ground, (i + 1)) > final_random)) k = 1; // one of the middle barriers too low
				if (i > 0 && i < point_numbers - 2 && (gsl_vector_get(parabola_ground, (i - 1)) < final_random || gsl_vector_get(parabola_ground, (i + 1)) < final_random)) k = 0; // one of the middle barriers OK
			} while (k == 1);
		} // otherwise it stays the same as for the ground state = excited state is compressed and shifted copy of the ground state

        gsl_vector_set (random_numbers_binning_excited, binning, stretch * final_random);
		// to see the distribution of potential in the excited state for all of ensemble.
        gsl_vector_set (random_excited, i, final_random);
		//this vector contains all barrier tops for this particular system/ensemble. bottoms are undefined

		if (i == 0) average_excited = gsl_vector_get(parabola_ground, 1); // height of left wall counted from the bottom of left well
		if (i > 0 && i < point_numbers - 2) average_excited = (gsl_vector_get(parabola_ground, (i - 1)) + gsl_vector_get(parabola_ground, (i + 1))) / 2;
		if (i == point_numbers - 1) average_excited = gsl_vector_get(parabola_ground, (point_numbers - 2)); //height of the right wall counted from the bottom of the right well
		//average of the bottoms of adjacent wells
		barrier_height_excited = gsl_vector_get(random_excited, i) + gsl_vector_get(parabola_ground, i) - average_excited;
		
		//By now all barriers are guaranteed to be higher than the bottoms of the corresponding neighboring wells.
		//
        gsl_vector_set (v_barrier_height_ground, i/2, barrier_height_ground); 
		gsl_vector_set (v_barrier_height_excited, i/2, stretch * barrier_height_excited); 
		//these vectors contain barrier heights adjusted to bottoms

		para_ground = gsl_vector_get(parabola_ground, i);
		// top's value before modification with random number, ground state
		para_excited = gsl_vector_get(parabola_excited, i);
		// top's value before modification with random number, excited state

		// tops of the barriers
		gsl_vector_set(some_parabola_ground, i, para_ground + gsl_vector_get(random_ground, i));
		gsl_vector_set(some_parabola_excited, i, para_excited + gsl_vector_get(random_excited, i));
		// randomness + original smooth baseline (parabola)
		
        ++binning;
		gsl_vector_set(excited_parabola_for_rate, i, stretch* (gsl_vector_get(some_parabola_excited, i))); 
		// let's not shift the excited state to higher energy just yet...
		gsl_vector_set (ground_parabola_for_rate, i, gsl_vector_get (some_parabola_ground, i));
        
    } // end of for loop
	
	//
	//SAVING EVERYTHING TO LARGE MATRIX. This matrix is called land_scape_file. But it seems that it is currently NOT saved into any actual file
	//
	j = 0;
    double ground_state, excited_state;
	for (i = 0; i < point_numbers ; i=i+1) 
	{
        gsl_matrix_set (land_scape_file, j, 0, i+1); // well or barrier number. Before I changed it, it probably was position
        gsl_matrix_set (land_scape_file, j, 1, gsl_vector_get (parabola_ground, j)); //Parabola + random bottoms; tops may not be defined
	
		gsl_matrix_set (land_scape_file, j, 2, stretch * gsl_vector_get (some_parabola_excited, j) + translation); //fully defined excited landscape
		// Mehdi: it should be calculated once for each ensemble since I wanted to check that correct points of excited state
		//are chosen for delta calculation.
		//
		//
        gsl_matrix_set (land_scape_file, j, 3 * ensemble+3, gsl_vector_get (random_ground, j)); //barrier tops and bottoms only; ground state, no baseline
		//
        ground_state = gsl_vector_get (parabola_ground, j) + gsl_vector_get (random_ground, j); //barrier tops and bottoms on the baseline
        gsl_matrix_set (land_scape_file, j, 3 * ensemble+4, ground_state);// full true ground state landscape
		//
        excited_state = stretch * (gsl_vector_get (parabola_excited, j)+gsl_vector_get (random_excited, j)) + translation;
		//Excited state = stretch * ground state + translation
        gsl_matrix_set (land_scape_file, j, 3 * ensemble+5, excited_state);
        ++j;
    } //end for

	
    gsl_vector_memcpy(parabola_ground, temp_parabola_ground);
	gsl_vector_memcpy(parabola_excited, temp_parabola_excited);
	// since initial randonmess-free parabola part is calculated once only, outside, it should be reset here, 
	//otherwise for the next system the modified parabola will be used instead of smooth parabola.
    gsl_vector_free (some_parabola_ground);
    gsl_vector_free (some_parabola_excited);
    gsl_vector_free (parabola_excited);
    gsl_vector_free (random_ground);
    gsl_vector_free (random_excited);
    gsl_vector_free (temp_parabola_ground);
    gsl_vector_free (temp_parabola_excited);
    gsl_vector_free (bottom_ground);
    gsl_vector_free (bottom_excited);
	//
    if (ensemble == ensemble_number)
	{
		// Since the data are appending to these files, they should be closed only after the last system.
    file.close();
    file1.close();
    } //end if

	// New 2021: AS HERE ARE TWO POSIBILITIES ABOUT GROUND/EXCITED STATE. WE HAVE TO DO THE CALCULATION TWICE HERE.
	// On a separate note, it looks like we are currently not writing anything to the file1, as all writing to the file1 was commented out
	
	//double w_square = 0;
	double delta = 0;
	double delta_square = 0;
	//asymmetry; the difference between the two neighboring well's minima / bottoms. 

	// 2021: JA Code starts here.
	
	gsl_vector* xcoordvector = gsl_vector_calloc(Hamilton_well_numbers); //Creating a length N vector for assigning x coordinates - long hi-res landscape
	
	//seems OK up to here
	
	int pointsperflat = (int)(Hamilton_well_numbers / point_numbers); //12 from 12.488
	//Calculates the number of points per flat line (barrier or well bottom); remainder is distributed between the outer top regions
	int leftoverpoints = Hamilton_well_numbers - (pointsperflat * point_numbers); //512-492=20 leftover points
	//Calculates number of extra points  that will be distributed to edge regions
	//(If there is an odd number of leftoverpoints, the extra point will go to the end.  
	int extrastartpoints = leftoverpoints / 2;
	int extraendpoints = leftoverpoints - extrastartpoints;
	double delta_x = width_well / pointsperflat;     //step of high resolution x-axis, 1/12 Angstroem
	int center_flat_point;  //center flat point determines the number of the point within the Middle well/barrier to place the x=0
	if (pointsperflat % 2 == 0) { //divisible by 2
		center_flat_point = (pointsperflat / 2) + 1; //This is the number of the point within the center barrier/well where the x=0 will be centered;
		//If pointsperflat is even, the well will be slightly off center. Say there are 10 points, the 6th point will be chosen. Which will have 5 to its left and 4 points to its right (And therefore a proportional discrepency in length)
	}
	else {
		center_flat_point = (pointsperflat + 1) / 2;
		//If pointsperflat is odd, the well will be centered. Say there are 11 points, this will give the 6th point with even spacing on either side.
	}

	double x_value;
	int number_of_full_flats_before_main = (point_numbers) / 2;
	int nth_point = extrastartpoints + (pointsperflat * ((point_numbers) / 2)) + center_flat_point;
	double xshift = delta_x * (nth_point - 1.0);
	x_value = (-xshift); //starting point of the next step, negative x

	for (int i = 0; i <= Hamilton_well_numbers - 1; i++) //creating the array of x-coordinates
	{
		gsl_vector_set(xcoordvector, i, x_value);
		x_value = x_value + delta_x;
	}

	//
	gsl_vector* potentialvector = gsl_vector_calloc(Hamilton_well_numbers); //allocating space for long array containing potential V
	gsl_vector* Vb_vector = gsl_vector_calloc(well_numbers);
	gsl_vector* Vr_vector = gsl_vector_calloc(well_numbers);
	gsl_vector* Vp_vector = gsl_vector_calloc(well_numbers);
	
		// stuff specific to ground state starts here.
	
	int t = 0;  ///big x vector counter (512/1024/2048)
	j = 0;  ///Medhi-point counter
	int b = 0;  ///TEST CODE print counter
	for (j = 0; j < point_numbers; j++) 
	{
		if (j == 0) { // beginning edge (barrier)
			while (t <= (extrastartpoints + pointsperflat - 1)) {
				gsl_vector_set(potentialvector, t, gsl_vector_get(ground_parabola_for_rate, j)); //must be still in MHz
				t = t + 1;
			}
		}

		if (j > 0 && j < (point_numbers - 1)) //segments other than end barriers
		{ 
			while (t <= ((extrastartpoints - 1) + (j + 1) * (pointsperflat))) {
				for (int k = 0; k < pointsperflat; k++) {
					gsl_vector_set(potentialvector, t, gsl_vector_get(ground_parabola_for_rate, j));
					t = t + 1;
				}
			}

		}
		if (j == (point_numbers - 1)) //end edge barrier
		{ 
			while (t <= Hamilton_well_numbers - 1) {
				gsl_vector_set(potentialvector, t, gsl_vector_get(ground_parabola_for_rate, j));
				t = t + 1;
			}

		}
		for (b = 0; b <= (Hamilton_well_numbers - 1); b++) 
		{

			//so far do nothing
		}
	}
			
		//creating the matrix for Miller-colbert algorithm
	gsl_matrix_complex* matrix1 = gsl_matrix_complex_calloc(Hamilton_well_numbers, Hamilton_well_numbers);

	//VZ: filling the Miller-colbert matrix...
	double scaling_constant = pow(h_planck, 2) / (pow(2 * pi * delta_x, 2) * 2 * mass); 
	// looks like it is in SI units...
	gsl_complex vec;
	double test;
	for (int i = 0; i <= Hamilton_well_numbers - 1; i++) {
		for (j = 0; j <= Hamilton_well_numbers - 1; j++) {
			if (j != i) {
				//gsl_complex vec;
				GSL_SET_COMPLEX(&vec, scaling_constant * pow(-1, j - i) * 2 / pow(j - i, 2), 0); //scaling constant in SI
				//looks like this is in cm-1
				gsl_matrix_complex_set(matrix1, i, j, vec);
			}
			if (j == i) {
				//gsl_complex vec;
				test = gsl_vector_get(potentialvector, j);
				GSL_SET_COMPLEX(&vec, scaling_constant * pow(-1, j - i)* pi * pi / 3 + test * 1e6 * h_planck, 0); // in Joules
				gsl_matrix_complex_set(matrix1, i, j, vec);
			}
		}
	}

	gsl_matrix_complex* eigenmatrix = gsl_matrix_complex_alloc(Hamilton_well_numbers, Hamilton_well_numbers);   
	//The corresponding eigenfunctions are stored in the COLUMNS of the matrix
	
	gsl_vector* eigenv = gsl_vector_calloc(Hamilton_well_numbers);          // Vector which contains eigenvalues / energies, in Joules
	gsl_eigen_hermv_workspace* eigenspace2 = gsl_eigen_hermv_alloc(Hamilton_well_numbers);
	gsl_eigen_hermv(matrix1, eigenv, eigenmatrix, eigenspace2);
	gsl_eigen_hermv_sort(eigenv, eigenmatrix, GSL_EIGEN_SORT_VAL_ASC);
	gsl_vector_memcpy(eigenvg, eigenv);
	double eigenfunction_integral, norm_eigenfunction_integral, eigenvector_norm_factor;
	//AG: Normalizing The Eigenfunctions
	//eigenfunctions should be normalized to integral of wavefunction squared =1 
	//calculate integral
	for (int i = 0; i < Hamilton_well_numbers; i++) {
		test = gsl_vector_get(eigenv, i);
		eigenfunction_integral = 0;
		 norm_eigenfunction_integral = 0;

		for (int j = 0; j < Hamilton_well_numbers; j++) {
			eigenfunction_integral += pow(GSL_REAL(gsl_matrix_complex_get(eigenmatrix, j, i)), 2) * delta_x;
		}
		//calculate normalization factor
		eigenvector_norm_factor = sqrt(1 / eigenfunction_integral);
		//cout << "Eigenfunction #" << i << " integral = " << eigenfunction_integral << endl;

		//normalize eigen matrix
		for (int k = 0; k < Hamilton_well_numbers; k++) {
			gsl_matrix_complex_set(eigenmatrix, k, i, gsl_complex_mul_real(gsl_matrix_complex_get(eigenmatrix, k, i), eigenvector_norm_factor));
		}

		//printing normalized eigenfunction integral values.
		for (int u = 0; u < Hamilton_well_numbers; u++) {
			norm_eigenfunction_integral += pow(GSL_REAL(gsl_matrix_complex_get(eigenmatrix, u, i)), 2) * delta_x;
		}
		//cout << "Normalized Eigenfunction #" << i << " integral = " << norm_eigenfunction_integral << endl;
	}
	//
	// 
	//VZ: creating matrix containing all coefficients c_n of decomposition of gaussians into stationary state wavefunctions
	int u;
	double sigma;
	double integral;
	double gauss_integral;
	double gauss_integral_normalized;
	double gauss_center;
	double normalization_factor;
	gsl_matrix* normalized_gauss_matrix = gsl_matrix_calloc(well_numbers, Hamilton_well_numbers);    //AG: defining a matrix containging "q" normalized gauss funtcions with N points each. 
	//VZ: the next variable is the width of the Gaussian wavefunction of the localized state (not eigenstate) expressed via well width - needs to be easily accessible
	int std_deviations_per_well = 4;
	// The 4 value here can be changed. With 4 std. deviations in each well, 95% of the area under the curve is within the well.
	sigma = pow(10, -10) / std_deviations_per_well; //well width 1 angstroem, in reality this value should be passed from other parts of the (large/Mehdi) program
	for (u = 0; u < well_numbers; u++) {   // u loop runs once for each well (q), to create however many gaussians   //SHOULD GO LESS THAN Q
		gauss_center = gsl_vector_get(xcoordvector, extrastartpoints + pointsperflat + (2 * u * pointsperflat) + center_flat_point - 1);
		cout << "Gauss Center = " << gauss_center << endl;
		//
		//VZ the above may need to be slightly modified, as Jing insisted on Gaussian centered on some of the values of x from the coordinate array, 
		// //and it is not really necessary.

		//creating individual Gaussian, each located in one particular well
		gsl_vector* gauss_vector = gsl_vector_calloc(Hamilton_well_numbers);         // generating vector to put in Gaussian values
		gsl_vector* product_vector = gsl_vector_calloc(Hamilton_well_numbers);      // generating vector to multiply gaussian values by stationary state wavefunction values
		integral = 0;
		gauss_integral = 0;
		gauss_integral_normalized = 0;      //AG: defining the integral squared value of the Normalized Gaussian Function.

		for (int i = 0; i < Hamilton_well_numbers; i++) {           // i loop is used for a specific gaussian 
			gsl_vector_set(gauss_vector, i, (1 / (sigma * pow(2 * pi, 1 / 2))) * exp((-pow((gsl_vector_get(xcoordvector, i) - gauss_center) / sigma, 2)) / 2));
			// VZ: proper normalization of Gaussian wavefunctions should be checked BEFORE multiplying Gaussian by stationary state wavefunction; if normalization is incorrect, it should e corrected first.
			// 2021: The multiplying part is commented out, I think that's Okay

			//AG: Taking the integral of the gaussian wavefunction squared
			gauss_integral = gauss_integral + pow(gsl_vector_get(gauss_vector, i), 2) * delta_x;
		}
		normalization_factor = sqrt(1 / gauss_integral);
		
		gsl_vector* gauss_vector_normalized = gsl_vector_calloc(Hamilton_well_numbers);

		for (j = 0; j < Hamilton_well_numbers; j++) {
			gsl_vector_set(gauss_vector_normalized, j, normalization_factor * gsl_vector_get(gauss_vector, j));

			gsl_matrix_set(normalized_gauss_matrix, u, j, gsl_vector_get(gauss_vector_normalized, j));

			//AG: Finding the integral squared value of the Normalized Gaussian Function.
			gauss_integral_normalized = gauss_integral_normalized + pow(gsl_vector_get(gauss_vector_normalized, j), 2) * delta_x;

		}
	}// end of cycle over u

	//gauss_vector_normalized is properly normalized...

	//AG: calculating integral of a product of each Gaussian with every stationary state wavefunction (each row with each column).
	//2021: Here is the multiplying part, AFTER the Normalization.
	//AG: Matrix muliplication qxN * NxN = qxN,
	//Normalized Gauss Matrix * Eigen Matrix * delta_x = Coefficient Matrix
	//
	double cn_value;
	for (i = 0; i < well_numbers; i++) 
	{
		for (j = 0; j < Hamilton_well_numbers; j++) {
			cn_value = 0;
			for ( k = 0; k < Hamilton_well_numbers; k++) {
				cn_value += gsl_matrix_get(normalized_gauss_matrix, i, k) * GSL_REAL(gsl_matrix_complex_get(eigenmatrix, k, j)) * delta_x;
			}
			gsl_matrix_set(coeff_matrixg, i, j, cn_value);
		}
	}
	//
	// Classical Attempt Frequency, in Hz:
	//
	gsl_matrix* classical_attempt_freq_matrix = gsl_matrix_calloc(well_numbers, Hamilton_well_numbers);
	gsl_vector* well_potential_U = gsl_vector_calloc(well_numbers);
	int counter = 0;
	double U, E, velocity; 
	for (int i = 1; i < point_numbers; i += 2)
	{
		test = gsl_vector_get(ground_parabola_for_rate, i);
		gsl_vector_set(well_potential_U, counter, test); // potentialvector is in MHz; i is always odd, so i-1 is even, so divisible
		counter += 1;
	}
	// E_kin=E-U (where U is the respective well bottom) and E_kin=mv^2/2.
	// E-U=mv^2/2 then v=sqrt(2(E-U)/m)
	for (i = 0; i < well_numbers; i++) {
		U = gsl_vector_get(well_potential_U, i) * h_planck*1e6; //in Joule
		for (j = 0; j < Hamilton_well_numbers; j++)
		{
			E = gsl_vector_get(eigenv, j); //eigenvalues are in Joules
			velocity = sqrt((2 / mass) * (E - U)); // 1610 m/s, not relativistic
			gsl_matrix_set(classical_attempt_freq_matrix_g, i, j, velocity / (2 * width_well)); // 8e12
		}
	}

	//AG: calculating tunneling probability using Garashchuk eq.13

	// NB: here I made some changes on 16/6/2021 to make it easier to track what is initial well and what is final well
	
	counter = 0;
	// getting Vb and Vr, Vp values .
	for (int i = 2; i < point_numbers - 1; i += 2) { // segment number 2 is the first barrier (0 is left wall, 1 is the first well
		gsl_vector_set(Vb_vector, counter, gsl_vector_get(potentialvector, i)); // barrier top
		//bottom of the well to the right of the barrier, product well for left to right tunneling
		gsl_vector_set(Vp_vector, counter, gsl_vector_get(potentialvector, (i + 1))); //bottom of the well to the right of the barrier, product well for left to right tunneling
		gsl_vector_set(Vr_vector, counter, gsl_vector_get(potentialvector, (i - 1))); //bottom of the well to the left of the barrier, reactant well for left to right tunneling
		//cout << counter << endl;
		counter += 1;
	}

	//calculating Transmission probablity for every eigenvector in each well (using all Vb Vr Vp values).
	//T(E) matrix is (q-1) rows times N columns, each row is one well and contains transmission probablity for each eigenvector energy (eigenvalue).
	double V_b, V_p, V_r, Energy, Length, h_bar, k_r, k_p, k_b, k_r2, k_p2, k_b2, transmission_probability_LtoR;
	for (int i = 0; i < well_numbers - 1; i++) { //cycle over wells
		V_b = gsl_vector_get(Vb_vector, i)*1e6*h_planck; //barrier height
		//double V_0 = gsl_vector_get(V0_vector, i); // bottom of the well to the right of the barrier. 
		V_p = gsl_vector_get(Vp_vector, i) * 1e6 * h_planck; // bottom of the well to the right of the barrier. 
		// for left to right tunneling, this would be "product well" using Garashchuk terminology
		V_r = gsl_vector_get(Vr_vector, i) * 1e6 * h_planck; // bottom of the well to the left of the barrier. 

		for (int j = 0; j < Hamilton_well_numbers; j++)
		{ //cycle over energy levels
			Energy = gsl_vector_get(eigenv, j); // current stationary state energy in Hz
			Length = width_well;
			//double k_r = sqrt(2 * mass * Energy); // reactant well, well to the left of the barrier for left to right transition
			//2021: k_r is reactant and k_p is product
			// //
			h_bar = 1.05 * pow(10, -34);
			k_r = sqrt(2 * mass * (Energy - V_r)) / h_bar; // reactant well, well to the left of the barrier for left to right transition
			k_b = sqrt(2 * mass * abs(V_b - Energy)) / h_bar;
			//double k_p = sqrt(2 * mass * (Energy - V_0)); //product well
			k_p = sqrt(2 * mass * (Energy - V_p)) / h_bar; //product well
			//New 2021: Some problems here, k_r and k_b are not right.
			//I changed the code from "Joule_E - V_r" to " V_r - Joule_E" and I got positive value. I think it's not right because v_r is the bottom of the well...
			//Maybe there are something wrong with the Joule_E?
			
			k_r2 = pow(k_r, 2);
			k_b2 = pow(k_b, 2);
			k_p2 = pow(k_p, 2);
			transmission_probability_LtoR = 0;

			if (Energy > V_b) { //eq 13 of Garashchuk
				transmission_probability_LtoR = (4 * k_r * k_b2 * k_p) / (pow((k_r + k_p), 2) * k_b2 + (k_r2 * k_p2 + k_b2 * (k_b2 - k_r2 - k_p2)) * pow(sin(k_b * Length), 2));
				//energy above the tip of the barrier
			}
			else {
				transmission_probability_LtoR = (4 * k_r * k_b2 * k_p) / (pow((k_r + k_p), 2) * k_b2 + (k_r2 * k_p2 + k_b2 * (k_b2 - k_r2 - k_p2)) * pow(sinh(k_b * Length), 2));
				//energy below the top of the barrier
			}

			// 2021: Professor ask me to do RtoL but I see no difference of the Garshchuk Function comaring to the LtoR one.
			// 2021: There ARE DIFFERENCES, q-1 means removing one well from calculation, in LtoR, the last well is removed, in RtoL the first well should be removed.
			// 2021: The differences is in classical_attempt_freq_matrix, starts from i+1 to num_wells

			gsl_matrix_set(TE_matrix_ground, i, j, transmission_probability_LtoR);
			//these are probabilities of one single act of tunneling and are not yet thermally averaged
			// 
		}
	//
	}
		gsl_vector_set_zero(potentialvector);
		//
		// //
		// 
	// End of ground state, beginning of the excited state 88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		// stuff specific to excited state starts here.88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		//
		//

		t = 0;  ///big x vector counter (512/1024/2048)
		j = 0;  ///Medhi-point counter
		b = 0;  ///TEST CODE print counter
		for (j = 0; j < point_numbers; j++) {
			if (j == 0) { // beginning edge (barrier)
				while (t <= (extrastartpoints + pointsperflat - 1)) {
					gsl_vector_set(potentialvector, t, gsl_vector_get(excited_parabola_for_rate, j));
					t = t + 1;
				}
			}

			if (j > 0 && j < (point_numbers - 1)) { //segments other than end barriers
				while (t <= ((extrastartpoints - 1) + (j + 1) * (pointsperflat))) {
					for (int k = 0; k < pointsperflat; k++) {
						gsl_vector_set(potentialvector, t, gsl_vector_get(excited_parabola_for_rate, j));
						t = t + 1;
					}
				}

			}
			if (j == (point_numbers - 1)) { //end edge barrier
				while (t <= Hamilton_well_numbers - 1) {
					gsl_vector_set(potentialvector, t, gsl_vector_get(excited_parabola_for_rate, j));
					t = t + 1;

				}

			}
			for (b = 0; b <= (Hamilton_well_numbers - 1); b++) {

				//so far do nothing
			}
		}
		for (int i = 0; i < Hamilton_well_numbers; i++)
			

		//creating the matrix for Miller-colbert algorithm
		gsl_matrix_complex_set_zero(matrix1);

		//VZ: filling the Miller-colbert matrix...
		// looks like it is in SI units...
		for (int i = 0; i <= Hamilton_well_numbers - 1; i++) {
			for (j = 0; j <= Hamilton_well_numbers - 1; j++) {
				if (j != i) {
					//gsl_complex vec;
					GSL_SET_COMPLEX(&vec, scaling_constant * pow(-1, j - i) * 2 / pow(j - i, 2), 0); //scaling constant in SI
					gsl_matrix_complex_set(matrix1, i, j, vec);
				}
				if (j == i) {
					test = gsl_vector_get(potentialvector, j);
					//gsl_complex vec;
					GSL_SET_COMPLEX(&vec, scaling_constant * pow(-1, j - i)*pi * pi / 3 + test * 1e6 * h_planck, 0); // in Hz
					gsl_matrix_complex_set(matrix1, i, j, vec);
				}
			}
		}

		gsl_matrix_complex_set_zero(eigenmatrix); 
		//The corresponding eigenfunctions are stored in the columns of the matrix

		//gsl_vector* eigenv = gsl_vector_calloc(Hamilton_well_numbers);          // Vector which contains eigenvalues / energies, in Joules
		//gsl_eigen_hermv_workspace* eigenspace2 = gsl_eigen_hermv_alloc(Hamilton_well_numbers);
		gsl_eigen_hermv(matrix1, eigenv, eigenmatrix, eigenspace2);
		gsl_eigen_hermv_sort(eigenv, eigenmatrix, GSL_EIGEN_SORT_VAL_ASC);
		gsl_vector_memcpy(eigenve, eigenv);
		//AG: Normalizing The Eigenfunctions
		//eigenfunctions should be normalized to integral of wavefunction squared =1 
		//calculate integral
		for (int i = 0; i < Hamilton_well_numbers; i++) {
			eigenfunction_integral = 0;
			norm_eigenfunction_integral = 0;

			for (int j = 0; j < Hamilton_well_numbers; j++) {
				eigenfunction_integral += pow(GSL_REAL(gsl_matrix_complex_get(eigenmatrix, j, i)), 2) * delta_x;
			}
			//calculate normalization factor
			eigenvector_norm_factor = sqrt(1 / eigenfunction_integral);
			//cout << "Eigenfunction #" << i << " integral = " << eigenfunction_integral << endl;

			//normalize eigen matrix
			for (int k = 0; k < Hamilton_well_numbers; k++) {
				gsl_matrix_complex_set(eigenmatrix, k, i, gsl_complex_mul_real(gsl_matrix_complex_get(eigenmatrix, k, i), eigenvector_norm_factor));
			}

			//printing normalized eigenfunction integral values.
			for (int u = 0; u < Hamilton_well_numbers; u++) {
				norm_eigenfunction_integral += pow(GSL_REAL(gsl_matrix_complex_get(eigenmatrix, u, i)), 2) * delta_x;
			}
			//cout << "Normalized Eigenfunction #" << i << " integral = " << norm_eigenfunction_integral << endl;
		}
		//
		// 
		//VZ: creating matrix containing all coefficients c_n of decomposition of gaussians into stationary state wavefunctions
		
		gsl_matrix_set_zero(normalized_gauss_matrix); 
		//AG: defining a matrix containging "q" normalized gauss funtcions with N points each. 
		//VZ: the next variable is the width of the Gaussian wavefunction of the localized state (not eigenstate) expressed via well width - needs to be easily accessible
		std_deviations_per_well = 4;
		// The 4 value here can be changed. With 4 std. deviations in each well, 95% of the area under the curve is within the well.
		sigma = pow(10, -10) / std_deviations_per_well; //well width 1 angstroem, in reality this value should be passed from other parts of the (large/Mehdi) program
		for (u = 0; u < well_numbers; u++) {   // u loop runs once for each well (q), to create however many gaussians   //SHOULD GO LESS THAN Q
			gauss_center = gsl_vector_get(xcoordvector, extrastartpoints + pointsperflat + (2 * u * pointsperflat) + center_flat_point - 1);
			cout << "Gauss Center = " << gauss_center << endl;
			//
			//VZ the above may need to be slightly modified, as Jing insisted on Gaussian centered on some of the values of x from the coordinate array, 
			// //and it is not really necessary.

			//creating individual Gaussian, each located in one particular well
			gsl_vector* gauss_vector = gsl_vector_calloc(Hamilton_well_numbers);         // generating vector to put in Gaussian values
			gsl_vector* product_vector = gsl_vector_calloc(Hamilton_well_numbers);      // generating vector to multiply gaussian values by stationary state wavefunction values
			integral = 0;
			gauss_integral = 0;
			gauss_integral_normalized = 0;      //AG: defining the integral squared value of the Normalized Gaussian Function.

			for (int i = 0; i < Hamilton_well_numbers; i++) {           // i loop is used for a specific gaussian 
				gsl_vector_set(gauss_vector, i, (1 / (sigma * pow(2 * pi, 1 / 2))) * exp((-pow((gsl_vector_get(xcoordvector, i) - gauss_center) / sigma, 2)) / 2));
				// VZ: proper normalization of Gaussian wavefunctions should be checked BEFORE multiplying Gaussian by stationary state wavefunction; if normalization is incorrect, it should e corrected first.
				// 2021: The multiplying part is commented out, I think that's Okay

				//AG: Taking the integral of the gaussian wavefunction squared
				gauss_integral = gauss_integral + pow(gsl_vector_get(gauss_vector, i), 2) * delta_x;
			}
			normalization_factor = sqrt(1 / gauss_integral);
			
			gsl_vector* gauss_vector_normalized = gsl_vector_calloc(Hamilton_well_numbers);

			for (int j = 0; j < Hamilton_well_numbers; j++) {
				gsl_vector_set(gauss_vector_normalized, j, normalization_factor * gsl_vector_get(gauss_vector, j));

				gsl_matrix_set(normalized_gauss_matrix, u, j, gsl_vector_get(gauss_vector_normalized, j));

				//AG: Finding the integral squared value of the Normalized Gaussian Function.
				gauss_integral_normalized = gauss_integral_normalized + pow(gsl_vector_get(gauss_vector_normalized, j), 2) * delta_x;

			}
		}// end of cycle over u

		//gauss_vector_normalized is properly normalized...

		//AG: calculating integral of a product of each Gaussian with every stationary state wavefunction (each row with each column).
		//2021: Here is the multiplying part, AFTER the Normalization.
		//AG: Matrix muliplication qxN * NxN = qxN,
		//Normalized Gauss Matrix * Eigen Matrix * delta_x = Coefficient Matrix
		//
		//
		for (int i = 0; i < well_numbers; i++)
		{
			for (int j = 0; j < Hamilton_well_numbers; j++) {
				double cn_value = 0;
				for (int k = 0; k < Hamilton_well_numbers; k++) {
					cn_value += gsl_matrix_get(normalized_gauss_matrix, i, k) * GSL_REAL(gsl_matrix_complex_get(eigenmatrix, k, j)) * delta_x;
				}
				gsl_matrix_set(coeff_matrixe, i, j, cn_value); // excited state c_n
			}
		}
		//
		// Classical Attempt Frequency, in Hz:
		//
		//
		counter = 0;
		for (int i = 1; i < point_numbers; i += 2)
		{
			gsl_vector_set(well_potential_U, counter, gsl_vector_get(excited_parabola_for_rate, i)); // potentialvector is in MHz
			counter += 1;
		}
		// E_kin=E-U (where U is the respective well bottom) and E_kin=mv^2/2.
		// E-U=mv^2/2 then v=sqrt(2(E-U)/m)
		for (int i = 0; i < well_numbers; i++) {
			U = gsl_vector_get(well_potential_U, i) * h_planck * 1e6; //in Joule
			for (int j = 0; j < Hamilton_well_numbers; j++)
			{
				 E = gsl_vector_get(eigenv, j);
				velocity = sqrt((2 / mass) * (E - U));
				test = velocity / (2 * width_well);
				gsl_matrix_set(classical_attempt_freq_matrix_e, i, j, test);
			}
		}

		//AG: calculating tunneling probability using Garashchuk eq.13

		// NB: here I made some changes on 16/6/2021 to make it easier to track what is initial well and what is final well

		counter = 0;
		// getting Vb and Vr, Vp values .
		for (int i = 2; i < point_numbers - 1; i += 2) { // segment number 2 is the first barrier (0 is left wall, 1 is the first well
			gsl_vector_set(Vb_vector, counter, gsl_vector_get(potentialvector, i)); // barrier top
			//bottom of the well to the right of the barrier, product well for left to right tunneling
			gsl_vector_set(Vp_vector, counter, gsl_vector_get(potentialvector, (i + 1))); //bottom of the well to the right of the barrier, product well for left to right tunneling
			gsl_vector_set(Vr_vector, counter, gsl_vector_get(potentialvector, (i - 1))); //bottom of the well to the left of the barrier, reactant well for left to right tunneling
			//cout << counter << endl;
			counter += 1;
		}

		//calculating Transmission probablity for every eigenvector in each well (using all Vb Vr Vp values).
		//T(E) matrix is (q-1) rows times N columns, each row is one well and contains transmission probablity for each eigenvector energy (eigenvalue).
		for (int i = 0; i < well_numbers - 1; i++) { //cycle over wells
			 V_b = gsl_vector_get(Vb_vector, i) * 1e6 * h_planck; //barrier height
			//double V_0 = gsl_vector_get(V0_vector, i); // bottom of the well to the right of the barrier. 
			V_p = gsl_vector_get(Vp_vector, i) * 1e6 * h_planck; // bottom of the well to the right of the barrier. 
			// for left to right tunneling, this would be "product well" using Garashchuk terminology
			 V_r = gsl_vector_get(Vr_vector, i) * 1e6 * h_planck; // bottom of the well to the left of the barrier. 

			for (int j = 0; j < Hamilton_well_numbers; j++)
			{ //cycle over energy levels
				Energy = gsl_vector_get(eigenv, j); // current stationary state energy in Hz
				Length = width_well;
				//double k_r = sqrt(2 * mass * Energy); // reactant well, well to the left of the barrier for left to right transition
				//2021: k_r is reactant and k_p is product
				// //
				k_r = sqrt(2 * mass * (Energy - V_r)) / h_bar; // reactant well, well to the left of the barrier for left to right transition
				k_b = sqrt(2 * mass * abs(V_b - Energy)) / h_bar;
				//double k_p = sqrt(2 * mass * (Energy - V_0)); //product well
				k_p = sqrt(2 * mass * (Energy - V_p)) / h_bar; //product well
				//New 2021: Some problems here, k_r and k_b are not right.
				//I changed the code from "Joule_E - V_r" to " V_r - Joule_E" and I got positive value. I think it's not right because v_r is the bottom of the well...
				//Maybe there are something wrong with the Joule_E?

				k_r2 = pow(k_r, 2);
				k_b2 = pow(k_b, 2);
				k_p2 = pow(k_p, 2);
				transmission_probability_LtoR = 0;

				if (Energy > V_b) { //eq 13 of Garashchuk
					transmission_probability_LtoR = (4 * k_r * k_b2 * k_p) / (pow((k_r + k_p), 2) * k_b2 + (k_r2 * k_p2 + k_b2 * (k_b2 - k_r2 - k_p2)) * pow(sin(k_b * Length), 2));
					//energy above the tip of the barrier
				}
				else {
					transmission_probability_LtoR = (4 * k_r * k_b2 * k_p) / (pow((k_r + k_p), 2) * k_b2 + (k_r2 * k_p2 + k_b2 * (k_b2 - k_r2 - k_p2)) * pow(sinh(k_b * Length), 2));
					//energy below the top of the barrier
				}
				gsl_matrix_set(TE_matrix_excited, i, j, transmission_probability_LtoR);
				//these are probabilities of one single act of tunneling and are not yet thermally averaged
				// 
		//
		//
			} // end of cycle over energies
		} //end of cycle over wells
} // end sub
