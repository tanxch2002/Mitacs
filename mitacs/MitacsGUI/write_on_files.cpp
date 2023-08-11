#include "pch.h"

#include "write_on_files.h"
//#include <string> now included in Global_variables.h
#include <iostream>
#include <fstream>
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <direct.h>
#include "Global_variables.h"
using namespace std;
#pragma warning(disable : 4996)


void SetDirectory()
{
    string folder_name;
    //folder_name = to_string(ensemble_number)+" system-cooling time = " + to_string(whole_cooling_time/(1e6*60*60)) + " h- number of burn = "+to_string(number_of_burn)+"- spontaneous recovery time = "+"- spontaneous recovery time = "+to_string(whole_recovery_time/(1e6*60*60))+" h";// this line does the same thing as below lines, however, I believe that doing it in several steps is more understandable.
    folder_name= "\\" + to_string(ensemble_number); // file startes with the number of molecules in the simulations
    if(ensemble_number==1)
        folder_name.append("sy");
    else    folder_name.append("sys");
    if(generate_energy_landscape==1)
        folder_name.append("-NewEl");// new energy landscape (EL) is produced.
    else
       folder_name.append ("-LoadEl");// energy landscape (EL) is loaded from somewhere else.
    if (minimum_frequency != maximum_frequency)
        folder_name.append ("-SMS");
    else
        folder_name.append ("-HB");
    if (burn_in_one_step == 1)
        folder_name.append ("-1Shot");
    if (mean_value_of_translation != burn_frequency)
        folder_name.append("-ShiftedBurning");// Burn wave length is not the same as absorption peak (which is determined by the man value of distribution of energy gap between the ground and excited states).
    if (curvature_of_parabola!=0){
        //folder_name.append("-Parabolic");
        folder_name.append("-curv=");
        folder_name.append(to_string(curvature_of_parabola/energy_conversion_factor));

    } //end if
    folder_name.append("-mu=");
    folder_name.append(to_string(mu/energy_conversion_factor));
    folder_name.append("-sig=");
    folder_name.append(to_string(sigma/energy_conversion_factor));
    folder_name.append("-stretch=");
    folder_name.append(to_string(stretch));
    folder_name.append("-mean=");
    folder_name.append(to_string(mean/energy_conversion_factor));
    folder_name.append("-dev=");
    folder_name.append(to_string(deviation/energy_conversion_factor));
    if (generate_energy_landscape == 0 && scaling_factor != 1)
	{// This is important when an energy landscape is loaded with some scaling factor other than 1.
        folder_name.append("-scaling=");
        folder_name.append(to_string(scaling_factor));
    } //end if

    if (bottom_decoupling == 1)
        folder_name.append("-decoupbott");
    else folder_name.append("-coupbott");
    if (barrier_decoupling == 1)
        folder_name.append("-decoupbarri");
    else folder_name.append("-coupbarri");

    if (minimum_frequency != maximum_frequency)
	{
        folder_name.append("-Scan#");
        folder_name.append(to_string(number_of_scan));
} //end if
    folder_name.append("-Tb=");
    folder_name.append(to_string(burn_temperature));
    if (starting_from_uniform_distribution == 1) folder_name.append("K-StarFromUniDist-");
    else{
        if (starting_from_equilibrium == 0) folder_name.append("K-StarFromNonEqu-");
        else folder_name.append("K-StarFromEqu-");
    } //end else
    folder_name.append("coolingt=" );
    folder_name.append(to_string(whole_cooling_time/(1e6*60*60)));
    folder_name.append("h-#ofburn=");
    folder_name.append(to_string(number_of_burn));
    folder_name.append("-SponRecovt=");
    folder_name.append(to_string(whole_recovery_time/(1e6*60*60)));
    folder_name.append("h-");
    folder_name.append(to_string(number_of_recovery_measurement));
    folder_name.append("steps");
    if(analytical_burn_calculation == 1 && recovery_between_two_acts_of_burn == 0)
        folder_name.append("-Analytical");
    else folder_name.append("-Multisteps");
    if(recovery_correction_during_frequency_scan == 1)
        folder_name.append("-ScanRecCorreON");
    else folder_name.append("-ScanRecCorreOFF");

    if(recovery_between_two_acts_of_burn == 1)
	{
        folder_name.append("-BurningRecCorrecON=");
        folder_name.append(to_string(time_interval_between_two_consecutive_burns));
        folder_name.append("s");

    } //end if
    else folder_name.append("-BurningRecCorrecOFF");
	//
	int n = folder_name.length();
	char char_array[500];
	strcpy(char_array, folder_name.c_str());
	if (_mkdir(char_array) == 0) cout << "New directory created";
		_chdir(char_array);
    //Mehdi deleted all old files then will write new files in the same folder.
	//not sure this is necessary...
    //QDirIterator files(QDir::current(), QDirIterator::Subdirectories);
    //while ( files.hasNext() )
        //fstream::remove( files.next() );

} //end of sub

void WriteParabola (gsl_vector * parabola, string filename)
{
	fstream new_file;
    new_file.open(filename, ios::out | ios::app);
	// the file name is chosen in the main() based on the ground or excited state information.it could be parabola or disordered straight line.
	if (!new_file) {
		cout << "Error opening the file" << filename;
	} //end if
    new_file << "\tX" << "\tBaseline"<< endl;
    int j = 0;
    for (double i = 0; i < 2 * interval; i += d_distance_of_two_wells)
	{
		new_file << i << "\t" << gsl_vector_get (parabola, j) << endl;
        ++j;
    } //end for
    new_file.close();
} //end sub

void WriteEnergyLandscape (gsl_matrix * land_scape_file) 
{
    string file_name;
    file_name.append("energy_landscape_");
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

    fstream file;
		    if (generate_energy_landscape == 1){
			file.open(file_name.append(".dat"), ios::out | ios::app);
    		if (!file) {
			cout << "Error opening the file" << file_name;
		} //end if
        //file.setRealNumberPrecision(45);//This determines the precision of data to be written in the file. 15 figures guarantees that there would not be a round-off error when this landscape will be used in another program.
    }//end if

			fstream file2;
			file.open(file_name.append(".txt"), ios::out | ios::app);// to view the landscape
	if (!file2) {
		cout << "Error opening the file" << file_name;
	} //end if
    
    file2 << "X" << "\tParabola" << "\tMiddleNumber";
    for ( int ensemble = 0; ensemble < ensemble_number; ++ensemble)
        file2 << "\tRandomNumbers(" << ensemble << ")" << "\tGroundState(" << ensemble << ")" << "\tExcitedState(" << ensemble << ")" ;
    file2 << endl;

    for ( int i = 0; i < point_numbers; i++)
    {
        for (int j = 0; j < 3 * ensemble_number + 3; j++)
		{
            if (generate_energy_landscape == 1)
                file << gsl_matrix_get (land_scape_file, i, j) << "\t";// the tab command in the last run of loop produces an extra blank column. To avoid of this column, one can use a condition to add tab to all column except for the last one, like: if (j != 3 * ensemble_number + 3) streamToWrite2 <<"\t"

            file2 << gsl_matrix_get (land_scape_file, i, j) << "\t";// the tab command in the last run of loop produces an extra blank column. To avoid of this column, one can use a condition to add tab to all column except for the last one, like: if (j != 3 * ensemble_number + 3) streamToWrite2 <<"\t"
        } //end for
        if (generate_energy_landscape == 1)
            file << "\n";
        file2 << "\n";
      } //end for
    file.close();
    file2.close();
} //end sub

/*writing the rate matrix in a file.*/
void WriteRateMatrix (gsl_matrix * rate_matrix)
{
	fstream file;
	file.open("rate_matrix.dat", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file rate_matrix.dat";
	}
	    for (int i = 0; i < well_numbers; i++)
        {
        for (int j = 0; j < well_numbers; j++)
            file << gsl_matrix_get (rate_matrix, i, j) << "\t";
        file << endl;}
    file.close();
}


/*Writing A0*/

void WriteA0 (gsl_vector * A0)
{
	fstream file;
	file.open("A0.dat", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file A0.dat";
	}
    for(int i = 0; i < well_numbers; i++) file << gsl_vector_get(A0, i) << endl;
    file.close();
}

/*Writing B, it is unnecessary and I could use WriteA0 function instead*/

void WriteB (gsl_vector * B)
{

	fstream file;
	file.open("B.dat", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file B.dat";
	}
    for(int i = 0; i < well_numbers; i++) file << gsl_vector_get(B, i) << endl;
    file.close();
}
/*writing U in a file, it is unnecessary and I could use WriteRateMatrix function instead.*/
void WriteU (gsl_matrix * U)
{
	fstream file;
	file.open("U.dat", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file U.dat";
	}
    for (int i = 0; i < well_numbers; i++)
	{
        for (int j = 0; j < well_numbers; j++)
            file<< gsl_matrix_get (U, i, j) << "\t";
        file << endl;
    }
    file.close();
}



/*Writing the rates of ensemble which is a huge matrix for large ensemble (18 Gb for ensemble = 24000) and it should be omitted for real program */

void WriteConcentrationAfterSimulationTime (gsl_matrix * rate_file, double number_of_row, double time_step, string filename)
{
    double time = 0;
    for (int i = 0; i < number_of_row + round_off; i++)
    {
        gsl_matrix_set (rate_file,i, 0, time/(1e9));// first column of rate_file.txt is time or X axis
        time += time_step;
	} //end for
    fstream file;
    file.open(filename, ios::out | ios::app);// the file name is chosen in the main() based on the ground or excited state information.
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    // these four lines produce the first row of rates_of _ensemble.txt which includes the names of each columns.
    file << "Time-S";
    for (int ensemble = 0 ; ensemble < ensemble_number; ensemble++)
        for (int i = 0; i < well_numbers; ++i)
            file << "\twell(" << i << "_" << ensemble <<")";
    file << "\n";
    for (int i = 0; i < number_of_row + round_off; ++i)
	{
        for (int j = 0; j < well_numbers * ensemble_number + 1; j++)
		{
            file << gsl_matrix_get (rate_file, i, j) << "\t";
        } //end for
        file << endl;
    } //end for
    file.close();
}

//*********** writting last row of rates for each well and corresponding ensemble or hole vesrus well_number******************//
  void WriteLastRowOfRateEquation (gsl_matrix * last_rates, string filename)
  {
	fstream file;
    file.open(filename, ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "Well number";
    for (int column_counter = 0; column_counter < ensemble_number; column_counter++)
        file << "\tensemble " << column_counter;
    file <<"\n";

    for (int row_counter = 0; row_counter < well_numbers; row_counter++)
    {
        file << row_counter;
        for (int column_counter = 0; column_counter < ensemble_number; column_counter++)
        {
            if (to_string(gsl_matrix_get(last_rates, row_counter, column_counter)) == "nan") {
                std::cout << "There is a problem in system " << column_counter << " in stage of " << filename << std::endl;
                //throw;
            } //end if
            file << "\t" << gsl_matrix_get (last_rates, row_counter, column_counter) ;
        } //end for
        file << "\n";

    } //end for
    file.close();

} //end sub



// writting total frequency in a file. The first well_numbers rows (e.g. random numbers= 111, 0-110) are frequency-probability for the first ensemble, the second well_numbers rows (111-221) are frequency-probability for second ensemble,...
void WriteWholeFrequencyOfEnsemble(gsl_matrix * total_frequency, string filename){

    fstream file;
    file.open(filename, ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file" << filename;
	}

    //file << "Frequency" << "\tProbability" << endl;//WARNING: I have to comment out this line later. July 29, 2014
    for (int i = 0; i < all_frequency; i++)
        file << gsl_matrix_get (total_frequency, i, 0) << "\t" << gsl_matrix_get (total_frequency, i, 1) << endl;


    file.close();

}

//as this function is not necessary to run each time I put some repetative codes (by comparison with Binning.cpp) inside of it.

void WriteDetailOfBinning(gsl_matrix * total_frequency, double bin_step, double smallest_frequency, string filename)
{
    double summation_bin;
    fstream file;// this file helps to follow the summation of different frequencies. The summation over the whole probability of first bin is the first value in the probability column in frequency_binned.txt and so on.
    file.open(filename, ios::out | ios::app);//"frequency_binned_details_probability.txt" or "frequency_binned_details_hole.txt"
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "bin_number" << "\tFrequency" << "\tProbability" << endl;
    for (int i = 0; i <= total_bins; i++)
    {// the less than or equal to signs are applied here to take into account the maximum frequency.
        summation_bin = 0;
        for (int j = 0; j < all_frequency; j++)
            if (gsl_matrix_get (total_frequency, j, 0) >= smallest_frequency + i * bin_step &&
                gsl_matrix_get (total_frequency, j, 0) < smallest_frequency + (i+1) * bin_step){
            summation_bin += gsl_matrix_get (total_frequency, j, 1);
            file << i << "\t" << gsl_matrix_get (total_frequency, j, 0) << "\t" << gsl_matrix_get (total_frequency, j, 1) << endl;// by studying of this file, one can sort the data of total_frequecny in qtiplot then bin it by him/herself and compare the result with this file.
        } //end if
    } //end for
    file.close();
} //end sub




void WriteBinnedFrequency (gsl_vector * bin, double bin_step, double smallest_frequency, string filename)
{
    fstream file;// the final binned result. This file is the final result of this program so far.
    file.open(filename, ios::out | ios::app);// "frequency_binned.txt" or "hole spectrum.txt"
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "Frequency_GHz" << "\tProbability_b" << endl;
    for (int i = 0; i <= total_bins; i++)
        file << (smallest_frequency + i * bin_step)/energy_conversion_factor << "\t" << gsl_vector_get (bin, i) << endl;// the binned value is devoted for the begining frequency of each bin. It could be done by considering binned value for the middle of each bin.//if energy is in MHz then by dividing it by energy_conversion_factor, it will be converted to GHz.

    file.close();
} //end sub

void WriteDistribution (gsl_vector * random_numbers_binning, gsl_vector * potential_binning, gsl_vector * lambda_binning, string filename)
{
    fstream file;
    file.open(filename, ios::out | ios::app);// the file name is chosen in the main() based on the ground or excited state information.
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "X" << "\tRandom_Numbers" << "\tPotential" << "\tLambda\n";
    for (int binning = 0; binning < (well_numbers-1) * ensemble_number; binning++)
		// as the random numbers are more than lambda and barriers (every ensemble has two random numbers) the remainign data will add to the file in the next loop.
        file << binning <<"\t" << gsl_vector_get (random_numbers_binning, binning) << "\t" << gsl_vector_get (potential_binning, binning) << "\t" << gsl_vector_get (lambda_binning, binning) << "\n";
    for (int binning = (well_numbers-1) * ensemble_number; binning < well_numbers * ensemble_number; binning++)
	{
        file << binning <<"\t" << gsl_vector_get (random_numbers_binning, binning) << "\n";
    }
    file.close();
}

//************************************writting resonant well number in SMS****************************//

void WriteResonantWell (gsl_matrix * time_delay_in_frequency_scan)
{
			fstream file;
			file.open("resonant_well_number_sms.txt", ios::out | ios::app);// it keeps resonant well number in each ensemble and in each scan step of frequency.
			if (!file) {
				cout << "Error opening the file resonant_well_number_sms.txt";
			}
    file << "BurnFrequency-GHz";
    for (int i = 0; i < ensemble_number; i++)
        file  << "\tResonant_well-" << i;
    file << endl;
        for (int i = 0; i < well_numbers; i++)
            file << gsl_matrix_get (time_delay_in_frequency_scan,i, 2)/energy_conversion_factor << "\t"
                    << gsl_matrix_get (time_delay_in_frequency_scan,i, 0) << endl;
    file.close();
} //end sub

//************************************writting the probabilities of each well after each scan in SMS****************************//
    void WriteProbabilitiesVsWellAfterEachScan (gsl_matrix * probabilities_in_each_scan, string filename)
		{
    fstream file;// This is probability vs well number (probabilities-vs-well_after_each_scan_sms)
    file.open(filename, ios::out | ios::app);// "frequency_binned.txt" or "hole spectrum.txt"
	if (!file) {
		cout << "Error opening the file"<< filename;
	}
    file << "Well number";
	for (int i = 0; i < ensemble_number; i++)
	{
		for (int j = 0; j < number_of_scan; j++)
			file << "\tScan_" << j << "_Syst_" << i;
		file << endl;
	} //end for
    for (int i = 0; i < well_numbers; i++)
	{
        file << i;
        for (int j = 0; j < number_of_scan * ensemble_number; j++)
            file << "\t" << gsl_matrix_get (probabilities_in_each_scan, i, j);
        file << endl;
    } //end for
    file.close();
} //end sub

//************************************writting the probabilities vs frequency of each well after each scan in SMS****************************//
void WriteProbabilitiesVsFrequencyAfterEachScan (gsl_vector * sort, gsl_matrix * probabilities_in_each_scan, string filename)
{
    fstream file;// This is probability vs ordered frequency (probabilities-vs-frequency_after_each_scan_sms)
    file.open(filename, ios::out | ios::app);// "frequency_binned.txt" or "hole spectrum.txt"
	if (!file) {
		cout << "Error opening the file"<< filename;
	}

    file << "Frequency";
	for (int i = 0; i < ensemble_number; i++)
	{
		for (int j = 0; j < number_of_scan; j++)
			file << "\tScan_" << j << "_Syst_" << i;
		file << endl;
	} //end for
    for (int i = 0; i < well_numbers; i++)
	{
        file << gsl_vector_get(sort,i) / energy_conversion_factor;// NB: energy_conversion_factor: Frequency should be in GHz. the first column (number 0) is ordered frequency
        for (int j = 0; j < number_of_scan * ensemble_number; j++)
            file << "\t" << gsl_matrix_get (probabilities_in_each_scan, i, j);
        file << endl;
    } //end for
    file.close();
} //end sub

//************************************writting the probabilities in resonant well after each act of burn for each system during each scan in SMS****************************//

void WriteProbabilitiesInResonantWell (gsl_matrix * prob_in_reson_well_after_each_act_burn)
{
			fstream file;
			file.open("probability_in_resonant_well_after_each_act_of_burn_for_ensmble.txt", ios::out | ios::app);// it keeps probability in resonant well after each act of burn in each system and in each scan.
			if (!file) {
				cout << "Error opening the file probability_in_resonant_well_after_each_act_of_burn_for_ensmble.txt";
			}
    file << "system";
    for (int system_counter = 0; system_counter < ensemble_number; system_counter++)
	{
        for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++)
		{
            for(int resonant_freq_counter = 0; resonant_freq_counter < well_numbers; resonant_freq_counter++)
			{
                for (int burn_counter = 0; burn_counter < number_of_burn; burn_counter++)
				{
                    file << "\tScan_" << scan_counter << "_ResonFreqNumber_" << resonant_freq_counter << "_BurnNumb_" << burn_counter;
                } //end for
            } //ened for
        } //end for
    } //end for

    file << endl;
    int prob_in_reson_well_after_each_burn_counter = 0;
    for (int system_counter = 0; system_counter < ensemble_number; system_counter++)
	{
        file << system_counter;
        prob_in_reson_well_after_each_burn_counter = 0;
        for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++)
		{
            for(int resonant_freq_counter = 0; resonant_freq_counter < well_numbers; resonant_freq_counter++)
			{
                for (int burn_counter = 0; burn_counter < number_of_burn; burn_counter++)
				{
                    file << "\t" << gsl_matrix_get(prob_in_reson_well_after_each_act_burn, system_counter, prob_in_reson_well_after_each_burn_counter);
                    ++prob_in_reson_well_after_each_burn_counter;
                }
            }
        }
    }
    file.close();

} //end sub


//************************************writting the averaged probabilities in resonant well after each scan for each system in SMS****************************//

void WriteAvgProbabilitiesInResonantWell (gsl_matrix * averaged_probability_in_resonant_well)
		{

			fstream file;
			file.open("averaged_probability_in_resonant_well_after_each_scan.txt", ios::out | ios::app);// it keeps the averaged probability in resonant wells after each scan in each system.
			if (!file) {
				cout << "Error opening the file averaged_probability_in_resonant_well_after_each_scan.txt";
			}
    file << "System";
    for (int system_counter = 0; system_counter < ensemble_number; system_counter++
        ){
        for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++)
        {// this loop goes over the whole band.
            //file << "\tScan_" << scan_counter;
            for (int resonant_freq_counter = 0; resonant_freq_counter < well_numbers; resonant_freq_counter++) //By considering the fact that by choosing correct value for frequency step all wells are in resonance with the laser once only, the number of well is the same as number of frequency.
                file << "\tScan_" << scan_counter << "_ResonFreqNumber_" << resonant_freq_counter;
        }
    }

    file << endl;
    int avg_prob_reson_well_counter = 0;
    for (int system_counter = 0; system_counter < ensemble_number; system_counter++)
    {
        file << system_counter;
        avg_prob_reson_well_counter = 0;
        for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++)
        {
            for (int resonant_freq_counter = 0; resonant_freq_counter < well_numbers; resonant_freq_counter++)
            {
                file <<"\t" << gsl_matrix_get(averaged_probability_in_resonant_well, system_counter, avg_prob_reson_well_counter);
                ++avg_prob_reson_well_counter;
            }

        }
        file << endl;
    }
    file.close();

}


//****It works only for one system**********writting the averaged probabilities in resonant well after each scan versus frequency in SMS****************************//
        void WriteAvgProbabilitiesInResonantWellVsFrequency (gsl_matrix * averaged_probability_in_resonant_well, gsl_vector * sorted_frequency){// in case that this program is used for burning purpose, the minimum_frequency is equal to maximum_frequency and then the averaged for starting_well is set only. NB. I could abonden this function for this situation since it is not useful.

			fstream file;
			file.open("averaged_probability_in_resonant_well_after_each_scan_vs_Frequency.txt", ios::out | ios::app);// it keeps the averaged probability in resonant wells after each scan versus frequency. Meaning that frequency on x-axis is sorted in ascending way and the change of averaged probability in resonant wells after each scan is written in the y-axis. If there is more than one system, then it needs some modifications.
			if (!file) {
				cout << "Error opening the file averaged_probability_in_resonant_well_after_each_scan_vs_Frequency.txt";
			}
    file << "Frequency";

    for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++)// this loop goes over the whole band.
        file << "\tScan_" << scan_counter;


    file << endl;

   
        for (int system_counter = 0; system_counter < ensemble_number; system_counter++)
        {
            for (int resonant_freq_counter = 0; resonant_freq_counter < well_numbers; resonant_freq_counter++)
            { //By considering the fact that by choosing correct value for frequency step all wells are in resonance with the laser once only, the number of well is the same as number of frequency.
                file << gsl_vector_get (sorted_frequency, resonant_freq_counter) / energy_conversion_factor;
                for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++)
                {// this loop goes over the whole band.
                    file  << "\t" << gsl_matrix_get(averaged_probability_in_resonant_well,system_counter,well_numbers*scan_counter+resonant_freq_counter);
                }

                file << endl;
            }

        }

    file.close();

}


void WriteBinnedAvgProProbabilitiesInResonantWellVsFrequency (gsl_matrix * averaged_probability_in_resonant_well, gsl_vector * sorted_frequency){


	fstream file;
	file.open("binned_averaged_probability_in_resonant_well_after_each_scan_vs_Frequency.txt", ios::out | ios::app);// it keeps the averaged probability in resonant wells after each scan versus frequency. Meaning that frequency on x-axis is sorted in ascending way and the change of averaged probability in resonant wells after each scan is written in the y-axis. If there is more than one system, then it needs some modifications.
	if (!file) {
		cout << "Error opening the file binned_averaged_probability_in_resonant_well_after_each_scan_vs_Frequency.txt";
	}
    file << "Binned_Frequency";

    for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++)// this loop goes over the whole band.
        file << "\tScan_" << scan_counter;


    file << endl;
    gsl_vector * binned_frequency = gsl_vector_calloc (averaged_probability_SMS_bin_number);// it keeps binned frequency.
    gsl_matrix * binned_probability = gsl_matrix_calloc (averaged_probability_SMS_bin_number, number_of_scan);// it keeps the probability of each binned frequency after each scan.
    double bin_step = (maximum_frequency - minimum_frequency) / averaged_probability_SMS_bin_number;
    double frequency, sum_frequency, avg_frequency, sum_probability, avg_probability;

    sum_probability = 0;
    avg_probability = 0;
    int number_of_well_in_bin = 0;// the number of frequencies (well) which their energy fits in a bin.
    gsl_matrix * bin_number = gsl_matrix_calloc (well_numbers,2);// the first column is bin number and the second column keeps the well number which happens to be in that bin.
    int bin_counter, previous_bin;
    for (int system_counter = 0; system_counter < ensemble_number; system_counter++)
    {
        for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++)
        {// this loop goes over the whole band.
            previous_bin = 0;
            for (int resonant_freq_counter = 0; resonant_freq_counter < well_numbers; resonant_freq_counter++)
            {
                frequency = gsl_vector_get (sorted_frequency, resonant_freq_counter);

                for (bin_counter = previous_bin; bin_counter < averaged_probability_SMS_bin_number; bin_counter++)
                {
                    if (frequency >= minimum_frequency + bin_counter*bin_step && frequency < minimum_frequency + (bin_counter+1)*bin_step){
                        gsl_matrix_set (bin_number, resonant_freq_counter, 0, bin_counter);
                        gsl_matrix_set (bin_number, resonant_freq_counter, 1, resonant_freq_counter);
                        sum_frequency += frequency;
                        sum_probability += gsl_matrix_get(averaged_probability_in_resonant_well,system_counter,well_numbers*scan_counter+resonant_freq_counter);
                        ++number_of_well_in_bin;
                        previous_bin = bin_counter;// next search for bin starts from this bin number.
                        break;// when the proper bin is found, there is no need to go over all bins.
                    }
                    number_of_well_in_bin = 0;
                    sum_frequency = 0;
                    avg_frequency = 0;
                    sum_probability = 0;
                    avg_probability = 0;
                }

                if(number_of_well_in_bin != 0){
                    avg_frequency = (sum_frequency/energy_conversion_factor) / number_of_well_in_bin;
                    gsl_vector_set(binned_frequency, bin_counter, avg_frequency);

                    avg_probability = sum_probability / number_of_well_in_bin;
                    gsl_matrix_set (binned_probability, bin_counter, scan_counter, avg_probability);
                }

            }
        }
    }


    for (bin_counter = 0; bin_counter < averaged_probability_SMS_bin_number; bin_counter++)
    {
        if (gsl_vector_get(binned_frequency, bin_counter) !=0 )
                file << gsl_vector_get(binned_frequency, bin_counter);
        for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++)
             if (gsl_vector_get(binned_frequency, bin_counter) !=0 )
                 file << "\t" << gsl_matrix_get (binned_probability, bin_counter, scan_counter);

        if (gsl_vector_get(binned_frequency, bin_counter) !=0 )
            file << endl;
    }


    gsl_vector_free (binned_frequency);
    gsl_matrix_free (binned_probability);
    gsl_matrix_free (bin_number);

    file.close();

}

//************************************writting the probabilities of each well after frequency scan****************************//
        //I beleive that I could some how use previous function, but I chose not. Since I want to have these two functions separately for the time being.
void WriteProbabilitiesAfterFrequencyScan (gsl_matrix * probabilities_in_each_frequency_step)
		{
			fstream file;
			file.open("probabilities_after_each_frequency_scan_sms.txt", ios::out | ios::app);
			if (!file) {
				cout << "Error opening the file probabilities_after_each_frequency_scan_sms.txt";
			}
    file << "Well number";
    
    for (int i = 0; i < ensemble_number; i++)
    {
        for (int k = 0; k < number_of_scan; k++)
        {
            for (int j = 0; j < well_numbers; j++)
            {
                file << "\tProb_Syst_" << i << "_Scan_" << k << "_Freq_" << j;//  So there is no negative value for wells which are not in resonance with the laser.
            }
        }

    }
    
    file << endl;
    int variable_burn_frequency_counter;

    for (int i = 0; i < well_numbers; ++i){// it scans over all wells
        variable_burn_frequency_counter = 0;
        file << i;
        for (int ii = 0; ii < ensemble_number; ii++)
        {
            for (int k = 0; k < number_of_scan; k++)
            {
                for (int j = 0; j < well_numbers; j++)
                {// It scans over all frequencies, however, in SMS-V05, the procedure of finding resonant well has been changed in such a way that in frequency scanning loop each well is burned once only so number_of_frequency_steps has been replaced with dimesnion.

                        file << "\t" << gsl_matrix_get (probabilities_in_each_frequency_step,i, variable_burn_frequency_counter);
                        ++variable_burn_frequency_counter;

                }
            }
        }
        file << endl;
    }
    file.close();

}






//************************************writting initial well number****************************//
        void WriteInitialWellNumber (gsl_vector * initial_well)
		{

			fstream file;
			file.open("initial_well_number.txt", ios::out | ios::app);// it keeps initial well number in each ensemble.
			if (!file) {
				cout << "Error opening the file initial_well_number.txt";
			}
    file << "ensemble" << "\tinitial_well" << endl;
    for (int i = 0; i < ensemble_number; i++)
        file << i << "\t" << gsl_vector_get (initial_well, i) << endl;//


}
/*writing the general information of each run in a file */
void WriteGeneralInformation (double calculation_time, double bin_step, double hole_burning_yield)
{

	fstream file;
	file.open("General_information.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file General_information.txt";
	}
    file << "Ensemeble number: " << ensemble_number << "\nReal Well numbers: " << well_numbers << "\nRandom numbers: " << well_numbers
            << "\nBurning frequency (Wavelength): " << burn_frequency << " GHz" << "\t(" << light_speed/burn_frequency << "  nm)" << "\nBurn Temperature: " << burn_temperature
            << " K" << "\nburn time: " << burn_time << " ns" << "\nmean value and deviation of potential barriers:\t"<< "mu = " << mu << "\t sigma = " << sigma << endl;// to write the basic information of the problem.

    file << "Total number of bins: " << total_bins << "\tThe width of bin: " << bin_step ;
    if (calculation_time < 60)
        file << "\nCalculation time: " << calculation_time << " seconds" << endl;
    else if (calculation_time < 3600)
        file << "\nCalculation time: " << calculation_time/60 << " minutes" << endl;
    else
        file << "\nCalculation time: " << calculation_time/3600 << " hours" << endl;

    file << "The avargae hole burning yield is " << hole_burning_yield << endl;
    file << "The cooling time is " << whole_cooling_time/(60*60*1e6) << " hours" << endl;
    file << "The number of burn is " << number_of_burn << " and total burning time is " << number_of_burn * time_interval_between_two_consecutive_burns
            << " sec or " << number_of_burn * time_interval_between_two_consecutive_burns/60 << " min." << endl;
    file << "The spontaneous recovery time is " << whole_recovery_time/(60*60*1e6) << " hours" << endl;
    file << "The curvature of parabola is " << curvature_of_parabola/energy_conversion_factor << " GHz" << endl;
    file.close();

}


//***********writing the probability of occupation of each well at different temperatures***********//
        void WriteRateAtDifferentTemp (gsl_matrix* output_at_different_T){
			fstream file;
			file.open("output_at_different_T.txt", ios::out | ios::app);
			if (!file) {
				cout << "Error opening the file output_at_different_T.txt";
			}
    file << "WellNumber";
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
        for (double temperature = initial_temperature; temperature >= final_temperature; temperature -= temperature_step) // writting the first row of this file
            file << "\tensemble_" << ensemble <<"_at_" << temperature;
    file << "\n";

    for ( int i = 0; i < well_numbers; i++)
    {
        file << i << "\t";
        for ( int j = 0; j< ensemble_number * (int ((initial_temperature - final_temperature)/temperature_step) + 1); j++)
            file << gsl_matrix_get (output_at_different_T, i, j) << "\t";
        file << "\n";
    }

    file.close();

}





//writing the rate elements to find out the reason that out_put_at_different_T has some negative elements.
void WriteRateElements (gsl_matrix* rate, double temperature){

    /*to see each element of rate matrix. Actually rate matrix is Q and different elements of following file is r01, r12,... while rate_matrix consists of a combination of those elements.*/
	fstream file;
	file.open("rate_elements.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file rate_elements.txt";
	}
    file << "Rate elements at T = " << temperature << " K" << endl;

    for (int i = 0; i < well_numbers-1; i++)// the first well does not have any tunneling rate so v_potential_of_well starts from i+1 or 1 for the first well. Also the tunneling rate for the last well is meaningless thus this counter works over 0 to well_numbers-1

        file << "r" << i << i+1 << " = " << gsl_matrix_get (rate,i, i+1) << "\tr" << i+1 << i << " = " << gsl_matrix_get (rate,i+1, i) << endl;
    file << "\n";

    //file.close();//NOTICE: Since I am useing append, I am not sure about closing this file.

}


//write the probability distribution after a long time or boltzman distribution.

void WriteLongtimeDist ( gsl_matrix * longtime_distribution/*or boltzman distribution*/, string filename){

    fstream file;
    file.open(filename, ios::out | ios::app);// the file name is chosen in the main()
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "Wellnumber" ;
    for (int column_counter = 0; column_counter < ensemble_number; column_counter++)
        file << "\tensemble " << column_counter;
    file <<"\n";

    for (int row_counter = 0; row_counter < well_numbers; row_counter++)
    {
        file << row_counter;
        for (int column_counter = 0; column_counter < ensemble_number; column_counter++)
        {
            /*double normalization_factor = 0;//TODO, check its application. I don't remember its purpose but it seems that it is useless now?!!
for (int i = 0; i < well_numbers; ++i)
    normalization_factor += gsl_matrix_get (longtime_distribution, i, column_counter);*/

            file << "\t" << gsl_matrix_get (longtime_distribution, row_counter, column_counter) ;
        }
        file << "\n";

    }
    file.close();

}

//write the initial well and the most probable well in terms of well condition (depth).


void WriteInitialwellsAndMostProbableWells (gsl_matrix * boltzman_distribution_BT, gsl_vector * initial_well, gsl_matrix * last_rates, gsl_vector * deepest_well){
	fstream file;
	file.open("InitialWell-MostProbableWell-Comparison.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file InitialWell-MostProbableWell-Comparison.txt";
	}
    file << "ensemble\tResonantWell\tProbabilityInRW\tMostProbableWell\tProbabilityInMPW_BT\tCoincidenceWellCriteria" << endl;//ProbabilityInRW: the pre-burn probability in resonant well. ProbabilityInMPW: the probability in the deepst well (the most probable one). CoincidenceWellCriteria: If the resonant well and the most probable one are the same this number is 0, otherwise -1.
    int mostProbableWell, resonantWell;

    for (int i = 0; i < ensemble_number; i++)
    {
        resonantWell = (int)gsl_vector_get(initial_well,i);
        mostProbableWell = (int)gsl_vector_get(deepest_well,i);
        file << i << "\t" << resonantWell << "\t" << gsl_matrix_get (last_rates, resonantWell, i)
                << "\t" << mostProbableWell << "\t" << gsl_matrix_get (boltzman_distribution_BT, mostProbableWell, i);
        if (mostProbableWell == resonantWell)file << "\t" << 0;
        else file << "\t" << -1;// I am intersted to know whether the most probable well is in resonant with laser or not. 0 shows that it is and -1 shows that it is not. Then I can count the number of system which the most probable well are/are not in resonant.
        file << endl;
    }

    file.close();

}






//write the discrepancy between different situations and equilibrium.

void WriteDiscrepancyOfEquilibrium (gsl_matrix * boltzman_distribution, gsl_matrix * output, string filename){
    fstream file;
    file.open(filename, ios::out | ios::app);// the file name is chosen in the main()
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "WellNumber";
    double difference;
    for (int ensemble = 0 ; ensemble < ensemble_number; ensemble++)
        file << "\tensemble-"<< ensemble;
    file << "\n";

    for (int i = 0; i < well_numbers; i++)
    {
        file << i ;

        for (int ensemble = 0 ; ensemble < ensemble_number; ensemble++)
        {
            difference = gsl_matrix_get(boltzman_distribution, i, ensemble)-gsl_matrix_get(output, i, ensemble);
            file << "\t" << std::abs(difference);//the absolute value of the difference between the probability of each well in each situation and boltzman for corresponding well.
        }
        file << "\n";
    }
    file.close();

}


//write the probability of resonant well (frequency)

void WriteProbabilityHeight (gsl_vector * probability_height, double post_burn_probability, string filename){
    fstream file;
    file.open(filename, ios::out | ios::app);// the file name is chosen in the main()
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "RecoveryTime(s)" << "\tProbabilityHight" << endl;
    file << 0 << "\t" << post_burn_probability << endl;// this is post burn probability at burn frequency.

    for (int i = 0; i < number_of_recovery_measurement; i++)
        file << (recovery_time * (i+1))/1e9 << "\t" << gsl_vector_get (probability_height, i) << endl;
    file.close();
}


// write the binned spectrum after each step of recovery
void WriteRecoveredSpectrum (gsl_matrix * recovered_probability, double smallest_frequency, double bin_step, string filename){

    fstream file;
    file.open(filename, ios::out | ios::app);// the file name is chosen in the main().
    if(!file) {
		cout << "Error opening the file" << filename;
	}
    file << "Frequency_b";

    for (int i = 0; i < number_of_recovery_measurement; i++)
        file << "\tRecoverNumb" << i;
    file << endl;

    for (int i = 0; i <= total_bins; i++)
    {
        file << smallest_frequency + i * bin_step;
        for (int j = 0; j < number_of_recovery_measurement; j++)
            file << "\t" << gsl_matrix_get (recovered_probability, i, j);

        file << "\n";

    }
    file.close();

}




// write the difference between pre burn (or even post bur) and equilibrium.

void WriteEquilibriumDifference (gsl_vector * equilibrium_ensemble, string filename){

    fstream file;
    file.open(filename, ios::out | ios::app);// the file name is chosen in the main().
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "Ensemble" << "\tAverage-diff" << endl;
    double averaged_difference = 0;
    for (int i = 0; i < ensemble_number; i++)
    {
        file << i << "\t" << gsl_vector_get (equilibrium_ensemble, i) << endl;
        averaged_difference += gsl_vector_get (equilibrium_ensemble, i);
    }

    file.close();

}




void WriteBurningYieldvsStartingWell (gsl_vector * initial_well, gsl_vector * hole_burning_yield_starting_well, gsl_vector * pre_burn_at_starting_well,
                                      gsl_matrix * energy_difference) {

	fstream file;
	file.open("initial well vs yield probability.txt", ios::out| ios::app);
	if (!file) {
		cout << "Error opening the file Initial well vs yield probability.txt";
	}
    file << "ensemble" << "\tinitial-well" << "\tBurning-Yield" << "\tPreBurn-Probability" << "\tFrequency-starting-well"
            "\tenergy_difference-b-e"<< endl;

    int starting_well;// this is not necessary, only for clarification.
    double difference;// this is not necessary, only for clarification.
    for (int i = 0; i < ensemble_number; i++)
    {
        starting_well = (int)gsl_vector_get (initial_well, i);
        difference = gsl_matrix_get(energy_difference, starting_well, i);
        file << i << "\t" << starting_well << "\t" << gsl_vector_get (hole_burning_yield_starting_well, i) << "\t" <<
                gsl_vector_get (pre_burn_at_starting_well, i)<< "\t" << difference <<
                "\t" << burn_frequency - difference << endl;
    }
    //file.close();//NOTICE: Since I am useing append, I am not sure about closing this file.

}



//*********writing the hole growth kinetics****************************//

        void WriteHoleGrowthKinetics (gsl_matrix * hole_growth_kinetics, gsl_vector * pre_burn_at_starting_well){

            //Hole_growth_kinetics is a matrix with number of burn steps rows and ensemble coulumns, and it contains A0 probabilities for resonant well
            //pre_burn_at_starting_well is a vector containing starting value of A0 in starting well. Not used as of 2022.
            // instead sum_resonant_probabilities_preburn is calculated outside, it is global variable.
			fstream file;
			file.open("HoleGrowthKinetics-time.txt", ios::out | ios::app);
			if (!file) {
				cout << "Error opening the file HoleGrowthKinetics-time.txt";
			}
    file << "Time-sec" << "\tDelta-F" << "\tProbability-RW" << endl;



	fstream file2;
	file2.open("HoleGrowthKinetics-dose.txt", ios::out | ios::app);
	if (!file2) {
		cout << "Error opening the file HoleGrowthKinetics-dose.txt";
	}
    //file2.setRealNumberPrecision(15);

    file2 << "dose" << "\tDelta-F" << "\tProbability" << "\tActsOfBurn "<< "\tDelta-F2" << endl;
    // delta-f is the result of (post-burn - pre-burn)/pre-burn

     // first data-lines in two files, corresponding to pre-burn situation
    //file << "Time-sec" << "\tDelta-F" << "\tProbability-RW"
    file << burn_time * 1e-6/100 << "\t" << 1 << "\t" << sum_resonant_probabilities_preburn << endl;
    //burn_time * 1e-6/100 is an aribitrary value for pre-burn.

    // burn_power is in mW
    file2 << (time_interval_between_two_consecutive_burns/100)*burn_power*1e-3/(pi*(beam_diameter/2)*(beam_diameter/2)) << "\t" << 1
            << "\t" << sum_resonant_probabilities_preburn << "\t" << 1e-6 /*small number instead of 0*/ << "\t" << 1 << endl;
    // dose is in J/cm^2. time_interval_between_two_consecutive_burns/100 is an arbitrary value considered as pre-burn dose since I don't want to use 0, 
    // as we like to report HGK on a logarithmic scale. 
    //time_interval_between_two_consecutive_burns is in second.

    double sum = 0, dose = 0, delta_fluorescence = 0;
    if (number_of_HGK_measurement != 1 && burn_in_one_step != 1)
    {
        for (int i = 0 ; i < number_of_HGK_measurement; i++)
        {
            sum = 0;

            for (int j = 0; j < ensemble_number; j++)
            { 
                sum += gsl_matrix_get (hole_growth_kinetics, i, j); // sum over molecules
            }

            delta_fluorescence = -(sum - sum_resonant_probabilities_preburn)/ sum_resonant_probabilities_preburn;
            // according to definition, this number should be negative since it is post-burn - pre-burn, 
            //however M would like tit to be positive since its drwaing makes more sense with positive number.
            file << time_interval_between_two_consecutive_burns * ((i+1)*HGK_measurement_frequency) << "\t" << 1-delta_fluorescence << "\t" << sum << endl;
            // NB: Mehdi had wrong time here, he had burn time which is excited state lifetime

            dose = time_interval_between_two_consecutive_burns*((i+1)*HGK_measurement_frequency)*burn_power*1e-3/(pi*(beam_diameter/2)*(beam_diameter/2));
            //burn_power should be in W, the conversion factor from mW to W is 1e-3.

            file2 << dose << "\t" << 1-delta_fluorescence << "\t" << sum << "\t" << (i+1)*HGK_measurement_frequency << "\t" << 1-delta_fluorescence << endl;

        }
    }
    else
    {// in case that HGK is measured once only, it should be the result of last burning not the first one. so instead of "i" in preious loop I substitute it by 1 and HGK_measurement_frequency with number_of_burn. For instance, number_of_burn = 10 and number_of_HGK_measurement = 10 means HGK_measurement_frequency = 1, so (i*HGK_measurement_frequency+1) = 9*1+1 = 10. To get the same result for number_of_HGK_measurement =  1 (HGK_measurement_frequency = 10), it should be 1*(number_of_burn) = 1*10
        // for time, by the same comparison I found following equation.

        sum = 0;// unnecessary since double sum = 0.Just to make sure that it does not have any value.

        for (int j = 0; j < ensemble_number; j++)
            sum += gsl_matrix_get (hole_growth_kinetics, 0, j);//hole_growth_kinetics has only 1 row in this situation.


        delta_fluorescence = -(sum - sum_resonant_probabilities_preburn)/ sum_resonant_probabilities_preburn;// according to definition, this number should be negative since it is post-burn - pre-burn, however I would like to find it positive since its drwaing makes more sense with positive number.
        file << burn_time*1e-6 * (number_of_burn) << "\t" << 1-delta_fluorescence << "\t" << sum << endl; // to convert time from micro-sec to sec, it should be multiplied by 1e-6.

        dose = (time_interval_between_two_consecutive_burns*number_of_burn)*burn_power*1e-3/(pi*(beam_diameter/2)*(beam_diameter/2));//burn_power should be in W, the conversion factor from mW to W is 1e-3.

        file2 << dose << "\t" << 1-delta_fluorescence << "\t" << sum << "\t" << number_of_burn << "\t" << 1-delta_fluorescence << endl;

    }

    file.close();//NOTICE: Since I am using append, I am not sure about closing this file.
    file2.close();//NOTICE: Since I am using append, I am not sure about closing this file.
}



//***********writing the probability of occupation of each well at different recovery time***********//
        void WriteRateAtDifferentTime (gsl_matrix* output_at_different_recovery_time, const double number_of_measurement, string filename){



    fstream file;
    file.open(filename, ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "WellNumber";
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
        for (int recovery_counter = 0; recovery_counter < number_of_measurement; ++recovery_counter) // writting the first row of this file
            file << "\tensemble_" << ensemble <<"_at_" << recovery_counter;
    file << "\n";

    for ( int i = 0; i < well_numbers; i++)
    {
        file << i << "\t";
        for ( int j = 0; j< ensemble_number * number_of_measurement; j++)
            file << gsl_matrix_get (output_at_different_recovery_time, i, j) << "\t";
        file << "\n";
    }


    file.close();

}





//*************writing the ratio of probability of each well after each recovery step with respect to the boltzamn distribution in corresponding well.

void WriteRatioProbability (gsl_matrix * longtime_distribution_BT, gsl_matrix * output_at_different_recovery_time, double post_burn_probability){
	fstream file;
	file.open("ratio.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file ratio.txt";
	}
    file << "Time";
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
        for (int i = 0; i < well_numbers; i++)
            file << "\tEnsemble" << ensemble << "_WellNumber" << i;
    file << endl;
    file << 0 << "\t" << post_burn_probability << endl;// this is post burn probability at burn frequency.

    for (int i = 0; i < number_of_recovery_measurement; i++)
    {
        file << (recovery_time * (i+1))/1e9 << "\t";

        for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
            for (int i = 0; i < well_numbers; ++i)
                file << gsl_matrix_get(output_at_different_recovery_time, i, ensemble)/gsl_matrix_get (longtime_distribution_BT,i ,ensemble) << "\t" ;

        file << endl;

    }
    file.close();

}

void WriteProbabilityRecoveredWell(gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well,
                                   gsl_matrix * last_rates_preburn, gsl_matrix *last_rates, const int number_of_measurement, string filename){


    fstream file;
    file.open(filename, ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "System" << "\tInitialWell" << "\tpreburnProb" << "\tpostburnProb";
    for (int recovery_counter = 0; recovery_counter < number_of_measurement; ++recovery_counter)
        file <<"\tProb_after_" << recovery_counter << "_recovery";
    file<< endl;
    int initialwell;
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        file << ensemble << "\t" << initialwell << "\t" << gsl_matrix_get(last_rates_preburn, initialwell, ensemble) << "\t" <<
                gsl_matrix_get(last_rates, initialwell, ensemble);
        for (int recovery_counter = 0; recovery_counter < number_of_measurement; recovery_counter++)
            file << "\t" << gsl_matrix_get (output_at_different_recovery_time, initialwell, ensemble*number_of_measurement+recovery_counter);
        file<< endl;

    }
    file.close();

}

// burned well probability is written for each system at different time
void WriteProbabilityInBurnedWell (gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well,
                                   gsl_matrix * last_rates_preburn, gsl_matrix *last_rates){

	fstream file;
	file.open("ProbabilitiesInBurnedWell.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file ProbabilitiesInBurnedWell.txt";
	}
    file << "RecoveryTime";
    int initialwell;
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        file << "\tW0-1-Syst-"<< initialwell-1 << "-" << ensemble << "\tW0-Syst-" << initialwell << "-" << ensemble << "\tW0+1-Syst-" << initialwell+1 << "-" << ensemble
                << "\tW11-Syst-" << ensemble;
    }
    file << endl;

    file << -1;//preburn
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        file << "\t" << gsl_matrix_get(last_rates_preburn, initialwell-1, ensemble)<< "\t" << gsl_matrix_get(last_rates_preburn, initialwell, ensemble)
                << "\t" << gsl_matrix_get(last_rates_preburn, initialwell+1, ensemble) << "\t" << gsl_matrix_get(last_rates_preburn, 11, ensemble);
    }
    file << endl;
    file << 0;/*initial time which is 0 with the probability of post burn in each burned well*/
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        file << "\t" << gsl_matrix_get(last_rates, initialwell-1, ensemble)<< "\t" << gsl_matrix_get(last_rates, initialwell, ensemble)
                << "\t" << gsl_matrix_get(last_rates, initialwell+1, ensemble) << "\t" << gsl_matrix_get(last_rates, 11, ensemble);
    }
    file << endl;
    for (int recovery_counter = 0; recovery_counter < number_of_recovery_measurement; recovery_counter++)
    {
        file << (recovery_time * (recovery_counter+1))/1e9 ;
        for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
        {
            initialwell = (int)gsl_vector_get(initial_well, ensemble);
            file << "\t" << gsl_matrix_get (output_at_different_recovery_time, initialwell-1, ensemble*number_of_recovery_measurement + recovery_counter)
                    << "\t" << gsl_matrix_get (output_at_different_recovery_time, initialwell, ensemble*number_of_recovery_measurement + recovery_counter)
                    << "\t" << gsl_matrix_get (output_at_different_recovery_time, initialwell+1, ensemble*number_of_recovery_measurement + recovery_counter)
                    << "\t" << gsl_matrix_get (output_at_different_recovery_time, 11, ensemble*number_of_recovery_measurement + recovery_counter);
        }
        file << endl;
    }
    file.close();

}

void WriteProbabilityInOnlyBurnedWell(gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well,
                                      gsl_matrix * last_rates_preburn, gsl_matrix *last_rates_postburn, const int number_of_measurement
                                      , const double time, string filename){
    fstream file;
    file.open(filename, ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file" << filename;
	}

    file << "RecoveryTime";
    int initialwell;
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        file << "\tW0-Syst-" << initialwell << "-" << ensemble;
    }

    file << endl;

    file << -1;//preburn
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        file << "\t" << gsl_matrix_get(last_rates_preburn, initialwell, ensemble);
    }
    file << endl;
    file << 0;/*initial time which is 0 with the probability of post burn in each burned well*/
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        file << "\t" << gsl_matrix_get(last_rates_postburn, initialwell, ensemble);
    }
    file << endl;
    for (int recovery_counter = 0; recovery_counter < number_of_measurement; recovery_counter++)
    {
        file << (time * (recovery_counter+1))/1e9 ;
        for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
        {
            initialwell = (int)gsl_vector_get(initial_well, ensemble);
            file << "\t" << gsl_matrix_get (output_at_different_recovery_time, initialwell, ensemble*number_of_measurement + recovery_counter);
        }
        file << endl;
    }

    file.close();

}

//writing the relative hole depth in a file. By now, I don't know a way to use some part of next function. Since it is used in thermocycling process as well, the result of relative hole depth will be overwritten which is not my goal.

void WriteRelativeHoleDepthRecovery (gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well,
                                     gsl_matrix * last_rates_preburn, gsl_matrix *last_rates_postburn){
	fstream file;
	file.open("Recovery-RelativeDepth.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file Recovery-RelativeDepth.txt";
	}
    file << "RecoveryTime-sec\t" << "Probability\t" << "PreburnSummation\t" << "PostburnSummation\t"
            << "HoleDepth\t" << "SummationOverResonantWellProb" << endl;
    int initialwell;
    double pre_burn_summation = 0;
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        pre_burn_summation += gsl_matrix_get(last_rates_preburn, initialwell, ensemble);//weight_factor * gsl_matrix_get(last_rates_preburn, initialwell, ensemble);
    }

    double post_burn_summation = 0;

    file << number_of_burn*burn_time/1e6 << "\t" << 1 << endl;// initial time is chosen post burn time in second (1e6 is for micro second) which is number_of_burn*burn_time/1e6. The depth of hole after burning is 1 and during recovery it reduces to smaller value (0, for 100% recovery).

    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        post_burn_summation += gsl_matrix_get(last_rates_postburn, initialwell, ensemble);
    }


    double hole_depth = pre_burn_summation - post_burn_summation;

    std::cout << "There is a " << (hole_depth / pre_burn_summation) * 100 << "% hole in this simulation.\n";

	fstream file2;
	file2.open("InitialWell-Recovery.txt", ios::out | ios::app);//debugging
	if (!file2) {
		cout << "Error opening the file InitialWell-Recovery.txt";
	}

file2 << "NoRecovery\tEnsemble\tInitialWell\tProbValue" << endl;
    double summation;
    double probability_in_burned_well;

    if (hole_depth == 0){
        for (int recovery_counter = 0; recovery_counter < number_of_recovery_measurement; recovery_counter++)
		{
            file <<(recovery_time * (recovery_counter+1))/1e6 ;//1e6 for micro second.
            file << "\t" << 1 << endl;// the hole depth in recovey remains constantly 1 if there is no hole.
        }
    }
	else
	{

        for (int recovery_counter = 0; recovery_counter < number_of_recovery_measurement; recovery_counter++)
		{
            file <<(recovery_time * (recovery_counter+1))/1e6 ;//1e6 for micro second.
            file2 << recovery_counter;
            summation = 0;
            for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
            {
                initialwell = (int)gsl_vector_get(initial_well, ensemble);
                probability_in_burned_well = gsl_matrix_get (output_at_different_recovery_time, initialwell, ensemble * number_of_recovery_measurement + recovery_counter);
                summation += probability_in_burned_well;//weight_factor * probability_in_burned_well;
                file2 << "\t" << ensemble << "\t" << initialwell << "\t" << ensemble * number_of_recovery_measurement + recovery_counter<< "\t"
                        << gsl_matrix_get (output_at_different_recovery_time, initialwell, ensemble * number_of_recovery_measurement + recovery_counter) << endl;//debugging
            }
            file << "\t" << (pre_burn_summation-summation) / hole_depth
                    << "\t"  << pre_burn_summation << "\t"  << post_burn_summation << "\t"  << hole_depth
                    << "\t"  << summation << endl;// the hole depth in recovey should be normalized to the hole depth just after burning.
        }
    }


    file.close();

}
//
//
#include "load_experimental_thermocycling_profile.h"
void WriteAveragedProbabilityInOnlyBurnedWell(gsl_matrix * output_at_different_recovery_time, gsl_vector * initial_well,
                                              gsl_matrix * last_rates_preburn, gsl_matrix *last_rates_postburn, const int number_of_measurement,
                                              const double time, string filename)
{

    fstream file;
    file.open(filename, ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    gsl_matrix * experimental_thermocycling_time_temperature = gsl_matrix_calloc(number_of_thermocycling_measurement,2);
	// first column keeps the cycling time (increasing T to maximal temperature and cooling it down to burn temperature), and the second column keeps the corresponing maximal cycling temperature.

    //LoadExperimentalThermocyclingProfile(experimental_thermocycling_time_temperature);
    // commented out in December 2022, 
    //TO DO : need to put another thermocycling profile resder here


    file << "RecoveryTime(sec)\t" ;// note this may not right for thermocycling, I think it should be replaced with temperature.

    if (filename == "AveragedProbabilityInOnlyBurnedWell_Thermocycling.txt")
        file << "Temperature" << "\tProbability" << endl;
    else  file << "Probability" << endl;
    file << burn_time * 1e-6/100;//preburn. burn_time * 1e-6/100 is an arbitrary number presenting the preburn situation
    int initialwell;
    double pre_burn_summation = 0;


    //double weight_factor;// weight_factor is the pre-burn probability in each system which determins the weight of that well in ensemble.

    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        //weight_factor = gsl_matrix_get(last_rates_preburn, initialwell, ensemble);
        pre_burn_summation += gsl_matrix_get(last_rates_preburn, initialwell, ensemble);
        //weight_factor * gsl_matrix_get(last_rates_preburn, initialwell, ensemble);
    }
    if (filename == "AveragedProbabilityInOnlyBurnedWell_Thermocycling.txt")file << "\t" << burn_temperature << "\t" << pre_burn_summation << endl;// the summation over preburn probabilities of burned wells in ensemble.
    else file << "\t" << pre_burn_summation << endl;
    // the summation over preburn probabilities of burned wells in ensemble.


    double post_burn_summation = 0;
    file << number_of_burn*burn_time;// initial time is chosen post burn time which is number_of_burn*burn_time

    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        //weight_factor = gsl_matrix_get(last_rates_preburn, initialwell, ensemble);
        post_burn_summation += gsl_matrix_get(last_rates_postburn, initialwell, ensemble);
        //weight_factor * gsl_matrix_get(last_rates_postburn, initialwell, ensemble);
    }

    if (filename == "AveragedProbabilityInOnlyBurnedWell_Thermocycling.txt")
        file << "\t" << burn_temperature << "\t" << post_burn_summation << endl;
    // the summation over preburn probabilities of burned wells in ensemble.
    else file << "\t" << post_burn_summation << endl;// the summation over preburn probabilities of burned wells in ensemble.




    double summation;
    double probability_in_burned_well;
double time_buck = 0;
    for (int recovery_counter = 0; recovery_counter < number_of_measurement; recovery_counter++)
    {
        if (filename == "AveragedProbabilityInOnlyBurnedWell_Thermocycling.txt")
        {
            time_buck+=gsl_matrix_get (experimental_thermocycling_time_temperature, recovery_counter, 0)/1e6;
            file << time_buck ;

         }
        else
            file << (time * (recovery_counter+1))/1e6 ;
        //in second. 1e6 for micro second.NB 1e6 should be substituted by another variable in constancts.h

        if (filename == "AveragedProbabilityInOnlyBurnedWell_Thermocycling.txt")
        {
            //double t = time * recovery_counter * 1e-6/60;// time in following equation should be in minute.
            file << "\t" << gsl_matrix_get (experimental_thermocycling_time_temperature, recovery_counter, 1);
        }
        summation = 0;
        for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
        {
            initialwell = (int)gsl_vector_get(initial_well, ensemble);
            //weight_factor = gsl_matrix_get(last_rates_preburn, initialwell, ensemble);
            probability_in_burned_well = gsl_matrix_get (output_at_different_recovery_time, initialwell, ensemble * number_of_measurement + recovery_counter);
            summation += probability_in_burned_well;//weight_factor * probability_in_burned_well;
        }
        file << "\t" << summation << endl;

    }
    file.close();

}


void WriteAveragedProbabilityInDeepestWell(gsl_matrix * output_at_different_recovery_time,
                                           gsl_matrix * last_rates_preburn, gsl_matrix *last_rates_postburn, const int number_of_measurement,
                                           const double time, gsl_vector * deepest_well, gsl_vector * initial_well, string filename){

    fstream file;
    file.open(filename, ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "RecoveryTime\t" << "ProbabilityRW\t" << "ProbabilityDW" << endl;
    file << burn_time * 1e-6/100;//preburn
    int resonantwell, deepestwell;

    double summation1 = 0, summation2 = 0, summation3 = 0, summation4 = 0;


    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        resonantwell = (int)gsl_vector_get(initial_well, ensemble);
        deepestwell = (int)gsl_vector_get(deepest_well, ensemble);

        summation1 += gsl_matrix_get(last_rates_preburn, resonantwell, ensemble);
        summation2 += gsl_matrix_get(last_rates_postburn, resonantwell, ensemble);
        summation3 += gsl_matrix_get(last_rates_preburn, deepestwell, ensemble);
        summation4 += gsl_matrix_get(last_rates_postburn, deepestwell, ensemble);
    }
    file << "\t" << summation1 << "\t" << summation3 << endl;// the summation over preburn probabilities of resonant and deepest wells in ensemble.


    file << 0;/*initial time which is 0 with the summation over postburn probabilities of burned wells in ensemble.*/

    file << "\t" << summation2 << "\t" << summation4 << endl;// the summation over postburn probabilities of resonant and deepest wells in ensemble.


    for (int recovery_counter = 0; recovery_counter < number_of_measurement; recovery_counter++)
    {
        file << (time * (recovery_counter+1))/1e6 ;
        summation2 = 0;
        summation4 = 0;
        double probability_in_deepest_well, probability_in_resonant_well;
        for (int ensemble = 0; ensemble < ensemble_number; ++ensemble){
            resonantwell = (int)gsl_vector_get(initial_well, ensemble);
            deepestwell = (int)gsl_vector_get(deepest_well, ensemble);

            probability_in_resonant_well = gsl_matrix_get (output_at_different_recovery_time, resonantwell, ensemble * number_of_measurement + recovery_counter);
            probability_in_deepest_well = gsl_matrix_get (output_at_different_recovery_time, deepestwell, ensemble * number_of_measurement + recovery_counter);
            summation2 += probability_in_resonant_well;
            summation4 += probability_in_deepest_well;
        }
        file << "\t" << summation2 << "\t" << summation4 << endl;
    }

    file.close();

}



void WriteProbabilitiesAtDifferentT (gsl_matrix * output_at_different_T, gsl_vector * initial_well){

	fstream file;
	file.open("ProbabilitiesInDifferentT.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file ProbabilitiesInDifferentT.txt";
	}
    file << "Temperature";

    int initialwell;
    for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
    {
        initialwell = (int)gsl_vector_get(initial_well, ensemble);
        file << "\tW0-1-Syst-"<< initialwell-1 << "-" << ensemble << "\tW0-Syst-" << initialwell << "-" << ensemble << "\tW0+1-Syst-" << initialwell+1 << "-" << ensemble
                << "\tW11-Syst-" << ensemble;
    }
    file << endl;

    int temperature_counter = 0;
    int number_of_cooling = (int ((initial_temperature - final_temperature)/temperature_step) + 1);
    int column_counter;
    for (double temperature = initial_temperature; temperature >= final_temperature; temperature -= temperature_step)
    {
        file << temperature;

        for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
        {
            initialwell = (int)gsl_vector_get(initial_well, ensemble);
            column_counter = number_of_cooling*ensemble + temperature_counter;
            file << "\t" << gsl_matrix_get(output_at_different_T, initialwell-1, column_counter)<< "\t" << gsl_matrix_get(output_at_different_T, initialwell, column_counter)
                    << "\t" << gsl_matrix_get(output_at_different_T, initialwell+1, column_counter) << "\t" << gsl_matrix_get(output_at_different_T, 11, column_counter);
        }
        file << endl;

        ++temperature_counter;
    }


    file.close();

}





//#include <iostream>
//using namespace std;
void WriteAveragedOutputAtdifferentT (gsl_matrix * output_at_different_T){

	fstream file;
	file.open("AveragedOutputAtDifferentT.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file AveragedOutputAtDifferentT.txt";
	}
    file << "WellNumber";

    for (double temperature = initial_temperature; temperature >= final_temperature; temperature -= temperature_step) // writting the first row of this file
        file << "\tAveraged_at_" << temperature;
    file << "\n";




    int number_of_cooling = (int ((initial_temperature - final_temperature)/temperature_step) + 1);
    int column_counter = 0;
    int column_number;
    double summation;
    for ( int i = 0; i < well_numbers; i++)
    {//for (double temperature = initial_temperature; temperature >= final_temperature; temperature -= temperature_step){
        file << i << "\t";

        for (double temperature = initial_temperature; temperature >= final_temperature; temperature -= temperature_step)
        {
            summation = 0;
            for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
            {

                column_number = number_of_cooling*ensemble + column_counter;
                summation +=  gsl_matrix_get(output_at_different_T, i, column_number);
                //cout << "output in well " << i << " in ensemble " << ensemble << " is " << gsl_matrix_get(output_at_different_T, i, column_number) << endl;

            }
            file << (summation / ensemble_number) << "\t";
            ++column_counter;
        }
        column_counter = 0;
        file << endl;
    }
    file.close();

}



//*************************************writing the translation energy in a file*****************//
        void WriteTranslationEnergy (double translation){

    string file_name;
    file_name.append("translation_energy_");// the file name is chosen in the main() based on the ground or excited state information.
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
fstream file;
file.open(file_name, ios::out | ios::app);
if (!file) {
	cout << "Error opening the file" << file_name;
}
    //file.setRealNumberPrecision(45);

    file << translation << endl;

    //file.close();//NOTICE: Since I am useing append, I am not sure about closing this file.
}



void WriteAllConstants ()
{ //this is a dangerous exercise. not sure why we are doing this

	fstream file;
	file.open("constants.h", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file constants.h";
	}

    file <<"#ifndef CONSTANTS_H" << endl;
    file <<"#define CONSTANTS_H" << endl;
    file <<"\n";
    file <<"#include <cmath>" << endl;
    file <<"\n";
    file <<"const double pi = " << pi <<";" << endl;
    file <<"const double round_off = " << round_off <<";// as simulation time and step are defined as double, the round off issue makes problem and the last step doesn't work." << endl;
    file <<"const double round_off_limit = " << round_off_limit << "; // I may use this limit to set some very small numbers to zero." << endl;
    file <<"const double accuracy_measure = "<< accuracy_measure << ";" << endl;
    file <<"const double extreme_small_number = "<< extreme_small_number<<";" << endl;
    file <<"const int balancing = "<< balancing <<";// balancing in solving eigensystem is on (1) or off (0)." << endl;
    file <<"const double energy_conversion_factor = "<< energy_conversion_factor <<";// to convert energy (frequency) from GHz to MHz (by multipilication).// All equation have been drived by assumption that energy is in GHz and time in nano second. However, later I decided to change the scales to MHz and microseconds." << endl;
    file <<"const double time_conversion_factor = "<< time_conversion_factor <<";// to convert time from nano second to micro second (by multipilication).NOTICE: I have to change it to a unit for example micro second or... . In this case I have to make sure that all time values are converted properly (in some cases I multiply time by 1e6 to convert second to microsec).// All equation have been drived by assumption that energy is in GHz and time in nano second. However, later I decided to change the scales to MHz and microseconds." << endl;
    file <<"\n";
    file <<"const double light_speed = "<<   light_speed <<";// m/s, the speed of light in vacuum" << endl;
    file <<"const double h_planck = "<<   h_planck <<";// J.S.  planck's constant" << endl;
    file <<"const double k_b =  20.8366455 * energy_conversion_factor ; // MHz, in GHz over Kelvin or GHz/K. k_b(J/K) / h (J.S) or 1.3806503*1e-23/6.626068*1e-34 => k_b = 2.08366455 * 1e10  in Hz/K or 2.08366455 * 1e1 GHz/K" << endl;
    file <<"\n";
    file <<"const double w_0 =  1e3 * energy_conversion_factor;//MHz. This coeficient is used in W and is equal to 10^12 Hz according to nonexponential hole burning in organic glass p. 4571. (7.6 is added after our discussion on Jan 14, 2011 according to fitting results on experimental data). I changed it to 10^12 Hz again on Feb 8, 2012 (look at my log book volume II)" << endl;
    //file <<"const double hopping_rate_pre_factor = "<<   w_0 <<"; // in GHz. For the time being I assumed that this is the same as w_0 , hoewever this is a pre factor of hopping rate which I have not found a reasonable value for it yet. Then, one could take it from 'U6+ in SrWO4 as a two-level system: a study of the thermal activation of spectral hole refilling, Journal of Luminescence 94–95 (2001) 683–686'" << endl;
    file <<"const double mass = "<<   mass <<";// decreased from 100 amu to 10 amu (look at my log book volume II, Feb 8, 2012) and one amu is 1.66* 10^-27 kg" << endl;
    file <<"\n";
    file <<"const double interval = "<<   interval <<";// if interval is defined as integer, the program does not work properly." << endl;
    file <<"const double curvature_of_parabola = "<<   curvature_of_parabola <<";//75740.5 * energy_conversion_factor;//75740.5 * energy_conversion_factor ;// MHz, the curvature of parabola. This value results in having smallest delta equals to 10 cm^-1 or ~300 GHz. This value should be chamged to have the same assymetry between two wells at the bottom of two wells if stretch is changing." << endl;
    file <<"//const double anharmonicity =  5 * mu;// The coefficient of anharmonicity term" << endl;
    file <<"const int well_numbers = "<<   well_numbers <<";// number of random numbers or teeth which will be produced in the landscape and it is correspond to well numbers. it should be chosen odd to make equal wells around the center." << endl;
    file <<"const int point_numbers =  well_numbers*2; /* the points which are needed for plotting the main curve (parabola). the distance between two adjacent points is 2*interval/(point_numbers-1)" << endl;
    file <<"                               (or (max-min)/(point_numbers-1) if the interval is not symmetrical with recpect to zero(. For example we" << endl;
    file <<"                               want to have 5 points for drawing parabola so each step should be one unit ([2*2]/[5-1]) to span -2 to 2 (our interval in this case) is five steps.*/" << endl;
    file <<"\n";
    file <<"\n";
    file <<"//coupling between the ground and excited state." << endl;
    file <<"const int bottom_decoupling = "<<   bottom_decoupling <<"; //if this value is 1 then a new set of random variable is used for the bottom of wells in the excited state; otherwise, the same value from the ground state is used (by multipling it by stretch)" << endl;
    file <<"const int barrier_decoupling = "<<   barrier_decoupling <<"; //if this value is 1 then a new set of random variable is used for the the barrier heights in the excited state; otherwise, the same value from the ground state is used (by multipling it by stretch)" << endl;
    file <<"\n";
    file <<"\n";
    file <<"\n";
    //file <<"const double burn_wavelength = "<<   burn_wavelength <<"; //nm" << endl;
    file <<"\n";
    file <<"//parameters for distribution of energy gap between the ground and excited states" << endl;
    file <<"const double inhomogeneous_bandwidth = "<<   inhomogeneous_bandwidth <<";// cm^-1. it is used for calculating the standard deviation of translation distribution." << endl;
    file <<"//const double standard_deviation_of_translation = "<<   standard_deviation_of_translation <<";//I changed the mean value from SDF peak to the burn frequency on Feb 13. producing random numbers with gaussian distribution to use in different landscape. 438794.499 and 5995.84916 GHz are corresponding to 14641 and 180 wavenumbers from the SDF peak and width of CP43, A-state, respectively, in the last paper. 180 comes from measuring inhomogeneous line broadening experiments." << endl;
    file <<"const double standard_deviation_of_translation =  (inhomogeneous_bandwidth*29.9792458/2.355)*energy_conversion_factor;// the relationship between bandwidth in cm^-1 and standard deviation is FWHM = 2.355 * sigma. It should be noted that 1 cm^-1 = 29.9792458 GHz then FWHM=180*29.97 = 5394.6 GHz." << endl;
    file <<"const double mean_value_of_translation = (light_speed/burn_wavelength) * energy_conversion_factor;//based on the last paper burn frequency was 14560 wavenumbers or 686.8131 nm and this formula converts it to the GHz (denominator has a 1e-9 for nm to m and the result shoudl be multiplied by 1e-9 to convert total frequency in GHz)." << endl;
    file <<"\n";
    file <<"const double burn_frequency = mean_value_of_translation;//frequency in MHz." << endl;
    file <<"const double burn_power = "<<   burn_power <<";//mW. 1microW is 1e-3, as an example. 5, 10, and 15 W are taken from Kohler's paper." << endl;
    file <<"const double beam_diameter = "<<   beam_diameter <<";//1/(2/sqrt(pi)); //cm. NB: In Kohler's paper the power is given in W/cm^2, so I have to conver the diameter in such a way that it gets rid of other constants in my equation (done on on September 12, 2012, in the second log book). In that calculation, area = pi r^2 = pi d^2/4. by assuming d = 2/sqrt(pi) the area would be 1cm^2. BTW sqrt(pi)/2 = 0.886 cm." << endl;
    file <<"const double absorption_cross_section_at_peak_RT = "<<   absorption_cross_section_at_peak_RT <<";//cm^2// September 12, 2012 (second log book)" << endl;
    file <<"const double homogeneous_width_at_RT = "<<   homogeneous_width_at_RT <<";//cm^-1 for Cp43" << endl;
    file <<"const dhomogeneous_width_at_BTouble  = "<<   homogeneous_width_at_BT <<";//cm^-1 for Cp43" << endl;
    file <<"const double phonon_coupling = "<<   phonon_coupling <<";// this is S or Huang–Rhys factor for Cp43." << endl;
    file <<"const double absorption_cross_section_at_peak_BT = absorption_cross_section_at_peak_RT * (homogeneous_width_at_RT/homogeneous_width_at_BT) * exp(-phonon_coupling);//=1.054e-12;//cm^2, see the notes on September 12, 2012 (in the second log book)." << endl;
    file <<"const double time_interval_between_two_consecutive_burns = 1*(beam_diameter*beam_diameter/(absorption_cross_section_at_peak_BT*burn_power*burn_wavelength))*1.56015e-13;//in seconds, see the notes on September 12, 2012 (in the second log book)." << endl;
    file <<"const double scan_time = "<<   scan_time <<";//second. For SMS, 2.6 second comes from the scan speed of 44cm^-1/s (Kohler's paper). Our range here is approximately between 434733-438215GHz (or 14505-14621 wavenumbers) so delta_t = 116/44 = 2.6 the time to scan the whole range. I assumed that laser stays there for 2.6/well_numbers as an approximate value." << endl;
    file <<"const double whole_burning_time = "<<   whole_burning_time <<";//scan_time/well_numbers ; // in second, for burning experiment it was 19min = 19 *60 seconds in one of our experiments." << endl;
    file <<"\n";
    file <<"//scan parameters:" << endl;
    file <<"const int number_of_scan = "<<   number_of_scan <<";" << endl;
    file <<"const double minimum_frequency = burn_frequency;//434733 * energy_conversion_factor; // MHz or 14505.5 cm^-1 or 689.4nm. This the minimum frequency that I found it postburn_spectrum.txt in a typlical hole burning calculation for flat landscape with 3000 systems." << endl;
    file <<"const double maximum_frequency = burn_frequency;//438215 * energy_conversion_factor; // MHz or 14621.7 cm^-1 or 683.9nm. This the maximum frequency that I found it postburn_spectrum.txt in a typlical hole burning calculation for flat landscape with 3000 systems." << endl;
    file <<"const double frequency_step = "<<   frequency_step <<";// MHz. I played with this value such that all wells are chosen once only." << endl;
    file <<"const int number_of_frequency_steps = ((maximum_frequency-minimum_frequency)/frequency_step)+1;" << endl;
    file <<"const double each_scan_time = scan_time * 1e6;// micro second. It means that 2.6 seconds takes to scan the whole band from minimum_frequency to maximum_frequency (see the note for whole_burning_time, or my notes in the second log book on April 9th, 2013)." << endl;
    file <<"const double time_interval_between_each_frequency_step = frequency_step * each_scan_time / (maximum_frequency - minimum_frequency);" << endl;
    file <<"const double scaling_factor = "<<   scaling_factor <<";// this factor scales the energy landscape (barrier heights and bottom of wells) in the ground satate which in turn affects the barrier heights in the excited state." << endl;
    file <<"const double scaling_parameters = "<<   scaling_parameters <<";// this is the same as scaling_factor but I want to use this parameter to change mu, sigma, mean, and deviation at the same time. The former factor is reserved for scaling an energy landscape when it is imported in the program." << endl;
    file <<"\n";
    file <<"const double mu = "<<   mu <<";// mu should be around 1000 wavenumbers (look at my log book volume II). This corresponds to the mean value of lambda close to 22 in the ground state." << endl;
    file <<"const double sigma = "<<   sigma <<";//mu/4.5;// I changed this parameter to find lambda distribution in the ground state close to 2 based on the last paper. When denominator is less than 5, random generator makes negative values for potentiel which is consequently become zero in lambda. This makes non gaussian distribution for lambda." << endl;
    file <<"const double stretch = "<<   stretch <<";//The stretch should not be exceeded 1. This resuls in the mean value of lambda close to 11 (10.94) and sigma almost 1.1 in the excited state. Changing this value affects the curvature if the energy difference between the bottom of two wells in the ground state at the bottom of parabola (in the ground state) should be a constant (~10 wavenumbers). see Curvature Calculation.nb" << endl;
    file <<"const int ensemble_number = "<<   ensemble_number <<";//" << endl;
    file <<"const int max_ensemble_for_writting = "<< max_ensemble_for_writing <<";// Some files were designed to check the validit" << endl;
    file <<"const double mean = "<<   mean <<";// the mean value for the bottom of wells" << endl;
    file <<"const double deviation = "<<   deviation <<";// the standard deviation value for the bottom of wells" << endl;
    file <<"const int generate_energy_landscape = "<<   generate_energy_landscape <<";//1 means that energy landscape is generated and 0 means that energy landscape (and Translation in case that bottom of wells are independent in two states) is/are imported." << endl;
    file <<"const double room_temperature = "<<   room_temperature <<";//Kelvin.// it could be the same as initial_temperature." << endl;
    file <<"\n";
    file <<"//parameters related to the cooling for finding the distribution before burning." << endl;
    file <<"const double extra_temperature = "<<   extra_temperature <<"; // I would like to add extra temperature to measure HGK at different temperatures while not changing the temperature step. It should be noted that the distribution at higher temperatures is close to each other." << endl;
    file <<"const double initial_temperature = room_temperature; //Kelvin." << endl;
    file <<"const double final_temperature = 5 + extra_temperature; // Kelvin." << endl;
    file <<"const int temperature_step = "<<   temperature_step <<"; // the step of reducing of temperature. NB. Changing of this parameter should not be done unless cooling time changes properly." << endl;
    file <<"const double cooling_time = "<<   cooling_time <<";//2*60*1e9 * time_conversion_factor;// In case of burning at different temperatures, the whole cooling time cannot be the same. Thus I assumed that the system cools down by 5 degrees every 2min.//  whole_cooling_time /(int ((initial_temperature - final_temperature)/temperature_step) + 1);// nano second. Taken time to cool down the sample in each step (temperature_step). It should be in such a way that the total amount of time to cool the sample down from initial T to burn T be almost an hour." << endl;
    file <<"const double whole_cooling_time = (int ((initial_temperature - final_temperature)/temperature_step) + 1) * cooling_time;//2*60*60*1e9 * time_conversion_factor;// 1 hour in nano second is 60*60*1e9 and in micro second is 60*60*1e6.  30*24*60*60*1e9 * time_conversion_factor;//" << endl;
    file <<"//const double cooling_time_step = cooling_time / 1.0; // in this way instead of using an aribitrary time step (1e-3) I am using such a time step which results in having fixed row number (e.g. 5001)." << endl;
    file <<"//NOTICE: There is a problem with the choice of time interval and step due to the round-off. For instance, for burn_time = 3.5 and burn_time_step = 0.1, A0_excited for last calculation is 0 because the loop time does not continue till the end, meaning simulation_time (in this example, burn_time). It is better to choose the time step as a multiplier of the total time." << endl;
    file <<"const double burn_time = "<<   burn_time <<"; // micro second. As rates are in GHz, the simulation time should be in ns and the maximum of it should be as equal as the life time of excited state." << endl;
    file <<"//const double burn_time_step = burn_time/1.0;" << endl;
    file <<"const int number_of_burn = "<<   number_of_burn <<";//int(whole_burning_time/time_interval_between_two_consecutive_burns);//20000;//corresponding to a 37% hole. //multiplication of this number by burn_time gives the total desired burn time, e.g. 1e6 * 3.5 ns = 3.5 ms." << endl;
    file <<"const int number_of_HGK_measurement = "<<   number_of_HGK_measurement <<";//0.5*number_of_burn;//NB: This number should be chosen in shuch a way that number_of_burn is divisible by it. The number of measurement to create the HGK curve" << endl;
    file <<"const int HGK_measurement_frequency = number_of_burn/number_of_HGK_measurement;//NB:There is a round off error when I use int(whole_burning_time/time_interval_between_two_consecutive_burns) which I did not spend time to catch it. since measuring the hole depth after each act of burn consumes a lot of memory, this parameter is defined. It means, how frequently the hole depth should be measured (e.g. every 100 acts of burn)" << endl;
    file <<"const int analytical_burn_calculation = "<<   analytical_burn_calculation <<";// 0 means calculating through multi-steps burning and 1 means calculating the result of post burn through analytical model." << endl;
    file <<"const int burn_in_one_step = "<<   burn_in_one_step <<";// in case that I don't need HGK, I can use the analytical function to calculate the result of burning in one step." << endl;
    file <<"const int recovery_between_two_acts_of_burn = "<<   recovery_between_two_acts_of_burn <<";//recovery within two consecutive absorption is activated if recovery_between_two_acts_of_burn = 1 otherwise it is not considered." << endl;
    file <<"const int recovery_correction_during_frequency_scan = "<<   recovery_correction_during_frequency_scan <<";// 1 means considering the recovery for the time that system is not in resonance with the laser during frequency scan." << endl;
    file <<"\n";
    file <<"const double whole_recovery_time = "<<  whole_recovery_time <<";//2*167*60*1e9 * time_conversion_factor;//1000*1e9 * time_conversion_factor;//167*60*1e9 * time_conversion_factor;//13*365*24*60*60*1e9 * time_conversion_factor;//time_interval_between_two_consecutive_burns*1e6;// an hour in nano second is 1 * 60*60*1e9." << endl;
    file <<"const int number_of_recovery_measurement = "<<   number_of_recovery_measurement <<";// it means that during the whole_recovery_time, we will measure the spectum number_of_recovery_measurement (e.g. 60) times." << endl;
    file <<"const double recovery_time = whole_recovery_time / number_of_recovery_measurement;// it means that calculation for each step will take recovery_time." << endl;
    file <<"//const double recovery_time_step = recovery_time / 1.0;" << endl;
    file <<"\n";
    file <<"const double number_of_longtime_measurement = "<<   number_of_longtime_measurement <<";" << endl;
    file <<"const double whole_long_time = "<< whole_long_time <<";// this time (13*365*24*60*60*1e9 = 13 years) is conidered as a longtime to calculate stationary state and then the recovery depth is compared with this distribution." << endl;
    file <<"const double long_time = whole_long_time / number_of_longtime_measurement;" << endl;
    file <<"//const double long_time_step = long_time / 1.0;" << endl;
    file <<"\n";
    file <<"const double burn_temperature = 5 + extra_temperature; //Kelvin." << endl;
    file <<"const int starting_from_equilibrium = "<<  starting_from_equilibrium <<";//It could be either 0 (non-equilibrium) or 1 (equilibrium). This parameter determins whether the pre burn distribution is the result of cooling cycle (starting_from_equilibrium = 0) or is Boltzman distribution at burn temperature (or equilibrium, starting_from_equilibrium = 1)." << endl;
    file <<"const int starting_from_uniform_distribution = "<<   starting_from_uniform_distribution <<";//It could be either 0 (non uniform initial distribution, it could be either non-equilibrium or the results of cooling) or 1 (uniform distribution, i.e. 1/22)." << endl;
    file <<"const int writing_details_of_burning_in_reson_well = "<<   writing_details_of_burning_in_reson_well <<"; // Writing of prob_in_reson_well_after_each_act_burn needs a great deal of memory. 0 means avoid doing that." << endl;
    file <<"\n";
    file <<"const double highest_experimental_cycling_temperature = "<<   highest_experimental_cycling_temperature <<";// This is the highest temperature during thermocycling experiment." << endl;
    file <<"const double highest_cycling_temperature = highest_experimental_cycling_temperature-burn_temperature;// Since the first temperature in thermocycling loop should be greater than burn temperature which results in having final temperature higher than what we have in practice (experimental_cycling_temperature)." << endl;
    file <<"//const int cycle_temperature_step = "<<   "NOT DEFINED" <<";" << endl;
    file <<"\n";
    file <<"const double whole_thermocycling_time = "<<   whole_thermocycling_time <<";// The measurment time for warming the sample up is about 250 minutes (this is half of the whole thermocycling experiment that we usually do in the lab). This estimation is not completely correct since it does not take the same time for warming up and cooling down. In the first approximation, I assumed that temperature profile is linear and the same for both warming and cooling. Also, the sample does not stay at each cycling maximum temperature for the same time, however, I belive that it is not a big issue since the rate at low temperature (which sample reaches to them faster) is low enough that there is no significant difference between 5 or 20 minutes." << endl;
    file <<"const int number_of_thermocycling_measurement = "<<   number_of_thermocycling_measurement <<";" << endl;
    file <<"//const double cycle_temperature_step = cycle_temperature-burn_temperature / number_of_thermocycling_measurement;" << endl;
    file <<"const double thermocycling_time = whole_thermocycling_time / number_of_thermocycling_measurement;// it means that calculation for each step will take thermocycling_time." << endl;
    file <<"//const double thermocycling_time_step = thermocycling_time / 1.0;" << endl;
    file <<"\n";
    file <<"\n";
    file <<"const double recovery_temperature = burn_temperature;//Kelvin." << endl;
    file <<"\n";
    file <<"\n";
    file <<"const int number_of_row_burning = "<<   number_of_row_burning <<";//burn_time/burn_time_step;// by omitting the loop in the SolveRateMatrix, there is no need to have this line. the number of row in the rate_file (or rates_of_ensemble.txt) is related to the integer answer of dividing burn_time by time_step" << endl;
    file <<"const int number_of_row_cooling = "<<   number_of_row_cooling <<";//cooling_time/cooling_time_step;// the number of row in the rates_file.txt is related to the integer answer of dividing cooling_time by cooling_time_step (this is the step numbers which takes the sample to be cooled down)" << endl;
    file <<"const int number_of_row_recovery = "<<  number_of_row_recovery <<";//recovery_time/recovery_time_step;" << endl;
    file <<"const int number_of_row_longtime = "<<   number_of_row_longtime <<";//long_time/long_time_step;" << endl;
    file <<"const int number_of_row_thermocycling = "<<   number_of_row_thermocycling <<";//thermocycling_time/thermocycling_time_step;" << endl;
    file <<"\n";
    file <<"\n";
    file <<"const int number_of_column = well_numbers * ensemble_number + 1;// the first column is for simulation time so the number of column should be increased by 1" << endl;
    file <<"\n";
    file <<"//well parameters" << endl;
    file <<"\n";
    file <<"\n";
    file <<"const double step = (2*interval)/(point_numbers-1); // in general, (max-min)/(point_numbers-1)" << endl;
    file <<"\n";
    file <<"const double step_well = step * (point_numbers/well_numbers);// as point numbers are double of well numbers here, the step_well is double of step" << endl;
    file <<"\n";
    file <<"//const double roughness = "<<   "NOT DEFINED" <<";// the random numbers can be multiplied by this number as it has been done for sinusoide." << endl;
    file <<"const int well_numbers = well_numbers-1; //There are as many state (or concentration elements [A], [B], ...) as (random numbers - 1). It should be mentioned as the first well is always the first point on the most left part of graph then having concentration beyond of this point is not meaningful. Thus, one well means one concentration instead of two." << endl;
    file <<"\n";
    file <<"const double d_distance_of_two_wells = (step_well) * 1e-9;// Nanometer. The distance between two adjacent wells is 2 times greater than step (or is equal to step_well) when point_numbers = 2 * well_numbers." << endl;
    file <<"\n";
    file <<"\n";
    file <<"const double lambda_coefficient = d_distance_of_two_wells * pi * sqrt(8 * mass * 1e6 / h_planck );// lambda_coefficient is the constant part of lambda, which is calculated in the createRateMatrix.h" << endl;
    file <<"\n";
    file <<"const int all_frequency = well_numbers * ensemble_number;// all possible frequencies which come out of the whole ensemble." << endl;
    file <<"const int total_bins = sqrt (all_frequency); // a fixed number. sqrt (all_frequency) is very big for number of bins, I changed it to ensemble number. Square-root choice is used in Excel and many others" << endl;
    file <<"const int averaged_probability_SMS_bin_number = "<<    averaged_probability_SMS_bin_number <<";" << endl;
    file <<"\n";
    file <<"#endif // CONSTANTS_H" << endl;
}




void WriteCutOffEnergy (double cutoff_energy_difference){

    string filename;
    filename = "cutoff_energy=" + to_string(cutoff_energy_difference/energy_conversion_factor) + " GHz.txt";
	fstream file;
	file.open(filename, ios::out | ios::app);

	if (!file) {
		cout << "Error opening the file" << filename;
	}
    file << "Cut-off energy difference is " << cutoff_energy_difference / energy_conversion_factor << " GHz." << endl;

    //file.close();//NOTICE: Since I am useing append, I am not sure about closing this file.
}


void WriteBarrierHeightandParabolaforRate (gsl_vector * parabola_for_rate, string filename1, gsl_vector * barrier_height, string filename2, int ensemble, int recovery_counter){

    fstream file1;
    file1.open(filename1, ios::out | ios::app);
	if (!file1) {
		cout << "Error opening the file" << filename1;
	}
    file1 << "ensemble\t" << ensemble << "\tRecoveryNumber\t" << recovery_counter << endl;

    for (unsigned int i=0; i < parabola_for_rate->size ; i++)
        file1 << gsl_vector_get (parabola_for_rate,i) << endl;




	fstream file2;
	file2.open("BarrierHeight.txt", ios::out | ios::app);
	if (!file2) {
		cout << "Error opening the file BarrierHeight.txt";
	}
    file2 << "ensemble\t" << ensemble << "\tRecoveryNumber\t" << recovery_counter << endl;

    for (unsigned int i=0; i < barrier_height->size ; i++)
        file2 << gsl_vector_get (barrier_height,i) << endl;

    //file1.close();//NOTICE: Since I am useing append, I am not sure about closing this file.
    //file2.close();//NOTICE: Since I am useing append, I am not sure about closing this file.

}
