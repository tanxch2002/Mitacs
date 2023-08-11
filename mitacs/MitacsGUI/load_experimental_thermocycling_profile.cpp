#include "pch.h"
#include "Global_variables.h"
#include "load_experimental_thermocycling_profile.h"
// This function loads time and temperature according to the experimental profile. 
//Time is converted to proper unit (e.g. micro second) and temperature is in kelvin.
void LoadExperimentalThermocyclingProfile (gsl_matrix * experimental_thermocycling_time_temperature)
{

double time_factor = 60*1e9 * time_conversion_factor;// it converts time from minute to microsecond.
double thermocycling_time [17]{15, 24, 18, 25, 22, 26, 24, 31, 23, 41, 25, 42, 26, 36, 31, 38, 64};
// measured time in experiment is in minutes.
double thermocycling_temperature[17]{7, 9, 11, 13, 15, 17, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49};
// maximal cycling temperature.
FILE* f;
fopen_s(&f, "thermocycling_profile.dat", "r");
gsl_matrix_fscanf(f, experimental_thermocycling_time_temperature); // file to matrix
fclose(f);

for (int i=0; i < number_of_thermocycling_measurement; ++i)
{
    thermocycling_time[i] = time_factor * thermocycling_time[i];
    gsl_matrix_set (experimental_thermocycling_time_temperature, i, 0, thermocycling_time[i]);
    gsl_matrix_set (experimental_thermocycling_time_temperature, i, 1, thermocycling_temperature[i]);
}

}
