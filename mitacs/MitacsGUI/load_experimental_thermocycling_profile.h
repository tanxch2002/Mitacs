#ifndef LOAD_EXPERIMENTAL_THERMOCYCLING_PROFILE_H
#define LOAD_EXPERIMENTAL_THERMOCYCLING_PROFILE_H
#include <gsl/gsl_matrix.h> // it is necessary for the rate matrix definition
#include "Global_variables.h"

void LoadExperimentalThermocyclingProfile (gsl_matrix * experimental_thermocycling_time_temperature);

#endif // LOAD_EXPERIMENTAL_THERMOCYCLING_PROFILE_H
