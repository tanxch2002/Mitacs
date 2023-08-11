#include "pch.h"
#include "check_situation.h"
//using namespace std;
void CheckSituation (double temperature, CString &situation)
{

    if (temperature == -1000) situation = "burn";
    else
        if (temperature == -100) situation = "recovery";
    else
        if (temperature == -10) situation = "thermocycling";
    else situation = "cooling";

}
