
#include "pch.h"
#include <cmath>
#include "Global_variables.h"
int to1D(int x, int y, int z) //x-system number (xMax=ensemble_number), y - well number, z - state number
{
	
	int a; 
	a = z * ensemble_number * well_numbers + y * ensemble_number + x; 
	return a;
}
// if it is done like this, it does the following (two molecules, four wells, 512 states)
// molecule0,well 0, state 0  ....0
// molecule0,well 0, state 1  ....8
// molecule0,well 0, state 2  ....16
//...
// molecule0,well 0, state 511  ....4088
// //
// molecule0,well y=1, state 0  ....2
// molecule0,well y=1, state 1  ....10
// molecule0,well y=1, state 2  ....18
///...
// molecule0,well y=1, state 511  ....4090
// //
// molecule0,well y=2, state 0  ....4
// molecule0,well y=2, state 1  ....12
// molecule0,well y=2, state 2  ....20
///...
// molecule0,well y=2, state 511  ....4092
// 
// molecule0,well y=3, state 0  ....6
// molecule0,well y=3, state 1  ....14
// molecule0,well y=3, state 2  ....22
///...
// molecule0,well y=3, state 511  ....4094
///
// second molecule starts here
// molecule1,well y=0, state 0  ....1
// molecule1,well y=0, state 1  ....9
// molecule1,well y=0, state 2  ....17
///...
// molecule1,well y=0, state 511  ....4089