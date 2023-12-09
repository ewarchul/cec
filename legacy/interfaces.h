/*
 * Interfaces for the particular CEC versions.
 * Almost each interface has signature as follows:
 *
 * `cec20XY_interface(char *, double *, double *, int, int, int)`
 *
 * where `datapath` is path to CEC's data, `x` is the input argument, 
 * `f` is the array for the output value, `nx` is the row number of given input vector
 * as a input, `mx` is its column number and `func_num` is index of
 * optimization function in benchmark.
 * The very one execption is CEC2021 which takes one additional parameter, i.e 
 * `suite` of benchmark. 
 * 
 * For more information check technical documentations stored in `doc`
 * directory. 
 */

#ifndef INTRERFACE_H
#define INTERFACES_H

#include "hybrid_funcs.h"
#include "complex_funcs.h"
#include "basic_funcs.h"
#include "utils.h"
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

void cec2013_interface(char *, double *, double *, int, int, int);
void cec2014_interface(char *, double *, double *, int, int, int);
void cec2015_interface(char *, double *, double *, int, int, int);
void cec2017_interface(char *, double *, double *, int, int, int);
void cec2019_interface(char *, double *, double *, int, int, int);
void cec2021_interface(char *, double *, double *, int, int, int, char *);

#endif 
