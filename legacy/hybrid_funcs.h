/* 
 * Hybrid functions
 *
 * The signature of each hybrid function is given below:
 *
 * `void cec20XY_hf0Z(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag)`
 *
 * where `x` and `f` stands for, respectively, input and output vectors, `nx`
 * is problem dimension and `Os`, `Mr`, `S` are vector with problem specific data. 
 * Last two parameteres are binary flags which decides if function will be
 * shifted (`s_flag`) or rotated (`r_flag`).
 *
 * More information about hybrid functions contain technical documentations of
 * CEC benchmarks.
 *
 */

#ifndef HYBRID_FUNCS_H
#define HYBRID_FUNCS_H

#include "basic_funcs.h"
#include <stdlib.h>
#include <stddef.h>
#include <math.h>


void cec2014_hf01(double *, double *, int, double *, double *, int *, int, int);
void cec2014_hf02(double *, double *, int, double *, double *, int *, int, int);
void cec2014_hf03(double *, double *, int, double *, double *, int *, int, int);
void cec2014_hf04(double *, double *, int, double *, double *, int *, int, int);
void cec2014_hf05(double *, double *, int, double *, double *, int *, int, int);
void cec2014_hf06(double *, double *, int, double *, double *, int *, int, int);

void cec2015_hf01(double *, double *, int, double *, double *, int *, int, int);
void cec2015_hf02(double *, double *, int, double *, double *, int *, int, int);
void cec2015_hf03(double *, double *, int, double *, double *, int *, int, int);

void cec2017_hf01(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf02(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf03(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf04(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf05(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf06(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf07(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf08(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf09(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf10(double *, double *, int, double *, double *, int *, int, int);

void cec2021_hf01(double *, double *, int, double *, double *, int *, int, int);
void cec2021_hf02(double *, double *, int, double *, double *, int *, int, int);
void cec2021_hf03(double *, double *, int, double *, double *, int *, int, int);

#endif
