/*
 * Complex functions
 *
 * The signature of each complex function is given below:
 *
 * `void cec20XY_hf0Z(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag)`
 *
 * where `x` and `f` stands for, respectively, input and output vectors, `nx`
 * is problem dimension and `Os`, `Mr`, are vector with problem specific data.
 * Last parameter, `r_flag`, is binary flags which decides if function will be
 * rotated.
 *
 * More information about complex functions contain technical documentations of
 * CEC benchmarks.
 *
 */

#ifndef COMPLEX_FUNCS_H
#define COMPLEX_FUNCS_H

#include "basic_funcs.h"
#include "hybrid_funcs.h"
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E 2.7182818284590452353602874713526625
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029
#endif

void cec2014_cf01(double *, double *, int, double *, double *, int);
void cec2014_cf02(double *, double *, int, double *, double *, int);
void cec2014_cf03(double *, double *, int, double *, double *, int);
void cec2014_cf04(double *, double *, int, double *, double *, int);
void cec2014_cf05(double *, double *, int, double *, double *, int);
void cec2014_cf06(double *, double *, int, double *, double *, int);
void cec2014_cf07(double *, double *, int, double *, double *, int *, int);
void cec2014_cf08(double *, double *, int, double *, double *, int *, int);

void cec2015_cf01(double *, double *, int, double *, double *, double *, int);
void cec2015_cf02(double *, double *, int, double *, double *, int *, double *,
                  int);
void cec2015_cf03(double *, double *, int, double *, double *, double *, int);
void cec2015_cf04(double *, double *, int, double *, double *, double *, int);
void cec2015_cf05(double *, double *, int, double *, double *, int *, double *,
                  int);
void cec2015_cf06(double *, double *, int, double *, double *, double *, int);
void cec2015_cf07(double *, double *, int, double *, double *, double *, int);

void cec2017_cf01(double *, double *, int, double *, double *, int);
void cec2017_cf02(double *, double *, int, double *, double *, int);
void cec2017_cf03(double *, double *, int, double *, double *, int);
void cec2017_cf04(double *, double *, int, double *, double *, int);
void cec2017_cf05(double *, double *, int, double *, double *, int);
void cec2017_cf06(double *, double *, int, double *, double *, int);
void cec2017_cf07(double *, double *, int, double *, double *, int);
void cec2017_cf08(double *, double *, int, double *, double *, int);
void cec2017_cf09(double *, double *, int, double *, double *, int *, int);
void cec2017_cf10(double *, double *, int, double *, double *, int *, int);

void cec2021_cf01(double *, double *, int, double *, double *, int);
void cec2021_cf02(double *, double *, int, double *, double *, int);
void cec2021_cf03(double *, double *, int, double *, double *, int);
void cec2021_cf01_s(double *, double *, int, double *, double *, int);
void cec2021_cf02_s(double *, double *, int, double *, double *, int);
void cec2021_cf03_s(double *, double *, int, double *, double *, int);
#endif
