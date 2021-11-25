/*
 * Affine transformations
 *
 * CEC's functions are usually transformed by written below affine
 * transformations:
 *
 * - rotation (`rotatefunc`)
 * - shift (`shiftfunc`)
 * - shift & rotate (`sr_func`)
 *
 * `asyfunc` and `osyfunc` are specific transformations for CEC2013. Check
 * technical raport for more information.
 *
 * `cf_cal` is helper function for complex functions.
 */

#ifndef AFFINE_TRANS_H
#define AFFINE_TRANS_H

#include "types.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E 2.7182818284590452353602874713526625
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029
#endif

double *rotate_modern(double *input, int problem_num, cec_state_t *state);
double *apply_transformation_rate(size_t dim, double input[dim], double t_rate);
double *shift_rotate_modern(double *input, int problem_num,
                            cec_state_t *state, cec_affine_transforms_t info);
double *shift_modern(double *input, int problem_num, cec_state_t *state);



void shiftfunc(double *, double *, int, double *);
void rotatefunc(double *, double *, int, double *);
void sr_func(double *, double *, int, double *, double *, double, int, int,
             double *);
void asyfunc(double *, double *x, int, double);
void oszfunc(double *, double *, int);
void cf_cal(double *, double *, int, double *, double *, double *, double *,
            int);
#endif
