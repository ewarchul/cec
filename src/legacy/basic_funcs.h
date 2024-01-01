/*
 * Basic functions (unimodal or multimodal)
 *
 * Each function takes at least two parameters: `x` given input and `f` output
 * vector. Additional parameters are the same as for hybrid and complex
 * function.
 *
 * Basic functions are building blocks for hybrid and complex function. 
 * More information about basic functions contain technical documentations of
 * CEC benchmarks.
 *
 */


#ifndef BASIC_FUNCS_H
#define BASIC_FUNCS_H

#include "affine_trans.h"
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E 2.7182818284590452353602874713526625
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029
#endif

void Lennard_Jones(double *, int, double *); 
void Hilbert(double *, int, double *);       
void Chebyshev(double *, int, double *);     
void sphere_func(double *, double *, int, double *, double *, int, int);
void ellips_func(double *, double *, int, double *, double *, int, int);
void bent_cigar_func(double *, double *, int, double *, double *, int, int);


void cec2017_bent_cigar_func(double *, double *, int, double *, double *, int, int);

void discus_func(double *, double *, int, double *, double *, int, int);
void dif_powers_func(double *, double *, int, double *, double *, int, int);
void rosenbrock_func(double *, double *, int, double *, double *, int, int);
void schaffer_F7_func(double *, double *, int, double *, double *, int, int, double *);
void ackley_func(double *, double *, int, double *, double *, int, int);
void rastrigin_func(double *, double *, int, double *, double *, int, int);
void weierstrass_func(double *, double *, int, double *, double *, int, int);
void griewank_func(double *, double *, int, double *, double *, int, int);
void schwefel_func(double *, double *, int, double *, double *, int, int);
void katsuura_func(double *, double *, int, double *, double *, int, int);
void bi_rastrigin_func(double *, double *, int, double *, double *, int, int);
void grie_rosen_func(double *, double *, int, double *, double *, int, int);
void escaffer6_func(double *, double *, int, double *, double *, int, int);
void step_rastrigin_func(double *, double *, int, double *, double *, int, int);
void happycat_func(double *, double *, int, double *, double *, int, int);
void hgbat_func(double *, double *, int, double *, double *, int, int);
void sum_diff_pow_func(double *, double *, int, double *, double *, int, int);
void zakharov_func(double *, double *, int, double *, double *, int, int);
void levy_func(double *, double *, int, double *, double *, int, int);
void dixon_price_func(double *, double *, int, double *, double *, int, int);

#endif
