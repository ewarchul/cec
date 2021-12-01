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

inline int cec_dimension_idx(size_t dim) {
  return dim == 10 ? 0 : dim == 30 ? 1 : dim == 50 ? 2 : dim == 100 ? 3 : -1;
}

inline double *shift_modern(double *input, int fn, cec_state_t *state) {
  int dim = state->dimension_;
  double *output = calloc(dim, sizeof(double));
  for (int i = 0; i < dim; ++i) {
    output[i] = input[i] - state->data_.shifts_[cec_dimension_idx(dim)]
                               .data_[fn - 1]
                               .data_[i];
  }
  return output;
}

inline double *shuffle_modern(size_t dim, int fn, double *input,
                              cec_state_t *state) {
  double *output = calloc(dim, sizeof(double));
  for (size_t i = 0; i < dim; ++i) {
    int index =
        state->data_.shuffles_[cec_dimension_idx(dim)].data_[fn - 1].data_[i] -
        1;
    output[i] = input[index];
  }
  return output;
}

inline double *rotate_modern(double *input, int fn, cec_state_t *state) {
  size_t dim = state->dimension_;
  double *output = calloc(dim, sizeof(double));
  for (size_t i = 0; i < dim; ++i) {
    output[i] = 0;
    for (size_t j = 0; j < dim; ++j) {
      output[i] =
          output[i] + input[j] * state->data_.rotates_[cec_dimension_idx(dim)]
                                     .data_[fn - 1]
                                     .data_[i * dim + j];
    }
  }
  return output;
}

inline double *apply_transformation_rate(size_t dim, double input[dim],
                                         double t_rate) {
  double *output = calloc(dim, sizeof(double));
  for (size_t i = 0; i < dim; ++i) {
    output[i] = input[i] * t_rate;
  }
  return output;
}
inline double *shift_rotate_modern(double *input, int problem_num,
                                   cec_state_t *state,
                                   cec_affine_transforms_t info) {
  int dim = state->dimension_;
  double *output = calloc(dim, sizeof(double));
  if (info.shift_) {
    if (info.rotate_) {
      output = shift_modern(input, problem_num, state);
      output = apply_transformation_rate(dim, output, info.transform_rate_);
      output = rotate_modern(output, problem_num, state);
    } else {
      output = shift_modern(input, problem_num, state);
      output = apply_transformation_rate(dim, output, info.transform_rate_);
    }
  } else {
    if (info.rotate_) {
      output = apply_transformation_rate(dim, output, info.transform_rate_);
      output = rotate_modern(output, problem_num, state);
    } else {
      output = apply_transformation_rate(dim, output, info.transform_rate_);
    }
  }

  return output;
}

void shiftfunc(double *, double *, int, double *);
void rotatefunc(double *, double *, int, double *);
void sr_func(double *, double *, int, double *, double *, double, int, int,
             double *);
void asyfunc(double *, double *x, int, double);
void oszfunc(double *, double *, int);
void cf_cal(double *, double *, int, double *, double *, double *, double *,
            int);
#endif
