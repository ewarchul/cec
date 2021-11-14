#include "affine_trans.h"

void shiftfunc(double *x, double *xshift, int nx, double *Os) {
  int i;
  for (i = 0; i < nx; i++) {
    xshift[i] = x[i] - Os[i];
  }
}

// TODO fix problem with dimension indices

double *shift_modern(size_t dim, double input[dim], int problem_num,
                     cec_state_t *state) {
  //double *output = calloc(10, sizeof(double));
  double *output = malloc(dim * sizeof(double));
  for (size_t i = 0; i < dim; ++i) {
    output[i] =
        input[i] - state->data_.shifts_[0].data_[problem_num - 1].data_[i];
  }
  return output;
}

void rotatefunc(double *x, double *xrot, int nx, double *Mr) {
  int i, j;
  for (i = 0; i < nx; i++) {
    xrot[i] = 0;
    for (j = 0; j < nx; j++) {
      xrot[i] = xrot[i] + x[j] * Mr[i * nx + j];
    }
  }
}

/*double *rotate_modern(size_t dim, double input[dim], int problem_num,*/
/*cec_state_t *state) {*/
/*double *output = calloc(dim, sizeof dim);*/
/*for (size_t i = 0; i < dim; ++i) {*/
/*for (size_t j = 0; j < dim; ++j) {*/
/*output[i] =*/
/*output[i] +*/
/*input[j] * state->rotate_values_.data_[dim][problem_num][i * dim + j];*/
/*}*/
/*}*/
/*return output;*/
/*}*/

/*double *apply_transformation_rate(size_t dim, double input[dim],*/
/*double t_rate) {*/
/*double *output = calloc(dim, sizeof dim);*/
/*for (size_t i = 0; i < dim; ++i) {*/
/*output[i] *= input[i] * t_rate;*/
/*}*/
/*return output;*/
/*}*/

/*double *shift_rotate_modern(size_t dim, double input[dim], int problem_num,*/
/*cec_state_t *state, cec_affine_transforms_t info) {*/
/*double *output = calloc(dim, sizeof dim);*/
/*if (info.shift_ == true) {*/
/*if (info.rotate_ == true) {*/
/*output = shift_modern(dim, input, problem_num, state);*/
/*output = apply_transformation_rate(dim, output, info.transform_rate_);*/
/*output = rotate_modern(dim, output, problem_num, state);*/
/*} else {*/
/*output = shift_modern(dim, input, problem_num, state);*/
/*output = apply_transformation_rate(dim, output, info.transform_rate_);*/
/*}*/
/*} else {*/
/*if (info.rotate_ == 1) {*/
/*output = apply_transformation_rate(dim, output, info.transform_rate_);*/
/*output = rotate_modern(dim, output, problem_num, state);*/
/*} else {*/
/*output = apply_transformation_rate(dim, output, info.transform_rate_);*/
/*}*/
/*}*/

/*return output;*/
/*}*/

void sr_func(double *x, double *sr_x, int nx, double *Os, double *Mr,
             double sh_rate, int s_flag, int r_flag, double *y) {
  int i;
  if (s_flag == 1) {
    if (r_flag == 1) {
      shiftfunc(x, y, nx, Os);
      for (i = 0; i < nx; i++) {
        y[i] = y[i] * sh_rate;
      }
      rotatefunc(y, sr_x, nx, Mr);
    } else {
      shiftfunc(x, sr_x, nx, Os);
      for (i = 0; i < nx; i++) {
        sr_x[i] = sr_x[i] * sh_rate;
      }
    }
  } else {

    if (r_flag == 1) {
      for (i = 0; i < nx; i++) {
        y[i] = x[i] * sh_rate;
      }
      rotatefunc(y, sr_x, nx, Mr);
    } else
      for (i = 0; i < nx; i++) {
        sr_x[i] = x[i] * sh_rate;
      }
  }
}

void cf_cal(double *x, double *f, int nx, double *Os, double *delta,
            double *bias, double *fit, int cf_num) {
  int i, j;
  double *w;
  double w_max = 0, w_sum = 0;
  w = (double *)malloc(cf_num * sizeof(double));
  for (i = 0; i < cf_num; i++) {
    fit[i] += bias[i];
    w[i] = 0;
    for (j = 0; j < nx; j++) {
      w[i] += pow(x[j] - Os[i * nx + j], 2.0);
    }
    if (w[i] != 0)
      w[i] = pow(1.0 / w[i], 0.5) * exp(-w[i] / 2.0 / nx / pow(delta[i], 2.0));
    else
      w[i] = INF;
    if (w[i] > w_max)
      w_max = w[i];
  }
  for (i = 0; i < cf_num; i++) {
    w_sum = w_sum + w[i];
  }
  if (w_max == 0) {
    for (i = 0; i < cf_num; i++)
      w[i] = 1;
    w_sum = cf_num;
  }
  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] = f[0] + w[i] / w_sum * fit[i];
  }
  free(w);
}
