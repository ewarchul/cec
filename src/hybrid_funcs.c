#include "hybrid_funcs.h"
#include "affine_trans.h"
#include "basic_funcs.h"
#include "types.h"
#include <string.h>

void cec2014_hf01(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 3;
  double fit[3];
  int G[3], G_nx[3];
  double Gp[3] = {0.3, 0.3, 0.4};
  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;
  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2014_hf02(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 3;
  double fit[3];
  int G[3], G_nx[3];
  double Gp[3] = {0.3, 0.3, 0.4};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  bent_cigar_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2014_hf03(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 4;
  double fit[4];
  int G[4], G_nx[4];
  double Gp[4] = {0.2, 0.2, 0.3, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  griewank_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  weierstrass_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2014_hf04(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) /* Hybrid Function 4 */
{
  int i, tmp, cf_num = 4;
  double fit[4];
  int G[4], G_nx[4];
  double Gp[4] = {0.2, 0.2, 0.3, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  discus_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  grie_rosen_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}
void cec2014_hf05(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) /* Hybrid Function 5 */
{
  int i, tmp, cf_num = 5;
  double fit[5];
  int G[5], G_nx[5];
  double Gp[5] = {0.1, 0.2, 0.2, 0.2, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2014_hf06(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) /* Hybrid Function 6 */
{
  int i, tmp, cf_num = 5;
  double fit[5];
  int G[5], G_nx[5];
  double Gp[5] = {0.1, 0.2, 0.2, 0.2, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  katsuura_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  happycat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  grie_rosen_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  ackley_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2015_hf01(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) /* Hybrid Function 1 */
{
  int i, tmp, cf_num = 3;
  double fit[3];
  int G[3], G_nx[3];
  double Gp[3] = {0.3, 0.3, 0.4};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;
  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2015_hf02(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) /* Hybrid Function 2 */
{
  int i, tmp, cf_num = 4;
  double fit[4];
  int G[4], G_nx[4];
  double Gp[4] = {0.2, 0.2, 0.3, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  griewank_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  weierstrass_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2015_hf03(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) /* Hybrid Function 3 */
{
  int i, tmp, cf_num = 5;
  double fit[5];
  int G[5], G_nx[5];
  double Gp[5] = {0.1, 0.2, 0.2, 0.2, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2017_hf01(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 3;
  double fit[3];
  int G[3], G_nx[3];
  double Gp[3] = {0.2, 0.4, 0.4};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;
  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                      */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  zakharov_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

typedef struct cec_shuffles_t cec_shuffles_t;
struct cec_shuffles_t {
  int shifts_[6];
  int partition_idx_[6];
};

cec_shuffles_t mk_shuffles(size_t dim, int fn_nums, double weights[fn_nums]) {
  int partition_idx[fn_nums];
  double tmp = 0;
  for (int i = 0; i < fn_nums - 1; ++i) {
    partition_idx[i] = ceil(weights[i] * dim);
    tmp += partition_idx[i];
  }
  partition_idx[fn_nums - 1] = dim - tmp;

  int shifts[fn_nums];
  shifts[0] = 0;
  for (int i = 1; i < fn_nums; ++i) {
    shifts[i] = shifts[i - 1] + partition_idx[i - 1];
  }

  cec_shuffles_t shuffles;
  memcpy(shuffles.shifts_, shifts, fn_nums * sizeof(int));
  memcpy(shuffles.partition_idx_, partition_idx, fn_nums * sizeof(int));

  return shuffles;
}

double cec2017_hf01_modern(size_t dim, int fn, double *input,
                           cec_state_t *state) {

  double weights[3] = {0.2, 0.4, 0.4};
  cec_shuffles_t shuffles = mk_shuffles(dim, 3, weights);

  cec_affine_transforms_t af_trans = {
      .shift_ = true, .rotate_ = true, .transform_rate_ = 1};
  double *shift_rotated = shift_rotate_modern(input, fn, state, af_trans);
  double *shuffled = shuffle_modern(dim, fn, shift_rotated, state);

  double y_0 = zakharov_func_modern(shuffles.partition_idx_[0],
                                    (shuffled + shuffles.shifts_[0]));
  double *tmp = apply_transformation_rate(10, shuffled, 2.048 / 100.0);

  double y_1 = rosenbrock_func_modern(shuffles.partition_idx_[1],
                                      (tmp + shuffles.shifts_[1]));

  double *tmp2 = apply_transformation_rate(10, shuffled, 5.12 / 100.0);
  double y_2 = rastrigin_func_modern(shuffles.partition_idx_[2],
                                     (tmp2 + shuffles.shifts_[2]));

  return y_0 + y_1 + y_2;
}

double cec2017_hf02_modern(size_t dim, int fn, double *input,
                           cec_state_t *state) {

  double weights[3] = {0.3, 0.3, 0.4};
  cec_shuffles_t shuffles = mk_shuffles(dim, 3, weights);

  cec_affine_transforms_t af_trans = {
      .shift_ = true, .rotate_ = true, .transform_rate_ = 1};
  double *shift_rotated = shift_rotate_modern(input, fn, state, af_trans);
  double *shuffled = shuffle_modern(dim, fn, shift_rotated, state);

  double y_0 = ellips_func_modern(shuffles.partition_idx_[0],
                                  (shuffled + shuffles.shifts_[0]));
  double *tmp = apply_transformation_rate(10, shuffled, 10.0);

  double y_1 = schwefel_func_modern(shuffles.partition_idx_[1],
                                    (tmp + shuffles.shifts_[1]));

  double y_2 = bent_cigar_func_modern(shuffles.partition_idx_[2],
                                      (shuffled + shuffles.shifts_[2]));

  return y_0 + y_1 + y_2;
}

void cec2017_hf02(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 3;
  double fit[3];
  int G[3], G_nx[3];
  double Gp[3] = {0.3, 0.3, 0.4};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                      */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  bent_cigar_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

double cec2017_hf03_modern(size_t dim, int fn, double *input,
                           cec_state_t *state) {

  double weights[3] = {0.3, 0.3, 0.4};
  cec_shuffles_t shuffles = mk_shuffles(dim, 3, weights);

  cec_affine_transforms_t af_trans = {
      .shift_ = true, .rotate_ = true, .transform_rate_ = 1};
  double *shift_rotated = shift_rotate_modern(input, fn, state, af_trans);
  double *shuffled = shuffle_modern(dim, fn, shift_rotated, state);

  double y_0 = bent_cigar_func_modern(shuffles.partition_idx_[0],
                                  (shuffled + shuffles.shifts_[0]));
  double *tmp = apply_transformation_rate(10, shuffled, 2.048 / 100.0);
  //TODO fix rosenbrock
  double y_1 = rosenbrock_func_modern(shuffles.partition_idx_[1],
                                    (tmp + shuffles.shifts_[1]));
  double y_2 = bi_rastrigin_func_modern(shuffles.partition_idx_[2],
                                      (shuffled + shuffles.shifts_[2]), state, fn);

  return y_0 + y_1 + y_2;
}

void cec2017_hf03(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 3;
  double fit[3];
  int G[3], G_nx[3];
  double Gp[3] = {0.3, 0.3, 0.4};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                      */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }

  i = 0;
  bent_cigar_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  bi_rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2017_hf04(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 4;
  double fit[4];
  int G[4], G_nx[4];
  double Gp[4] = {0.2, 0.2, 0.2, 0.4};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                      */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  ackley_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  schaffer_F7_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0, y);
  i = 3;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2017_hf05(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 4;
  double fit[4];
  int G[4], G_nx[4];
  double Gp[4] = {0.2, 0.2, 0.3, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                      */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;

  bent_cigar_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}
void cec2017_hf06(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 4;
  double fit[4];
  int G[4], G_nx[4];
  double Gp[4] = {0.2, 0.2, 0.3, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                      */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2017_hf07(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 5;
  double fit[5];
  int G[5], G_nx[5];
  double Gp[5] = {0.1, 0.2, 0.2, 0.2, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                      */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  katsuura_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  ackley_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  grie_rosen_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2017_hf08(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 5;
  double fit[5];
  int G[5], G_nx[5];
  double Gp[5] = {0.2, 0.2, 0.2, 0.2, 0.2};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                      */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }

  i = 0;
  ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  ackley_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  discus_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2017_hf09(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 5;
  double fit[5];
  int G[5], G_nx[5];
  double Gp[5] = {0.2, 0.2, 0.2, 0.2, 0.2};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  bent_cigar_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  grie_rosen_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  weierstrass_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2017_hf10(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) {
  int i, tmp, cf_num = 6;
  double fit[6];
  int G[6], G_nx[6];
  double Gp[6] = {0.1, 0.1, 0.2, 0.2, 0.2, 0.2};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                      */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }

  i = 0;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  katsuura_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  ackley_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 5;
  schaffer_F7_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0, y);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2021_hf01(
    double *x, double *f, int nx, double *Os, double *Mr, int *S, int s_flag,
    int r_flag) /* Hybrid Function 1  /F17 Hybrid Function 1 in cec2014*/
{
  int i, tmp, cf_num = 3;
  double fit[3];
  int G[3], G_nx[3];
  double Gp[3] = {0.3, 0.3, 0.4};

  tmp = 0;
  for (i = 1; i < cf_num; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  // G_nx[cf_num-1]=nx-tmp;
  G_nx[0] = nx - tmp;
  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  i = 1;
  rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  i = 2;
  ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2021_hf02(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) /* Hybrid Function 5 */
{
  int i, tmp, cf_num = 4;
  double fit[4];
  int G[4], G_nx[4];
  double Gp[4] = {0.2, 0.2, 0.3, 0.3};

  if (nx == 5) {

    G_nx[0] = 1;
    G_nx[1] = 1;
    G_nx[2] = 1;
    G_nx[3] = 2;

  } else {
    tmp = 0;
    for (i = 1; i < cf_num; i++) {
      G_nx[i] = ceil(Gp[i] * nx);
      tmp += G_nx[i];
    }
    G_nx[0] = nx - tmp;
  }

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2021_hf03(double *x, double *f, int nx, double *Os, double *Mr, int *S,
                  int s_flag, int r_flag) // hf5 in cec2014 F21 in cec2014
{
  int i, tmp, cf_num = 5;
  double fit[5];
  int G[5], G_nx[5];
  double Gp[5] = {0.1, 0.2, 0.2, 0.2, 0.3};

  if (nx == 5) // deal with D=6**04/01/2020
  {
    G_nx[0] = 1;
    G_nx[1] = 1;
    G_nx[2] = 1;
    G_nx[3] = 1;
    G_nx[4] = 1;

  } // deal with D=5**04/01/2020
  else {
    tmp = 0;
    for (i = 1; i < cf_num; i++) {
      G_nx[i] = ceil(Gp[i] * nx);
      tmp += G_nx[i];
    }
    G_nx[0] = nx - tmp;
  }

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}
