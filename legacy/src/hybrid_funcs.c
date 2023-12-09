#include "hybrid_funcs.h"

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
