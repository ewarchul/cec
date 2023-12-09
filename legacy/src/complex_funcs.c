#include "complex_funcs.h"

void cec2014_cf01(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 1 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 30, 40, 50};
  double bias[5] = {0, 100, 200, 300, 400};

  i = 0;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+4;
  i = 1;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 2;
  bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+30;
  i = 3;
  discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 4;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, 0);
  fit[i] = 10000 * fit[i] / 1e+10;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf02(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 2 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {20, 20, 20};
  double bias[3] = {0, 100, 200};

  i = 0;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, 0);
  i = 1;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 2;
  hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf03(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 3 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  double bias[3] = {0, 100, 200};
  i = 0;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 4e+3;
  i = 1;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+10;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf04(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 10, 10, 10, 10};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 4e+3;
  i = 1;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+10;
  i = 3;
  weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 400;
  i = 4;
  griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf05(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 10, 10, 20, 20};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1000;
  i = 1;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  i = 2;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 3;
  weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 400;
  i = 4;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf06(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 6 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 30, 40, 50};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  grie_rosen_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 1;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  i = 2;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 3;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 2e+7;
  i = 4;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf07(double *x, double *f, int nx, double *Os, double *Mr, int *SS,
                  int r_flag) /* Composition Function 7 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  double bias[3] = {0, 100, 200};
  i = 0;
  cec2014_hf01(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 1;
  cec2014_hf02(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 2;
  cec2014_hf03(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf08(double *x, double *f, int nx, double *Os, double *Mr, int *SS,
                  int r_flag) /* Composition Function 8 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  double bias[3] = {0, 100, 200};
  i = 0;
  cec2014_hf04(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 1;
  cec2014_hf05(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 2;
  cec2014_hf06(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2015_cf01(double *x, double *f, int nx, double *Os, double *Mr,
                  double *bias, int r_flag) /* Composition Function 1 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {20, 20, 20};

  i = 0;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, 0);
  i = 1;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 2;
  hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2015_cf02(double *x, double *f, int nx, double *Os, double *Mr, int *SS,
                  double *bias, int r_flag) /* Composition Function 2 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  i = 0;
  cec2015_hf01(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 1;
  cec2015_hf02(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 2;
  cec2015_hf03(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2015_cf03(double *x, double *f, int nx, double *Os, double *Mr,
                  double *bias, int r_flag) /* Composition Function 3 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 10, 10, 20, 20};
  i = 0;
  hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1000.0;
  i = 1;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+3;
  i = 2;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 4.0e+3;
  i = 3;
  weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 400.0;
  i = 4;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+10;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2015_cf04(double *x, double *f, int nx, double *Os, double *Mr,
                  double *bias, int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 20, 30, 30};
  i = 0;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 4.0e+3;
  i = 1;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+3;
  i = 2;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+10;
  i = 3;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 4;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+3;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2015_cf05(double *x, double *f, int nx, double *Os, double *Mr, int *SS,
                  double *bias, int r_flag) /* Composition Function 5 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 10, 10, 20, 20};
  i = 0;
  cec2015_hf03(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 1;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+3;
  i = 2;
  cec2015_hf01(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 3;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 4.0e+3;
  i = 4;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2015_cf06(double *x, double *f, int nx, double *Os, double *Mr,
                  double *bias, int r_flag) /* Composition Function 6 */
{
  int i, cf_num = 7;
  double fit[7];
  double delta[7] = {10, 20, 30, 40, 50, 50, 50};
  i = 0;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+3;
  i = 1;
  grie_rosen_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 4.0e+3;
  i = 2;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 4.0e+3;
  i = 3;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 4;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+10;
  i = 5;
  bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+10;
  i = 6;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000.0 * fit[i] / 1.0e+3;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2015_cf07(double *x, double *f, int nx, double *Os, double *Mr,
                  double *bias, int r_flag) /* Composition Function 7 */
{
  int i, cf_num = 10;
  double fit[10];
  double delta[10] = {10, 10, 20, 20, 30, 30, 40, 40, 50, 50};
  i = 0;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 100.0 * fit[i] / 1.0e+3;
  i = 1;
  weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 100.0 * fit[i] / 400.0;
  i = 2;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 100.0 * fit[i] / 1.0e+3;
  i = 3;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 100.0 * fit[i] / 4.0e+3;
  i = 4;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 100.0 * fit[i] / 1.0e+5;
  i = 5;
  hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 100.0 * fit[i] / 1000.0;
  i = 6;
  katsuura_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 100.0 * fit[i] / 1.0e+7;
  i = 7;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 8;
  grie_rosen_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 100.0 * fit[i] / 4.0e+3;
  i = 9;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 100.0 * fit[i] / 1.0e+5;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2017_cf01(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) {
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 20, 30};
  double bias[3] = {0, 100, 200};

  i = 0;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 1;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 2;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);

  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2017_cf02(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) {
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 20, 30};
  double bias[3] = {0, 100, 200};

  i = 0;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 1;
  griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 2;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2017_cf03(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) {
  int i, cf_num = 4;
  double fit[4];
  double delta[4] = {10, 20, 30, 40};
  double bias[4] = {0, 100, 200, 300};

  i = 0;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 1;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 2;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 3;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}
void cec2017_cf04(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) {
  int i, cf_num = 4;
  double fit[4];
  double delta[4] = {10, 20, 30, 40};
  double bias[4] = {0, 100, 200, 300};

  i = 0;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 1;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 2;
  griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 3;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2017_cf05(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) {
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 30, 40, 50};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  i = 1;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 3;
  discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 4;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2017_cf06(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) {
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 20, 30, 40};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 2e+7;
  i = 1;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 2;
  griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 3;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 4;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2017_cf07(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) {
  int i, cf_num = 6;
  double fit[6];
  double delta[6] = {10, 20, 30, 40, 50, 60};
  double bias[6] = {0, 100, 200, 300, 400, 500};
  i = 0;
  hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1000;
  i = 1;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  i = 2;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 3;
  bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+30;
  i = 4;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 5;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 2e+7;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2017_cf08(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) {
  int i, cf_num = 6;
  double fit[6];
  double delta[6] = {10, 20, 30, 40, 50, 60};
  double bias[6] = {0, 100, 200, 300, 400, 500};
  i = 0;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 1;
  griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 2;
  discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 3;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 4;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 5;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 2e+7;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2017_cf09(double *x, double *f, int nx, double *Os, double *Mr, int *SS,
                  int r_flag) {
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  double bias[3] = {0, 100, 200};
  i = 0;
  cec2017_hf05(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 1;
  cec2017_hf06(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 2;
  cec2017_hf07(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2017_cf10(double *x, double *f, int nx, double *Os, double *Mr, int *SS,
                  int r_flag) {
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  double bias[3] = {0, 100, 200};
  i = 0;
  cec2017_hf05(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 1;
  cec2017_hf08(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  i = 2;
  cec2017_hf09(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1,
               r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2021_cf01(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 2 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 20, 30};
  double bias[3] = {0, 0, 0};

  i = 0;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 1;
  griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 2;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2021_cf02(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 3 */
{
  int i, cf_num = 4;
  double fit[4];
  double delta[4] = {10, 20, 30, 40};
  double bias[4] = {0, 0, 0, 0};

  i = 0;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 1;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 2;
  griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 3;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2021_cf03(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 30, 40, 50};
  double bias[5] = {0, 0, 0, 0, 0};
  i = 0;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  i = 1;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 3;
  discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 4;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2021_cf01_s(double *x, double *f, int nx, double *Os, double *Mr,
                    int r_flag) /* Composition Function 2 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 20, 30};
  double bias[3] = {0, 100, 200};

  i = 0;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 1;
  griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 2;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2021_cf02_s(double *x, double *f, int nx, double *Os, double *Mr,
                    int r_flag) /* Composition Function 3 */
{
  int i, cf_num = 4;
  double fit[4];
  double delta[4] = {10, 20, 30, 40};
  double bias[4] = {0, 100, 200, 300};

  i = 0;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 1;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 2;
  griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 3;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2021_cf03_s(double *x, double *f, int nx, double *Os, double *Mr,
                    int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 30, 40, 50};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  i = 1;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  i = 3;
  discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 4;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}
