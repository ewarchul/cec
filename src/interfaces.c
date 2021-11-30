#include "interfaces.h"
#include "affine_trans.h"
#include "basic_funcs.h"
#include "types.h"
#include "utils.h"

numeric cec_eval(int fn, numeric *input, cec_state_t *state) {
  numeric output = 1;
  switch (state->version_) {
  case CEC_2017:
    output = cec_interface_2017(fn, input, state);
    break;
  }
  return output;
}

numeric cec_interface_2017(int fn, numeric *input, cec_state_t *state) {
  numeric output = 0;
  double opt_vals[30] = {100,  200,  300,  400,  500,  600,  700,  800,
                         900,  1000, 1100, 1200, 1300, 1400, 1500, 1600,
                         1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400,
                         2500, 2600, 2700, 2800, 2900, 3000};
  cec_affine_transforms_t affine_transforms = {
      .rotate_ = true, .shift_ = true, .transform_rate_ = 1};
  double *shiftrot = shift_rotate_modern(input, fn, state, affine_transforms);
  switch (fn) {
  case 1: {
    output = bent_cigar_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 2: {
    output = sum_diff_pow_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 3: {
    output = zakharov_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 4: {
    output = rosenbrock_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 5: {
    output = rastrigin_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 6: {
    output = schaffer_F7_func_modern(state->dimension_, input, shiftrot);
    break;
  }
  case 7: {
    output = bent_cigar_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 8: {
    output = step_rastrigin_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 9: {
    output = levy_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 10: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 11: {
    output = cec2017_hf01_modern(state->dimension_, shiftrot, state);
    break;
  }
  case 12: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 13: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 14: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 15: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 16: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 17: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 18: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 19: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  case 20: {
    output = schwefel_func_modern(state->dimension_, shiftrot);
    break;
  }
  }
  return output + opt_vals[fn];
}

CecData cd = {
    .prevDimension = 0,
    .prevFunction = 0,
    .dataLoaded = 0,
};

void cec2014_interface(char *datapath, double *x, double *f, int nx, int mx,
                       int func_num) {
  if (!(nx == 2 || nx == 10 || nx == 20 || nx == 30 || nx == 50 || nx == 100)) {
    perror("Error: Test functions are only defined for D = 2, 10, 20, 30, 50, "
           "100.");
  }
  if (nx == 2 && ((func_num >= 17 && func_num <= 22) ||
                  (func_num >= 29 && func_num <= 30))) {
    perror("Error: hf0{1..6}, cf0{7..8} are NOT defined for D=2.");
  }

  int shuffleFlag =
      ((func_num >= 17 && func_num <= 22) || (func_num == 29 || func_num == 30))
          ? 1
          : 0;
  if (cd.dataLoaded == 1) {
    if ((cd.prevDimension != nx) || (cd.prevFunction != func_num)) {
      cd.dataLoaded = 0;
    }
  }
  if (!cd.dataLoaded) {
    free(cd.M);
    free(cd.OShift);
    if (shuffleFlag) {
      free(cd.SS);
      loadShuffleData(&cd, datapath, nx, func_num, 2014);
    }
    loadMatrixData(&cd, datapath, nx, func_num, 2014);
    loadOShiftData(&cd, datapath, nx, func_num, 2014);
    cd.prevFunction = func_num;
    cd.prevDimension = nx;
    cd.dataLoaded = 1;
  }
  for (int i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      ellips_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      bent_cigar_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 200.0;
      break;
    case 3:
      discus_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 300.0;
      break;
    case 4:
      rosenbrock_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 400.0;
      break;
    case 5:
      ackley_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 500.0;
      break;
    case 6:
      weierstrass_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 600.0;
      break;
    case 7:
      griewank_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 700.0;
      break;
    case 8:
      rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 0);
      f[i] += 800.0;
      break;
    case 9:
      rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 900.0;
      break;
    case 10:
      schwefel_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 0);
      f[i] += 1000.0;
      break;
    case 11:
      schwefel_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1100.0;
      break;
    case 12:
      katsuura_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1200.0;
      break;
    case 13:
      happycat_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1300.0;
      break;
    case 14:
      hgbat_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1400.0;
      break;
    case 15:
      grie_rosen_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1500.0;
      break;
    case 16:
      escaffer6_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1600.0;
      break;
    case 17:
      cec2014_hf01(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1700.0;
      break;
    case 18:
      cec2014_hf02(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1800.0;
      break;
    case 19:
      cec2014_hf03(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1900.0;
      break;
    case 20:
      cec2014_hf04(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 2000.0;
      break;
    case 21:
      cec2014_hf05(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 2100.0;
      break;
    case 22:
      cec2014_hf06(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 2200.0;
      break;
    case 23:
      cec2014_cf01(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2300.0;
      break;
    case 24:
      cec2014_cf02(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2400.0;
      break;
    case 25:
      cec2014_cf03(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2500.0;
      break;
    case 26:
      cec2014_cf04(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2600.0;
      break;
    case 27:
      cec2014_cf05(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2700.0;
      break;
    case 28:
      cec2014_cf06(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2800.0;
      break;
    case 29:
      cec2014_cf07(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1);
      f[i] += 2900.0;
      break;
    case 30:
      cec2014_cf08(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1);
      f[i] += 3000.0;
      break;
    default:
      perror("Error: There are only 30 test functions in this test suite "
             "[CEC2014]!");
      f[i] = 0.0;
      break;
    }
  }
}

void cec2015_interface(char *datapath, double *x, double *f, int nx, int mx,
                       int func_num) {
  if (!(nx == 2 || nx == 10 || nx == 30 || nx == 50 || nx == 100)) {
    perror("Error: Test functions are only defined for D = 2, 10, 20, 30, 50, "
           "100.");
  }
  if (nx == 2 && ((func_num >= 6 && func_num <= 8) || (func_num == 10) ||
                  (func_num == 13))) {
    perror("Error: hf0{1..3}, cf0{2..5} are NOT defined for D=2.");
  }

  int cf_nums[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 5, 7, 10};
  int bShuffle[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0};
  int biasFlag = cf_nums[func_num] > 1 ? 1 : 0;
  int shuffleFlag = bShuffle[func_num] == 1 ? 1 : 0;

  if (cd.dataLoaded == 1) {
    if ((cd.prevDimension != nx) || (cd.prevFunction != func_num)) {
      cd.dataLoaded = 0;
    }
  }
  if (!cd.dataLoaded) {
    free(cd.M);
    free(cd.OShift);
    if (shuffleFlag) {
      free(cd.SS);
      loadShuffleData(&cd, datapath, nx, func_num, 2015);
    }
    if (biasFlag) {
      free(cd.bias);
      loadBiasData(&cd, datapath, func_num);
    }
    loadMatrixData(&cd, datapath, nx, func_num, 2015);
    loadOShiftData_(&cd, datapath, nx, func_num);
    cd.prevFunction = func_num;
    cd.prevDimension = nx;
    cd.dataLoaded = 1;
  }

  for (int i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      ellips_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      bent_cigar_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 200.0;
      break;
    case 3:
      ackley_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 300.0;
      break;
    case 4:
      rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 400.0;
      break;
    case 5:
      schwefel_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 500.0;
      break;
    case 6:
      cec2015_hf01(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 600.0;
      break;
    case 7:
      cec2015_hf02(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 700.0;
      break;
    case 8:
      cec2015_hf03(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 800.0;
      break;
    case 9:
      cec2015_cf01(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.bias, 1);
      f[i] += 900.0;
      break;
    case 10:
      cec2015_cf02(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, cd.bias, 1);
      f[i] += 1000.0;
      break;
    case 11:
      cec2015_cf03(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.bias, 1);
      f[i] += 1100.0;
      break;
    case 12:
      cec2015_cf04(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.bias, 1);
      f[i] += 1200.0;
      break;
    case 13:
      cec2015_cf05(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, cd.bias, 1);
      f[i] += 1300.0;
      break;
    case 14:
      cec2015_cf06(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.bias, 1);
      f[i] += 1400.0;
      break;
    case 15:
      cec2015_cf07(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.bias, 1);
      f[i] += 1500.0;
      break;
    default:
      perror("Error: There are only 15 test functions in this test suite! "
             "[CEC2015-LB]");
      f[i] = 0.0;
      break;
    }
  }
}

void cec2017_interface(char *datapath, double *x, double *f, int nx, int mx,
                       int func_num) {
  if (!(nx == 2 || nx == 10 || nx == 20 || nx == 30 || nx == 50 || nx == 100)) {
    perror("Error: Test functions are only defined for D = 2, 10, 20, 30, 50, "
           "100.");
  }
  if (nx == 2 && ((func_num >= 17 && func_num <= 22) ||
                  (func_num >= 29 && func_num <= 30))) {
    perror("Error: hf0{1..6}, cf0{7..8} are NOT defined for D=2.");
  }

  int shuffleFlag =
      ((func_num >= 11 && func_num <= 20) || (func_num == 29 || func_num == 30))
          ? 1
          : 0;

  if (cd.dataLoaded == 1) {
    if ((cd.prevDimension != nx) || (cd.prevFunction != func_num)) {
      cd.dataLoaded = 0;
    }
  }

  if (!cd.dataLoaded) {
    free(cd.M);
    free(cd.OShift);
    if (shuffleFlag) {
      free(cd.SS);
      loadShuffleData(&cd, datapath, nx, func_num, 2017);
    }
    loadMatrixData(&cd, datapath, nx, func_num, 2017);
    loadOShiftData(&cd, datapath, nx, func_num, 2017);
    cd.prevFunction = func_num;
    cd.prevDimension = nx;
    cd.dataLoaded = 1;
  }
  //  schaffer_F7_func depended on global state a bit too much...
  double *schafferF7_hack = malloc(nx * sizeof(double));
  for (int i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      bent_cigar_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      sum_diff_pow_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 200.0;
      break;
    case 3:
      zakharov_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 300.0;
      break;
    case 4:
      rosenbrock_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 400.0;
      break;
    case 5:
      rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 500.0;
      break;
    case 6:
      schaffer_F7_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1,
                       schafferF7_hack);
      f[i] += 600.0;
      break;
    case 7:
      bi_rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 700.0;
      break;
    case 8:
      step_rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 800.0;
      break;
    case 9:
      levy_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 900.0;
      break;
    case 10:
      schwefel_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1000.0;
      break;
    case 11:
      cec2017_hf01(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1100.0;
      break;
    case 12:
      cec2017_hf02(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1200.0;
      break;
    case 13:
      cec2017_hf03(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1300.0;
      break;
    case 14:
      cec2017_hf04(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1400.0;
      break;
    case 15:
      cec2017_hf05(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1500.0;
      break;
    case 16:
      cec2017_hf06(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1600.0;
      break;
    case 17:
      cec2017_hf07(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1700.0;
      break;
    case 18:
      cec2017_hf08(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1800.0;
      break;
    case 19:
      cec2017_hf09(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 1900.0;
      break;
    case 20:
      cec2017_hf10(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      f[i] += 2000.0;
      break;
    case 21:
      cec2017_cf01(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2100.0;
      break;
    case 22:
      cec2017_cf02(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2200.0;
      break;
    case 23:
      cec2017_cf03(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2300.0;
      break;
    case 24:
      cec2017_cf04(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2400.0;
      break;
    case 25:
      cec2017_cf05(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2500.0;
      break;
    case 26:
      cec2017_cf06(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2600.0;
      break;
    case 27:
      cec2017_cf07(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2700.0;
      break;
    case 28:
      cec2017_cf08(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      f[i] += 2800.0;
      break;
    case 29:
      cec2017_cf09(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1);
      f[i] += 2900.0;
      break;
    case 30:
      cec2017_cf10(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1);
      f[i] += 3000.0;
      break;
    default:
      perror("\nError: There are only 30 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
  free(schafferF7_hack);
}

void cec2019_interface(char *datapath, double *x, double *f, int nx, int mx,
                       int func_num) {
  if (!(nx == 2 || nx == 9 || nx == 10 || nx == 16 || nx == 18)) {
    perror("Error: Test functions are only defined for D=10, 9, 16, 18 \n\
          F1 is defined on D=9 \n F2 is defined on D=16 \n\
          F3 is defined on D=18 \n F4-F10 are defined on D=10.");
  }
  int externalDataFlag = func_num > 3 ? 1 : 0;

  if (cd.dataLoaded == 1 && externalDataFlag) {
    if ((cd.prevDimension != nx) || (cd.prevFunction != func_num)) {
      cd.dataLoaded = 0;
    }
  }

  if (!cd.dataLoaded && externalDataFlag) {
    free(cd.M);
    free(cd.OShift);
    loadMatrixData(&cd, datapath, nx, func_num, 2019);
    loadOShiftData(&cd, datapath, nx, func_num, 2019);
    cd.prevFunction = func_num;
    cd.prevDimension = nx;
    cd.dataLoaded = 1;
  }

  for (int i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      Chebyshev(&x[i * nx], nx, &f[i]);
      f[i] += 1.0;
      break;
    case 2:
      Hilbert(&x[i * nx], nx, &f[i]);
      f[i] += 1.0;
      break;
    case 3:
      Lennard_Jones(&x[i * nx], nx, &f[i]);
      f[i] += 1.0;
      break;
    case 4:
      rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1.0;
      break;
    case 5:
      griewank_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1.0;
      break;
    case 6:
      weierstrass_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1.0;
      break;
    case 7:
      schwefel_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1.0;
      break;
    case 8:
      escaffer6_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1.0;
      break;
    case 9:
      happycat_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1.0;
      break;
    case 10:
      ackley_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      f[i] += 1.0;
      break;
    default:
      perror("Error: There are only 10 test functions in this test suite! "
             "[CEC2019]");
      f[i] = 0.0;
      break;
    }
  }
}

void cec2021_interface(char *datapath, double *x, double *f, int nx, int mx,
                       int func_num, char *suite) {
  if (!(nx == 10 || nx == 20)) {
    perror("Error: Test functions are only defined for D = 10, 20.");
  }
  if (func_num < 1 || func_num > 10) {
    perror("Error: Test function is not defined");
  }

  int const shuffleFlag = (func_num >= 5 && func_num <= 7) ? 1 : 0;

  if (cd.dataLoaded == 1) {
    if ((cd.prevDimension != nx) || (cd.prevFunction != func_num)) {
      cd.dataLoaded = 0;
    }
  }

  if (!cd.dataLoaded) {
    free(cd.M);
    free(cd.OShift);
    if (shuffleFlag) {
      free(cd.SS);
      loadShuffleData(&cd, datapath, nx, func_num, 2021);
    }
    loadMatrixDataSuite(&cd, datapath, nx, func_num, suite);
    loadOShiftDataSuite(&cd, datapath, nx, func_num, suite);
    cd.prevFunction = func_num;
    cd.prevDimension = nx;
    cd.dataLoaded = 1;
  }

  int shiftFlag = 0;
  if (!strcmp(suite, "bias_shift") || !strcmp(suite, "shift") ||
      !strcmp(suite, "bias_shift_rot") || !strcmp(suite, "shift_rot")) {
    shiftFlag = 1;
  }

  int biasFlag = 0;
  if (!strcmp(suite, "bias") || !strcmp(suite, "bias_shift") ||
      !strcmp(suite, "bias_rot") || !strcmp(suite, "bias_shift_rot")) {
    biasFlag = 1;
  }

  for (int i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      bent_cigar_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      break;
    case 2:
      schwefel_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      break;
    case 3:
      bi_rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      break;
    case 4:
      grie_rosen_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      break;
    case 6:
      cec2021_hf02(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      break;
    case 7:
      cec2021_hf03(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1, 1);
      break;
    case 8:
      shiftFlag ? cec2021_cf01_s(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1)
                : cec2021_cf01(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      break;
    case 9:
      shiftFlag ? cec2021_cf02_s(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1)
                : cec2021_cf02(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      break;
    case 10:
      shiftFlag ? cec2021_cf03_s(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1)
                : cec2021_cf03(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
      break;
    default:
      perror("Error: There are only 10 test functions in this test suite!");
      f[i] = 0.0;
      break;
    }
  }
  f[0] += getFunctionBias(biasFlag, func_num);
}
