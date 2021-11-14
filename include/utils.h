#ifndef UTILS_H
#define UTILS_H

/*
 * Helper function for getting right value of bias in CEC 2021 benchmark.
 */

#include "types.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


cec_state_t cec_mk_state(cec_version_t version, cec_suite_t suite);
cec_benchmark_info_t cec_mk_benchmark_info(cec_version_t version);

typedef struct CecData CecData;
struct CecData {
  int dataLoaded;
  int prevFunction;
  int prevDimension;
  double *M;
  double *OShift;
  double *bias;
  int *SS;
};
double getFunctionBias(const int, const int);
void loadMatrixData(CecData *, char *, int, int, int);
void loadMatrixDataSuite(CecData *, char *, int, int, char *);
void loadOShiftData(CecData *, char *, int, int, int);
void loadOShiftData_(CecData *, char *, int, int);
void loadOShiftDataSuite(CecData *, char *, int, int, char *);
void loadShuffleData(CecData *, char *, int, int, int);
void loadBiasData(CecData *, char *, int);

#endif
