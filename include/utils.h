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

#define BUFFER_SIZE 1024
char error_buf[BUFFER_SIZE];
char error_str[BUFFER_SIZE];

#define log_error(MSG, ...)                                                    \
  {                                                                            \
    snprintf(error_str, (BUFFER_SIZE - 1), "[ERROR] (%s:%s:%i) ", __FILE__,    \
             __func__, __LINE__);                                              \
    snprintf(error_buf, (BUFFER_SIZE - 1), MSG, ##__VA_ARGS__);                \
    strcat(error_str, error_buf);                                              \
    puts(error_str);                                                           \
  }

cec_state_t cec_init(cec_version_t version, cec_suite_t suite);
cec_external_data_t cec_read_external_data(cec_version_t version, int dim,
                                           int problem, cec_affineT_type_t tag);
cec_problem_data_t cec_load_problems_data(cec_version_t version,
                                          int problem_nums, int dim,
                                          cec_affineT_type_t tag);
cec_problem_data_t *cec_load_affine_data(cec_version_t version,
                                         cec_benchmark_info_t info,
                                         cec_affineT_type_t tag);

cec_benchmark_data_t cec_load_benchmark_data(cec_version_t version);
void cec_print_state(cec_state_t state);
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
