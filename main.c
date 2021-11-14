#include "cec.h"
#include "utils.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define REP 1

CecData cd = {
    .prevDimension = 0,
    .prevFunction = 0,
    .dataLoaded = 0,
};

int main() {
  cec_load_benchmark_data(CEC_2017);
//  double input[10] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
 // double result = sphere_func_modern(10, input);

 // cec_external_data_t file_data = cec_read_external_data(CEC_2017, 10, 1);
  //printf("size of loaded data = %zu\n", file_data.size_);
//  cec_problem_data_t benchmark_data = cec_load_problems_data(CEC_2017, 30, 10);
 // int dims[4] = { 10 };
//  cec_load_affine_data(CEC_2017, 30, 1, dims);
//  printf("size of loaded data = %zu\n", benchmark_data.size_);
 // printf("dim 10 -> data[problem = 30][0] = %LF\n", benchmark_data.data_[29].data_[0]);
 // printf("Result = %f\n", result);
  return 0;
}

/*int main() {
  int year[] = { 2014, 2015, 2017, 2019, 2021 };
  char dataPath[50];
  int fn_max[] = { 30, 15, 30, 10, 10 };
  int dims[1] = {10};
  for (int y = 0; y < 2; ++y) {
    for (int d = 0; d < 1; ++d) {
      sprintf(dataPath, "../data/cec%d", year[y]);
      for (int fn = 1; fn < fn_max[y] + 1; ++fn) {
        for (int i = 0; i < REP; i++) {
          double *output = malloc(dims[d] * sizeof(double));
          switch (year[y]) {
          case 2014:
            cec2014_interface(dataPath, input, output, dims[d], 1, fn);
            break;
          case 2015:
            cec2015_interface(dataPath, input, output, dims[d], 1, fn);
            break;
          case 2017:
            cec2017_interface(dataPath, input, output, dims[d], 1, fn);
            break;
          case 2019:
            cec2019_interface(dataPath, input, output, dims[d], 1, fn);
            break;
          case 2021:
            cec2021_interface(dataPath, input, output, dims[d], 1, fn,
                              "bias_shift_rot");
            break;
          }
          free(output);
        }
      }
    }
  }
  free(cd.M);
  free(cd.OShift);
  free(cd.SS);
  free(cd.bias);
  return 0;
}*/
