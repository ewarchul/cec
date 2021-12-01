#include "affine_trans.h"
#include "basic_funcs.h"
#include "cec.h"
#include "hybrid_funcs.h"
#include "interfaces.h"
#include "utils.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define REP 1
#define BENCH_ITERS 10000
#define FN 11

extern CecData cd;

int main() {
  double input[10] = {[0 ... 9] = 1};
  double output[10];
  cec_state_t state = cec_mk_state(CEC_2017, 10, NONE);

  clock_t start_time = clock();
  double result = 0;
  for (int i = 0; i < BENCH_ITERS; ++i) {
    result = cec_eval(FN, input, &state);
  }
  clock_t end_time = clock();
  double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

  clock_t start_time_old = clock();
  for (int i = 0; i < BENCH_ITERS; ++i) {
    cec2017_interface("../data/cec2017", input, output, 10, 1, FN);
  }
  double elapsed_time_old = (double)(clock() - start_time_old) / CLOCKS_PER_SEC;
  printf("[CEC OLD] Done in %lf secs\n", elapsed_time_old);


  printf("[CEC NEW] Done in %lf secs\n", elapsed_time);

  printf("speedup [pct] = %lf\n", (elapsed_time_old / elapsed_time));
  printf("Solution: NEW -> %lf =?= OLD -> %lf\n", result, output[0]);

  return 0;
}
