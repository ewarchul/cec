#include "basic_funcs.h"
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
  cec_state_t state = cec_mk_state(CEC_2017, NONE);

  double input[10] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  double *shifted_input = shift_modern(10, input, 1, &state);

  double result = bent_cigar_func_modern(10, shifted_input);

  printf("%lf\n", result);
  return 0;
}
