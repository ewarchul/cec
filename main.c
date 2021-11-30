#include "affine_trans.h"
#include "basic_funcs.h"
#include "cec.h"
#include "utils.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define REP 1

int main() {
  cec_state_t state = cec_mk_state(CEC_2017, 10, NONE);
  double input[10] = {[0 ... 9] = 1};
  double result = cec_eval(11, input, &state);
  printf("%lf\n", result);
  return 0;
}
