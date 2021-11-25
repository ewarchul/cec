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
  cec_state_t state = cec_mk_state(CEC_2017, 30, NONE);

  numeric input[30] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                       0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                       0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  cec_affine_transforms_t transforms = {
      .bias_ = false,
      .rotate_ = true,
      .shift_ = true,
      .shift_rotate_ = false,
      .shuffle_ = false,
      .transform_rate_ = 1,
  };

  double *sr_input = shift_rotate_modern(input, 1, &state, transforms);
  double result = bent_cigar_func_modern(30, sr_input);
  printf("%lf\n", result);

  return 0;
}
