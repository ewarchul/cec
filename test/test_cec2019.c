#define UNIT_TESTING 1

#include "cec.h"
#include "interfaces.h"
#include "cecs/cec2019.h"
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <cmocka.h>

#define REPEATS 5000

CecData cd = {
  .prevDimension = 0,
  .prevFunction = 0,
  .dataLoaded = 0,
};


inline static void fillRandomArray(size_t n, double* input) {
  for(size_t i = 0; i < n; ++i) {
    input[i] = rand();
  }
}


static void test_cec2019_equiv(void **state) {
  (void)state;
  int dims[1] = {10};
  int fnNum = 10;
  char *dataPath = "../data/cec2019";
  int row = 1;
  for (int d = 0; d < 1; ++d) {
    for (int f = 1; f < fnNum+1; ++f) {
      double *original_output = malloc(dims[d] * sizeof(double));
      double *rewritten_output = malloc(dims[d] * sizeof(double));
      double *input = malloc(dims[d] * sizeof(double));
      for(int i = 0; i < REPEATS; ++i) {
      fillRandomArray(dims[d], input);
      cec2019_interface(dataPath, input, rewritten_output, dims[d], row, f);  
      cec19_test_func(input, original_output, dims[d], row, f);
      assert_true(original_output[0] == rewritten_output[0]);
      }

      free(rewritten_output);
      free(original_output);
      free(input);
    }
  }
}



int main() {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_cec2019_equiv)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
