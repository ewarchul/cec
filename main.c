#include <cec/cec.h>
#include <cec/session.h>
#include <cec/util.h>

#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdio.h>

int main(void) {
  u8   inputs_num = 1;
  u8   dim        = 10;
  f64* input      = init_generic_matrix(dim, inputs_num, 5);
  print_matrix_f64(input, dim, inputs_num);

  cec_session* session = create_cec_session(CEC2014, NULL);

  struct timeval st, et;
  gettimeofday(&st, NULL);
  f64* output = cec(session, 1, dim, inputs_num, input);
  gettimeofday(&et, NULL);
  int elapsed = ((et.tv_sec - st.tv_sec) * 1000000) + (et.tv_usec - st.tv_usec);
  printf("Elapsed: %d us\n", elapsed);
  print_vec_f64(output, dim);
  free(input);
  destroy_cec_session(session);

  return 0;
}
