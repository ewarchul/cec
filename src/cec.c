#include "cec/cec.h"
#include "benchmark/cec2014.h"
#include "cec/session.h"

f64* cec(cec_session* session, u8 fn_num, u8 dim, u8 input_num, f64* inputs) {
  mdspan_f64 inputs_s = {dim, input_num, inputs};
  return cec2014(session, fn_num, inputs_s);
}
