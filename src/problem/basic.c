#include "basic.h"

f64 sphere_fn(span_f64 input) {
  f64 result = 0.0;
  for (u8 i = 0; i < input.size; ++i) {
    result += input.data[i] * input.data[i];
  }
  return result;
}
