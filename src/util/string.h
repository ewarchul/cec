#pragma once

#include "cec/types.h"
#include <stdio.h>

static str concat(const str lhs, const str rhs) {
  u8    size        = lhs.len + rhs.len;
  byte* result_data = malloc(size);
  snprintf(result_data, size, "%s%s", lhs.data, rhs.data);
  return (str){.data = result_data, .len = size};
}

static str to_str(int input) {
  byte* parsed = malloc(32);
  sprintf(parsed, "%d", input);
  return (str){.data = parsed, .len = 32};
}
