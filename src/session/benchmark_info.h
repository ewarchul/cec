#pragma once

#include "cec/types.h"

typedef struct benchmark_info benchmark_info;
struct benchmark_info {
  u8*   dimensions;
  u8    dim_num;
  u8    problem_num;
  char* storage_dir;
};

benchmark_info mk_benchmark_info(cec_version version);
