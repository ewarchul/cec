#include "benchmark_info.h"
#include "cec/allocator.h"

PRIVATE u8* get_dimensions(cec_version version) {
  u8* dimensions = generic_alloc(u8, 100);
  switch (version) {
    case CEC2014: {
      dimensions[0] = 10;
      dimensions[1] = 20;
      dimensions[2] = 30;
      dimensions[3] = 50;
      dimensions[4] = 100;
      break;
    }
  }
  return dimensions;
}

PRIVATE u8 get_dimensions_num(cec_version version) {
  switch (version) {
    case CEC2014:
      return 5;
  }
  abort();
}

PRIVATE char* get_storage_dir(cec_version version) {
  switch (version) {
    case CEC2014:
      return "../data/cec2014";
  }
  abort();
}

PRIVATE u8 get_problem_num(cec_version version) {
  switch (version) {
    case CEC2014:
      return 30;
  }
  abort();
}

benchmark_info mk_benchmark_info(cec_version version) {
  return (benchmark_info){
      .dimensions  = get_dimensions(version),
      .dim_num     = get_dimensions_num(version),
      .storage_dir = get_storage_dir(version),
      .problem_num = get_problem_num(version),
  };
}
