#ifndef TYPES_H
#define TYPES_H

#include <stdbool.h>
#include <stddef.h>

#define MAX_PATH 1024

typedef struct fpath fpath;
struct fpath {
  char ptr[MAX_PATH];
};

typedef enum {
  SHIFT_SUITE,
  ROTATE_SUITE,
  BIAS_SUITE,
  BIAS_ROT_SUITE,
  SHIFT_ROTATE_SUITE,
  BIAS_SHIFT_ROTATE_SUITE
} cec_suite_t;

typedef enum {
  CEC_2013 = 2013,
  CEC_2014 = 2014,
  CEC_2015 = 2015,
  CEC_2017 = 2017,
  CEC_2019 = 2019,
  CEC_2021 = 2021
} cec_version_t;

typedef enum {
  SHIFT,
  ROT,
  SHUFFLE,
  BIAS
} cec_affineT_type_t;

typedef long double numeric;

typedef struct cec_interface_t cec_interface_t;

typedef struct cec_problems_range_t cec_problems_range_t;
struct cec_problems_range_t {
  unsigned min_;
  unsigned max_;
};

typedef struct cec_affine_transforms_t cec_affine_transforms_t;
struct cec_affine_transforms_t {
  bool shift_rotate_;
  bool shift_;
  bool rotate_;
  bool shuffle_;
  bool bias_;
  double transform_rate_;
};

typedef struct cec_external_data_t cec_external_data_t;
struct cec_external_data_t {
  size_t size_;
  numeric *data_;
  bool valid_;
};

typedef struct cec_problem_data_t cec_problem_data_t;
struct cec_problem_data_t {
  int problem_nums_;
  size_t size_;
  cec_external_data_t *data_;
};

typedef struct cec_benchmark_data_t cec_benchmark_data_t;
struct cec_benchmark_data_t {
  cec_problem_data_t *rotates_;
  cec_problem_data_t *shuffles_;
  cec_problem_data_t *shifts_;
  cec_problem_data_t *biases_;
};

typedef struct cec_state_t cec_state_t;
struct cec_state_t {
  int dimensions_;
  cec_version_t version_;
  cec_suite_t suite_;
  cec_benchmark_data_t data_;
};

typedef double (*function_t)(size_t n, double input[n]);

typedef struct cec_benchmark_info_t cec_benchmark_info_t;
struct cec_benchmark_info_t {
  int problem_nums_;
  size_t dimension_nums_;
  int *dimensions_;
};


#endif
