#include "io.h"

cec_external_data_t cec_read_external_data(cec_version_t version, int dim,
                                           int problem,
                                           cec_affineT_type_t tag) {
  char filepath[256];
  switch (tag) {
  case SHIFT:
    sprintf(filepath, "../data/cec%d/shift_data_%d.txt", version, problem);
    break;
  case ROT:
    sprintf(filepath, "../data/cec%d/M_%d_D%d.txt", version, problem, dim);
    break;
  case SHUFFLE:
    sprintf(filepath, "../data/cec%d/shuffle_data_%d_D%d.txt", version, problem,
            dim);
    break;
  case BIAS:
    sprintf(filepath, "../data/cec%d/bias_%d.txt", version, problem);
    break;
  }
  FILE *fhandler = fopen(filepath, "r");
  if (fhandler == NULL) {
    log_error("Error: Cannot open input file for reading. DIMENSION = %d, "
              "PROBLEM = %d\n",
              dim, problem);
    exit(EXIT_FAILURE);
  }
  size_t size = dim * dim * 10;
  numeric *data = malloc(size * sizeof(numeric));
  for (size_t i = 0; i < size; ++i) {
    if (fscanf(fhandler, "%lf", &data[i]) == -1) {
      break;
    }
  }
  fclose(fhandler);
  cec_external_data_t result = {.data_ = data, .size_ = size, .valid_ = true};
  return result;
}

cec_problem_data_t cec_load_problems_data(cec_version_t version,
                                          int problem_nums, int dim,
                                          cec_affineT_type_t tag) {
  cec_external_data_t *data =
      malloc(problem_nums * sizeof(cec_external_data_t));

  size_t total_size = 0;
  for (int p = 0; p < problem_nums; ++p) {
    data[p] = cec_read_external_data(version, dim, p + 1, tag);
    total_size += data[p].size_;
  }

  cec_problem_data_t result = {
      .problem_nums_ = problem_nums, .data_ = data, .size_ = total_size};

  return result;
}

cec_problem_data_t *cec_load_affine_data(cec_version_t version,
                                         cec_benchmark_info_t info,
                                         cec_affineT_type_t tag) {

  cec_problem_data_t *result =
      malloc(info.dimension_nums_ * sizeof(cec_problem_data_t));

  for (size_t d = 0; d < info.dimension_nums_; ++d) {
    result[d] = cec_load_problems_data(version, info.problem_nums_,
                                       info.dimensions_[d], tag);
  }
  return result;
}

cec_benchmark_data_t cec_load_benchmark_data(cec_version_t version) {
  cec_benchmark_info_t info = cec_mk_benchmark_info(version);

  cec_problem_data_t *rotates = cec_load_affine_data(version, info, ROT);
  cec_problem_data_t *shifts = cec_load_affine_data(version, info, SHIFT);
  cec_problem_data_t *shuffles = cec_load_affine_data(version, info, SHUFFLE);

  cec_benchmark_data_t result = {
      .rotates_ = rotates,
      .shifts_ = shifts,
      .shuffles_ = shuffles,
  };
  return result;
}
