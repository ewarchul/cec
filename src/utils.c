#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

#define dsize 32

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
    sprintf(filepath, "../data/cec%d/shuffle_data_%d_D%d.txt", version, problem, dim);
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
    if (fscanf(fhandler, "%LF", &data[i]) == -1) {
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

double getFunctionBias(const int biasFlag, const int fnNumber) {
  double bias = 0.0;
  double fnBiasDict[10] = {100.0,  1100.0, 700.0,  1900.0, 1700.0,
                           1600.0, 2100.0, 2200.0, 2400.0, 2500.0};
  if (biasFlag) {
    bias = fnBiasDict[fnNumber - 1];
  } else {
    bias = 0.0;
  }
  return bias;
}

void loadMatrixData(CecData *cd, char *dataPath, int dim, int fn,
                    int cecVersion) {
  int funcTreshold, coeff = 0;
  if (cecVersion == 2014) {
    funcTreshold = 23;
    coeff = 10;
  } else if (cecVersion == 2015) {
    int cf_nums[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 5, 7, 10};
    funcTreshold = -1;
    coeff = cf_nums[fn];
  } else if (cecVersion == 2017) {
    funcTreshold = 20;
    coeff = 10;
  } else if (cecVersion == 2019) {
    funcTreshold = 100;
    coeff = 1;
  } else if (cecVersion == 2021) {
    funcTreshold = 7;
    coeff = 10;
  } else {
    funcTreshold = -1;
    coeff = -1;
  }
  char FileName[256];
  sprintf(FileName, "%s/M_%d_D%d.txt", dataPath, fn, dim);
  FILE *fptMData = fopen(FileName, "r");
  if (fptMData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int MatrixSize = fn < funcTreshold ? dim * dim : dim * dim * coeff;
  cd->M = malloc(MatrixSize * sizeof(double));
  if (cd->M == NULL) {
    perror("Error: there is insufficient memory available!");
  } else {
    for (int i = 0; i < MatrixSize; ++i) {
      if (fscanf(fptMData, "%lf", &cd->M[i]) == -1) {
        break;
      }
    }
  }
  fclose(fptMData);
}

void loadMatrixDataSuite(CecData *cd, char *dataPath, int dim, int fn,
                         char *suite) {
  int funcTreshold = 7;
  int coeff = 10;
  char FileName[256];
  if (!strcmp(suite, "basic") || !strcmp(suite, "bias") ||
      !strcmp(suite, "bias_shift") || !strcmp(suite, "shift")) {
    sprintf(FileName, "%s/M_%d_D%d_nr.txt", dataPath, fn, dim);
  } else {
    sprintf(FileName, "%s/M_%d_D%d.txt", dataPath, fn, dim);
  }
  FILE *fptMData = fopen(FileName, "r");
  if (fptMData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int MatrixSize = fn < funcTreshold ? dim * dim : dim * dim * coeff;
  cd->M = malloc(MatrixSize * sizeof(double));
  if (cd->M == NULL) {
    perror("Error: there is insufficient memory available!");
  } else {
    for (int i = 0; i < MatrixSize; ++i) {
      if (fscanf(fptMData, "%lf", &cd->M[i]) == -1) {
        break;
      }
    }
  }
  fclose(fptMData);
}

void loadOShiftDataSuite(CecData *cd, char *dataPath, int dim, int fn,
                         char *suite) {
  char FileName[256];
  int funcTreshold = 7;
  int coeff = 10;
  if (!strcmp(suite, "basic") || !strcmp(suite, "rot") ||
      !strcmp(suite, "bias") || !(strcmp(suite, "bias_rot"))) {
    sprintf(FileName, "%s/shift_data_%d_ns.txt", dataPath, fn);
  } else {
    sprintf(FileName, "%s/shift_data_%d.txt", dataPath, fn);
  }
  FILE *fptOShiftData = fopen(FileName, "r");
  if (fptOShiftData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int OShiftSize = fn < funcTreshold ? dim : coeff * dim;
  cd->OShift = malloc(OShiftSize * sizeof(double));
  if (cd->OShift == NULL) {
    perror("Error: there is insufficient memory available!");
  }

  if (fn < funcTreshold) {
    for (int i = 0; i < OShiftSize; ++i) {
      if (fscanf(fptOShiftData, "%lf", &cd->OShift[i]) == -1) {
        break;
      }
    }
  } else {
    for (int i = 0; i < coeff - 1; i++) {
      for (int j = 0; j < dim; j++) {
        int count = fscanf(fptOShiftData, "%lf", &cd->OShift[i * dim + j]);
        if (count == -1) {
          break;
        }
      }
      int count = fscanf(fptOShiftData, "%*[^\n]%*c");
      if (count == -1) {
        break;
      }
    }
    for (int j = 0; j < dim; j++) {
      if (fscanf(fptOShiftData, "%lf", &cd->OShift[(coeff - 1) * dim + j]) ==
          -1) {
        break;
      }
    }
  }
  fclose(fptOShiftData);
}

void loadOShiftData(CecData *cd, char *dataPath, int dim, int fn,
                    int cecVersion) {
  char FileName[256];
  int funcTreshold = 0;
  int coeff = 0;
  if (cecVersion == 2014) {
    funcTreshold = 23;
    coeff = 10;
  } else if (cecVersion == 2015) {
    int coeffs[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 5, 7, 10};
    coeff = coeffs[fn];
  } else if (cecVersion == 2017) {
    funcTreshold = 20;
    coeff = 10;
  } else if (cecVersion == 2019) {
    funcTreshold = 100;
    coeff = 1;
  } else {
    funcTreshold = -1;
    coeff = -1;
  }
  sprintf(FileName, "%s/shift_data_%d.txt", dataPath, fn);
  FILE *fptOShiftData = fopen(FileName, "r");
  if (fptOShiftData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int OShiftSize = fn < funcTreshold ? dim : coeff * dim;
  cd->OShift = malloc(OShiftSize * sizeof(double));
  if (cd->OShift == NULL) {
    perror("Error: there is insufficient memory available!");
  }

  if (fn < funcTreshold) {
    for (int i = 0; i < OShiftSize; ++i) {
      if (fscanf(fptOShiftData, "%lf", &cd->OShift[i]) == -1) {
        break;
      }
    }
  } else {
    for (int i = 0; i < coeff - 1; i++) {
      for (int j = 0; j < dim; j++) {
        int count = fscanf(fptOShiftData, "%lf", &cd->OShift[i * dim + j]);
        if (count == -1) {
          break;
        }
      }
      int count = fscanf(fptOShiftData, "%*[^\n]%*c");
      if (count == -1) {
        break;
      }
    }
    for (int j = 0; j < dim; j++) {
      if (fscanf(fptOShiftData, "%lf", &cd->OShift[(coeff - 1) * dim + j]) ==
          -1) {
        break;
      }
    }
  }
  fclose(fptOShiftData);
}

void loadOShiftData_(CecData *cd, char *dataPath, int dim, int fn) {
  int coeffs[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 5, 7, 10};
  int coeff = coeffs[fn];
  char fileName[256];
  char tmpchar;
  sprintf(fileName, "%s/shift_data_%d.txt", dataPath, fn);
  FILE *fpt = fopen(fileName, "r");
  if (fpt == NULL) {
    perror("Cannot open input file for reading");
  }
  cd->OShift = malloc(coeff * dim * sizeof(double));
  if (cd->OShift == NULL) {
    perror("Error: there is insufficient memory available!");
  }
  for (int i = 0; i < dim * coeff; i++) {
    if (fscanf(fpt, "%lf", &cd->OShift[i]) == -1) {
      break;
    }
    if (coeff > 1 && ((i + 1) % dim) == 0) {
      if (fscanf(fpt, "%c", &tmpchar) == -1) {
        break;
      }
      while (tmpchar != '\n') {
        if (fscanf(fpt, "%c", &tmpchar) == -1) {
          break;
        }
      }
    }
  }
  fclose(fpt);
}

void loadShuffleData(CecData *cd, char *dataPath, int dim, int fn,
                     int cecVersion) {
  int coeff = 0;
  int shuffleFlag = 0;
  if (cecVersion == 2014 || cecVersion == 2017 || cecVersion == 2019 ||
      cecVersion == 2021) {
    coeff = 10;
  } else if (cecVersion == 2015) {
    int cf_nums[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 5, 7, 10};
    coeff = cf_nums[fn];
  }
  if (cecVersion == 2014) {
    shuffleFlag = (fn >= 17 && fn <= 22) ? 1 : 0;
  } else if (cecVersion == 2015) {
    shuffleFlag = 0;
  } else if (cecVersion == 2017) {
    shuffleFlag = (fn >= 11 && fn <= 20) ? 1 : 0;
  } else if (cecVersion == 2021) {
    shuffleFlag = (fn >= 5 && fn <= 7) ? 1 : 0;
  }

  char FileName[256];
  sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", dataPath, fn, dim);
  FILE *fptShuffleData = fopen(FileName, "r");
  if (fptShuffleData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int ShuffleSize = shuffleFlag ? dim : coeff * dim;
  cd->SS = malloc(ShuffleSize * sizeof(int));
  for (int i = 0; i < ShuffleSize; ++i) {
    if (fscanf(fptShuffleData, "%d", &cd->SS[i]) == -1) {
      break;
    }
  }
  fclose(fptShuffleData);
}

void loadBiasData(CecData *cd, char *dataPath, int fn) {
  int coeffs[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 5, 7, 10};
  int coeff = coeffs[fn];
  char fileName[256];
  sprintf(fileName, "%s/bias_%d.txt", dataPath, fn);
  FILE *fptBiasData = fopen(fileName, "r");
  if (fptBiasData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  cd->bias = malloc(coeff * sizeof(double));
  if (cd->bias == NULL) {
    perror("Error: there is insufficient memory available!");
  } else {
    for (int i = 0; i < coeff; ++i) {
      if (fscanf(fptBiasData, "%lf", &cd->bias[i]) == -1) {
        perror("Cannot read bias matrix data.");
      }
    }
  }
  fclose(fptBiasData);
}

cec_benchmark_info_t cec_mk_benchmark_info(cec_version_t version) {
  cec_benchmark_info_t _;
  switch (version) {
  case CEC_2013: {
    int n = 10;
    int fns = 28;
    int dims[10] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    cec_benchmark_info_t info = {
        .dimension_nums_ = n, .dimensions_ = dims, .problem_nums_ = fns};
    return info;
  }
  case CEC_2014: {
    int n = 1;
    int fns = 20;
    int dims[1] = {10};
    cec_benchmark_info_t info = {
        .dimension_nums_ = n, .dimensions_ = dims, .problem_nums_ = fns};
    return info;
  }
  case CEC_2015: {
    int n = 1;
    int fns = 30;
    int dims[10] = {10};
    cec_benchmark_info_t info = {
        .dimension_nums_ = n, .dimensions_ = dims, .problem_nums_ = fns};
    return info;
  }
  case CEC_2017: {
    int n = 4;
    int fns = 30;
    int *dims = malloc(n * sizeof(int));
    dims[0] = 10;
    dims[1] = 30;
    dims[2] = 50;
    dims[3] = 100;
    cec_benchmark_info_t info = {
        .dimension_nums_ = n, .dimensions_ = dims, .problem_nums_ = fns};
    return info;
  }
  case CEC_2019: {
    int n = 1;
    int fns = 30;
    int dims[1] = {10};
    cec_benchmark_info_t info = {
        .dimension_nums_ = n, .dimensions_ = dims, .problem_nums_ = fns};
    return info;
  }
  case CEC_2021: {
    int n = 2;
    int fns = 10;
    int dims[2] = {10, 20};
    cec_benchmark_info_t info = {
        .dimension_nums_ = n, .dimensions_ = dims, .problem_nums_ = fns};
    return info;
  }
  }
  return _;
}
