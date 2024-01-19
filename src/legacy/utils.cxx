#include "utils.h"

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
