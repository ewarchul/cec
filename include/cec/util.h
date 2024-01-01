#pragma once

#include "types.h"

#include <stdlib.h>
#include <stdio.h>

static f64* generic_alloc_matrix(u8 ncol, u8 nrow) {
  return malloc(sizeof(f64) * ncol * nrow);
}

static f64* init_generic_matrix(u8 ncol, u8 nrow, f64 init_value) {
  f64* output = generic_alloc_matrix(ncol, nrow);
  for (u8 row = 0; row < nrow; ++row) {
    for (u8 col = 0; col < ncol; ++col) {
      output[row * ncol + col] = init_value;
    }
  }
  return output;
}

static void generic_free_matrix(u8 ncol, f64** matrix) {
  for (u8 col = 0; col < ncol; ++col) {
    free(matrix[col]);
  }
}

static void print_vec_f64(f64* vec, u8 nrow) {
  for (u8 row = 0; row < nrow; ++row) {
    printf("vec[%i] = %lf\n", row, vec[row]);
  }
}

static void print_matrix_f64(f64* matrix, u8 ncol, u8 nrow) {
  for (u8 row = 0; row < nrow; ++row) {
    for (u8 col = 0; col < ncol; ++col) {
      printf("vec[row=%i][col=%i] = %F\n", row, col, matrix[row * nrow + col]);
    }
  }
}
