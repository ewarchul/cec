#include "trans.h"
#include <string.h>

f64* shift(span_f64 input, span_f64 shift_vec, arena_allocator* allocator) {
  f64* output = arena_alloc(allocator, f64, input.size);
  for (u8 col = 0; col < input.size; ++col) {
    output[col] = input.data[col] - shift_vec.data[col];
  }
  return output;
}

f64* rotate(span_f64 input, mdspan_f64 rotate_matrix,
            arena_allocator* allocator) {
  f64* output = arena_alloc(allocator, f64, input.size);
  for (u8 col = 0; col < rotate_matrix.ncol; ++col) {
    for (u8 row = 0; row < rotate_matrix.nrow; ++row) {
      output[row * rotate_matrix.nrow + col] +=
          input.data[row] * rotate_matrix.data[row * rotate_matrix.nrow + col];
    }
  }
  return output;
}

f64* scale(span_f64 input, f64 scale_mul, arena_allocator* allocator) {
  f64* output = arena_alloc(allocator, f64, input.size);
  for (u8 col = 0; col < input.size; ++col) {
    output[col] = input.data[col] * scale_mul;
  }
  return output;
}

f64* scale_shift_rotate(span_f64 input, cec_ctx* ctx,
                        arena_allocator* allocator) {
  (void)ctx;
  f64* output = arena_alloc(allocator, f64, input.size);
  memcpy(output, input.data, input.size * sizeof(f64));
  u8                    dim         = input.size;
  span_f64              shift_vec_s = {dim, ctx->shift_table};
  mdspan_f64            rot_mat_s   = {dim, dim, ctx->rotate_table};
  f64                   scale_mul   = ctx->scale_mul;
  affine_transform_mask mask        = ctx->mask;
  if (mask.shift) {
    output = shift(SPAN_F64{input.size, output}, shift_vec_s, allocator);
  }
  if (mask.scale) {
    output = scale(SPAN_F64{input.size, output}, scale_mul, allocator);
  }
  if (mask.rotate) {
    output = rotate(SPAN_F64{input.size, output}, rot_mat_s, allocator);
  }

  return output;
}
