#pragma once

#include "cec/types.h"
#include "cec/allocator.h"
#include "cec/context.h"

f64* shift(span_f64 input, span_f64 shift_vec, arena_allocator* allocator);
f64* rotate(span_f64 input, mdspan_f64 rotate_matrix, arena_allocator* allocator);
f64* scale(span_f64 input, f64 scale_mul, arena_allocator* allocator);
f64* scale_shift_rotate(span_f64 input, cec_ctx* ctx, arena_allocator* allocator);
