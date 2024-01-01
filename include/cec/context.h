#pragma once

#include "types.h"

typedef struct cec_ctx cec_ctx;
struct cec_ctx {
  u8                    fn_max;
  f64*                  shift_table;
  f64*                  rotate_table;
  f64                   scale_mul;
  affine_transform_mask mask;
};

cec_ctx* mk_ctx(byte* datapath, u8 dimension, u8 fn_max);
void     free_cec_ctx(cec_ctx* ctx);
