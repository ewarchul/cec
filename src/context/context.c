#include "cec/context.h"
#include "cec/allocator.h"
#include "util/string.h"
#include "io/io.h"

cec_ctx* mk_ctx(byte* datapath, u8 dimension, u8 fn_max) {
  cec_ctx* ctx       = generic_alloc(cec_ctx, fn_max);
  ctx->fn_max        = fn_max;
  char json_path[32] = "/affine_info.json";
  str  ext_path      = concat((str){datapath, 32}, (str){json_path, 32});
  for (u8 fn = 0; fn < fn_max; ++fn) {
    ctx[fn].rotate_table = read_table_data(
        datapath,
        (table_info){.dim = dimension, .fn_num = fn + 1, .type = rotate_table});
    ctx[fn].shift_table = read_table_data(
        datapath,
        (table_info){.dim = dimension, .fn_num = fn + 1, .type = shift_table});
    ctx[fn].scale_mul = read_scale_mul(ext_path.data, fn + 1);
    ctx[fn].mask      = read_affine_mask(ext_path.data, fn + 1);
  }
  free(ext_path.data);
  return ctx;
}

void free_cec_ctx(cec_ctx* ctx) {
  for (u8 i = 0; i < ctx[i].fn_max; ++i) {
    free(ctx[i].rotate_table);
    free(ctx[i].shift_table);
  }
  free(ctx);
}
