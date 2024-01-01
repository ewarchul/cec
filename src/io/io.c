#include "io.h"
#include "consts.h"
#include "util/json.h"
#include "util/string.h"
#include "cec/allocator.h"
#include <json-c/json.h>

PRIVATE byte* get_table_format(table_type t) {
  switch (t) {
    case rotate_table:
      return "%s/rot/M_%d_D%d.txt";
    case shift_table:
      return "%s/shift/shift_data_%d.txt";
    case shuffle_table:
      return "%s/shuffle/shuffle_data_%d_D%d.txt";
  }
  perror("Unknown table type. Aborting...");
  abort();
}

PRIVATE filename_t get_table_path(byte* datapath, table_info info) {
  filename_t fname;
  byte*      format = get_table_format(info.type);
  sprintf(fname.name, format, datapath, info.fn_num, info.dim);
  return fname;
}

f64* read_table_data(byte* datapath, table_info info) {
  filename_t table_path = get_table_path(datapath, info);
  f64*       table    = generic_alloc(f64, info.dim * info.dim * CEC_SIZE_MUL);
  FILE*      fhandler = fopen(table_path.name, "r");
  if (fhandler == NULL) {
    perror("Failed to open matrix file");
    exit(EXIT_FAILURE);
  }
  for (u32 row = 0; row < info.dim * CEC_SIZE_MUL; ++row) {
    for (u32 col = 0; col < info.dim; ++col) {
      fscanf(fhandler, "%lf", &table[row * info.dim + col]);
    }
  }
  fclose(fhandler);
  return table;
}

f64 read_scale_mul(char* filepath, u8 fn_num) {
  json_object* affine_info = read_json(filepath);
  str          fn_str      = to_str(fn_num);
  json_object* problem     = get_key_js(affine_info, fn_str.data);
  f64          scale_mul   = get_key(f64, problem, "scale_mul");
  json_object_put(affine_info);
  free(fn_str.data);
  return scale_mul;
}

affine_transform_mask read_affine_mask(const char* filepath, u8 fn_num) {
  json_object*          affine_info = read_json(filepath);
  str                   fn_str      = to_str(fn_num);
  json_object*          problem     = get_key_js(affine_info, fn_str.data);
  json_object*          affine_mask = get_key_js(problem, "affine_mask");
  affine_transform_mask mask        = {
             .shift  = get_key(bool, affine_mask, "shift"),
             .rotate = get_key(bool, affine_mask, "rotate"),
             .scale  = get_key(bool, affine_mask, "scale"),
  };
  json_object_put(affine_info);
  free(fn_str.data);
  return mask;
}
