#pragma once

#include "cec/types.h"

f64*                  read_table_data(byte* datapath, table_info info);
f64                   read_scale_mul(char* filepath, u8 fn_num);
affine_transform_mask read_affine_mask(const char* filepath, u8 fn_num);
