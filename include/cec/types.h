#pragma once

#include <sys/types.h>
#include <stdbool.h>

#define SPAN_F64   (span_f64)
#define MDSPAN_F64 (mdspan_f64)

#define PRIVATE    static

typedef enum { CEC2014 } cec_version;
typedef enum { rotate_table, shift_table, shuffle_table } table_type;

typedef u_int8_t  u8;
typedef u_int32_t u32;
typedef u_int64_t u64;

typedef int32_t i32;
typedef float   f32;
typedef double  f64;
typedef char    byte;

typedef struct affine_transform_mask affine_transform_mask;
struct affine_transform_mask {
  bool shift;
  bool rotate;
  bool scale;
};

typedef struct filename_t filename_t;
struct filename_t {
  byte name[128];
};

typedef struct table_info table_info;
struct table_info {
  table_type type;
  u8         dim;
  u8         fn_num;
};

typedef struct span_f64 span_f64;
struct span_f64 {
  u8   size;
  f64* data;
};

typedef struct mdspan_f64 mdspan_f64;
struct mdspan_f64 {
  u8   ncol;
  u8   nrow;
  f64* data;
};

typedef f64 (*optim_problem)(span_f64 input);

typedef struct str str;
struct str {
  byte* data;
  u8    len;
};
