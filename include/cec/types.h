#pragma once

#include <cstdint>
#include <mdspan/mdspan.hpp>

enum class cec_version { CEC2014 };
enum class table_type { rotate_table, shift_table, shuffle_table };

using u8  = uint8_t;
using u32 = uint32_t;
using u64 = uint64_t;

using i32 = int32_t;
using i64 = int64_t;
using f32 = float;
using f64 = double;

struct geom_transform_info {
  bool do_shift_;
  bool do_rotate_;
  bool do_scale_;
  f64  scale_mul_;
};

struct table_info {
  table_type type;
  u8         dim;
  u8         fn_num;
};

using real = f64;

template <class ValueType>
using span_2d = Kokkos::mdspan<ValueType, Kokkos::dextents<size_t, 2>,
                                   Kokkos::layout_right>;
