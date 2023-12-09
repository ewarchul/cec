#pragma once

#include "cecxx/types.hpp"
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/view/transform.hpp>

namespace cecxx::affine {
using namespace cecxx::types;

template <typename Matrix> struct rotate {
  rotate(Matrix rotate_matrix) : rotate_{std::move(rotate_matrix)} {}
  auto operator()(numeric_range auto &&xs) {
    return rotate_ | ranges::views::transform([&](auto &&row) {
             return ranges::inner_product(row, xs, 0);
           });
  }

  Matrix rotate_{};
};

} // namespace cecxx::affine
