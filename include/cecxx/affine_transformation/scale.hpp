#pragma once

#include "cecxx/types.hpp"
#include <range/v3/view/transform.hpp>

namespace cecxx::affine {
using namespace cecxx::types;

struct scale {
  scale(double scale_mul) : scale_mul_{scale_mul} {}
  auto operator()(numeric_range auto &&xs) {
    return xs | ranges::views::transform(
                    [this](auto &&x) { return x * scale_mul_; });
  }

  double scale_mul_{};
};

} // namespace cecxx::affine
