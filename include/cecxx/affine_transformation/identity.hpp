#pragma once
#include "cecxx/types.hpp"

namespace cecxx::affine {
using namespace cecxx::types;

struct identity {
  constexpr auto &&operator()(types::numeric_range auto &&t) {
    return (decltype(t) &&)t;
  }
};

struct dummy_transformation {
  dummy_transformation(double mm) : mul_{mm} {}
  constexpr auto operator()(numeric_range auto &&t) {
    t[0] *= mul_;
    return t;
  }

  double mul_;
};
} // namespace cecxx::affine
