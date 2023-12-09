#pragma once
#include "cecxx/types.hpp"

namespace cecxx::affine {
using namespace cecxx::types;

struct identity {
  constexpr auto &&operator()(types::numeric_range auto &&t) {
    return (decltype(t) &&)t;
  }
};
} // namespace cecxx::affine
