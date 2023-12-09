#pragma once

#include "cecxx/affine_transformation/identity.hpp"
#include <functional>
#include <numeric>
#include <ranges>

namespace cecxx::unimodal {
using namespace cecxx::types;

template <numeric_range R, affine_transformation<R> T = affine::identity>
[[nodiscard]] constexpr auto sphere(R &&xs, T &&p = T{}) {
  const auto xs_square =
      std::invoke(p, std::forward<decltype(xs)>(xs)) |
      std::ranges::views::transform([](auto &&x) { return x * x; });
  return std::reduce(xs_square.begin(), xs_square.end());
}

} // namespace cecxx::unimodal
