#pragma once

#include "cecxx/affine_transformation/identity.hpp"
#include "cecxx/types.hpp"
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/transform.hpp>

namespace cecxx::unimodal {
using namespace cecxx::types;
namespace r = ranges;
namespace rv = ranges::views;

template <numeric_range R, typename AffineFun = affine::identity>
[[nodiscard]] constexpr auto sphere(R &&xs, AffineFun &&fn = AffineFun{}) {
  const auto affined = fn(std::forward<decltype(xs)>(xs));
  const auto ys = affined | rv::transform([](auto &&x) { return x * x; });
  return r::accumulate(ys, 0.0, r::plus{});
}

} // namespace cecxx::unimodal
