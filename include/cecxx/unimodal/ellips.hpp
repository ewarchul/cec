#pragma once

#include "cecxx/affine_transformation/identity.hpp"
#include "cecxx/types.hpp"
#include <cmath>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>

namespace cecxx::unimodal {
using namespace cecxx::types;
namespace r = ranges;
namespace rv = ranges::views;

template <numeric_range R, typename AffineFun = affine::identity>
[[nodiscard]] constexpr auto ellips(R &&xs, AffineFun &&fn = AffineFun{}) {
  const auto affined = fn(std::forward<decltype(xs)>(xs));
  const auto n = r::size(xs);
  const auto ys = affined | rv::enumerate | rv::transform([n](auto &&zipped) {
                    const auto &[i, x] = zipped;
                    return std::pow(10.0, 6.0 * i / (n - 1)) * x * x;
                  });
  return r::accumulate(ys, 0.0, r::plus{});
}

} // namespace cecxx::unimodal
