#pragma once

#include "cecxx/affine_transformation/identity.hpp"
#include "cecxx/types.hpp"
#include <cmath>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/view/transform.hpp>

namespace cecxx::unimodal {
using namespace cecxx::types;
namespace r = ranges;
namespace rv = ranges::views;

template <numeric_range R, typename AffineFun = affine::identity>
[[nodiscard]] constexpr auto bent_cigar(R &&xs, AffineFun &&fn = AffineFun{}) {
  const auto affined = fn(std::forward<decltype(xs)>(xs));
  const auto head = affined | rv::take(1) | rv::transform([](auto &&x) { return x * x; });
  const auto tail = affined | rv::drop(1) | rv::transform([](auto &&x) {
                      return std::pow(10.0, 6.0) * x * x;
                    });
  return r::accumulate(rv::concat(head, tail), 0.0, r::plus{});
}

} // namespace cecxx::unimodal
