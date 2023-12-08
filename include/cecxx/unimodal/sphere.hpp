#pragma once

#include "cecxx/types.hpp"
#include <functional>
#include <numeric>
#include <ranges>
#include <span>
#include <type_traits>

namespace cecxx::problems::unimodal {

[[nodiscard]] constexpr auto sphere(const std::ranges::range auto &x) {
  using T = typename std::decay_t<decltype(*x.begin())>;
  return std::accumulate(std::begin(x), std::end(x), types::as<T>(0),
                         std::plus{});
}

} // namespace cecxx::problems::unimodal
