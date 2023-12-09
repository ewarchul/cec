#pragma once

#include <concepts>
#include <ranges>
#include <type_traits>

namespace cecxx::types {

template <typename R>
consteval auto as(auto x)
  requires(std::convertible_to<decltype(x), R>)
{
  return static_cast<R>(x);
}

template <typename T>
concept numeric_range = std::ranges::range<T> and
                        (std::integral<std::ranges::range_value_t<T>> or
                         std::floating_point<std::ranges::range_value_t<T>>);
} // namespace cecxx::types
