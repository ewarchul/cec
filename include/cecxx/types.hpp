#pragma once

#include <concepts>
namespace cecxx::types {

template <typename R>
consteval auto as(auto x)
  requires(std::convertible_to<decltype(x), R>)
{
  return static_cast<R>(x);
}
} // namespace cecxx::types
