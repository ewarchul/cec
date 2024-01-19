#include "include/cec/types.h"
#include <cstdlib>
#include <cec/cec.h>
#include <fmt/core.h>
#include <vector>

template <typename... T> auto println(fmt::format_string<T...>&& format, T&&... t) -> void {
  fmt::print("{}\n", fmt::format(std::forward<fmt::format_string<T...>>(format), std::forward<T>(t)...));
}

auto main() -> int {
  // const auto cec14 = cec{cec_version};
  // const auto input = ...;
  // const auto fn_num = ...;
  // const auto result = evaluate(cec14, fn_num, input);

  auto cec_session = cec{cec_version::CEC2014};

  return EXIT_SUCCESS;
}
