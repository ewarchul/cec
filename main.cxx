#include <fmt/core.h>
#include <fmt/chrono.h>

#include <Eigen/Eigen>
#include <chrono>
#include <cec/cec.hpp>
#include <cstdlib>
#include <iostream>

#include "include/cec/types.h"

template <typename... T> auto println(fmt::format_string<T...>&& format, T&&... t) -> void {
  fmt::print("{}\n", fmt::format(std::forward<fmt::format_string<T...>>(format), std::forward<T>(t)...));
}

auto main() -> int {
  // initialize input for evaluation 
  constexpr int nrow{10};
  constexpr int ncol{2};
  Eigen::MatrixXd input = Eigen::MatrixXd::Ones(nrow, ncol);

  // create view for this input
  auto spanmd = span_2d<real>{input.data(), nrow, ncol};

  // create cec session 
  auto cec = cec_session{cec_version::CEC2014};

  //evaluate 
  auto start = std::chrono::system_clock::now();
  const auto result = cec.eval(1, spanmd);
  auto end = std::chrono::system_clock::now();
  println("[duration = {}] Result -> {}", end-start, result[0]);

  return EXIT_SUCCESS;
}
