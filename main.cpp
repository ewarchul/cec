#include "cecxx/unimodal/sphere.hpp"
#include <cecxx/unimodal/Unimodal.hpp>
#include <iostream>
#include <vector>

auto main() -> int {
  const auto x = std::vector<double>{};
  auto y = cecxx::problems::unimodal::sphere(x);

  std::cout << y << std::endl;

  return 0;
}
