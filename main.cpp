#include "cecxx/unimodal/sphere.hpp"
#include <cecxx/unimodal/Unimodal.hpp>
#include <iostream>
#include <vector>

auto main() -> int {
  auto x = std::vector<double>{1, 2, 3};
  auto affine = cecxx::affine::dummy_transformation{1.1};
  auto y = cecxx::unimodal::sphere(x, affine);

  std::cout << y << std::endl;

  return 0;
}
