#include "cecxx/affine_transformation/scale.hpp"
#include "cecxx/unimodal/bent_cigar.hpp"
#include <cecxx/unimodal/Unimodal.hpp>
#include <cecxx/affine_transformation/AffineTransformation.hpp>
#include <iostream>
#include <vector>

auto main() -> int {
  auto scaler_2 = cecxx::affine::scale{2};

  const auto x = std::vector<double>{2, 2, 2};
  auto y = cecxx::unimodal::sphere(x, scaler_2);

  int xx[] = {2, 2, 2};
  auto yy = cecxx::unimodal::sphere(xx, scaler_2);

  std::cout << y << std::endl;
  std::cout << yy << std::endl;

  return 0;
}
