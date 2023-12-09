#include "cecxx/unimodal/bent_cigar.hpp"
#include <cecxx/unimodal/Unimodal.hpp>
#include <cecxx/affine_transformation/AffineTransformation.hpp>
#include <iostream>
#include <vector>

auto main() -> int {
  auto x = std::vector<double>{2, 2, 2};
  auto y = cecxx::unimodal::bent_cigar(x);

  std::cout << y << std::endl;

  return 0;
}
