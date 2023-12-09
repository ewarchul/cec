#include "cecxx/affine_transformation/rotate.hpp"
#include <cecxx/affine_transformation/AffineTransformation.hpp>
#include <cecxx/unimodal/Unimodal.hpp>
#include <iostream>
#include <vector>


auto main() -> int {
  auto scaler_2 = cecxx::affine::scale{2};
  auto rot_matrix = std::vector{
    std::vector{1, 2, 3},
    std::vector{4, 5, 6},
    std::vector{7, 8, 9}
  };
  //  int rot_matrix_c[2][3] = {{1, 2, 3}, {4, 5, 6}};
  auto rotator = cecxx::affine::rotate{rot_matrix};
  const auto x = std::vector<double>{1, 2, 3};
  auto xr = rotator(x);
  std::cout << "rotated: " << xr << std::endl;
  auto y = cecxx::unimodal::sphere(x, scaler_2);
  int xx[] = {2, 2, 2};
  auto yy = cecxx::unimodal::sphere(xx, scaler_2);
  std::cout << y << std::endl;
  std::cout << yy << std::endl;

  return 0;
}
