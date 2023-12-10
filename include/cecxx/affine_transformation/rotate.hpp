#pragma once

#include "cecxx/types.hpp"
#include <iostream>
#include <mdspan/mdspan.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>

namespace cecxx::affine {
using namespace cecxx::types;

template <typename Matrix> struct rotate {
  rotate(Matrix rotate_matrix) : rotate_{std::move(rotate_matrix)} {}
  auto operator()(numeric_range auto &&xs) {
    auto span_mat = Kokkos::mdspan{
        rotate_.data(), Kokkos::extents{rotate_.rows(), rotate_.cols()}};
    for (int i{0}; i < rotate_.rows(); ++i) {
      for (int j{0}; j < rotate_.cols(); ++j) {
        std::cout << span_mat[i, j];
      }
      std::cout << std::endl;
    }

    auto rows_num =
        ranges::views::iota(0, static_cast<int>(span_mat.extent(0)));
    return rows_num | ranges::views::transform([span_mat, &xs](auto &&row_num) {
             auto mat_row =
                 Kokkos::submdspan(span_mat, row_num, Kokkos::full_extent);
             auto mat_row_span =
                 std::span{mat_row.data_handle(), mat_row.size()};
             std::cout << "mat row span -> " << mat_row_span[0]
                       << mat_row_span[1] << mat_row_span[2] << std::endl;
             return ranges::inner_product(mat_row_span, xs, 0);
           });
  }

  Matrix rotate_{};
};

} // namespace cecxx::affine
