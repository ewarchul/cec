#pragma once

#include <Eigen/Eigen>
#include <fmt/format.h>

#include <array>
#include <cec/specification.hpp>
#include <cstdio>
#include <experimental/__p0009_bits/full_extent_t.hpp>
#include <experimental/__p2630_bits/submdspan.hpp>
#include <iostream>
#include <nlohmann/json.hpp>
#include <range/v3/all.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/zip.hpp>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "types.h"

namespace fs = std::filesystem;
namespace rv = ranges::views;

using table_t = std::vector<std::vector<real>>;

using fn_idx_t = u8;
using dim_t = u8;
using fn_geom_data_t = std::unordered_map<fn_idx_t, table_t>;
using dim_fn_geom_data_t = std::unordered_map<dim_t, fn_geom_data_t>;
using trans_geom_mask_t = std::unordered_map<fn_idx_t, geom_transform_info>;

using eval_fn_t = std::function<real(std::span<const real>)>;
using benchmark_t = std::unordered_map<fn_idx_t, eval_fn_t>;

inline auto sphere_fn(std::span<const real> input) -> real {
  real result = 0.0;
  for (double i : input) {
    result += i * i;
  }
  return result;
}

const auto cec_2014_problems = benchmark_t{{1, sphere_fn}};

inline auto dispatch_problem_set(cec_version ver) -> benchmark_t {
  switch (ver) {
  case cec_version::CEC2014:
    return cec_2014_problems;
  }
  throw std::runtime_error("Unreachable");
}

template <auto TableType> constexpr auto get_table_path(const u8 fn_num, const u8 dimension) -> std::string {
  switch (TableType) {
  case table_type::rotate_table:
    return fmt::format("rot/M_{}_D{}.txt", fn_num, dimension);
  case table_type::shift_table:
    return fmt::format("shift/shift_data_{}.txt", fn_num);
  case table_type::shuffle_table:
    return fmt::format("shuffle/shuffle_data_{}_D{}.txt", fn_num, dimension);
  }
  throw std::runtime_error("Unreachable");
}

template <auto TableType> auto load_cec_table(const fs::path& data_storage, const u8 fn_num, const u8 dimension) -> table_t {
  const auto table_path = data_storage / get_table_path<TableType>(fn_num, dimension);
  const auto fd = std::fopen(table_path.c_str(), "r");
  if (fd == nullptr) {
    throw std::runtime_error(fmt::format("Cannot open file [{}].", table_path.string()));
  }
  table_t output{};
  output.resize(dimension, std::vector<f64>(dimension));
  for (auto row{0}; row < dimension; ++row) {
    for (auto col{0}; col < dimension; ++col) {
      std::fscanf(fd, "%lf", &output[row][col]);
    }
  }
  std::fclose(fd);
  return output;
}

template <auto TableType> auto load_fn_data(const cec_spec& spec) -> fn_geom_data_t {
  fn_geom_data_t output{};
  auto fns = rv::closed_iota(1u, spec.total_fn_num_);
  for (const auto& fn : fns) {
    output[fn] = load_cec_table<TableType>(spec.data_storage_, fn, 100);
  }
  return output;
}

template <auto TableType> auto load_dim_fn_data(const cec_spec& spec) -> dim_fn_geom_data_t {
  dim_fn_geom_data_t output{};
  auto zipped = rv::cartesian_product(spec.dimensions_, rv::closed_iota(1u, spec.total_fn_num_));
  for (const auto& [dim, fn] : zipped) {
    output[dim][fn] = load_cec_table<TableType>(spec.data_storage_, fn, dim);
  }
  return output;
}

inline auto to_geom_info(const nlohmann::json& js) {
  return geom_transform_info{.do_shift_ = js["affine_mask"]["shift"].get<bool>(),
                             .do_rotate_ = js["affine_mask"]["rotate"].get<bool>(),
                             .do_scale_ = js["affine_mask"]["scale"].get<bool>(),
                             .scale_mul_ = js["scale_mul"].get<double>()};
}

template <auto Size> auto load_json(const fs::path& filepath) -> nlohmann::json {
  auto fd = std::fopen(filepath.c_str(), "r");
  if (fd == nullptr) {
    throw std::runtime_error(fmt::format("Cannot open JSON file [{}].", filepath.string()));
  }
  std::array<char, Size> buf{};
  std::fread(buf.data(), Size, Size, fd);
  if (std::ferror(fd)) {
    throw std::runtime_error(fmt::format("IO error occured during reading of JSON file [{}]", filepath.string()));
  }
  return nlohmann::json::parse(std::begin(buf), std::end(buf));
}

inline auto load_mask(const cec_spec& spec) -> trans_geom_mask_t {
  const auto geom_masks = load_json<1024>(spec.data_storage_ / "affine_info.json");
  auto mapped = geom_masks.get<std::unordered_map<std::string, nlohmann::json>>();
  trans_geom_mask_t output{};
  for (const auto& [prob, geom_info] : mapped) {
    output.insert({std::stoi(std::string{prob}), to_geom_info(geom_info)});
  }
  return output;
}

struct cec_ctx {
  cec_ctx(const cec_spec& spec) {
    shift_ = load_fn_data<table_type::shift_table>(spec);
    rotate_ = load_dim_fn_data<table_type::rotate_table>(spec);
    mask_ = load_mask(spec);
  }
  fn_geom_data_t shift_;
  dim_fn_geom_data_t rotate_;
  trans_geom_mask_t mask_;
};

auto shift(auto& input, std::span<const real> shift_vec) -> void {
  for (u8 col = 0; col < input.size(); ++col) {
    input[col] = input[col] - shift_vec[col];
  }
}

// auto rotate(std::span<const real> input, span_2d<const real> rot_mat) {
//   std::vector<real> rotated(input.size());
//   for (u8 row = 0; row < rot_mat.extent(0); ++row) {
//     for (u8 col = 0; col < rot_mat.extent(1); ++col) {
//       rotated[row] += input[row] * rot_mat[row, col];
//     }
//   }
//   return rotated;
// }

auto scale(auto& input, const real scale_mul) -> void {
  for (u8 col = 0; col < input.size(); ++col) {
    input[col] = input[col] * scale_mul;
  }
}

auto apply_geom_transforms(auto&& input, u8 fn_index, const cec_ctx& ctx) {
  const auto dim = input.size();
  auto trans_info = ctx.mask_.at(fn_index);
  auto scale_mul = ctx.mask_.at(fn_index).scale_mul_;
  auto shift_vec = std::span<const real>{ctx.shift_.at(fn_index).data()->data(), 1000};

  if (trans_info.do_scale_) [[likely]] {
    scale(input, trans_info.scale_mul_);
  }

  if (trans_info.do_shift_) [[likely]] {
    shift(input, shift_vec); 
  }

  return input;
}

// basic_fn := x -> geom_trans -> fn
// hybrid_fn := x ->  geom_trans -> [fn] 
// complex_fn := x -> [[basic_fn | hybrid_fn]] -> complex_trans

class cec_session {
 public:
  cec_session(cec_version ver) : spec_{ver}, ctx_{spec_}, fns_{dispatch_problem_set(ver)} {}
  auto eval(const fn_idx_t fn, span_2d<const real> input) -> std::vector<real> {
    std::vector<real> output(input.extent(0));
    for (auto col{0u}; col < input.extent(1); ++col) {
      auto row = Kokkos::submdspan(input, Kokkos::full_extent, col);
      Eigen::VectorXd vec = Eigen::VectorXd::Ones(row.size());
      Eigen::VectorXd transformed = apply_geom_transforms(std::move(vec), fn, ctx_);
      output[col] = fns_[fn](std::span<const real>{row.data_handle(), row.size()});
    }
    return output;
  }

 private:
  cec_spec spec_;
  cec_ctx ctx_;
  benchmark_t fns_;
};
