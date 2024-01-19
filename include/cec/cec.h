#pragma once

#include <array>
#include <fmt/format.h>

#include <cec/specification.hpp>
#include <cstdio>
#include <nlohmann/json.hpp>
#include <range/v3/all.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/zip.hpp>
#include <unordered_map>
#include <vector>

#include "types.h"

namespace fs = std::filesystem;
namespace rv = ranges::views;

using fn_index = u8;
using table = std::vector<std::vector<f64>>;
using geom_data_dict = std::unordered_map<fn_index, table>;
using transformation_mask_dict = std::unordered_map<fn_index, geom_transform_info>;

template <auto TableType> constexpr auto get_table_path(const u8 fn_num, const u8 dimension) -> std::string {
  switch (TableType) {
  case table_type::rotate_table:
    return fmt::format("rot/M_{}_D{}.txt", fn_num, dimension);
  case table_type::shift_table:
    return fmt::format("shift/shift_data_{}.txt", fn_num);
  case table_type::shuffle_table:
    return fmt::format("shuffle/shuffle_data_{}_D{}.txt", fn_num, dimension);
  }
}

template <auto TableType>
auto load_table(const fs::path &data_storage, const u8 dimension, const u8 fn_num) -> table {
  const auto table_path = data_storage / get_table_path<TableType>(fn_num, dimension);
  const auto fd = std::fopen(table_path.c_str(), "r");
  if (fd == nullptr) {
    throw std::runtime_error(fmt::format("Cannot open file [{}].", table_path.string()));
  }
  table output{};
  output.resize(dimension, std::vector<f64>(dimension));
  for (auto row{0}; row < dimension; ++row) {
    for (auto col{0}; col < dimension; ++col) {
      std::fscanf(fd, "%lf", &output[row][col]);
    }
  }
  std::fclose(fd);
  return output;
}

template <auto TableType> auto load_cec_table(const cec_spec &spec) -> geom_data_dict {
  geom_data_dict output{};
  auto zipped = rv::cartesian_product(spec.dimensions_, rv::closed_iota(1u, spec.total_fn_num_));
  for (const auto &[dim, prob] : zipped) {
    output[prob] = load_table<TableType>(spec.data_storage_, dim, prob);
  }
  return output;
}

inline auto to_geom_info(const nlohmann::json &js) {
  return geom_transform_info{.do_shift_ = js["affine_mask"]["shift"].get<bool>(),
                             .do_rotate_ = js["affine_mask"]["rotate"].get<bool>(),
                             .do_scale_ = js["affine_mask"]["scale"].get<bool>(),
                             .scale_mul_ = js["scale_mul"].get<double>()};
}

template <auto Size> auto load_json(const fs::path &filepath) -> nlohmann::json {
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

inline auto load_mask(const cec_spec &spec) -> transformation_mask_dict {
  const auto geom_masks = load_json<1024>(spec.data_storage_ / "affine_info.json");
  auto mapped = geom_masks.get<std::unordered_map<std::string, nlohmann::json>>();
  transformation_mask_dict output{};
  for (const auto &[prob, geom_info] : mapped) {
    output.insert({std::stoi(std::string{prob}), to_geom_info(geom_info)});
  }
  return output;
}

struct cec_ctx {
  cec_ctx(const cec_spec &spec) {
    shift_ = load_cec_table<table_type::shift_table>(spec);
    rotate_ = load_cec_table<table_type::rotate_table>(spec);
    mask_ = load_mask(spec);
  }
  geom_data_dict shift_;
  geom_data_dict rotate_;
  transformation_mask_dict mask_;
};

class cec {
public:
  cec(cec_version version) : spec_{version}, ctx_{spec_} {}
  auto eval(const fn_index fn, matrix_span<real> input);

private:
  cec_spec spec_;
  cec_ctx ctx_;
};
