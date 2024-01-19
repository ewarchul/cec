#pragma once

#include <cec/types.h>
#include <filesystem>
#include <vector>

struct cec_spec {
  cec_spec(cec_version ver) {
    switch (ver) {
      case cec_version::CEC2014: {
        total_fn_num_ = 30;
        data_storage_ = "../data/cec2014";
        dimensions_   = {10, 30, 50, 100};
      } break;
    }
  }

  std::vector<u8>       dimensions_;
  u8                    total_fn_num_;
  std::filesystem::path data_storage_;
};
