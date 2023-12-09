#include "cecxx/unimodal/Unimodal.hpp"
#include <boost/ut.hpp>

#ifdef COMPILE_EIGEN_TESTS
#include <Eigen/Core>
#endif

namespace ut = boost::ut;

boost::ut::suite<"unimodal functions"> errors = [] {
  using namespace ut;
  "sphere_stl_vec"_test = [&] {
    std::vector<int> in{1, 2, 3};
    ut::expect(cecxx::unimodal::sphere(in) == 14);
  };

  "sphere_c_array"_test = [&] {
    int in[] = {1, 2, 3};
    ut::expect(cecxx::unimodal::sphere(std::span(in)) == 14);
  };

#ifdef COMPILE_EIGEN_TESTS
  "sphere_eigen"_test = [&] {
    Eigen::VectorXd in(3);
    in << 1, 2, 3;
    ut::expect(cecxx::unimodal::sphere(in) == 14);
  };
#endif 


  "bent_cigar"_test = [&] {
    std::vector<int> in{2, 2, 2};
    ut::expect(cecxx::unimodal::bent_cigar(in) == 8000004);
  };


};

