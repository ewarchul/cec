cmake_minimum_required(VERSION 3.25)
project(cecxx.tests)

find_package(ut 2.0 REQUIRED)
if(CECXX_TEST_EIGEN)
  find_package(Eigen3 REQUIRED NO_MODULE)
endif(CECXX_TEST_EIGEN)

enable_testing()

add_executable(cecxx.tests main.cpp unimodal.test.cpp)
target_link_libraries(cecxx.tests PRIVATE Boost::ut cecxx)
if(CECXX_TEST_EIGEN)
  add_definitions(-DCOMPILE_EIGEN_TESTS)
  target_link_libraries(cecxx.tests PRIVATE Eigen3::Eigen)
endif(CECXX_TEST_EIGEN)


