cmake_minimum_required(VERSION 3.25)
set(PROJECT_NAME cecxx)

project(${PROJECT_NAME} VERSION 0.0.1 DESCRIPTION "Implementation of CECs benchmarks in C++23" LANGUAGES CXX)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20 CACHE STRING "Default value for CXX_STANDARD property of targets.")
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(CECXX_BUILD_TESTS "Build the tests" ON)
option(CECXX_TEST_EIGEN "Tests function with Eigen3 types" ON)

add_library(cecxx INTERFACE)
target_include_directories(cecxx INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>)
target_compile_features(cecxx INTERFACE cxx_std_20)


add_executable(main main.cpp)
target_link_libraries(main cecxx)
set_target_properties(main PROPERTIES CXX_STANDARD 23)

message(STATUS "tests -> ${CECXX_BUILD_TESTS}")
if(CECXX_BUILD_TESTS) 
  add_subdirectory(test)
endif(CECXX_BUILD_TESTS)
