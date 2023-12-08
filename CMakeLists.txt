cmake_minimum_required(VERSION 3.1)
set(PROJECT_NAME cecxx)

project(${PROJECT_NAME} VERSION 0.0.1 DESCRIPTION "Implementation of CECs benchmarks in C++20" LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(cecxx INTERFACE)
target_include_directories(cecxx INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>)
target_compile_features(cecxx INTERFACE cxx_std_20)


add_executable(main main.cpp)
target_link_libraries(main cecxx)
set_target_properties(main PROPERTIES CXX_STANDARD 20)







