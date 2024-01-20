cmake_minimum_required(VERSION 3.25)

project(cec LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 23 CACHE STRING "Default value for CXX_STANDARD property of targets.")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


find_package(fmt REQUIRED)
find_package(range-v3 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(mdspan REQUIRED)

add_library(cec INTERFACE ${CEC_SOURCES})
target_include_directories(cec INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)


add_executable(main main.cxx)
target_compile_features(main PRIVATE cxx_std_23)
target_link_libraries(main cec Eigen3::Eigen std::mdspan fmt::fmt nlohmann_json::nlohmann_json range-v3::range-v3)
target_compile_options(main PUBLIC
  -O0
  #  -Werror
  -Wall
  -Wextra
  -Wempty-body
  -Wformat-security
  -Wformat
  -Winit-self
  -Wunused-variable
  -Wpedantic
  -Wno-unused-function
)





