cmake_minimum_required(VERSION 3.25)

project(cec LANGUAGES CXX VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 23 CACHE STRING "Default value for CXX_STANDARD property of targets.")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


set(CEC_SOURCES 
  # src/io/io.cxx
  # src/util/allocator.cxx
  # src/context/context.cxx
  # src/geometry/trans.cxx
  # src/session/session.cxx
  # src/session/benchmark_info.cxx
  # src/problem/basic.cxx
  # src/benchmark/cec2014.cxx
  src/cec.cxx
)

find_package(fmt REQUIRED)
find_package(range-v3 REQUIRED)
find_package(nlohmann_json REQUIRED)

add_library(cec SHARED ${CEC_SOURCES})
target_include_directories(cec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src
)
target_link_libraries(cec fmt::fmt nlohmann_json::nlohmann_json range-v3::range-v3)

target_compile_options(cec PRIVATE
  -O0
  #  -Werror
  -Wall
  -Wextra
  -Wempty-body
  -Wformat-security
  -Wformat
  -Winit-self
  -Wunused-variable
  #  -Wpedantic
  -Wno-unused-function
)

add_executable(main main.cxx)
target_link_libraries(main cec)
#target_link_options(main PRIVATE -fsanitize=address)
target_compile_options(main PRIVATE -O3)




