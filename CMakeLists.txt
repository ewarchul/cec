cmake_minimum_required(VERSION 3.25)

project(cec LANGUAGES C VERSION 0.0.1)
set(CMAKE_C_STANDARD 17 CACHE STRING "Default value for C_STANDARD property of targets.")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(json-c CONFIG)

set(CEC_SOURCES 
  src/io/io.c
  src/util/allocator.c
  src/context/context.c
  src/geometry/trans.c
  src/session/session.c
  src/session/benchmark_info.c
  src/problem/basic.c
  src/benchmark/cec2014.c
  src/cec.c
)

add_library(cec SHARED ${CEC_SOURCES})
target_include_directories(cec PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src
)
target_link_libraries(cec json-c::json-c)

target_compile_options(cec PRIVATE
  -O3
  -Werror
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

add_executable(main main.c)
target_link_libraries(main cec)
target_link_options(main PRIVATE -fsanitize=address)
target_compile_options(main PRIVATE -O3)




