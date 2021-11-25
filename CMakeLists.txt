cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME cec)
project(${PROJECT_NAME} VERSION 1.0
             DESCRIPTION "Re-implementation of CECs benchmarks in C"
             LANGUAGES C
)

set(CMAKE_C_STANDARD 11)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(EXECUTABLE ${PROJECT_NAME})
set(TEST ${PROJECT_NAME}-test)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${PROJECT_SOURCE_SIR}/build)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)


include(FetchContent)
FetchContent_Declare(
  Tau
  URL https://github.com/jasmcaus/tau/archive/dev.zip
)
FetchContent_MakeAvailable(Tau)
enable_testing()

set(EXTRA_LIBS m)

include_directories(include)

file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/src/*.c")
file(GLOB_RECURSE INC_FILES "${PROJECT_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE TEST_FILES "${PROJECT_SOURCE_DIR}/tests/*.c")


add_executable(${EXECUTABLE} main.c ${SRC_FILES} ${INC_FILES})
target_compile_options(${EXECUTABLE}
  PUBLIC
  -Wall
  -fpic
  -fstack-protector-strong
  -g
  #-O3
  -Wextra
  -pedantic
)

target_link_libraries(${EXECUTABLE}  ${EXTRA_LIBS} -fno-omit-frame-pointer)


add_executable(${TEST} ${TEST_FILES} ${INC_FILES} ${SRC_FILES})
target_compile_options(${TEST}
  PUBLIC
  -g
)
target_link_libraries(${TEST} Tau ${EXTRA_LIBS})




