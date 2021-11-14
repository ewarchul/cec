cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME cec)


project(${PROJECT_NAME} VERSION 1.0
             DESCRIPTION "Re-implementation of CECs benchmarks in C"
             LANGUAGES C
)

include(CTest)
enable_testing()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(EXECUTABLE ${PROJECT_NAME})
set(TEST ${PROJECT_NAME}-test)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${PROJECT_SOURCE_SIR}/build)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

set(EXTRA_LIBS m)


find_path(CMOCKA_INCLUDE_DIR
  NAMES
  cmocka.h
  PATHS
  $(CMOCKA_INCLUDE_DIR)/include
)

find_library(CMOCKA_LIBRARY
  NAMES
  cmocka cmocka_shared
  PATHS
  ${CMOCKA_ROOT_DIR}/include
)

if (CMOCKA_LIBRARY)
  set(CMOCKA_LIBRARIES
    ${CMOCKA_LIBRARIES}
    ${CMOCKA_LIBRARY}
  )
endif (CMOCKA_LIBRARY)


include_directories(${CMOCKA_INCLUDE_DIR} include)

file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/src/*.c")
file(GLOB_RECURSE INC_FILES "${PROJECT_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE ORIGNAL_CEC_SRC "${PROJECT_SOURCE_DIR}/test/cecs/*.c")
file(GLOB_RECURSE ORIGNAL_CEC_INC "${PROJECT_SOURCE_DIR}/test/cecs/*.h")


add_executable(${EXECUTABLE} main.c ${SRC_FILES} ${INC_FILES})
target_compile_options(${EXECUTABLE}
  PUBLIC
  -Wall
  -fpic
  -fstack-protector-strong
  -O3
  -Wextra
  -pedantic
  -Werror
)

target_link_libraries(${EXECUTABLE} ${CMOCKA_LIBRARIES} ${EXTRA_LIBS} -fno-omit-frame-pointer)

list(APPEND tests_names "test_cec2014")
list(APPEND cecs_ver "2014")

list(APPEND tests_names "test_cec2015")
list(APPEND cecs_ver "2015")

list(APPEND tests_names "test_cec2017")
list(APPEND cecs_ver "2017")

list(APPEND tests_names "test_cec2019")
list(APPEND cecs_ver "2019")

list(APPEND tests_names "test_cec2021")
list(APPEND cecs_ver "2021")

list(LENGTH tests_names count)
math(EXPR count "${count} -1")
foreach(i RANGE ${count})
  list(GET tests_names ${i} test_name)
  list(GET cecs_ver ${i} cec_ver)
  add_executable(
    ${test_name}
    ${PROJECT_SOURCE_DIR}/test/${test_name}.c
    ${PROJECT_SOURCE_DIR}/test/cecs/cec${cec_ver}.c
    ${SRC_FILES}
    ${INC_FILES}
  )
target_link_libraries(${test_name} ${CMOCKA_LIBRARIES} ${EXTRA_LIBS})
target_compile_options(${test_name} PUBLIC -Wall -O3 -Wextra -pedantic)
add_test(${test_name} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name}")
endforeach()



