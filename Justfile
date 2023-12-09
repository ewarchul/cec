compiler := "g++"
cpu_cores := "15"

alias c := clean
alias i := init
alias b := build
alias r := run
alias t := test

clean: 
  rm -rf build-*
  mkdir build-{{compiler}}

init:
  cmake -B build-{{compiler}} -S . -DCMAKE_CXX_COMPILER={{compiler}}

build: init
  cmake --build build-{{compiler}} -j {{cpu_cores}}

run:
  ./build-{{compiler}}/main

test:
  ./build-{{compiler}}/test/cecxx.tests
