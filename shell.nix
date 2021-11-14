{ pkgs ? import <nixpkgs> {} }:
  pkgs.mkShell rec {
    buildInputs = with pkgs; [
      gcc
      cmake
      cmocka
      ninja
      gdb
      valgrind
    ];
  }
