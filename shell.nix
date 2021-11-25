{ pkgs ? import <nixpkgs> {} }:
  pkgs.mkShell rec {
    buildInputs = with pkgs; [
      gcc
      cmake
      ninja
      gdb
      valgrind
    ];
  }
