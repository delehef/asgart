{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = [
    pkgs.cmake pkgs.cargo pkgs.rustfmt
  ];
}
