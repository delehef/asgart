{ pkgs, lib, config, inputs, ... }:

{
  cachix.enable = false;

  # https://devenv.sh/packages/
  packages = [ pkgs.git pkgs.figlet pkgs.cmake ];

  enterShell = ''
  figlet "ASGART loaded"
  '';

  languages.rust.enable = true;
}
