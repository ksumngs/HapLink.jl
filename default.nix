{ pkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/2f82431c7fdfa641f9816011286a2fa2c489eedb.tar.gz") {} }:

pkgs.mkShell {
  buildInputs = [
    pkgs.julia
    pkgs.pre-commit
    pkgs.nodejs
  ];

}
