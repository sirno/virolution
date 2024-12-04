{
  description = "Virolution flake.";

  inputs = {

    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";

    flake-utils.url = "github:numtide/flake-utils";

    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };

  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      rust-overlay,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
      in
      {

        packages.default = pkgs.rustPlatform.buildRustPackage rec {

          name = "virolution";
          version = "0.5.0-dev";
          src = ./.;

          nativeBuildInputs = with pkgs; [
            rust-bin.nightly.latest.default
          ];

          cargoLock = {
            lockFile = ./Cargo.lock;
          };

        };

        devShell.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            rust-bin.nightly.latest.default
          ];
        };

      }
    );

}
