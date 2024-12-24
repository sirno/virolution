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
        overlays = [
          (import rust-overlay)
        ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
        rust-select = pkgs.rust-bin.selectLatestNightlyWith (
          toolchain:
          toolchain.default.override {
            extensions = [
              "rust-analyzer"
              "rust-src"
            ];
          }
        );
        rustPlatform = pkgs.makeRustPlatform {
          cargo = rust-select;
          rustc = rust-select;
        };
      in
      {
        checks = {

          # run tests for sequential version
          sequential = rustPlatform.buildRustPackage rec {
            name = "virolution-sequential-test";
            src = ./.;

            buildType = "debug";
            cargoLock = {
              lockFile = ./Cargo.lock;
            };
            buildInputs = with pkgs; [
              rust-select
            ];

            doInstall = false;
          };

          # run tests for parallel version
          parallel = rustPlatform.buildRustPackage rec {
            name = "virolution-parallel-test";
            buildFeatures = [ "parallel" ];
            src = ./.;

            buildType = "debug";
            cargoLock = {
              lockFile = ./Cargo.lock;
            };
            buildInputs = with pkgs; [
              rust-select
            ];

            doInstall = false;
          };
        };

        packages.default = rustPlatform.buildRustPackage rec {

          name = "virolution";
          version = "0.5.0-dev";
          src = ./.;

          nativeBuildInputs = with pkgs; [
            rust-select
          ];

          cargoLock = {
            lockFile = ./Cargo.lock;
          };

        };

        packages.sequential = rustPlatform.buildRustPackage rec {

          name = "virolution-sequential";
          version = "0.5.0-dev";
          src = ./.;

          nativeBuildInputs = with pkgs; [
            rust-select
          ];

          cargoLock = {
            lockFile = ./Cargo.lock;
          };

        };

        packages.parallel = rustPlatform.buildRustPackage rec {

          name = "virolution-parallel";
          version = "0.5.0-dev";
          src = ./.;

          buildFeatures = [ "parallel" ];

          nativeBuildInputs = with pkgs; [
            rust-select
          ];

          cargoLock = {
            lockFile = ./Cargo.lock;
          };

        };

        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            rust-select
            hyperfine
          ];
        };

      }
    );

}
