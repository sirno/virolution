{
  description = "Virolution flake.";

  inputs = {

    nixpkgs.url = "github:nixos/nixpkgs/release-24.11";

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

        packageDefaults = {
          name = "virolution-sequential";
          version = "0.5.0-dev";
          src = pkgs.lib.fileset.toSource {
            root = ./.;
            fileset = pkgs.lib.fileset.unions [
              ./Cargo.toml
              ./Cargo.lock
              ./build.rs
              ./examples
              ./macros
              ./src
              ./tests
            ];
          };

          nativeBuildInputs = with pkgs; [
            rust-select
          ];

          cargoLock = {
            lockFile = ./Cargo.lock;
          };
        };

      in
      {
        apps = {
          benchmark = {
            type = "app";
            program = "${self.packages.${system}.default-benchmark}/bin/bench.sh";
          };
          parallel-benchmark = {
            type = "app";
            program = "${self.packages.${system}.parallel-benchmark}/bin/bench.sh parallel";
          };
        };
        checks = {
          # run tests for sequential version
          sequential = rustPlatform.buildRustPackage (
            packageDefaults
            // {
              name = "virolution-sequential-test";
              buildType = "debug";
              doInstall = false;
            }
          );

          # run tests for parallel version
          parallel = rustPlatform.buildRustPackage (
            packageDefaults
            // {
              name = "virolution-parallel-test";
              buildFeatures = [ "parallel" ];
              buildType = "debug";
              doInstall = false;
            }
          );
        };

        packages.default = rustPlatform.buildRustPackage (
          packageDefaults
          // {
            doCheck = false;
          }
        );

        packages.parallel = rustPlatform.buildRustPackage (
          packageDefaults
          // {
            name = "virolution-parallel";
            buildFeatures = [ "parallel" ];
            doCheck = false;
          }
        );

        packages.default-benchmark = pkgs.stdenv.mkDerivation {
          name = "virolution-sequential-benchmark";
          src = packageDefaults.src;

          buildInputs = with pkgs; [
            self.packages.${system}.default
            hyperfine
            pkgs.makeWrapper
          ];

          installPhase = ''
            mkdir -p $out/bin

            # Copy the script
            cp $src/tests/benchmark.sh $out/bin/benchmark.sh
            chmod +x $out/bin/benchmark.sh

            # Wrap the script so that virolution (the parallel one) is on the PATH
            makeWrapper $out/bin/benchmark.sh $out/bin/bench.sh \
              --prefix PATH : "${self.packages.${system}.default}/bin"
          '';
        };

        packages.parallel-benchmark = pkgs.stdenv.mkDerivation rec {
          name = "virolution-parallel-benchmark";
          src = packageDefaults.src;

          buildInputs = [
            self.packages.${system}.parallel
            pkgs.hyperfine
            pkgs.makeWrapper
          ];

          installPhase = ''
            mkdir -p $out/bin

            # Copy the script
            cp $src/tests/benchmark.sh $out/bin/benchmark.sh
            chmod +x $out/bin/benchmark.sh

            # Wrap the script so that virolution (the parallel one) is on the PATH
            makeWrapper $out/bin/benchmark.sh $out/bin/bench.sh \
              --prefix PATH : "${self.packages.${system}.parallel}/bin"
          '';
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
