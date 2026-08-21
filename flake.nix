{
  description = "SnapWave — fast, implicit, unstructured-grid short wave solver";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      systems = [
        "x86_64-linux"
        # "aarch64-linux"
        # "x86_64-darwin"
        # "aarch64-darwin"
      ];
      forAllSystems = nixpkgs.lib.genAttrs systems;

      # Keep only what is needed to build the legacy Fortran oracle:
      # Makefile, Fortran/C sources and the bundled third_party code.
      # Excludes Windows project files, docs, testcases, the prebuilt zip and
      # the Windows-only netcdf binaries.
      legacySource = pkgs:
        pkgs.lib.cleanSourceWith {
          src = self;
          filter =
            path: type:
            let
              rel = pkgs.lib.removePrefix (toString self + "/") (toString path);
              top = builtins.head (pkgs.lib.splitString "/" rel);
            in
            builtins.elem top [ "Makefile" "src" "utils_lgpl" "third_party_open" ]
            && rel != "third_party_open/netcdf";
        };

      # Sources for the pure-Rust build (plan.md Phase 12): only Rust sources
      # and the manifest. No Fortran/C/NetCDF input is compiled or linked.
      cargoSource = pkgs:
        pkgs.lib.cleanSourceWith {
          src = self;
          filter =
            path: type:
            let
              rel = pkgs.lib.removePrefix (toString self + "/") (toString path);
              top = builtins.head (pkgs.lib.splitString "/" rel);
            in
            builtins.elem top [ "Cargo.toml" "Cargo.lock" "src_rust" ];
        };

      # The production binary: pure Rust, built by Cargo with no foreign
      # toolchain (plan.md Phase 12).
      mkSnapwave = pkgs:
        pkgs.rustPlatform.buildRustPackage {
          pname = "snapwave";
          version = "unstable-${self.shortRev or "dirty"}";

          src = cargoSource pkgs;
          cargoLock.lockFile = ./Cargo.lock;

          # Integration tests need the testcases/ tree and a live oracle;
          # they run via `cargo test` in the dev shell, not in the sandbox.
          doCheck = false;

          meta = {
            description = "Fast, implicit, unstructured-grid short wave solver";
            mainProgram = "snapwave";
            license = nixpkgs.lib.licenses.lgpl21;
            platforms = nixpkgs.lib.platforms.unix;
          };
        };

      # The legacy stand-alone Makefile build (argument-less Fortran
      # binary), kept as the numerical oracle: point SNAPWAVE_ORACLE at it
      # for wrapper-vs-oracle regression runs.
      mkSnapwaveLegacy = pkgs:
        pkgs.stdenv.mkDerivation (finalAttrs: {
          pname = "snapwave-legacy";
          version = "unstable-${self.shortRev or "dirty"}";

          src = legacySource pkgs;

          nativeBuildInputs = [ pkgs.gfortran ];
          buildInputs = [
            pkgs.netcdffortran
            pkgs.netcdf
            pkgs.hdf5
            pkgs.zlib
            pkgs.curl
          ];

          # Serial build: the Makefile relies on compilation order for the
          # Fortran module files (e.g. m_ec_triangle needs precision.mod),
          # which is not fully parallel-safe under make -j.
          enableParallelBuilding = false;

          makeFlags = [
            "FC=${pkgs.gfortran}/bin/gfortran"
            "CC=${pkgs.stdenv.cc}/bin/cc"
            # Absolute paths so the build does not rely on PATH lookups
            "NF_CONFIG=${pkgs.lib.getDev pkgs.netcdffortran}/bin/nf-config"
          ];

          installPhase = ''
            runHook preInstall
            install -Dm755 SnapWave/lnx64/bin/snapwave $out/bin/snapwave
            runHook postInstall
          '';

          meta = {
            description = "SnapWave legacy stand-alone Fortran build (numerical oracle)";
            mainProgram = "snapwave";
            license = nixpkgs.lib.licenses.lgpl21;
            platforms = nixpkgs.lib.platforms.unix;
          };
        });
    in
    {
      packages = forAllSystems (
        system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
        in
        {
          snapwave = mkSnapwave pkgs;
          snapwave-legacy = mkSnapwaveLegacy pkgs;
          default = mkSnapwave pkgs;
        }
      );

      checks = forAllSystems (
        system:
        let pkgs = nixpkgs.legacyPackages.${system};
        in
        {
          # Smoke test: run the "linear shoaling & refraction" coarse
          # testcase (2 timesteps, small grid) through the pure-Rust wrapper
          # and verify that it exits cleanly and produces valid NetCDF
          # map/history output.
          smoke-test = pkgs.stdenvNoCC.mkDerivation {
            name = "snapwave-smoke-test";
            dontUnpack = true;
            nativeBuildInputs = [ pkgs.netcdf ]; # ncdump

            buildCommand = ''
              cp -r "${self}/testcases/31_linear_shoaling_refraction" testcase
              chmod -R u+w testcase
              cd testcase/run/coarse

              # The checked-in testcase uses Windows path separators
              sed -i 's|\\|/|g' SnapWave.inp

              # No `mkdir -p ../../output` here: the wrapper owns
              # output-directory policy and creates missing directories
              # itself (plan.md Phase 5).

              ${self.packages.${system}.snapwave}/bin/snapwave --verbose SnapWave.inp > logfile.txt 2>&1

              test -f ../../output/shoalref_coarse_neu_map.nc
              test -f ../../output/shoalref_coarse_neu_his.nc
              ncdump -h ../../output/shoalref_coarse_neu_map.nc > /dev/null
              ncdump -h ../../output/shoalref_coarse_neu_his.nc > /dev/null

              mkdir -p $out
              cp logfile.txt $out/
            '';
          };
        }
      );

      devShells = forAllSystems (
        system:
        let pkgs = nixpkgs.legacyPackages.${system};
        in
        {
          default = pkgs.mkShell {
            packages =
              [
                pkgs.gfortran # Fortran compiler (Makefile oracle default FC)
                (pkgs.lib.getDev pkgs.netcdffortran) # nf-config, needed by the Makefile oracle
                pkgs.hdf5
                pkgs.zlib
                pkgs.curl
                pkgs.netcdf # ncdump & friends to inspect output files
                pkgs.gnumake
                pkgs.rustc # Rust toolchain for the pure-Rust wrapper
                pkgs.cargo
              ]
              ++ nixpkgs.lib.optionals
                pkgs.stdenv.hostPlatform.isLinux
                [ pkgs.gdb ];

            shellHook = ''
              echo "SnapWave dev shell"
              echo "  make                    release build (Fortran oracle) -> SnapWave/lnx64/bin/snapwave"
              echo "  make DEBUG=1            debug build (-g -O0 -fcheck=all -fbacktrace)"
              echo "  make clean              remove build artefacts"
              echo "  cargo build             pure-Rust wrapper (no Fortran/C/NetCDF)"
              echo "  cargo test              build + run tests (smoke, regression, CLI)"
              echo "  nix flake check         build + smoke-test the Rust wrapper"
              echo "  nix build .#snapwave-legacy  build the Fortran oracle via Nix"
            '';
          };
        }
      );
    };
}
