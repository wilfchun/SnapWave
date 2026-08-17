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

      # Keep only what is needed to build: Makefile, Fortran/C sources and the
      # bundled third_party code. Excludes Windows project files, docs,
      # testcases, the prebuilt zip and the Windows-only netcdf binaries.
      buildSource = pkgs:
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

      mkSnapwave = pkgs:
        pkgs.stdenv.mkDerivation (finalAttrs: {
          pname = "snapwave";
          version = "unstable-${self.shortRev or "dirty"}";

          src = buildSource pkgs;

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
            description = "Fast, implicit, unstructured-grid short wave solver";
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
          snapwave = mkSnapwave pkgs;
        in
        {
          inherit snapwave;
          default = snapwave;
        }
      );

      checks = forAllSystems (
        system:
        let pkgs = nixpkgs.legacyPackages.${system};
        in
        {
          # Smoke test: run the "linear shoaling & refraction" coarse testcase
          # (2 timesteps, small grid) and verify that the solver exits cleanly
          # and produces valid NetCDF map/history output.
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

              mkdir -p ../../output

              ${self.packages.${system}.snapwave}/bin/snapwave > logfile.txt

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
                pkgs.gfortran # Fortran compiler (Makefile default FC)
                (pkgs.lib.getDev pkgs.netcdffortran) # nf-config, needed by the Makefile
                pkgs.hdf5
                pkgs.zlib
                pkgs.curl
                pkgs.netcdf # ncdump & friends to inspect output files
                pkgs.gnumake
              ]
              ++ nixpkgs.lib.optionals
                pkgs.stdenv.hostPlatform.isLinux
                [ pkgs.gdb ];

            shellHook = ''
              echo "SnapWave dev shell"
              echo "  make             release build  -> SnapWave/lnx64/bin/snapwave"
              echo "  make DEBUG=1     debug build (-g -O0 -fcheck=all -fbacktrace)"
              echo "  make clean       remove build artefacts"
              echo "  nix flake check  build + run the testcase smoke test"
            '';
          };
        }
      );
    };
}
