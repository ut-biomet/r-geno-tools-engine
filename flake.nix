{
  description = "Flake for a R environment";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
  }:
    flake-utils.lib.eachDefaultSystem (system:
        let
          pkgs = import nixpkgs {inherit system;};
          Rpkgs = pkgs;
          # # for specific R version:
          # Rpkgs =
          #   (import (pkgs.fetchFromGitHub {
          #       # Update to get other R version
          #       name = "nixpkgs-R-4.3.1";
          #       url = "https://github.com/NixOS/nixpkgs/";
          #       ref = "refs/heads/nixpkgs-unstable";
          #       rev = "976fa3369d722e76f37c77493d99829540d43845";
          #     })
          #     {inherit system;})
          #   .pkgs;
          R-with-my-packages = Rpkgs.rWrapper.override {
            packages = with Rpkgs.rPackages; [
              # list necessary R packages here
              gaston
              digest
              jsonlite
              manhattanly
              R6
              plotly
              networkD3
              ellipse
              /*
              developement packages
              */
              testthat
            ];
          };
        in {

          packages = {
          };

          devShells.default =
            pkgs.mkShell {
              LOCALE_ARCHIVE = if "${system}" == "x86_64-linux" then "${pkgs.glibcLocalesUtf8}/lib/locale/locale-archive" else "";
              R_LIBS_USER = "''"; # to no use users' installed R packages
              nativeBuildInputs = [pkgs.bashInteractive];
              buildInputs = [
                R-with-my-packages
              ];
            };

          apps = {
            # example of "apps" that could be run with `nix run .\#<name> -- --args`
            testsRepo = let
              testsRepo = pkgs.writeShellApplication {
                name = "testsRepo";
                text = ''
                  echo "simple unit tests:"
                  Rscript --vanilla ${./tests/testthat.R}
                '';
              };
            in {
              type = "app";
              program = "${testsRepo}/bin/testsRepo";
            };
          };
        }
        );
}

