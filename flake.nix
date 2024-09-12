{
  description = "Flake for a R environment";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    # nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
  }:
    flake-utils.lib.eachDefaultSystem (
      system: let
        pkgs = import nixpkgs {inherit system;};

        vegan_withoutX = pkgs.rPackages.vegan.override {
          requireX = false;
        };
        vcfR_deps = with pkgs.rPackages; [
          ape
          dplyr
          magrittr
          memuse
          pinfsc50
          Rcpp
          stringr
          tibble
          viridisLite
          vegan_withoutX
        ];
        vcfR_withoutX = pkgs.rPackages.vcfR.overrideAttrs (oldAttrs: {
          propagatedBuildInputs = vcfR_deps;
          nativeBuildInputs = vcfR_deps;
        });

        R-packages = with pkgs.rPackages; [
          argparse
          R6
          digest
          ellipse
          gaston
          jsonlite
          manhattanly
          networkD3
          qqman
          plotly
          vcfR_withoutX
          pandoc
          Rd2md
          roxygen2

          sommer

          caret
          RAINBOWR

          (
            pkgs.rPackages.buildRPackage {
              name = "breedSimulatR";
              src = pkgs.fetchFromGitHub {
                owner = "ut-biomet";
                repo = "breedSimulatR";
                rev = "21fc8eb7f2e83f685a3eb99ca4bc611dee652ddd";
                sha256 = "sha256-JZvjTqlj4LK3tvLrrb2oVBgwTKU0Nntur6He2tQveCc=";
              };
              propagatedBuildInputs = with pkgs.rPackages; [data_table R6 vcfR_withoutX];
              nativeBuildInputs = with pkgs.rPackages; [data_table R6 vcfR_withoutX];
            }
          )

          /*
          developement packages
          */
          testthat
          AGHmatrix

          languageserver
          styler
        ];
      in rec {
        packages.r-geno-tools-engine = pkgs.callPackage ./nix_pkgs/default.nix {
          inherit pkgs;
        };
        packages.default = packages.r-geno-tools-engine;

        devShells.default = pkgs.mkShell {
          LOCALE_ARCHIVE =
            if "${system}" == "x86_64-linux"
            then "${pkgs.glibcLocalesUtf8}/lib/locale/locale-archive"
            else "";
          R_LIBS_USER = "''"; # to no use users' installed R packages
          R_PROFILE_USER = "''"; # to disable`.Rprofile` files (eg. when the project already use `renv`)
          nativeBuildInputs = [pkgs.bashInteractive];
          buildInputs = [
            (pkgs.rWrapper.override {packages = R-packages;})
            (pkgs.rstudioWrapper.override {packages = R-packages;})
            pkgs.pandoc
            pkgs.python3
            (
              (pkgs.callPackage ./nix_pkgs/default.nix {
                inherit pkgs;
              })
              .overrideAttrs (oldAttrs: {
                doCheck = false;
              })
            )
          ];
        };

        apps = {
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

          writeDoc = let
            writeDoc = pkgs.writeShellApplication {
              name = "writeDoc";
              text = ''
                Rscript --vanilla -e "source('tools/tools.R'); writeDoc()"
              '';
            };
          in {
            type = "app";
            program = "${writeDoc}/bin/writeDoc";
          };

          createResultExample = let
            createResultExample = pkgs.writeShellApplication {
              name = "createResultExample";
              text = ''
                Rscript --vanilla -e "source('tools/tools.R'); createResultExample()"
              '';
            };
          in {
            type = "app";
            program = "${createResultExample}/bin/createResultExample";
          };
        };
      }
    );
}
