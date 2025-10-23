{
  lib,
  pkgs,
  runCommand,
  ...
}:
let
  # To undo the effect of `requireX = true`
  vegan_withoutX = pkgs.rPackages.vegan.override { requireX = false; };
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

  engines_R_deps = with pkgs.rPackages; [
    argparse
    R6
    digest
    ellipse
    gaston
    jsonlite
    manhattanly
    networkD3
    plotly
    qqman
    vcfR_withoutX
    sommer
    caret
    RAINBOWR

    (pkgs.rPackages.buildRPackage {
      name = "breedSimulatR";
      src = pkgs.fetchFromGitHub {
        owner = "ut-biomet";
        repo = "breedSimulatR";
        rev = "d5f485f2c8154a4a0af661e176fbaf09376f4c01";
        sha256 = "sha256-PBgDW3ZopeDkaMXDyT+CeGEixEuhGJrSGb2eLOJ1O2U=";
      };
      propagatedBuildInputs = with pkgs.rPackages; [
        data_table
        R6
        vcfR_withoutX
      ];
      nativeBuildInputs = with pkgs.rPackages; [
        data_table
        R6
        vcfR_withoutX
      ];
    })
  ];

  tests_R_deps = with pkgs.rPackages; [
    testthat
    AGHmatrix
  ];

  R_with_packages = pkgs.rWrapper.override { packages = engines_R_deps; };
  R_with_packages_test = pkgs.rWrapper.override { packages = engines_R_deps ++ tests_R_deps; };
in
pkgs.stdenv.mkDerivation (finalAttrs: rec {
  pname = "r-geno-tools-engine";
  version = "v1.3.0";

  src = pkgs.lib.sources.cleanSource ../.;

  # To skip the Makefile `make deps`
  # dontBuild = true;
  doCheck = true;

  propagatedBuildInputs = [
    R_with_packages
    pkgs.python3
    pkgs.pandoc
    pkgs.fontconfig
  ];

  nativeBuildInputs = [
    R_with_packages_test
    pkgs.makeWrapper
  ];

  patchPhase = ''
    runHook prePatch

    substituteInPlace ./r-geno-tools-engine.R --replace-fail '#!/usr/bin/env Rscript' '#!${R_with_packages}/bin/Rscript --vanilla'

    runHook postPatch
  '';

  buildPhase = ''

    mkdir -p build/src
    cp -r ./src/* build/src/.

    mkdir -p build/bin
    cp ./r-geno-tools-engine.R build/bin/r-geno-tools-engine.R
    wrapProgram $(pwd)/build/bin/r-geno-tools-engine.R \
      --set PATH ${lib.makeBinPath (propagatedBuildInputs ++ [ pkgs.coreutils ])} \
      --set FONTCONFIG_FILE ${pkgs.fontconfig.out}/etc/fonts/fonts.conf \
      --set FONTCONFIG_PATH ${pkgs.fontconfig.out}/etc/fonts/ \
      --set RGENOROOT $(pwd)/build \
      --set R_LIBS_USER "\"\""

    mkdir build/test_data
    cp ./data/geno/testMarkerData01.vcf.gz build/test_data/.
    cp ./data/pheno/testPhenoData01.csv build/test_data/.
    cp ./data/pedigree/testPedData_char.csv build/test_data/.
    cp ./data/geno/breedGame_geno.vcf.gz build/test_data/.
    cp ./data/geno/breedGame_phasedGeno.vcf.gz build/test_data/.
    cp ./data/pedigree/breedGame_pedigree.csv build/test_data/.
    cp ./data/crossingTable/breedGame_small_crossTable.csv build/test_data/.
    cp ./data/SNPcoordinates/breedingGame_SNPcoord.csv build/test_data/.
    cp ./data/markerEffects/breedGame_markerEffects.csv build/test_data/.
    cp ./data/markerEffects/breedGame_markerEffects_2traits.json build/test_data/.
    cp ./data/genomic_selection/geno_G1.vcf.gz build/test_data/.
    cp ./data/genomic_selection/geno_G2.vcf.gz build/test_data/.
    cp ./data/genomic_selection/pheno_test.csv build/test_data/.
    cp ./data/genomic_selection/pheno_train.csv build/test_data/.
  '';

  checkPhase = ''
    export XDG_CACHE_HOME="$(mktemp -d)"
    export FONTCONFIG_FILE=${pkgs.fontconfig.out}/etc/fonts/fonts.conf
    export FONTCONFIG_PATH=${pkgs.fontconfig.out}/etc/fonts/

    RGENOROOT="build"
    ROOT_DATA_FILES="."
    Rscript --vanilla ./tests/testthat.R
  '';

  installPhase = ''
    runHook preInstall

    substituteInPlace ./build/bin/r-geno-tools-engine.R --replace-fail $(pwd)/build $out
    cp -r ./build $out
    mv $out/bin/r-geno-tools-engine.R $out/bin/r-geno-tools-engine

    runHook postInstall
  '';

  passthru.tests =
    runCommand "r-geno-tools-engine-test" { nativeBuildInputs = [ finalAttrs.finalPackage ]; }
      ''
        # Necessary for the derivation to be successful
        mkdir $out
        mkdir $out/logs
        r-geno-tools-engine > $out/logs/help.txt

        ROOT_DATA_FILES=${finalAttrs.finalPackage}/test_data

        echo "Test: r-geno-tools-engine gwas..."
        r-geno-tools-engine gwas \
          --genoFile "$ROOT_DATA_FILES/testMarkerData01.vcf.gz" \
          --phenoFile "$ROOT_DATA_FILES/testPhenoData01.csv" \
          --trait "Flowering.time.at.Arkansas" \
          --test "score" \
          --response "quantitative" \
          --thresh-maf 0.05 \
          --thresh-callrate 0.95 \
          --outFile "$out/gwasRes_score.json" | tee $out/logs/gwas.log

        echo "Test: r-geno-tools-engine gwas-manplot..."
        r-geno-tools-engine gwas-manplot \
          --gwasFile "$out/gwasRes_score.json" \
          --adj-method "bonferroni" \
          --thresh-p 0.05 \
          --filter-nPoints 3000 \
          --outFile "$out/manPlot.html"  | tee $out/logs/manPlot_html.log

        echo "Test: r-geno-tools-engine gwas-manplot..."
        r-geno-tools-engine gwas-manplot \
          --gwasFile "$out/gwasRes_score.json" \
          --adj-method "bonferroni" \
          --thresh-p 0.05 \
          --filter-nPoints 3000 \
          --no-interactive \
          --outFile "$out/manPlot.png" | tee $out/logs/manPlot_html.log

        echo "Test: r-geno-tools-engine gwas-adjresults..."
        r-geno-tools-engine gwas-adjresults \
          --gwasFile "$out/gwasRes_score.json" \
          --adj-method "bonferroni" \
          --filter-nPoints 3000 \
          --outFile "$out/gwasRes_adj.json" | tee $out/logs/gwasRes_adj.log

        echo "Test: r-geno-tools-engine ldplot..."
        r-geno-tools-engine ldplot \
          --genoFile "$ROOT_DATA_FILES/testMarkerData01.vcf.gz" \
          --from 42 \
          --to 62 \
          --outFile "$out/ldplot.png" | tee $out/logs/ldplot.log

        echo "Test: r-geno-tools-engine relmat-ped..."
        r-geno-tools-engine relmat-ped \
          --pedFile "$ROOT_DATA_FILES/breedGame_pedigree.csv" \
          --outFile "$out/pedRelMat.json" | tee $out/logs/pedRelMat.log

        echo "Test: r-geno-tools-engine relmat-geno..."
        r-geno-tools-engine relmat-geno \
          --genoFile "$ROOT_DATA_FILES/breedGame_geno.vcf.gz" \
          --outFile "$out/genoRelMat.json" | tee $out/logs/genoRelMat.log

        echo "Test: r-geno-tools-engine relmat-combined..."
        r-geno-tools-engine relmat-combined \
          --ped-relmatFile $out/pedRelMat.json \
          --geno-relmatFile $out/genoRelMat.json \
          --combine-method Martini \
          --tau 1.3 \
          --omega 0.7 \
          --outFile "$out/combinedRelMat.json" | tee $out/logs/combinedRelMat.log

        echo "Test: r-geno-tools-engine relmat-heatmap..."
        r-geno-tools-engine relmat-heatmap \
          --relmatFile "$out/genoRelMat.json" \
          --no-interactive \
          --outFile "$out/relMat_heatmap.png" | tee $out/logs/relMat_heatmap.log

        echo "Test: r-geno-tools-engine pedNetwork..."
        r-geno-tools-engine pedNetwork \
          --pedFile "$ROOT_DATA_FILES/breedGame_pedigree.csv" \
          --outFile "$out/pedNet.html" | tee $out/logs/pedNet.log

        echo "Test: r-geno-tools-engine crossing-simulation..."
        r-geno-tools-engine crossing-simulation \
          --genoFile "$ROOT_DATA_FILES/breedGame_phasedGeno.vcf.gz" \
          --crossTableFile "$ROOT_DATA_FILES/breedGame_small_crossTable.csv" \
          --SNPcoordFile "$ROOT_DATA_FILES/breedingGame_SNPcoord.csv" \
          --outFile "$out/crossSim.vcf.gz" | tee $out/logs/crossSim.log

        echo "Test: r-geno-tools-engine progeny-blup-calculation 1..."
        r-geno-tools-engine progeny-blup-calculation \
          --genoFile "$ROOT_DATA_FILES/breedGame_phasedGeno.vcf.gz" \
          --crossTableFile "$ROOT_DATA_FILES/breedGame_small_crossTable.csv" \
          --SNPcoordFile "$ROOT_DATA_FILES/breedingGame_SNPcoord.csv" \
          --markerEffectsFile "$ROOT_DATA_FILES/breedGame_markerEffects.csv" \
          --outFile "$out/progBlups.json" | tee $out/logs/progBlups.log

        echo "Test: r-geno-tools-engine progeny-blup-calculation 2..."
        r-geno-tools-engine progeny-blup-calculation \
          --genoFile "$ROOT_DATA_FILES/breedGame_phasedGeno.vcf.gz" \
          --crossTableFile "$ROOT_DATA_FILES/breedGame_small_crossTable.csv" \
          --SNPcoordFile "$ROOT_DATA_FILES/breedingGame_SNPcoord.csv" \
          --markerEffectsFile "$ROOT_DATA_FILES/breedGame_markerEffects_2traits.json" \
          --outFile "$out/progBlups_2traits.json" | tee $out/logs/progBlups_2traits.log

        echo "Test: r-geno-tools-engine progeny-blup-plot (with 1 trait)..."
        r-geno-tools-engine progeny-blup-plot \
          --progeniesBlupFile "$out/progBlups.json" \
          --y-axis-name "Phenotypic trait" \
          --error-bar-interval 0.95 \
          --outFile "$out/blupPlot_1trait.html" | tee $out/logs/blupPlot_1trait.log

        echo "Test: r-geno-tools-engine progeny-blup-plot (with 2 traits)..."
        r-geno-tools-engine progeny-blup-plot \
          --progeniesBlupFile "$out/progBlups_2traits.json" \
          --y-axis-name "Phenotypic trait" \
          --error-bar-interval 0.95 \
          --trait "trait1" \
          --outFile "$out/blupPlot_1trait_with_multi-trait-input.html" | tee $out/logs/blupPlot_1trait_with_multi-trait-input.log

        echo "Test: r-geno-tools-engine progeny-blup-plot-2-traits..."
        r-geno-tools-engine progeny-blup-plot-2-traits \
          --progeniesBlupFile "$out/progBlups_2traits.json" \
          --x-trait "trait1" \
          --y-trait "trait2" \
          --confidence-level 0.95 \
          --outFile "$out/blupPlot_2traits.html" | tee $out/logs/blupPlot_2traits.log

        echo "Test: DONE"
      '';
})
