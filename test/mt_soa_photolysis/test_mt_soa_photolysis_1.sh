#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part.spec
../../extract_aero_time out/mt_soa_photolysis_0001
../../extract_gas out/mt_soa_photolysis_0001

# AHOM (column 9) and its photolysis product AMTN1 (column 32) against the
# analytic solution for 1% of the fixed photolysis frequency of NO2; the
# tolerance covers the Monte Carlo sampling of 5000 particles
../../numeric_diff --by col --min-col 9 --max-col 9 --rel-tol 0.05 ref_mt_soa_photolysis_aero.txt out/mt_soa_photolysis_0001_aero_time.txt
../../numeric_diff --by col --min-col 32 --max-col 32 --rel-tol 0.05 ref_mt_soa_photolysis_aero.txt out/mt_soa_photolysis_0001_aero_time.txt
# formic acid (column 30) released to the gas phase
../../numeric_diff --by col --min-col 30 --max-col 30 --rel-tol 0.05 ref_mt_soa_photolysis_gas.txt out/mt_soa_photolysis_0001_gas.txt
