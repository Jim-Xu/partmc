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
../../extract_gas out/hpld_photolysis_0001

# HPLD (column 151) against the exponential decay at 0.055 times the
# photolysis frequency of NO2
../../numeric_diff --by col --min-col 151 --max-col 151 --rel-tol 1e-4 ref_hpld_photolysis_gas.txt out/hpld_photolysis_0001_gas.txt
