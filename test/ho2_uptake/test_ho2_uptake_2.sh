#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_zero.spec
../../extract_gas out/ho2_uptake_zero_0001

../../numeric_diff --by col --min-col 16 --max-col 16 --rel-tol 0.05 ref_ho2_uptake_zero_gas.txt out/ho2_uptake_zero_0001_gas.txt
