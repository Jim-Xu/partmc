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
../../extract_gas out/monoterpene_oh_0001

# API (column 95) against the solution with OH in steady state: the OH from
# the photolysis of H2O2 is consumed by API. The tolerance is 1% of the API
# that reacts (11.6 of 1000 ppb)
../../numeric_diff --by col --min-col 95 --max-col 95 --rel-tol 1e-4 ref_monoterpene_oh_gas.txt out/monoterpene_oh_0001_gas.txt
# H2O2 (column 17); the reference leaves out the H2O2 formed from HO2
../../numeric_diff --by col --min-col 17 --max-col 17 --rel-tol 1e-3 ref_monoterpene_oh_gas.txt out/monoterpene_oh_0001_gas.txt
