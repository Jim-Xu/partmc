#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

# fixed photolysis rates: the factor 0.5 must give the same result as
# rates that are half as large
../../partmc run_fixed_factor.spec
../../partmc run_fixed_half.spec
../../partmc run_fixed_full.spec
../../extract_gas out/fixed_factor_0001
../../extract_gas out/fixed_half_0001
../../extract_gas out/fixed_full_0001
../../numeric_diff --by col --rel-tol 1e-12 out/fixed_half_0001_gas.txt out/fixed_factor_0001_gas.txt

# and halving the rates must change the result
if ../../numeric_diff --by col --rel-tol 1e-6 out/fixed_full_0001_gas.txt out/fixed_half_0001_gas.txt ; then
    echo "halving the photolysis rates had no effect"
    exit 1
fi
