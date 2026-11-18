#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

# MOSAIC's own acidity of these ammonium bisulfate particles is higher than
# the ISORROPIA II acidity, so the IEPOX uptake must differ
../../partmc run_mosaic.spec
../../extract_gas out/mosaic_0001
if ../../numeric_diff --by col --rel-tol 1e-3 out/isorropia_0001_gas.txt out/mosaic_0001_gas.txt ; then
    echo "the acidity option had no effect on the IEPOX uptake"
    exit 1
fi
