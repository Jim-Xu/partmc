#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_gas out/mosaic_0001
# The gas phase depends on the aerosol surface area (reactive uptake of
# organic gases), and hence on the randomly sampled number of
# particles, so the tolerance is much larger than without uptake.
../../numeric_diff --by row --rel-tol 5e-3 ref_mosaic_0001_gas.txt out/mosaic_0001_gas.txt
