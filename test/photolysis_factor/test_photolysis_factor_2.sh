#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

# clear-sky photolysis: at noon with the factor 0 the chemistry must be
# the same as at night
../../partmc run_solar_dark.spec
../../partmc run_solar_night.spec
../../partmc run_solar_noon.spec
../../extract_gas out/solar_dark_0001
../../extract_gas out/solar_night_0001
../../extract_gas out/solar_noon_0001
../../numeric_diff --by col --rel-tol 1e-12 out/solar_night_0001_gas.txt out/solar_dark_0001_gas.txt

# and the sun must matter without the factor
if ../../numeric_diff --by col --rel-tol 1e-6 out/solar_noon_0001_gas.txt out/solar_night_0001_gas.txt ; then
    echo "the noon run without the factor did not differ from the night run"
    exit 1
fi
