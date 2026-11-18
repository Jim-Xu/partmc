#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part_norl.spec
../../extract_gas out/norl_0001
../../extract_aero_time out/norl_0001

# gases: passive tracers, exact up to rounding
../../numeric_diff --by col --rel-tol 1e-10 ref_norl_gas.txt out/norl_0001_gas.txt
# particles: Monte Carlo sampling of 10000 particles
../../numeric_diff --by col --rel-tol 0.05 ref_norl_aero.txt out/norl_0001_aero_time.txt
