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
../../extract_gas out/beta_pinene_0001

# BPIN (column 149) against the analytic solution of BPIN + O3 with the IUPAC
# rate constant
../../numeric_diff --by col --min-col 149 --max-col 149 --rel-tol 1e-4 ref_beta_pinene_gas.txt out/beta_pinene_0001_gas.txt
# O3 (column 12): the HO2 formed from the OH of the reaction consumes some more
# O3 than the analytic solution (0.06 ppb of 100 ppb)
../../numeric_diff --by col --min-col 12 --max-col 12 --rel-tol 1e-3 ref_beta_pinene_gas.txt out/beta_pinene_0001_gas.txt
