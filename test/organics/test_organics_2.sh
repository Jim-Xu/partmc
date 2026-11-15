#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

# species data file without the O:C and OM:OC columns
../../partmc run_part_7col.spec
../../test_organics_process out/organics_7col
../../numeric_diff --by elem --abs-tol 1e-30 ref_o_c_7col.txt out/organics_7col_o_c.txt
