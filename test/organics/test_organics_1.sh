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
../../test_organics_process out/organics
../../numeric_diff --by elem --rel-tol 1e-10 ref_o_c.txt out/organics_o_c.txt
