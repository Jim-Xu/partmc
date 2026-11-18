#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

# without the iepox_acidity line MOSAIC must use the ISORROPIA II acidity
../../partmc run_default.spec
../../partmc run_isorropia.spec
../../extract_gas out/default_0001
../../extract_gas out/isorropia_0001
../../numeric_diff --by col --rel-tol 1e-12 out/isorropia_0001_gas.txt out/default_0001_gas.txt
