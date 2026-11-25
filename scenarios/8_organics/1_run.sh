#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

mkdir -p out

# The input files of the box case (organics.spec, the initial, background
# and emission files, the profiles) are in the repository. They are
# written by make_case.py: to run another case, edit the parameters there
# (or choose a variant: ./make_case.py forest) and rerun it before running.
# PartMC must be built with the MOSAIC library that has the CRACMM3
# organic aerosol.

../../build/partmc organics.spec

# Now run ./2_process.sh to process the data
