#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

mkdir -p out

# The emission files gas_emit.dat and aero_emit_comp_*.dat are in the
# repository. To change the emission scalings, edit make_emissions.py
# and regenerate them with ./make_emissions.py before running.

../../build/partmc organics.spec

# Now run ./2_process.sh to process the data
