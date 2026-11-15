#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

# The data should have already been processed by ./2_process.sh

./plot_organic_mass.py
./plot_o_c.py
./plot_diam_o_c_dist.py

# Now view the PDF files in the out/ directory
