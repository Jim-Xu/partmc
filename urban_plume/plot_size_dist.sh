#!/bin/sh

echo Plotting size distribution
gnuplot -persist <<ENDPLOT
set logscale
set xlabel "radius (m)"
set ylabel "aerosol volume density (#/m^3)"
set title "URBAN_PLUME test case"
plot "out/urban_plume_summary_aero_binned.d" index 0 using 1:3 title "VL/SO4_a" with lines
set terminal postscript eps
set output "out/urban_plume_size_dist.eps"
replot
ENDPLOT
epstopdf out/urban_plume_size_dist.eps
