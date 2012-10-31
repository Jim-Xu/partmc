
# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set logscale y
set xlabel "dimensionless time"
set ylabel "normalized number concentration"

set xrange [0.1:1000]
set yrange [1e-6:1]

plot "out/part_brown_free_df_2_4_0001_dimless_t_series.txt" using 1:2 title "Df = 2.4, PartMC", \
     "ref_free_df_2_4_dimless_time_regrid.txt" using 1:2 with lines title "Df = 2.4, Ref"
