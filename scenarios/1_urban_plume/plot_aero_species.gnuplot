# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set key top left

set title "Aerosol bulk mass concentrations (with coag)"

set xrange [0:24]
set xtics 3

set xlabel "time / h"
set ylabel "mass concentration / (ug/m^3)"

#    column  1: time (s)
#    column  2: aerosol number concentration (#/m^3)
#    column  3: aerosol mass concentration (kg/m^3)
#    column  4: aerosol SO4 concentration (kg/m^3)
#    column  5: aerosol NO3 concentration (kg/m^3)
#    column  6: aerosol Cl concentration (kg/m^3)
#    column  7: aerosol NH4 concentration (kg/m^3)
#    column  8: aerosol MSA concentration (kg/m^3)
#    column  9: aerosol AHOM concentration (kg/m^3)
#    column 10: aerosol AELHOM concentration (kg/m^3)
#    column 11: aerosol AOP3 concentration (kg/m^3)
#    column 12: aerosol AROCN2ALK concentration (kg/m^3)
#    column 13: aerosol AROCN1ALK concentration (kg/m^3)
#    column 14: aerosol AROCP0ALK concentration (kg/m^3)
#    column 15: aerosol AROCP1ALK concentration (kg/m^3)
#    column 16: aerosol AROCP2ALK concentration (kg/m^3)
#    column 17: aerosol AROCP3ALK concentration (kg/m^3)
#    column 18: aerosol AROCN2OXY2 concentration (kg/m^3)
#    column 19: aerosol AROCN2OXY4 concentration (kg/m^3)
#    column 20: aerosol AROCN2OXY8 concentration (kg/m^3)
#    column 21: aerosol AROCN1OXY1 concentration (kg/m^3)
#    column 22: aerosol AROCN1OXY3 concentration (kg/m^3)
#    column 23: aerosol AROCN1OXY6 concentration (kg/m^3)
#    column 24: aerosol AROCP0OXY2 concentration (kg/m^3)
#    column 25: aerosol AROCP0OXY4 concentration (kg/m^3)
#    column 26: aerosol AROCP1OXY1 concentration (kg/m^3)
#    column 27: aerosol AROCP1OXY3 concentration (kg/m^3)
#    column 28: aerosol AROCP2OXY2 concentration (kg/m^3)
#    column 29: aerosol AROCP3OXY2 concentration (kg/m^3)
#    column 30: aerosol ATRPN concentration (kg/m^3)
#    column 31: aerosol AHONIT concentration (kg/m^3)
#    column 32: aerosol AMTN1 concentration (kg/m^3)
#    column 33: aerosol CO3 concentration (kg/m^3)
#    column 34: aerosol Na concentration (kg/m^3)
#    column 35: aerosol Ca concentration (kg/m^3)
#    column 36: aerosol OIN concentration (kg/m^3)
#    column 37: aerosol OC concentration (kg/m^3)
#    column 38: aerosol BC concentration (kg/m^3)
#    column 39: aerosol AISO3NOS concentration (kg/m^3)
#    column 40: aerosol AISO3OS concentration (kg/m^3)
#    column 41: aerosol AORGC concentration (kg/m^3)
#    column 42: aerosol ASOAT concentration (kg/m^3)
#    column 43: aerosol AISO4 concentration (kg/m^3)
#    column 44: aerosol AISO5 concentration (kg/m^3)
#    column 45: aerosol AGLYOLIG concentration (kg/m^3)
#    column 46: aerosol H2O concentration (kg/m^3)

set ytics nomirror
set y2tics

set multiplot layout 2,1

plot "out/urban_plume_0001_aero_time.txt" using ($1/3600):($5*1e9) axes x1y1 with lines title "NO3", \
     "out/urban_plume_0001_aero_time.txt" using ($1/3600):($37*1e9) axes x1y1 with lines title "OC", \
     "out/urban_plume_0001_aero_time.txt" using ($1/3600):($7*1e9) axes x1y1 with lines title "NH4"

plot "out/urban_plume_0001_aero_time.txt" using ($1/3600):($4*1e9) axes x1y1 with lines title "SO4", \
     "out/urban_plume_0001_aero_time.txt" using ($1/3600):($38*1e9) axes x1y1 with lines title "BC", \
     "out/urban_plume_0001_aero_time.txt" using ($1/3600):((sum [i=9:32] column(i) + sum [i=39:45] column(i))*1e9) axes x1y1 with lines title "SOA" # AHOM to AMTN1 and AISO3NOS to AGLYOLIG (31 species)

unset multiplot
