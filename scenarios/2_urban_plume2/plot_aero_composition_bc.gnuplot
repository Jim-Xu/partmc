# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xrange [0.01:1]
set yrange [0:80]

set logscale x

set xlabel "diameter / um"
set ylabel "BC dry mass fraction / %"

#    column  1: particle ID number
#    column  2: number concentration (m^{-3})
#    column  3: particle diameter (m)
#    column  4: particle total mass (kg)
#    column  5: SO4 mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column  6: NO3 mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column  7: Cl mass (kg) - density = 0.2200E+04 (kg/m^3)
#    column  8: NH4 mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column  9: MSA mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 10: AHOM mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 11: AELHOM mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 12: AOP3 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 13: AROCN2ALK mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 14: AROCN1ALK mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 15: AROCP0ALK mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 16: AROCP1ALK mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 17: AROCP2ALK mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 18: AROCP3ALK mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 19: AROCN2OXY2 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 20: AROCN2OXY4 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 21: AROCN2OXY8 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 22: AROCN1OXY1 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 23: AROCN1OXY3 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 24: AROCN1OXY6 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 25: AROCP0OXY2 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 26: AROCP0OXY4 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 27: AROCP1OXY1 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 28: AROCP1OXY3 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 29: AROCP2OXY2 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 30: AROCP3OXY2 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 31: ATRPN mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 32: AHONIT mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 33: AMTN1 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 34: CO3 mass (kg) - density = 0.2600E+04 (kg/m^3)
#    column 35: Na mass (kg) - density = 0.2200E+04 (kg/m^3)
#    column 36: Ca mass (kg) - density = 0.2600E+04 (kg/m^3)
#    column 37: OIN mass (kg) - density = 0.2600E+04 (kg/m^3)
#    column 38: OC mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 39: BC mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 40: AISO3NOS mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 41: AISO3OS mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 42: AORGC mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 43: ASOAT mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 44: AISO4 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 45: AISO5 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 46: AGLYOLIG mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 47: H2O mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 48: frozen flag

set multiplot layout 2,2

set title "BC composition at 0 hours"
plot "out/urban_plume2_0001_00000001_aero_particles.txt" using ($3*1e6):(($39/($4-$47))*100) with points notitle
set title "BC composition at 5 hours"
plot "out/urban_plume2_0001_00000006_aero_particles.txt" using ($3*1e6):(($39/($4-$47))*100) with points notitle
set title "BC composition at 7 hours"
plot "out/urban_plume2_0001_00000008_aero_particles.txt" using ($3*1e6):(($39/($4-$47))*100) with points notitle
set title "BC composition at 24 hours"
plot "out/urban_plume2_0001_00000025_aero_particles.txt" using ($3*1e6):(($39/($4-$47))*100) with points notitle

unset multiplot
