#!/usr/bin/env python3
"""Reference solution for the HPLD photolysis test case.

Writes ref_hpld_photolysis_gas.txt in the column layout of extract_gas (time
and one column per gas species of gas_data.dat), with the mixing ratio of
HPLD and zeros in all other columns. HPLD, the hydroperoxy aldehyde from the
isomerization of the isoprene peroxy radical, photolyzes in MOSAIC at 0.055
times the photolysis frequency of NO2 (rates.dat) and has no other loss here
(its reaction with NO3 needs NOx):

    HPLD + hv -> OH + ISOPRD        0.055 J_NO2

so HPLD = HPLD_0 exp(-0.055 J_NO2 t).
"""
import math

T_MAX, T_OUT = 3600.0, 300.0          # s (run_part.spec)
HPLD_0 = 10.0                         # ppb (gas_init.dat)
J_NO2 = 5.0e-3                        # s-1 (rates.dat)
RATIO = 0.055                         # J_HPLD / J_NO2 in MOSAIC

def main():
    species = [l.split()[0] for l in open('gas_data.dat') if l.strip() and not l.startswith('#')]
    with open('ref_hpld_photolysis_gas.txt', 'w') as out:
        for i in range(int(round(T_MAX / T_OUT)) + 1):
            c = {'HPLD': HPLD_0 * math.exp(-RATIO * J_NO2 * i * T_OUT)}
            out.write(''.join('%30.15e' % x for x in [i * T_OUT] + [c.get(s, 0.0) for s in species]) + '\n')

if __name__ == '__main__':
    main()
