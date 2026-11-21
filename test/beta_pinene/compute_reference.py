#!/usr/bin/env python3
"""Reference solution for the beta-pinene test case.

Writes ref_beta_pinene_gas.txt in the column layout of extract_gas (time and
one column per gas species of gas_data.dat), with the analytic mixing ratios
of O3 and BPIN and zeros in all other columns. In the dark and without other
reactive gases the only reaction is

    BPIN + O3 -> products,    k = 1.39e-15 exp(-1280 / T) cm3 molecule-1 s-1

(IUPAC preferred value, datasheet Ox_VOC19; Cox et al., 2020, Atmos. Chem.
Phys. 20, 13497), which consumes one O3 per BPIN, so that with the constant
difference D = [O3] - [BPIN] and E = exp(-k D t):

    [BPIN](t) = [BPIN]_0 D E / ([O3]_0 - [BPIN]_0 E),   [O3] = [BPIN] + D.
"""
import math

T, P = 298.0, 1e5                     # K, Pa (temp.dat, pressure.dat)
T_MAX, T_OUT = 3600.0, 600.0          # s (run_part.spec)
O3_0, BPIN_0 = 100.0, 1.0             # ppb (gas_init.dat)
K_BOLTZMANN = 1.380649e-23            # J K-1

def main():
    species = [l.split()[0] for l in open('gas_data.dat') if l.strip() and not l.startswith('#')]
    air = P / (K_BOLTZMANN * T) * 1e-6                # molecules cm-3
    k = 1.39e-15 * math.exp(-1280.0 / T)              # cm3 molecule-1 s-1
    d = (O3_0 - BPIN_0) * 1e-9 * air                  # molecules cm-3
    with open('ref_beta_pinene_gas.txt', 'w') as f:
        for i in range(int(round(T_MAX / T_OUT)) + 1):
            t = i * T_OUT
            e = math.exp(-k * d * t)
            bpin = BPIN_0 * (O3_0 - BPIN_0) * e / (O3_0 - BPIN_0 * e)
            c = {'BPIN': bpin, 'O3': bpin + O3_0 - BPIN_0}
            f.write(''.join('%30.15e' % x for x in [t] + [c.get(s, 0.0) for s in species]) + '\n')

if __name__ == '__main__':
    main()
