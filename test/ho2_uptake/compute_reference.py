#!/usr/bin/env python3
"""Reference solutions for the HO2 uptake test case.

Writes ref_ho2_uptake_<case>_gas.txt for the cases default, zero and one
(uptake coefficients 0.2, 0 and 1) in the column layout of extract_gas (time
and one column per gas species of gas_data.dat), with the mixing ratio of HO2
and zeros in all other columns. HO2 is lost on the particles at the frequency

    k = N pi D^2 / (D / (2 Dg) + 4 / (c gamma))

for N particles per volume of the diameter D, with the mean molecular speed
c = sqrt(8 R T / (pi M)) of HO2 and the gas-phase diffusivity
Dg = 8.55e-6 m^2/s (T / 273.15 K)^1.75 (101325 Pa / p) that CMAQ uses for its
uptake reactions, so HO2 = HO2_0 exp(-k t). The reaction of HO2 with itself
is negligible at the mixing ratio of the test (0.1% in 600 s).
"""
import math

T_MAX, T_OUT = 600.0, 60.0            # s (run_*.spec)
HO2_0 = 1.0e-5                        # ppb (gas_init.dat)
TEMP, PRES = 298.0, 1.0e5             # K, Pa (temp.dat, pressure.dat)
NUM_CONC, DIAM = 1.0e9, 2.0e-7        # m-3, m (aero_init_dist.dat)
GAMMA = {'default': 0.2, 'zero': 0.0, 'one': 1.0}
R_GAS, M_HO2 = 8.314472, 33.01e-3     # J/(mol K), kg/mol

def frequency(gamma):
    if gamma <= 0.0:
        return 0.0
    diffusivity = 0.0855e-4 * (TEMP / 273.15) ** 1.75 * 101325.0 / PRES
    speed = math.sqrt(8.0 * R_GAS * TEMP / (math.pi * M_HO2))
    return NUM_CONC * math.pi * DIAM ** 2 / (DIAM / (2.0 * diffusivity) + 4.0 / (speed * gamma))

def main():
    species = [l.split()[0] for l in open('gas_data.dat') if l.strip() and not l.startswith('#')]
    for case, gamma in GAMMA.items():
        k = frequency(gamma)
        print('%-8s gamma %.1f  k = %.4e s-1' % (case, gamma, k))
        with open('ref_ho2_uptake_%s_gas.txt' % case, 'w') as out:
            for i in range(int(round(T_MAX / T_OUT)) + 1):
                c = {'HO2': HO2_0 * math.exp(-k * i * T_OUT)}
                out.write(''.join('%30.15e' % x for x in [i * T_OUT] + [c.get(s, 0.0) for s in species]) + '\n')

if __name__ == '__main__':
    main()
