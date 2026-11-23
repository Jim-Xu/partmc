#!/usr/bin/env python3
"""Reference solution for the monoterpene + OH test case.

Writes ref_monoterpene_oh_gas.txt in the column layout of extract_gas (time
and one column per gas species of gas_data.dat), with the mixing ratios of
API and H2O2 and zeros in all other columns. OH comes from the photolysis of
H2O2 at the fixed frequency J (rates.dat) and reacts with API or with H2O2:

    H2O2 + hv -> 2 OH            J
    API + OH  -> products        kA = 12.1e-12 exp(444 / T)
    H2O2 + OH -> HO2 + H2O       kH = 2.9e-12 exp(-160 / T)

(rate constants of MOSAIC, cm3 molecule-1 s-1). OH is in steady state, so with
the fraction f = kA [API] / (kA [API] + kH [H2O2]) of the OH that reacts with
API:

    d[API]/dt  = -2 J [H2O2] f
    d[H2O2]/dt = -J [H2O2] - 2 J [H2O2] (1 - f)

Left out: the reactions of OH with the products of API, which take 0.3% of the
OH at the end of the run, and the H2O2 formed from HO2 (0.25 ppb in the run).
The two equations are integrated with the classical Runge-Kutta method.
"""
import math

T, P = 298.0, 1e5                     # K, Pa (temp.dat, pressure.dat)
T_MAX, T_OUT = 600.0, 60.0            # s (run_part.spec)
H2O2_0, API_0 = 1000.0, 1000.0        # ppb (gas_init.dat)
J = 1.0e-5                            # s-1 (rates.dat)
K_BOLTZMANN = 1.380649e-23            # J K-1

def main():
    species = [l.split()[0] for l in open('gas_data.dat') if l.strip() and not l.startswith('#')]
    ka = 12.1e-12 * math.exp(444.0 / T)
    kh = 2.9e-12 * math.exp(-160.0 / T)
    def rhs(y):
        api, h2o2 = y
        f = ka * api / (ka * api + kh * h2o2)      # the ratio does not depend on the unit
        return (-2.0 * J * h2o2 * f, -J * h2o2 - 2.0 * J * h2o2 * (1.0 - f))
    n_sub = 600; h = T_OUT / n_sub
    y = (API_0, H2O2_0)
    with open('ref_monoterpene_oh_gas.txt', 'w') as out:
        for i in range(int(round(T_MAX / T_OUT)) + 1):
            c = {'API': y[0], 'H2O2': y[1]}
            out.write(''.join('%30.15e' % x for x in [i * T_OUT] + [c.get(s, 0.0) for s in species]) + '\n')
            for _ in range(n_sub):
                k1 = rhs(y)
                k2 = rhs(tuple(a + 0.5 * h * b for a, b in zip(y, k1)))
                k3 = rhs(tuple(a + 0.5 * h * b for a, b in zip(y, k2)))
                k4 = rhs(tuple(a + h * b for a, b in zip(y, k3)))
                y = tuple(a + h / 6.0 * (b + 2 * c2 + 2 * c3 + d) for a, b, c2, c3, d in zip(y, k1, k2, k3, k4))

if __name__ == '__main__':
    main()
