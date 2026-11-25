#!/usr/bin/env python3
"""Reference solution for the test of the photolysis of the monoterpene SOA.

Writes ref_mt_soa_photolysis_aero.txt and ref_mt_soa_photolysis_gas.txt in
the column layouts of extract_aero_time and extract_gas, with the analytic
values of AHOM, AMTN1 and HCOOH and zeros in all other columns. With the
fixed photolysis frequency J of NO2 (rates.dat) MOSAIC photolyzes the
highly oxygenated molecules of the monoterpene SOA in the particles at the
rate k = 0.01 J (SOAphot01 of CRACMM3, CMAQ 6.0):

    AHOM -> 0.8 AMTN1 + 0.2 HCOOH by moles,    y = 0.8 * 262 / 250 by mass,

so that with the initial mass concentration M0 of AHOM

    AHOM(t) = M0 exp(-k t),   AMTN1(t) = y M0 (1 - exp(-k t)),
    HCOOH(t) = 0.2 M0 (1 - exp(-k t)) / (0.250 kg/mol)  (mol m-3).

At 270 K both aerosol species stay in the particles (less than 0.1% in the
gas phase).
"""
import math

T, P, R = 270.0, 1e5, 8.314472        # K, Pa (temp.dat, pressure.dat), J mol-1 K-1
T_MAX, T_OUT = 3600.0, 600.0          # s (run_part.spec)
J_NO2 = 0.008                         # s-1 (rates.dat)
J_FACTOR, Y_MASS, Y_HCOOH = 0.01, 0.8 * 262.0 / 250.0, 0.2
NUM_CONC, DIAM = 1e9, 2e-7            # m-3, m (aero_init_dist.dat)
MASS_FRAC = dict(SO4=0.364, NH4=0.136, AHOM=0.5)     # aero_init_comp.dat

def main():
    gas = [l.split()[0] for l in open('gas_data.dat') if l.strip() and not l.startswith('#')]
    aero = [l.split() for l in open('aero_data.dat') if l.strip() and not l.startswith('#')]
    names = [a[0] for a in aero]
    dens = {a[0]: float(a[1]) for a in aero}
    mw = {a[0]: float(a[3].replace('d', 'e')) for a in aero}
    # mass of AHOM in a particle of the given volume and mass fractions
    vol = math.pi / 6.0 * DIAM**3
    m0 = NUM_CONC * vol * MASS_FRAC['AHOM'] / sum(f / dens[s] for s, f in MASS_FRAC.items())
    k = J_FACTOR * J_NO2
    air = P / (R * T)                                 # mol m-3
    fmt = lambda v: ''.join('%30.15e' % x for x in v) + '\n'
    with open('ref_mt_soa_photolysis_aero.txt', 'w') as fa, open('ref_mt_soa_photolysis_gas.txt', 'w') as fg:
        for i in range(int(round(T_MAX / T_OUT)) + 1):
            t = i * T_OUT
            gone = 1.0 - math.exp(-k * t)
            row = [t, 0.0, 0.0] + [0.0] * len(names)
            row[3 + names.index('AHOM')] = m0 * (1.0 - gone)
            row[3 + names.index('AMTN1')] = Y_MASS * m0 * gone
            fa.write(fmt(row))
            row = [t] + [0.0] * len(gas)
            row[1 + gas.index('HCOOH')] = Y_HCOOH * m0 * gone / mw['AHOM'] / air * 1e9
            fg.write(fmt(row))
    print('columns: AHOM %d, AMTN1 %d (aerosol); HCOOH %d (gas)' % (
        4 + names.index('AHOM'), 4 + names.index('AMTN1'), 2 + gas.index('HCOOH')))
    print('AHOM at the start %.4e kg m-3; photolyzed after %g s: %.2f%%' % (m0, T_MAX, 100 * (1 - math.exp(-k * T_MAX))))

if __name__ == '__main__':
    main()
