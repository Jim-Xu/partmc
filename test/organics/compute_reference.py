#!/usr/bin/env python

# Compute the reference values of the organics test from the input
# files: for each of the three monodisperse populations the O:C, the
# OM:OC, the organic mass and the organic carbon mass of one particle.
# Writes ref_o_c.txt and ref_o_c_7col.txt.

import numpy as np

# species A, B, C, D
rho = np.array([1800.0, 1000.0, 1400.0, 1800.0])  # density (kg/m^3)
o_c = np.array([0.0, 0.2, 0.8, 0.0])              # O:C (1)
om_oc = np.array([0.0, 1.4, 2.0, 0.0])            # OM:OC (1), 0 = not organic

# mass fractions and diameters (m) of the populations
w = np.array([[0.0, 0.6, 0.4, 0.0],
              [0.5, 0.25, 0.0, 0.25],
              [1.0, 0.0, 0.0, 0.0]])
diam = np.array([1e-7, 2e-7, 4e-7])

def reference(om_oc, o_c):
    organic = om_oc > 0
    lines = []
    for (w_pop, d) in zip(w, diam):
        # mass of each species in one particle (kg)
        mass = np.pi / 6 * d**3 / np.sum(w_pop / rho) * w_pop
        org_mass = np.sum(mass[organic])
        carbon = mass[organic] / om_oc[organic]
        carbon_mass = np.sum(carbon)
        if carbon_mass > 0:
            pop_o_c = np.sum(carbon * o_c[organic]) / carbon_mass
            pop_om_oc = org_mass / carbon_mass
        else:
            pop_o_c = 0.0
            pop_om_oc = 0.0
        print("O:C = %.10f, OM:OC = %.10f, organic mass = %.10e kg, "
              "organic carbon mass = %.10e kg"
              % (pop_o_c, pop_om_oc, org_mass, carbon_mass))
        # smallest and largest value of each quantity are the same
        values = [pop_o_c, pop_o_c, pop_om_oc, pop_om_oc,
                  org_mass, org_mass, carbon_mass, carbon_mass]
        lines.append("".join("%30.20e" % v for v in values))
    return "\n".join(lines) + "\n"

print("species data with O:C and OM:OC")
open("ref_o_c.txt", "w").write(reference(om_oc, o_c))
print("species data without O:C and OM:OC: nothing is organic carbon")
open("ref_o_c_7col.txt", "w").write(reference(np.zeros(4), np.zeros(4)))
