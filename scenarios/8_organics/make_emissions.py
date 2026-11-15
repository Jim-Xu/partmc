#!/usr/bin/env python

# Generate the emission files of the 8_organics scenario (urban plume
# case):
#
#   gas_emit.dat              gas emissions, from gas_emit_base.dat plus
#                             the precursors of the organic aerosol scheme
#   aero_emit_comp_*.dat      composition of the emitted combustion
#                             particles, with semivolatile primary
#                             organic matter instead of OC
#
# Usage: ./make_emissions.py
#
# All scaling factors are collected below, see 0_README.md for where
# they come from. Edit them and rerun this script to change the
# emissions.

import math, re

######################################################################
# gas-phase precursors: (new species, reference species, molar ratio)

GAS_SCALING = [
    # shadow precursors of the aromatic SOA: same emissions as the
    # CBM-Z species (SOAXYL without the naphthalene part of XYL)
    ("SOATOL", "TOL", 1.0),
    ("SOAXYL", "XYL", 0.998),
    # naphthalene and long-chain alkanes from the lumped species of an
    # inventory that lacks them (CMAQv5.1 VOC emission updates)
    ("NAPH", "XYL", 0.002),
    ("SOAALK", "PAR", 0.108),
    # benzene: urban emission ratio to CO in Los Angeles, 1.30 pptv/ppbv
    # (Borbon et al., 2013, Table 2)
    ("BENZ", "CO", 1.3e-3),
    # terpenes relative to isoprene, from the global annual emissions
    # of MEGAN2.1 (Guenther et al., 2012, Table 6):
    # isoprene 535 Tg, alpha-pinene 66.1 Tg, all other monoterpenes
    # 96.2 Tg, sesquiterpenes 29.0 Tg
    ("API", "ISOP", (66.1 / 136.24) / (535.0 / 68.12)),
    ("LIM", "ISOP", (96.2 / 136.24) / (535.0 / 68.12)),
    ("SESQ", "ISOP", (29.0 / 204.36) / (535.0 / 68.12)),
]

######################################################################
# primary organic matter (POA)

# mass fractions of the volatility bins of the emitted POA
# (CMAQ v5.5, CMAQ_Control_DESID_cb6r5_ae7_aq.nml; Murphy et al., 2017)
POA_SPLIT = [
    ("ALVPO1", 0.09),
    ("ASVPO1", 0.09),
    ("ASVPO2", 0.14),
    ("ASVPO3", 0.18),
    ("AIVPO1", 0.50),
]

# mass of the precursor of potential combustion SOA (PCVOC) emitted per
# mass of POA (same sources), and its molecular weight (kg/mol)
PCVOC_PER_POA = 6.579
PCVOC_MOLEC_WEIGHT = 170e-3

# composition of the emitted combustion particles (mass fractions), as
# in the 1_urban_plume scenario; "POA" is split according to POA_SPLIT
AERO_EMIT_COMP = {
    "aero_emit_comp_cooking.dat": [("POA", 1.0)],
    "aero_emit_comp_diesel.dat": [("POA", 0.3), ("BC", 0.7)],
    "aero_emit_comp_gasol.dat": [("POA", 0.8), ("BC", 0.2)],
}

######################################################################

def read_data_lines(filename):
    """Return the non-comment lines of a file as lists of words."""
    lines = []
    for line in open(filename):
        words = line.split("#")[0].split()
        if words:
            lines.append(words)
    return lines

def read_densities(filename):
    """Return {species: density (kg/m^3)} from an aero_data file."""
    return dict((w[0], float(w[1])) for w in read_data_lines(filename))

def read_modes(filename):
    """Return the log-normal modes of an aero_dist file as dicts."""
    modes = []
    for words in read_data_lines(filename):
        if words[0] == "mode_name":
            modes.append({})
        modes[-1][words[0]] = words[1]
    return modes

def expand_comp(comp):
    """Replace the entry "POA" of a composition by the volatility bins."""
    species = []
    for (name, frac) in comp:
        if name == "POA":
            species.extend((s, frac * f) for (s, f) in POA_SPLIT)
        else:
            species.append((name, frac))
    return species

def poa_mass_flux(dist_filename, densities):
    """POA mass emitted per area and time (kg m^{-2} s^{-1}) by all
    modes of the emission distribution, for a unit rate scaling."""
    flux = 0.0
    for mode in read_modes(dist_filename):
        if mode["mass_frac"] not in AERO_EMIT_COMP:
            continue
        assert mode["mode_type"] == "log_normal"
        assert mode["diam_type"] == "geometric"
        species = expand_comp(AERO_EMIT_COMP[mode["mass_frac"]])
        total = sum(f for (s, f) in species)
        density = total / sum(f / densities[s] for (s, f) in species)
        poa_frac = sum(f for (s, f) in species
                       if s in dict(POA_SPLIT)) / total
        ln_sigma = float(mode["log10_geom_std_dev"]) * math.log(10.0)
        mean_vol = math.pi / 6.0 * float(mode["geom_mean_diam"])**3 \
            * math.exp(4.5 * ln_sigma**2)
        mode_flux = float(mode["num_conc"]) * mean_vol * density * poa_frac
        print("  mode %-10s POA mass flux %.3e kg m^{-2} s^{-1}"
              % (mode["mode_name"], mode_flux))
        flux += mode_flux
    return flux

def rate_at(times, rates, t):
    """Rate scaling in effect at time t (PartMC uses the entry of the
    last time that is not after t)."""
    rate = rates[0]
    for (time, r) in zip(times, rates):
        if time <= t:
            rate = r
    return rate

def write_aero_comps():
    for (filename, comp) in sorted(AERO_EMIT_COMP.items()):
        with open(filename, "w") as f:
            f.write("# generated by make_emissions.py\n")
            f.write("#       mass fractions\n")
            for (name, frac) in expand_comp(comp):
                f.write("%-14s %g\n" % (name, frac))
        print("Wrote %s" % filename)

def write_gas_emit():
    base = read_data_lines("gas_emit_base.dat")
    assert base[0][0] == "time" and base[1][0] == "rate"
    times = [float(x) for x in base[0][1:]]
    rates = [float(x) for x in base[1][1:]]
    emit = dict((w[0], [float(x) for x in w[1:]]) for w in base[2:])
    order = [w[0] for w in base[2:]]
    for (name, values) in emit.items():
        assert len(values) == len(times), name

    for (new, ref, ratio) in GAS_SCALING:
        assert new not in emit
        emit[new] = [ratio * x for x in emit[ref]]
        order.append(new)

    # PCVOC follows the POA of the particle emissions, which have their
    # own times and rate scalings
    aero = dict((w[0], w[1:]) for w in read_data_lines("aero_emit.dat"))
    assert len(set(aero["dist"])) == 1
    aero_times = [float(x) for x in aero["time"]]
    aero_rates = [float(x) for x in aero["rate"]]
    print("POA emissions of %s:" % aero["dist"][0])
    flux = poa_mass_flux(aero["dist"][0], read_densities("aero_data.dat"))
    print("  total      POA mass flux %.3e kg m^{-2} s^{-1}" % flux)
    emit["PCVOC"] = []
    for (time, rate) in zip(times, rates):
        aero_rate = rate_at(aero_times, aero_rates, time)
        if aero_rate > 0.0:
            # the value is multiplied by the gas rate scaling in PartMC
            assert rate > 0.0
            emit["PCVOC"].append(PCVOC_PER_POA * flux * aero_rate
                                 / PCVOC_MOLEC_WEIGHT / rate)
        else:
            emit["PCVOC"].append(0.0)
    order.append("PCVOC")

    with open("gas_emit.dat", "w") as f:
        f.write("# generated by make_emissions.py from gas_emit_base.dat\n")
        f.write("# time (s)\n")
        f.write("# rate = scaling parameter\n")
        f.write("# emissions (mol m^{-2} s^{-1})\n")
        f.write("%-8s" % "time" + "".join("%12g" % x for x in times) + "\n")
        f.write("%-8s" % "rate" + "".join("%12g" % x for x in rates) + "\n")
        for name in order:
            f.write("%-8s" % name
                    + "".join("%12.3E" % x for x in emit[name]) + "\n")
    print("Wrote gas_emit.dat")

if __name__ == "__main__":
    assert abs(sum(f for (s, f) in POA_SPLIT) - 1.0) < 1e-12
    write_aero_comps()
    write_gas_emit()
