#!/usr/bin/env python

# Generate the emission files of the 8_organics scenario (urban plume
# case) for the organic aerosol of CRACMM3:
#
#   gas_emit.dat              gas emissions, from gas_emit_base.dat plus
#                             the precursors of the organic aerosol
#   aero_emit_comp_*.dat      composition of the emitted combustion
#                             particles, with the semivolatile primary
#                             organic matter of CRACMM3 instead of OC
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
    # precursors of the aromatic SOA (the species TOL and XYL of
    # CRACMM3): same emissions as the CBM-Z species (SOAXYL without the
    # naphthalene part of XYL)
    ("SOATOL", "TOL", 1.0),
    ("SOAXYL", "XYL", 0.998),
    # naphthalene and long-chain alkanes from the lumped species of an
    # inventory that lacks them (CMAQv5.1 VOC emission updates); the
    # long-chain alkanes are the species HC10 of CRACMM3
    ("NAPH", "XYL", 0.002),
    ("HC10", "PAR", 0.108),
    # benzene: urban emission ratio to CO in Los Angeles, 1.30 pptv/ppbv
    # (Borbon et al., 2013, Table 2)
    ("BENZ", "CO", 1.3e-3),
    # terpenes relative to isoprene, from the global annual emissions
    # of MEGAN2.1 (Guenther et al., 2012, Table 6):
    # isoprene 535 Tg, alpha-pinene 66.1 Tg, beta-pinene 18.9 Tg and
    # camphene 4.0 Tg (both react slowly with ozone; BPIN), all other
    # monoterpenes 73.3 Tg, sesquiterpenes 29.0 Tg
    ("API", "ISOP", (66.1 / 136.24) / (535.0 / 68.12)),
    ("BPIN", "ISOP", (22.9 / 136.24) / (535.0 / 68.12)),
    ("LIM", "ISOP", (73.3 / 136.24) / (535.0 / 68.12)),
    ("SESQ", "ISOP", (29.0 / 204.36) / (535.0 / 68.12)),
]

######################################################################
# primary organic matter (POA)

# volatility profiles of the emitted POA by source, from the emission
# control files of CMAQ for CRACMM (CMAQ_Control_DESID_cracmm1_aq.nml of
# v5.4, CMAQ_Control_DESID_cracmm2.nml of v5.5): mass emitted per mass
# of POA of the inventory as (species without the phase prefix, part
# that CMAQ emits in the particles, part that CMAQ emits as vapor).
#   cooking:  meat cooking profile of Woody et al. (2016), the default
#             of CMAQ for sources without a profile of their own, with
#             an O:C of 0.1 to 0.2 (Mohr et al., 2009)
#   diesel, gasoline: Lu et al. (2020). The profile of gasoline exhaust
#             adds up to 1.4: the vapors that a filter does not collect
POA_PROFILE = {
    "cooking": [("ROCN2OXY2", 0.35, 0.0), ("ROCP0OXY2", 0.175, 0.175),
                ("ROCP1OXY1", 0.009, 0.091), ("ROCP2OXY2", 0.0, 0.100),
                ("ROCP3OXY2", 0.0, 0.100)],
    "diesel": [("ROCN1ALK", 0.187, 0.0), ("ROCP0ALK", 0.043, 0.043),
               ("ROCP1ALK", 0.0, 0.364), ("ROCP2ALK", 0.0, 0.363)],
    "gasoline": [("ROCN1ALK", 0.228, 0.0), ("ROCP0ALK", 0.054, 0.054),
                 ("ROCP1ALK", 0.0, 0.515), ("ROCP2ALK", 0.0, 0.550)],
}

# composition of the emitted combustion particles (mass fractions), as
# in the 1_urban_plume scenario; "POA" is split according to the
# profile of the source
AERO_EMIT_COMP = {
    "aero_emit_comp_cooking.dat": ("cooking", [("POA", 1.0)]),
    "aero_emit_comp_diesel.dat": ("diesel", [("POA", 0.3), ("BC", 0.7)]),
    "aero_emit_comp_gasol.dat": ("gasoline", [("POA", 0.8), ("BC", 0.2)]),
}

######################################################################
# organic gases of intermediate volatility (IVOC) of gasoline exhaust:
# mass emitted per mass of non-methane organic gases (Lu et al., 2020;
# the same CMAQ files). The non-methane organic gases are the organic
# species of gas_emit_base.dat below, with the molecular weights (g/mol)
# of the carbon bond species in CMAQ (GC_cb6r5_ae7_aq.nml). The profile
# of diesel exhaust is not used: the inventory does not tell the diesel
# part of the organic gases.
IVOC_PER_NMOG = [("VROCP3ALK", 0.0032), ("VROCP4ALK", 0.0034),
                 ("VROCP5ALK", 0.0019), ("VROCP6ALK", 0.0064),
                 ("VROCP5ARO", 0.0060), ("VROCP6ARO", 0.0253)]
NMOG_MOLEC_WEIGHT = {"PAR": 14.0, "OLET": 42.1, "OLEI": 56.1, "ETH": 28.0,
                     "TOL": 92.1, "XYL": 106.2, "HCHO": 30.0, "ALD2": 44.0,
                     "AONE": 58.1, "CH3OH": 32.0, "ANOL": 46.1}

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

def read_molec_weights(filename):
    """Return {species: molecular weight (kg/mol)} from an aero_data file."""
    return dict((w[0], float(w[3].replace("d", "e")))
                for w in read_data_lines(filename))

def split_poa(source):
    """Split the POA of a source between the particles and the gas phase.

    The particles carry the mass of the POA of the inventory: the
    species of the profile in the order of increasing volatility, as
    far as that mass goes. The rest of the profile is emitted as vapor.
    Returns ([(aerosol species, mass per mass of POA)],
    [(vapor, mass per mass of POA)])."""
    particle, vapor, room = [], [], 1.0
    for (name, in_particle, as_vapor) in POA_PROFILE[source]:
        total = in_particle + as_vapor
        assert total >= in_particle > -1e-12
        part = min(total, room)
        room -= part
        if part > 1e-12:
            particle.append(("A" + name, part))
        if total - part > 1e-12:
            vapor.append(("V" + name, total - part))
    assert abs(room) < 1e-12, "the profile of %s is below 1" % source
    return particle, vapor

def expand_comp(source, comp):
    """Replace the entry "POA" of a composition by the species of the
    volatility profile of the source."""
    species = []
    for (name, frac) in comp:
        if name == "POA":
            species.extend((s, frac * f) for (s, f) in split_poa(source)[0])
        else:
            species.append((name, frac))
    return species

def poa_mass_flux(dist_filename, densities):
    """POA mass emitted per area and time (kg m^{-2} s^{-1}) by the
    modes of the emission distribution, for a unit rate scaling, as
    {source: flux}."""
    flux = {}
    for mode in read_modes(dist_filename):
        if mode["mass_frac"] not in AERO_EMIT_COMP:
            continue
        assert mode["mode_type"] == "log_normal"
        assert mode["diam_type"] == "geometric"
        (source, comp) = AERO_EMIT_COMP[mode["mass_frac"]]
        species = expand_comp(source, comp)
        total = sum(f for (s, f) in species)
        density = total / sum(f / densities[s] for (s, f) in species)
        poa_frac = dict(comp)["POA"] / total
        ln_sigma = float(mode["log10_geom_std_dev"]) * math.log(10.0)
        mean_vol = math.pi / 6.0 * float(mode["geom_mean_diam"])**3 \
            * math.exp(4.5 * ln_sigma**2)
        mode_flux = float(mode["num_conc"]) * mean_vol * density * poa_frac
        print("  mode %-10s POA mass flux %.3e kg m^{-2} s^{-1}"
              % (mode["mode_name"], mode_flux))
        flux[source] = flux.get(source, 0.0) + mode_flux
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
        (source, species) = comp
        with open(filename, "w") as f:
            f.write("# generated by make_emissions.py\n")
            f.write("#       mass fractions\n")
            for (name, frac) in expand_comp(source, species):
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

    # organic gases of intermediate volatility, from the mass of the
    # non-methane organic gases
    molec_weight = read_molec_weights("aero_data.dat")
    molec_weight.update({"VROCP4ALK": 254.50e-3, "VROCP5ALK": 198.39e-3,
                         "VROCP6ALK": 184.37e-3, "VROCP5ARO": 190.33e-3,
                         "VROCP6ARO": 176.30e-3})    # GC_cracmm3.nml of CMAQ
    nmog = [sum(1e-3 * mw * emit[name][i]
                for (name, mw) in NMOG_MOLEC_WEIGHT.items())
            for i in range(len(times))]              # kg m^{-2} s^{-1}
    for (name, frac) in IVOC_PER_NMOG:
        mw = molec_weight.get(name, molec_weight.get("A" + name[1:]))
        emit[name] = [frac * x / mw for x in nmog]
        order.append(name)

    # the part of the POA profiles that the particles do not carry is
    # emitted as vapor. It follows the POA of the particle emissions,
    # which have their own times and rate scalings
    aero = dict((w[0], w[1:]) for w in read_data_lines("aero_emit.dat"))
    assert len(set(aero["dist"])) == 1
    aero_times = [float(x) for x in aero["time"]]
    aero_rates = [float(x) for x in aero["rate"]]
    print("POA emissions of %s:" % aero["dist"][0])
    flux = poa_mass_flux(aero["dist"][0], read_densities("aero_data.dat"))
    for source in sorted(flux):
        for (name, per_poa) in split_poa(source)[1]:
            mw = molec_weight["A" + name[1:]]
            values = []
            for (time, rate) in zip(times, rates):
                aero_rate = rate_at(aero_times, aero_rates, time)
                if aero_rate > 0.0:
                    # the value is multiplied by the gas rate scaling in
                    # PartMC
                    assert rate > 0.0
                    values.append(per_poa * flux[source] * aero_rate
                                  / mw / rate)
                else:
                    values.append(0.0)
            if name in emit:
                emit[name] = [a + b for (a, b) in zip(emit[name], values)]
            else:
                emit[name] = values
                order.append(name)
            print("  %s of %s exhaust as vapor: %.3f per mass of POA"
                  % (name, source, per_poa))

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
    write_aero_comps()
    write_gas_emit()
