#!/usr/bin/env python
#
# Write the input files of the 8_organics box case.
#
# Usage: ./make_case.py [case]        (default: urban_forest)
#
# The case is a box: a well-mixed boundary layer of constant depth,
# defined by
#   - its environment: place and date (they set the solar photolysis of
#     MOSAIC and the light of the biogenic emissions), temperature,
#     initial relative humidity (the water vapor is then constant),
#     pressure and depth;
#   - its initial air: the regional air plus the emissions of a city,
#     an increment of CO with the gases and particles emitted with it;
#   - the background air, with which the box exchanges at a constant
#     rate: the regional air;
#   - the biogenic emissions of the land beneath (MEGAN2.1). Nothing
#     else is emitted;
#   - the dry deposition of particles and gases to that land.
# All parameters are below (see 0_README.md for where they come from);
# the variants in CASES change some of them. Other cases are set up by
# adding a variant or by changing the parameters.
#
# Writes organics.spec, temp.dat, pres.dat, height.dat, gas_init.dat,
# gas_back.dat, gas_emit.dat, aero_init_dist.dat, aero_back.dat,
# aero_back_dist.dat, aero_emit.dat, aero_emit_dist.dat, the
# compositions aero_comp_*.dat and, with deposition, drydep_params.dat
# and gas_deposition.dat. The species files gas_data.dat and
# aero_data.dat are not changed.

import math, os, sys

######################################################################
# the case

CASE = {
    # place and date: southeastern United States, 21 June, start at
    # 06:00 local time (CDT, UTC - 5), two days
    "latitude": 33.0,                  # degrees
    "longitude": -86.0,                # degrees (west negative)
    "start_day": 172,                  # day of the year
    "start_local_hour": 6.0,           # h
    "utc_offset": -5.0,                # h
    "hours": 48,
    # boundary layer
    "height": 1000.0,                  # m
    "pressure": 1.0e5,                 # Pa
    "temp_min": 294.0,                 # K, at 06:00 local time
    "temp_max": 303.0,                 # K, at 15:00 local time
    "rel_humidity": 0.90,              # at the start (1)
    # exchange with the background (regional) air: e-folding time of
    # one day
    "exchange_rate": 1.0 / 86400.0,    # s^{-1}
    # urban increment of CO in the initial air
    "delta_co": 200.0,                 # ppb
    # land cover for the biogenic emissions: fractions of the plant
    # functional types of MEGAN2.1 (Guenther et al., 2012, Table 3):
    # 7 broadleaf deciduous temperate tree, 1 needleleaf evergreen
    # temperate tree, 13 cool C3 grass
    "land_cover": {7: 0.50, 1: 0.25, 13: 0.25},
    # PPFD transmission of the atmosphere (ratio of the PPFD above the
    # canopy to that at the top of the atmosphere), the value of the
    # standard conditions of MEGAN (Guenther et al., 2006)
    "ppfd_transmission": 0.6,
    # dry deposition of particles and gases to the forest (see below)
    "deposition": True,
    # wind speed at 20 m for the particle dry deposition: with the
    # roughness length of the forest a friction velocity of 0.28 m/s
    "wind_speed": 2.0,                 # m/s
    # particle-resolved run
    "n_part": 1000,
    "n_repeat": 3,
    "rand_init": 7,
}

CASES = {
    "urban_forest": {},                    # the default
    "forest": {"delta_co": 0.0},           # without the city
    "urban": {"land_cover": {}},           # without the biogenic emissions
    "urban_forest_nodep": {"deposition": False},   # without deposition
}

######################################################################
# regional air (initial air without the urban increment, and the
# background air): idealized round values for the summertime rural
# boundary layer in the eastern United States (ppb). Not listed: 0.

REGIONAL_GAS = {
    "O3": 40.0, "CO": 100.0, "CH4": 1900.0, "C2H6": 1.0,
    "NO": 0.05, "NO2": 0.3, "HNO3": 0.5, "PAN": 0.2, "ONIT": 0.1,
    "SO2": 0.5, "NH3": 0.5,
    "H2O2": 1.0, "CH3OOH": 0.5, "HCHO": 1.0, "ALD2": 0.2, "AONE": 1.5,
    "CH3OH": 3.0, "RCOOH": 0.2,
    "PAR": 1.0, "ETH": 0.2, "OLET": 0.02, "TOL": 0.05, "XYL": 0.02,
    "SOATOL": 0.05, "SOAXYL": 0.02, "BENZ": 0.05,
}

# regional particles: one aged accumulation mode of partly neutralized
# sulfate (molar NH4/SO4 below) and aged organic matter, represented by
# the nonvolatile AORGC (O:C 0.667), which no process of the box forms
REGIONAL_AERO_MASS = {"SO4": 2.0, "AORGC": 2.0}     # ug m^{-3}
REGIONAL_NH4_PER_SO4 = 1.5                           # mol/mol
REGIONAL_MODE = (1.2e-7, 0.20)  # geom. mean diameter (m), log10 geom. std. dev.

######################################################################
# urban increment: emission ratios to CO (mol/mol)

# Los Angeles 2010 (Borbon et al., 2013, Table 2, ppt per ppb of CO),
# lumped into the species of CBM-Z:
#   TOL = toluene 3.18 + ethylbenzene 0.57 + n-propylbenzene 0.11
#         + isopropylbenzene 0.03
#   XYL = m+p-xylene 1.79 + o-xylene 0.67 + 1,3,5-, 1,2,3- and
#         1,2,4-trimethylbenzene 0.31 + 0.24 + 0.62
#         + 1-ethyl-2-methylbenzene 0.12 + 1-ethyl-3,4-dimethylbenzene 0.35
#   PAR = 4 x butanes (4.42 + 3.08) + 5 x pentanes (3.26 + 8.69)
#         + 6 x n-hexane 1.13 + 9, 10 and 11 x n-nonane, n-decane and
#         undecane (0.22, 0.18, 0.29) + 1.5 x propane 11.2
#   OLET = propene 3.74 + 1-butene 0.34 + 2 x 1,3-butadiene 0.35
#          + styrene 0.22
#   OLEI = cis- and trans-2-butene 0.09 + 0.10
#   ALD2 = acetaldehyde 5.42 + propanal 1.13 + butanal 0.25
#   ANOL = ethanol 48.01 + isopropanol 10.60
#   ISOP = isoprene of the traffic 0.30
_PAR = (4 * (4.42 + 3.08) + 5 * (3.26 + 8.69) + 6 * 1.13
        + 9 * 0.22 + 10 * 0.18 + 11 * 0.29 + 1.5 * 11.2)
URBAN_RATIO = dict((k, 1e-3 * v) for (k, v) in {
    "TOL": 3.89, "XYL": 4.10, "BENZ": 1.30, "ETH": 10.35, "C2H6": 18.4,
    "PAR": _PAR, "OLET": 3.74 + 0.34 + 2 * 0.35 + 0.22, "OLEI": 0.19,
    "ALD2": 5.42 + 1.13 + 0.25, "AONE": 11.78, "CH3OH": 21.72,
    "ANOL": 48.01 + 10.60, "ISOP": 0.30}.items())
# NOx, SO2, NH3 and formaldehyde, which Borbon et al. do not give: the
# ratios of the daily emission totals of the 1_urban_plume scenario
# (Riemer et al., 2009), which this scenario used until 2026
URBAN_RATIO.update({"NO": 0.1093, "NO2": 0.00575, "SO2": 0.00862,
                    "NH3": 0.0210, "HCHO": 0.00575})
# precursors of the organic aerosol of CRACMM3: toluene and xylenes
# (TOL and XYL of CRACMM3, tracers here, so that TOL and XYL of CBM-Z
# keep their chemistry); naphthalene and the long-chain alkanes from
# the lumped species, as in the CMAQv5.1 VOC emission updates (Pouliot
# and Pye) for inventories without them
URBAN_DERIVED = [("SOATOL", "TOL", 1.0), ("SOAXYL", "XYL", 0.998),
                 ("NAPH", "XYL", 0.002), ("HC10", "PAR", 0.108)]
# organic gases of intermediate volatility of gasoline exhaust: mass per
# mass of non-methane organic gases (Lu et al., 2020; emission control
# files of CMAQ for CRACMM), with the molecular weights (g/mol) of the
# carbon bond species of CMAQ for the organic gases above
IVOC_PER_NMOG = [("VROCP3ALK", 0.0032), ("VROCP4ALK", 0.0034),
                 ("VROCP5ALK", 0.0019), ("VROCP6ALK", 0.0064),
                 ("VROCP5ARO", 0.0060), ("VROCP6ARO", 0.0253)]
NMOG_MOLEC_WEIGHT = {"PAR": 14.0, "OLET": 42.1, "OLEI": 56.1, "ETH": 28.0,
                     "TOL": 92.1, "XYL": 106.2, "HCHO": 30.0, "ALD2": 44.0,
                     "AONE": 58.1, "CH3OH": 32.0, "ANOL": 46.1}
# molecular weights of the vapors without a particle-phase species
# (g/mol, GC_cracmm3.nml of CMAQ); the others take that of their
# particle-phase species in aero_data.dat
VAPOR_MOLEC_WEIGHT = {"VROCP4ALK": 254.50, "VROCP5ALK": 198.39,
                      "VROCP6ALK": 184.37, "VROCP5ARO": 190.33,
                      "VROCP6ARO": 176.30}

# primary organic matter: 9.1 ug m^{-3} per ppm of CO, the ratio of the
# former inventory of this scenario (Los Angeles 2010: 9.5 ug m^{-3}
# ppm^{-1} of hydrocarbon-like and cooking OA, Hayes et al., 2015),
# divided among the sources as in that inventory (mass fractions)
POA_PER_CO = 9.1                        # ug m^{-3} ppm^{-1}
POA_SOURCES = {"cooking": 0.436, "diesel": 0.325, "gasoline": 0.239}
# particles of the sources (1_urban_plume, Riemer et al., 2009): geom.
# mean diameter (m), log10 geom. std. dev., mass fractions ("POA" is
# split with the volatility profile of the source)
URBAN_MODES = {
    "cooking": (8.64e-8, 0.28, [("POA", 1.0)]),
    "diesel": (5.0e-8, 0.24, [("POA", 0.3), ("BC", 0.7)]),
    "gasoline": (5.0e-8, 0.24, [("POA", 0.8), ("BC", 0.2)]),
}
# volatility profiles of the POA by source (emission control files of
# CMAQ for CRACMM: CMAQ_Control_DESID_cracmm1_aq.nml of v5.4,
# CMAQ_Control_DESID_cracmm2.nml of v5.5): mass per mass of POA of the
# inventory, as (species without the phase prefix, part that CMAQ emits
# in the particles, part that CMAQ emits as vapor).
#   cooking: meat cooking profile of Woody et al. (2016), the default of
#            CMAQ for sources without a profile of their own
#   diesel, gasoline: Lu et al. (2020); the profile of gasoline exhaust
#            adds up to 1.4, the vapors that a filter does not collect
POA_PROFILE = {
    "cooking": [("ROCN2OXY2", 0.35, 0.0), ("ROCP0OXY2", 0.175, 0.175),
                ("ROCP1OXY1", 0.009, 0.091), ("ROCP2OXY2", 0.0, 0.100),
                ("ROCP3OXY2", 0.0, 0.100)],
    "diesel": [("ROCN1ALK", 0.187, 0.0), ("ROCP0ALK", 0.043, 0.043),
               ("ROCP1ALK", 0.0, 0.364), ("ROCP2ALK", 0.0, 0.363)],
    "gasoline": [("ROCN1ALK", 0.228, 0.0), ("ROCP0ALK", 0.054, 0.054),
                 ("ROCP1ALK", 0.0, 0.515), ("ROCP2ALK", 0.0, 0.550)],
}

######################################################################
# biogenic emissions, MEGAN2.1 (Guenther et al., 2012)

# emission factors at standard conditions (ug m^{-2} h^{-1}) of the
# plant functional types 1 to 15 (Table 2)
MEGAN_EF = {
    "isoprene": [600, 3000, 1, 7000, 10000, 7000, 10000, 11000, 2000,
                 4000, 4000, 1600, 800, 200, 1],
    "myrcene": [70, 70, 60, 80, 30, 80, 30, 30, 30, 50, 30,
                0.3, 0.3, 0.3, 0.3],
    "sabinene": [70, 70, 40, 80, 50, 80, 50, 50, 50, 70, 50,
                 0.7, 0.7, 0.7, 0.7],
    "limonene": [100, 100, 130, 80, 80, 80, 80, 80, 60, 100, 60,
                 0.7, 0.7, 0.7, 0.7],
    "3-carene": [160, 160, 80, 40, 30, 40, 30, 30, 30, 100, 30,
                 0.3, 0.3, 0.3, 0.3],
    "t-beta-ocimene": [70, 70, 60, 150, 120, 150, 120, 120, 90, 150, 90,
                       2, 2, 2, 2],
    "beta-pinene": [300, 300, 200, 120, 130, 120, 130, 130, 100, 150, 100,
                    1.5, 1.5, 1.5, 1.5],
    "alpha-pinene": [500, 500, 510, 600, 400, 600, 400, 400, 200, 300, 200,
                     2, 2, 2, 2],
    "other monoterpenes": [180, 180, 170, 150, 150, 150, 150, 150, 110,
                           200, 110, 5, 5, 5, 5],
    "alpha-farnesene": [40, 40, 40, 60, 40, 60, 40, 40, 40, 40, 40,
                        3, 3, 3, 4],
    "beta-caryophyllene": [80, 80, 80, 60, 40, 60, 40, 40, 50, 50, 50,
                           1, 1, 1, 4],
    "other sesquiterpenes": [120, 120, 120, 120, 100, 120, 100, 100, 100,
                             100, 100, 2, 2, 2, 2],
}
# temperature coefficient beta (K^{-1}) and light-dependent fraction LDF
# (Table 4)
MEGAN_PARAM = {
    "isoprene": (0.13, 1.0), "myrcene": (0.1, 0.6), "sabinene": (0.1, 0.6),
    "limonene": (0.1, 0.2), "3-carene": (0.1, 0.2),
    "t-beta-ocimene": (0.1, 0.8), "beta-pinene": (0.1, 0.2),
    "alpha-pinene": (0.1, 0.6), "other monoterpenes": (0.1, 0.4),
    "alpha-farnesene": (0.17, 0.5), "beta-caryophyllene": (0.17, 0.5),
    "other sesquiterpenes": (0.17, 0.5),
}
# model species: as in CRACMM (Pye et al., 2023), monoterpenes with one
# double bond are API and those with two or more LIM; beta-pinene has a
# species of its own in the CBM-Z of this MOSAIC (BPIN; CRACMM3 counts
# it as API). The class "other monoterpenes" (mostly terpinenes,
# terpinolene, phellandrenes, with two double bonds) is counted as LIM
MEGAN_SPECIES = {
    "isoprene": "ISOP", "alpha-pinene": "API", "3-carene": "API",
    "sabinene": "API", "beta-pinene": "BPIN", "limonene": "LIM",
    "myrcene": "LIM", "t-beta-ocimene": "LIM", "other monoterpenes": "LIM",
    "alpha-farnesene": "SESQ", "beta-caryophyllene": "SESQ",
    "other sesquiterpenes": "SESQ",
}
MOLEC_WEIGHT_BIO = {"ISOP": 68.12, "API": 136.23, "BPIN": 136.23,
                    "LIM": 136.23, "SESQ": 204.35}    # g/mol
EMIT_STEP = 1800.0               # s, time step of the emission file

######################################################################
# dry deposition

# particles: the loss function "drydep" of PartMC, the size-resolved
# scheme of Zhang et al. (2001) with its own coefficients and the
# parameters of its land use category 5 (mixed broadleaf and needleleaf
# trees) in the seasonal category 1 (midsummer with lush vegetation),
# Table 3; neutral stability, reference height 20 m. In the order of the
# drydep_params file ("u_mean" is the wind speed of the case).
DRYDEP_PARAMS = [("z_ref", 20.0), ("u_mean", None), ("z_rough", 1.15),
                 ("A", 5.0e-3), ("alpha", 0.8), ("eps_0", 3.0),
                 ("gamma", 0.56), ("C_B", 1.0), ("C_IN", 0.5), ("C_IM", 1.0),
                 ("nu", 2.0), ("beta", 2.0)]
# gases: daytime (10:00 to 15:00) dry deposition velocities (cm/s)
# measured above a forest in the southeastern United States by Nguyen et
# al. (2015, Fig. 3), for the species of CBM-Z and of the organic aerosol
# block that correspond to the measured compounds: nitric acid 3.8,
# hydrogen peroxide 5.2, the isoprene hydroxy hydroperoxides and
# epoxydiols 2.5 (ISPX, ISHP, IEPOX), the hydroperoxy aldehydes 2.4
# (HPLD), formic acid 1.0, the hydroxy nitrates of isoprene 1.5 (ISON)
# and of methacrolein and methyl vinyl ketone 1.5 (ONIT). By analogy:
# IPX (dihydroxy dihydroperoxides) and ROOH (the other hydroperoxides;
# measured 0.8 to 4.1) like the hydroxy hydroperoxides, RCOOH like
# formic acid, INALD (nitrooxy aldehydes) like the hydroxy nitrates, the
# monoterpene nitrates VTRPN and VHONIT like the measured monoterpene
# nitrooxy hydroperoxides 0.8
GAS_VD_DAY = {"HNO3": 3.8, "H2O2": 5.2, "ISPX": 2.5, "ISHP": 2.5,
              "IEPOX": 2.5, "IPX": 2.5, "ROOH": 2.5, "HPLD": 2.4,
              "HCOOH": 1.0, "RCOOH": 1.0, "ONIT": 1.5, "ISON": 1.5,
              "INALD": 1.5, "VTRPN": 0.8, "VHONIT": 0.8}
# the oxygenated vapors of the organic aerosol deposit like the
# hydroperoxides and epoxydiols: their effective Henry constants are
# large (Hodzic et al., 2014; Knote et al., 2015). CMAQ deposits these
# vapors too, but not VHOM and VELHOM (GC_cracmm3 namelist of CMAQ 6.0);
# the alkane-like vapors are hardly soluble and are not deposited here
VAPOR_VD_DAY = 2.5
DEPOSITED_VAPORS = ["VOP3", "VROCN2OXY2", "VROCN2OXY4", "VROCN2OXY8",
                    "VROCN1OXY1", "VROCN1OXY3", "VROCN1OXY6", "VROCP0OXY2",
                    "VROCP0OXY4", "VROCP1OXY1", "VROCP1OXY3", "VROCP2OXY2",
                    "VROCP3OXY2", "VMTN1"]
# diurnal cycle of the gas deposition velocities: the sine of the solar
# elevation relative to its mean from 10:00 to 15:00 local time, at
# least this fraction (at night). The measured velocities are close to
# zero at night (Nguyen et al., 2015, Fig. 3), within their error bars
# below a tenth of the daytime mean
GAS_VD_NIGHT_FRACTION = 0.1

######################################################################

R_GAS = 8.314462618              # J mol^{-1} K^{-1}

def read_data_lines(filename):
    """Return the non-comment lines of a file as lists of words."""
    lines = []
    for line in open(filename):
        words = line.split("#")[0].split()
        if words:
            lines.append(words)
    return lines

def read_aero_data(filename):
    """Return {species: (density (kg/m^3), molecular weight (g/mol))}."""
    return dict((w[0], (float(w[1]), 1e3 * float(w[3].replace("d", "e"))))
                for w in read_data_lines(filename))

def split_poa(source):
    """Split the POA of a source between the particles and the gas phase.

    The particles carry the mass of the POA of the inventory: the
    species of the profile in the order of increasing volatility, as
    far as that mass goes. The rest of the profile is vapor. Returns
    ([(aerosol species, mass per mass of POA)],
    [(vapor, mass per mass of POA)])."""
    particle, vapor, room = [], [], 1.0
    for (name, in_particle, as_vapor) in POA_PROFILE[source]:
        total = in_particle + as_vapor
        part = min(total, room)
        room -= part
        if part > 1e-12:
            particle.append(("A" + name, part))
        if total - part > 1e-12:
            vapor.append(("V" + name, total - part))
    assert abs(room) < 1e-12, "the profile of %s is below 1" % source
    return particle, vapor

def local_hour(case, t):
    return case["start_local_hour"] + t / 3600.0

def temperature(case, t):
    """Temperature (K): half cosines from the minimum at 06:00 to the
    maximum at 15:00 local time and back."""
    h = (local_hour(case, t) - 6.0) % 24.0
    tmin, tmax = case["temp_min"], case["temp_max"]
    if h <= 9.0:
        return tmin + (tmax - tmin) * 0.5 * (1.0 - math.cos(math.pi * h / 9.0))
    return tmax - (tmax - tmin) * 0.5 * (1.0 - math.cos(math.pi * (h - 9.0)
                                                        / 15.0))

def sin_solar_elevation(case, t):
    """Sine of the solar elevation, as MOSAIC computes it
    (gas/solarzenithangle.f90, with the time of PartMC's mosaic.F90)."""
    start_time = 3600.0 * (case["start_local_hour"] - case["utc_offset"])
    tsec = case["start_day"] * 86400.0 + start_time + t
    tmar21 = (79 * 24 + 12) * 3600.0
    if tsec >= tmar21:
        tsec -= tmar21
    else:
        tsec += ((365 - 79) * 24 - 12) * 3600.0
    dec = 0.4092797 * math.sin(1.992385e-7 * tsec)
    lat = math.radians(case["latitude"])
    lon = math.radians(case["longitude"])
    return (math.sin(lat) * math.sin(dec) + math.cos(lat) * math.cos(dec)
            * math.cos(lon + 7.272205e-5 * tsec))

def ppfd_toa(case):
    """PPFD at the top of the atmosphere (umol m^{-2} s^{-1}),
    Guenther et al. (2006), Eq. (13)."""
    return 3000.0 + 99.0 * math.cos(2.0 * 3.14 * (case["start_day"] - 10)
                                    / 365.0)

def ppfd(case, t):
    """PPFD above the canopy (umol m^{-2} s^{-1})."""
    return (case["ppfd_transmission"] * ppfd_toa(case)
            * max(sin_solar_elevation(case, t), 0.0))

def daily_means(case):
    """Mean PPFD above the canopy and mean temperature of the first day."""
    times = [60.0 * (i + 0.5) for i in range(24 * 60)]
    return (sum(ppfd(case, t) for t in times) / len(times),
            sum(temperature(case, t) for t in times) / len(times))

def gamma_light_temp(case, t, p_daily, t_daily):
    """Canopy-scale response of the light-dependent emissions to light
    and temperature, the parameterized canopy environment emission
    activity of MEGAN (PCEEA, Guenther et al., 2006, Eqs. 11, 12 and 14,
    with LAI 5, for which Eq. 15 gives 1). It is 1 at the standard
    conditions (solar elevation 60 degrees, PPFD transmission 0.6,
    303 K, daily means 400 umol m^{-2} s^{-1} and 297 K)."""
    sin_a = sin_solar_elevation(case, t)
    if sin_a <= 0.0:
        return 0.0
    phi = ppfd(case, t) / (sin_a * ppfd_toa(case))
    gamma_p = sin_a * (2.46 * (1.0 + 0.0005 * (p_daily - 400.0)) * phi
                       - 0.9 * phi**2)
    t_opt = 313.0 + 0.6 * (t_daily - 297.0)
    e_opt = 1.75 * math.exp(0.08 * (t_daily - 297.0))
    x = (1.0 / t_opt - 1.0 / temperature(case, t)) / 0.00831
    ct1, ct2 = 80.0, 200.0
    gamma_t = e_opt * ct2 * math.exp(ct1 * x) \
        / (ct2 - ct1 * (1.0 - math.exp(ct2 * x)))
    return max(gamma_p, 0.0) * gamma_t

def biogenic_emissions(case):
    """Emission times (s) and {species: [mol m^{-2} s^{-1}]}: the
    light-dependent fraction LDF of each compound class follows the
    canopy-scale response of isoprene to light and temperature (PCEEA),
    the rest the exponential temperature response of MEGAN2.1 (Eq. 11),
    1 at the standard temperature of 303 K. Leaf age, soil moisture and
    CO2 factors are 1."""
    times = [EMIT_STEP * i for i in range(int(case["hours"] * 3600
                                              / EMIT_STEP))]
    emit = {}
    if not case["land_cover"]:
        return times, emit
    p_daily, t_daily = daily_means(case)
    for (cls, factors) in MEGAN_EF.items():
        ef = sum(frac * factors[pft - 1]
                 for (pft, frac) in case["land_cover"].items())
        beta, ldf = MEGAN_PARAM[cls]
        spec = MEGAN_SPECIES[cls]
        values = emit.setdefault(spec, [0.0] * len(times))
        for (i, t) in enumerate(times):
            tm = t + 0.5 * EMIT_STEP
            gamma = ((1.0 - ldf) * math.exp(beta * (temperature(case, tm)
                                                    - 303.0))
                     + ldf * gamma_light_temp(case, tm, p_daily, t_daily))
            values[i] += ef * gamma * 1e-6 / MOLEC_WEIGHT_BIO[spec] / 3600.0
    return times, emit

def gas_deposition_profile(case):
    """Times (s), rate scalings (the diurnal cycle) and {species:
    daytime velocity (m/s)} of the gas dry deposition."""
    times = [EMIT_STEP * i for i in range(int(case["hours"] * 3600
                                              / EMIT_STEP))]
    t10 = (10.0 - case["start_local_hour"]) * 3600.0
    ref = [max(sin_solar_elevation(case, t10 + 60.0 * (i + 0.5)), 0.0)
           for i in range(5 * 60)]
    s_ref = sum(ref) / len(ref)
    rates = [max(max(sin_solar_elevation(case, t + 0.5 * EMIT_STEP), 0.0)
                 / s_ref, GAS_VD_NIGHT_FRACTION) for t in times]
    vd = dict((k, 1e-2 * v) for (k, v) in GAS_VD_DAY.items())
    vd.update((k, 1e-2 * VAPOR_VD_DAY) for k in DEPOSITED_VAPORS)
    return times, rates, vd

def mode_text(name, comp_file, num_conc, diam, log10_std):
    return ("mode_name %s\nmass_frac %s\ndiam_type geometric\n"
            "mode_type log_normal\nnum_conc %.6e\ngeom_mean_diam %g\n"
            "log10_geom_std_dev %g\n\n"
            % (name, comp_file, num_conc, diam, log10_std))

def mode_num_conc(masses, diam, log10_std, aero_data):
    """Number concentration (m^{-3}) of a log-normal mode with the
    species masses (ug m^{-3})."""
    total = sum(masses.values())
    density = total / sum(m / aero_data[s][0] for (s, m) in masses.items())
    ln_sigma = log10_std * math.log(10.0)
    mean_vol = math.pi / 6.0 * diam**3 * math.exp(4.5 * ln_sigma**2)
    return total * 1e-9 / (density * mean_vol)

def write_comp(filename, masses):
    total = sum(masses.values())
    with open(filename, "w") as f:
        f.write("# generated by make_case.py\n#       mass fractions\n")
        for (name, mass) in sorted(masses.items()):
            f.write("%-14s %.6g\n" % (name, mass / total))

def write(filename, text):
    with open(filename, "w") as f:
        f.write(text)

SPEC = """run_type particle               # particle-resolved run
output_prefix out/organics      # prefix of output files
n_repeat %(n_repeat)d                      # number of Monte Carlo repeats
n_part %(n_part)d                     # total number of particles
restart no                      # whether to restart from saved state (yes/no)
do_select_weighting no          # whether to select weighting explicitly (yes/no)

t_max %(t_max)d                    # total simulation time (s)
del_t 60                        # timestep (s)
t_output 3600                   # output interval (0 disables) (s)
t_progress 3600                 # progress printing interval (0 disables) (s)

do_camp_chem no                 # whether to use CAMP for chemistry
do_tchem no                     # whether to use TChem for chemistry

gas_data gas_data.dat           # file containing gas data
gas_init gas_init.dat           # initial gas concentrations

aerosol_data aero_data.dat      # file containing aerosol data
do_fractal no                   # whether to do fractal treatment
aerosol_init aero_init_dist.dat # aerosol initial condition file

temp_profile temp.dat           # temperature profile file
pressure_profile pres.dat       # pressure profile file
height_profile height.dat       # height profile file
gas_emissions gas_emit.dat      # gas emissions file
gas_background gas_back.dat     # background gas concentrations file
aero_emissions aero_emit.dat    # aerosol emissions file
aero_background aero_back.dat   # aerosol background file
%(loss_lines)s
rel_humidity %(rel_humidity)g               # initial relative humidity (1)
latitude %(latitude)g                     # latitude (degrees, -90 to 90)
longitude %(longitude)g                   # longitude (degrees, -180 to 180)
altitude 0                      # altitude (m)
start_time %(start_time)d                # start time (s since 00:00 UTC)
start_day %(start_day)d                   # start day of year (UTC)

do_coagulation yes              # whether to do coagulation (yes/no)
coag_kernel brown               # coagulation kernel
do_condensation no              # whether to do condensation (yes/no)
do_mosaic yes                   # whether to do MOSAIC (yes/no)
do_optical yes                  # whether to compute optical props (yes/no)
do_fixed_photolysis no          # whether to use fixed photolysis rates (yes/no)
do_nucleation no                # whether to do nucleation (yes/no)
do_immersion_freezing no        # whether to do freezing (yes/no)

rand_init %(rand_init)d                     # random initialization (0 to use time)
allow_doubling yes              # whether to allow doubling (yes/no)
allow_halving yes               # whether to allow halving (yes/no)
record_removals yes             # whether to record particle removals (yes/no)
do_parallel no                  # whether to run in parallel (yes/no)
"""

def main(name):
    case = dict(CASE)
    case.update(CASES[name])
    aero_data = read_aero_data("aero_data.dat")
    gas_species = [w[0] for w in read_data_lines("gas_data.dat")]
    t_max = int(case["hours"] * 3600)
    start_time = int(round(3600.0 * (case["start_local_hour"]
                                     - case["utc_offset"]))) % 86400
    print("case %s" % name)

    # environment
    times = range(0, t_max + 1, 3600)
    write("temp.dat", "# generated by make_case.py\n# time (s)\n"
          "# temp (K)\n" + "time " + " ".join("%d" % t for t in times)
          + "\ntemp " + " ".join("%.3f" % temperature(case, t)
                                 for t in times) + "\n")
    write("pres.dat", "# time (s)\n# pressure (Pa)\ntime 0\n"
          "pressure %g\n" % case["pressure"])
    write("height.dat", "# time (s)\n# height (m)\ntime 0\n"
          "height %g\n" % case["height"])
    t0 = temperature(case, 0.0)
    cair = case["pressure"] / (R_GAS * t0)                 # mol m^{-3}
    ppb_of_mass = lambda mass, mw: mass * 1e-6 / mw / cair * 1e9

    # initial gases: regional air + urban increment
    dco = case["delta_co"]
    gas = dict(REGIONAL_GAS)
    ratio = dict(URBAN_RATIO)
    for (new, ref, r) in URBAN_DERIVED:
        ratio[new] = r * ratio[ref]
    nmog = sum(ratio.get(k, 0.0) * mw
               for (k, mw) in NMOG_MOLEC_WEIGHT.items())  # g per mol of CO
    for (k, frac) in IVOC_PER_NMOG:
        mw = VAPOR_MOLEC_WEIGHT.get(k, aero_data.get("A" + k[1:], (0, 0))[1])
        ratio[k] = frac * nmog / mw
    for (k, r) in ratio.items():
        gas[k] = gas.get(k, 0.0) + r * dco
    gas["CO"] = REGIONAL_GAS["CO"] + dco

    # initial particles: regional mode + urban modes
    so4 = REGIONAL_AERO_MASS["SO4"]
    regional = dict(REGIONAL_AERO_MASS)
    regional["NH4"] = (REGIONAL_NH4_PER_SO4 * so4 / aero_data["SO4"][1]
                       * aero_data["NH4"][1])
    write_comp("aero_comp_regional.dat", regional)
    diam, log10_std = REGIONAL_MODE
    num_regional = mode_num_conc(regional, diam, log10_std, aero_data)
    init_modes = mode_text("regional", "aero_comp_regional.dat",
                           num_regional, diam, log10_std)
    print("regional particles: %s ug m^-3, %.0f cm^-3"
          % (", ".join("%s %.3g" % kv for kv in sorted(regional.items())),
             num_regional * 1e-6))
    poa_total = POA_PER_CO * dco * 1e-3                    # ug m^{-3}
    for (source, (diam, log10_std, comp)) in sorted(URBAN_MODES.items()):
        if dco <= 0.0:
            break
        poa = poa_total * POA_SOURCES[source]
        poa_frac = dict(comp)["POA"]
        particle, vapor = split_poa(source)
        masses = {}
        for (s, frac) in comp:
            if s == "POA":
                for (sp, f) in particle:
                    masses[sp] = masses.get(sp, 0.0) + poa * f
            else:
                masses[s] = masses.get(s, 0.0) + poa * frac / poa_frac
        for (sp, f) in vapor:
            gas[sp] = gas.get(sp, 0.0) + ppb_of_mass(poa * f,
                                                     aero_data["A" + sp[1:]][1])
        comp_file = "aero_comp_%s.dat" % source
        write_comp(comp_file, masses)
        num = mode_num_conc(masses, diam, log10_std, aero_data)
        init_modes += mode_text(source, comp_file, num, diam, log10_std)
        print("%s particles: POA %.3g ug m^-3, total %.3g ug m^-3, "
              "%.0f cm^-3%s" % (source, poa, sum(masses.values()),
                                num * 1e-6,
                                "; vapor " + ", ".join(
                                    "%s %.3g ug m^-3" % (sp, poa * f)
                                    for (sp, f) in vapor) if vapor else ""))
    write("aero_init_dist.dat", "# generated by make_case.py\n" + init_modes)

    for k in gas:
        assert k in gas_species, "unknown gas species %s" % k
    write("gas_init.dat", "# generated by make_case.py\n"
          "# species    initial mixing ratio (ppb)\n"
          + "".join("%-12s %.6e\n" % (k, gas[k]) for k in gas_species
                    if gas.get(k, 0.0) > 0.0))

    # background air: the regional air
    rate = case["exchange_rate"]
    write("gas_back.dat", "# generated by make_case.py\n# time (s)\n"
          "# rate (s^{-1})\n# concentrations (ppb)\n"
          "time %-12s 0\nrate %-12s %.6e\n" % ("", "", rate)
          + "".join("%-17s %.6e\n" % (k, REGIONAL_GAS[k]) for k in gas_species
                    if REGIONAL_GAS.get(k, 0.0) > 0.0))
    write("aero_back.dat", "# generated by make_case.py\n# time (s)\n"
          "# rate (s^{-1})\n# aerosol distribution filename\n"
          "time 0\nrate %.6e\ndist aero_back_dist.dat\n" % rate)
    diam, log10_std = REGIONAL_MODE
    write("aero_back_dist.dat", "# generated by make_case.py\n"
          + mode_text("regional", "aero_comp_regional.dat", num_regional,
                      diam, log10_std))

    # emissions: biogenic gases only
    write("aero_emit.dat", "# generated by make_case.py: no particle "
          "emissions\ntime 0\nrate 0\ndist aero_emit_dist.dat\n")
    write("aero_emit_dist.dat", "# no particle emissions\n")
    times, emit = biogenic_emissions(case)
    if emit:
        order = ["ISOP", "API", "BPIN", "LIM", "SESQ"]
        write("gas_emit.dat", "# generated by make_case.py: biogenic "
              "emissions (MEGAN2.1)\n# time (s)\n# rate = scaling "
              "parameter\n# emissions (mol m^{-2} s^{-1})\n"
              + "%-8s" % "time" + "".join("%11d" % t for t in times) + "\n"
              + "%-8s" % "rate" + "".join("%11g" % 1 for t in times) + "\n"
              + "".join("%-8s" % s + "".join("%11.3E" % x for x in emit[s])
                        + "\n" for s in order))
        for s in order:
            day = [x * MOLEC_WEIGHT_BIO[s] * 1e6 * 3600.0
                   for x in emit[s][:int(86400 / EMIT_STEP)]]
            print("%-5s emission: daily mean %7.1f, maximum %7.1f, "
                  "minimum %6.1f ug m^-2 h^-1"
                  % (s, sum(day) / len(day), max(day), min(day)))
    else:
        write("gas_emit.dat", "# generated by make_case.py: no gas "
              "emissions\ntime 0\nrate 0\nNO 0.0\n")

    # dry deposition
    if case["deposition"]:
        params = [(k, case["wind_speed"] if k == "u_mean" else v)
                  for (k, v) in DRYDEP_PARAMS]
        write("drydep_params.dat", "# generated by make_case.py: particle "
              "dry deposition, Zhang et al. (2001), land use category 5, "
              "seasonal category 1\n"
              + "".join("%-8s %g\n" % kv for kv in params))
        times_d, rates_d, vd = gas_deposition_profile(case)
        for k in vd:
            assert k in gas_species, "unknown gas species %s" % k
        write("gas_deposition.dat", "# generated by make_case.py: gas dry "
              "deposition\n# time (s)\n# rate = diurnal scaling\n"
              "# dry deposition velocities (m/s)\n"
              + "%-11s" % "time" + "".join("%10d" % t for t in times_d) + "\n"
              + "%-11s" % "rate" + "".join("%10.4f" % r for r in rates_d)
              + "\n" + "".join("%-11s" % k + "%10.3E" % vd[k] * len(times_d)
                               + "\n" for k in gas_species if k in vd))
        loss_lines = ("loss_function drydep            # loss function "
                      "specification\ndrydep_params drydep_params.dat "
                      "# particle dry deposition parameters\ngas_deposition "
                      "gas_deposition.dat # gas dry deposition velocities\n")
        params = dict(params)
        print("dry deposition: friction velocity %.2f m/s; gas velocities "
              "x %.2f at noon, x %.2f at night (daily mean %.2f)"
              % (0.4 * params["u_mean"] / math.log(params["z_ref"]
                                                   / params["z_rough"]),
                 max(rates_d), min(rates_d), sum(rates_d[:48]) / 48.0))
    else:
        for f in ("drydep_params.dat", "gas_deposition.dat"):
            if os.path.exists(f):
                os.remove(f)
        loss_lines = ("loss_function none              # loss function "
                      "specification\n")

    write("organics.spec", SPEC % dict(case, t_max=t_max,
                                       start_time=start_time,
                                       loss_lines=loss_lines))
    print("initial gases (ppb): " + ", ".join(
        "%s %.3g" % (k, gas[k]) for k in ["CO", "NO", "NO2", "O3", "SO2",
                                          "NH3", "TOL", "XYL", "PAR", "HC10",
                                          "VROCP6ARO", "VROCP2ALK"]
        if k in gas))
    print("temperature %.1f to %.1f K, initial RH %.2f; noon PPFD %.0f "
          "umol m^-2 s^-1" % (min(temperature(case, t) for t in times),
                              max(temperature(case, t) for t in times),
                              case["rel_humidity"],
                              max(ppfd(case, t) for t in times)))

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "urban_forest")
