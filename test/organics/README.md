
# Organic matter test case

This tests the per-particle properties of the organic matter: the atomic oxygen-to-carbon ratio O:C, the ratio of organic mass to organic carbon mass OM:OC, the organic mass and the organic carbon mass. They are computed from the species properties O:C and OM:OC in the last two columns of `aero_data.dat`.

The population consists of three monodisperse subpopulations with the species A (inorganic), B and C (organic with different O:C and OM:OC) and D (not organic, like black carbon):

1. organic particles with B and C,
2. mixed particles with A, B and D,
3. inorganic particles with A only, for which O:C and OM:OC are zero.

There is no coagulation and no chemistry, so all particles of a subpopulation stay identical, and the reference values in `ref_o_c.txt` follow directly from the input files (see `compute_reference.py`).

The second test uses a species data file without the O:C and OM:OC columns (`aero_data_7col.dat`). Then no species counts as organic carbon and all values are zero.
