# Acidity of the IEPOX uptake (MOSAIC)

This tests the optional spec line `iepox_acidity` (after `do_fixed_photolysis`), which chooses the acidity used by MOSAIC for the acid-catalyzed uptake of isoprene epoxydiols (IEPOX -> AISO3NOS, AISO3OS; IEPOX is a gas species of the CRACMM3 organic aerosol block): `isorropia` (the default) for one ISORROPIA II calculation on the totals of the gases and all particles, as CMAQ does for each of its modes, or `mosaic` for MOSAIC's own acidity of each particle. Each run is 1 h in the dark (polar night) with 1 ppb IEPOX over wet ammonium bisulfate particles at 70% RH, with flat weighting and no coagulation.

- Test 1: the run without the `iepox_acidity` line must be identical to the run with `iepox_acidity isorropia`.
- Test 2: the run with `iepox_acidity mosaic` must differ from it (MOSAIC's own pH of these particles is lower than the ISORROPIA II pH, so the uptake is faster).

The uptake rate itself is checked against the CMAQ formula with the ISORROPIA II acidity by the uptake test of the VBS port (vbs_work/tests/uptake_test.py, check 5; the formula is the same for CRACMM3), the split into AISO3NOS and AISO3OS by vbs_work/cracmm/upt_test.py.
