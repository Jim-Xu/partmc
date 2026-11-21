# beta-pinene test case (MOSAIC)

This tests the gas species `BPIN` (beta-pinene) of the MOSAIC interface. beta-pinene has its own species because it reacts with ozone at about a fifth of the rate of alpha-pinene (`API`), while its reaction with the nitrate radical forms the monoterpene nitrates; its products are those of the other monoterpenes.

The run is 1 h in the dark (polar night) at 298 K with 100 ppb O3 and 1 ppb BPIN and no other reactive gases, so the only reaction is BPIN + O3, which consumes one O3 per BPIN. Test 1 compares the mixing ratios of O3 and BPIN with the analytic solution for the IUPAC rate constant, 1.39e-15 exp(-1280/T) cm3 molecule-1 s-1 (Cox et al., 2020, Atmos. Chem. Phys. 20, 13497, datasheet Ox_VOC19); `compute_reference.py` writes the reference in the column layout of `extract_gas`, with zeros for all other species. BPIN falls to 0.847 ppb (alpha-pinene would fall to 0.46 ppb).

The rate constants and products of all BPIN reactions are checked term by term by the rate test of the VBS port (vbs_work/tests/rate_test.py).
