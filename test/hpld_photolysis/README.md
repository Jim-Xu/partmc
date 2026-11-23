# HPLD photolysis test case (MOSAIC)

This tests the gas species `HPLD` of the MOSAIC interface, the hydroperoxy aldehyde that the isoprene peroxy radical forms by isomerization (CMAQ cb6r5_ae7_aq, CB6r5: ISO2 = HO2 + HPLD, R155). MOSAIC photolyzes it to OH and `ISOPRD` (R163) at 0.055 times the photolysis frequency of NO2, which is the rate of the CMAQ table HPALD (quantum yield 1) at the surface.

The run is 1 h at 298 K with 10 ppb HPLD and a fixed photolysis frequency of NO2 of 5e-3 s-1 as the only photolysis. Without NOx HPLD has no other loss, so it decays exponentially at 2.75e-4 s-1, to 37% after 1 h. Test 1 compares the mixing ratio of HPLD with this solution; `compute_reference.py` writes the reference file.

The rate constants and products of the isoprene reactions that were changed with HPLD (ISOPP + HO2, the isomerization of ISOPP, ISPX + OH, HPLD + NO3) are checked term by term by the rate test of the VBS port (vbs_work/tests/rate_test.py).
