# Monoterpene + OH test case (MOSAIC)

This tests that the reactions of the monoterpenes with OH in the MOSAIC interface consume OH. In the MOSAIC library of 2012 these were tracer reactions: the monoterpene was lost, OH was not. They now form the gas-phase products of the monoterpene reactions of CMAQ cb6r5_ae7_aq (CB6r5: TERP + OH, R172), like the reactions with O3 and NO3.

The run is 10 min at 298 K with 1000 ppb H2O2 and 1000 ppb alpha-pinene (`API`) and a fixed photolysis frequency of H2O2 of 1e-5 s-1 as the only source of OH (20 ppt s-1). Almost all OH reacts with API (97%; the rest with H2O2), so API decreases by the amount of OH formed, 11.6 ppb. Test 1 compares the mixing ratios of API and H2O2 with this solution; `compute_reference.py` writes the reference file. If API did not consume OH, OH would be limited by its reaction with H2O2 alone and 31% of the API would react instead of 1.2%.

The rate constants and products of all monoterpene reactions are checked term by term by the rate test of the VBS port (vbs_work/tests/rate_test.py).
