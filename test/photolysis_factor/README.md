# Photolysis factor test case

This tests the optional scenario line `photolysis_factor_profile`, a time profile of a factor that multiplies all photolysis rates of MOSAIC (for example the ratio of the actual to the clear-sky photolysis under clouds or haze). Every run is 2 h of MOSAIC gas and aerosol chemistry with the initial air of the `mosaic` test case, no emissions, no exchange with the background and no coagulation, so the runs differ only in their photolysis.

- Test 1, fixed photolysis rates (`do_fixed_photolysis yes`): the factor 0.5 (`factor_half.dat`) with the rates of `rates.dat` must give the same result as the rates of `rates_half.dat`, which are exact decimal halves; the full rates must give a different result.
- Test 2, clear-sky photolysis from the solar zenith angle: at noon with the factor 0 (`factor_zero.dat`) the chemistry must be the same as at midnight; the noon run without the factor must differ. The initial air forms no methylglyoxal or IEPOX, whose particle uptake depends on day and night, so the photolysis is the only difference between noon and midnight.

Both comparisons should agree to rounding (relative tolerance 1e-12).
