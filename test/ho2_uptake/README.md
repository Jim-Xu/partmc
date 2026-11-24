# HO2 uptake test case (MOSAIC)

This tests the loss of HO2 on the particles in the MOSAIC interface, the heterogeneous reaction `HET_HO2` of the CRACMM3 mechanism of CMAQ 6.0: HO2 is lost on the wet surface of the particles without products, with the uptake coefficient 0.2 (Ivatt et al., Nature Geosci. 15, 536, 2022) and limited by gas-phase diffusion. PartMC computes the loss frequency of every time step from the wet diameters of its particles (`mosaic_ho2_uptake_freq()` in `src/mosaic.F90`) and MOSAIC applies it as a first-order loss in its gas-phase chemistry. The optional spec file line `ho2_uptake_coeff` sets the uptake coefficient; 0 switches the loss off.

The runs are 10 min in the dark at 298 K and 1000 hPa with 1e-5 ppb of HO2 as the only gas and 1e9 m-3 solid ammonium sulfate particles of 200 nm diameter (relative humidity 30%), so that HO2 has no other loss than the particles and decays exponentially: at 2.26e-3 s-1 with the default uptake coefficient 0.2 (`run_default.spec`, without the line; test 1), not at all with 0 (`run_zero.spec`, test 2), and at 6.59e-3 s-1 with 1 (`run_one.spec`, test 3), where gas-phase diffusion is half of the resistance (18% with 0.2). `compute_reference.py` writes the reference files.

The number of particles of a run is a Poisson sample (10000 on average), so the tolerance of the comparison is 5%.
