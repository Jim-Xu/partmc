# Residual layer and gas dry deposition test case

This tests the optional scenario lines `gas_deposition` (first-order loss of the gases at the rate v_d / h) and `residual_layer` (exchange of the mixed layer with a residual layer above it when the mixing height changes).

The mixing height is 1000 m, falls to 200 m between 1 and 2 h (evening), stays at 200 m until 4 h and grows to 1500 m between 4 and 6 h (morning). There is no chemistry and no coagulation, and in tests 1 and 2 no exchange with the background other than by entrainment, so every quantity follows from mixing, emission and deposition:

- gases: NO and HCl are in the initial air only, NO2 only in the background, NH3 is emitted during the night (2-4 h), and HNO3 is in the initial air and deposits while the mixing height is constant (velocities and scaling in `gas_depo.dat`);
- particles: three mono-disperse single-species populations, SI in the initial air, SB in the background, SE emitted during the night.

With the residual layer (`run_part_rl.spec`, test 1) the air that leaves the mixed layer in the evening keeps its composition and is entrained again in the morning, before any background air above the old mixing height of 1000 m: the morning mixed layer has (200 m of night air + 800 m of residual-layer air + 500 m of background air) / 1500 m. Without it (`run_part_norl.spec`, test 2) the growing mixed layer entrains only background air. Test 3 (`run_part_rldil.spec`) repeats test 1 with exchange with the background air at 1/(4 h) (`gas_back_dil.dat`, `aero_back_dil.dat`): the residual layer mixes with the background air at the same rate as the mixed layer, so that the residual-layer air entrained in the morning has partly been replaced by background air. `compute_reference.py` computes the three references with the time stepping of PartMC (profiles evaluated at the end of each 60 s step); the gases must agree to rounding, the particle concentrations within the Monte Carlo error.
