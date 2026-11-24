
# Organic aerosol scenario

This scenario is the test bed for the organic aerosol scheme of MOSAIC: every formation pathway of the scheme has emissions, and the processing program writes the organic aerosol diagnostics. It is independent of the other scenarios.

The case that is set up here is an urban plume. It takes the meteorology, the initial and background conditions, the emissions of inorganic gases and CBM-Z organic gases and the particle size distributions from the `1_urban_plume` scenario and adds the emissions of the precursors of the organic aerosol. Other cases (chamber experiments, biogenic or biomass-burning environments) can be set up in the same way by replacing the input files.

The organic aerosol scheme follows the `cb6r5_ae7_aq` mechanism (AERO7) of CMAQ: SOA from aromatics, long alkanes, monoterpenes (including organic nitrates), isoprene (including the reactive uptake of the epoxydiol and of methylglyoxal) and sesquiterpenes, particle-phase oligomerization, semivolatile primary organic aerosol (POA) with OH aging, and potential SOA from combustion emissions. Unlike in CMAQ, the potential combustion SOA has a fresh class `APCSO` (O:C 0.51) and an aged class `APCSO2` (O:C 1.28); the fresh class turns into the aged one by OH (1.03e-11 cm³ s⁻¹ × [OH], mass conserving), so that the O:C of urban SOA rises with photochemical age as in the SIMPLE parameterization of Hodzic and Jimenez (2011), with the calibration factor of Hayes et al. (2015). The monoterpene nitrates are also split, into a stable species (`MTNO3`, `AMTNO3`) and a hydrolyzable one (`MTNO3H`, `AMTNO3H`, 12% of the nitrate formed): CMAQ hydrolyzes all particulate monoterpene nitrate with a lifetime of 3 h, but only about a tenth of the nitrates from the oxidation by the nitrate radical hydrolyzes (Boyd et al., 2015; Takeuchi and Ng, 2019). The monoterpene SOA of low volatility (`AMT1`, `AMT2`) photolyzes in the particles as the monoterpene SOA of the CRACMM3 mechanism of CMAQ 6.0: at 1% of the photolysis frequency of NO2 it turns into the more oxidized `AMTN1` (O:C 1.0; 84% of the mass) and formic acid, so that the O:C of the monoterpene SOA rises from that of fresh chamber SOA (about 0.4) during the day. The O:C and OM:OC of the organic species in `aero_data.dat` are those of CMAQ, where the O:C of the species without a surrogate molecule follows from OM:OC with a relation for molecules of carbon, hydrogen and oxygen (Simon and Bhave, 2012; Pye et al., 2017). For the species with nitrogen or sulfur that relation counts the nitrogen or sulfur as oxygen, so their O:C is set from the molecular composition instead: 0.50 for the monoterpene nitrates (C10H17NO5; CMAQ 0.587) and 1.04 for the isoprene epoxydiol SOA `AISO3` (60% 2-methyltetrols with O:C 0.8 and 40% organosulfate with O:C 1.4, the two species into which CMAQ itself splits it; CMAQ 1.307). These are molecular ratios: an aerosol mass spectrometer assigns the oxygen of the NO2 and SO3 groups to inorganic nitrate and sulfate (Farmer et al., 2010) and reports about 0.30 and 0.80 for these species, which matters when the O:C of a run is compared with such measurements.

In the gas phase, the reactions of the monoterpenes (`API`, `LIM`, `BPIN`) with OH, O3 and NO3 consume the oxidant and form the gas-phase products of the monoterpene reactions of CB6r5 (TERP and APIN in `cb6r5_ae7_aq`), written with the species of CBM-Z. In the MOSAIC library of 2012 these were tracer reactions: OH was not consumed and no radicals were formed, so that OH was too high wherever monoterpenes are an important sink of it. One gas species was added for this, `XO2N`, the organic peroxy radical that forms an organic nitrate with NO instead of NO2. A `gas_data.dat` without `XO2N` still works, but this radical is then lost at the end of every time step.

The isoprene peroxy radical (`ISOPP`) also follows CB6r5: its reaction with HO2 forms the hydroxy hydroperoxides `ISPX` with 0.12 OH and HO2 (in CBM-Z it formed the generic hydroperoxide `ROOH`), and it isomerizes to HO2 and the hydroperoxy aldehydes, the gas species `HPLD`, which photolyze to OH at 0.055 times the photolysis frequency of NO2. The reaction of `ISPX` with OH, which forms the epoxydiols `EPOX`, consumes OH and returns 0.933 OH. Under low NO these reactions return most of the radicals that the oxidation of isoprene consumed in CBM-Z. The three isoprene peroxy radicals (`ISOPP`, `ISOPN`, `ISOPO2`) react with NO at 2.7e-12 exp(350/T) cm3 s-1 and `ISOPP` reacts with HO2 at 2.12e-13 exp(1300/T) cm3 s-1, the measured rate constants (Wennberg et al., Chem. Rev. 118, 3337, 2018), instead of the CBM-Z values (4.0e-12 and 1.7e-13 exp(1300/T)); together with the isomerization these set how much of the isoprene ends as `ISPX` and `EPOX`. A `gas_data.dat` needs `ISPX` and `HPLD`: without them the products of `ISOPP` + HO2 and of the isomerization are lost at the end of every time step.

Three photolysis frequencies also follow CB6r5 (ratios to that of NO2, from the cross sections and quantum yields of `cb6r5_ae7_aq`): methylglyoxal `MGLY` 1.41e-2, the organic nitrates `ONIT` 3.2e-4 and the isoprene products `ISOPRD` (methacrolein, methyl vinyl ketone) 2.87e-3. The CBM-Z values were 2.0, 0.31 and 0.05 times these. The value for `ISOPRD` is the one of the CMAQ table, which has a quantum yield of 0.05 at all wavelengths; the quantum yields of methacrolein and methyl vinyl ketone are small and uncertain, and this radical source with them. Where isoprene dominates, the faster photolysis of `ISOPRD` adds about a tenth to the primary source of radicals; in an urban plume the slower photolysis of `MGLY` lowers HO2 by about a tenth and the ozone maximum by a few percent. The organic aerosol does not change beyond the noise of the particle sampling, except for the SOA from methylglyoxal, which follows its gas concentration. The photolysis of `ISOPRD` also has the products of CB6r5 (R161), written with the species of CBM-Z: 1.10 HO2 + 0.50 XO2 + 0.34 CH3O2 + 0.208 C2O3 + 0.26 HCHO + 0.24 OLET + 0.24 PAR + 0.17 AONE + 0.128 ALD2, where ALD2 stands for the glycolaldehyde of CB6r5; CBM-Z formed 0.97 C2O3 and 0.7 CH3O2, that is more peroxyacetyl nitrate and less HO2.

Three further reactions of the radicals follow CMAQ. The peroxyacetyl radical reacts with HO2 as in CB6r5 (R57: 3.14e-12 exp(580/T) cm3 s-1, 0.5 OH + 0.5 CH3O2 + 0.13 O3 + 0.5 `RCOOH`); in CBM-Z this reaction ended the radical chain (0.4 O3 + 0.4 `RCOOH`). The reaction of the isoprene peroxy radical `ISOPP` with the other organic peroxy radicals forms 0.8 HO2, 0.598 HCHO and 0.072 `XO2` besides `ISOPRD` (CB6r5 R154); in CBM-Z it formed `ISOPRD` only. O(1D) is quenched by O2 and N2 at 3.3e-11 exp(55/T) and 2.15e-11 exp(110/T) cm3 s-1 (the values of the CRACMM3 mechanism of CMAQ 6.0, whose sum for air is the O1D + M of CB6r5) and reacts with water at 2.14e-10 cm3 s-1 (CB6r5 R11), so that an O(1D) atom forms 12% less OH than with the rate constants of CBM-Z. Together these lower HO2 by 10 to 15% and the ozone maximum by 4% in the urban plume of this scenario; where isoprene dominates, OH is 3% lower in the afternoon.

## Emissions

The emission files `gas_emit.dat` and `aero_emit_comp_*.dat` are generated by `make_emissions.py` from `gas_emit_base.dat` (the gas emissions of `1_urban_plume`) and the particle emissions in `aero_emit_dist.dat`. All scaling factors are at the top of that script. To change them, edit the script and run `./make_emissions.py`.

| Species | Emission | Source of the scaling |
|---|---|---|
| `SOATOL` | 1 × `TOL` | Shadow precursor of the toluene SOA, reacts at the CB6r5 rate without changing the CBM-Z chemistry. |
| `SOAXYL` | 0.998 × `XYL` | Shadow precursor of the xylene SOA; xylene without naphthalene, see `NAPH`. |
| `NAPH` | 0.002 × `XYL` | CMAQv5.1 VOC emission updates (Pouliot and Pye), estimate for inventories without naphthalene, based on the 2011 NEI. |
| `SOAALK` | 0.108 × `PAR` | Same source, estimate for inventories without the long-chain alkane tracer. |
| `BENZ` | 1.3e-3 × `CO` | Urban emission ratio of benzene to CO in Los Angeles, 1.30 pptv/ppbv (Borbon et al., 2013, Table 2; Paris: 1.07 at the ground, 2.50 airborne). |
| `API` | 0.0618 × `ISOP` | Global annual emissions of MEGAN2.1 (Guenther et al., 2012, Table 6): α-pinene 66.1 Tg, isoprene 535 Tg. |
| `BPIN` | 0.0214 × `ISOP` | Same table: β-pinene 18.9 Tg and camphene 4.0 Tg. β-pinene has its own species because it reacts with ozone at a fifth of the rate of α-pinene (IUPAC rate constants), and its reaction with the nitrate radical forms the monoterpene nitrates; camphene, which reacts even more slowly with ozone, is emitted with it. |
| `LIM` | 0.0685 × `ISOP` | Same table: all other monoterpenes, 73.3 Tg. |
| `SESQ` | 0.0181 × `ISOP` | Same table: all sesquiterpenes, 29.0 Tg. |
| `PCVOC` | 6.579 × POA (mass) | CMAQ v5.5 emission control file for `cb6r5_ae7_aq`; Murphy et al. (2017). |
| POA in particles | `ALVPO1` 0.09, `ASVPO1` 0.09, `ASVPO2` 0.14, `ASVPO3` 0.18, `AIVPO1` 0.50 (mass fractions) | Same sources. |

All ratios are molar unless noted otherwise.

Things to be aware of:

* The primary organic matter of the cooking, diesel and gasoline particles is emitted entirely in the particle phase, with the volatility distribution above. It then evaporates according to the gas-particle partitioning, so most of it is in the gas phase shortly after emission. CMAQ instead emits the more volatile bins directly as vapors; the total is the same.
* The size distributions and the organic mass fractions of the emitted particles are those of `1_urban_plume`. The POA species have a density of 1400 kg/m^3 and `OC` has 1000 kg/m^3, so the emitted organic mass is larger than in `1_urban_plume` (by 40% for the cooking particles, less for the others).
* The initial and background particles contain the nonvolatile `OC`, as in `1_urban_plume`.
* The terpene emissions follow the diurnal profile of isoprene. Real monoterpene emissions continue at night, so the nighttime chemistry with the nitrate radical is weaker here than in a forested environment.
* The global ratios of terpenes to isoprene are not specific to any region.

## Running

    ./1_run.sh
    ./2_process.sh
    ./3_plot.sh

`organics_process.F90` writes `out/organics_process.nc` with the time series (averaged over the repeats) of

* `mass_conc_<group>`: mass concentration of the species groups `poa`, `oxidized_poa`, `soa_combustion`, `soa_anthro`, `soa_monoterpene`, `soa_isoprene`, `soa_sesquiterpene`, `soa_bio_oligomer`, `soa_glyoxal`,
* `tot_org_mass_conc`, `tot_sec_org_mass_conc`: all organic matter and secondary organic matter (everything except the group `poa`),
* `tot_o_c`, `tot_om_oc`: carbon-weighted O:C and OM:OC of all organic matter; `sec_o_c`: O:C of the secondary organic matter,

and one file per output time with the two-dimensional number distributions of dry diameter against the per-particle O:C (`diam_o_c_dist`) and against the secondary organic fraction of the dry mass (`diam_sec_frac_dist`).

## References

A. Borbon et al. (2013) Emission ratios of anthropogenic volatile organic compounds in northern mid-latitude megacities: Observations versus emission inventories in Los Angeles and Paris, _J. Geophys. Res. Atmos._ 118, 2041-2057, DOI: <http://dx.doi.org/10.1002/jgrd.50059>.

D. K. Farmer et al. (2010) Response of an aerosol mass spectrometer to organonitrates and organosulfates and implications for atmospheric chemistry, _Proc. Natl. Acad. Sci. USA_ 107, 6670-6675, DOI: <http://dx.doi.org/10.1073/pnas.0912340107>.

A. B. Guenther et al. (2012) The Model of Emissions of Gases and Aerosols from Nature version 2.1 (MEGAN2.1): an extended and updated framework for modeling biogenic emissions, _Geosci. Model Dev._ 5, 1471-1492, DOI: <http://dx.doi.org/10.5194/gmd-5-1471-2012>.

P. L. Hayes et al. (2015) Modeling the formation and aging of secondary organic aerosols in Los Angeles during CalNex 2010, _Atmos. Chem. Phys._ 15, 5773-5801, DOI: <http://dx.doi.org/10.5194/acp-15-5773-2015>.

A. Hodzic and J. L. Jimenez (2011) Modeling anthropogenically controlled secondary organic aerosols in a megacity: a simplified framework for global and climate models, _Geosci. Model Dev._ 4, 901-917, DOI: <http://dx.doi.org/10.5194/gmd-4-901-2011>.

B. N. Murphy et al. (2017) Semivolatile POA and parameterized total combustion SOA in CMAQv5.2: impacts on source strength and partitioning, _Atmos. Chem. Phys._ 17, 11107-11133, DOI: <http://dx.doi.org/10.5194/acp-17-11107-2017>.

G. Pouliot and H. O. T. Pye, CMAQv5.1 VOC emission updates, CMAQ wiki, <https://cmascenter.org/cmaq/wiki/airqualitymodeling.org/index.php/CMAQv5.1_VOC_Emission_Updates>.

H. O. T. Pye et al. (2017) On the implications of aerosol liquid water and phase separation for organic aerosol mass, _Atmos. Chem. Phys._ 17, 343-369, DOI: <http://dx.doi.org/10.5194/acp-17-343-2017>.

N. Riemer, M. West, R. A. Zaveri, and R. C. Easter (2009) Simulating the evolution of soot mixing state with a particle-resolved aerosol model, _J. Geophys. Res._ 114(D09202), DOI: <http://dx.doi.org/10.1029/2008JD011073>.

H. Simon and P. V. Bhave (2012) Simulating the degree of oxidation in atmospheric organic particles, _Environ. Sci. Technol._ 46, 331-339, DOI: <http://dx.doi.org/10.1021/es202361w>.
