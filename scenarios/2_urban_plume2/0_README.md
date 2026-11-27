
# Urban plume 2 scenario

R. A. Zaveri, J. C. Barnard, R. C. Easter, N. Riemer, and M. West, (2010) Particle-resolved simulation of aerosol size, composition, mixing state, and the associated optical and cloud condensation nuclei activation properties in an evolving urban plume, _J. Geophys. Res._ 115(D17210), DOI: <http://dx.doi.org/10.1029/2009JD013616>.

This scenario is a slight tweak of `1_urban_plume`. It scales back the NOx emissions, extends the simulation duration to 48 hours, and adds black carbon aerosols to the background.

The gas and aerosol species are those of MOSAIC with the organic aerosol of the CRACMM3 mechanism of CMAQ (see `scenarios/8_organics`). The organic matter of the initial and background particles, the SORGAM species `API1` in the paper, is the aged regional organic matter `AORGC`. As in `1_urban_plume`, the primary organic aerosol is emitted as the nonvolatile species `OC`, and the SOA of the aromatics forms from the tracer species `SOATOL` and `SOAXYL`, which are emitted with toluene and the xylenes.
