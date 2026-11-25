#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import netCDF4, numpy

# (variable name in the process file, legend label)
groups = [
    ("poa", "POA (alkane-like ROC, OC)"),
    ("soa_roc", "oxygenated ROC SOA"),
    ("soa_monoterpene", "monoterpene SOA"),
    ("soa_isoprene", "isoprene SOA"),
    ("soa_glyoxal", "methylglyoxal SOA"),
    ("soa_cloud", "AORGC (cloud SOA; here regional OA)"),
]

(figure, axes) = mpl_helper.make_fig(figure_width=7, right_margin=2.6)

ncf = netCDF4.Dataset("out/organics_process.nc")
time = numpy.array(ncf.variables["time"][:]) / 3600
mass_concs = [numpy.array(ncf.variables["mass_conc_" + name][:]) * 1e9
              for (name, label) in groups]

axes.stackplot(time, mass_concs, labels=[label for (name, label) in groups])
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"organic mass conc. / $\rm \mu g\ m^{-3}$")
axes.set_xlim(time.min(), time.max())
handles, labels = axes.get_legend_handles_labels()
axes.legend(handles[::-1], labels[::-1], loc='center left',
            bbox_to_anchor=(1.0, 0.5))
axes.grid(True)

out_filename = "out/organics_organic_mass.pdf"
print("Writing %s" % out_filename)
figure.savefig(out_filename)
