#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import netCDF4, numpy

for (filename, index) in mpl_helper.get_filename_list('out/', r'organics_([0-9]+)_process\.nc'):
    (figure, axes, cbar_axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1, colorbar=True)

    ncf = netCDF4.Dataset(filename)
    diam_edges = numpy.array(ncf.variables["diam_edges"][:]) * 1e6
    o_c_edges = numpy.array(ncf.variables["o_c_edges"][:])
    diam_o_c_dist = numpy.array(ncf.variables["diam_o_c_dist"][:]) * 1e-6

    p = axes.pcolor(diam_edges, o_c_edges, diam_o_c_dist,
                    norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e6), linewidths = 0.1)

    axes.set_xscale("log")
    axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \mu m$")
    axes.set_xlim(1e-2, 1e0)

    axes.set_yscale("linear")
    axes.set_ylabel(r"O:C of the organic matter")
    axes.set_ylim(0, 1.2)

    axes.grid(True)
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                           orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number conc. $n(D_{\rm dry},{\rm O{:}C})$ / $\rm cm^{-3}$")

    out_filename = "out/organics_diam_o_c_dist_%s.pdf" % index
    print("Writing %s" % out_filename)
    figure.savefig(out_filename)
    matplotlib.pyplot.close(figure)
