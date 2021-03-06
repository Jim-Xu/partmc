#!/usr/bin/env python
# Copyright (C) 2007-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import config
import config_filelist
import config_matplotlib
import os, sys, math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
sys.path.append("../../tool")
import partmc

fig_base_dir = "figs"
data_base_dir = "data"
data_type = "diam_bc_num"

value_min = 9.8 # 80
value_max = 1e6

x_axis = partmc.log_grid(min = config.diameter_axis_min,
                         max = config.diameter_axis_max,
                         n_bin = config.num_diameter_bins)
y_axis = partmc.linear_grid(min = config.bc_axis_min,
                            max = config.bc_axis_max,
                            n_bin = config.num_bc_bins)

def make_plot(value, out_filename):
    (figure, axes, cbar_axes) = config_matplotlib.make_fig(colorbar = True, right_margin = 0.9)

    axes.grid(True)
    axes.grid(True, which = 'minor')
    axes.minorticks_on()
    axes.set_xscale('log')

    axes.set_xbound(config.diameter_axis_min, config.diameter_axis_max)
    axes.set_ybound(config.bc_axis_min, config.bc_axis_max)

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    yaxis.set_minor_locator(matplotlib.ticker.MaxNLocator(8))

    axes.set_xlabel(r"dry diameter $D\ (\rm\mu m)$")
    axes.set_ylabel(r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$")

    axes.set_xbound(x_axis.min, x_axis.max)
    axes.set_ybound(y_axis.min, y_axis.max)
    
    mask = np.ma.make_mask(value <= 0.0)
    value_pos = np.ma.array(value, mask=mask)
    p = axes.pcolor(x_axis.edges(), y_axis.edges(), value_pos.transpose(),
                    norm = matplotlib.colors.LogNorm(vmin=value_min, vmax=value_max),
                    cmap=matplotlib.cm.jet, linewidths = 0.1)
    figure.colorbar(p, cax = cbar_axes, format = matplotlib.ticker.LogFormatterMathtext())
    cbar_axes.set_ylabel(r"number conc. $(\rm cm^{-3})$")
    figure.savefig(out_filename)

if __name__ == "__main__":
    for run in config_filelist.runs:
        data_dir = os.path.join(data_base_dir, run["name"])
        fig_dir = os.path.join(fig_base_dir, run["name"])
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)
        for loop in run["loops"]:
            for index in loop["indices"]:
                data_name = "%s_%04d_%08d" % (data_type, loop["num"], index["num"])
                print run["name"] + " " + data_name
                data_filename = os.path.join(data_dir, data_name + ".txt")
                value = np.loadtxt(data_filename)
                fig_filename = os.path.join(fig_dir, data_name + ".pdf")
                make_plot(value, fig_filename)
                plt.close('all')
