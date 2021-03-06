#!/usr/bin/env python

import os, sys
import numpy as np

sys.path.append("../../tool")
import mpl_helper
import matplotlib
import config

num_conc_min = 4e0
num_conc_max = 4e5
mass_conc_min = 1e-5
mass_conc_max = 4e1

x_array = np.loadtxt("data/2d_bc_10K_wei+1_12_x_values.txt") * 1e6
y_array = np.loadtxt("data/2d_bc_10K_wei+1_12_y_values.txt") * 100

line1 = np.zeros([len(x_array)])
line1[:] = 70

def min_pos(a):
    ma = np.ma.masked_less_equal(a, 0)
    return ma.min()

num_avg1 = np.loadtxt("data/2d_bc_10K_wei+1_12_hist_average_num.txt") / 1e6
num_avg2 = np.loadtxt("data/2d_bc_10K_wei-1_12_hist_average_num.txt") / 1e6
num_avg3 = np.loadtxt("data/2d_bc_10K_wei-4_12_hist_average_num.txt") / 1e6
print "num range: ", min_pos(num_avg1), num_avg1.max()
print "num range: ", min_pos(num_avg2), num_avg2.max()
print "num range: ", min_pos(num_avg3), num_avg3.max()

mass_avg1 = np.loadtxt("data/2d_bc_10K_wei+1_12_hist_average_mass.txt") * 1e9
mass_avg2 = np.loadtxt("data/2d_bc_10K_wei-1_12_hist_average_mass.txt") * 1e9
mass_avg3 = np.loadtxt("data/2d_bc_10K_wei-4_12_hist_average_mass.txt") * 1e9
print "mass range: ", min_pos(mass_avg1), mass_avg1.max()
print "mass range: ", min_pos(mass_avg2), mass_avg2.max()
print "mass range: ", min_pos(mass_avg3), mass_avg3.max()

(figure, axes_array, cbar_axes_array) \
    = mpl_helper.make_fig_array(3,2, figure_width=config.figure_width_double, 
                                top_margin=1, bottom_margin=0.45,
                                left_margin=1.07, right_margin=0.65,
                                vert_sep=0.3, horiz_sep=0.3,
                                colorbar="shared", colorbar_location="top")

axes = axes_array[2][0]
cbar_axes = cbar_axes_array[0]
p = axes.pcolor(x_array, y_array, num_avg1.transpose(), 
                norm = matplotlib.colors.LogNorm(vmin=num_conc_min, vmax=num_conc_max),
                linewidths=0.1, edgecolors='face')
axes.plot(x_array,line1, 'k')
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"BC mass frac. $w_{\rm BC}$ / \%")
axes.set_ylim(0, 80)
axes.set_xlim(5e-3, 5)
axes.grid(True)
axes.text(-0.3, 0.5, r'$\alpha = +1$', horizontalalignment='center',
           verticalalignment='center', transform=axes.transAxes,
           rotation='vertical', bbox=dict(edgecolor='black', facecolor='white',
                                          boxstyle="round,pad=0.5"))
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"number conc. $n(D_{\rm dry},w_{\rm BC})$ / $\rm cm^{-3}$")

axes = axes_array[1][0]
axes.pcolor(x_array, y_array, num_avg2.transpose(),
            norm = matplotlib.colors.LogNorm(vmin=num_conc_min, vmax=num_conc_max),
            linewidths=0.1, edgecolors='face')
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"BC mass frac. $w_{\rm BC}$ / \%")
axes.set_ylim(0, 80)
axes.set_xlim(5e-3, 5)
axes.grid(True)
axes.text(-0.3, 0.5, r'$\alpha = -1$', horizontalalignment='center',
           verticalalignment='center', transform=axes.transAxes,
           rotation='vertical', bbox=dict(edgecolor='black', facecolor='white',
                                          boxstyle="round,pad=0.5"))

axes = axes_array[0][0]
axes.pcolor(x_array, y_array, num_avg3.transpose(),
            norm = matplotlib.colors.LogNorm(vmin=num_conc_min, vmax=num_conc_max),
            linewidths=0.1, edgecolors='face')
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"BC mass frac. $w_{\rm BC}$ / \%")
axes.set_ylim(0, 80)
axes.set_xlim(5e-3, 5)
axes.grid(True)
axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \upmu m$")
axes.text(-0.3, 0.5, r'$\alpha = -4$', horizontalalignment='center',
           verticalalignment='center', transform=axes.transAxes,
           rotation='vertical', bbox=dict(edgecolor='black', facecolor='white',
                                          boxstyle="round,pad=0.5"))

axes = axes_array[2][1]
cbar_axes = cbar_axes_array[1]
p = axes.pcolor(x_array, y_array, mass_avg1.transpose(),
                norm = matplotlib.colors.LogNorm(vmin=mass_conc_min, vmax=mass_conc_max),
                linewidths=0.1, edgecolors='face')
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylim(0, 80)
axes.set_xlim(5e-3, 5)
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar_axes.set_xlabel(r"mass conc. $m^{\rm BC}(D_{\rm dry},w_{\rm BC})$ / ($\rm \upmu g \ m^{-3}$)")
cbar.set_ticks([1e-5, 1e-3, 1e-1, 1e1])

axes = axes_array[1][1]
axes.pcolor(x_array, y_array, mass_avg2.transpose(),
            norm = matplotlib.colors.LogNorm(vmin=mass_conc_min, vmax=mass_conc_max),
            linewidths=0.1, edgecolors='face')
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylim(0, 80)
axes.set_xlim(5e-3, 5)
axes.grid(True)

axes = axes_array[0][1]
axes.pcolor(x_array, y_array, mass_avg3.transpose(),
            norm = matplotlib.colors.LogNorm(vmin=mass_conc_min, vmax=mass_conc_max),
            linewidths=0.1, edgecolors='face')
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylim(0, 80)
axes.set_xlim(5e-3, 5)
axes.grid(True)
axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \upmu m$")

mpl_helper.remove_fig_array_axes(axes_array)

figure.savefig("figs/num_mass_bc_diameter_2d.pdf")

