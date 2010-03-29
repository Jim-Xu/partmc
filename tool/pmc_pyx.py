#!/usr/bin/env python
# Copyright (C) 2007-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *

color_list = [color.hsb(0/3.0, 1, 1),
	      color.hsb(1/3.0, 1, 1),
	      color.hsb(2/3.0, 1, 1),
	      color.hsb(1/6.0, 1, 1),
	      color.hsb(3/6.0, 1, 1),
	      color.hsb(5/6.0, 1, 1),
	      color.hsb(1/12.0, 1, 1),
	      color.hsb(3/12.0, 1, 1),
	      color.hsb(5/12.0, 1, 1),
	      color.hsb(7/12.0, 1, 1),
	      color.hsb(9/12.0, 1, 1),
	      color.hsb(11/12.0, 1, 1),
	      ]

class listpalette(color.gradient):

    def __init__(self, colorlist):
        self.colorclass = colorlist[0][1].__class__
        self.colorlist = colorlist

    def getcolor(self, param):
	for i in range(len(self.colorlist)):
	    if self.colorlist[i][0] >= param:
		break
	else:
	    raise ValueError
	if i == 0:
	    i = 1
	# list[i-1] < param < list[i]
	alpha = (param - self.colorlist[i-1][0]) \
	    / (self.colorlist[i][0] - self.colorlist[i-1][0])
        colordict = {}
        for key in self.colorlist[0][1].color.keys():
            colordict[key] = alpha * self.colorlist[i][1].color[key] \
		+ (1 - alpha) * self.colorlist[i-1][1].color[key]
        return self.colorclass(**colordict)

rainbow_palette = listpalette([[0, color.rgb(0, 0, 1)],
			       [0.3, color.rgb(0, 1, 1)],
			       [0.5, color.rgb(0, 1, 0)],
			       [0.7, color.rgb(1, 1, 0)],
			       [1, color.rgb(1, 0, 0)]])

grid_painter = graph.axis.painter.regular(gridattrs = [style.linestyle.dotted])

aerosol_species_tex = {
    "SO4": "SO$_4$",
    "NO3": "NO$_3$",
    "Cl": "Cl",
    "NH4": "NH$_4$",
    "MSA": "MSA",
    "ARO1": "ARO1",
    "ARO2": "ARO2",
    "ALK1": "ALK1",
    "OLE1": "OLE1",
    "API1": "API1",
    "API2": "API2",
    "LIM1": "LIM1",
    "LIM2": "LIM2",
    "CO3": "CO$_3$",
    "Na": "Na",
    "Ca": "Ca",
    "OIN": "OIN",
    "OC": "Organic carbon",
    "BC": "Black carbon",
    "H2O": "H$_2$O",
    }

gas_species_tex = {
    "H2SO4": "H$_2$SO$_4$",
    "HNO3": "HNO$_3$",
    "HCl": "HC$\ell$",
    "NH3": "NH$_3$",
    "NO2": "NO$_2$",
    "NO3": "NO$_3$",
    "N2O5": "N$_2$O$_5$",
    "HNO4": "HNO$_4$",
    "O3": "O$_3$",
    "O1D": "O$_1$D",
    "O3P": "O$_3$P",
    "HO2": "HO$_2$",
    "H2O2": "H$_2$O$_2$",
    "SO2": "SO$_2$",
    "CH4": "CH$_4$",
    "C2H6": "C$_2$H$_6$",
    "CH3O2": "CH$_3$O$_2$",
    "CH3OH": "CH$_3$OH",
    "CH3OOH": "CH$_3$OOH",
    "C2O3": "C$_2$O$_3$",
    "CH3SO2H": "CH$_3$SO$_2$H",
    "CH3SCH2OO": "CH$_3$SCH$_2$OO",
    "CH3SO2": "CH$_3$SO$_2$",
    "CH3SO3": "CH$_3$SO$_3$",
    "CH3SO2OO": "CH$_3$SO$_2$OO",
    "CH3SO2CH2OO": "CH$_3$SO$_2$CH$_2$OO",
    }

def tex_species(species):
    if species in aerosol_species_tex.keys():
	return aerosol_species_tex[species]
    if species in gas_species_tex.keys():
	return gas_species_tex[species]
    return species
