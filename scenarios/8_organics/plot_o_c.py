#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import netCDF4, numpy

(figure, axes) = mpl_helper.make_fig(top_margin=0.6, right_margin=0.8)
axes2 = axes.twinx()

ncf = netCDF4.Dataset("out/organics_process.nc")
time = numpy.array(ncf.variables["time"][:]) / 3600
tot_o_c = numpy.array(ncf.variables["tot_o_c"][:])
tot_o_c_err = numpy.array(ncf.variables["tot_o_c_ci_offset"][:])
sec_o_c = numpy.array(ncf.variables["sec_o_c"][:])
sec_o_c_err = numpy.array(ncf.variables["sec_o_c_ci_offset"][:])
tot_om_oc = numpy.array(ncf.variables["tot_om_oc"][:])
tot_om_oc_err = numpy.array(ncf.variables["tot_om_oc_ci_offset"][:])

# there is no secondary organic matter at the start
sec = (sec_o_c > 0)

axes.errorbar(time, tot_o_c, tot_o_c_err, fmt="b-", label=r'O:C, all organics')
axes.errorbar(time[sec], sec_o_c[sec], sec_o_c_err[sec], fmt="g-",
              label=r'O:C, secondary')
axes2.errorbar(time, tot_om_oc, tot_om_oc_err, fmt="r-", label=r'OM:OC')
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"O:C")
axes2.set_ylabel(r"OM:OC")
axes.set_ylim(bottom=0)
axes.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.35, 1.0))
axes2.legend(loc='lower center', bbox_to_anchor=(0.9, 1.0))
axes.grid(True)

out_filename = "out/organics_o_c.pdf"
print("Writing %s" % out_filename)
figure.savefig(out_filename)
