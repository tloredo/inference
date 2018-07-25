"""
Use the Kaplan-Meier estimator to infer the luminosity function of
X-ray emitting stars in the open cluster NGC 2516, based on detections
and nondetections of X-ray counterparts to stellar optical counterparts
using the EPIC camera on XMM-Newton.

The data used are log10(Lx) for Lx in erg/s/cm2 (detections and upper limits).

Deep X-ray survey of the young open cluster NGC 2516 with XMM-Newton
I. Pillitteri1, G. Micela2, F. Damiani2, and S. Sciortino2
astro-ph/0601177, subm. to A&Ap
"""

import csv, re
from numpy import *
from pylab import *
from KM import KM

class Detection(object):

    # regexp to match 'val$_{minus}^{plus}$'
    Lx_re = re.compile(r'(?P<val>.*)\$_\{(?P<minus>.*)\}\^\{(?P<plus>.*)\}\$')

    # regexp to match 'val $\pm$ err'
    fx_re = re.compile(r'(?P<val>.*)\$\\pm\$(?P<err>.*)')

    def __init__(self, record):
        """
        Read and store detection data from a record read from the
        detection table.
        """
        self.x_id = int(record[0])
        self.opt_id = int(record[1])
        self.dist = float(record[2])
        self.Vmag = float(record[3])

        # Parse f_x \pm err in TeX
        m = self.fx_re.match(record[6]) 
        if m is None:
            raise ValueError, 'Bad f_x entry in table!'
        self.fx = float(m.group('val'))
        self.fx_err = float(m.group('err'))

        # Parse Log Lx -err +err in TeX
        m = self.Lx_re.match(record[7]) 
        if m is None:
            raise ValueError, 'Bad log(L_x) entry in table!'
        self.logLx = float(m.group('val'))
        self.logLx_merr = float(m.group('minus'))
        self.logLx_perr = float(m.group('plus'))

        
class Nondetection(object):

    def __init__(self, record):
        """
        Read and store nondetection data from a record read from the
        nondetection table.
        """
        self.opt_id = int(record[0])
        self.Vmag = float(record[1])
        self.logLx_lim = float(record[7])


# Gather detection data.
ifile = file('tab_b1.tex', 'r')
reader = csv.reader(ifile, delimiter='&')
dtxns = []
for record in reader:
    dtxn = Detection(record)
    dtxns.append(dtxn)
ifile.close()

# Gather nondetection data.
ifile = file('tab_c1.tex', 'r')
reader = csv.reader(ifile, delimiter='&')
nondtxns = []
# First row is special: has \leq:
record = reader.next()
# regexp for '$\leq$ lim'.
lim_re = re.compile(r'\s*\$\\leq\$(?P<lim>.*)')
m = lim_re.match(record[7])
record[7] = m.group('lim')
nondtxns = [ Nondetection(record) ]
# Get the data in the rest of the rows.
for record in reader:
    nondtxn = Nondetection(record)
    nondtxns.append(nondtxn)
ifile.close()


# Make arrays just of the log(L) estimates and limits.
logLx_vals = array([ d.logLx for d in dtxns ])
logLx_lims = array([ n.logLx_lim for n in nondtxns ])

# Find the KM estimate.
kmest = KM(logLx_vals, logLx_lims)
