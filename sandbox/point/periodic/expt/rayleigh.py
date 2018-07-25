from numpy import *
import string as s
import _rayleigh
# from histo import histo
try:
	import macfs
	ismac = 1
except:
	ismac = 0

def rd_pcfdat (fname, data):
    """Read data from pcfreq data files."""
    
    infile = open(fname, 'r')
    lines = infile.readlines()[1:]
    for line in lines:
        data.append(eval(string.split(line)[0]))

data = []
rd_pcfdat('7565.dat', data)
data = N.array(data)
n = len(data)
print n, ' events'

# Compute the Rayleigh statistic.
ofname = '7565.ray'
ofile = open(ofname, 'w')
if ismac:
	ofileFS = macfs.FSSpec(ofname)
	ofileFS.SetCreatorType('R*ch','TEXT')

flo = 750.
fhi = 770.
nf = 200
df = (fhi-flo) / (nf-1.)
Rmax = 0.
for i in range(0, 2):
    f = flo + i*df
    w = 2*pi*f
    S = 0
    C = 0
#    wd = w*data
#    S = N.sum(N.sin(wd))
#    C = N.sum(N.cos(wd))
#    P = (S*S + C*C)/n
    S, C, P = raystat.raystat(w, data)
#    phi = atan(S/C)/(2*pi)
#    R = sqrt(P)
#    Rmax = max(Rmax, R)
#    line = 'i.f.R.P.phi %5i %#10.5g %#10.4g %#10.4g %#10.4g\n'
#    line = line % (i, f, R, P, phi)
    line = 'i.f.P %5i %#10.5g %#10.4g\n'
    line = line % (i, f, P)
    ofile.write(line)

ofile.write('\n')
#print 'Max R, P = ', Rmax, Rmax*Rmax

ofile.close()

wlo = 2*pi*flo
whi = 2*pi*fhi
#P = raystat.rsgrid(nf, wlo, whi, data)
#print P[0:5]

import cephes
def log_I0(x):
	return x + log(cephes.i0e(x))
k = .01
S = 1.e1
print 'cephes:  ', log_I0(k*S) - 200*log_I0(k)
print 'raystat: ', raystat.lml_wk(k,200,S)
