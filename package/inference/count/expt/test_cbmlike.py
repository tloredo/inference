from _cbmlike import *
from numpy import *
from scipy.special import hyp1f1
from scipy.stats import poisson
from ioutils import cprint
from numutils import lin_stepper


n_on, t_on = 6, 1.
n_off, t_off = 5, 1.
cut = 1.e-6

c = sigctp(n_on, t_on, n_off, t_off)
offset = slml_offset(n_on, t_on, n_off, t_off)

print 'c length:', len(c)
print c
print 'c sum =', sum(c)
print 'offset =', offset

# This used to be exposed from _cbmlike; use scipy's instead.
def ppois(i, mu):
    return poisson.pmf(i, mu)

def psig2(s,t_on, c):
    ps = 0.
    st = s*t_on
    for i in range(len(c)):
        ps += c[i] * ppois(i, st)
    return t_on*ps

print 'Comparing C coef and recursive algorithms:'
print 's   coeff  exp(recurs)  ratio   nterms  coeff/scipy'
for s in [.5, 1., 1.5, 3.]:
    ps = psigc(s, t_on, c)
    ps2 = psig2(s, t_on, c)
    slml, nt, err = slmlike(s, n_on, t_on, n_off, t_off, offset, cut)
    ls = exp(slml)
    cprint(s, ps, ls, ls/ps, nt, ps/ps2)

print

print 's  coeff  exp(recurs)   ratio'
for s in lin_stepper(0, 5, 11):
    ps = psigc(s, t_on, c)
    slml, nt, err = slmlike(s, n_on, t_on, n_off, t_off, offset, cut)
    ls = exp(slml)
    cprint(s, ps, ls, ls/ps)

print
print 'Vector test:'

t_off = 1.
t_on = 2.
n_on = array([1, 20, 100])
n_off = array([1, 5, 20])
s = array([0., 5., 30.])
sh = 0.5*s
s2 = 2*s

print 's near peak:', s
offset = vslml_offset(n_on, t_on, n_off, t_off)
slml = vslmlike(sh, n_on, t_on, n_off, t_off, offset, cut)
print 'half peak:', slml
slml = vslmlike(s, n_on, t_on, n_off, t_off, offset, cut)
print 'near peak:', slml
slml = vslmlike(s2, n_on, t_on, n_off, t_off, offset, cut)
print 'dbl peak: ', slml
s = s+5
slml = vslmlike(s, n_on, t_on, n_off, t_off, offset, cut)
print 'shift 5:  ', slml
