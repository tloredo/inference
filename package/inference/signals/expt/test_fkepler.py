from numpy import *
from inference.signals import fkepler

# t2EA, t2TA, vt2EA, vt2TA already tested
# Check consistency of "store" versions with these.

T_max = 100.
nt = 50
epochs = T_max*random.random(nt)  # random epochs in [0,T_max)
epochs.sort()

tau = 31.
e = .7
Mp = .3*pi

fkepler.settol(1.e-4)
fkepler.setup_Mp(tau, e, Mp)

cs_EA1 = fkepler.vt2EA(epochs)
cs_EA2 = zeros((nt,2), float)
fkepler.t2EA_store(epochs, cs_EA2)
deltas = cs_EA2 - cs_EA1
print 'min, max EA deltas:', deltas.min(0), deltas.max(0)

cs_TA1 = fkepler.vt2TA(epochs)
cs_TA2 = zeros((nt,2), float)
fkepler.t2TA_store(epochs, cs_TA2)
deltas = cs_TA2 - cs_TA1
print 'min, max TA deltas:', deltas.min(0), deltas.max(0)

cs_EA3 = zeros((nt,2), float)
cs_TA3 = zeros((nt,2), float)
fkepler.t2anoms(epochs, cs_EA3, cs_TA3)
dea = cs_EA3 - cs_EA1
dta = cs_TA3 - cs_TA1
print 'min, max EA deltas:', dea.min(0), dea.max(0)
print 'min, max TA deltas:', dta.min(0), dta.max(0)
