from numpy import pi, sqrt, ones
from numpy import array, sin, cos, fromfunction
from inference.gauss import vecba

#   data = array([ [0,0], [1,0], [2,0], [3,0], [4,0] ], float)
data = fromfunction(lambda i,j: 1.*i, (10000,2))
data = .05*data
data[:,1] = 2. + cos(0.5*pi*data[:,0])
print "# of data: ", len(data)
#   print data

def set_wt (w, times): return(w*times)
def cos_wt (w, times, wt): return cos(wt)
def sin_wt (w, times, wt): return sin(wt)
def const (w, times, wt): return ones(len(times), float)

print "Running test case of vecba:"
# import profile
# profile.run("test = vecba( (cos_wt, sin_wt, const), data, sigma=1., setup=set_wt )")
test = vecba.BA( (cos_wt, sin_wt, const), data, sigma=1., setup=set_wt )
print "Marg: ", test.marg_stats(0.53*pi)
print "Amps: ", test.amplitudes()
R = test.residuals()
rms = sqrt(sum(R*R))
print "RMS residual: ", rms
#   print "Resids: ", test.residuals()

def LoopMarg():
    for i in range(1000):
        if i%2 == 0:
            w = .52*pi
        else:
            w = .53*pi
        t = test.marg_stats(w)
        if (i+1)%50 == 0: print i+1

#LoopMarg()

#import leaktest
#leaktest.eval_test("LoopMarg()",globals(),locals())

class DataBasis:
    """
    Data and basis functions for constant plus sinusoid.
    """
    
    def __init__(self, data, sigmas):
        self.times = data[:,0].copy()
        self.smpls = data[:,1].copy()
        self.sigs = sigmas  # scalar or array matching samples
        self.std_smpls = self.smpls/self.sigs
        self.n_data = len(self.times)
        self.n_basis = 3
        self.ones = ones(self.n_data, float)
        self.std_basis = (self.cos_wt, self.sin_wt, self.const)
    
    def set_nonlin(self, w):
        self.w = w
        self.wtimes = w*self.times
    
    def cos_wt(self):
        return cos(self.wtimes)/self.sigs
    
    def sin_wt(self):
        return sin(self.wtimes)/self.sigs
    
    def const(self):
        return self.ones


db = DataBasis(data, 1.)
test2 = vecba.BAObj(db)
print "Marg: ", test2.marg_stats(0.53*pi)
print "Amps: ", test2.amplitudes()
R = test2.residuals()
rms = sqrt(sum(R*R))
print "RMS residual: ", rms
