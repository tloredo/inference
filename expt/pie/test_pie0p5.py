"""
2021-05-31 Modified from pie_test0p4 for PIE v0.5, updated for Py-3
"""

from matplotlib.pyplot import *
from numpy import *
from inference.pie import *
from inference.pie.inference_base import Inference
from inference.pie.signalmodel import StopStepping
from inference.pie.logger import log_debug


ion()


################################################################################
#  Testing, 1 2 3...
#

# log_debug()

class Test(SignalModel):
    # def __init__(self):
    #    ParameterizedModel.__init__(self)

    x = RealParam(1.2)
    y = RealParam(3.4, '2nd param')


a = array([1.])


class EmptyPreds(PredictorSet):
    d1 = SampledGaussian(a, a, a)


class Mix(ChisqrInference, Test, EmptyPreds):
    pass


if 0:
    fph = RealParamHandler()
    f = RealParamValue(1.23, fph)
    print(f, 2*f, f-1)
    print(f.addOne())

a = Mix()
print(a.x, a.y)
b = Mix()
b.x = 5.6
b.y = 7.8
print('a:', a.x, a.y)
print('b:', b.x, b.y)
# print 'a dict:', a.__dict__
# print 'b dict:', b.__dict__
# print 'Test dict:', Test.__dict__
# print 'Mix dict:', Mix.__dict__

# Manually step through a grid to verify step order.
if 0:
    a.x.step(1,3,5)
    a.y.logstep(1,1000,4)

    def do_grid():
        while True:
            try:
                a.next_step()
                print(a.step_nums(), a.x, a.y)
            except StopStepping:
                print('Completed the grid!')
                break
    do_grid()
    # a.x.vary(1., .1, (1,5))

print('Testing with signal(x) = p*x + q...')
from numpy import array
locns1 = array([1,2,3],float)
vals1 = array([3,5,7],float)  # for p=2, q=1
sigs1 = array([1,116],float)  # bad dimen, to test checking
locns2 = array([3,4,5],float)
vals2 = array([7,9,11],float)
sigs2 = array([2,2,2],float)


class MyModel(SignalModel):
    # def __init__(self):
    #    ParameterizedModel.__init__(self)

    p = RealParam(1.2)
    q = RealParam(3.4, '2nd param')

    def signal(self, x):
        return self.p*x + self.q

    def log_prior(self):
        # print self.p, self.q
        return 0.


class MyData(PredictorSet):
    d1 = SampledGaussian(locns1, vals1, sigs2)
    d2 = SampledGaussian(locns2, vals2, sigs2)


class MyBayes(BayesianInference, MyModel, MyData):
    pass


class MyChi(ChisqrInference, MyModel, MyData):
    pass


print('Instantiating MyInf...')
m = MyBayes()

print('params: ', m.params)
print('predictors: ', m.predictors)

print('Varying p & q to find (2,1):')
m.p.vary(2.5, .1, (-5, 5.))
m.q.vary(1.5, .1, (-5., 5.))
m.map_bounds()
print(m.fit())
print()

print('Info & covar matrix:')
info = m.obs_info()
covar = linalg.inv(info)
print(info)
print(covar)
print()

print('g2 = Step p and q by 5 & 3:')
m.p.range(0, 4)
m.p.step(0, 4, 50)
m.q.step(-1, 3, 30)
g2 = m.fit()
print(g2)

print('pg = Step p with q fixed:')
m.q = 1.
pg = m.fit()
print(pg)

print('qg = Step q with p fixed:')
m.p = 2.
m.q.step(-1.5, 3.5, 30)
qg = m.fit()
print(qg)

print('qpro = q profile:')
m.p.vary(2., .1, (0, 4.))
qpro = m.fit()
print(qpro)

print('ppro = p profile:')
m.q.vary(3., .1, (-10., 10.))
m.p.step(1, 3, 25)
ppro = m.fit()
print(ppro)
