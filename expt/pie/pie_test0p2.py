from pie import *
from pie.parameterizedmodel import StopStepping
from pie.chisqr import SampledChisqrPredictor
from Numeric import array, Float

################################################################################
#  Testing, 1 2 3...
#

class test(ParameterizedModel):
    #def __init__(self):
    #    ParameterizedModel.__init__(self)

    x = RealParam(1.2)
    y = RealParam(3.4, '2nd param')

if 0:
    fph = RealParamHandler()
    f = RealParamValue(1.23, fph)
    print f, 2*f, f-1
    print f.addOne()

a = test()
print a.x, a.y, a.x.addOne(), a.y.addOne()
b = test()
b.x = 5.6
b.y = 7.8
print 'a:', a.x, a.y, a.x.addOne(), a.y.addOne()
print 'b:', b.x, b.y, b.x.addOne(), b.y.addOne()
print 'a dict:', a.__dict__
print 'b dict:', b.__dict__
print 'test dict:', test.__dict__
a.x.step(1,3,5)
a.y.logStep(1,1000,4)
def doGrid():
    while True:
        try:
            a.nextStep()
            print a.stepNums()
        except StopStepping:
            print 'Completed the grid!'
            break
doGrid()
# a.x.vary(1., .1, (1,5))

from Numeric import array
locns = array([1,2,3],Float)
vals = array([10,20,30],Float)
sigs = array([4,5,6],Float)
scp = SampledChisqrPredictor()
scp.setData(locns, vals, sigs)
#locns = array([(1,2),(1,3),(2,3)],Float)
#scp.setData(locns, vals, sigs)
all = array([(1,10,4),(2,20,5),(3,30,6)],Float)
#all = [((1,2),10,4),((1,3),20,5),((2,3),30,6)]
#all = [(1,10,4),(2,20,5),(3,30,6)]
scp.setData(all)

class MyModel(ParameterizedModel):
    #def __init__(self):
    #    ParameterizedModel.__init__(self)

    p = RealParam(1.2)
    q = RealParam(3.4, '2nd param')

    def signal(self, x):
        return self.p*self.q*x

class Mix(MyModel, SampledChisqrPredictor):
    pass

m = Mix()
m.setData(all)
print 'chisqr =', m.chisqr()
m.setData(locns, vals, sigs)
print 'chisqr =', m.chisqr()

class MyInf(MyModel, SampledChisqrInference):
    pass

inf = MyInf()
inf.setData(locns, vals, sigs)
inf.p.step(1,3,5)
inf.q.logStep(1,1000,4)
results = inf.doGrid()
inf.q.vary(5.,.1,(1.,1000))
