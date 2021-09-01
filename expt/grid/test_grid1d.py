from numpy import *
from inference.grid.grid1d import PeakParams, LabeledMaxList, BasicInterval, \
    RefinedIntervalList1D, BasicGrid1D, Grid1D
from inference.utils.ioutils import safe_open

if 0:
    print('Test of PeakParams:')
    x = [0, 1., 3.]
    y = [-8., -2., -2.]
    print(PeakParams(x, y))

if 0:
    print("Testing LML:")
    lml = LabeledMaxList(3)
    lml.update(3.,"three")
    print(lml)
    lml.update(3.1,"three.one")
    print(lml)
    lml.update(1.,"one")
    print(lml)
    lml.update(10.,"ten")
    print(lml)
    lml.update(20.,"twenty")
    print(lml)
    lml.update(4.,"four")
    print(lml)
    from copy import copy
    lml2 = copy(lml.list)
    lml.limitRange(3.)
    print(lml)
    print(lml2)

if 0:
    print("Testing BasicInterval:")
    bi = BasicInterval(5., 10., 6)
    print(bi.locate(5.))
    print(bi)
    bi = BasicInterval(1., 100., 5, stype="log")
    print(bi.locate(10.))
    print(bi)

if 0:
    print("Testing RefinedIntervalList1D:")
    cil = RefinedIntervalList1D(5., 10., 6)
    # cil.merge()
    print(cil)
    cil.update(7., 8., 5)
    print(cil)
    print("# of pts: ", cil.npts)
    print("cum: ", cil.cumpts)
    cil.update(6., 9., 16)
    cil.update(7., 8., 11)
    print(cil)
    print("# of pts: ", cil.npts)
    print("cum: ", cil.cumpts)
#    cil.merge()
#    print cil
    print(list(map(str,cil.nonbase())))
    print(cil.locate(7.5))
    print(cil.nlocate(6))

if 1:
    print("Testing BasicGrid1D:")
    from numpy import log
    lml = LabeledMaxList(4)
    sg1 = BasicGrid1D(10, 1000, 200, lml=lml)
#    def f(x): return -((x-300.)*(x-300.))/2.e4
#    def f(x): return log(x)

    def f(x):
        lx = zeros((len(x),2), float)
        lx[:,0] = log(x)
        lx[:,1] = x
        return lx
    sg1.vecfunc(f,0)
#    print sg1
    print(sg1.trapzd(0), sg1.logtrapzd(0))
    print(sg1.lml, sg1.point(2))

if 1:
    print("Testing Grid1D:")
    rg = Grid1D(0., 100., 101, sigpts=101, nsig=5)

    def bexp(x):
        if x > -300.:
            return exp(x)
        else:
            return 1.e-300

    def f(x): return log(bexp(-(x-50.)**2/2.) + .2*bexp(-(x-72.)**2/2.))
#    def f(x): return -(x-50.234)**2/2., bexp(-(x-50.234)**2/2.)
#    def f(x): return log( 2*x )
    rg.firstpass(f,0)
    print(rg.base.lml)
# Use the next with log(2*x) & print the contributions to logtrapzd
#    rg.lml = rg.base.lml = [(2., (1, 0., 10., 2.)), (1., (2, 0., 70., 2.))]
    print("Base trapzd =    ", rg.logtrapzd())
    rg.refine()
    print("Refined trapzd = ", rg.logtrapzd())
    print(rg.ilist)
    print(rg.point(98))
    ofile = safe_open('rg.dat')
    print(rg.totpts)
    for i in range(rg.totpts):
        point = rg.point(i)
        line = "%5i  %#10.4g  %#10.4g\n" % (i, point[0], point[2])
        ofile.write(line)
    ofile.close()

    t = 1.2*sqrt(2*pi)
    print('true value = ', repr(log(t)))
    from numpy import arange
    fsum = 0.
    d = 1.
    for x in arange(0., 100., d):
        fsum += 0.5*d*(bexp(f(x)) + bexp(f(x+d)))
    print(log(fsum), fsum, (fsum-t)/t)

    fsum = 0.
    d = 1.
    # print arange(0.,45.,d)
    for x in arange(0., 45., d):
        fsum += 0.5*d*(bexp(f(x)) + bexp(f(x+d)))
    d = .02
    # print arange(45., 55., d)
    for x in arange(45., 55., d):
        fsum += 0.5*d*(bexp(f(x)) + bexp(f(x+d)))
    d = 1.
    for x in arange(55., 67., d):
        fsum += 0.5*d*(bexp(f(x)) + bexp(f(x+d)))
    d = .02
    for x in arange(67., 77., d):
        fsum += 0.5*d*(bexp(f(x)) + bexp(f(x+d)))
    d = 1.
    for x in arange(77., 100., d):
        fsum += 0.5*d*(bexp(f(x)) + bexp(f(x+d)))
    print(log(fsum), fsum, (fsum-t)/t)

    print('------- log test -------')
    rg = Grid1D(1., 100., 100, stype="log")
    def f(x): return 1/x
    rg.firstpass(f)
    print(rg.base.lml)
    print("Base trapzd =    ", rg.trapzd())
    rg.refine()
    print("Refined trapzd = ", rg.trapzd())
    print('True value = ', log(100.))
    print(rg.ilist)
    nlast = len(rg)-1
    print('Point', nlast, ':', rg.point(nlast))
    ofile = safe_open('logrg.dat')
    print(rg.totpts)
    for i in range(rg.totpts):
        point = rg.point(i)
        line = "%5i  %#10.4g  %#10.4g\n" % (i, point[0], point[2])
        ofile.write(line)
    ofile.close()
    print('****** Done! *******')
