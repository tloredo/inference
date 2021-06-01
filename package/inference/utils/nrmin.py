from numpy import array, identity, log, exp, clip


class MinBracket:
    """Bracket a minimum of a function of one variable."""

    def __init__(self, f, x1, x2):
        """Pass the function and two starting points that
        need not bracket the minimum."""
        GOLD = 1.618034
        TINY = 1.e-20
        MAXMAG = 100.
        self.x1 = x1
        self.f1 = f(x1)
        self.xm = x2
        self.fm = f(x2)
        if (self.fm > self.f1):
            self.x1, self.xm = self.xm, self.x1
            self.f1, self.fm = self.fm, self.f1
        self.x2 = self.xm + GOLD*(self.xm - self.x1)
        self.f2 = f(self.x2)
        while 1:
            if (self.fm < self.f2):
                break
            d1 = (self.xm - self.x1)
            d2 = (self.xm - self.x2)
            r = d1*(self.fm - self.f2)
            q = d2*(self.fm - self.f1)
            u = 2.*(q - r)
            if (abs(u) < TINY):
                if (u >= 0):
                    u = TINY
                else:
                    u = -TINY
            u = self.xm - (d2*q - d1*r)/u
            ulim = self.xm - MAXMAG*d2
            if ((self.xm-u)*(u-self.x2) > 0.):
                fu = f(u)
                if (fu < self.f2):
                    self.x1, self.xm = self.xm, u
                    self.f1, self.fm = self.fm, fu
                    break
                elif (fu > self.fm):
                    self.x2, self.f2 = u, fu
                    break
                u = self.x2 - GOLD*d2
                fu = f(u)
            elif ((self.x2-u)*(u-ulim) > 0.):
                fu = f(u)
                if (fu < self.f2):
                    self.xm, self.x2 = self.x2, u
                    self.fm, self.f2 = self.f2, fu
                    u = self.x2 - GOLD*d2
                    fu = f(u)
            elif ((u-ulim)*(ulim-self.x2) >= 0.):
                u = ulim
                fu = f(u)
            else:
                u = self.x2 - GOLD*d2
                fu = f(u)
            self.x1, self.xm, self.x2 = self.xm, self.x2, u
            self.f1, self.fm, self.f2 = self.fm, self.f2, fu
#           print '  Bracketing ', iter, u, fu
        if (self.x2 < self.x1):
            self.x1, self.x2 = self.x2, self.x1
            self.f1, self.f2 = self.f2, self.f1


# IterErr = "Too many iterations!"  # exceptions once could be strings...

class IterErr(Exception):
    pass


def brent(f, bracket, tol=1.e-6, maxit=100):
    """Find a minimum of f(x) using Brent's method."""
    CGOLD = 0.3819660
    EPS0 = 1.e-10
    a = bracket.x1
    b = bracket.x2
    v = bracket.xm
    w, x = v, v
    fv = f(v)
    fw, fx = fv, fv
    e = 0.
    iter = 0
    while 1:
        iter = iter + 1
        if (iter > maxit):
            raise IterErr('Too many iterations!')
        xm = 0.5*(a+b)
        tol1 = tol*abs(x) + EPS0
        tol2 = 2.*tol1
        if (abs(x-xm) <= (tol2-0.5*(b-a))):
            return (x, fx)
        if (abs(e) > tol1):
            xw = x - w
            xv = x - v
            r = xw*(fx-fv)
            q = xv*(fx-fw)
            p = xv*q - xw*r
            q = 2.*(q-r)
            if (q > 0.):
                p = -p
            q = abs(q)
            etemp = e
            e = d
            if (abs(p) >= abs(0.5*q*etemp) or p <= q*(a-x) or p >= q*(b-x)):
                if (x >= xm):
                    e = a-x
                else:
                    e = b-x
                d = CGOLD*e
            else:
                d = p/q
                u = x+d
                if (u-a < tol2 or b-u < tol2):
                    if (xm-x < 0.):
                        d = -tol1
                    else:
                        d = tol1
        else:
            if (x >= xm):
                e = a-x
            else:
                e = b-x
            d = CGOLD*e
        if (abs(d) >= tol1):
            u = x+d
        elif (d < 0):
            u = x-tol1
        else:
            u = x+tol1
        fu = f(u)
#       print '  Brent', iter, fu
        if (fu <= fx):
            if (u >= x):
                a = x
            else:
                b = x
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
        else:
            if (u < x):
                a = u
            else:
                b = u
            if (fu <= fw or w == x):
                v, fv = w, fw
                w, fw = u, fu
            elif (fu <= fv or v == x or v == w):
                v, fv = u, fu


def linmin(func, pt, drxn, tol=1.e-4):
    """Minimize f along drxn, starting at pt."""
    def f1dim(x, f=func, p=pt, d=drxn):
        pp = p + x*d
        return f(pp)
    brak = MinBracket(f1dim, 0., 1.)
    l, f = brent(f1dim, brak, tol=tol)
    drxn = l*drxn
    pt = pt + drxn
    return (pt, drxn, f)


def powell(f, p, drxns=None, ftol=1.e-6, maxit=200, nlog=None, shuffle=False):
    """
    Minimize a multivariate function using Powell's direction set method.
        f = function of a 1-d array argument
        p = starting pt
        drxns = list of direction vectors (i.e., *rows* are directions)

        shuffle = reorder the initial directions

    Returns:
        Minimizing point (array)
        Function value at the minimum
        Number of iterations
    """
    n = len(p)
    if drxns is None:
        drxns = identity(n)
        drxns = array(drxns, float)
    if shuffle:
        vec = drxns[0].copy()  # must copy!
        for i in range(n-1):
            drxns[i] = drxns[i+1]
        drxns[-1] = vec
    fret = f(p)
    p0 = p
    iter = 0
    while 1:
        iter = iter + 1
        # print('Powell', iter)
        fp = fret
        ibig, delta = 0, 0.
        for i in range(n):
            vec = drxns[i]
            fptt = fret
            p, vec, fret = linmin(f, p, vec)
            if abs(fptt-fret) > delta:
                delta = abs(fptt-fret)
                ibig = i
        if 2.*abs(fp-fret) <= ftol*(abs(fp)+abs(fret)):
            return (p, fret, iter)
        if iter == maxit:
            raise IterErr('Too many iterations!')
        if nlog is not None and iter % nlog == 0:
            print('Powell iteration', iter, fret)
        ptt = 2.*p - p0
        vec = p - p0
        p0 = p
        fptt = f(ptt)
#       print iter, tuple(ptt), fptt
        if fptt >= fp:
            continue
        t = 2.*(fp - 2.*fret + fptt) * (fp-fret-delta)**2 - delta*(fp-fptt)**2
        if t >= 0.:
            continue
        p, vec, fret = linmin(f, p, vec)
        drxns[ibig] = drxns[-1]
        drxns[-1] = vec


if __name__ == '__main__':

    def test_brent():
        '''
        Test brent on a smooth function with minumum at 2.
        '''
        def func(x): return -exp(-(x-2.)**2)

        fbrak = MinBracket(func,0.,.1)
        print("bracket: ", fbrak.x1, fbrak.xm, fbrak.x2)
        print("         ", fbrak.f1, fbrak.fm, fbrak.f2)
        print(brent(func, fbrak))

    def test_linmin():
        '''
        Test linmin on a 2-D function with min at x0=2, x1=3.
        '''
        def func(x):
            return -exp(-(x[0]-2.)**2 - (x[1]-3.)**2)
        pt = array([1., 2.], float)
        drxn = array([.1, .1])
        print(linmin(func, pt, drxn))

    # test_linmin()

    def test_powell():
        def func(x):
            # return exp((.1*x[0]-2.)**2 + (.1*x[1]-3.)**2 + abs(x[0]*x[1]))
            return exp((.1*x[0]-2.)**2 + (.1*x[1]-3.)**2)
        pt = array([1., 2.], float)
        print(powell(func, pt))

    # test_powell()
