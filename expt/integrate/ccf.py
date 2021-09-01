"""
Fast Clenshaw-Curtis quadrature rules using Waldvogel's DFT-based algorithm

See:

Waldvogel, J. (2006) "Fast Construction of the Fejer and
Clenshaw-Curtis Quadrature Rules," BIT Num. Math., 43, 195-202

Created 31 Mar 2009 by Tom Loredo
"""

from numpy import *


def ccrule(order, dft=None):
    """
    Calculate nodes and weights for a Clenshaw-Curtis rule over [-1,1].
    
    Note that the order is the degree of the polynomial that is integrated
    exactly; the number of nodes is order+1.  This is a closed rule, with nodes
    at the boundaries.
    
    If dft=True, DFTs of the weights for the CC rule and Fejer's 2nd rule
    are additionally returned.
    
    The algorithm uses an FFT of size order; for very large order, setting
    order to a power of 2 will significantly improve performance.
    """
    npts = order+1
    # The nodes are formally defined as x=cos(theta) for equally spaced
    # theta; reverse them so they increase in x (the rule is symmetric).
    nodes = cos(pi*arange(npts, dtype=float)/order)[::-1]
    # Be careful to use float here or later when dividing by odds!
    odds = arange(1, order, 2, dtype=float)
    l = len(odds)
    m = order-l
    v0 = zeros(npts)
    v0[:l] = 2./odds/(odds-2)
    v0[l] = 1./odds[-1]
    # v2 is the DFT of the weights for Fejer's 2nd rule.
    v2 = -v0[:-1] - v0[-1:0:-1]  # (up to penult) - (last to 2nd)
    g0 = -ones(order)
    g0[l] += order
    g0[m] += order
    # v2+g is the DFT of the weights for Clenshaw-Curtis.
    g = g0/(order**2 - 1 + mod(order,2))
    wts = empty(npts)
    wts[:-1] = fft.ifft(v2+g).real
    wts[-1] = wts[0]
    # Make nodes contiguous before returning.
    if dft:
        return ascontiguousarray(nodes), wts, v2+g, v2
    else:
        return ascontiguousarray(nodes), wts


def ccrule2(order):
    """
    Calculate nodes and weights for a Clenshaw-Curtis rule over [-1,1].
    
    Note that the order is the degree of the polynomial that is integrated
    exactly; the number of nodes is order+1.  This is a closed rule, with nodes
    at the boundaries.
    """
    npts = order+1
    # The nodes are formally defined as x=cos(theta) for equally spaced
    # theta; reverse them so they increase in x (the rule is symmetric).
    nodes = cos(pi*arange(npts, dtype=float)/order)[::-1]


# N=N1-1; bma=b-a;
# c=zeros(N1,2);
# c(1:2:N1,1)=(2./[1 1-(2:2:N).^2 ])'; c(2,2)=1;
# f=real(ifft([c(1:N1,:);c(N:-1:2,:)]));
# w=bma*([f(1,1); 2*f(2:N,1); f(N1,1)])/2;
# x=0.5*((b+a)+N*bma*f(1:N1,2));


def f1rule(order, dft=None):
    """
    Calculate nodes and weights for Fejer's 1st rule over [-1,1].
    
    Note that the order is the degree of the polynomial that is integrated
    exactly; the number of nodes is order+1.  This is an open rule, with
    no nodes at the boundaries.
    
    If dft=True, DFTs of the weights are additionally returned.

    The algorithm uses an FFT of size order+1; for very large order, choosing
    order+1 equal to a power of 2 will significantly improve performance.
    """
    npts = order+1
    # The nodes are formally defined as x=cos(theta) for equally spaced
    # theta; reverse them so they increase in x (the rule is symmetric).
    kvals = arange(npts, 0, -1, dtype=float)-0.5  # [npts-1/2, ..., 1/2]
    nodes = cos(pi*kvals/npts)
    l = npts / 2  # int divide
    m = npts-l
    K = arange(m, dtype=float)
    v0 = zeros(npts+1, dtype=complex)
    v0[:m] = 2*exp(1j*pi*K/npts)/(1-4*K**2)
    # v1 is the DFT of the weights for Fejer's 1st rule.
    v1 = v0[:-1] + v0[-1:0:-1].conj()  # (up to penult) - (last to 2nd)
    wts = fft.ifft(v1).real
    if dft:
        return nodes, wts, v1
    else:
        return nodes, wts


class ClenshawCurtis:
    """
    Clenshaw-Curtis quadrature with fast calculation of the quadrature rule
    """
    
    def __init__(self, order, l=None, u=None, err=False):
        """
        Define a Clenshaw-Curtis quadrature rule of specified order.
        
        The number of nodes in the rule is order+1; it integrates polynomials 
        of degree <= order exactly.  It is a closed rule; the nodes include
        the boundary points.
        
        The default range is [-1, 1]; if provided, [u, l] is the range.
        
        If err=True and order is even, an order/2 rule is defined to use
        for error estimation.
        """
        self.order = order
        self.npts = order+1
        # Store the raw weights for [-1,1].
        self._nodes, self._wts = ccrule(order)
        # Store weights for the specified interval.
        self.set_range(l, u)
        if err:
            raise NotImplementedError

    def set_range(self, l, u):
        """
        Set the range for the quadrature.
        """
        if l is None and u is None:
            self.l = -1.
            self.u = 1.
            self.nodes = self._nodes
            self.wts = self._wts
        elif l is None or u is None:
            raise ValueError('Specify both integration limits!')
        else:
            self.l = l
            self.u = u
            hw = 0.5*(u-l)
            # Be careful with the range; if one endpoint is much smaller than
            # the other, it can get lost in roundoff.
            self.nodes = hw*self._nodes + 0.5*(l+u)
            self.nodes[0] = l
            self.nodes[-1] = u
            self.wts = hw*self._wts

    def quad(self, f, l=None, u=None):
        """
        Return the Clenshaw-Curtis quadrature of f(x).
        """
        if l is None and u is None:
            if callable(f):
                vals = f(self.nodes)
            else:
                vals = f
            return sum(self.wts*vals)
        else:  # f must be callable when l, u given
            hw = 0.5*(u-l)
            nodes = hw*self._nodes + 0.5*(l+u)
            wts = hw*self._wts
            vals = f(nodes)
            return sum(wts*vals)

    def _explicit(self):
        l = int(floor(self.order/2))
        v = 1000*ones(self.order)  # just to check all elements get touched
        for k in range(l):
            v[k] = 2/(1-4.*k**2)
        v[l] = (self.order-3.)/(2*l-1) -1
        for k in range(1, int((self.order-1)/2.)+1):
            v[self.order-k] = v[k]
        return v


class Fejer1:
    """
    Fejer's 1st rule quadrature with fast calculation of the  rule via FFT
    """
    
    def __init__(self, order, l=None, u=None):
        """
        Define a Fejer's 1st quadrature rule of specified order.
        
        The number of nodes in the rule is order+1; it integrates polynomials 
        of degree <= order exactly.  It is an open rule; the boundaries are
        not among the nodes.
        
        The default range is [-1, 1]; if provided, [u, l] is the range.
        """
        self.order = order
        self.npts = order+1
        # Store the raw weights for [-1,1].
        self._nodes, self._wts = f1rule(order)
        # Store weights for the specified interval.
        self.set_range(l, u)

    def set_range(self, l, u):
        """
        Set the range for the quadrature.
        """
        if l is None and u is None:
            self.l = -1.
            self.u = 1.
            self.nodes = self._nodes
            self.wts = self._wts
        elif l is None or u is None:
            raise ValueError('Specify both integration limits!')
        else:
            self.l = l
            self.u = u
            hw = 0.5*(u-l)
            self.nodes = hw*self._nodes + 0.5*(l+u)
            self.wts = hw*self._wts

    def quad(self, f, l=None, u=None):
        """
        Return the Clenshaw-Curtis quadrature of f(x).
        """
        if l is None and u is None:
            if callable(f):
                vals = f(self.nodes)
            else:
                vals = f
            return sum(self.wts*vals)
        else:  # f must be callable when l, u given
            hw = 0.5*(u-l)
            nodes = hw*self._nodes + 0.5*(l+u)
            wts = hw*self._wts
            vals = f(nodes)
            return sum(wts*vals)


if __name__ == '__main__':

    def f(x):
        """
        Integrates to 2 over [-1,1]; to interval size over any
        other interval.
        """
        return 1

    def g(x):
        """
        Integrates to 2 over [-1,1]; to u**3 - l**3 over [l,u].
        """
        return 3*x**2

    def h(x):
        """
        Integrates to 4 over [-1,1] (odd term vanishes); to  101000000887.6
        over [2, 10].
        """
        return 11*x**10 + x**9 + 3*x**2
    
    cc8 = ClenshawCurtis(8)
    print 'Testing ClenshawCurtis:'
    print 'Should be 2, 2, 4:\n', cc8.quad(f), cc8.quad(g), cc8.quad(h)
    
    # An order-8 example tabulated by John Burkardt
    # http://people.sc.fsu.edu/~burkardt/f_src/clenshaw_curtis/clenshaw_curtis.html
    bcc8_nodes = array([   -1.0000000000,
   -0.9238795325,
   -0.7071067812,
   -0.3826834324,
    0.0000000000,
    0.3826834324,
    0.7071067812,
    0.9238795325,
    1.0000000000 ])
    bcc8_wts = array([    0.0158730159,
    0.1462186492,
    0.2793650794,
    0.3617178587,
    0.3936507937,
    0.3617178587,
    0.2793650794,
    0.1462186492,
    0.0158730159 ])
    print 'bcc8 quad:', sum(bcc8_wts*f(bcc8_nodes)), sum(bcc8_wts*g(bcc8_nodes)),\
        sum(bcc8_wts*h(bcc8_nodes))
    cc = ClenshawCurtis(10)
    print 'Exact for 10th degree:', cc.quad(h)
    cc7 = ClenshawCurtis(7)
    cc7.set_range(2, 10)
    print 'Should be 8, 992, 100999998841.6:'
    print '  Order 7: ', cc7.quad(f), cc7.quad(g), cc7.quad(h)
    cc10 = ClenshawCurtis(10, 2, 10)
    print '  Order 10:', cc10.quad(f), cc10.quad(g), cc10.quad(h, 2, 10)

    print 'Comparing 5th order explicit and fast DFT vectors:'
    n, w, dft, v2 = ccrule(10, True)
    v = cc._explicit()
    print v2-v
    print
    
    print 'Testing Fejer1:'
    f1 = Fejer1(7)
    print 'Should be 2, 2, 4:\n', f1.quad(f), f1.quad(g), f1.quad(h)
    f1 = Fejer1(10)
    print 'Exact for 10th degree:', f1.quad(h)
