"""
Fast Clenshaw-Curtis quadrature rules using Waldvogel's DFT-based algorithm

See:

Waldvogel, J. (2006) "Fast Construction of the Fejer and
Clenshaw-Curtis Quadrature Rules," BIT Num. Math., 43, 195-202

Created 5 Apr 2009 by Tom Loredo
"""

from numpy import *


class Quad:
    """
    A simple quadrature rule
    """
    
    def __init__(self, *args):
        """
        Define a quadrature rule from its nodes and weights.
        
        For closed rules (where the nodes include the boundaries), the
        signature may be:
        
            QuadRule(nodes, wts)
        
        where `nodes` and `wts` are arrays containing the (ordered)
        nodes and weights for the rule.
        
        For open rules (and optionally for closed rules), the signature is
        
            QuadRule(l, u, nodes, wts)

        where `a` and `b` are the integration limits, and `nodes` and `wts` 
        are arrays containing the (ordered) nodes and weights for the rule.
        """
        if len(args) == 2:
            self.type = 'closed'
            self.nodes = args[0]
            self.wts = args[1]
            self.l = self.nodes[0]
            self.u = self.nodes[-1]
        elif len(args) == 4:
            self.l = args[0]
            self.u = args[1]
            self.nodes = args[2]
            self.wts = args[3]
            if self.l==self.nodes[0] and self.u==self.nodes[-1]:
                self.type = 'closed'
            else:
                self.type = 'open'
        else:
            raise ValueError('Specify (nodes, wts) or (a, b, nodes, wts)!')
        self.npts = len(self.nodes)
        if len(self.wts) != self.npts:
            raise ValueError('Mismatch in length of nodes and wts!')

    def quad(self, f):
        """
        Evaluate the quadrature given a function or array of function values.
        
        If f is callable, find the array of values f(self.nodes) and
        evaluate the quadrature.  Otherwise, f must be a vector of function
        values at the nodes, used to evaluate the quadrature.
        """
        if callable(f):
            fvals = f(self.nodes)
        elif len(f) == self.npts:
            fvals = f
        else:
            raise ValueError('Argument must be callable or array of values!')
        return sum(self.wts*fvals)


class PLMapQuad(Quad):
    """
    A quadrature rule with nodes and weights adjusted to absorb a power-law factor
    """
    
    def __init__(self, indx, *args):
        """
        Define a quadrature rule that absorbs a power-law factor.
        
        Define a quadrature rule that absorbs a power-law factor in the
        integrand using nodes and weights from a basic rule.  The basic
        rule will be applied using a transformed ordinate whose Jacobian
        is the inverse of the specified power law, "flattening" the
        integrand and hopefully improving quadrature accuracy.
        
        `indx` specifies the index for the power law factor, defined
        for ordinate x according to x**(-indx), i.e., indx is the
        *negative* log slope of the power law factor.
        
        The remaining arguments define the basic quadrature rule that will
        be transformed.  For closed rules (where the nodes include the 
        boundaries), the signature may be:
        
            PLMapQuad(indx, nodes, wts)
        
        where `nodes` and `wts` are arrays containing the (ordered)
        nodes and weights for a rule defined over x.
        
        For open rules (and optionally for closed rules), the signature is
        
            PLMapQuad(indx, l, u, nodes, wts)

        where `l` and `u` are the integration limits, and `nodes` and `wts` 
        are arrays containing the (ordered) nodes and weights for a rule
        defined over x.
        
        The transformed rule will span [l,u] in x, but will have different
        nodes and weights.
        """
        self.indx = indx
        Quad.__init__(self, *args)
        # Save the base rule before transforming.
        self._nodes = self.nodes
        self._wts = self.wts
        delta_x = self.u - self.l
        # For indx near 1, use the indx=1 logarithmic case.
        if abs(indx-1.) < 0.2:
            # Transform the rule to span the new y coordinate.
            self.y_l = 0.
            self.y_u = log(self.u/self.l)
            scale = self.y_u / delta_x
            self.y_nodes = (self._nodes - self.l) * scale
            self.y_wts = self._wts * scale
            # Transform the y rule to an x rule.
            self.nodes = self.l * exp(self.y_nodes)
            self.wts = self.y_wts * self.nodes
        else:
            i1 = 1. - indx
            lpow = self.l**i1
            # Transform the rule to span the new y coordinate.
            self.y_l = 0.
            self.y_u = lpow * ((self.u/self.l)**i1 - 1.) / i1
            scale = self.y_u / delta_x
            self.y_nodes = (self._nodes - self.l) * scale
            self.y_wts = self._wts * scale
            # Transform the y rule to an x rule.
            logn = log(self.l) + log(1. + i1*self.y_nodes/lpow)/(1.-indx)
            self.nodes = exp(logn)
            log_fac = indx * logn
            self.wts = self.y_wts * exp(log_fac)

    def quad(self, f):
        """
        Evaluate the quadrature given a function or array of function values.
        
        If f is callable, find the array of values f(self.nodes) and
        evaluate the quadrature.  Otherwise, f must be a vector of function
        values at the nodes, used to evaluate the quadrature.
        """
        if callable(f):
            fvals = f(self.nodes)
        elif len(f) == self.npts:
            fvals = f
        else:
            raise ValueError('Argument must be callable or array of values!')
        return sum(self.wts*fvals)


class ExpMapQuad(Quad):
    """
    A quadrature rule adjusted to absorb an exponential factor in the integrand
    """
    
    def __init__(self, alpha, *args):
        """
        Define a quadrature rule that absorbs an exponential factor.
        
        Define a quadrature rule that absorbs an exponential factor in the
        integrand using nodes and weights from a basic rule.  The basic
        rule will be applied using a transformed ordinate whose Jacobian
        is the inverse of the exponential, "flattening" the
        integrand and hopefully improving quadrature accuracy.
        
        `alpha` specifies the exponential slope, defined for ordinate x
        according to exp(alpha*x).
        
        The remaining arguments define the basic quadrature rule that will
        be transformed.  For closed rules (where the nodes include the 
        boundaries), the signature may be:
        
            ExpMapQuad(alpha, nodes, wts)
        
        where `nodes` and `wts` are arrays containing the (ordered)
        nodes and weights for a rule defined over x.
        
        For open rules (and optionally for closed rules), the signature is
        
            ExpMapQuad(alpha, l, u, nodes, wts)

        where `l` and `u` are the integration limits, and `nodes` and `wts` 
        are arrays containing the (ordered) nodes and weights for a rule
        defined over x.
        
        The transformed rule will span [l,u] in x, but will have different
        nodes and weights.
        """
        self.alpha = alpha
        Quad.__init__(self, *args)
        # Save the base rule before transforming.
        self._nodes = self.nodes
        self._wts = self.wts
        delta_x = self.u - self.l
        # Transform the rule to span the new y coordinate.
        self.y_l = exp(alpha*self.l)/alpha
        self.y_u = exp(alpha*self.u)/alpha
        scale = (self.y_u - self.y_l) / delta_x
        self.scale = scale
        # Since y is exponential in x, be careful about roundoff when
        # calculating y nodes; offset from the boundary closest to 0 since
        # it may be exponentially small.
        if abs(self.y_l) < abs(self.y_u):
            self.y_nodes = self.y_l + (self._nodes - self.l) * scale
        else:
            self.y_nodes = self.y_u - (self.u - self._nodes) * scale
        self.y_wts = self._wts * scale
        # Transform the y rule to an x rule.
        self.nodes = log(alpha*self.y_nodes)/alpha
        self.wts = self.y_wts / self.y_nodes / alpha

    def quad(self, f):
        """
        Evaluate the quadrature given a function or array of function values.
        
        If f is callable, find the array of values f(self.nodes) and
        evaluate the quadrature.  Otherwise, f must be a vector of function
        values at the nodes, used to evaluate the quadrature.
        """
        if callable(f):
            fvals = f(self.nodes)
        elif len(f) == self.npts:
            fvals = f
        else:
            raise ValueError('Argument must be callable or array of values!')
        return sum(self.wts*fvals)


class CompositeQuad:
    """
    Composite quadrature rule built over contiguous intervals
    """

    @staticmethod
    def isquad(obj):
        """
        Return True if obj has a QuadRule interface.
        """
        if isinstance(obj, Quad):
            return True
        try:
            assert hasattr(obj, 'l')
            assert hasattr(obj, 'u')
            assert hasattr(obj, 'npts')
            assert hasattr(obj, 'nodes')
            assert hasattr(obj, 'wts')
            return True
        except AssertionError:
            return False

    def __init__(self, *args):
        """
        Define a composite quadrature rule from a collection of rules.
        
        Each argument should describe a basic quadrature rule.  There are
        three possible formats:
        
          * The argument may be a QuadRule instance.
        
          * For a closed rule, the argument may be a 2-tuple of the
            form (nodes, wts), where `nodes` and `wts` are arrays of
            values for the nodes and weights of the rule.
          
          * For an open or closed rule, the argument may be a 4-tuple of the
            form (a, b, nodes, wts), where `a` and `b` are the limits of
            integration, and `nodes` and `wts` are as above.
        
        The rules must be contiguous.
        """
        # Collect all the rules as QuadRule instances.
        self.rules = []
        last = None
        for i, arg in enumerate(args):
            if CompositeQuad.isquad(arg):
                if last:
                    if last.u != arg.l:
                        raise ValueError(\
                            'Rule %i not contiguous with previous!' % (i+1))
                self.rules.append(arg)
                last = arg
            else:
                rule = Quad(*arg)
                self.rules.append(rule)
                last = rule
        # Make an array of all the distinct nodes.
        prev = self.rules[0]
        distinct = [ prev.nodes ]
        self.starts = [0]  # keeps track of starting indices for rules
        self.npts = prev.npts
        for rule in self.rules[1:]:
            if prev.nodes[-1] != rule.nodes[0]:
                distinct.append(rule.nodes)
                self.starts.append(self.npts)
                self.npts += rule.npts
            else:
                distinct.append(rule.nodes[1:])
                self.starts.append(self.npts-1)
                self.npts += rule.npts - 1
            prev = rule
        self.nodes = concatenate(distinct)
        self.factor = None
        self.factors = None
    
    def set_factor(self, factor):
        """
        Set an integrand factor that will be constant for repeated evaluations.
        
        Consider a set of quadratures of integrands that may be written as
        f(x)*g(x), where the factor g(x) will be the same for each of the
        quadratures.  To save work, g(x) can be computed just once on the
        nodes using this method.

        Note that for such cases (inner) product quadrature rules exist,
        potentially evaluating f() and g() on *different* nodes, with
        higher order exactness than regular quadrature rules.
        """
        self.factor = factor
        self.factors = factor(self.nodes)

    def reset_factor(self, factor):
        """
        Reset the composite quadrature to not use a pre-calculated integrand factor.
        """
        self.factor = None
        self.factors = None

    def quad(self, f):
        """
        Evaluate the quadrature of the callable f.
        """
        fvals = f(self.nodes)
        self.fvals = fvals
        if self.factor:
            self.ivals = self.factor * self.fvals
        else:
            self.ivals = self.fvals
        result = 0.
        for rule, start in zip(self.rules, self.starts):
            result += rule.quad(self.ivals[start:start+rule.npts])
        return result


if __name__ == '__main__':
    from ccf import ClenshawCurtis, Fejer1
    cc = ClenshawCurtis(6, 3, 5)
    f1 = Fejer1(5, 5, 8)
    # 4-part composite rule over [0, 10] with CC and Fejer1 pieces:
    cq = CompositeQuad(ClenshawCurtis(4, 0, 3),
            (cc.nodes, cc.wts),
            (f1.l, f1.u, f1.nodes, f1.wts),
            Fejer1(7, 8., 10) )
    
    def f(x):
        """
        Integrates to 10 over [0,10].
        """
        try:
            return ones(len(x), float)
        except TypeError:
            return 1.

    def g(x):
        """
        Integrates to 1000 over [0,10].
        """
        return 3*x**2

    def h(x):
        """
        Integrates to 1.e11 + 1.e9 + 1000 over [0,10].
        """
        return 11*x**10 + x**9 + 3*x**2

    print 'Single order-10 CC & F1 quads for f=1, g=3*x**2, h(x)=10th degree:'
    print '(Exact results:  10   1000   101000001000)'
    cc10 = ClenshawCurtis(10, 0, 10)
    print cc10.quad(f), cc10.quad(g), cc10.quad(h)
    f10 = Fejer1(10, 0, 10)
    print f10.quad(f), f10.quad(g), f10.quad(h)
    print 'Composite quad for above cases:'
    print cq.quad(f), cq.quad(g), cq.quad(h)
    print

    def plm3(x):
        """
        Power-law integrand; integral = -1/x**2
        """
        return 2 * x**(-3.)

    def plm2(x):
        """
        Power-law integrand; integral = -1/x.
        """
        return x**(-2.)

    def plm1(x):
        """
        Power-law integrand; integral = ln(x).
        """
        return x**(-1.)

    def pl25(x):
        """
        Power-law integrand; integral = x**3.5.
        """
        return x**2.5 / 3.5

    cc = ClenshawCurtis(6, 1, 10)
    pl3 = PLMapQuad(3, cc.nodes, cc.wts)
    pl1 = PLMapQuad(1, cc.nodes, cc.wts)
    
    def roll(x):
        """
        Rolling power law, -1.5 for x<<10, -2.5 for x>>10.
        """
        return (x/10)**(-1.5) / (1. + x/10)
    
    l, u = 1., 10.
    cc = ClenshawCurtis(10, l, u)
    cc100 = ClenshawCurtis(100, l, u)
    rl, ru = roll(l), roll(u)
    indx = -log(ru/rl) / log(u/l)
    pl = PLMapQuad(indx, cc.nodes, cc.wts)
    print '10-pt, 100-pt CC and 10-pt CC with PL map for rolling power law:'
    print cc.quad(roll)
    print cc100.quad(roll)
    print pl.quad(roll)
    print
    
    def aexp(x):
        return exp(-.5*x)

    def rexp(x):
        """
        Rolling exponential, -.6 at 1, -.4 at 100
        """
        return exp((-.5 + .002*(x-50))*x)
    
    def rexp2(x):
        """
        Rolling exponential, -.6 at 1, -.4 at 10
        """
        return exp((-.5 + .02*(x-5))*x)
    
    eq = ExpMapQuad(-.5, cc.nodes, cc.wts)
    print 'ExpMapQuad (exact, quad):', (exp(-.5*100)-exp(-.5*1))/(-.5), eq.quad(aexp)
    