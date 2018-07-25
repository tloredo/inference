"""
Make cubature rules

Generate Python code implementing cubature rules from Ronald Cool's
online Encyclopedia of Cubature Formulas.  Using textual rule
descriptions copied-and-pasted from the web site, this script
processes the descriptions and generates Python code that takes a
function as an argument and returns the cubature of the function.

The generated code is in "_cub2d_auto.py".
"""

import re, pprint
from string import Template
from numpy import array, sum, sin, cos, pi, zeros, linspace, array2string
from numpy.core import arrayprint

def stripped_lines(text):
    """
    Return a generator that supplies stripped lines from a block of text,
    ignoring blank or comment ('#') lines, and stripping any trailing commas.
    """
    line_iter = iter(text.split('\n'))
    while 1:
        next = line_iter.next().strip()
        if next:
            if next[0] == '#':
                continue
            elif next[-1] == ',':
                yield next[:-1]
            else:
                yield next
        else:
            continue


def event_handler(regexp):
    """
    Decorator associating a regexp with an event handler.  The
    regexp will be compiled first.
    """
    def decorator(func):
        func.event_re = regexp
        func.event_cre = re.compile(regexp)
        return func
    return decorator

# This isn't quite a FSM, since there aren't really ordered transitions
# (except within the generator and weight processing).

class Processor(object):

    def __init__(self):
        self.handlers = []
    
    def add_handler(self, handler):
        """
        If the handler's regexp matches the current line, if there are no
        named groups in the match, call handler(line_gen, result), otherwise
        call handler(line_gen, result, **dct) where dct is the groupdict of
        the match, and result is the result of the previous processing.
        """
        if hasattr(handler, 'event_re'):
            self.handlers.append(handler)
        else:
            raise ValueError('Object is not an event handler!')
    
    def add_handlers(self, hlist):
        for h in hlist:
            self.add_handler(h)
    
    def process(self, line_gen, result=None):
        while True:
            try:
                line = line_gen.next()
            except StopIteration:
                break
            # print line
            for handler in self.handlers:
                m = handler.event_cre.match(line)
                # print ' ', handler.event_re, '==>', m
                if m:
                    dct = m.groupdict()
                    if dct:
                        handler(line_gen, result, **dct)
                    else:
                        handler(line_gen, result)
                    break
            if not m:
                raise ValueError('Unrecognized key: %s' % line)


# These are for experimenting at an interactive prompt:
#
# Match a keyword at the start of the line, with arg matching subsequent
# text after any whitespace, and ignoring trailing whitespace.
region_re = re.compile(r'^Region:\s*(?P<arg>.*?)\s*$')
dimension_re = re.compile(r'^Dimension:\s*(?P<arg>.*?)\s*$')
degree_re = re.compile(r'^Degree:\s*(?P<arg>.*?)\s*$')
points_re = re.compile(r'^Points:\s*(?P<arg>.*?)\s*$')
points_circ_re = re.compile(r'^Points:\s*(?P<arg>.*?)\s*\[\s*(?P<geom>.*?)\s*\]\s*$')
struct_re = re.compile(r'^Structure:\s*(?P<arg>.*?)\s*$')
rulestruct_re = re.compile(r'^Rule struct:\s*(?P<arg>.*?)\s*$')

# Match a keyword at the start of the line, with arg matching subsequent
# subsequent text in brackets after any whitespace, and ignoring trailing 
# whitespace.
gen_re = re.compile(r'^Generator:\s*\[\s*(?P<arg>.*?)\s*\]\s*$')

# Match a keyword at the start of the line.
wt_re = re.compile(r'^Corresponding weight:\s*$')
ID_re = re.compile(r'^Rule\s+(?P<arg>.*?),.*?$')


# Define event handlers for rule elements, using the above regexps.

@event_handler(r'^Rule\s+(?P<arg>.*?),.*?$')
def id_handler(line_gen, result, arg):
    result.id = arg
    result.ident = arg.replace('-', '_') # rule id as Python identifier

@event_handler(r'^Region:\s*(?P<arg>.*?)\s*$')
def region_handler(line_gen, result, arg):
    result.region = arg

@event_handler(r'^Dimension:\s*(?P<arg>.*?)\s*$')
def dimension_handler(line_gen, result, arg):
    result.dimen = int(arg)

@event_handler(r'^Degree:\s*(?P<arg>.*?)\s*$')
def degree_handler(line_gen, result, arg):
    result.degree = int(arg)

# Note we should try to match points_circ before points.
@event_handler(r'^Points:\s*(?P<arg>.*?)\s*\[\s*(?P<geom>.*?)\s*\]\s*$')
def points_circ_handler(line_gen, result, arg, geom):
    result.npts = int(arg)
    result.circ_geom = geom

@event_handler(r'^Points:\s*(?P<arg>.*?)\s*$')
def points_handler(line_gen, result, arg):
    result.npts = int(arg)

@event_handler(r'^Structure:\s*(?P<arg>.*?)\s*$')
def struct_handler(line_gen, result, arg):
    result.struct = arg
    if arg == 'Circular symmetric':
        caption = line_gen.next()
        if caption.split()[0] != 'i':
            raise ValueError('Invalid circular symmetric rule spec!')
        npts = 0
        result.circ_spec = []
        while npts < result.npts:
            row = line_gen.next().split()
            n = int(row[1])
            r, w = float(row[2]), float(row[3])
            result.circ_spec.append((n, r, w))
            npts += n

@event_handler(r'^Rule struct:\s*(?P<arg>.*?)\s*$')
def rulestruct_handler(line_gen, result, arg):
    result.rule_struct = array([int(s) for s in arg.split()])
    result.ngens = result.rule_struct.sum()
    result.gen_types = []
    result.gens = []
    result.gen_wts = []

@event_handler(r'^Generator:\s*\[\s*(?P<arg>.*?)\s*\]\s*$')
def generator_handler(line_gen, result, arg):
    result.gen_types.append(arg)
    result.gens.append( eval(line_gen.next()) )

@event_handler(r'^Corresponding weight:\s*$')
def weight_handler(line_gen, result):
    result.gen_wts.append( eval(line_gen.next()) )


def xcomb(items, n):
    """
    Generator returning all combinations of n objects selected
    (without replacement, all permutations) from the sequence items.
    
    This is Ulrich Hoffmann's recipe from the ASPN Python Cookbook:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/190465
    """
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xcomb(items[:i]+items[i+1:],n-1):
                yield [items[i]]+cc

def xsubsets(items, n):
    """
    Generator returning all *unique* combinations of n objects selected
    from the sequence items (i.e., different permutations of the same
    set of selected items are not returned).
    
    This is Ulrich Hoffmann's recipe from the ASPN Python Cookbook:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/190465
    """
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xsubsets(items[i+1:],n-1):
                yield [items[i]]+cc

def xperm(items):
    """
    Generator returning all permutations of objects from the sequence items.
    
    This is Ulrich Hoffmann's recipe from the ASPN Python Cookbook:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/190465
    
    For a non-lazy algorithm, use numpy.random.permutation.
    """
    return xcomb(items, len(items))

def flip_signs(seq):
    """
    Generator taking a sequence of non-negative values and returning
    all unique sequences with components of all signs.
    
    This algorithm is adapted from part of Fortran function FLSM in Alan
    Genz's ADAPT.  We go through the sequence from start to end, flipping
    signs, and yielding a value whenever an element becomes negative.  After
    a yield, the sign-flipping starts again from the beginning.
    
    E.g., for input [1, 1]:
    [1, 1] -> yield (on start)
    [-1, 1] -> yield
    [1, 1] -> [1, -1] -> yield
    [-1, -1] -> yield
    """
    yield seq
    i = 0
    while i < len(seq):
        if seq[i] != 0.:
            seq[i] = -seq[i]
        if seq[i] < 0.:
            yield seq
            i = 0
        else:
            i += 1

class RuleMaker(object):

    ruleProcessor = Processor()
    ruleProcessor.add_handlers([
        id_handler, region_handler, dimension_handler, degree_handler,
        points_circ_handler, points_handler, struct_handler, rulestruct_handler,
        generator_handler, weight_handler
        ])

    rule_template = Template("""
${id}_absc = array(
    ${absc} )

${id}_wts = array(
    ${wts} )

def ${name}(func):
    cub = 0.
    for absc, wt in zip(${id}_absc, ${id}_wts):
        cub += wt * func(*absc)
    return cub

""")


    def __init__(self, text):
        # Set attributes via processing the text description.
        self.ruleProcessor.process(stripped_lines(text), self)
        if self.struct == 'Fully symmetric':
            self._construct_fully_sym()
        elif self.struct == 'Circular symmetric':
            self._construct_circ_sym()
        else:
            raise ValueError('Unrecognized rule structure!')

    def _construct_fully_sym(self):
        self.abscissas = []
        self.weights = []
        for type, gen, wt in zip(self.gen_types, self.gens, self.gen_wts):
            if type == 'Origin':
                self.abscissas.append(array(gen))
                self.weights.append(wt)
            elif type == 'Cornerpoints of the unit-cube':
                flipper = flip_signs(list(gen))
                for coords in flipper:
                    self.abscissas.append(array(coords))
                    self.weights.append(wt)
            elif type == 'Fully symmetric':
                # Some coordinate values may be repeated, so we'll use a 
                # dict to detect permutation uniqueness.
                d = {}
                # Go through all permutations of the generator coordinates.
                perms = xperm(list(gen))
                for perm in perms:
                    t = tuple(perm)
                    if d.has_key(t):
                        continue
                    d[t] = None  # placeholder to mark the permutation
                    # Now go through all distinct sign flips.
                    flipper = flip_signs(list(perm))
                    for coords in flipper:
                        self.abscissas.append(array(coords))
                        self.weights.append(wt)
        print 'Constructed fully symmetric rule with', len(self.abscissas),\
            'points.'
        if len(self.abscissas) != self.npts:
            raise RuntimeError('Mismatch with expected # of points in rule!')
        self.abscissas = array(self.abscissas)
        self.weights = array(self.weights)

    def _construct_circ_sym(self):
        if self.dimen != 2:
            raise ValueError('Circularly symmetric rules valid only for 2-D!')
        self.abscissas = []
        self.weights = []
        for n, r, wt in self.circ_spec:
            if n == 1:
                if r != 0.:
                    raise ValueError('Nonzero radius in 1-pt circular rule!')
                self.abscissas.append(zeros((self.dimen), float))
                self.weights.append(wt)
            else:
                # Note we need a grid over n+1 angles; the 1st and last
                # points correspond so we drop the last.
                angles = linspace(0., 2*pi, n+1)[:-1]
                xvals, yvals = r*cos(angles), r*sin(angles)
                for x, y in zip(xvals, yvals):
                    self.abscissas.append(array([x,y]))
                    self.weights.append(wt)
        print 'Constructed circularly symmetric rule with', \
            len(self.abscissas), 'points.'
        if len(self.abscissas) != self.npts:
            raise RuntimeError('Mismatch with expected # of points in rule!')
        self.abscissas = array(self.abscissas)
        self.weights = array(self.weights)

    def cubature(self, func):
        sum = 0
        for absc, wt in zip(self.abscissas, self.weights):
            sum += wt * func(*absc)
        return sum

    def write_rule(self, ofile):
        """
        Write Python code implementing the rule to ofile.
        """
        # name = Template('cub_${dimen}d_${id}')
        # name.substitute({'id' : self.ident, 'dimen' : self.dimen})
        # Write 16 digits after the decimal - what the Encyc provided.
        # This is trickier than it should be!
        # numpy's array2string dumbly handles precision requests, losing
        # precision by insisting on %f format.  Temporarily replace its 
        # float formatter with something that uses %e.
        tmp = arrayprint._floatFormat
        arrayprint._floatFormat = float_e_format
        absc = array2string(self.abscissas, 80, separator=', ', prefix='    ')
        wts = array2string(self.weights, 80, separator=', ', prefix='    ')
        arrayprint._floatFormat = tmp
        ofile.write(self.rule_template.substitute(
            {'name' : self.ident, 'id' : self.ident,
             'absc' : absc, 'wts' : wts}))

def float_e_format(data, precision, suppress_small, sign = 0):
    return '%25.16e'

# Rule descriptions, cut-and-pasted from Cool's Encyclopedia of Cubature
# Formulas web site:

fully_sym_desc = [
"""Rule e2r2-3-4a, quality P
Region: Entire space with w = exp(r^2) 
Dimension: 2 
Degree: 3 
Points: 4 
Structure: Fully symmetric 
Rule struct: 0 1 0 0 
Generator: [ Fully symmetric ] 
( 1., 0., ) 
Corresponding weight: 
0.78539816339744830,
""",
"""
Rule e2-r2-3-4b, quality P
Region: Entire space with w = exp(r^2) 
Dimension: 2 
Degree: 3 
Points: 4 
Structure: Fully symmetric 
Rule struct: 0 0 1 0 
Generator: [ Fully symmetric ] 
( 0.70710678118654752, 0.70710678118654752, ) 
Corresponding weight: 
0.78539816339744830,
""",
"""
Rule e2r2-5-9a, quality P
Region: Entire space with w = exp(r^2) 
Dimension: 2 
Degree: 5 
Points: 9 
Structure: Fully symmetric 
Rule struct: 1 1 1 0 
Generator: [ Cornerpoints of the unit-cube ] 
( 1., 1., ) 
Corresponding weight: 
0.19634954084936207,
Generator: [ Origin ] 
( 0., 0., ) 
Corresponding weight: 
1.5707963267948966,

Generator: [ Fully symmetric ] 
( 1.4142135623730950, 0., ) 
Corresponding weight: 
0.19634954084936207,
""",
"""
Rule e2r2-5-9b, quality P
Region: Entire space with w = exp(r^2) 
Dimension: 2 
Degree: 5 
Points: 9 
Structure: Fully symmetric 
Rule struct: 1 0 0 1 
Generator: [ Origin ] 
( 0., 0., ) 
Corresponding weight: 
1.5707963267948966,
Generator: [ Fully symmetric ] 
( 1.3065629648763765, 0.54119610014619698, )
Corresponding weight:
0.19634954084936207,
""",
"""
Rule e2r2-7-12, quality P, minimum N
Region: Entire space with w = exp(r^2) 
Dimension: 2 
Degree: 7 
Points: 12 
Structure: Fully symmetric 
Rule struct: 0 1 2 0 
Generator: [ Fully symmetric ]
( 1.7320508075688772, 0., )
Corresponding weight:
0.087266462599716478,

Generator: [ Fully symmetric ]
( 0.53523313465963489, 0.53523313465963489, )
Corresponding weight:
0.66127983844512042,

Generator: [ Fully symmetric ]
( 1.4012585384440735, 1.4012585384440735, )
Corresponding weight:
0.036851862352611409,
""",
"""
Rule e2r2-9-21, quality P
Region: Entire space with w = exp(r^2) 
Dimension: 2 
Degree: 9 
Points: 21 
Structure: Fully symmetric 
Rule struct: 1 0 1 2 
Generator: [ Origin ]
( 0., 0., )
Corresponding weight:
1.0471975511965977,

Generator: [ Fully symmetric ]
( 1.5381890013208515, 1.5381890013208515, )
Corresponding weight:
0.012238016879947463,

Generator: [ Fully symmetric ]
( 0.43091398228826897, 1.0403183802565386, )
Corresponding weight:
0.24426215416421333,

Generator: [ Fully symmetric ]
( 0.54093734825794375, 2.1069972930282899, )
Corresponding weight:
0.011418225194962370,
""",
"""
Rule e2r2-11-28a, quality P
Region: Entire space with w = exp(r^2) 
Dimension: 2 
Degree: 11 
Points: 28 
Structure: Fully symmetric 
Rule struct: 0 3 0 2 
Generator: [ Fully symmetric ]
( 2.7578163962570077, 0., )
Corresponding weight:
#10 ^ -4 x 8.1766458176754177,
8.1766458176754177E-4,

Generator: [ Fully symmetric ]
( 1.7320508075688772, 0., )
Corresponding weight:
0.043633231299858239,

Generator: [ Fully symmetric ]
( 0.62805153015975584, 0., )
Corresponding weight:
0.53732552144981741,

Generator: [ Fully symmetric ]
( 1.2247448713915890, 2.1213203435596425, )
Corresponding weight:
#10 ^ -3 x 3.6361026083215199,
3.6361026083215199E-3,

Generator: [ Fully symmetric ]
( 0.70710678118654752, 1.2247448713915890, )
Corresponding weight:
0.098174770424681038,
""",
"""
Rule e2r2-13-36, quality P
Region: Entire space with w = exp(r^2) 
Dimension: 2 
Degree: 13 
Points: 36 
Structure: Fully symmetric 
Rule struct: 0 3 2 2 
Generator: [ Fully symmetric ]
( 2.3589322619806809, 0., )
Corresponding weight:
#10 ^ -3 x 4.0553465034375560,
4.0553465034375560E-3,

Generator: [ Fully symmetric ]
( 1.3764866796963506, 0., )
Corresponding weight:
0.11962714138576389,

Generator: [ Fully symmetric ]
( 0.54608441525854036, 0., )
Corresponding weight:
0.44548945821552661,

Generator: [ Fully symmetric ]
( 1.9855313399178849, 1.9855313399178849, )
Corresponding weight:
#10 ^ -4 x 5.0912352191625667,
5.0912352191625667E-4,

Generator: [ Fully symmetric ]
( 0.84624998844801990, 0.84624998844801990, )
Corresponding weight:
0.18135069243713553,

Generator: [ Fully symmetric ]
( 0.99625342617276321, 3.0304360653685791, )
Corresponding weight:
#10 ^ -5 x 7.5744694043053037,
7.5744694043053037E-5,

Generator: [ Fully symmetric ]
( 0.93182402277619985, 1.7792467389053190, )
Corresponding weight:
0.017107455972791175,
"""
]

circ_sym_desc = [
"""
Rule e2r2-9-19a, quality P
Region: Entire space with w = exp(r^2)
Dimension: 2 
Degree: 9 
Points: 19 [ 6-gon ] 
Structure: Circular symmetric
i	r(i)	w(i)
1	6	2.E0 	2.4543692606170259E-2
2	6	1.1024690870412086E0 	3.2447320274634297E-1
3	6	2.3752273089667399E0 	5.5031089588349583E-3
4	1	0.E0 	1.0144726277217040E0
""",
"""
Rule e2r2-11-28c, quality P
Region: Entire space with w = exp(r^2)
Dimension: 2 
Degree: 11 
Points: 28 [ 4-gon ] 
Structure: Circular symmetric 
i	r(i)	w(i)
1	8	2.2360679774997896E0 	7.5398223686155037E-3
2	4	1.3765979318520745E0 	1.5137267598023796E-1
3	4	2.6825150931410771E0 	1.0294743039068897E-3
4	4	6.1791116665497929E-1 	5.2680403237854616E-1
5	4	1.5261898960020960E0 	9.0730267409198283E-2
6	4	2.9220203379272203E0 	3.8206858832799763E-4
""",
"""
Rule e2r2-13-41, quality P
Region: Entire space with w = exp(r^2)
Dimension: 2 
Degree: 13 
Points: 41 [ 10-gon ] 
Structure: Circular symmetric
i	r(i)	w(i)
1	10	2.4494897427831780E0 	2.4240684055476799E-3
2	10	9.2897653516474092E-1 	1.9404444552734862E-1
3	10	1.7162505238049738E0 	4.4486553170629952E-2
4	10	3.0727002354041088E0 	1.1853582819051242E-4
5	1	0.E0 	7.3085662427262551E-1
"""
]


fully_sym_rules = {}
for desc in fully_sym_desc:
    rule = RuleMaker(desc)
    fully_sym_rules[rule.ident] = rule
print 'Created', len(fully_sym_rules.keys()), 'fully-symmetric cubature rules.'

circ_sym_rules = {}
for desc in circ_sym_desc:
    rule = RuleMaker(desc)
    circ_sym_rules[rule.ident] = rule
print 'Created', len(circ_sym_rules.keys()), 'circular-symmetric cubature rules.'

# Write Python code implementing the rules to a file.

ofname = '_cub2d_auto.py'
ofile = open(ofname, 'w')
ofile.write(
"""# NOTE:  This file is machine-generated; do not edit!
#
# If changes are necessary, modify make_ECF_rules.py (in integrate/expt in the 
# developer distribution) and re-run it to produce a new version of this file.

from numpy import array

""")

# Write out an __all__ that accesses only functions and the rule_dict.
names = fully_sym_rules.keys() + circ_sym_rules.keys()
all = names + ['rule_dict']
ofile.write('__all__ = %s\n\n' % pprint.pformat(all))

# Write out the rules.
for rule in fully_sym_rules.values():
    rule.write_rule(ofile)

for rule in circ_sym_rules.values():
    rule.write_rule(ofile)

# Write out a dictionary of name:function to simplify access.
ofile.write('names = %s\n\n' % pprint.pformat(names))
ofile.write("""
rule_dict = {}
for name in names:
    rule_dict[name] = ( eval(name), eval(name+'_absc'), eval(name+'_wts') )
""")

ofile.close()
print 'Wrote rules as a Python module to "%s".' % ofname
