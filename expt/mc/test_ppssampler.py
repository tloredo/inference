
from numpy import zeros, array, rand, linspace, random, Float, ones, arange, sum
from _ppssampler import set_rng_state, get_rng_state, local_irand, local_rand
from _ppssampler import equalprob, equalprobi
from population import *
from messages import restart, elapsedcpu

# Check passing RandomKit state back and forth with numpy.
print '*** Checking RandomKit state maintenance ***'
state0 = random.get_state()
id, key, pos = state0
set_rng_state(id, key, pos)
print 'Should be same:', rand(), local_rand()
for i in range(100):
    local_rand()
print 'Should be dif: ', rand(), local_rand()
random.set_state(get_rng_state())
print 'Should be same:', rand(), local_rand()
print

# Check equal-weight samplers.
print '*** Checking equiprobability samplers ***'
pool = zeros(10)  # Workspace
samp = equalprobi(5,10,pool)
print 'equalprobi:', samp

# Check equalprob with a list input.
pool = range(10,20)
samp = equalprob(5,pool)
print 'equalprob:', samp
# Check the original pool is unchanged.
print pool

# Check it again with an array.
pool = array(pool)
samp = equalprob(5,pool)
print 'equalprob:', samp
print pool
print

# Check PPS.
print '*** Checking Population samplers ***'
p5 = Population(['a', 'b', 'c', 'd', 'e'])
wts = array([1., .5, .2, 2.5, 1., .5, .7, .8, .7])
p1 = Population(weights=wts)

print 'With replacement (5-table):'
prob2 = wts[2]/sum(wts)
prob3 = wts[3]/sum(wts)
sample = p1.sample(10000)
n2 = sample.count(2)
n3 = sample.count(3)
print 'Frac expected, found of 2:', prob2, n2/10000.
print 'Frac expected, found of 3:', prob3, n3/10000.

print 'Without replacement (Sampford PPS):'

def test_nsamp(pop, nsamp, nrep):
    probs = pop.weights/sum(pop.weights)
    counts = zeros(pop.npopn,Float)
    restart()
    for i in range(nrep):
        if pop.npopn>10000 and (i+1)%10==0:
            print 'Repetition', i+1, 'of', nrep,'...'
        sample = pop.subset_pps(nsamp)
        for n in range(pop.npopn):
            if n in sample: counts[n] += 1
    print 'CPU time:', elapsedcpu()
    counts /= nsamp
    for n,c in enumerate(counts):
        print 'Frac expected, found of %i:  %f  %f' % (n, probs[n], c/nrep)
        if n>=20:
            print '...'
            break

def test_nsamp5(pop, nsamp, nrep):
    probs = pop.weights/sum(pop.weights)
    counts = zeros(pop.npopn,Float)
    restart()
    for i in range(nrep):
        if pop.npopn>10000 and (i+1)%10==0:
            print 'Repetition', i+1, 'of', nrep,'...'
        sample = pop.subset_pps5(nsamp)
        for n in range(pop.npopn):
            if n in sample: counts[n] += 1
    print 'CPU time:', elapsedcpu()
    counts /= nsamp
    for n,c in enumerate(counts):
        print 'Frac expected, found of %i:  %f  %f' % (n, probs[n], c/nrep)
        if n>=20:
            print '...'
            break

print 'Direct PPS:'
test_nsamp(p1, 3, 10000)
print '\n5-table PPS:'
test_nsamp5(p1, 3, 10000)

# 1-d population with same weights.
p1d = Population1D(linspace(10,90,9), wts)

# A large population.
wts = rand(10000)
large = Population(weights=wts)
