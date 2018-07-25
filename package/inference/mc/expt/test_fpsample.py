#! /usr/bin/env python

from numpy import zeros, array, rand, linspace, random
from _fpsample import equalprob, equalprobi
from _fpsample import prepwts, wtratios, SampfordPPS
from _fpsample import set_rng_state, get_rng_state, local_irand, local_rand

pool = zeros(10)
samp = equalprobi(5,10,pool)
print 'equalprobi:', samp

pool = range(10,20)
samp = equalprob(5,pool)
print 'equalprob:', samp
print pool

pool = array(pool)
samp = equalprob(5,pool)
print 'equalprob:', samp
print pool
print

wts = rand(6)
ind, sorted, cum = prepwts(wts)
print wts
print ind
print wts.take(ind)
print cum
print

wts = linspace(.1, 1., 1000)
ind, sorted, cumwt = prepwts(wts)
nsamp = 10
cumrat = wtratios(nsamp, sorted, cumwt[-1])
ntry, samp = SampfordPPS(nsamp, cumwt, cumrat)

state0 = random.get_state()
id, key, pos = state0
set_rng_state(key, pos)
print rand(), local_rand()
