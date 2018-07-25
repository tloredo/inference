import os, pickle
import os.path as path
from numpy import random

__all__ = ['save_rng', 'restore_rng']

MT_id = 'MT19937'  # NumPy's RNG as of 1.0.3

def save_rng(fname='.numpy-rng-state'):
    """
    Save the state of NumPy's RNG to a file in the CWD ('.numpy-rng-state' by
    default).  Backup the previous two saved states if present:  If the RNG
    state file exists (from a previous save), rename the previous one with a
    '.1' suffix.  If a '.1' file exists, rename it with a '.2' suffix.
    """
    state = random.get_state()
    if state[0] == MT_id:
        id, state = state[0], (state[1], state[2]) # (ID, (key, pos)) for MT
    else:
        raise RuntimeError('numpy.random using unrecognized RNG type!')
    if path.exists(fname):
        fname1 = fname + '.1'
        if path.exists(fname1):
            fname2 = fname + '.2'
            os.rename(fname1, fname2)
        os.rename(fname, fname1)
    ofile = open(fname, 'w')
    pickle.dump((id, state), ofile)
    ofile.close()
    
def restore_rng(fname='.numpy-rng-state', notify=True):
    """
    Restore the state of NumPy's RNG from the contents of a file in the CWD
    if the file exists; otherwise use (and save) the default initialization.
    The default file name is '.numpy-rng-state'.
    """
    if os.access(fname, os.R_OK | os.W_OK):
        rng_file = open(fname, 'r')
        id, state = pickle.load(rng_file)
        rng_file.close()
        if id == MT_id:
            # Note key is numpy,uint32 -> need repr() to see value.
            if notify:
                print('Recovered RNG state:  %s [%s %s ...] %i' %\
            	    (id, repr(state[0][0]), repr(state[0][1]), state[1]))
            random.set_state((id, state[0], state[1]))
        else:
            raise ValueError('Invalid ID for RNG in %s!' % fname)
    else:
        print('No accessible RNG status file; using (and saving) default initialization.')
        save_rng(fname)
