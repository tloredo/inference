from numpy import sin, cos, zeros
import fkepler

fkepler.settol(1.e-5)

def v_rad(K, tau, e, long, T0, v0, times):
    """
    Calculate the radial (line-of-sight) velocity for a Keplerian orbit
    as a function of orbital parameters, COM velocity, and time.  The
    parameters are:
      K = velocity amplitude
      tau = period
      e = eccentricity
      long = longitude of periastron (radians)
      T0 = time of periastron crossing
      v0 = COM velocity
    The times argument can be either a single float (with the velocity
    returned as a float) or an array of floats (with the velocity
    returned as a corresponding array).

    tau, T0 and times should be in the same units (typically days).
    K can use a different unit of time (typically m/s).
    """
    fkepler.setup_Tp(tau, e, T0)
    if type(times)==float:
        c, s = fkepler.t2TA(times)  # Get cos, sin of true anomaly
    else:
        c, s = fkepler.vt2TA(times)
    return v0 + K*cos(long)*(e+c) + K*sin(long)*s
    
    
class Orbit:
    """
    Keplerian orbit object, with methods for calculating various orbit
	observables as a function of time:
	    EA_cs = eccentric anomaly (cosine & sine)
		TA_cs = true anomaly (cosine & sine)
		v_rad = radial (line-of-sight) velocity
		K1, K2 = terms of v_rad independent of (K, long), useful in
		         model fitting
	The observed radial velocity is
		v_rad = v0 + K*cos(long)*K1 - K*sin(long)*K2
	with v0 the COM velocity.
	
	For most efficient use if K1 and K2 are needed, create the Kepler instance,
		kep = Kepler(K, tau, e, long, T0)
	Then, if velocities are needed at a set of times, first calculate the
	true anomalies (ignore the return value; it's stored in the instance):
		kep.true_anom(times)
	Then call K1, K2 or v_rad *without* time arguments (the stored anomaly
	values will be used):
		K1 = kep.K1()
		K2 = kep.K2()
		v_rad = kep.v_rad()
	"""
    
    def __init__(self, K, tau, e, long, T0, v0):
        """Constructor, setting the velocity amplitude, period, eccentricity, 
		longitude of periastron, and time of periastron crossing for the orbit."""
        self.K, self.tau, self.e, self.long, self.T0 = K, tau, e, long, T0
        self.v0 = v0
        self.Kclong, self.Kslong = K*cos(long), K*sin(long)
        
    def ecc_anom(self, times):
        """Calculate the cosine & sine of the eccentric anomaly for a vector
		of observing times."""
        fkepler.setup_Tp(tau, e, T0)
        self.cs_EA = fkepler.vt2EA(times)
        return self.cs_EA
        
    def true_anom(self, times):
        """Calculate the cosine & sine of the true anomaly for a vector
		of observing times."""
        fkepler.setup_Tp(tau, e, T0)
        self.cs_TA = fkepler.vt2TA(times)
        return self.cs_TA
        
    def K1(self, times=None):
        """The quantity [e + cos(nu)] for nu=true anomaly."""
        if times is not None:
            self.true_anom(times)
        return self.e + self.cs_TA[:,0]
        
    def K2(self, times=None):
        """The quantity sin(nu) for nu=true anomaly."""
        if times is not None:
            self.true_anom(times)
        return self.cs_TA[:,1]
        
    def v_rad(self, times=None):
        """The line-of-sight velocity vs. time (including COM velocity)."""
        if times is not None:
            self.true_anom(times)
        K1 = self.K1()
        K2 = self.K2()
        return self.v0 + self.Kclong*K1 - self.Kslong*K2
        
    def phase_vel(self, n):
        """Return two arrays with data for plotting a velocity curve vs. phase.
		The first array will have n phases in [-.5, 1.5], the corresponding 
		velocities will be in the second array."""
        dt = 2./(n-1)
        phases = zeros(n, float)
        for i in range(n):
            phases[i] = -.5 + i*dt
        times = self.tau*phases
        v_rad = self.v_rad(times)
        return phases, v_rad
    
    def writeVelCurve(self, n, ofile=None, label=None):
        """Create/write a string with data for plotting a velocity curve vs. phase.
		The curve is written as n phases in [-.5, 1.5] and the corresponding 
		velocities."""
        dt = 2./(n-1)
        phases = zeros(n, float)
        for i in range(n):
            phases[i] = -.5 + i*dt
        times = self.tau*phases
        v_rad = self.v_rad(times)
        if label:
            s = label + '  %-6.4f  %-6.4f\n'
        else:
            s = '%-6.4f  %-6.4f\n'
        out = ''
        for i in range(n):
            out += s % (phases[i], v_rad[i])
            # print phases[i], times[i], self.cs_TA[i], v_rad[i]
        if ofile:
            ofile.write(out)
            ofile.write('\n')
        else:
            return out
            
            
if __name__ == '__main__':
    # Write a text file with data for a v_rad vs. phase plot;
    # seperately, use matplotlib to plot such data on the screen.
    import ioutils
    from scipy import array, pi
    v0 = 15.
    K = 10.
    tau = 17.1
    e = 0.6
    long = 0.4*pi
    T0 = 37.5
    orbit = Orbit(K, tau, e, long, T0, v0)
    ofile = ioutils.openText('kepler.dat')
    orbit.writeVelCurve(200, ofile, 'vcurve')
    ofile.close()
    
    import pylab as pl
    phase, vel = orbit.phase_vel(200)
    pl.plot(phase, vel)
    pl.xlabel('Orbital Phase')
    pl.ylabel('Radial Velocity (km/s)')
    pl.savefig('kepler.eps')
    pl.show()
