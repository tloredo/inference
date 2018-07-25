#FFLAGS = -u
FC = g77
FFLAGS = -fno-automatic

IO = /Users/loredo/fort/util/ioutil.a
MISC = /Users/loredo/fort/util/dpnr.a

#-----------------------------------------------------------------------
# Simulate arrival times:

events: events.o $(MISC) $(IO)
	$(FC) $(FFLAGS) -o events events.o $(MISC) $(IO)

#-----------------------------------------------------------------------
# Piecewise-constant model:

pcfreq: pcfreq.o pc.o peaks.o $(MISC) $(IO)
	$(FC) $(FFLAGS) -o pcfreq pcfreq.o pc.o peaks.o $(MISC) $(IO)

pcfreq2: pcfreq2.o pc.o peaks.o $(MISC) $(IO)
	$(FC) $(FFLAGS) -o pcfreq2 pcfreq2.o pc.o peaks.o $(MISC) $(IO)

pcfreqd: pcfreqd.o pcd.o peaks.o $(MISC) $(IO)
	$(FC) $(FFLAGS) -o pcfreqd pcfreqd.o pcd.o peaks.o $(MISC) $(IO)

shape: shape.o pc.o $(MISC) $(IO)
	$(FC) $(FFLAGS) -o shape shape.o pc.o $(MISC) $(IO)

pcf-test: pcf-test.o pcd.o peaks.o $(MISC) $(IO)
	$(FC) $(FFLAGS) -o pcf-test pcf-test.o pcd.o peaks.o $(MISC) $(IO)


#-----------------------------------------------------------------------
# Python support:

_pcperiodic.o: _pcperiodic.f
	$(FC) $(FFLAGS) _pcperiodic.f

_pcperiodic.pyf: _pcperiodic.f
	f2py -m _pcperiodic -h _pcperiodic.pyf skip: factln gammln : _pcperiodic.f

_pcperiodic.so: _pcperiodic.pyf _pcperiodic.f
	f2py -c _pcperiodic.pyf _pcperiodic.f

