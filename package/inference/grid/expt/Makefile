#FFLAGS = -u
FC = g77
FFLAGS = -fno-automatic
	
gridquad.o:  gridquad.f
	$(FC) $(FFLAGS) -c gridquad.f

gridquad.pyf:  gridquad.f
	f2py --overwrite-signature -m gridquad -h gridquad.pyf gridquad.f

gridquad.so:  gridquad.f gridquad.pyf
	f2py -c gridquad.pyf gridquad.f

