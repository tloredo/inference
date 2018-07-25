#FFLAGS = -u
FC = g77
#FFLAGS = -fno-automatic

fisher_coinc.o:  fisher_coinc.f
	$(FC) $(FFLAGS) -c fisher_coinc.f

fisher_coinc.pyf:  fisher_coinc.f
	f2py --overwrite-signature -m fisher_coinc -h fisher_coinc.pyf fisher_coinc.f

fisher_coinc.so:  fisher_coinc.f fisher_coinc.pyf
	f2py -c fisher_coinc.pyf fisher_coinc.f
