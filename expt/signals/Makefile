FC = g77
FFLAGS= -O -fno-automatic

# Test compilation of subroutine collection:
_cosmo.o:  _cosmo.f

# Build default .pyf file:
_cosmo.pyf: _cosmo.f
	f2py -m _cosmo _cosmo.f skip: dd_l_ig lbt_ig age_ig qromb trapzd polint : \
	    -h _cosmo.pyf --overwrite-signature

# Build the lib:
_cosmo.so: _cosmo.pyf _cosmo.f
	f2py -c _cosmo.pyf _cosmo.f
