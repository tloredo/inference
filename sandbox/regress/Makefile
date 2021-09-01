FC = g77
FFLAGS= -O -fno-automatic

# Test compilation of subroutine collection from Matt's distribution:
_regress.o:  _regress.f

# Build default .pyf file:
_regress.pyf: _regress.f
	f2py -m _regress _regress.f skip: class nclass : -h _regress.pyf --overwrite-signature

# Build the lib:
_regress.so: _regress.pyf _regress.f
	f2py -c _regress.pyf _regress.f
