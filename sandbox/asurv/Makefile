FC = g77
FFLAGS= -O -fno-automatic

# Test compilation of Eric's ASURV distribution:
asurv_lib.o:  asurv_lib.f

asurv:  asurv_main.f  asurv_lib.o
	g77 -o asurv asurv_main.f asurv_lib.o

# Build default .pyf file:
_asurvkm.pyf: _asurvkm.f
	f2py -m _asurvkm _asurvkm.f only: kmestm kmdif plestm quart : -h _asurvkm.pyf --overwrite-signature

# Build the KM est lib:
_asurvkm.so: _asurvkm.pyf _asurvkm.f
	f2py -c _asurvkm.pyf _asurvkm.f
