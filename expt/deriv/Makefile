# FC= g77
# FFLAGS= -O -fno-automatic

FC= gfortran
# -fno-automatic saves local vars (vs. using stack for small vars)
# -O is a gcc optimization option
FFLAGS= -fno-automatic -O

_obsinfo.o:  _obsinfo.f
	$(FC) -c _obsinfo.f

_obsinfo-raw.pyf:
	f2py -m _obsinfo _obsinfo.f skip: rid_ex : -h _obsinfo-raw.pyf --overwrite-signature

_obsinfo.so:  _obsinfo.pyf _obsinfo.f
	f2py -c _obsinfo.pyf _obsinfo.f

test_ninfo:  test_ninfo.f
	$(FC) $(FFLAGS) -o test_ninfo test_ninfo.f