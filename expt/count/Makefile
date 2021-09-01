#FFLAGS = -u
FC = gfortran
FFLAGS = -fno-automatic

io = /home/astrosun/loredo/fort/ioutil.a
nr = /home/astrosun/loredo/fort/dpnr.a

sdens: sdens.o on_off.o $(nr) $(io)
	f77 -o sdens sdens.o on_off.o $(nr) $(io)

signif: signif.o on_off.o $(nr) $(io)
	f77 -o signif signif.o on_off.o $(nr) $(io)

cbmpdf: cbmpdf.f _cbmlike.o $(nr) $(io)
	f77 -u -o cbmpdf cbmpdf.f cbmlike.o $(nr) $(io)

_cbmlike.o:  _cbmlike.f
	$(FC) $(FFLAGS) -c cbmlike.f

_cbmlike.pyf:  _cbmlike.f
	f2py --overwrite-signature -m _cbmlike -h _cbmlike.pyf _cbmlike.f\
	    skip: ppois lppois gammln :
	echo ''
	echo '*** Edit _cbmlike.pyf to add docstring blocks ***'

_cbmlike.so:  _cbmlike.f _cbmlike.pyf
	f2py -c _cbmlike.pyf _cbmlike.f

clean:
	rm -f _cbmlike.so _cbmlike.o