CFLAGS = -g -I. -I/home/laplace/include/python2.2 -I/home/laplace/include/python2.2/Numeric

#...............................................................................
fkeplermodule.so:  fkeplermodule.o
	ld -G fkeplermodule.o -o $@

fkeplermodule.o:  fkeplermodule.c

#...............................................................................
mkarraymodule.so:  mkarraymodule.o
	ld -G mkarraymodule.o -o $@

mkarraymodule.o:  mkarraymodule.c

#...............................................................................
vbamodule.so:  vbamodule.o
	ld -G vbamodule.o -o $@

vbamodule.o:  vbamodule.c
