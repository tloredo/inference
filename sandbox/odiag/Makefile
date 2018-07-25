FFLAGS = -u -C
# FFLAGS = -u

IO = /home/lagrange/loredo/fort/ioutil.a
NR = /home/lagrange/loredo/fort/dpnr.a

#-----------------------------------------------------------------------
# Test nxtarc:
t_nxt:  t_nxt.o output.o $(NR) $(IO)
	f77 -o t_nxt t_nxt.o output.o $(NR) $(IO)

