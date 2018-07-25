c=======================================================================
c	tnolike.f:	Subprograms for TNO Distribution (TNOD)
c			program that read in survey data and
c			define the TNO dist'n likelihood function.
c
c	Created 23 Mar 98 by TJL
c
c=======================================================================
c-----------------------------------------------------------------------
	block data tnolike_bd

c---  	Block data subprogram initializing models and data.

	include 'tnolike.cb'
	include 'tnodata.cb'
	data curmod/0/
	data npars/2, 1, 3, 3, 4, 2, 3/
	data nsurveys/0/
	data lcut/1.e-5/
	data s_area(0)/1.d0/

	end

c-----------------------------------------------------------------------
	subroutine reset_data

c---	"RESET DATA"
c
c	Reset global variables to correspond to NO data present.

c+++  Globals:
	include 'tnodata.cb'

	nsurveys = 0

	return
	end

c-----------------------------------------------------------------------
	subroutine set_quad (n, mlo, mhi, eps)

c---	"SET QUADrature rules"
c
c	Define abscissas and weights for a Gauss-Hermite
c	quadrature over the uncertain object magnitude,
c	and a range and accuracy for Romberg quadrature over the
c	detection efficiency.

c+++  Arguments:
	integer n
	real*8 mlo, mhi, eps

c+++  Globals:
	include 'quad.cb'

	if (n .gt. nmqmax) call quit('Too many m quad points!')
	nmq = n
	if (nmq .eq. 1) then
	    mq_absc(1) = 0.
	    mq_wts(1) = 1.
	else
	    call gauher(mq_absc, mq_wts, nmq)
	endif

	m_lo = mlo
	m_hi = mhi
	delm = mhi - mlo
	rq_eps = eps

	return
	end

c-----------------------------------------------------------------------
	subroutine read_survey (lun, fname, sname, n)

c---	"READ SURVEY"
c
c	Read in the data defining a survey from file fname, 
c	using logical unit lun.  Return # of objects in survey, n.

c+++  Arguments:
	integer lun, n
	character*(*) fname, sname

c+++  Globals:
	include 'tnodata.cb'
	include 'quad.cb'
	real*8 pars(epmax), opars(2)
	common /eparcb/ pars, opars

c+++  Locals:
	real*8 one
	parameter (one = 1.)
	character*20 type
	character*72 args
	integer i
	real*8 sys
	logical ok

c+++  Functions:
	real*8 sl_scl
	logical rdpfnd, match

c---  First make sure we can handle another one!
	if (nsurveys .ge. nsmax) call quit('Too many surveys!')

c---  Read in the next survey.
	nsurveys = nsurveys + 1
	if (nsurveys .gt. nsmax)
     *      call quit('Too many surveys!')
	call opar (lun,fname)
	call rdpw('name',s_name(nsurveys))
	sname = s_name(nsurveys)

	call rdpw('survey.type',type)
	call lowify(type)
	if (match('list',type)) then
	    s_type(nsurveys) = 1
	else if (match('number',type)) then
	    s_type(nsurveys) = 2
	else if (match('on-off',type)) then
	    s_type(nsurveys) = 3
	else
	    call quit('Unknown survey type!')
	endif

	call rdpw('efficiency.type',type)
	call lowify(type)
	if (match('pritchett',type)) then
	    eta_type(nsurveys) = 1
	    call rdpda('mc.alpha',pars,2)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    if (eta_par(1,nsurveys) .lt. m_lo) 
     *          call quit('Eta mc below quadrature limit!')
	    if (eta_par(1,nsurveys) .gt. m_hi) 
     *          call quit('Eta mc above quadrature limit!')
	else if (match('heaviside',type)) then
	    eta_type(nsurveys) = 2
	    call rdpda('A.mcut',pars,2)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    if (eta_par(2,nsurveys) .lt. m_lo) 
     *          call quit('Eta cutoff below quadrature limit!')
	    if (eta_par(2,nsurveys) .gt. m_hi) 
     *          call quit('Eta cutoff above quadrature limit!')
	else if (match('tanh',type)) then
	    eta_type(nsurveys) = 3
	    call rdpda('A.mc.dm',pars,3)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    eta_par(3,nsurveys) = pars(3)
	    if (eta_par(2,nsurveys) .lt. m_lo) 
     *          call quit('Eta mc below quadrature limit!')
	    if (eta_par(2,nsurveys) .gt. m_hi) 
     *          call quit('Eta mc above quadrature limit!')
	else if (match('topc',type)) then
	    eta_type(nsurveys) = 4
	    call rdpda('A.mc.dm',pars,3)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    eta_par(3,nsurveys) = pars(3)
	    if (eta_par(2,nsurveys) .lt. m_lo) 
     *          call quit('Eta mc below quadrature limit!')
	    if (eta_par(2,nsurveys) .gt. m_hi) 
     *          call quit('Eta mc above quadrature limit!')
	else if (match('linear',type)) then
	    eta_type(nsurveys) = 5
	    call rdpda('m1.m2',pars,2)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    eta_par(3,nsurveys) = pars(2) - pars(1)
	    if (eta_par(1,nsurveys) .lt. m_lo) 
     *          call quit('Eta break below quadrature limit!')
	    if (eta_par(2,nsurveys) .gt. m_hi) 
     *          call quit('Eta cutoff above quadrature limit!')
	else if (match('tanh_cut',type)) then
	    eta_type(nsurveys) = 6
	    call rdpda('A.mc.dm.cut.w',pars,5)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    eta_par(3,nsurveys) = pars(3)
	    eta_par(4,nsurveys) = pars(4)
	    eta_par(5,nsurveys) = pars(5)
	    if (eta_par(2,nsurveys) .lt. m_lo) 
     *          call quit('Eta mc below quadrature limit!')
	    if (eta_par(2,nsurveys) .gt. m_hi) 
     *          call quit('Eta mc above quadrature limit!')
	else if (match('tanhp1',type)) then
	    eta_type(nsurveys) = 7
	    call rdpda('A.mc.dm',pars,3)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    eta_par(3,nsurveys) = pars(3)
	    if (eta_par(2,nsurveys) .lt. m_lo) 
     *          call quit('Eta mc below quadrature limit!')
	    if (eta_par(2,nsurveys) .gt. m_hi) 
     *          call quit('Eta mc above quadrature limit!')
	else if (match('dbltanh',type)) then
	    eta_type(nsurveys) = 8
	    call rdpda('A.mc.dm.dm2',pars,4)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    eta_par(3,nsurveys) = pars(3)
	    eta_par(4,nsurveys) = pars(4)
	    if (eta_par(2,nsurveys) .lt. m_lo) 
     *          call quit('Eta mc below quadrature limit!')
	    if (eta_par(2,nsurveys) .gt. m_hi) 
     *          call quit('Eta mc above quadrature limit!')
	else if (match('erfc2',type)) then
	    eta_type(nsurveys) = 9
	    call rdpda('A.mc.dm',pars,3)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    eta_par(3,nsurveys) = pars(3)
	    if (eta_par(2,nsurveys) .lt. m_lo) 
     *          call quit('Eta mc below quadrature limit!')
	    if (eta_par(2,nsurveys) .gt. m_hi) 
     *          call quit('Eta mc above quadrature limit!')
	else if (match('erfc-rt2',type)) then
	    eta_type(nsurveys) = 10
	    call rdpda('A.mc.dm',pars,3)
	    eta_par(1,nsurveys) = pars(1)
	    eta_par(2,nsurveys) = pars(2)
	    eta_par(3,nsurveys) = pars(3)
	    if (eta_par(2,nsurveys) .lt. m_lo) 
     *          call quit('Eta mc below quadrature limit!')
	    if (eta_par(2,nsurveys) .gt. m_hi) 
     *          call quit('Eta mc above quadrature limit!')
	else
	    call quit('Unknown efficiency type!')
	endif

	call rdpd('sq.degrees',s_area(nsurveys))


c---  Read in object data for "list" type surveys.
	if (s_type(nsurveys) .eq. 1) then
	    call rdpi('num.objects',nobjects(nsurveys))
	    call rdpd('sys.err',sys)
	    do 20 i=1, nobjects(nsurveys)
	        call rdpda('object',opars,2)
	        if (.not. rdpfnd()) call quit('Missing object data!')
	        obj_m(i,nsurveys) = opars(1)
	        obj_msig(i,nsurveys) = sqrt(opars(2)**2+sys**2)
20	    continue

c---  Read in number seen for "number" type surveys.
	else if (s_type(nsurveys) .eq. 2) then
	    call rdpi('num.objects',nobjects(nsurveys))
	    if (nobjects(nsurveys) .gt. nomax)
     *          call quit('Too many objects!')

c---  Read in bin data for "on-off" type surveys.  Calculate
c---  log scale factors.  Also do V to R conversion.
	else if (s_type(nsurveys) .eq. 3) then
	    call rdpi('num.bins',nobjects(nsurveys))
	    call rdpd('bin.size',dm(nsurveys))
	    call rdpd('V-R',sys)
	    if (nobjects(nsurveys) .gt. nbmax) 
     *          call quit('Too many data bins!')
	    ok = .true.
	    do 60 i=1, nobjects(nsurveys)
	        call fndcmd(lun, 'bin', args, ok)
	        call rddarg(args, bin_m(i,nsurveys), 1, ok)
	        call rdiarg(args, nsb(i,nsurveys), 1, ok)
	        call rdiarg(args, nb(i,nsurveys), 1, ok)
	        if (.not. ok) call quit('Problem with bin data!')
	        bin_m(i,nsurveys) = bin_m(i,nsurveys) - sys
	        bin_scl(i,nsurveys) = sl_scl(nsb(i,nsurveys), one,
     *                 nb(i,nsurveys), one)
60	    continue
	endif

	n = nobjects(nsurveys)
	call cpar

	return
	end

c-----------------------------------------------------------------------
	function eta (m)

c---	"ETA"
c
c	The detection efficiency at red magnitude m for the
c	current survey.
c
c	This functional form was borrowed from the Baksan
c	neutrino detector group.

c+++  Arguments:
	real*8 eta, m

c+++  Globals:
	include 'tnodata.cb'

c+++  Locals:
	real*8 rt2
	parameter (rt2 = 1.41421356237d0)
	real*8 erfcc
	real*8 one
	parameter (one = 1.)
	real*8 arg, arg2

c---  Default "perfect" survey if ns=0.
	if (ns .eq. 0) then
	    eta = one
	    return
	endif

c---  Pritchett: 
	if (eta_type(ns) .eq. 1) then
	    arg = eta_par(2,ns)*(m-eta_par(1,ns))
	    eta = 0.5 * ( one - arg / sqrt(one + arg**2) )

c---  Heaviside:
	else if (eta_type(ns) .eq. 2) then
	    if (m .gt. eta_par(2,ns)) then
	        eta = 0.
	    else
	        eta = eta_par(1,ns)
	    endif

c---  tanh:  Hyperbolic tangent.
	else if (eta_type(ns) .eq. 3) then
	    arg = (m - eta_par(2,ns)) / eta_par(3,ns)
	    arg = max(arg, -200.)
	    arg = min(arg, 200.)
	    eta = eta_par(1,ns) * (one - tanh(arg))

c---  topc:  Top of Cauchy variance integrand.
	else if (eta_type(ns) .eq. 4) then
	    if (m .ge. eta_par(2,ns)) then
	        eta = 0.
	    else
	        arg = (eta_par(2,ns) - m) / eta_par(3,ns)
	        arg = arg**2
	        eta = eta_par(1,ns) * arg / (one + arg)
	    endif

c---  Linear dropoff:
	else if (eta_type(ns) .eq. 5) then
	    if (m .gt. eta_par(2,ns)) then
	        eta = 0.
	    else if (m .gt. eta_par(1,ns)) then
	        eta = 1. - (m - eta_par(1,ns))/eta_par(3,ns)
	    else
	        eta = 1.
	    endif

c---  tanh_cut:  Hyperbolic tangent with smooth cutoff.
	else if (eta_type(ns) .eq. 6) then
	    arg = (m - eta_par(2,ns)) / eta_par(3,ns)
	    arg = max(arg, -200.)
	    arg = min(arg, 200.)
	    eta = eta_par(1,ns) * (one - tanh(arg))
	    arg = (m - eta_par(4,ns)) / (rt2*eta_par(5,ns))
	    arg = max(arg, -200.)
	    arg = min(arg, 200.)
	    eta = eta * 0.5 * erfcc(arg)

c---  tanh:  Hyperbolic tangent plus one.
	else if (eta_type(ns) .eq. 7) then
	    arg = (eta_par(2,ns) - m) / eta_par(3,ns)
	    arg = max(arg, -200.)
	    arg = min(arg, 200.)
	    eta = eta_par(1,ns) * (one + tanh(arg))

c---  tanh:  Double hyperbolic tangent.
	else if (eta_type(ns) .eq. 8) then
	    arg = (m - eta_par(2,ns)) / eta_par(3,ns)
	    arg = max(arg, -200.)
	    arg = min(arg, 200.)
	    arg2 = (m - eta_par(2,ns)) / eta_par(4,ns)
	    eta = 0.25 * eta_par(1,ns) * (one - tanh(arg))
     *            * (one - tanh(arg2))

c---  erfc2:  Complementary error function, double-width.
	else if (eta_type(ns) .eq. 9) then
	    arg = 0.5 * (m - eta_par(2,ns)) / eta_par(3,ns)
	    arg = max(arg, -200.)
	    arg = min(arg, 200.)
	    eta = 0.5 * eta_par(1,ns) * erfcc(arg)

c---  erfc-rt2:  Complementary error function, rt2-width.
	else if (eta_type(ns) .eq. 10) then
	    arg = (m - eta_par(2,ns)) / (rt2*eta_par(3,ns))
	    arg = max(arg, -200.)
	    arg = min(arg, 200.)
	    eta = 0.5 * eta_par(1,ns) * erfcc(arg)

	endif

	return
	end

c-----------------------------------------------------------------------
	subroutine ljliked (ljl, dljl, nd)

c---	"Log Joint LIKElihood and Derivatives"
c
c	Returns the log joint likelihood & its 1st derivs.
c	(Presently not implemented.)

c+++  Arguments:
	integer nd
	real*8 ljl, dljl(nd)

	call quit('Derivatives not yet implemented!')

	end

c-----------------------------------------------------------------------
	function ljlike ()

c---	"Log Joint LIKElihood"
c
c	Returns the log likelihood for all surveys.

c+++  Arguments:
	real*8 ljlike

c+++  Globals:
	include 'tnodata.cb'
	include 'tnolike.cb'

c+++  Locals:
	real*8 one, log0
	parameter (one = 1., log0 = -1.e10)
	integer i, nt
	real*8 nexp, m1, m2, sll

c+++  Functions:
	real*8 N_exp, N_exp_rng, olike, sllike

c---  First make sure the model is set.
	call setmod

c---  Just loop thru the surveys, adding the log likelihoods
c---  (we're presuming independence here!).
	ljlike = 0.
	do 100 ns=1, nsurveys

c>>>  For list data, use the Poisson point process form.
	    if (s_type(ns) .eq. 1) then
c	        print *, par(1), par(2)
	        ljlike = ljlike - N_exp()
c	        print *, par(1), par(2), N_exp()
c	        write(9,'(a,1pg12.4)') 'Nexp = ', N_exp()
	        do 20 i=1, nobjects(ns)
	            ljlike = ljlike + log(olike(i))
20		continue

c>>>  For number data, use the Poisson counting probability.
	    else if (s_type(ns) .eq. 2) then
	        nexp = N_exp()
	        ljlike = ljlike + nobjects(ns)*log(nexp) - nexp

c>>>  For on/off data, use the Poisson on/off marginal likelihood.
	    else if (s_type(ns) .eq. 3) then
	        do 60 i=1, nobjects(ns)
	            m1 = bin_m(i,ns) - 0.5*dm(ns)
	            m2 = bin_m(i,ns) + 0.5*dm(ns)
	            nexp = N_exp_rng(m1,m2)
	            sll = sllike(nexp, nsb(i,ns), one,
     *                  nb(i,ns), one, bin_scl(i,ns), lcut, nt)
	            if (nt .eq. 0) then
	                ljlike = ljlike + log0
	            else
	                ljlike = ljlike + sll
	            endif
60		continue
	    endif

100	continue

	return
	end

c-----------------------------------------------------------------------
	function olike (i)

c---	"Object LIKElihood"
c
c	Returns the likelihood factor associated with object i
c	in the current survey, presumed to be a "list" survey.

c+++  Arguments:
	real*8 olike
	integer i

c+++  Globals:
	include 'tnodata.cb'
	include 'quad.cb'

c+++  Locals:
	real*8 rtpi, rt2
	parameter (rtpi = 1.772453850905516d0, rt2= 1.414213562373095d0)
	real*8 m
	integer j

c+++  Functions:
	real*8 sigma

c---  For nmq = 1 or zero uncertainty, it's easy!
	if (nmq .eq. 1 .or. (obj_msig(i,ns) .le. 0.)) then
	    olike = sigma(obj_m(i,ns))
	    return
	endif

c---  Otherwise, it's still pretty easy; just sum the quadrature,
c---  noting that the Gaussian is included in the wts.  The rt2
c---  factors come from the def'ns of gauher absc/wts.
	olike = 0.
	do 20 j=1, nmq
	    m = obj_msig(i,ns)*rt2*mq_absc(j) + obj_m(i,ns)
	    olike = olike + mq_wts(j)*sigma(m)
c	    write(9,'(i10,4(1pg12.4))') j, m, mq_wts(j), sigma(m),
c     *            mq_wts(j)*sigma(m)/rtpi
20	continue
	olike = olike / rtpi
c	            write(9,'(i5,1pg12.4)') i, olike

	return
	end

c-----------------------------------------------------------------------
	function N_exp ()

c---	"Number EXPected"
c
c	Returns the expected number of objects for the current survey,
c	found by integrating the efficiency and the model
c	object density, and multiplying by the observing area.

c+++  Arguments:
	real*8 N_exp

c+++  Globals:
	include 'tnodata.cb'
	include 'quad.cb'

c+++  Locals:
	real*8 rt2
	parameter (rt2 = 1.41421356237d0)
	real*8 m1, m2, m3, m4

c+++  Functions:
	real*8 N_exp_rng

c---  Build the integral out of pieces whose size depends on
c---  characteristics of the efficiency function.

c---  Pritchett:
	if (eta_type(ns) .eq. 1) then
	    N_exp = N_exp_rng(m_lo, eta_par(1,ns)) +
     *              N_exp_rng(eta_par(1,ns), m_hi)

c---  Heaviside:
	else if (eta_type(ns) .eq. 2) then
	    N_exp = N_exp_rng(m_lo, eta_par(2,ns))

c---  Tanh:
	else if (eta_type(ns) .eq. 3) then
	    m1 = eta_par(2,ns) - 1.5*eta_par(3,ns)
	    m2 = eta_par(2,ns) + 1.5*eta_par(3,ns)
	    N_exp = N_exp_rng(m_lo, m1) +
     *              N_exp_rng(m1, m2) +
     *              N_exp_rng(m2, m_hi)

c---  Topc:
	else if (eta_type(ns) .eq. 4) then
	    m1 = eta_par(2,ns) - 3.*eta_par(3,ns)
	    N_exp = N_exp_rng(m_lo, m1) +
     *              N_exp_rng(m1, eta_par(2,ns))

c---  Linear:
	else if (eta_type(ns) .eq. 5) then
	    N_exp = N_exp_rng(m_lo, eta_par(1,ns)) +
     *              N_exp_rng(eta_par(1,ns), eta_par(2,ns))


c---  Tanh with cutoff:
c***  Note that we cut the integral off at 6*w above the multiplicative
c     cutoff, to avoid integration problems.
	else if (eta_type(ns) .eq. 6) then
	    if (eta_par(4,ns) .lt. eta_par(2,ns)) then
	        m1 = eta_par(4,ns) - 1.5*eta_par(5,ns)
	        m2 = eta_par(4,ns) + 1.5*eta_par(5,ns)
	        m3 = max(eta_par(2,ns) + 1.5*eta_par(3,ns), m2)
	        m3 = min(m3, m_hi)
	    else
	        m1 = eta_par(2,ns) - 1.5*eta_par(3,ns)
	        m2 = eta_par(2,ns) + 1.5*eta_par(3,ns)
	        m3 = max(eta_par(4,ns) + 1.5*eta_par(5,ns), m2)
	        m3 = min(m3, m_hi)
	    endif
	    m4 = max(eta_par(4,ns) + 6.*eta_par(5,ns), m3)
	    m4 = min(m4, m_hi)
	    N_exp = N_exp_rng(m_lo, m1) +
     *              N_exp_rng(m1, m2) +
     *              N_exp_rng(m2, m3) +
     *              N_exp_rng(m3, m4)

c---  Tanh plus one:
	else if (eta_type(ns) .eq. 7) then
	    m1 = eta_par(2,ns) - 1.5*eta_par(3,ns)
	    m2 = eta_par(2,ns) + 1.5*eta_par(3,ns)
	    N_exp = N_exp_rng(m_lo, m1) +
     *              N_exp_rng(m1, m2) +
     *              N_exp_rng(m2, m_hi)

c---  Double tanh:
	else if (eta_type(ns) .eq. 8) then
	    m1 = eta_par(2,ns) - 1.5*eta_par(3,ns)
	    m2 = eta_par(2,ns) + 1.5*eta_par(3,ns)
	    N_exp = N_exp_rng(m_lo, m1) +
     *              N_exp_rng(m1, m2) +
     *              N_exp_rng(m2, m_hi)

c---  Erfc, 2*width:
	else if (eta_type(ns) .eq. 9) then
c***  Note that we cut the integral off at 6*w above the multiplicative
c     cutoff, to avoid integration problems for bright surveys
c     that use erfc2 (e.g., Larsen).
	    m1 = eta_par(2,ns) - 1.5*eta_par(3,ns)/rt2
	    m2 = eta_par(2,ns) + 1.5*eta_par(3,ns)/rt2
	    m3 = min(eta_par(2,ns) + 6.*eta_par(3,ns), m_hi)
c	    print *, 'lo-1: ', m_lo, m1, N_exp_rng(m_lo, m1)
c	    print *, '1-2:  ', m1, m2, N_exp_rng(m1, m2)
c	    print *, '2-3:  ', m2, m3, N_exp_rng(m2, m3)
	    N_exp = N_exp_rng(m_lo, m1) +
     *              N_exp_rng(m1, m2) +
     *              N_exp_rng(m2, m3)

c---  Erfc, rt2*width:
	else if (eta_type(ns) .eq. 10) then
	    m1 = eta_par(2,ns) - 1.5*eta_par(3,ns)
	    m2 = eta_par(2,ns) + 1.5*eta_par(3,ns)
	    m3 = min(eta_par(2,ns) + 6.*eta_par(3,ns), m_hi)
c	    print *, 'lo-1: ', m_lo, m1, N_exp_rng(m_lo, m1)
c	    print *, '1-2:  ', m1, m2, N_exp_rng(m1, m2)
c	    print *, '2-3:  ', m2, m3, N_exp_rng(m2, m3)
	    N_exp = N_exp_rng(m_lo, m1) +
     *              N_exp_rng(m1, m2) +
     *              N_exp_rng(m2, m3)

	else
	    pause 'No quadrature defined for this eta choice!'
	endif

	return
	end

c-----------------------------------------------------------------------
	function N_exp_rng (m1, m2)

c---	"Number EXPected in m RaNGe"
c
c	Returns the expected number of objects for the current survey,
c	found by integrating the efficiency and the model
c	object density, and multiplying by the observing area.

c+++  Arguments:
	real*8 N_exp_rng, m1, m2

c+++  Globals:
	include 'tnodata.cb'
	include 'quad.cb'

c+++  Locals:
	real*8 es1, es2, y1, y2

c+++  Functions:
	real*8 sigma, eta, N_exp_ignd, N_exp_tignd
	external N_exp_ignd, N_exp_tignd

c---  Since sigma(m) is typically exponential, transform
c---  the integral to facilitate integration of an exponential.
c---  Start by finding the slope (set it to zero if eta*sigma
c---  vanishes at an endpoint).
	es1 = sigma(m1)*eta(m1)
	es2 = sigma(m2)*eta(m2)
	if (es1 .eq. 0. .and. es2 .eq. 0.) then
	    print *, 'N_exp_rng vanishes:'
	    print *, m1, sigma(m1), eta(m1)
	    print *, m2, sigma(m2), eta(m2)
	    N_exp_rng = 0.
	    return
	endif
	if (es1 .le. 0. .or. es2 .le. 0.) then
	    beta = 0.
	else
	    beta = log(es2/es1) / (m2 - m1)
	endif
c	write(9,'(5(1pg12.4))') m1, m2, es1, es2, beta

c---  If the slope is at or near 0 or large just integrate normally.
	if (abs(beta) .le. 0.1 .or. abs(beta) .ge. 20.) then
	    call qromb(N_exp_ignd, m1, m2, N_exp_rng, rq_eps)

c---  Otherwise, transform away the dominant exp'l part.
	else
	    y1 = exp(beta*m1)
	    y2 = exp(beta*m2)
	    call qromb(N_exp_tignd, y1, y2, N_exp_rng, rq_eps)
	    N_exp_rng = N_exp_rng / beta
	endif

c---  Don't forget the solid angle factor!
	N_exp_rng = s_area(ns) * N_exp_rng

c	print *, 'Range: ', m1, m2, beta, N_exp_rng
c	write(9,'(4(1pg12.4), i5)') m1, m2, beta, N_exp_rng
	return
	end

c-----------------------------------------------------------------------
	function N_exp_ignd (m)

c---	"Number EXPected InteGraND"
c
c	Efficiency * object density for the current survey.

c+++  Arguments:
	real*8 N_exp_ignd, m

c+++  Globals:
	include 'quad.cb'

c+++  Functions:
	real*8 sigma, eta

	N_exp_ignd = eta(m) * sigma(m)

	return
	end

c-----------------------------------------------------------------------
	function N_exp_tignd (y)

c---	"Number EXPected Transformed InteGraND"
c
c	Efficiency * object density for the current survey,
c	as a function of an exponential coordinate

c+++  Arguments:
	real*8 N_exp_tignd, y

c+++  Globals:
	include 'quad.cb'

c+++  Locals:
	real*8 m

c+++  Functions:
	real*8 sigma, eta

	m = log(y) / beta
	N_exp_tignd = eta(m) * sigma(m) / y

	return
	end

c=======================================================================
c	Code below defines the model sigma(m) function.
c=======================================================================
c-----------------------------------------------------------------------
	subroutine set_model (n, np)

c---	"SET MODEL"
c
c	Set the model to use for the TNO mag dist'n, and
c	return the # of params.
c
c	n	Model
c	1	Simple exponential (m_0, alpha)

c+++  Arguments:
	integer n, np

c+++  Globals:
	include 'tnolike.cb'

	if (n .lt. 1 .or. n .gt. nmodels) call quit('Invalid model!')

	curmod = n
	np = npars(n)
	curnp = np

	return
	end

c-----------------------------------------------------------------------
	subroutine set_params (params, np)

c---	"SET PARAMeterS"
c
c	Set the parameters of the current model.
c
c	This also converts alpha to base e.

c+++  Arguments:
	integer np
	real*8 params(np)

c+++  Locals:
	real*8 l10
	parameter (l10 = 2.302585092994046d0)
	integer i
	real*8 e

c+++  Globals:
	include 'tnolike.cb'

	do 20 i=1, curnp
	    par(i) = params(i)
20	continue

c---  The l10 factors below adjust base-10 exponents to be
c---  base-e exponents.
	if (curmod .eq. 1) then
	    par(2) = l10*par(2)
	else if (curmod .eq. 3) then
	    e = exp(-par(2)*par(3))
	    par(1) = 10.**par(1) * par(2)*e / (1. - e)
	else if (curmod .eq. 4) then
	    e = exp(-par(2)*par(3))
	    par(1) = par(1) * par(2)*e / (1. - e)
	else if (curmod .eq. 5) then
	    dpl_c = 10.**((par(4)-par(3))*(par(2)-23.))
	    par(1) = (1.+dpl_c)*par(1)
	    par(3) = l10*par(3)
	    par(4) = l10*par(4)
c-- Integral form (incorrect!).
c	    par(1) = l10*(1.+dpl_c)*par(1)
c	    par(3) = l10*par(3)
c	    par(4) = l10*par(4)
	else if (curmod .eq. 6) then
	    par(2) = l10*par(2)
	else if (curmod .eq. 7) then
	    par(2) = l10*par(2)
	    par(3) = l10*par(3)
        endif

	return
	end

c-----------------------------------------------------------------------
	function sigma (m)

c---	"SIGMA"
c
c	Return the object density per unit area per unit m.
c	The area unit is whatever the units of the input
c	survey areas are.

c+++  Arguments:
	real*8 sigma, m

c+++  Globals:
	include 'tnolike.cb'

c+++  Locals:
	real*8 dm

c---  1: Simple exponential (m_0, alpha).
	if (curmod .eq. 1) then
	    sigma = par(2) * exp( par(2) * (m - par(1)) )

c---  2: Constant (A).
	else if (curmod .eq. 2) then
	    sigma = par(1)

c---  3,4: Exponential, scaled below cutoff (logN or N, alpha, m_c).
	else if (curmod .eq. 3 .or. curmod .eq. 4) then
	    sigma = par(1) * exp(par(2)*m)

c---  5: Double power law, as in Bernstein et al. (2004) (Sig23, Req, a1, a2).
	else if (curmod .eq. 5) then
	    dm = m - 23.
	    sigma = par(1) / (exp(-par(3)*dm) + dpl_c*exp(-par(4)*dm))
c-- Incorrect integral form.
c	    e1 = exp(-par(3)*dm)
c	    e2 = exp(-par(4)*dm)
c	    sigma = par(1) * (par(3)*e1 + dpl_c*par(4)*e2) /
c     *                       (e1 + dpl_c*e2)**2

c---  6: Simple exponential (sig23, alpha).
	else if (curmod .eq. 6) then
	    sigma = par(1) * exp( par(2) * (m - 23.) )

c---  7: Rolling power law, as in Bernstein et al. (2004) (Sig23, a1, a2).
	else if (curmod .eq. 7) then
	    dm = m - 23.
	    sigma = par(1) * exp(par(2)*dm + par(3)*dm**2)
	endif

	return
	end

c-----------------------------------------------------------------------
	subroutine sigma_d (m, s, ds, nds)

c---	"SIGMA with Derivatives"
c
c	Return the object density per unit area per unit m.
c	The area unit is whatever the units of the input
c	survey areas are.

c+++  Arguments:
	integer nds
	real*8 m, s, ds(*)

c+++  Globals:
	include 'tnolike.cb'

c---  1: Simple exponential (m_0, alpha)
	if (curmod .eq. 1) then
	    s = par(2) * exp( par(2) * (m - par(1)) )
	    ds(1) = - par(2) * s
	    ds(2) = (m - par(1) - 1./par(2)) * s

c---  2: Constant (A).
	else if (curmod .eq. 2) then
	    s = par(1)
	    ds(1) = 1.

c---  3: Exponential, scaled below cutoff (logN, alpha, m_c).
	else if (curmod .eq. 3) then
	    call quit('Derivs not implemented for this model!')

c---  4: Exponential, scaled below cutoff (N, alpha, m_c).
	else if (curmod .eq. 4) then
	    call quit('Derivs not implemented for this model!')

c---  5: Double power law (Sig23, Req, alpha1, alpha2).
	else if (curmod .eq. 5) then
	    call quit('Derivs not implemented for this model!')

c---  6: Simple exponential (sig23, alpha)
	else if (curmod .eq. 6) then
	    call quit('Derivs not implemented for this model!')

	endif

	return
	end

c-----------------------------------------------------------------------
	subroutine wmodel (lun, name, n, lo, hi, npts)

c---	"Write MODEL"
c
c	Writes the predicted & observed cumulative distribution for
c	survey n to the indicated file.
c
c	The cumulative is calculated as piecewise exponential in m.

c+++  Arguments:
	integer lun, n, npts
	character*(*) name
	real*8 lo, hi

c+++  Globals:
	include 'tnodata.cb'
	include 'tnolike.cb'
	include 'quad.cb'

c+++  Locals:
	real*8 one, log0, tiny
	parameter (one = 1., log0 = -1.e10, tiny = 1.e-30)
	integer maxpts
	parameter (maxpts = 200)
	real*8 l10
	parameter (l10 = 2.302585092994046d0)

	real*8 mvals(maxpts), ndens(maxpts)
	real*8 del, w, Ntot, cum, frxn, alpha, ratio
	real*8 lambda, d, dt, dcum, dcump, rno, pks, mks
	logical logp
	integer i

	real*8 dvals(nomax)
	integer indx(nomax), no

c+++  Functions:
	real*8 eta, sigma, N_exp, N_exp_rng, probks

c---  Check npts.
	if (abs(npts) .gt. maxpts) call quit('Too many pts for wmodel!')

c---  Make sure the survey and model are set, and init step sizes.
	ns = n
	call setmod
	call setstp (lo, hi, npts, logp, del)

c---  Fill the grid.
	do 100 i=1, abs(npts)
	    call stpinf(i, lo, del, logp, mvals(i), w)
	    ndens(i) = s_area(ns) * eta(mvals(i)) * sigma(mvals(i))
100	continue

c---  Open the output file.
	call openf(lun, name)

c---  Calculate cumulative; write results on the fly.
c---  Output:   m  ndens  cum   frxn   slope   eta
	cum = 0.
	Ntot = N_exp()
	write(lun,'(a,1pg12.4,/)') 'Ntot  ', Ntot
	do 200 i=1, abs(npts)
	    if (i .eq. 1) then
	        if (ndens(i)*ndens(i+1) .le. 0.) then
	            alpha = 0.
	        else
	            ratio = ndens(i+1) / ndens(i)
	            alpha = log(ratio) / (mvals(i+1)-mvals(i))
	        endif
	        w = mvals(i+1)-mvals(i)
	    else
	        if (ndens(i)*ndens(i-1) .le. 0.) then
	            alpha = 0.
	        else
	            ratio = ndens(i) / ndens(i-1)
	            alpha = log(ratio) / (mvals(i)-mvals(i-1))
	            if (ratio .gt. 1.e-8) then
	                cum = cum + ndens(i)*
     *                     (exp(alpha*(mvals(i)-mvals(i-1)))-1.)/alpha
	            endif
	        endif
	        w = mvals(i)-mvals(i-1)
	    endif
	    frxn = cum / Ntot
	    alpha = log10(sigma(mvals(i)+.05*w) / sigma(mvals(i)-.05*w))/(.1*w)
	    write(lun,'(6hpred  ,6(1pg12.4))') mvals(i), ndens(i), cum,
     *           frxn, alpha, eta(mvals(i))
200	continue

c---  Write out the data values (at their ML estimates) as a
c---  cumulative histogram, after copying and sorting.
	no = nobjects(ns)
	do 300 i=1, no
	    dvals(i) = obj_m(i,ns)
300	continue

	call indexx (no, dvals, indx)

	write(lun,*)
	cum = 0.
	write(lun,'(5hobs  ,2(1pg12.4))') dvals(indx(1)), cum
	do 320 i=1, no-1
	    cum = cum + 1.
	    write(lun,'(5hobs  ,2(1pg12.4))') dvals(indx(i)), cum
	    write(lun,'(5hobs  ,2(1pg12.4))') dvals(indx(i+1)), cum
320	continue
	cum = cum + 1.
	write(lun,'(5hobs  ,2(1pg12.4))') dvals(indx(no)), cum
	write(lun,'(5hobs  ,2(1pg12.4))') 
     *     2*dvals(indx(no))-dvals(indx(no-1)), cum

c---  Do the KS test.
	dcum = 0.
	dcump = 0.
	d = 0.
	rno = no
	do 400 i=1, no
	    dcum = i / rno
	    cum = N_exp_rng(m_lo,dvals(indx(i))) / Ntot
	    dt = max(abs(dcump-cum),abs(dcum-cum))
	    if (dt .gt. d) then
	        d = dt
	        mks = dvals(indx(i))
	    endif
	    dcump = dcum
400	continue
	lambda = sqrt(rno)*d
	pks = probks(lambda)
	write(lun,'(a,4(1pg12.4))') 'm.d.l.p  ', mks, d, lambda, pks

	close(lun)
	return
	end

c-----------------------------------------------------------------------
      FUNCTION probks(alam)
      REAL*8 probks,alam,EPS1,EPS2
      PARAMETER (EPS1=0.001, EPS2=1.e-8)
      INTEGER j
      REAL*8 a2,fac,term,termbf
      a2=-2.*alam**2
      fac=2.
      probks=0.
      termbf=0.
      do 11 j=1,100
        term=fac*exp(a2*j**2)
        probks=probks+term
        if(abs(term).le.EPS1*termbf.or.abs(term).le.EPS2*probks)return
        fac=-fac
        termbf=abs(term)
11    continue
      probks=1.
      return
      END

c-----------------------------------------------------------------------
	subroutine wbins (lun, name, n, lo, hi, nbins)

c---	"Write MODEL"
c
c	Writes the predicted & observed number of objects
c	in bins.

c+++  Arguments:
	integer lun, n, nbins
	character*(*) name
	real*8 lo, hi

c+++  Globals:
	include 'tnodata.cb'
	include 'tnolike.cb'
	include 'quad.cb'

c+++  Locals:
	integer maxbins
	parameter (maxbins = 50)

	real*8 npred(maxbins)
	real*8 del, w, m, blo, bhi, fac, chi2
	real*8 ocum, etab, blo1, ncum, sig
	logical logb
	integer nobs(maxbins), i, b

c+++  Functions:
	real*8 N_exp_rng, eta

c---  Check npts.
	if (abs(nbins) .gt. maxbins) 
     *      call quit('Too many bins for wbins!')

c---  Make sure the survey and model are set, and set bin info.
	ns = n
	call setmod
	call setbin (lo, hi, nbins, logb, del)
	do 20 i=1, abs(nbins)
	    nobs(i) = 0
20	continue

c---  Bin the objects.
	do 100 i=1, nobjects(ns)
	    call binnum (obj_m(i,ns), lo, del, logb, b)
	    if (b .ge. 1 .and. b .le. abs(nbins)) nobs(b) = nobs(b)+1
100	continue

c---  Calculate predictions and write output.
	call openf(lun, name)
	chi2 = 0.
	do 200 i=1, abs(nbins)
	    call bininf(i, lo, del, logb, m, w)
	    blo = m - 0.5*w
	    bhi = m + 0.5*w
	    npred(i) = N_exp_rng(blo, bhi)
	    chi2 = chi2 + (nobs(i)-npred(i))**2/npred(i)
	    write(lun,210) m, w/2., npred(i), sqrt(npred(i)), nobs(i)
200	continue
210	format('bin ',4(1pg12.4), i4)
	write(lun,'(a,1pg12.4,i5)') 'chisq  ',chi2, abs(nbins)

c---  Write out a rough predicted underlying diff'l dist'n.
	write(lun,*)
	do 300 i=1, abs(nbins)
	    call bininf(i, lo, del, logb, m, w)
	    blo = m - 0.5*w
	    bhi = m + 0.5*w
	    etab = (eta(blo) + 4.*eta(m) + eta(bhi)) / 6.
	    fac = 1./(etab * s_area(ns) * w)
	    write(lun,310) m, w/2., fac*npred(i), fac*sqrt(npred(i)), 
     *         fac*nobs(i)
300	continue
310	format('diff ',5(1pg12.4))

c---  Write out a rough predicted underlying cum dist'n.
c---  NOTE:  We use the ideal survey for ncum.
	write(lun,*)
	ocum = 0.
	do 400 i=1, abs(nbins)
	    call bininf(i, lo, del, logb, m, w)
	    blo = m - 0.5*w
	    if (i .eq. 1) blo1 = blo
	    bhi = m + 0.5*w
	    etab = (eta(blo) + 4.*eta(m) + eta(bhi)) / 6.
	    fac = 1./(etab * s_area(ns))
	    ocum = ocum + fac*nobs(i)
	    npred(i) = N_exp_rng(m_lo, bhi)
	    ns = 0
	    ncum = N_exp_rng(m_lo, bhi)
	    ns = n
	    sig = sqrt(npred(i)) * ncum / npred(i)
	    write(lun,410) bhi, ncum, sig, ocum
400	continue
410	format('cum  ',6(1pg12.4))

	close(lun)
	return
	end
