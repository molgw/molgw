module m_cg

use m_vars
implicit none

private::vdflt,parck,itsum,v2norm,rmdcon,lvmul,livmul,litvmu,ltvmul,lupdat,dbdog,dotprd,stopx,vaxpy,reldst,assst,wzbfgs
private::vcopy,vscopy,vvmulp

public:: deflt,sumit

contains

subroutine deflt ( alg, iv, liv, lv, v )
!
!! DEFLT: supply default values to IV and V.
!
!  Discussion:
!
!   ALG = 1 means regression constants.
!   ALG = 2 means general unconstrained optimization constants.
!
  integer liv
  integer lv

  integer alg
  integer iv(liv)
  real(dp) v(lv)
  !external vdflt
  integer miv, mv
  integer miniv(2), minv(2)
  integer algsav, covprt, covreq, dtype, hc, ierr, inith, inits
  integer ipivot, ivneed, lastiv, lastv, lmat, mxfcal, mxiter
  integer nfcov, ngcov, nvdflt, outlev, parprt, parsav, perm
  integer prunit, qrtyp, rdreq, rmat, solprt, statpr, vneed
  integer vsave, x0prt

  parameter (algsav=51, covprt=14, covreq=15, dtype=16, hc=71 )
  parameter (ierr=75, inith=25, inits=25, ipivot=76, ivneed=3 )
  parameter (lastiv=44, lastv=45, lmat=42, mxfcal=17, mxiter=18 )
  parameter (nfcov=52, ngcov=53, nvdflt=50, outlev=19, parprt=20 )
  parameter (parsav=49, perm=58, prunit=21, qrtyp=80, rdreq=57 )
  parameter (rmat=78, solprt=22, statpr=23, vneed=4, vsave=60 )
  parameter (x0prt=24)

  data miniv(1)/80/, miniv(2)/59/, minv(1)/98/, minv(2)/71/

  if ( alg < 1 .or. 2 < alg ) then
    iv(1) = 67
    return
  end if

  miv = miniv(alg)

  if ( liv < miv ) then
    iv(1) = 15
    return
  end if

  mv = minv(alg)

  if ( lv < mv ) then
    iv(1) = 16
    return
  end if

  call vdflt(alg, lv, v)
  iv(1) = 12
  iv(algsav) = alg
  iv(ivneed) = 0
  iv(lastiv) = miv
  iv(lastv) = mv
  iv(lmat) = mv + 1
  iv(mxfcal) = 1000
  iv(mxiter) = 150
  iv(outlev) = 1
  iv(parprt) = 1
  iv(perm) = miv + 1
! change  
  iv(prunit) = 2
! change  
  iv(solprt) = 1
  iv(statpr) = 1
  iv(vneed) = 0
  iv(x0prt) = 1
!
!  General optimization values.
!
  if ( 2 <= alg ) then

    iv(dtype) = 0
    iv(inith) = 1
    iv(nfcov) = 0
    iv(ngcov) = 0
    iv(nvdflt) = 25
    iv(parsav) = 47
!
!  Regression values.
!
  else

    iv(covprt) = 3
    iv(covreq) = 1
    iv(dtype) = 1
    iv(hc) = 0
    iv(ierr) = 0
    iv(inits) = 0
    iv(ipivot) = 0
    iv(nvdflt) = 32
    iv(parsav) = 67
    iv(qrtyp) = 1
    iv(rdreq) = 3
    iv(rmat) = 0
    iv(vsave) = 58

  end if

  return
end subroutine deflt

subroutine vdflt ( alg, lv, v )
!
!! VDFLT supplies default values to V.
!
!  Discussion:
!
!    alg = 1 means regression constants.
!    alg = 2 means general unconstrained optimization constants.
!
  implicit none

  integer alg, lv
  real(dp) v(lv)
!  real(dp) rmdcon
  real(dp) machep, mepcrt, one, sqteps, three
  integer afctol, bias, cosmin, decfac, delta0, dfac, dinit, dltfdc
  integer dltfdj, dtinit, d0init, epslon, eta0, fuzz, huberc
  integer incfac, lmax0, lmaxs, phmnfc, phmxfc, rdfcmn, rdfcmx
  integer rfctol, rlimit, rsptol, sctol, sigmin, tuner1, tuner2
  integer tuner3, tuner4, tuner5, xctol, xftol

  parameter (one=1.d+0, three=3.d+0)

  parameter (afctol=31, bias=43, cosmin=47, decfac=22, delta0=44 )
  parameter ( dfac=41, dinit=38, dltfdc=42, dltfdj=43, dtinit=39 )
  parameter ( d0init=40, epslon=19, eta0=42, fuzz=45, huberc=48 )
  parameter ( incfac=23, lmax0=35, lmaxs=36, phmnfc=20, phmxfc=21 )
  parameter ( rdfcmn=24, rdfcmx=25, rfctol=32, rlimit=46, rsptol=49 )
  parameter ( sctol=37, sigmin=50, tuner1=26, tuner2=27, tuner3=28 )
  parameter ( tuner4=29, tuner5=30, xctol=33, xftol=34)

  machep = rmdcon(3)
  v(afctol) = 1.d-20

  if ( machep > 1.d-10 ) then
    v(afctol) = machep**2
  end if

  v(decfac) = 0.5d+0
  sqteps = rmdcon(4)
  v(dfac) = 0.6d+0
  v(delta0) = sqteps
  v(dtinit) = 1.d-6
  mepcrt = machep ** (one/three)
  v(d0init) = 1.d+0
  v(epslon) = 0.1d+0
  v(incfac) = 2.d+0
  v(lmax0) = 1.d+0
  v(lmaxs) = 1.d+0
  v(phmnfc) = -0.1d+0
  v(phmxfc) = 0.1d+0
  v(rdfcmn) = 0.1d+0
  v(rdfcmx) = 4.d+0
  v(rfctol) = max (1.d-10, mepcrt**2)
  v(sctol) = v(rfctol)
  v(tuner1) = 0.1d+0
  v(tuner2) = 1.d-4
  v(tuner3) = 0.75d+0
  v(tuner4) = 0.5d+0
  v(tuner5) = 0.75d+0
  v(xctol) = sqteps
  v(xftol) = 1.d+2 * machep

  if ( alg < 2 ) then
    v(cosmin) = max (1.d-6, 1.d+2 * machep)
    v(dinit) = 0.d+0
    v(dltfdc) = mepcrt
    v(dltfdj) = sqteps
    v(fuzz) = 1.5d+0
    v(huberc) = 0.7d+0
    v(rlimit) = rmdcon(5)
    v(rsptol) = 1.d-3
    v(sigmin) = 1.d-4
  else
    v(bias) = 0.8d+0
    v(dinit) = -1.0d+0
    v(eta0) = 1.0d+3 * machep
  end if

  return
end subroutine vdflt

subroutine sumit ( d, fx, g, iv, liv, lv, n, v, x)
!
!  SUMIT carries out unconstrained minimization iterations for SUMSL.
!
!  Discussion:
!
!    The routine uses double-dogleg/BFGS steps.
!
!    parameters iv, n, v, and x are the same as the corresponding
!    ones to sumsl (which see), except that v can be shorter (since
!    the part of v that sumsl uses for storing g is not needed).
!    moreover, compared with sumsl, iv(1) may have the two additional
!    output values 1 and 2, which are explained below, as is the use
!    of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
!    output value from sumsl (and smsno), is not referenced by
!    sumit or the subroutines it calls.
!
!    fx and g need not have been initialized when sumit is called
!    with iv(1) = 12, 13, or 14.
!
! iv(1) = 1 means the caller should set fx to f(x), the function value
!             at x, and call sumit again, having changed none of the
!             other parameters.  an exception occurs if f(x) cannot be
!             (e.g. if overflow would occur), which may happen because
!             of an oversized step.  in this case the caller should set
!             iv(toobig) = iv(2) to 1, which will cause sumit to ig-
!             nore fx and try a smaller step.  the parameter nf that
!             sumsl passes to calcf (for possible use by calcg) is a
!             copy of iv(nfcall) = iv(6).
! iv(1) = 2 means the caller should set g to g(x), the gradient vector
!             of f at x, and call sumit again, having changed none of
!             the other parameters except possibly the scale vector d
!             when iv(dtype) = 0.  the parameter nf that sumsl passes
!             to calcg is iv(nfgcal) = iv(7).  if g(x) cannot be
!             evaluated, then the caller may set iv(nfgcal) to 0, in
!             which case sumit will return with iv(1) = 65.
!
!  Parameters:
!
! d.... scale vector.
! fx... function value.
! g.... gradient vector.
! iv... integer value array.
! liv.. length of iv (at least 60).
! lv... length of v (at least 71 + n*(n+13)/2).
! n.... number of variables (components in x and g).
! v.... floating-point value array.
! x.... vector of parameters to be optimized.
!
  integer liv
  integer lv
  integer n

  integer iv(liv)
  real(dp) d(n)
  real(dp) fx
  real(dp) g(n)
  real(dp) v(lv)
  real(dp) x(n)
  integer dg1, g01, i, k, l, lstgst, nwtst1, step1
  integer        temp1, w, x01, z
  real(dp) t
  real(dp) half, negone, one, onep2, zero
  !logical stopx
  !real(dp) dotprd, reldst, v2norm
  integer cnvcod, dg, dgnorm, dinit, dstnrm, dst0, f, f0, fdif
  integer gthg, gtstep, g0, incfac, inith, irc, kagqt, lmat, lmax0
  integer lmaxs, mode, model, mxfcal, mxiter, nextv, nfcall, nfgcal
  integer ngcall, niter, nreduc, nwtstp, preduc, radfac, radinc
  integer radius, rad0, reldx, restor, step, stglim, stlstg, toobig
  integer tuner4, tuner5, vneed, xirc, x0

  parameter (cnvcod=55, dg=37, g0=48, inith=25, irc=29, kagqt=33 )
  parameter ( mode=35, model=5, mxfcal=17, mxiter=18, nfcall=6 )
  parameter ( nfgcal=7, ngcall=30, niter=31, nwtstp=34, radinc=8 )
  parameter ( restor=9, step=40, stglim=11, stlstg=41, toobig=2 )
  parameter ( vneed=4, xirc=13, x0=43)

  parameter (dgnorm=1, dinit=38, dstnrm=2, dst0=3, f=10, f0=13 )
  parameter ( fdif=11, gthg=44, gtstep=4, incfac=23, lmat=42 )
  parameter ( lmax0=35, lmaxs=36, nextv=47, nreduc=6, preduc=7 )
  parameter ( radfac=16, radius=8, rad0=9, reldx=17, tuner4=29 )
  parameter ( tuner5=30)

  parameter (half=0.5d+0, negone=-1.d+0, one=1.d+0, onep2=1.2d+0, zero=0.d+0)
!
  i = iv(1)
  if (i == 1) go to 50
  if (i == 2) go to 60
!
!   check validity of iv and v input values
!
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  if (iv(1) == 12 .or. iv(1) == 13) then
    iv(vneed) = iv(vneed) + n*(n+13)/2
  end if
  call parck(2, d, iv, liv, lv, n, v)
  i = iv(1) - 2
  if (i > 12) then
    return
  end if
  go to (180, 180, 180, 180, 180, 180, 120, 90, 120, 10, 10, 20), i
!
!   storage allocation
!
10    l = iv(lmat)
  iv(x0) = l + n*(n+1)/2
  iv(step) = iv(x0) + n
  iv(stlstg) = iv(step) + n
  iv(g0) = iv(stlstg) + n
  iv(nwtstp) = iv(g0) + n
  iv(dg) = iv(nwtstp) + n
  iv(nextv) = iv(dg) + n
  if (iv(1) /= 13) go to 20
     iv(1) = 14
     return
!
!   initialization
!
 20   iv(niter) = 0
  iv(nfcall) = 1
  iv(ngcall) = 1
  iv(nfgcal) = 1
  iv(mode) = -1
  iv(model) = 1
  iv(stglim) = 1
  iv(toobig) = 0
  iv(cnvcod) = 0
  iv(radinc) = 0
  v(rad0) = 0.0D+00
  if (v(dinit) >= 0.0D+00) call vscopy(n, d, v(dinit))
  if (iv(inith) /= 1) go to 40
!
!  set the initial hessian approximation to diag(d)**-2
!
     l = iv(lmat)
     call vscopy(n*(n+1)/2, v(l), zero)
     k = l - 1

     do i = 1, n
       k = k + i
       t = d(i)
       if (t <= 0.0D+00) t = one
       v(k) = t
     end do
!
!  compute initial function value
!
 40   iv(1) = 1
  return

 50   v(f) = fx
  if (iv(mode) >= 0) go to 180
  iv(1) = 2
  if (iv(toobig) == 0) then
    return
  end if
     iv(1) = 63
     go to 300
!
!   make sure gradient could be computed
!
 60   if (iv(nfgcal) /= 0) go to 70
     iv(1) = 65
     go to 300

 70   dg1 = iv(dg)
  call vvmulp(n, v(dg1), g, d, -1)
  v(dgnorm) = v2norm(n, v(dg1))

  if (iv(cnvcod) /= 0) go to 290
  if (iv(mode) == 0) go to 250
!
!   allow first step to have scaled 2-norm at most v(lmax0)
!
  v(radius) = v(lmax0)

  iv(mode) = 0
!
!  main loop
!
!   print iteration summary, check iteration limit
!
 80   call itsum(d, g, iv, liv, lv, n, v, x)
 90   k = iv(niter)
  if (k < iv(mxiter)) go to 100
     iv(1) = 10
     go to 300
!
!   update radius
!
 100  iv(niter) = k + 1
  if(k>0)v(radius) = v(radfac) * v(dstnrm)
!
!   initialize for start of next iteration
!
  g01 = iv(g0)
  x01 = iv(x0)
  v(f0) = v(f)
  iv(irc) = 4
  iv(kagqt) = -1
!
!      copy x to x0, g to g0
!
  call vcopy(n, v(x01), x)
  call vcopy(n, v(g01), g)
!
!  Check STOPX and function evaluation limit
!
 110  if ( .not. stopx ( ) ) go to 130
     iv(1) = 11
     go to 140
!
!  Come here when restarting after func. eval. limit or STOPX.
!
 120  if (v(f) >= v(f0)) go to 130
     v(radfac) = one
     k = iv(niter)
     go to 100

 130  if (iv(nfcall) < iv(mxfcal)) go to 150
     iv(1) = 9
 140     if (v(f) >= v(f0)) go to 300
!
!  in case of STOPX or function evaluation limit with
!  improved v(f), evaluate the gradient at x.
!
          iv(cnvcod) = iv(1)
          go to 240
!
!  Compute candidate step
!
 150  step1 = iv(step)
  dg1 = iv(dg)
  nwtst1 = iv(nwtstp)
  if (iv(kagqt) >= 0) go to 160
     l = iv(lmat)
     call livmul(n, v(nwtst1), v(l), g)
     v(nreduc) = half * dotprd(n, v(nwtst1), v(nwtst1))
     call litvmu(n, v(nwtst1), v(l), v(nwtst1))
     call vvmulp(n, v(step1), v(nwtst1), d, 1)
     v(dst0) = v2norm(n, v(step1))
     call vvmulp(n, v(dg1), v(dg1), d, -1)
     call ltvmul(n, v(step1), v(l), v(dg1))
     v(gthg) = v2norm(n, v(step1))
     iv(kagqt) = 0
 160  call dbdog(v(dg1), lv, n, v(nwtst1), v(step1), v)
  if (iv(irc) == 6) go to 180
!
!   check whether evaluating f(x0 + step) looks worthwhile
!
  if (v(dstnrm) <= 0.0D+00) go to 180
  if (iv(irc) /= 5) go to 170
  if (v(radfac) <= one) go to 170
  if (v(preduc) <= onep2 * v(fdif)) go to 180
!
!  Compute f(x0 + step)
!
 170  x01 = iv(x0)
  step1 = iv(step)
  call vaxpy(n, x, one, v(step1), v(x01))
  iv(nfcall) = iv(nfcall) + 1
  iv(1) = 1
  iv(toobig) = 0
  return
!
!  Assess candidate step.
!
 180  x01 = iv(x0)
  v(reldx) = reldst(n, d, x, v(x01))
  call assst(iv, liv, lv, v)
  step1 = iv(step)
  lstgst = iv(stlstg)
  if (iv(restor) == 1) call vcopy(n, x, v(x01))
  if (iv(restor) == 2) call vcopy(n, v(lstgst), v(step1))
  if (iv(restor) /= 3) go to 190
     call vcopy(n, v(step1), v(lstgst))
     call vaxpy(n, x, one, v(step1), v(x01))
     v(reldx) = reldst(n, d, x, v(x01))

 190  k = iv(irc)
  go to (200,230,230,230,200,210,220,220,220,220,220,220,280,250), k
!
!      recompute step with changed radius
!
 200     v(radius) = v(radfac) * v(dstnrm)
     go to 110
!
!   compute step of length v(lmaxs) for singular convergence test.
!
 210  v(radius) = v(lmaxs)
  go to 150
!
!   convergence or false convergence
!
 220  iv(cnvcod) = k - 4
  if (v(f) >= v(f0)) go to 290
     if (iv(xirc) == 14) go to 290
          iv(xirc) = 14
!
!  Process acceptable step.
!
 230  if (iv(irc) /= 3) go to 240
     step1 = iv(step)
     temp1 = iv(stlstg)
!
!      set  temp1 = hessian * step  for use in gradient tests
!
     l = iv(lmat)
     call ltvmul(n, v(temp1), v(l), v(step1))
     call lvmul(n, v(temp1), v(l), v(temp1))
!
!   compute gradient
!
 240  iv(ngcall) = iv(ngcall) + 1
  iv(1) = 2
  return
!
!   initializations -- g0 = g - g0, etc.
!
 250  g01 = iv(g0)
  call vaxpy(n, v(g01), negone, v(g01), g)
  step1 = iv(step)
  temp1 = iv(stlstg)
  if (iv(irc) /= 3) go to 270
!
!   set v(radfac) by gradient tests
!
!  Set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))
!
     call vaxpy(n, v(temp1), negone, v(g01), v(temp1))
     call vvmulp(n, v(temp1), v(temp1), d, -1)
!
!  Do gradient tests
!
     if (v2norm(n, v(temp1)) <= v(dgnorm) * v(tuner4)) then
       go to 260
     end if

     if (dotprd(n, g, v(step1)) >= v(gtstep) * v(tuner5))  then
       go to 270
     end if

 260               v(radfac) = v(incfac)
!
!   update h, loop
!
 270  w = iv(nwtstp)
  z = iv(x0)
  l = iv(lmat)
  call wzbfgs(v(l), n, v(step1), v(w), v(g01), v(z))
!
!  Use the n-vectors starting at v(step1) and v(g01) for scratch.
!
  call lupdat(v(temp1), v(step1), v(l), v(g01), v(l), n, v(w), v(z))
  iv(1) = 2
  go to 80
!
!   misc. details
!
!   bad parameters to assess
!
 280  iv(1) = 64
  go to 300
!
!  Print summary of final iteration and other requested items
!
 290  iv(1) = iv(cnvcod)
  iv(cnvcod) = 0
 300  call itsum(d, g, iv, liv, lv, n, v, x)

  return
end subroutine sumit 

subroutine parck ( alg, d, iv, liv, lv, n, v )
!
!  PARCK checks parameters, prints changed values.
!
!  Discussion:
!
!    alg = 1 for regression, alg = 2 for general unconstrained opt.
!
  integer alg, liv, lv, n
  integer iv(liv)
  real(dp) d(n), v(lv)
!  real(dp) rmdcon
  integer max0
  integer i, ii, iv1, j, k, l, m, miv1, miv2, ndfalt, parsv1, pu
  integer ijmp, jlim(2), miniv(2), ndflt(2)
  character*1 varnm(2), sh(2)
  character*4 cngd(3), dflt(3), vn(2,34), which(3)
  real(dp) big, machep, tiny, vk, vm(34), vx(34)
  integer algsav, dinit, dtype, dtype0, epslon, inits, ivneed
  integer lastiv, lastv, lmat, nextiv, nextv, nvdflt, oldn
  integer parprt, parsav, perm, prunit, vneed

  parameter (algsav=51, dinit=38, dtype=16, dtype0=54, epslon=19 )
  parameter ( inits=25, ivneed=3, lastiv=44, lastv=45, lmat=42 )
  parameter ( nextiv=46, nextv=47, nvdflt=50, oldn=38, parprt=20 )
  parameter ( parsav=49, perm=58, prunit=21, vneed=4)
  save big, machep, tiny

  data big/0.d+0/, machep/-1.d+0/, tiny/1.d+0/

     data vn(1,1),vn(2,1)/'epsl','on..'/
     data vn(1,2),vn(2,2)/'phmn','fc..'/
     data vn(1,3),vn(2,3)/'phmx','fc..'/
     data vn(1,4),vn(2,4)/'decf','ac..'/
     data vn(1,5),vn(2,5)/'incf','ac..'/
     data vn(1,6),vn(2,6)/'rdfc','mn..'/
     data vn(1,7),vn(2,7)/'rdfc','mx..'/
     data vn(1,8),vn(2,8)/'tune','r1..'/
     data vn(1,9),vn(2,9)/'tune','r2..'/
     data vn(1,10),vn(2,10)/'tune','r3..'/
     data vn(1,11),vn(2,11)/'tune','r4..'/
     data vn(1,12),vn(2,12)/'tune','r5..'/
     data vn(1,13),vn(2,13)/'afct','ol..'/
     data vn(1,14),vn(2,14)/'rfct','ol..'/
     data vn(1,15),vn(2,15)/'xcto','l...'/
     data vn(1,16),vn(2,16)/'xfto','l...'/
     data vn(1,17),vn(2,17)/'lmax','0...'/
     data vn(1,18),vn(2,18)/'lmax','s...'/
     data vn(1,19),vn(2,19)/'scto','l...'/
     data vn(1,20),vn(2,20)/'dini','t...'/
     data vn(1,21),vn(2,21)/'dtin','it..'/
     data vn(1,22),vn(2,22)/'d0in','it..'/
     data vn(1,23),vn(2,23)/'dfac','....'/
     data vn(1,24),vn(2,24)/'dltf','dc..'/
     data vn(1,25),vn(2,25)/'dltf','dj..'/
     data vn(1,26),vn(2,26)/'delt','a0..'/
     data vn(1,27),vn(2,27)/'fuzz','....'/
     data vn(1,28),vn(2,28)/'rlim','it..'/
     data vn(1,29),vn(2,29)/'cosm','in..'/
     data vn(1,30),vn(2,30)/'hube','rc..'/
     data vn(1,31),vn(2,31)/'rspt','ol..'/
     data vn(1,32),vn(2,32)/'sigm','in..'/
     data vn(1,33),vn(2,33)/'eta0','....'/
     data vn(1,34),vn(2,34)/'bias','....'/

  data vm(1)/1.0d-3/, vm(2)/-0.99d+0/, vm(3)/1.0d-3/, vm(4)/1.0d-2/
  data vm(5)/1.2d+0/, vm(6)/1.d-2/, vm(7)/1.2d+0/, vm(8)/0.d+0/
  data vm(9)/0.d+0/, vm(10)/1.d-3/, vm(11)/-1.d+0/, vm(13)/0.d+0/
  data vm(15)/0.d+0/, vm(16)/0.d+0/, vm(19)/0.d+0/, vm(20)/-10.d+0/
  data vm(21)/0.d+0/, vm(22)/0.d+0/, vm(23)/0.d+0/, vm(27)/1.01d+0/
  data vm(28)/1.d+10/, vm(30)/0.d+0/, vm(31)/0.d+0/, vm(32)/0.d+0/
  data vm(34)/0.d+0/

  data vx(1)/0.9d+0/, vx(2)/-1.d-3/, vx(3)/1.d+1/, vx(4)/0.8d+0/
  data vx(5)/1.d+2/, vx(6)/0.8d+0/, vx(7)/1.d+2/, vx(8)/0.5d+0/
  data vx(9)/0.5d+0/, vx(10)/1.d+0/, vx(11)/1.d+0/, vx(14)/0.1d+0/
  data vx(15)/1.d+0/, vx(16)/1.d+0/, vx(19)/1.d+0/, vx(23)/1.d+0/
  data vx(24)/1.d+0/, vx(25)/1.d+0/, vx(26)/1.d+0/, vx(27)/1.d+10/
  data vx(29)/1.d+0/, vx(31)/1.d+0/, vx(32)/1.d+0/, vx(33)/1.d+0/
  data vx(34)/1.d+0/

  data varnm(1)/'p'/, varnm(2)/'n'/, sh(1)/'s'/, sh(2)/'h'/
  data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/
  data dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/
  data ijmp/33/, jlim(1)/0/, jlim(2)/24/, ndflt(1)/32/, ndflt(2)/25/
  data miniv(1)/80/, miniv(2)/59/

  pu = 0
  if (prunit <= liv) pu = iv(prunit)
  if (alg < 1 .or. alg > 2) go to 340
  if (iv(1) == 0) call deflt(alg, iv, liv, lv, v)
  iv1 = iv(1)
  if (iv1 /= 13 .and. iv1 /= 12) go to 10
  miv1 = miniv(alg)
  if (perm <= liv) miv1 = max0(miv1, iv(perm) - 1)
  if (ivneed <= liv) miv2 = miv1 + max0(iv(ivneed), 0)
  if (lastiv <= liv) iv(lastiv) = miv2
  if (liv < miv1) go to 300
  iv(ivneed) = 0
  iv(lastv) = max0(iv(vneed), 0) + iv(lmat) - 1
  iv(vneed) = 0
  if (liv < miv2) go to 300
  if (lv < iv(lastv)) go to 320
 10   if (alg == iv(algsav)) go to 30
  if (pu /= 0) write(pu,20) alg, iv(algsav)
 20 format(/39h the first parameter to deflt should be,i3,12h rather than,i3)
     iv(1) = 82
     return
 30   if (iv1 < 12 .or. iv1 > 14) go to 60
     if (n >= 1) go to 50
          iv(1) = 81
          if (pu == 0) then
            return
          end if
          write(pu,40) varnm(alg), n
 40           format(/8h /// bad,a1,2h =,i5)
          return
 50      if (iv1 /= 14) iv(nextiv) = iv(perm)
     if (iv1 /= 14) iv(nextv) = iv(lmat)
     if (iv1 == 13) then
       return
     end if
     k = iv(parsav) - epslon
     call vdflt(alg, lv-k, v(k+1))
     iv(dtype0) = 2 - alg
     iv(oldn) = n
     which(1) = dflt(1)
     which(2) = dflt(2)
     which(3) = dflt(3)
     go to 110
 60   if (n == iv(oldn)) go to 80
     iv(1) = 17
     if (pu == 0) then
       return
     end if
     write(pu,70) varnm(alg), iv(oldn), n
 70      format(/5h /// ,1a1,14h changed from ,i5,4h to ,i5)
     return

 80   if (iv1 <= 11 .and. iv1 >= 1) go to 100
     iv(1) = 80
     if (pu /= 0) write(pu,90) iv1
 90      format(/13h ///  iv(1) =,i5,28h should be between 0 and 14.)
     return

 100  which(1) = cngd(1)
  which(2) = cngd(2)
  which(3) = cngd(3)

 110  if (iv1 == 14) iv1 = 12
  if (big > tiny) go to 120
     tiny = rmdcon(1)
     machep = rmdcon(3)
     big = rmdcon(6)
     vm(12) = machep
     vx(12) = big
     vx(13) = big
     vm(14) = machep
     vm(17) = tiny
     vx(17) = big
     vm(18) = tiny
     vx(18) = big
     vx(20) = big
     vx(21) = big
     vx(22) = big
     vm(24) = machep
     vm(25) = machep
     vm(26) = machep
     vx(28) = rmdcon(5)
     vm(29) = machep
     vx(30) = big
     vm(33) = machep
 120  m = 0
  i = 1
  j = jlim(alg)
  k = epslon
  ndfalt = ndflt(alg)

  do l = 1, ndfalt
    vk = v(k)
    if (vk >= vm(i) .and. vk <= vx(i)) go to 140
      m = k
      if (pu /= 0) write(pu,130) vn(1,i), vn(2,i), k, vk,vm(i), vx(i)
 130  format(/6h ///  ,2a4,5h.. v(,i2,3h) =,d11.3,7h should, &
      11h be between,d11.3,4h and,d11.3)
 140  k = k + 1
     i = i + 1
     if (i == j) i = ijmp
  end do

  if (iv(nvdflt) == ndfalt) go to 170
     iv(1) = 51
     if (pu == 0) then
       return
     end if
     write(pu,160) iv(nvdflt), ndfalt
 160     format(/13h iv(nvdflt) =,i5,13h rather than ,i5)
     return
 170  if ((iv(dtype) > 0 .or. v(dinit) > 0.0D+00) .and. iv1 == 12) then
             go to 200
  end if

  do i = 1, n
     if (d(i) > 0.0D+00) go to 190
          m = 18
          if (pu /= 0) write(pu,180) i, d(i)
 180     format(/8h ///  d(,i3,3h) =,d11.3,19h should be positive)
 190     continue
  end do

 200  if (m == 0) go to 210
     iv(1) = m
     return

 210  if (pu == 0 .or. iv(parprt) == 0) then
        return
      end if
  if (iv1 /= 12 .or. iv(inits) == alg-1) go to 230
     m = 1
     write(pu,220) sh(alg), iv(inits)
 220 format(/22h nondefault values..../5h init,a1,14h      iv(25) =,i3)
 230  if (iv(dtype) == iv(dtype0)) go to 250
     if (m == 0) write(pu,260) which
     m = 1
     write(pu,240) iv(dtype)
 240     format(20h dtype      iv(16) =,i3)
 250  i = 1
  j = jlim(alg)
  k = epslon
  l = iv(parsav)
  ndfalt = ndflt(alg)

  do ii = 1, ndfalt
     if (v(k) == v(l)) go to 280
          if (m == 0) write(pu,260) which
 260          format(/1h ,3a4,9halues..../)
          m = 1
          write(pu,270) vn(1,i), vn(2,i), k, v(k)
 270          format(1x,2a4,5h.. v(,i2,3h) =,d15.7)
 280     k = k + 1
     l = l + 1
     i = i + 1
     if (i == j) i = ijmp
  end do

  iv(dtype0) = iv(dtype)
  parsv1 = iv(parsav)
  call vcopy(iv(nvdflt), v(parsv1), v(epslon))
  return

 300  iv(1) = 15
  if (pu == 0) then
    return
  end if
  write(pu,310) liv, miv2
 310  format(/10h /// liv =,i5,17h must be at least,i5)
  if (liv < miv1) then
    return
  end if
  if (lv < iv(lastv)) go to 320
  return

 320  iv(1) = 16
  if (pu == 0) then
    return
  end if
  write(pu,330) lv, iv(lastv)
 330  format(/9h /// lv =,i5,17h must be at least,i5)
  return

 340  iv(1) = 67
  if (pu == 0) then
    return
  end if
  write(pu,350) alg
 350  format(/10h /// alg =,i5,15h must be 1 or 2)

  return
end subroutine parck

subroutine itsum ( d, g, iv, liv, lv, p, v, x )

!
! ITSUM prints an iteration summary.
!
  integer liv
  integer lv
  integer p

  real(dp) d(p)
  real(dp) g(p)
  integer iv(liv)
  real(dp) v(lv)
  real(dp) x(p)
  integer alg, i, iv1, m, nf, ng, ol, pu
  character*4 model1(6), model2(6)
  real(dp) nreldf, oldf, preldf, reldf
  integer algsav, dstnrm, f, fdif, f0, needhd, nfcall, nfcov, ngcov
  integer ngcall, niter, nreduc, outlev, preduc, prntit, prunit
  integer reldx, solprt, statpr, stppar, sused, x0prt

  parameter (algsav=51, needhd=36, nfcall=6, nfcov=52, ngcall=30 )
  parameter ( ngcov=53, niter=31, outlev=19, prntit=39, prunit=21 )
  parameter ( solprt=22, statpr=23, sused=64, x0prt=24)
  parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6, preduc=7, reldx=17 )
  parameter ( stppar=5)

  data model1/'    ','    ','    ','    ','  g ','  s '/
  data model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/

  pu = iv(prunit)

  if ( pu == 0 ) then
    return
  end if

  iv1 = iv(1)
  if (iv1 > 62) iv1 = iv1 - 51
  ol = iv(outlev)
  alg = iv(algsav)
  if (iv1 < 2 .or. iv1 > 15) go to 370
  if (iv1 >= 12) go to 120
  if (iv1 == 2 .and. iv(niter) == 0) go to 390
  if (ol == 0) go to 120
  if (iv1 >= 10 .and. iv(prntit) == 0) go to 120
  if (iv1 > 2) go to 10
     iv(prntit) = iv(prntit) + 1
     if (iv(prntit) < iabs(ol)) then
       return
     end if
 10   nf = iv(nfcall) - iabs(iv(nfcov))
  iv(prntit) = 0
  reldf = 0.0D+00
  preldf = 0.0D+00
  oldf = max (abs(v(f0)), abs(v(f)))
  if (oldf <= 0.0D+00) go to 20
     reldf = v(fdif) / oldf
     preldf = v(preduc) / oldf
 20   if (ol > 0) go to 60
!
!  print short summary line
!
     if (iv(needhd) == 1 .and. alg == 1) write(pu,30)
 30   format(/'   it   nf',6x,'f',7x,'reldf',3x,'preldf',3x,'reldx', &
     &  2x,'model  stppar')
     if (iv(needhd) == 1 .and. alg == 2) write(pu,40)
 40  format(/'    it   nf',7x,'f',8x,'reldf',4x,'preldf',4x,'reldx   stppar')
     iv(needhd) = 0
     if (alg == 2) go to 50
     m = iv(sused)
     write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
       model1(m), model2(m), v(stppar)
     go to 120

 50  write(pu,110) iv(niter), nf, v(f), reldf, preldf, v(reldx), v(stppar)
     go to 120
!
!  print long summary line
!
 60   if (iv(needhd) == 1 .and. alg == 1) write(pu,70)
 70   format(/11h    it   nf,6x,1hf,7x,5hreldf,3x,6hpreldf,3x,5hreldx, &
        2x,13hmodel  stppar,2x,6hd*step,2x,7hnpreldf)
  if (iv(needhd) == 1 .and. alg == 2) write(pu,80)
 80   format(/11h    it   nf,7x,1hf,8x,5hreldf,4x,6hpreldf,4x,5hreldx, &
        3x,6hstppar,3x,6hd*step,3x,7hnpreldf)
  iv(needhd) = 0
  nreldf = 0.0D+00
  if (oldf > 0.0D+00) nreldf = v(nreduc) / oldf
  if (alg == 2) go to 90
  m = iv(sused)
  write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
              model1(m), model2(m), v(stppar), v(dstnrm), nreldf
  go to 120

 90   write(pu,110) iv(niter), nf, v(f), reldf, preldf, &
             v(reldx), v(stppar), v(dstnrm), nreldf
 100  format(i6,i5,d10.3,2d9.2,d8.1,a3,a4,2d8.1,d9.2)
 110  format(i6,i5,d11.3,2d10.2,3d9.1,d10.2)

 120  if (iv(statpr) < 0) go to 430
  go to (999, 999, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, &
    330, 350, 520), iv1

 130  write(pu,140)
 140  format(/' x-convergence')
  go to 430

 150  write ( pu, 160 )
 160  format(/'relative function convergence')
  go to 430

 170  write(pu,180)
 180  format(/'x- and relative function convergence')
  go to 430

 190  write(pu,200)
 200  format(/'Absolute function convergence.')
  go to 430

 210  write(pu,220)
 220  format(/'Singular convergence.')
  go to 430

 230  write(pu,240)
 240  format(/'False convergence.')
  go to 430

 250  write(pu,260)
 260  format(/'Function evaluation limit.')
  go to 430

 270  write(pu,280)
 280  format(/'Iteration limit.')
  go to 430

 290  write(pu,300)
 300  format(/'STOPX')
  go to 430

 310  write(pu,320)
 320  format(/'Initial f(x) cannot be computed.')

  go to 390

 330  write(pu,340)
 340  format(/'Bad parameters to assess.')
  go to 999

 350  write(pu,360)
 360  format(/'Gradient could not be computed.')
  if (iv(niter) > 0) go to 480
  go to 390

 370  write(pu,380) iv(1)
 380  format(/'iv(1) =',i5)
  go to 999
!
!   initial call on itsum
!
 390  if (iv(x0prt) /= 0) write(pu,400) (i, x(i), d(i), i = 1, p)
 400  format(/23h     i     initial x(i),8x,4hd(i)//(1x,i5,d17.6,d14.3))
!     the following are to avoid undefined variables when the
!     function evaluation limit is 1...
!
  v(dstnrm) = 0.0D+00
  v(fdif) = 0.0D+00
  v(nreduc) = 0.0D+00
  v(preduc) = 0.0D+00
  v(reldx)  = 0.0D+00
  if (iv1 >= 12) go to 999
  iv(needhd) = 0
  iv(prntit) = 0
  if (ol == 0) go to 999
  if (ol < 0 .and. alg == 1) write(pu,30)
  if (ol < 0 .and. alg == 2) write(pu,40)
  if (ol > 0 .and. alg == 1) write(pu,70)
  if (ol > 0 .and. alg == 2) write(pu,80)
  if (alg == 1) write(pu,410) v(f)
  if (alg == 2) write(pu,420) v(f)
 410  format(/11h     0    1,d10.3)
!365  format(/11h     0    1,e11.3)
 420  format(/11h     0    1,d11.3)
  go to 999
!
!  Print various information requested on solution
!
 430  iv(needhd) = 1
  if (iv(statpr) == 0) go to 480
     oldf = max (abs(v(f0)), abs(v(f)))
     preldf = 0.0D+00
     nreldf = 0.0D+00
     if (oldf <= 0.0D+00) go to 440
          preldf = v(preduc) / oldf
          nreldf = v(nreduc) / oldf
 440     nf = iv(nfcall) - iv(nfcov)
     ng = iv(ngcall) - iv(ngcov)
     write(pu,450) v(f), v(reldx), nf, ng, preldf, nreldf
 450  format(/9h function,d17.6,8h   reldx,d17.3/12h func. evals, &
    i8,9x,11hgrad. evals,i8/7h preldf,d16.3,6x,7hnpreldf,d15.3)

     if (iv(nfcov) > 0) write(pu,460) iv(nfcov)
 460     format(/1x,i4,50h extra func. evals for covariance and diagnostics.)
     if (iv(ngcov) > 0) write(pu,470) iv(ngcov)
 470     format(1x,i4,50h extra grad. evals for covariance and diagnostics.)

 480  if (iv(solprt) == 0) go to 999
     iv(needhd) = 1
     write(pu,490)
 490  format(/22h     i      final x(i),8x,4hd(i),10x,4hg(i)/)
     do i = 1, p
          write(pu,510) i, x(i), d(i), g(i)
     end do
 510     format(1x,i5,d16.6,2d14.3)
  go to 999

 520  write(pu,530)
 530  format(/'Inconsistent dimensions.')
 999  continue

  return
end subroutine itsum

function v2norm ( p, x ) result(vv2norm)
!
!! V2NORM returns the 2-norm of the p-vector X.
!
!  Discussion:
!
!    The routine tries to avoid underflow.
!
!  Parameters:
!
  integer p

  real(dp) x(p)
  integer i, j
  real(dp) r, scale
  real(dp), save :: sqteta = 0.0D+00
  real(dp) t, xi
!  real(dp) rmdcon
  real(dp) vv2norm

  vv2norm = 0.0D+00

  if (p <= 0 ) then
    return
  end if

  if ( all ( x(1:p) == 0.0D+00 ) ) then
    return
  end if

  scale = 0.0D+00
  do i = 1, p
    if ( x(i) /= 0.0D+00 ) then
      scale = abs(x(i))
      exit
    end if
  end do

  if ( scale == 0.0D+00 ) then
    return
  end if

  if ( p <= i ) then
    vv2norm = scale
    return
  end if

  t = 1.0D+00
  if ( sqteta == 0.0D+00 ) then
    sqteta = rmdcon(2)
  end if
!
!  sqteta is (slightly larger than) the square root of the
!  smallest positive floating point number on the machine.
!  the tests involving sqteta are done to prevent underflows.
!
  j = i + 1
  do i = j, p
    xi = abs(x(i))
    if (xi <= scale) then
      r = xi / scale
      if (r > sqteta) t = t + r*r
    else
      r = scale / xi
      if (r <= sqteta) r = 0.0D+00
      t = 1.0D+00  +  t * r*r
      scale = xi
    end if
  end do

  vv2norm = scale * sqrt(t)

end function v2norm

function rmdcon ( k ) result(mcons)
!
!! RMDCON returns machine dependent constants.
!
!  Discussion:
!
!    Comments below contain data statements for various machines.
!    To convert to another machine, place a c in column 1 of the
!    data statement line(s) that correspond to the current machine
!    and remove the c from column 1 of the data statement line(s)
!    that correspond to the new machine.
!
!    the constant returned depends on k...
!
!         k = 1... smallest pos. eta such that -eta exists.
!         k = 2... square root of eta.
!         k = 3... unit roundoff = smallest pos. no. machep such
!                  that 1 + machep > 1 .and. 1 - machep < 1.
!         k = 4... square root of machep.
!         k = 5... square root of big (see k = 6).
!         k = 6... largest machine no. big such that -big exists.
!
  integer k
  real(dp) mcons
  real(dp) big, eta, machep
  integer bigi(4), etai(4), machei(4)
  equivalence (big,bigi(1)), (eta,etai(1)), (machep,machei(1))
!
!  ibm 360, ibm 370, or xerox
!
!     data big/z7fffffffffffffff/, eta/z0010000000000000/,
!    1     machep/z3410000000000000/
!
!  data general
!
!     data big/0.7237005577d+76/, eta/0.5397605347d-78/,
!    1     machep/2.22044605d-16/
!
!  dec 11
!
!     data big/1.7d+38/, eta/2.938735878d-39/, machep/2.775557562d-17/
!
!  hp3000
!
!     data big/1.157920892d+77/, eta/8.636168556d-78/,
!    1     machep/5.551115124d-17/
!
!  honeywell
!
!     data big/1.69d+38/, eta/5.9d-39/, machep/2.1680435d-19/
!
!  dec10
!
!     data big/"377777100000000000000000/,
!    1     eta/"002400400000000000000000/,
!    2     machep/"104400000000000000000000/
!
!  burroughs
!
!     data big/o0777777777777777,o7777777777777777/,
!    1     eta/o1771000000000000,o7770000000000000/,
!    2     machep/o1451000000000000,o0000000000000000/
!
!  control data
!
!     data big/37767777777777777777b,37167777777777777777b/,
!    1     eta/00014000000000000000b,00000000000000000000b/,
!    2     machep/15614000000000000000b,15010000000000000000b/
!
!  prime
!
!     data big/1.0d+9786/, eta/1.0d-9860/, machep/1.4210855d-14/
!
!  univac
!
!     data big/8.988d+307/, eta/1.2d-308/, machep/1.734723476d-18/
!
!  vax
!
  data big/1.7d+38/, eta/2.939d-39/, machep/1.3877788d-17/
!
!  cray 1
!
!     data bigi(1)/577767777777777777777b/,
!    1     bigi(2)/000007777777777777776b/,
!    2     etai(1)/200004000000000000000b/,
!    3     etai(2)/000000000000000000000b/,
!    4     machei(1)/377224000000000000000b/,
!    5     machei(2)/000000000000000000000b/
!
!  port library -- requires more than just a data statement...
!
!     external d1mach
!     real(dp) d1mach, zero
!     data big/0.d+0/, eta/0.d+0/, machep/0.d+0/, zero/0.d+0/
!     if (big > 0.0D+00) go to 1
!        big = d1mach(2)
!        eta = d1mach(1)
!        machep = d1mach(4)
!1    continue
!
! end of port
!
!  body -
!
  go to (10, 20, 30, 40, 50, 60), k

 10   mcons = eta
  return

 20   mcons = sqrt(256.d+0*eta)/16.d+0
  return

 30   mcons = machep
  return

 40   mcons = sqrt(machep)
  return

 50   mcons = sqrt(big/256.d+0)*16.d+0
  return

 60   mcons = big

end function rmdcon

subroutine vcopy ( p, y, x )
!
!! VCOPY sets y = x.
!
!  Discussion:
!
!    x and y are p-vectors
!
  implicit none

  integer p
  real(dp) x(p)
  real(dp) y(p)

  y(1:p) = x(1:p)

  return
end subroutine vcopy

subroutine vscopy ( p, y, s )
!
!! VSCOPY sets the vector Y to scalar S.
!
  implicit none
  integer p
  real(dp) s
  real(dp) y(p)

  y(1:p) = s

  return
end subroutine vscopy

subroutine vvmulp ( n, x, y, z, k )
!
!! VVMULP sets x(i) = y(i) * z(i)**k, 1 <= i <= n (for k = 1 or -1)
!
  implicit none
  integer n
  integer k
  real(dp) x(n)
  real(dp) y(n)
  real(dp) z(n)

  if ( k < 0 ) then
    x(1:n) = y(1:n) / z(1:n)
  else
    x(1:n) = y(1:n) * z(1:n)
  end if

  return
end subroutine vvmulp

subroutine lvmul ( n, x, l, y )
!
!! LVMUL computes x = L * y.
!
!  Discussion:
!
!    L  is an  n x n  lower triangular matrix stored compactly by rows.
!    x and y may occupy the same storage.
!
  integer n

  real(dp) x(n), l(*), y(n)
!     dimension l(n*(n+1)/2)
  integer i, ii, ij, i0, j, np1
  real(dp) t

  np1 = n + 1
  i0 = n*(n+1)/2

  do ii = 1, n
    i = np1 - ii
    i0 = i0 - i
    t = 0.0D+00
    do j = 1, i
      ij = i0 + j
      t = t + l(ij)*y(j)
    end do
    x(i) = t
  end do

  return
end subroutine lvmul

subroutine livmul ( n, x, l, y )
!
!! LIVMUL solves L * x = y.
!
!  Discussion:
!
!    L is an  n x n  lower triangular
!    matrix stored compactly by rows.  x and y may occupy the same
!    storage.
!
  integer n

  real(dp) x(n), l(*), y(n)
!  external dotprd
!  real(dp) dotprd
  integer i, j, k
  real(dp) t

  do k = 1, n
    if (y(k) /= 0.0D+00 ) go to 20
    x(k) = 0.0D+00
  end do

  return

20 continue

  j = k*(k+1)/2
  x(k) = y(k) / l(j)

  if (k >= n) then
    return
  end if

  k = k + 1

  do i = k, n
     t = dotprd(i-1, l(j+1), x)
     j = j + i
     x(i) = (y(i) - t)/l(j)
  end do

  return
end subroutine livmul

subroutine litvmu ( n, x, l, y )
!
!! LITVMU solves L' * x = y.
!
!  Discussion:
!
!    L is an  n x n  lower triangular
!    matrix stored compactly by rows.  x and y may occupy the same
!    storage.
!
  integer n

  real(dp) l(*)
  real(dp) x(n)
  real(dp) y(n)
  integer i, ii, ij, i0, j
  real(dp) xi

  x(1:n) = y(1:n)

  i0 = n*(n+1)/2

  do ii = 1, n
    i = n+1 - ii
    xi = x(i)/l(i0)
    x(i) = xi
    if ( i <= 1 ) then
      exit
    end if
    i0 = i0 - i
    if ( xi /= 0.0D+00 ) then
      do j = 1, i-1
        ij = i0 + j
        x(j) = x(j) - xi*l(ij)
      end do
    end if
  end do

  return
end subroutine litvmu

subroutine ltvmul ( n, x, l, y )
!
!! LTVMUL computes  x = (l**t)*y.
!
!  Discussion:
!
!    L is an  n x n  lower triangular matrix stored compactly by rows.
!    x and y may occupy the same storage.
!
  integer n
  real(dp) x(n), l(*), y(n)
!     dimension l(n*(n+1)/2)
  integer i, ij, i0, j
  real(dp) yi

  i0 = 0
  do i = 1, n
    yi = y(i)
    x(i) = 0.0D+00
    do j = 1, i
      ij = i0 + j
      x(j) = x(j) + yi * l(ij)
    end do
    i0 = i0 + i
  end do

  return
end subroutine ltvmul

subroutine lupdat ( beta, gamma, l, lambda, lplus, n, w, z )
!
!! LUPDAT computes lplus = secant update of L.
!
!  Discussion:
!
!    this routine updates the cholesky factor  l  of a symmetric
!    positive definite matrix to which a secant update is being
!    applied -- it computes a cholesky factor  lplus  of
!    l * (i + z*w**t) * (i + w*z**t) * l**t.  it is assumed that  w
!    and  z  have been chosen so that the updated matrix is strictly
!    positive definite.
!
!    this code uses recurrence 3 of ref. 1 (with d(j) = 1 for all j)
!    to compute  lplus  of the form  l * (i + z*w**t) * q,  where  q
!    is an orthogonal matrix that makes the result lower triangular.
!    lplus may have some negative diagonal elements.
!
!  Reference:
!
!    D Goldfarb,
!    Factorized Variable Metric Methods for Unconstrained Optimization,
!    Mathematics of Computation,
!    Volume 30, pages 796-811, 1976.
!
!  Parameters:
!
!   beta = scratch vector.
!  gamma = scratch vector.
!      l (input) lower triangular matrix, stored rowwise.
! lambda = scratch vector.
!  lplus (output) lower triangular matrix, stored rowwise, which may
!             occupy the same storage as  l.
!      n (input) length of vector parameters and order of matrices.
!      w (input, destroyed on output) right singular vector of rank 1
!             correction to  l.
!      z (input, destroyed on output) left singular vector of rank 1
!             correction to  l.
!
  integer n
  real(dp) beta(n), gamma(n), l(*), lambda(n), lplus(*), w(n), z(n)
!     dimension l(n*(n+1)/2), lplus(n*(n+1)/2)
!
  integer i, ij, j, jj, jp1, k, nm1
  integer np1
  real(dp) a, b, bj, eta, gj, lj, lij, ljj, nu, s, theta, wj, zj
  real(dp) one

  parameter (one=1.d+0 )

  nu = one
  eta = 0.0D+00
  if (n <= 1) go to 30
  nm1 = n - 1
!
!  temporarily store s(j) = sum over k = j+1 to n of w(k)**2 in
!  lambda(j).
!
  s = 0.0D+00
  do i = 1, nm1
     j = n - i
     s = s + w(j+1)**2
     lambda(j) = s
  end do
!
!  compute lambda, gamma, and beta by goldfarb*s recurrence 3.
!
  do 20 j = 1, nm1
     wj = w(j)
     a = nu*z(j) - eta*wj
     theta = one + a*wj
     s = a*lambda(j)
     lj = sqrt(theta**2 + a*s)
     if (theta > 0.0D+00) lj = -lj
     lambda(j) = lj
     b = theta*wj + s
     gamma(j) = b * nu / lj
     beta(j) = (a - b*eta) / lj
     nu = -nu / lj
     eta = -(eta + (a**2)/(theta - lj)) / lj
 20      continue
 30   lambda(n) = one + (nu*z(n) - eta*w(n))*w(n)
!
!  update l, gradually overwriting  w  and  z  with  l*w  and  l*z.
!
  np1 = n + 1
  jj = n * (n + 1) / 2

  do k = 1, n

     j = np1 - k
     lj = lambda(j)
     ljj = l(jj)
     lplus(jj) = lj * ljj
     wj = w(j)
     w(j) = ljj * wj
     zj = z(j)
     z(j) = ljj * zj
     if (k == 1) go to 50
     bj = beta(j)
     gj = gamma(j)
     ij = jj + j
     jp1 = j + 1

     do i = jp1, n
          lij = l(ij)
          lplus(ij) = lj*lij + bj*w(i) + gj*z(i)
          w(i) = w(i) + lij*wj
          z(i) = z(i) + lij*zj
          ij = ij + i
     end do

 50      jj = jj - j

  end do

  return
end subroutine lupdat

subroutine dbdog ( dig, lv, n, nwtstp, step, v )
!
!! DBDOG: compute a double dogleg step.
!
!  Discussion:
!
!    This subroutine computes a candidate step (for use in an
!    unconstrained minimization code) by the double dogleg algorithm of
!    dennis and mei (ref. 1), which is a variation on powell*s dogleg
!    scheme (ref. 2, p. 95).
!
!    let  g  and  h  be the current gradient and hessian approxima-
!    tion respectively and let d be the current scale vector.  this
!    routine assumes dig = diag(d)**-2 * g  and  nwtstp = h**-1 * g.
!    the step computed is the same one would get by replacing g and h
!    by  diag(d)**-1 * g  and  diag(d)**-1 * h * diag(d)**-1,
!    computing step, and translating step back to the original
!    variables, i.e., premultiplying it by diag(d)**-1.
!
!  Reference:
!
!    John Dennis, Howell Mei,
!    Two New Unconstrained Optimization Algorithms Which Use
!    Function and Gradient Values,
!    Journal of Optimization Theory and Applications,
!    Volume 28, pages 453-482, 1979.
!
!    M J D Powell,
!    A Hybrid Method for Non-linear Equations,
!    in Numerical Methods for Non-linear Equations,
!    edited by Philip Rabinowitz,
!    Gordon and Breach, London, 1970.
!
!  Parameters:
!
!    dig (input) diag(d)**-2 * g -- see algorithm notes.
!      g (input) the current gradient vector.
!     lv (input) length of v.
!      n (input) number of components in  dig, g, nwtstp,  and  step.
! nwtstp (input) negative newton step -- see algorithm notes.
!   step (output) the computed step.
!      v (i/o) values array, the following components of which are
!             used here...
! v(bias)   (input) bias for relaxed newton step, which is v(bias) of
!             the way from the full newton to the fully relaxed newton
!             step.  recommended value = 0.8 .
! v(dgnorm) (input) 2-norm of diag(d)**-1 * g -- see algorithm notes.
! v(dstnrm) (output) 2-norm of diag(d) * step, which is v(radius)
!             unless v(stppar) = 0 -- see algorithm notes.
! v(dst0) (input) 2-norm of diag(d) * nwtstp -- see algorithm notes.
! v(grdfac) (output) the coefficient of  dig  in the step returned --
!             step(i) = v(grdfac)*dig(i) + v(nwtfac)*nwtstp(i).
! v(gthg)   (input) square-root of (dig**t) * (hessian) * dig -- see
!             algorithm notes.
! v(gtstep) (output) inner product between g and step.
! v(nreduc) (output) function reduction predicted for the full newton
!             step.
! v(nwtfac) (output) the coefficient of  nwtstp  in the step returned --
!             see v(grdfac) above.
! v(preduc) (output) function reduction predicted for the step returned.
! v(radius) (input) the trust region radius.  d times the step returned
!             has 2-norm v(radius) unless v(stppar) = 0.
! v(stppar) (output) code telling how step was computed... 0 means a
!             full newton step.  between 0 and 1 means v(stppar) of the
!             way from the newton to the relaxed newton step.  between
!             1 and 2 means a true double dogleg step, v(stppar) - 1 of
!             the way from the relaxed newton to the Cauchy step.
!             greater than 2 means 1 / (v(stppar) - 1) times the Cauchy
!             step.
!
  integer lv
  integer n

  real(dp) dig(n), nwtstp(n), step(n), v(lv)
!  external dotprd, v2norm
!  real(dp) dotprd, v2norm
  real(dp) cfact, cnorm, ctrnwt, ghinvg, femnsq, gnorm
  real(dp) nwtnrm, relax, rlambd, t, t1, t2
  real(dp) half, two
  integer bias, dgnorm, dstnrm, dst0, grdfac, gthg, gtstep
  integer nreduc, nwtfac, preduc, radius, stppar
  parameter (half=0.5d+0, two=2.d+0)
  parameter (bias=43, dgnorm=1, dstnrm=2, dst0=3, grdfac=45 )
  parameter ( gthg=44, gtstep=4, nreduc=6, nwtfac=46, preduc=7 )
  parameter ( radius=8, stppar=5)

  nwtnrm = v(dst0)
  rlambd = 1.0D+00
  if (nwtnrm > 0.0D+00 ) rlambd = v(radius) / nwtnrm
  gnorm = v(dgnorm)
  ghinvg = two * v(nreduc)
  v(grdfac) = 0.0D+00
  v(nwtfac) = 0.0D+00
  if (rlambd < 1.0D+00 ) go to 30
!
!  The Newton step is inside the trust region.
!
     v(stppar) = 0.0D+00
     v(dstnrm) = nwtnrm
     v(gtstep) = -ghinvg
     v(preduc) = v(nreduc)
     v(nwtfac) = -1.0D+00
     step(1:n) = -nwtstp(1:n)
     return

 30   v(dstnrm) = v(radius)
  cfact = (gnorm / v(gthg))**2
!
!  Cauchy step = -cfact * g.
!
  cnorm = gnorm * cfact
  relax = 1.0D+00 - v(bias) * ( 1.0D+00 - gnorm*cnorm/ghinvg)
  if (rlambd < relax) go to 50
!
!  Step is between relaxed Newton and full Newton steps.
!
     v(stppar) = 1.0D+00 -  (rlambd - relax) / ( 1.0D+00 - relax)
     t = -rlambd
     v(gtstep) = t * ghinvg
     v(preduc) = rlambd * ( 1.0D+00 - half*rlambd) * ghinvg
     v(nwtfac) = t
     step(1:n) = t * nwtstp(1:n)
     return

 50   if (cnorm < v(radius)) go to 70
!
!  The Cauchy step lies outside the trust region --
!  step = scaled Cauchy step.
!
     t = -v(radius) / gnorm
     v(grdfac) = t
     v(stppar) = 1.0D+00  +  cnorm / v(radius)
     v(gtstep) = -v(radius) * gnorm
  v(preduc) = v(radius)*(gnorm - half*v(radius)*(v(gthg)/gnorm)**2)
     step(1:n) = t * dig(1:n)
     return
!
!  Compute dogleg step between Cauchy and relaxed Newton
!  femur = relaxed newton step minus Cauchy step.
!
 70   ctrnwt = cfact * relax * ghinvg / gnorm
!
!  ctrnwt = inner product of Cauchy and relaxed Newton steps,
!  scaled by gnorm**-1.
!
  t1 = ctrnwt - gnorm*cfact**2
!
!  t1 = inner prod. of femur and Cauchy step, scaled by gnorm**-1.
!
  t2 = v(radius)*(v(radius)/gnorm) - gnorm*cfact**2
  t = relax * nwtnrm
  femnsq = (t/gnorm)*t - ctrnwt - t1
!
!  femnsq = square of 2-norm of femur, scaled by gnorm**-1.
!
  t = t2 / (t1 + sqrt(t1**2 + femnsq*t2))
!
!  Dogleg step  =  Cauchy step  +  t * femur.
!
  t1 = (t - 1.0D+00 ) * cfact
  v(grdfac) = t1
  t2 = -t * relax
  v(nwtfac) = t2
  v(stppar) = two - t
  v(gtstep) = t1*gnorm**2 + t2*ghinvg
  v(preduc) = -t1*gnorm * ((t2 + 1.0D+00 )*gnorm) &
                  - t2 * ( 1.0D+00 + half*t2)*ghinvg &
                   - half * (v(gthg)*t1)**2

  step(1:n) = t1 * dig(1:n) + t2 * nwtstp(1:n)

  return
end subroutine dbdog

function dotprd ( p, x, y ) result(dtpr)
!
!! DOTPRD returns the inner product of vectors X and Y.
!
  integer p
  integer i
  real(dp) dtpr
!  real(dp) rmdcon
  real(dp), save :: sqteta = 0.0D+00
  real(dp) t
  real(dp) x(p)
  real(dp) y(p)

  dtpr = 0.0D+00

  if ( sqteta == 0.0D+00 ) then
    sqteta = rmdcon(2)
  end if

  do i = 1, p

    t = max ( abs ( x(i) ), abs ( y(i) ) )
    if ( t > 1.0D+00 ) go to 10
    if (t < sqteta) go to 20
    t = (x(i)/sqteta)*y(i)
    if (abs(t) < sqteta) go to 20
 10   dtpr = dtpr + x(i)*y(i)

 20 continue

  end do

end function dotprd

function stopx ( ) result(stpx)
!
!! STOPX checks to see if the BREAK key has been pressed.
!
!  Discussion:
!
!     this function may serve as the stopx (asynchronous interruption)
!     function for the nl2sol (nonlinear least-squares) package at
!     those installations which do not wish to implement a
!     dynamic stopx.
!
!     at installations where the nl2sol system is used
!     interactively, this dummy stopx should be replaced by a
!     function that returns .true. if and only if the interrupt
!     (break) key has been pressed since the last call on stopx.
!
  logical stpx

  stpx = .false.

end function stopx

subroutine vaxpy ( p, w, a, x, y )
!
!! VAXPY sets w = a*x + y.
!
!  Discussion:
!
!    w, x, y = p-vectors, a = scalar
!
!  Parameters:
!
  implicit none

  integer p

  real(dp) a
  real(dp) w(p)
  real(dp) x(p)
  real(dp) y(p)

  w(1:p) = a * x(1:p) + y(1:p)

  return
end subroutine vaxpy

function reldst ( p, d, x, x0 ) result(rreldst)
!
!! RELDST computes the relative difference between X and X0.
!
  integer p

  real(dp) rreldst
  real(dp) d(p), x(p), x0(p)
  integer i
  real(dp) emax, t, xmax

  emax = 0.0D+00
  xmax = 0.0D+00

  do i = 1, p
    t = abs(d(i) * (x(i) - x0(i)))
    if (emax < t) emax = t
    t = d(i) * (abs(x(i)) + abs(x0(i)))
    if (xmax < t) xmax = t
  end do

  rreldst = 0.0D+00
  if ( xmax > 0.0D+00 ) rreldst = emax / xmax

end function reldst

subroutine assst ( iv, liv, lv, v )
!
!! ASSST assesses a candidate step.
!
!  Discussion:
!
!    This subroutine is called by an unconstrained minimization
!    routine to assess the next candidate step.  it may recommend one
!    of several courses of action, such as accepting the step, recom-
!    puting it using the same or a new quadratic model, or halting due
!    to convergence or false convergence.  See the return code listing
!    below.
!
!  Reference:
!
!    John Dennis, David Gay, Roy Welsch,
!    An Adaptive Nonlinear Least-squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, 1981.
!
!    M J D Powell,
!    A Fortran Subroutine for Solving Systems of Nonlinear Algebraic Equations,
!    in Numerical Methods for Nonlinear Algebraic Equations,
!    edited by Philip Rabinowitz,
!    Gordon and Breach, London, 1970.
!
!  Parameters:
!
!  iv (i/o) integer parameter and scratch vector -- see description
!             below of iv values referenced.
!
! liv (in)  length of iv array.
!
!  lv (in)  length of v array.
!
!   v (i/o) real parameter and scratch vector -- see description
!             below of v values referenced.
!
!   iv values referenced
!
!    iv(irc) (i/o) on input for the first step tried in a new iteration,
!             iv(irc) should be set to 3 or 4 (the value to which it is
!             set when step is definitely to be accepted).  on input
!             after step has been recomputed, iv(irc) should be
!             unchanged since the previous return of assst.
!                on output, iv(irc) is a return code having one of the
!             following values...
!                  1 = switch models or try smaller step.
!                  2 = switch models or accept step.
!                  3 = accept step and determine v(radfac) by gradient
!                       tests.
!                  4 = accept step, v(radfac) has been determined.
!                  5 = recompute step (using the same model).
!                  6 = recompute step with radius = v(lmaxs) but do not
!                       evaulate the objective function.
!                  7 = x-convergence (see v(xctol)).
!                  8 = relative function convergence (see v(rfctol)).
!                  9 = both x- and relative function convergence.
!                 10 = absolute function convergence (see v(afctol)).
!                 11 = singular convergence (see v(lmaxs)).
!                 12 = false convergence (see v(xftol)).
!                 13 = iv(irc) was out of range on input.
!             return code i has precdence over i+1 for i = 9, 10, 11.
! iv(mlstgd) (i/o) saved value of iv(model).
!  iv(model) (i/o) on input, iv(model) should be an integer identifying
!             the current quadratic model of the objective function.
!             if a previous step yielded a better function reduction,
!             then iv(model) will be set to iv(mlstgd) on output.
! iv(nfcall) (in)  invocation count for the objective function.
! iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
!             function reduction this iteration.  iv(nfgcal) remains
!             unchanged until a function reduction is obtained.
! iv(radinc) (i/o) the number of radius increases (or minus the number
!             of decreases) so far this iteration.
! iv(restor) (out) set to 1 if v(f) has been restored and x should be
!             restored to its initial value, to 2 if x should be saved,
!             to 3 if x should be restored from the saved value, and to
!             0 otherwise.
!  iv(stage) (i/o) count of the number of models tried so far in the
!             current iteration.
! iv(stglim) (in)  maximum number of models to consider.
! iv(switch) (out) set to 0 unless a new model is being tried and it
!             gives a smaller function value than the previous model,
!             in which case assst sets iv(switch) = 1.
! iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
!             overflow).
!   iv(xirc) (i/o) value that iv(irc) would have in the absence of
!             convergence, false convergence, and oversized steps.
!
!   v values referenced
!
! v(afctol) (in)  absolute function convergence tolerance.  if the
!             absolute value of the current function value v(f) is less
!             than v(afctol), then assst returns with iv(irc) = 10.
! v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
!             nonzero.
! v(dstnrm) (in)  the 2-norm of d*step.
! v(dstsav) (i/o) value of v(dstnrm) on saved step.
!   v(dst0) (in)  the 2-norm of d times the newton step (when defined,
!             i.e., for v(nreduc) >= 0).
!      v(f) (i/o) on both input and output, v(f) is the objective func-
!             tion value at x.  if x is restored to a previous value,
!             then v(f) is restored to the corresponding value.
!   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
!             value of v(f) if an earlier step gave a bigger function
!             decrease, and for the input value of v(f) otherwise).
! v(flstgd) (i/o) saved value of v(f).
!     v(f0) (in)  objective function value at start of iteration.
! v(gtslst) (i/o) value of v(gtstep) on saved step.
! v(gtstep) (in)  inner product between step and gradient.
! v(incfac) (in)  minimum factor by which to increase radius.
!  v(lmaxs) (in)  maximum reasonable step size (and initial step bound).
!             if the actual function decrease is no more than twice
!             what was predicted, if a return with iv(irc) = 7, 8, 9,
!             or 10 does not occur, if v(dstnrm) > v(lmaxs), and if
!             v(preduc) <= v(sctol) * abs(v(f0)), then assst re-
!             turns with iv(irc) = 11.  if so doing appears worthwhile,
!             then assst repeats this test with v(preduc) computed for
!             a step of length v(lmaxs) (by a return with iv(irc) = 6).
! v(nreduc) (i/o)  function reduction predicted by quadratic model for
!             newton step.  if assst is called with iv(irc) = 6, i.e.,
!             if v(preduc) has been computed with radius = v(lmaxs) for
!             use in the singular convervence test, then v(nreduc) is
!             set to -v(preduc) before the latter is restored.
! v(plstgd) (i/o) value of v(preduc) on saved step.
! v(preduc) (i/o) function reduction predicted by quadratic model for
!             current step.
! v(radfac) (out) factor to be used in determining the new radius,
!             which should be v(radfac)*dst, where  dst  is either the
!             output value of v(dstnrm) or the 2-norm of
!             diag(newd)*step  for the output value of step and the
!             updated version, newd, of the scale vector d.  for
!             iv(irc) = 3, v(radfac) = 1.0D+00 is returned.
! v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
!             value of v(dstnrm) -- suggested value = 0.1.
! v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
!  v(reldx) (in) scaled relative change in x caused by step, computed
!             (e.g.) by function  reldst  as
!                 max (d(i)*abs(x(i)-x0(i)), 1 <= i <= p) /
!                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p).
! v(rfctol) (in)  relative function convergence tolerance.  if the
!             actual function reduction is at most twice what was pre-
!             dicted and  v(nreduc) <= v(rfctol)*abs(v(f0)),  then
!             assst returns with iv(irc) = 8 or 9.
! v(stppar) (in)  marquardt parameter -- 0 means full newton step.
! v(tuner1) (in)  tuning constant used to decide if the function
!             reduction was much less than expected.  suggested
!             value = 0.1.
! v(tuner2) (in)  tuning constant used to decide if the function
!             reduction was large enough to accept step.  suggested
!             value = 10**-4.
! v(tuner3) (in)  tuning constant used to decide if the radius
!             should be increased.  suggested value = 0.75.
!  v(xctol) (in)  x-convergence criterion.  if step is a newton step
!             (v(stppar) = 0) having v(reldx) <= v(xctol) and giving
!             at most twice the predicted function decrease, then
!             assst returns iv(irc) = 7 or 9.
!  v(xftol) (in)  false convergence tolerance.  if step gave no or only
!             a small function decrease and v(reldx) <= v(xftol),
!             then assst returns with iv(irc) = 12.
!
!  notes
!
!   application and usage restrictions
!
!        this routine is called as part of the nl2sol (nonlinear
!     least-squares) package.  it may be used in any unconstrained
!     minimization solver that uses dogleg, goldfeld-quandt-trotter,
!     or levenberg-marquardt steps.
!
!   algorithm notes
!
!        see (1) for further discussion of the assessing and model
!     switching strategies.  while nl2sol considers only two models,
!     assst is designed to handle any number of models.
!
!   usage notes
!
!        on the first call of an iteration, only the i/o variables
!     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and
!     v(preduc) need have been initialized.  between calls, no i/o
!     values execpt step, x, iv(model), v(f) and the stopping toler-
!     ances should be changed.
!        after a return for convergence or false convergence, one can
!     change the stopping tolerances and call assst again, in which
!     case the stopping tests will be repeated.
!
!   history
!
!        john dennis designed much of this routine, starting with
!     ideas in (2). roy welsch suggested the model switching strategy.
!        david gay and stephen peters cast this subroutine into a more
!     portable form (winter 1977), and david gay cast it into its
!     present form (fall 1978).
!
  integer liv
  integer lv

  integer iv(liv)
  real(dp) v(lv)
  logical goodx
  integer i, nfc
  real(dp) emax, emaxs, gts, rfac1, xmax
  real(dp) half, one, onep2, two
  integer afctol, decfac, dstnrm, dstsav, dst0, f, fdif, flstgd, f0
  integer gtslst, gtstep, incfac, irc, lmaxs, mlstgd, model, nfcall
  integer nfgcal, nreduc, plstgd, preduc, radfac, radinc, rdfcmn
  integer rdfcmx, reldx, restor, rfctol, sctol, stage, stglim
  integer stppar, switch, toobig, tuner1, tuner2, tuner3, xctol
  integer xftol, xirc

  parameter ( half=0.5d+0, one=1.d+0, onep2=1.2d+0, two=2.d+0)
  parameter ( irc=29, mlstgd=32, model=5, nfcall=6, nfgcal=7 )
  parameter ( radinc=8, restor=9, stage=10, stglim=11, switch=12 )
  parameter ( toobig=2, xirc=13)
  parameter (afctol=31, decfac=22, dstnrm=2, dst0=3, dstsav=18 )
  parameter (f=10, fdif=11, flstgd=12, f0=13, gtslst=14, gtstep=4 )
  parameter (incfac=23, lmaxs=36, nreduc=6, plstgd=15, preduc=7 )
  parameter (radfac=16, rdfcmn=24, rdfcmx=25, reldx=17, rfctol=32 )
  parameter (sctol=37, stppar=5, tuner1=26, tuner2=27, tuner3=28 )
  parameter (xctol=33, xftol=34)

  nfc = iv(nfcall)
  iv(switch) = 0
  iv(restor) = 0
  rfac1 = one
  goodx = .true.
  i = iv(irc)

  if (i >= 1 .and. i <= 12) then
        go to (20,30,10,10,40,280,220,220,220,220,220,170), i
  end if

  iv(irc) = 13
  return
!
!  Initialize for new iteration.
!
 10   iv(stage) = 1
  iv(radinc) = 0
  v(flstgd) = v(f0)
  if (iv(toobig) == 0) go to 110
     iv(stage) = -1
     iv(xirc) = i
     go to 60
!
!  Step was recomputed with new model or smaller radius
!  first decide which
!
 20   if (iv(model) /= iv(mlstgd)) go to 30
!
!  Old model retained, smaller radius tried
!  do not consider any more new models this iteration
!
     iv(stage) = iv(stglim)
     iv(radinc) = -1
     go to 110
!
!  A new model is being tried.  decide whether to keep it.
!
 30   iv(stage) = iv(stage) + 1
!
!  Now we add the possibiltiy that step was recomputed with
!  the same model, perhaps because of an oversized step.
!
 40   if (iv(stage) > 0) go to 50
!
!  Step was recomputed because it was too big.
!
     if (iv(toobig) /= 0) go to 60
!
!  Restore iv(stage) and pick up where we left off.
!
     iv(stage) = -iv(stage)
     i = iv(xirc)
     go to (20, 30, 110, 110, 70), i

 50   if (iv(toobig) == 0) go to 70
!
!  Handle oversize step
!
  if (iv(radinc) > 0) go to 80
     iv(stage) = -iv(stage)
     iv(xirc) = iv(irc)

 60      v(radfac) = v(decfac)
     iv(radinc) = iv(radinc) - 1
     iv(irc) = 5
     iv(restor) = 1
     return

 70   if (v(f) < v(flstgd)) go to 110
!
!  The new step is a loser.  restore old model.
!
  if (iv(model) == iv(mlstgd)) go to 80
     iv(model) = iv(mlstgd)
     iv(switch) = 1
!
!  Restore step, etc. only if a previous step decreased v(f).
!
 80   if (v(flstgd) >= v(f0)) go to 110
     iv(restor) = 1
     v(f) = v(flstgd)
     v(preduc) = v(plstgd)
     v(gtstep) = v(gtslst)
     if (iv(switch) == 0) rfac1 = v(dstnrm) / v(dstsav)
     v(dstnrm) = v(dstsav)
     nfc = iv(nfgcal)
     goodx = .false.

 110  v(fdif) = v(f0) - v(f)
  if (v(fdif) > v(tuner2) * v(preduc)) go to 140
  if(iv(radinc)>0) go to 140
!
!         no (or only a trivial) function decrease
!         so try new model or smaller radius
!
     if (v(f) < v(f0)) go to 120
          iv(mlstgd) = iv(model)
          v(flstgd) = v(f)
          v(f) = v(f0)
          iv(restor) = 1
          go to 130
 120     iv(nfgcal) = nfc
 130     iv(irc) = 1
     if (iv(stage) < iv(stglim)) go to 160
          iv(irc) = 5
          iv(radinc) = iv(radinc) - 1
          go to 160
!
!  Nontrivial function decrease achieved
!
 140  iv(nfgcal) = nfc
  rfac1 = 1.0D+00
  v(dstsav) = v(dstnrm)
  if (v(fdif) > v(preduc)*v(tuner1)) go to 190
!
!  Decrease was much less than predicted -- either change models
!  or accept step with decreased radius.
!
  if (iv(stage) >= iv(stglim)) go to 150
!
!  Consider switching models
!
     iv(irc) = 2
     go to 160
!
!  Accept step with decreased radius
!
 150  iv(irc) = 4
!
!   set v(radfac) to fletcher*s decrease factor
!
 160  iv(xirc) = iv(irc)
  emax = v(gtstep) + v(fdif)
  v(radfac) = half * rfac1

  if (emax < v(gtstep)) then
    v(radfac) = rfac1 * max (v(rdfcmn),half * v(gtstep)/emax)
  end if
!
!  Do false convergence test
!
 170  if (v(reldx) <= v(xftol)) go to 180
     iv(irc) = iv(xirc)
     if (v(f) < v(f0)) go to 200
          go to 230

 180  iv(irc) = 12
  go to 240
!
!  Handle good function decrease
!
 190  if (v(fdif) < (-v(tuner3) * v(gtstep))) go to 210
!
!  Increasing radius looks worthwhile.  see if we just
!  recomputed step with a decreased radius or restored step
!  after recomputing it with a larger radius.
!
  if (iv(radinc) < 0) go to 210
  if (iv(restor) == 1) go to 210
!
!  We did not.  try a longer step unless this was a newton step.
!
     v(radfac) = v(rdfcmx)
     gts = v(gtstep)
     if (v(fdif) < (half/v(radfac) - 1.0D+00 ) * gts) then
       v(radfac) = max (v(incfac), half*gts/(gts + v(fdif)))
     end if
     iv(irc) = 4
     if (v(stppar) == 0.0D+00 ) go to 230
     if (v(dst0) >= 0.0D+00 .and. (v(dst0) < two*v(dstnrm) &
              .or. v(nreduc) < onep2*v(fdif)))  then
       go to 230
     end if
!
!  Step was not a newton step.  recompute it with a larger radius.
!
          iv(irc) = 5
          iv(radinc) = iv(radinc) + 1
!
!  Save values corresponding to good step
!
 200  v(flstgd) = v(f)
  iv(mlstgd) = iv(model)
  if (iv(restor) /= 1) iv(restor) = 2
  v(dstsav) = v(dstnrm)
  iv(nfgcal) = nfc
  v(plstgd) = v(preduc)
  v(gtslst) = v(gtstep)
  go to 230
!
!  Accept step with radius unchanged.
!
 210  v(radfac) = 1.0D+00
  iv(irc) = 3
  go to 230
!
!  Come here for a restart after convergence.
!
 220  iv(irc) = iv(xirc)
  if (v(dstsav) >= 0.0D+00 ) go to 240
     iv(irc) = 12
     go to 240
!
!  Perform convergence tests.
!
 230  iv(xirc) = iv(irc)
 240  if (iv(restor) == 1 .and. v(flstgd) < v(f0)) iv(restor) = 3
  if (abs(v(f)) < v(afctol)) iv(irc) = 10

  if (half * v(fdif) > v(preduc)) then
    return
  end if

  emax = v(rfctol) * abs(v(f0))
  emaxs = v(sctol) * abs(v(f0))
  if (v(dstnrm) > v(lmaxs) .and. v(preduc) <= emaxs) then
    iv(irc) = 11
  end if
  if (v(dst0) < 0.0D+00 ) go to 250
  i = 0

  if ((v(nreduc) > 0.0D+00 .and. v(nreduc) <= emax) .or. &
      (v(nreduc) == 0.0D+00 .and. v(preduc) == 0.0D+00 )) then
    i = 2
  end if

  if (v(stppar) == 0.0D+00 .and. v(reldx) <= v(xctol) .and. goodx) then
    i = i + 1
  end if

  if (i > 0) iv(irc) = i + 6
!
!  Consider recomputing step of length v(lmaxs) for singular
!  convergence test.
!
 250  if (iv(irc) > 5 .and. iv(irc) /= 12) then
     return
  end if

  if (v(dstnrm) > v(lmaxs)) go to 260
     if (v(preduc) >= emaxs) then
       return
     end if
          if (v(dst0) <= 0.0D+00 ) go to 270
               if (half * v(dst0) <= v(lmaxs)) then
                 return
               end if
                    go to 270
 260  if (half * v(dstnrm) <= v(lmaxs)) then
        return
      end if
  xmax = v(lmaxs) / v(dstnrm)
  if (xmax * (two - xmax) * v(preduc) >= emaxs) then
    return
  end if
 270  if (v(nreduc) < 0.0D+00 ) go to 290
!
!   recompute v(preduc) for use in singular convergence test
!
  v(gtslst) = v(gtstep)
  v(dstsav) = v(dstnrm)
  if (iv(irc) == 12) v(dstsav) = -v(dstsav)
  v(plstgd) = v(preduc)
  i = iv(restor)
  iv(restor) = 2
  if (i == 3) iv(restor) = 0
  iv(irc) = 6
  return
!
!  Perform singular convergence test with recomputed v(preduc)
!
 280  v(gtstep) = v(gtslst)
  v(dstnrm) = abs(v(dstsav))
  iv(irc) = iv(xirc)
  if (v(dstsav) <= 0.0D+00 ) iv(irc) = 12
  v(nreduc) = -v(preduc)
  v(preduc) = v(plstgd)
  iv(restor) = 3

 290  if (-v(nreduc) <= v(rfctol) * abs(v(f0))) iv(irc) = 11

  return
end subroutine assst

subroutine wzbfgs ( l, n, s, w, y, z )
!
!! WZBFGS compute Y and Z for LUPDAT corresponding to BFGS update.
!
!  Discussion:
!
!    When S is computed in certain ways, for example by GQTSTP or
!    DBLDOG, it is possible to save N**2/2 operations since L'*S
!    or L*L'*S is then known.
!
!    If the BFGS update to L*L' would reduce its determinant to
!    less than EPS times its old value, then this routine in effect
!    replaces Y by THETA*Y + (1-THETA)*L*L'*S, where THETA
!    (between 0 and 1) is chosen to make the reduction factor = EPS.
!
!  Parameters:
!
!    l (i/o) cholesky factor of hessian, a lower triang. matrix stored
!             compactly by rows.
!
!    n (input) order of  l  and length of  s,  w,  y,  z.
!
!    s (input) the step just taken.
!
!    w (output) right singular vector of rank 1 correction to l.
!
!    y (input) change in gradients corresponding to s.
!
!    z (output) left singular vector of rank 1 correction to l.
!
  implicit none

  integer n

!  real(dp) dotprd
  real(dp) cs
  real(dp) cy
  real(dp), parameter :: eps = 0.1D+00
  real(dp) epsrt
  real(dp) l(n*(n+1)/2)
  real(dp) s(n)
  real(dp) shs
  real(dp) theta
  real(dp) w(n)
  real(dp) y(n)
  real(dp) ys
  real(dp) z(n)

  call ltvmul ( n, w, l, s )
  shs = dotprd ( n, w, w )
  ys = dotprd ( n, y, s )

  if ( ys < eps * shs ) then
    theta = ( 1.0D+00 - eps ) * shs / ( shs - ys )
    epsrt = sqrt ( eps )
    cy = theta / ( shs * epsrt )
    cs = ( 1.0D+00 + ( theta - 1.0D+00 ) / epsrt ) / shs
  else
    cy = 1.0D+00 / ( sqrt ( ys ) * sqrt ( shs ) )
    cs = 1.0D+00 / shs
  end if

  call livmul ( n, z, l, y )

  z(1:n) = cy * z(1:n) - cs * w(1:n)

  return
end subroutine wzbfgs

end module m_cg
