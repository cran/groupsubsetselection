MODULE lsq

!  Module for unconstrained linear least-squares calculations.
!  The algorithm is suitable for updating LS calculations as more
!  data are added.   This is sometimes called recursive estimation.
!  Only one dependent variable is allowed.
!  Based upon Applied Statistics algorithm AS 274.
!  Translation from Fortran 77 to Fortran 90 by Alan Miller.
!  A function, VARPRD, has been added for calculating the variances
!  of predicted values, and this uses a subroutine BKSUB2.

!  Version 1.14, 19 August 2002 - ELF90 compatible version
!  Author: Alan Miller
!  e-mail : amiller @ bigpond.net.au
!  WWW-pages: http://www.ozemail.com.au/~milleraj
!             http://users.bigpond.net.au/amiller/

!  Bug fixes:
!  1. In REGCF a call to TOLSET has been added in case the user had
!     not set tolerances.
!  2. In SING, each time a singularity is detected, unless it is in the
!     variables in the last position, INCLUD is called.   INCLUD assumes
!     that a new observation is being added and increments the number of
!     cases, NOBS.   The line:  nobs = nobs - 1 has been added.
!  3. row_ptr was left out of the DEALLOCATE statement in routine startup
!     in version 1.07.
!  4. In COV, now calls SS if rss_set = .FALSE.  29 August 1997
!  5. In TOLSET, correction to accomodate negative values of D.  19 August 2002

!  Other changes:
!  1. Array row_ptr added 18 July 1997.   This points to the first element
!     stored in each row thus saving a small amount of time needed to
!     calculate its position.
!  2. Optional parameter, EPS, added to routine TOLSET, so that the user
!     can specify the accuracy of the input data.
!  3. Cosmetic change of lsq_kind to dp (`Double precision')
!  4. Change to routine SING to use row_ptr rather than calculate the position
!     of first elements in each row.

!  The PUBLIC variables are:
!  dp       = a KIND parameter for the floating-point quantities calculated
!             in this module.   See the more detailed explanation below.
!             This KIND parameter should be used for all floating-point
!             arguments passed to routines in this module.

!  nobs    = the number of observations processed to date.
!  ncol    = the total number of variables, including one for the constant,
!            if a constant is being fitted.
!  r_dim   = the dimension of array r = ncol*(ncol-1)/2
!  vorder  = an integer vector storing the current order of the variables
!            in the QR-factorization.   The initial order is 0, 1, 2, ...
!            if a constant is being fitted, or 1, 2, ... otherwise.
!  initialized = a logical variable which indicates whether space has
!                been allocated for various arrays.
!  tol_set = a logical variable which is set when subroutine TOLSET has
!            been called to calculate tolerances for use in testing for
!            singularities.
!  rss_set = a logical variable indicating whether residual sums of squares
!            are available and usable.
!  d()     = array of row multipliers for the Cholesky factorization.
!            The factorization is X = Q.sqrt(D).R where Q is an ortho-
!            normal matrix which is NOT stored, D is a diagonal matrix
!            whose diagonal elements are stored in array d, and R is an
!            upper-triangular matrix with 1's as its diagonal elements.
!  rhs()   = vector of RHS projections (after scaling by sqrt(D)).
!            Thus Q'y = sqrt(D).rhs
!  r()     = the upper-triangular matrix R.   The upper triangle only,
!            excluding the implicit 1's on the diagonal, are stored by
!            rows.
!  tol()   = array of tolerances used in testing for singularities.
!  rss()   = array of residual sums of squares.   rss(i) is the residual
!            sum of squares with the first i variables in the model.
!            By changing the order of variables, the residual sums of
!            squares can be found for all possible subsets of the variables.
!            The residual sum of squares with NO variables in the model,
!            that is the total sum of squares of the y-values, can be
!            calculated as rss(1) + d(1)*rhs(1)^2.   If the first variable
!            is a constant, then rss(1) is the sum of squares of
!            (y - ybar) where ybar is the average value of y.
!  sserr   = residual sum of squares with all of the variables included.
!  row_ptr() = array of indices of first elements in each row of R.
!
!--------------------------------------------------------------------------

!     General declarations

IMPLICIT NONE

INTEGER, SAVE                :: nobs, ncol, r_dim
INTEGER, ALLOCATABLE, SAVE   :: vorder(:), row_ptr(:), corder(:)
LOGICAL, SAVE                :: initialized = .false.,                  &
                                tol_set = .false., rss_set = .false.

! Note. dp is being set to give at least 12 decimal digit
!       representation of floating point numbers.   This should be adequate
!       for most problems except the fitting of polynomials.   dp is
!       being set so that the same code can be run on PCs and Unix systems,
!       which will usually represent floating-point numbers in `double
!       precision', and other systems with larger word lengths which will
!       give similar accuracy in `single precision'.

INTEGER, PARAMETER           :: dp = SELECTED_REAL_KIND(12,60)
REAL (dp), ALLOCATABLE, SAVE :: d(:), rhs(:), r(:), tol(:), rss(:), &
				wtslb(:) ! Weights lower bounds for each weight.
INTEGER :: ng  !Yi
integer, allocatable, save :: gorder(:), g_ptr(:), gnv(:)
REAL (dp), SAVE              :: zero = 0.0_dp, one = 1.0_dp, vsmall
REAL (dp), SAVE              :: sserr, toly
INTEGER offset,nforce,option,dbfid!debug
PUBLIC                       :: dp, nobs, ncol, r_dim, vorder, row_ptr, &
                                initialized, tol_set, rss_set,          &
                                d, rhs, r, tol, rss, sserr, offset, nforce, option, &
                                ng, gnv, g_ptr, gorder, corder,dbfid ! Added by Yi Guo for GSS

PRIVATE                      :: zero, one, vsmall


CONTAINS

SUBROUTINE startup(nvar, fit_const)

!     Allocates dimensions for arrays and initializes to zero
!     The calling program must set nvar = the number of variables, and
!     fit_const = .true. if a constant is to be included in the model,
!     otherwise fit_const = .false.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)  :: nvar
LOGICAL, INTENT(IN)  :: fit_const

!     Local variable
INTEGER   :: i

vsmall = 10. * TINY(zero)

nobs = 0
IF (fit_const) THEN
  ncol = nvar + 1
  offset=nforce+1
ELSE
  ncol = nvar
  offset=nforce
END IF

IF (initialized) DEALLOCATE(d, rhs, r, tol, rss, vorder, row_ptr, gorder, g_ptr) !Yi: add gorder and g_ptr
r_dim = ncol * (ncol - 1)/2

!ALLOCATE( d(ncol), rhs(ncol), r(r_dim), tol(ncol), rss(ncol), vorder(ncol), row_ptr(ncol) )    !Yi
ALLOCATE( d(ncol), rhs(ncol), r(r_dim), tol(ncol), rss(ncol), vorder(ncol), row_ptr(ncol),gorder(ng), g_ptr(ng+1))    !Yi
! Yi
! Get g_ptr
g_ptr = 1
gorder(1) =1
if (ng>=2) then
    do i = 2, ng
        g_ptr(i) = g_ptr(i-1) + gnv(i-1)
        gorder(i) = i
    end do
end if
g_ptr(ng+1) = g_ptr(ng) + gnv(ng)
g_ptr = g_ptr + nforce
! Get g_ptr done.

d = zero
rhs = zero
r = zero
sserr = zero

IF (fit_const) THEN
  DO i = 1, ncol
    vorder(i) = i-1
  END DO
ELSE
  DO i = 1, ncol
    vorder(i) = i
  END DO
END IF ! (fit_const)

! row_ptr(i) is the position of element R(i,i+1) in array r().

row_ptr(1) = 1
DO i = 2, ncol-1
  row_ptr(i) = row_ptr(i-1) + ncol - i + 1
END DO
row_ptr(ncol) = 0

initialized = .true.
tol_set = .false.
rss_set = .false.

RETURN
END SUBROUTINE startup




SUBROUTINE includ(weight, xrow, yelem)

!     ALGORITHM AS75.1  APPL. STATIST. (1974) VOL.23, NO. 3

!     Calling this routine updates D, R, RHS and SSERR by the
!     inclusion of xrow, yelem with the specified weight.

!     *** WARNING  Array XROW is overwritten.

!     N.B. As this routine will be called many times in most applications,
!          checks have been eliminated.
!
!--------------------------------------------------------------------------


IMPLICIT NONE
REAL (dp),INTENT(IN)                    :: weight, yelem
REAL (dp), DIMENSION(:), INTENT(IN OUT) :: xrow

!     Local variables

INTEGER     :: i, k, nextr
REAL (dp)   :: w, y, xi, di, wxi, dpi, cbar, sbar, xk

nobs = nobs + 1
w = weight
y = yelem
rss_set = .false.
nextr = 1
DO i = 1, ncol

!     Skip unnecessary transformations.   Test on exact zeroes must be
!     used or stability can be destroyed.

  IF (ABS(w) < vsmall) RETURN
  xi = xrow(i)
  IF (ABS(xi) < vsmall) THEN
    nextr = nextr + ncol - i
  ELSE
    di = d(i)
    wxi = w * xi
    dpi = di + wxi*xi
    cbar = di / dpi
    sbar = wxi / dpi
    w = cbar * w
    d(i) = dpi
    DO k = i+1, ncol
      xk = xrow(k)
      xrow(k) = xk - xi * r(nextr)
      r(nextr) = cbar * r(nextr) + sbar * xk
      nextr = nextr + 1
    END DO
    xk = y
    y = xk - xi * rhs(i)
    rhs(i) = cbar * rhs(i) + sbar * xk
  END IF
END DO ! i = 1, ncol

!     Y * SQRT(W) is now equal to the Brown, Durbin & Evans recursive
!     residual.

sserr = sserr + w * y * y

RETURN
END SUBROUTINE includ



SUBROUTINE regcf(beta, nreq, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Modified version of AS75.4 to calculate regression coefficients
!     for the first NREQ variables, given an orthogonal reduction from
!     AS75.1.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)                  :: nreq
INTEGER, INTENT(OUT)                 :: ifault
REAL (dp), DIMENSION(:), INTENT(OUT) :: beta

!     Local variables

INTEGER   :: i, j, nextr

!     Some checks.

ifault = 0
IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
IF (ifault /= 0) RETURN

IF (.NOT. tol_set) CALL tolset()

DO i = nreq, 1, -1
  IF (SQRT(d(i)) < tol(i)) THEN
    beta(i) = zero
    d(i) = zero
    ifault = -i
  ELSE
    beta(i) = rhs(i)
    nextr = row_ptr(i)
    DO j = i+1, nreq
      beta(i) = beta(i) - r(nextr) * beta(j)
      nextr = nextr + 1
    END DO ! j = i+1, nreq
  END IF
END DO ! i = nreq, 1, -1

RETURN
END SUBROUTINE regcf



SUBROUTINE tolset(eps)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Sets up array TOL for testing for zeroes in an orthogonal
!     reduction formed using AS75.1.

REAL (dp), INTENT(IN), OPTIONAL :: eps

!     Unless the argument eps is set, it is assumed that the input data are
!     recorded to full machine accuracy.   This is often not the case.
!     If, for instance, the data are recorded to `single precision' of about
!     6-7 significant decimal digits, then singularities will not be detected.
!     It is suggested that in this case eps should be set equal to
!     10.0 * EPSILON(1.0)
!     If the data are recorded to say 4 significant decimals, then eps should
!     be set to 1.0E-03
!     The above comments apply to the predictor variables, not to the
!     dependent variable.

!     Correction - 19 August 2002
!     When negative weights are used, it is possible for an alement of D
!     to be negative.

!     Local variables.
!
!--------------------------------------------------------------------------

!     Local variables

INTEGER    :: col, row, pos
REAL (dp)  :: eps1, ten = 10.0, total, work(ncol)

!     EPS is a machine-dependent constant.

IF (PRESENT(eps)) THEN
  eps1 = MAX(ABS(eps), ten * EPSILON(ten))
ELSE
  eps1 = ten * EPSILON(ten)
END IF

!     Set tol(i) = sum of absolute values in column I of R after
!     scaling each element by the square root of its row multiplier,
!     multiplied by EPS1.

work = SQRT(ABS(d))
DO col = 1, ncol
  pos = col - 1
  total = work(col)
  DO row = 1, col-1
    total = total + ABS(r(pos)) * work(row)
    pos = pos + ncol - row - 1
  END DO
  tol(col) = eps1 * total
END DO

tol_set = .TRUE.
RETURN
END SUBROUTINE tolset




SUBROUTINE sing(lindep, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Checks for singularities, reports, and adjusts orthogonal
!     reductions produced by AS75.1.

!     Correction - 19 August 2002
!     When negative weights are used, it is possible for an alement of D
!     to be negative.

!     Auxiliary routines called: INCLUD, TOLSET
!
!--------------------------------------------------------------------------

INTEGER, INTENT(OUT)                :: ifault
LOGICAL, DIMENSION(:), INTENT(OUT)  :: lindep

!     Local variables

REAL (dp)  :: temp, x(ncol), work(ncol), y, weight
INTEGER    :: pos, row, pos2

ifault = 0

work = SQRT(ABS(d))
IF (.NOT. tol_set) CALL tolset()

DO row = 1, ncol
  temp = tol(row)
  pos = row_ptr(row)         ! pos = location of first element in row

!     If diagonal element is near zero, set it to zero, set appropriate
!     element of LINDEP, and use INCLUD to augment the projections in
!     the lower rows of the orthogonalization.

  lindep(row) = .FALSE.
  IF (work(row) <= temp) THEN
    lindep(row) = .TRUE.
    ifault = ifault - 1
    IF (row < ncol) THEN
      pos2 = pos + ncol - row - 1
      x = zero
      x(row+1:ncol) = r(pos:pos2)
      y = rhs(row)
      weight = d(row)
      r(pos:pos2) = zero
      d(row) = zero
      rhs(row) = zero
      CALL includ(weight, x, y)
                             ! INCLUD automatically increases the number
                             ! of cases each time it is called.
      nobs = nobs - 1
    ELSE
      sserr = sserr + d(row) * rhs(row)**2
    END IF ! (row < ncol)
  END IF ! (work(row) <= temp)
END DO ! row = 1, ncol

RETURN
END SUBROUTINE sing



SUBROUTINE ss()

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculates partial residual sums of squares from an orthogonal
!     reduction from AS75.1.
!
!--------------------------------------------------------------------------

!     Local variables

INTEGER    :: i
REAL (dp)  :: total

total = sserr
rss(ncol) = sserr
DO i = ncol, 2, -1
  total = total + d(i) * rhs(i)**2
  rss(i-1) = total
END DO

rss_set = .TRUE.
RETURN
END SUBROUTINE ss



SUBROUTINE cov(nreq, var, covmat, dimcov, sterr, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculate covariance matrix for regression coefficients for the
!     first nreq variables, from an orthogonal reduction produced from
!     AS75.1.

!     Auxiliary routine called: INV
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                   :: nreq, dimcov
INTEGER, INTENT(OUT)                  :: ifault
REAL (dp), INTENT(OUT)                :: var
REAL (dp), DIMENSION(:), INTENT(OUT)  :: covmat, sterr

!     Local variables.

INTEGER                :: dim_rinv, pos, row, start, pos2, col, pos1, k
REAL (dp)              :: total
REAL (dp), ALLOCATABLE :: rinv(:)

!     Check that dimension of array covmat is adequate.

IF (dimcov < nreq*(nreq+1)/2) THEN
  ifault = 1
  RETURN
END IF

!     Check for small or zero multipliers on the diagonal.

ifault = 0
DO row = 1, nreq
  IF (ABS(d(row)) < vsmall) ifault = -row
END DO
IF (ifault /= 0) RETURN

!     Calculate estimate of the residual variance.

IF (nobs > nreq) THEN
  IF (.NOT. rss_set) CALL ss()
  var = rss(nreq) / (nobs - nreq)
ELSE
  ifault = 2
  RETURN
END IF

dim_rinv = nreq*(nreq-1)/2
ALLOCATE ( rinv(dim_rinv) )

CALL INV(nreq, rinv)
pos = 1
start = 1
DO row = 1, nreq
  pos2 = start
  DO col = row, nreq
    pos1 = start + col - row
    IF (row == col) THEN
      total = one / d(col)
    ELSE
      total = rinv(pos1-1) / d(col)
    END IF
    DO K = col+1, nreq
      total = total + rinv(pos1) * rinv(pos2) / d(k)
      pos1 = pos1 + 1
      pos2 = pos2 + 1
    END DO ! K = col+1, nreq
    covmat(pos) = total * var
    IF (row == col) sterr(row) = SQRT(covmat(pos))
    pos = pos + 1
  END DO ! col = row, nreq
  start = start + nreq - row
END DO ! row = 1, nreq

DEALLOCATE(rinv)
RETURN
END SUBROUTINE cov



SUBROUTINE inv(nreq, rinv)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Invert first nreq rows and columns of Cholesky factorization
!     produced by AS 75.1.
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: nreq
REAL (dp), DIMENSION(:), INTENT(OUT) :: rinv

!     Local variables.

INTEGER    :: pos, row, col, start, k, pos1, pos2
REAL (dp)  :: total

!     Invert R ignoring row multipliers, from the bottom up.

pos = nreq * (nreq-1)/2
DO row = nreq-1, 1, -1
  start = row_ptr(row)
  DO col = nreq, row+1, -1
    pos1 = start
    pos2 = pos
    total = zero
    DO k = row+1, col-1
      pos2 = pos2 + nreq - k
      total = total - r(pos1) * rinv(pos2)
      pos1 = pos1 + 1
    END DO ! k = row+1, col-1
    rinv(pos) = total - r(pos1)
    pos = pos - 1
  END DO ! col = nreq, row+1, -1
END DO ! row = nreq-1, 1, -1

RETURN
END SUBROUTINE inv



SUBROUTINE partial_corr(in, cormat, dimc, ycorr, ifault)

!     Replaces subroutines PCORR and COR of:
!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculate partial correlations after the variables in rows
!     1, 2, ..., IN have been forced into the regression.
!     If IN = 1, and the first row of R represents a constant in the
!     model, then the usual simple correlations are returned.

!     If IN = 0, the value returned in array CORMAT for the correlation
!     of variables Xi & Xj is:
!       sum ( Xi.Xj ) / Sqrt ( sum (Xi^2) . sum (Xj^2) )

!     On return, array CORMAT contains the upper triangle of the matrix of
!     partial correlations stored by rows, excluding the 1's on the diagonal.
!     e.g. if IN = 2, the consecutive elements returned are:
!     (3,4) (3,5) ... (3,ncol), (4,5) (4,6) ... (4,ncol), etc.
!     Array YCORR stores the partial correlations with the Y-variable
!     starting with YCORR(IN+1) = partial correlation with the variable in
!     position (IN+1).
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: in, dimc
INTEGER, INTENT(OUT)                 :: ifault
REAL (dp), DIMENSION(:), INTENT(OUT) :: cormat, ycorr

!     Local variables.

INTEGER    :: base_pos, pos, row, col, col1, col2, pos1, pos2
REAL (dp)  :: rms(in+1:ncol), sumxx, sumxy, sumyy, work(in+1:ncol)

!     Some checks.

ifault = 0
IF (in < 0 .OR. in > ncol-1) ifault = ifault + 4
IF (dimc < (ncol-in)*(ncol-in-1)/2) ifault = ifault + 8
IF (ifault /= 0) RETURN

!     Base position for calculating positions of elements in row (IN+1) of R.

base_pos = in*ncol - (in+1)*(in+2)/2

!     Calculate 1/RMS of elements in columns from IN to (ncol-1).

IF (d(in+1) > zero) rms(in+1) = one / SQRT(d(in+1))
DO col = in+2, ncol
  pos = base_pos + col
  sumxx = d(col)
  DO row = in+1, col-1
    sumxx = sumxx + d(row) * r(pos)**2
    pos = pos + ncol - row - 1
  END DO ! row = in+1, col-1
  IF (sumxx > zero) THEN
    rms(col) = one / SQRT(sumxx)
  ELSE
    rms(col) = zero
    ifault = -col
  END IF ! (sumxx > zero)
END DO ! col = in+1, ncol-1

!     Calculate 1/RMS for the Y-variable

sumyy = sserr
DO row = in+1, ncol
  sumyy = sumyy + d(row) * rhs(row)**2
END DO ! row = in+1, ncol
IF (sumyy > zero) sumyy = one / SQRT(sumyy)

!     Calculate sums of cross-products.
!     These are obtained by taking dot products of pairs of columns of R,
!     but with the product for each row multiplied by the row multiplier
!     in array D.

pos = 1
DO col1 = in+1, ncol
  sumxy = zero
  work(col1+1:ncol) = zero
  pos1 = base_pos + col1
  DO row = in+1, col1-1
    pos2 = pos1 + 1
    DO col2 = col1+1, ncol
      work(col2) = work(col2) + d(row) * r(pos1) * r(pos2)
      pos2 = pos2 + 1
    END DO ! col2 = col1+1, ncol
    sumxy = sumxy + d(row) * r(pos1) * rhs(row)
    pos1 = pos1 + ncol - row - 1
  END DO ! row = in+1, col1-1

!     Row COL1 has an implicit 1 as its first element (in column COL1)

  pos2 = pos1 + 1
  DO col2 = col1+1, ncol
    work(col2) = work(col2) + d(col1) * r(pos2)
    pos2 = pos2 + 1
    cormat(pos) = work(col2) * rms(col1) * rms(col2)
    pos = pos + 1
  END DO ! col2 = col1+1, ncol
  sumxy = sumxy + d(col1) * rhs(col1)
  ycorr(col1) = sumxy * rms(col1) * sumyy
END DO ! col1 = in+1, ncol-1

ycorr(1:in) = zero

RETURN
END SUBROUTINE partial_corr




SUBROUTINE vmove(from, to, ifault)

!     ALGORITHM AS274 APPL. STATIST. (1992) VOL.41, NO. 2

!     Move variable from position FROM to position TO in an
!     orthogonal reduction produced by AS75.1.
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)    :: from, to
INTEGER, INTENT(OUT)   :: ifault

!     Local variables

REAL (dp)  :: d1, d2, x, d1new, d2new, cbar, sbar, y, tmplb
INTEGER    :: m, first, last, inc, m1, m2, mp1, col, pos, row

!     Check input parameters

ifault = 0
IF (from < 1 .OR. from > ncol) ifault = ifault + 4
IF (to < 1 .OR. to > ncol) ifault = ifault + 8
IF (ifault /= 0) RETURN

IF (from == to) RETURN

IF (.NOT. rss_set) CALL ss()

IF (from < to) THEN
  first = from
  last = to - 1
  inc = 1
ELSE
  first = from - 1
  last = to
  inc = -1
END IF

DO m = first, last, inc

!     Find addresses of first elements of R in rows M and (M+1).

  m1 = row_ptr(m)
  m2 = row_ptr(m+1)
  mp1 = m + 1
  d1 = d(m)
  d2 = d(mp1)

!     Special cases.

  IF (d1 < vsmall .AND. d2 < vsmall) GO TO 40
  x = r(m1)
  IF (ABS(x) * SQRT(d1) < tol(mp1)) THEN
    x = zero
  END IF
  IF (d1 < vsmall .OR. ABS(x) < vsmall) THEN
    d(m) = d2
    d(mp1) = d1
    r(m1) = zero
    DO col = m+2, ncol
      m1 = m1 + 1
      x = r(m1)
      r(m1) = r(m2)
      r(m2) = x
      m2 = m2 + 1
    END DO ! col = m+2, ncol
    x = rhs(m)
    rhs(m) = rhs(mp1)
    rhs(mp1) = x
    GO TO 40
  ELSE IF (d2 < vsmall) THEN
    d(m) = d1 * x**2
    r(m1) = one / x
    r(m1+1:m1+ncol-m-1) = r(m1+1:m1+ncol-m-1) / x
    rhs(m) = rhs(m) / x
    GO TO 40
  END IF

!     Planar rotation in regular case.

  d1new = d2 + d1*x**2
  cbar = d2 / d1new
  sbar = x * d1 / d1new
  d2new = d1 * cbar
  d(m) = d1new
  d(mp1) = d2new
  r(m1) = sbar
  DO col = m+2, ncol
    m1 = m1 + 1
    y = r(m1)
    r(m1) = cbar*r(m2) + sbar*y
    r(m2) = y - x*r(m2)
    m2 = m2 + 1
  END DO ! col = m+2, ncol
  y = rhs(m)
  rhs(m) = cbar*rhs(mp1) + sbar*y
  rhs(mp1) = y - x*rhs(mp1)

!     Swap columns M and (M+1) down to row (M-1).

  40 pos = m
  DO row = 1, m-1
    x = r(pos)
    r(pos) = r(pos-1)
    r(pos-1) = x
    pos = pos + ncol - row - 1
  END DO ! row = 1, m-1

!     Adjust variable order (VORDER), the tolerances (TOL) and
!     the vector of residual sums of squares (RSS).
! code for moving constraints indicator vector
    m1 = corder(m)
    corder(m) = corder(mp1)
    corder(mp1) = m1
! code for moving constraints indicator vector

! code for moving constraints vector
    tmplb = wtslb(m)
    wtslb(m) = wtslb(mp1)
    wtslb(mp1) = tmplb
! code for moving constraints vector
  m1 = vorder(m)
  vorder(m) = vorder(mp1)
  vorder(mp1) = m1
  x = tol(m)
  tol(m) = tol(mp1)
  tol(mp1) = x
  rss(m) = rss(mp1) + d(mp1) * rhs(mp1)**2
END DO

RETURN
END SUBROUTINE vmove



SUBROUTINE reordr(list, n, pos1, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Re-order the variables in an orthogonal reduction produced by
!     AS75.1 so that the N variables in LIST start at position POS1,
!     though will not necessarily be in the same order as in LIST.
!     Any variables in VORDER before position POS1 are not moved.

!     Auxiliary routine called: VMOVE
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)               :: n, pos1
INTEGER, DIMENSION(:), INTENT(IN) :: list
INTEGER, INTENT(OUT)              :: ifault

!     Local variables.

INTEGER    :: next, i, l, j

!     Check N.

ifault = 0
IF (n < 1 .OR. n > ncol+1-pos1) ifault = ifault + 4
IF (ifault /= 0) RETURN

!     Work through VORDER finding variables which are in LIST.

next = pos1
i = pos1
10 l = vorder(i)
DO j = 1, n
  IF (l == list(j)) GO TO 40
END DO
30 i = i + 1
IF (i <= ncol) GO TO 10

!     If this point is reached, one or more variables in LIST has not
!     been found.

ifault = 8
RETURN

!     Variable L is in LIST; move it up to position NEXT if it is not
!     already there.

40 IF (i > next) CALL vmove(i, next, ifault)
next = next + 1
IF (next < n+pos1) GO TO 30

RETURN
END SUBROUTINE reordr



SUBROUTINE hdiag(xrow, nreq, hii, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
!
!                         -1           -1
! The hat matrix H = x(X'X) x' = x(R'DR) x' = z'Dz

!              -1
! where z = x'R

! Here we only calculate the diagonal element hii corresponding to one
! row (xrow).   The variance of the i-th least-squares residual is (1 - hii).
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: nreq
INTEGER, INTENT(OUT)                 :: ifault
REAL (dp), DIMENSION(:), INTENT(IN)  :: xrow
REAL (dp), INTENT(OUT)               :: hii

!     Local variables

INTEGER    :: col, row, pos
REAL (dp)  :: total, wk(ncol)

!     Some checks

ifault = 0
IF (nreq > ncol) ifault = ifault + 4
IF (ifault /= 0) RETURN

!     The elements of xrow.inv(R).sqrt(D) are calculated and stored in WK.

hii = zero
DO col = 1, nreq
  IF (SQRT(d(col)) <= tol(col)) THEN
    wk(col) = zero
  ELSE
    pos = col - 1
    total = xrow(col)
    DO row = 1, col-1
      total = total - wk(row)*r(pos)
      pos = pos + ncol - row - 1
    END DO ! row = 1, col-1
    wk(col) = total
    hii = hii + total**2 / d(col)
  END IF
END DO ! col = 1, nreq

RETURN
END SUBROUTINE hdiag



FUNCTION varprd(x, nreq) RESULT(fn_val)

!     Calculate the variance of x'b where b consists of the first nreq
!     least-squares regression coefficients.
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: nreq
REAL (dp), DIMENSION(:), INTENT(IN)  :: x
REAL (dp)                            :: fn_val

!     Local variables

INTEGER    :: ifault, row
REAL (dp)  :: var, wk(nreq)

!     Check input parameter values

fn_val = zero
ifault = 0
IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
IF (nobs <= nreq) ifault = ifault + 8
IF (ifault /= 0) THEN
  ! WRITE(*, '(1x, a, i4)') 'Error in function VARPRD: ifault =', ifault
  RETURN
END IF

!     Calculate the residual variance estimate.

var = sserr / (nobs - nreq)

!     Variance of x'b = var.x'(inv R)(inv D)(inv R')x
!     First call BKSUB2 to calculate (inv R')x by back-substitution.

CALL BKSUB2(x, wk, nreq)
DO row = 1, nreq
  IF(d(row) > tol(row)) fn_val = fn_val + wk(row)**2 / d(row)
END DO

fn_val = fn_val * var

RETURN
END FUNCTION varprd



SUBROUTINE bksub2(x, b, nreq)

!     Solve x = R'b for b given x, using only the first nreq rows and
!     columns of R, and only the first nreq elements of R.
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: nreq
REAL (dp), DIMENSION(:), INTENT(IN)  :: x
REAL (dp), DIMENSION(:), INTENT(OUT) :: b

!     Local variables

INTEGER    :: pos, row, col
REAL (dp)  :: temp

!     Solve by back-substitution, starting from the top.

DO row = 1, nreq
  pos = row - 1
  temp = x(row)
  DO col = 1, row-1
    temp = temp - r(pos)*b(col)
    pos = pos + ncol - col - 1
  END DO
  b(row) = temp
END DO

RETURN
END SUBROUTINE bksub2


SUBROUTINE gmove(first, last, ifault)
! Revised from vmove.

!     Move group of variables from position FROM to position TO in an
!     orthogonal reduction produced by AS75.1.
!
! It calls vmove to acutally move the variables.
! Yi Guo@CMIS, 8/1/2010
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)    :: first, last
INTEGER, INTENT(OUT)   :: ifault

!     Local variables
INTEGER    :: m, nvg, new_first, new_last, g_ind

!     Check input parameters

ifault = 0
IF (first < 1 .OR. first > ncol) ifault = ifault + 4
IF (last < 1 .OR. last > ncol) ifault = ifault + 8
IF (ifault /= 0) RETURN

IF (first == last) RETURN

IF (.NOT. rss_set) CALL ss()

IF (first > last) THEN
  ifault = -999
  return
END IF
! ng: number of groups
! gnv: contains the number of variables in indexed group
! g_ptr: the start position of variables of indexed group.
!       Notice this indices are always in descent order. But the groups indices will change after move. For example, the sequence of variables in groups
!       is 1122233334455566 and g_ptr is 1 3 6 10 12 15 17(one more for dummy variable). After moving group 2 to group 5, it becomes 1133334455522266.So
!       g_ptr is 1 3 7 9 12 15 17.The group in second position is group 3. g_ptr(2) indicates the start position of the group at second position, not
!       the start position of group 2.
! gorder: the order of groups.

g_ind = gorder(first)
nvg = gnv(g_ind)
new_first = g_ptr(first)
new_last = g_ptr(last+1)
DO m =1, nvg
    new_first = g_ptr(first)
    new_last = new_last - 1
    call vmove(new_first,new_last,ifault)
END DO

! RSS has been updated in vmove, so don't have to worry here.
! Update gorder, g_ptr

gorder(first:(last-1)) = gorder((first+1):last)
g_ptr(first:(last-1)) = g_ptr((first+1):last) - nvg
g_ptr(last) = g_ptr(last-1) + gnv(gorder((last-1)))
gorder(last) = g_ind
RETURN
END SUBROUTINE gmove

SUBROUTINE qrgupdate(weight, xrow, yelem, n, wd, wr, wrhs, wsserr)
! Original comments from includ
!     ALGORITHM AS75.1  APPL. STATIST. (1974) VOL.23, NO. 3

!     Calling this routine updates D, R, RHS and SSERR by the
!     inclusion of xrow, yelem with the specified weight.

!     *** WARNING  Array XROW is overwritten.

!     N.B. As this routine will be called many times in most applications,
!          checks have been eliminated.
!
!--------------------------------------------------------------------------
!
! Yi's comments from here:
! This program is to get RSS reduction for a group of variables using QR
! update. It works for gadd1 to add a group in the inner most loop, which
! is designed to fully utilize the current QR decomposition. It is revised
! from includ which updates the public variables: r, rhs, d and does it row
! by row. qrupdate updates only local variables but the updating is the same.
! Yi Guo @CMIS, 03/01/2012

IMPLICIT NONE
REAL (dp),INTENT(IN)                    :: weight, yelem
integer, intent (in) :: n !the number of variables
REAL (dp), DIMENSION(:), INTENT(IN OUT) :: xrow, wd, wr, wrhs
real (dp), intent(in out) :: wsserr

!     Local variables

INTEGER     :: i, k, nextr
REAL (dp)   :: w, y, xi, di, wxi, dpi, cbar, sbar, xk

w = weight
y = yelem
rss_set = .false.
nextr = 1
DO i = 1, n

!     Skip unnecessary transformations.   Test on exact zeroes must be
!     used or stability can be destroyed.

  IF (ABS(w) < vsmall) RETURN
  xi = xrow(i)
  IF (ABS(xi) < vsmall) THEN
    nextr = nextr + n - i
  ELSE
    di = wd(i)
    wxi = w * xi
    dpi = di + wxi*xi
    cbar = di / dpi
    sbar = wxi / dpi
    w = cbar * w
    wd(i) = dpi
    DO k = i+1, n
      xk = xrow(k)
      xrow(k) = xk - xi * wr(nextr)
      wr(nextr) = cbar * wr(nextr) + sbar * xk
      nextr = nextr + 1
    END DO
    xk = y
    y = xk - xi * wrhs(i)
    wrhs(i) = cbar * wrhs(i) + sbar * xk
  END IF
END DO ! i = 1, ncol

!     Y * SQRT(W) is now equal to the Brown, Durbin & Evans recursive
!     residual.

wsserr = wsserr + w * y * y

RETURN
END SUBROUTINE qrgupdate

SUBROUTINE gaddposweights(beta, nreq, poswtsindicator, wr, wrow_ptr, wrhs, ifault)
! Original comments from regcf subroutine.
!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Modified version of AS75.4 to calculate regression coefficients
!     for the first NREQ variables, given an orthogonal reduction from
!     AS75.1.
!
!--------------------------------------------------------------------------
!
! Yi's comments
! Revised from regcf. It tests whether all the non-forced variables have positive weights.
! It will be used in gadd1 subroutine. If the weights are not positive, it is meaningless
! to report it at all. It will increase the computational time for sure, but it will save
! time in following operations.
!
! Effected subroutines: report (if this function is in action)
! New constraints vector is introduced to adapt to different positive weights constraints pattern.
! Yi Guo@CMIS, 05/01/2012
!

IMPLICIT NONE
INTEGER, INTENT(IN)                  :: nreq, poswtsindicator(:)
INTEGER, INTENT(OUT)                 :: ifault
REAL (dp), DIMENSION(nreq), intent (out) :: beta
real (dp), intent (in) :: wr(:), wrhs(:)
integer, intent (in) :: wrow_ptr(:)

!     Local variables

INTEGER   :: i, j, nextr

!     Some checks.

ifault = 0
IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
IF (ifault /= 0) RETURN

!IF (.NOT. tol_set) CALL tolset()

DO i = nreq, 1, -1
!  IF (SQRT(d(i)) < tol(i)) THEN
!    beta(i) = zero
!    d(i) = zero
!    ifault = -i
!  ELSE
    beta(i) = wrhs(i)
    nextr = wrow_ptr(i)
    DO j = i+1, nreq
      beta(i) = beta(i) - wr(nextr) * beta(j)
      nextr = nextr + 1
    END DO ! j = i+1, nreq
!  END IF
  if (poswtsindicator(i)==1) then
      if (beta(i) < 0) then
        ifault = -10 * nreq
        return
      end if
  end if
END DO ! i = nreq, 1, -1

RETURN
END SUBROUTINE gaddposweights

SUBROUTINE validweights(beta, nreq, ifault)
! Original comments from regcf subroutine.
!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Modified version of AS75.4 to calculate regression coefficients
!     for the first NREQ variables, given an orthogonal reduction from
!     AS75.1.
!
!--------------------------------------------------------------------------
!
! Yi's comments
! Revised from regcf. It test whether all the non-forced variables have valid weights.
! It will be used in report subroutine. If the weights are not valid, it is meaningless
! to report it at all. It will increase the computational time for sure, but it will save
! the time in following operations.
!
! Effected subroutines: report (if this function is put in action)
! Yi Guo@CMIS, 06/01/2012
!

IMPLICIT NONE
INTEGER, INTENT(IN)                  :: nreq
INTEGER, INTENT(OUT)                 :: ifault
REAL (dp), DIMENSION(nreq), intent (out) :: beta

!     Local variables

INTEGER   :: i, j, nextr

!     Some checks.

ifault = 0
IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
IF (ifault /= 0) RETURN

IF (.NOT. tol_set) CALL tolset()

DO i = nreq, 1, -1
  IF (SQRT(d(i)) < tol(i)) THEN
    beta(i) = zero
    d(i) = zero
    ifault = -i
  ELSE
    beta(i) = rhs(i)
    nextr = row_ptr(i)
    DO j = i+1, nreq
      beta(i) = beta(i) - r(nextr) * beta(j)
      nextr = nextr + 1
    END DO ! j = i+1, nreq
  END IF
  if (corder(i) ==1) then
      !if (beta(i) < 0) then
      if (beta(i) < wtslb(i)) then ! Change it to weights constraint test.
        ifault = -10 * nreq
        return
      end if
  end if
END DO ! i = nreq, 1, -1

RETURN
END SUBROUTINE validweights

SUBROUTINE gaddwtsconstraint(beta, nreq, wtsconindicator, wtslowbounds,wr, wrow_ptr, wrhs, ifault)
! Revised from gaddposweights. It tests whether all the non-forced variables have weights
! violation, i.e. wts(i) < epsilon(i) where epsilon is a vector of real numbers showing the
! lower bound of that weight. If weight is not good, then it will return a negative number.
! Otherwise return is zero which means all weights are fine.
!

IMPLICIT NONE
INTEGER, INTENT(IN)                  :: nreq, wtsconindicator(:)
REAL (dp), INTENT(IN)                :: wtslowbounds(:)
INTEGER, INTENT(OUT)                 :: ifault
REAL (dp), DIMENSION(nreq), intent (out) :: beta
real (dp), intent (in) :: wr(:), wrhs(:)
integer, intent (in) :: wrow_ptr(:)

!     Local variables

INTEGER   :: i, j, nextr

!     Some checks.

ifault = 0
IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
IF (ifault /= 0) RETURN

!IF (.NOT. tol_set) CALL tolset()

DO i = nreq, 1, -1
!  IF (SQRT(d(i)) < tol(i)) THEN
!    beta(i) = zero
!    d(i) = zero
!    ifault = -i
!  ELSE
    beta(i) = wrhs(i)
    nextr = wrow_ptr(i)
    DO j = i+1, nreq
      beta(i) = beta(i) - wr(nextr) * beta(j)
      nextr = nextr + 1
    END DO ! j = i+1, nreq
!  END IF
  if (wtsconindicator(i)==1) then
      if (beta(i) < wtslowbounds(i)) then
        ifault = -10 * nreq
        return
      end if
!      write(dbfid,*) 'nreq:',nreq
!      write(dbfid,*) 'beta(',i,')',beta(i)
  end if
END DO ! i = nreq, 1, -1

RETURN
END SUBROUTINE gaddwtsconstraint

END MODULE lsq
MODULE find_subsets

! A module of routines for finding and recording best-fitting subsets of
! regression variables

!   Version 1.10, 17 February 2004
!   Author: Alan Miller
!           formerly of CSIRO Division of Mathematical & Information Sciences
!   Phone: (+61) 3 9592-5085
!   e-mail: amiller @ bigpond.net.au
!   WWW-page: http://users.bigpond.net.au/amiller/

! 17 Feb 2004 Correction to subroutine EFROYM for the case in which all the
!             variables are selected. Thanks to David Jones.
! 12 Nov 1999 Made changes to routines exadd1 & add2 to prevent the calculation
!             of negative residual sums of squares which could occur in cases
!             in which the true RSS is zero.   Routine seq2 changed to avoid
!             cycling.
! 24 May 2000 Changed lsq_kind to dp (cosmetic change only)
! 4 June 2000 Added routine random_pick which picks a random set of variables.
! 29 Aug 2002 Set value of size in subroutine EFROYM when IER /= 0.

USE lsq

IMPLICIT NONE

INTEGER, SAVE                 :: max_size, nbest, lopt_dim1
REAL (dp), ALLOCATABLE, SAVE  :: bound(:), ress(:,:), wts(:,:)  !wts is for weights of valid groups variables
INTEGER, ALLOCATABLE, SAVE    :: lopt(:,:), wts_ptr(:), vwts(:,:), rank_ind(:,:), nvwts(:,:)!Locators of weights and variables, weights code

! Add wts to store regression coefficients for valid group variables. The following describes the structure of it.
! wts is a 2-D array with nbest columns and (nvar_max+1)*nvar_max/2a*M + nvar_max*nforce rows where M is the number of variables in largest group.
! wts will be initialized to be zero. It should be treated in the same way as ress, lopt etc. in terms of memory management.
! Effected subroutines (functions): init_subsets, report, validweights, modual find_sub

CONTAINS

SUBROUTINE init_subsets(nvar_max, fit_const, nvar)

INTEGER, INTENT(IN)           :: nvar_max
LOGICAL, INTENT(IN)           :: fit_const
INTEGER, INTENT(IN), OPTIONAL :: nvar

!     Local variables

INTEGER    :: i, ier, m, idx !Yi weights code
REAL (dp)  :: eps = 1.E-14
LOGICAL    :: lindep(ncol)

!     The LSQ module has probably already been initialized, but just in case ..

IF (.NOT. initialized) THEN
  IF (PRESENT(nvar)) CALL startup(nvar, fit_const)
END IF

IF (fit_const) THEN
  max_size = nvar_max + 1
ELSE
  max_size = nvar_max
END IF

lopt_dim1 = max_size * (max_size + 1) / 2
IF (ALLOCATED(bound)) DEALLOCATE(bound, ress, lopt, wts, wts_ptr,vwts, rank_ind, nvwts) ! Yi, weights code
ALLOCATE (bound(max_size), ress(max_size,nbest), lopt(lopt_dim1,nbest))

! Yi, weights code
call maximum(gnv,ng, m, idx)
allocate (wts(lopt_dim1*m+nforce*max_size,nbest),wts_ptr(max_size+1), &
           vwts(lopt_dim1*m+nforce*max_size,nbest),rank_ind(max_size,nbest),&
           nvwts(max_size,nbest))
rank_ind = 0 ! New code for init rank_ind = 0
wts_ptr = 1
nvwts = 0
wts = 0
vwts = 0
do i = 2, max_size + 1
    wts_ptr(i) = wts_ptr(i-1) + (i-1) * m + nforce
end do
! Yi, weights code

bound = HUGE(eps)
ress  = HUGE(eps)
lopt  = 0

CALL tolset(eps)
CALL sing(lindep, ier)

CALL ss()
DO i = 1, max_size
  !CALL report(i, rss(i))   !Yi
  CALL report(i, rss(g_ptr(i+1)-1)) !Yi
END DO

RETURN
END SUBROUTINE init_subsets


SUBROUTINE add1(first, last, ss, smax, jmax, ier)

! Calculate the reduction in residual sum of squares when one variable,
! selected from those in positions FIRST .. LAST, is added in position FIRST,
! given that the variables in positions 1 .. FIRST-1 (if any) are already
! included.

INTEGER, INTENT(IN)     :: first, last
INTEGER, INTENT(OUT)    :: jmax, ier
REAL (dp), INTENT(OUT)  :: ss(:), smax

!     Local variables

INTEGER    :: j, inc, pos, row, col,ok
REAL (dp)  :: zero = 0.0_dp, diag, dy, ssqx, sxx(ncol), sxy(ncol)

!  AAG Code to skip combinations we don't want
INTEGER, ALLOCATABLE :: ulist(:)
ALLOCATE (ulist(first))
ulist(1:first-1)=vorder(1:first-1)


!     Check call arguments

jmax = 0
smax = zero
ier = 0
IF (first > ncol) ier = 1
IF (last < first) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

!     Accumulate sums of squares & products from row FIRST

sxx(first:last) = zero
sxy(first:last) = zero
inc = ncol - last
pos = row_ptr(first)
DO row = first, last
  diag = d(row)
  dy = diag * rhs(row)
  sxx(row) = sxx(row) + diag
  sxy(row) = sxy(row) + dy
  DO col = row+1, last
    sxx(col) = sxx(col) + diag * r(pos)**2
    sxy(col) = sxy(col) + dy * r(pos)
    pos = pos + 1
  END DO
  pos = pos + inc
END DO

!     Incremental sum of squares for a variable = SXY * SXY / SXX.
!     Calculate whenever sqrt(SXX) > TOL for that variable.

DO j = first, last
  ssqx = sxx(j)
  ulist(first)=j										!   AAG
  call legal(ulist,ok,1)								!   AAG
  !IF (SQRT(ssqx) > tol(j) .and. ok == 1) THEN			!   AAG
  IF (SQRT(ssqx) > tol(j)) THEN
    !write(*,*),ok,'---',ulist
    ss(j) = sxy(j)**2 / sxx(j)
    IF (ss(j) > smax) THEN
      smax = ss(j)
      jmax = j
    END IF
  ELSE
    ss(j) = zero
  END IF
END DO

RETURN
END SUBROUTINE add1


SUBROUTINE add2(first, last, smax, j1, ier)
! Origianl comments
!     Calculate the maximum reduction in residual sum of squares when 2
!     variables, selected from those in positions FIRST .. LAST, are
!     added, given that the variables in positions 1 .. FIRST-1 (if
!     any) are already included.    J1, J2 are the positions of the two
!     best variables.   N.B. J2 < J1.
!
! Comments from Yi
! To incorporate adding 2 variables a time, it is a similar implementation of exadd1. It updates the NBEST subsets of ivar variables found
! from a call to subroutine add2.
! Yi Guo@CMIS, 12/1/2010
!
! All variables are in group (suppose second PC is included). So if there are 60 classes, it is still indexed by classes, use
! 2 * (i-1) + 1/2 to index the class and its second PC (i means i-th class). This subroutine will then test RSS reduction from
! first class to last class with 2nd PC. It is different from original add2 where arbitary mix of 2 variables are considered.
! It just tests from first to last "bicolumnly". Therefore it is simper than original version. In this case, j2 is no longer
! needed, so just j1 stays.  To make it consistent with the old version, keep j2 here which is j1 + 1. Other arguments are the same.
! N.B. the variables is arranged in this way: x1 p2x1 x2 p2x2 x3 p2x3 .... xn p2xn where p2- means the second PC of xi.

! Remove j2 from the scene to keep cran happy. j2 is simply j1+1 according to quotes above.

INTEGER, INTENT(IN)     :: first, last
INTEGER, INTENT(OUT)    :: j1, ier ! removed j2 here.
REAL (dp), INTENT(OUT)  :: smax
!     Local variables

INTEGER    :: i, row, ltemp, m, j2
REAL (dp)  :: zero = 0.0_dp, temp, det, two = 2.0, xx1, xx2, x1x2, yx1, yx2
real (dp) :: ss(nbest) !contains the nbest RSS reduction in order to update boundaries and nbest mixtures it finds so far
integer :: idx(nbest)  !the corresponding indices of variables of ss
real (dp), allocatable :: beta(:),wrhs(:), wtsvars(:), subwts(:,:) ! weights code
integer :: fposwts, wtsdim, maxval, maxidx !weights code
integer, allocatable :: wrow_ptr(:), subvars(:,:), subnvars(:), subrank(:), subloc(:) !Locators of weights and variables, weights code !weights
real (dp) :: rtemp

!     Check call arguments

smax = zero
j1 = 0
ier = 0
IF (first > ncol) ier = 1
IF (last <= first) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN
ss = -1
idx = 0

! Yi, weights code
wtsdim = first*2+nforce
allocate (subwts(wtsdim,nbest), subvars(wtsdim,nbest),&
           subrank(nbest), subnvars(nbest),subloc(nbest))
subrank = 1
subnvars = 0
! Yi, weights code

!! The following 3 lines are for debugging. Output data and read those in R to compare results.
!OPEN(10, FILE='data.txt')
!write(10,*),nforce,ncol,r,d,rhs,ng,g_ptr,first,last,vorder,gorder,row_ptr
!close(10)
!! Debugging code ends here.

! Start my own one from scratch
do i = g_ptr(first), g_ptr(last), 2
    xx1 = d(i)
    x1x2 = d(i)*r(row_ptr(i))
    xx2 = x1x2*r(row_ptr(i)) + d(i+1)   !xx2 = d(i)*r(row_ptr(i))**2 + d(i+1)
    yx1 = d(i)*rhs(i)
    yx2 = d(i)*r(row_ptr(i))*rhs(i) + d(i+1)*rhs(i+1)
    do row = g_ptr(first),i-1
        xx1 = xx1 + d(row)*r(row_ptr(row)+i-row-1)**2
        xx2 = xx2 + d(row)*r(row_ptr(row)+i-row)**2
        x1x2 = x1x2 + d(row)*r(row_ptr(row)+i-row-1)*r(row_ptr(row)+i-row)
        yx1 = yx1 + d(row)*r(row_ptr(row)+i-row-1)*rhs(row)
        yx2 = yx2 + d(row)*r(row_ptr(row)+i-row)*rhs(row)
    end do
!     Calculate reduction in RSS for pair i, i+1.
!     The sum of squares & cross-products are in:
!              ( xx1   x1x2)      ( yx1)
!              ( x1x2  xx2 )      ( yx2)

    det = MAX( (xx1 * xx2 - x1x2**2), zero)
    temp = SQRT(det)
    IF (temp < tol(i)*SQRT(xx2) .OR.temp < tol(i+1)*SQRT(xx1)) CYCLE
    temp = (xx2*yx1**2 - two*x1x2*yx2*yx1 + xx1*yx2**2) / det

    ! weights code
    if (allocated(beta)) deallocate(beta)
    allocate(beta(g_ptr(first)+1))
    beta(g_ptr(first):g_ptr(first)+1) = (/ (xx2*yx1 - x1x2*yx2)/det, (xx1*yx2 - x1x2*yx1)/det /)

    if (beta(g_ptr(first))*corder(i) <0 .or. beta(g_ptr(first)+1)*corder(i+1)<0) cycle ! If the coefficient of 1st PC is negative, discard it.
    if (allocated(wrhs)) deallocate(wrhs)
    allocate (wrhs(g_ptr(first)-1))
    do m = 1, g_ptr(first)-1
        rtemp = 0
        do j2 = 1, 2
            !wrhs(m) = rhs(m) - r(row_ptr(m)+i-1-m) * beta(g_ptr(first)) + r(row_ptr(m)+i-m) * beta(g_ptr(first)+1)
        rtemp = rtemp + r(row_ptr(m)+i-m-2+j2) * beta(g_ptr(first) + j2 - 1)
        end do
        wrhs(m) = rhs(m) - rtemp
    end do

    call gaddposweights(beta, g_ptr(first)-1, corder(1:g_ptr(first)-1), r, row_ptr, wrhs, fposwts) ! Get the weights of groups before first.
    if (fposwts < -2) cycle ! If negative weights found, discard it.

    if (allocated(wtsvars)) deallocate(wtsvars)
    allocate (wtsvars(2+g_ptr(first)-1))
    wtsvars =(/ vorder(1:g_ptr(first)-1), vorder(i:i+1) /) ! Assemble the variables.
    ! weights code

    ! Update ss and idx which are nbest RSS reduction and corresponding indices resp.
    do j2=1,nbest
        if (temp>ss(j2)) then
            if (j2 < nbest) then
                ss(j2+1:nbest) = ss(j2:nbest-1)
                idx(j2+1:nbest) = idx(j2:nbest-1)
                subloc(j2+1:nbest) = subloc(j2:nbest-1)
            end if
            ss(j2) = temp
            idx(j2) = gorder((i-nforce-1)/2+1)
            call maximum(subrank,nbest,maxval,maxidx)
            subwts(1:2+g_ptr(first)-1,maxidx) = beta
            subvars(1:2+g_ptr(first)-1,maxidx) = wtsvars
            where(subrank >= j2) subrank = subrank + 1
            subrank(maxidx) = j2
            subnvars(maxidx) = 2+g_ptr(first)-1 - nforce
            subloc(j2) = maxidx
            exit
        end if
    end do
    IF (temp > smax) THEN
      smax = temp
      j1 = i
    end if
end do
j1 = (j1-1)/2+1
j2 = j1 + 1

! Update the ress, lop, bound etc.
ltemp = gorder(first)
do i= 1, nbest
    ! It is possible that less than nbest combinations are found, for example first = 5, last = 6, nbest =4.
    if (ss(i)>0) then
        gorder(first) = idx(i)
        !call report(first,rss(g_ptr(first)-1)-ss(i)) ! update ! Commented for weights
        call gaddreport(first,rss(g_ptr(first)-1)-ss(i), subwts(:,subloc(i)), subvars(:,subloc(i)), subnvars(subloc(i))) !Weights update
    end if
end do
gorder(first) = ltemp

RETURN
END SUBROUTINE add2


SUBROUTINE bakwrd(first, last, ier)

!     Backward elimination from variables in positions FIRST .. LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.
!     On exit, the array VORDER contains the numbers of the variables
!     in the order in which they were deleted.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: pos, jmin, i
REAL (dp)  :: ss(last), smin

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

!     For POS = LAST, ..., FIRST+1 call DROP1 to find best variable to
!     find which variable to drop next.

DO pos = last, first+1, -1
  CALL drop1(first, pos, ss, smin, jmin, ier)
  CALL exdrop1(first, pos, ss, smin, jmin)
  IF (jmin > 0 .AND. jmin < pos) THEN
    CALL vmove(jmin, pos, ier)
    IF (nbest > 0) THEN
      DO i = jmin, pos-1
        CALL report(i, rss(i))
      END DO
    END IF
  END IF
END DO

RETURN
END SUBROUTINE bakwrd


SUBROUTINE drop1(first, last, ss, smin, jmin, ier)

! Calculate the increase in the residual sum of squares when the variable in
! position J is dropped from the model (i.e. moved to position LAST),
! for J = FIRST, ..., LAST-1.

INTEGER, INTENT(IN)     :: first, last
INTEGER, INTENT(OUT)    :: jmin, ier
REAL (dp), INTENT(OUT)  :: ss(:), smin

!     Local variables

INTEGER    :: j, pos1, inc, pos, row, col, i
REAL (dp)  :: large = HUGE(1.0_dp), zero = 0.0_dp, d1, rhs1, d2, x, wk(last), &
              vsmall = TINY(1.0_dp)

!     Check call arguments

jmin = 0
smin = large
ier = 0
IF (first > ncol) ier = 1
IF (last < first) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

!     POS1 = position of first element of row FIRST in r.

pos1 = row_ptr(first)
inc = ncol - last

!     Start of outer cycle for the variable to be dropped.

DO j = first, last
  d1 = d(j)
  IF (SQRT(d1) < tol(j)) THEN
    ss(j) = zero
    smin = zero
    jmin = j
    GO TO 50
  END IF
  rhs1 = rhs(j)
  IF (j == last) GO TO 40

!     Copy row J of R into WK.

  pos = pos1
  DO i = j+1, last
    wk(i) = r(pos)
    pos = pos + 1
  END DO
  pos = pos + inc

!     Lower the variable past each row.

  DO row = j+1, last
    x = wk(row)
    d2 = d(row)
    IF (ABS(x) * SQRT(d1) < tol(row) .OR. d2 < vsmall) THEN
      pos = pos + ncol - row
      CYCLE
    END IF
    d1 = d1 * d2 / (d2 + d1 * x**2)
    DO col = row+1, last
      wk(col) = wk(col) - x * r(pos)
      pos = pos + 1
    END DO
    rhs1 = rhs1 - x * rhs(row)
    pos = pos + inc
  END DO
  40 ss(j) = rhs1 * d1 * rhs1
  IF (ss(j) < smin) THEN
    jmin = j
    smin = ss(j)
  END IF

!     Update position of first element in row of r.

  50 IF (j < last) pos1 = pos1 + ncol - j
END DO

RETURN
END SUBROUTINE drop1


SUBROUTINE efroym(first, last, fin, fout, size, ier, lout)

!     Efroymson's stepwise regression from variables in positions FIRST,
!     ..., LAST.  If FIRST > 1, variables in positions prior to this are
!     forced in.  If LAST < ncol, variables in positions after this are
!     forced out.

!     A report is written to unit LOUT if LOUT >= 0.

INTEGER, INTENT(IN)    :: first, last, lout
INTEGER, INTENT(OUT)   :: size, ier
REAL (dp), INTENT(IN)  :: fin, fout

!     Local variables

INTEGER    :: jmax, jmin, i
REAL (dp)  :: one = 1.0, eps, zero = 0.0, ss(last), smax, base, var, f, smin

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (fin < fout .OR. fin <= zero) ier = ier + 256
IF (nobs <= ncol) ier = ier + 512
IF (ier /= 0) THEN
  size = 0
  RETURN
END IF

!     EPS approximates the smallest quantity such that the calculated value of
!     (1 + EPS) is > 1.   It is used to test for a perfect fit (RSS = 0).

eps = EPSILON(one)

!     SIZE = number of variables in the current subset

size = first - 1

!     Find the best variable to add next

20 CALL add1(size+1, last, ss, smax, jmax, ier)
IF (nbest > 0) CALL exadd1(size+1, smax, jmax, ss, last)

!     Calculate 'F-to-enter' value

IF (size > 0) THEN
  base = rss(size)
ELSE
  base = rss(1) + ss(1)
END IF
var = (base - smax) / (nobs - size - 1)
IF (var < eps*base) THEN
  ier = -1
  f = zero
ELSE
  f = smax / var
END IF
!IF (lout >= 0) WRITE(lout, 900) vorder(jmax), f
900 FORMAT(' Best variable to add:  ', i4, '  F-to-enter = ', f10.2)

!     Exit if F < FIN or IER < 0 (perfect fit)

IF (f < fin .OR. ier < 0) RETURN

!     Add the variable to the subset (in position FIRST).

!IF (lout >= 0) WRITE(lout, '(50x, "Variable added")')
size = size + 1
IF (jmax > first) CALL vmove(jmax, first, ier)
DO i = first, MIN(jmax-1, max_size)
  CALL report(i, rss(i))
END DO

!     See whether a variable entered earlier can be deleted now.

30 IF (size <= first) GO TO 20
CALL drop1(first+1, size, ss, smin, jmin, ier)
CALL exdrop1(first+1, size, ss, smin, jmin)
var = rss(size) / (nobs - size)
f = smin / var
!IF (lout >= 0) WRITE(lout, 910) vorder(jmin), f
910 FORMAT(' Best variable to drop: ', i4, '  F-to-drop  = ', f10.2)

IF (f < fout) THEN
  !IF (lout >= 0) WRITE(lout, '(50x, "Variable dropped")')
  CALL vmove(jmin, size, ier)
  IF (nbest > 0) THEN
    DO i = jmin, size-1
      CALL report(i, rss(i))
    END DO
  END IF
  size = size - 1
  GO TO 30
END IF

IF (size >= last) RETURN
GO TO 20
END SUBROUTINE efroym


SUBROUTINE exadd1(ivar, smax, jmax, ss, last)

!     Update the NBEST subsets of IVAR variables found from a call
!     to subroutine ADD1.

INTEGER, INTENT(IN)    :: ivar, jmax, last
REAL (dp), INTENT(IN)  :: smax, ss(:)

!     Local variables

REAL (dp)  :: zero = 0.0_dp, ssbase, sm, temp, wk(last)
INTEGER    :: i, j, ltemp, jm

IF (jmax == 0) RETURN
IF (ivar <= 0) RETURN
IF (ivar > max_size) RETURN
ltemp = vorder(ivar)
jm = jmax
sm = smax
IF (ivar > 1) ssbase = rss(ivar-1)
IF (ivar == 1) ssbase = rss(ivar) + ss(1)
wk(ivar:last) = ss(ivar:last)

DO i = 1, nbest
  temp = MAX(ssbase - sm, zero)
  IF (temp >= bound(ivar)) EXIT
  vorder(ivar) = vorder(jm)
  IF (jm == ivar) vorder(ivar) = ltemp
  CALL report(ivar, temp)
  IF (i >= nbest) EXIT
  wk(jm) = zero
  sm = zero
  jm = 0
  DO j = ivar, last
    IF (wk(j) <= sm) CYCLE
    jm = j
    sm = wk(j)
  END DO
  IF (jm == 0) EXIT
END DO

!     Restore VORDER(IVAR)

vorder(ivar) = ltemp

RETURN
END SUBROUTINE exadd1


SUBROUTINE exdrop1(first, last, ss, smin, jmin)
! Record any new subsets of (LAST-1) variables found from a call to DROP1

INTEGER, INTENT(IN)    :: first, last, jmin
REAL (dp), INTENT(IN)  :: ss(:), smin

! Local variables
INTEGER    :: list(1:last), i
REAL (dp)  :: rss_last, ssq

IF (jmin == 0 .OR. last < 1 .OR. last-1 > max_size) RETURN

rss_last = rss(last)
IF (rss_last + smin > bound(last-1)) RETURN

list = vorder(1:last)
DO i = first, last-1
  vorder(i:last-1) = list(i+1:last)
  ssq = rss_last + ss(i)
  CALL report(last-1, ssq)
  vorder(i) = list(i)
END DO

RETURN
END SUBROUTINE exdrop1


SUBROUTINE forwrd(first, last, ier)

!     Forward selection from variables in positions FIRST .. LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.
!     On exit, the array VORDER contains the numbers of the variables
!     in the order in which they were added.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: pos, jmax
REAL (dp)  :: ss(last), smax

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

!     For POS = FIRST .. max_size, call ADD1 to find best variable to put
!     into position POS.

DO pos = first, max_size
  CALL add1(pos, last, ss, smax, jmax, ier)
  IF (nbest > 0) CALL exadd1(pos, smax, jmax, ss, last)

!     Move the best variable to position POS.

  IF (jmax > pos) CALL vmove(jmax, pos, ier)
END DO

RETURN
END SUBROUTINE forwrd

SUBROUTINE report(nv, ssq)

!     Update record of the best NBEST subsets of NV variables, if
!     necessary, using SSQ.

INTEGER, INTENT(IN)    :: nv
REAL (dp), INTENT(IN)  :: ssq
real (dp), allocatable :: beta(:)

!     Local variables
INTEGER    :: rank, pos1, j, list(nv), ok, ipos, idx, m, temp(nbest)
REAL (dp)  :: under1 = 0.99999_dp, above1 = 1.00001_dp



!   AAG Mod 21 Sept 2009   ----    Test to see that the combination of basis functions is legal
INTEGER, ALLOCATABLE :: sublist(:)
ALLOCATE (sublist(nv))
!sublist=vorder(1:nv)    !Yi
sublist=gorder(1:nv)    !Yi
!CALL legal(sublist,ok,0) !Commented by Yi Guo
!write(*,*),ok,offset,size(sublist)
!if(ok == 0)RETURN  !Commented by Yi Guo

!     If residual sum of squares (SSQ) for the new subset > the
!     appropriate bound, return.

IF(nv > max_size) RETURN
IF(ssq >= bound(nv)) RETURN
pos1 = (nv*(nv-1))/2 + 1

! New stuff to incorporate weights constraint test
! Invalid weight? Yes, discard it; No, good and move on.
allocate (beta(g_ptr(nv+1)-1))
call validweights(beta, g_ptr(nv+1)-1,ipos)
if (ipos < -(g_ptr(nv+1)-1)) return
! Weights constraint test ends

!     Find rank of the new subset
DO rank = 1, nbest
  IF(ssq < ress(nv,rank)*above1) THEN
    !list = vorder(1:nv)    !Yi
    list = gorder(1:nv) !Yi
    CALL shell(list, nv)

!     Check list of variables if ssq is almost equal to ress(nv,rank) -
!     to avoid including the same subset twice.

    IF (ssq > ress(nv,rank)*under1) THEN
      IF (same_vars(list, lopt(pos1:,rank), nv)) RETURN
    END IF

!     Record the new subset, and move the others down one place.

    DO j = nbest-1, rank, -1
      ress(nv,j+1) = ress(nv,j)
      lopt(pos1:pos1+nv-1, j+1) = lopt(pos1:pos1+nv-1, j)
    END DO
    ress(nv,rank) = ssq
    lopt(pos1:pos1+nv-1, rank) = list(1:nv)
    bound(nv) = ress(nv,nbest)
    ! Code for weights
    ! Record the weights and correponding variables
    temp = rank_ind(nv,1:nbest)
    m = 0 ! New code for init rank_ind = 0
    idx = 0 ! New code for init rank_ind = 0
    ! New code for init rank_ind = 0
    do j = 1, nbest
        if (temp(j) == 0) then
            idx = j
            m = nbest + 1
            exit
        end if
    end do
    ! New code for init rank_ind = 0
    if (idx==0) then ! New code for init rank_ind = 0
    call maximum(temp,nbest,m,idx)
    end if ! New code for init rank_ind = 0

    wts(wts_ptr(nv):wts_ptr(nv)+g_ptr(nv+1)-2,idx) = beta
    vwts(wts_ptr(nv):wts_ptr(nv)+g_ptr(nv+1)-2,idx) = vorder(1:g_ptr(nv+1)-1)
    where(temp >= rank .and. temp /= 0) temp = temp + 1 ! New code for init rank_ind = 0
    temp(idx) = rank
    rank_ind(nv,1:nbest) = temp
    nvwts(nv,idx) = g_ptr(nv+1)-1 - nforce
    RETURN
  END IF
END DO

RETURN
END SUBROUTINE report


SUBROUTINE originalreport(nv, ssq)

!     Update record of the best NBEST subsets of NV variables, if
!     necessary, using SSQ.

INTEGER, INTENT(IN)    :: nv
REAL (dp), INTENT(IN)  :: ssq

!     Local variables
INTEGER    :: rank, pos1, j, list(nv), ok
REAL (dp)  :: under1 = 0.99999_dp, above1 = 1.00001_dp



!   AAG Mod 21 Sept 2009   ----    Test to see that the combination of basis functions is legal
INTEGER, ALLOCATABLE :: sublist(:)
ALLOCATE (sublist(nv))
!sublist=vorder(1:nv)    !Yi
sublist=gorder(1:nv)    !Yi
!CALL legal(sublist,ok,0) !Commented by Yi Guo
!write(*,*),ok,offset,size(sublist)
!if(ok == 0)RETURN  !Commented by Yi Guo


!     If residual sum of squares (SSQ) for the new subset > the
!     appropriate bound, return.

IF(nv > max_size) RETURN
IF(ssq >= bound(nv)) RETURN
pos1 = (nv*(nv-1))/2 + 1

!     Find rank of the new subset

DO rank = 1, nbest
  IF(ssq < ress(nv,rank)*above1) THEN
    !list = vorder(1:nv)    !Yi
    list = gorder(1:nv) !Yi
    CALL shell(list, nv)

!     Check list of variables if ssq is almost equal to ress(nv,rank) -
!     to avoid including the same subset twice.

    IF (ssq > ress(nv,rank)*under1) THEN
      IF (same_vars(list, lopt(pos1:,rank), nv)) RETURN
    END IF

!     Record the new subset, and move the others down one place.

    DO j = nbest-1, rank, -1
      ress(nv,j+1) = ress(nv,j)
      lopt(pos1:pos1+nv-1, j+1) = lopt(pos1:pos1+nv-1, j)
    END DO
    ress(nv,rank) = ssq
    lopt(pos1:pos1+nv-1, rank) = list(1:nv)
    bound(nv) = ress(nv,nbest)
    RETURN
  END IF
END DO

RETURN
END SUBROUTINE originalreport

!!!!!!!!!! Subroutine REPORT goes here

SUBROUTINE shell(l, n)

!      Perform a SHELL-sort on integer array L, sorting into increasing order.

!      Latest revision - 5 July 1995

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN OUT) :: l(:)

!     Local variables
INTEGER   :: start, finish, temp, new, i1, i2, incr, it

incr = n
DO
  incr = incr/3
  IF (incr == 2*(incr/2)) incr = incr + 1
  DO start = 1, incr
    finish = n

!      TEMP contains the element being compared; IT holds its current
!      location.   It is compared with the elements in locations
!      IT+INCR, IT+2.INCR, ... until a larger element is found.   All
!      smaller elements move INCR locations towards the start.   After
!      each time through the sequence, the FINISH is decreased by INCR
!      until FINISH <= INCR.

    20 i1 = start
    temp = l(i1)
    it = i1

!      I2 = location of element new to be compared with TEMP.
!      Test I2 <= FINISH.

    DO
      i2 = i1 + incr
      IF (i2 > finish) THEN
        IF (i1 > it) l(i1) = temp
        finish = finish - incr
        EXIT
      END IF
      new = l(i2)

!     If TEMP > NEW, move NEW to lower-numbered position.

      IF (temp > new) THEN
        l(i1) = new
        i1 = i2
        CYCLE
      END IF

!     TEMP <= NEW so do not swap.
!     Use NEW as the next TEMP.

      IF (i1 > it) l(i1) = temp
      i1 = i2
      temp = new
      it = i1

!     Repeat until FINISH <= INCR.
    END DO

    IF (finish > incr) GO TO 20
  END DO

!      Repeat until INCR = 1.

  IF (incr <= 1) RETURN
END DO

RETURN
END SUBROUTINE shell



FUNCTION same_vars(list1, list2, n) RESULT(same)

LOGICAL              :: same
INTEGER, INTENT(IN)  :: n, list1(:), list2(:)

same = ALL(list1(1:n) == list2(1:n))

RETURN
END FUNCTION same_vars



SUBROUTINE seq2(first, last, ier)

! Sequential replacement algorithm applied to the variables in positions
! FIRST, ..., LAST.   2 variables at a time are added or replaced.
! If FIRST > 1, variables in positions prior to this are forced in.
! If LAST < NP, variables in positions after this are left out.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER  :: nv, nsize

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0 .OR. nbest <= 0) RETURN

nv = MIN(max_size, last-1)

!     Outer loop; SIZE = current size of subset being considered.

DO nsize = first+1, nv
  CALL replace2(first, last, nsize)
END DO

RETURN
END SUBROUTINE seq2



SUBROUTINE replace2(first, last, nsize)
! Replace 2 variables at a time from those in positions first, ..., nsize
! with 2 from positions nsize, .., last - if they reduce the RSS.

INTEGER, INTENT(IN)  :: first, last, nsize

! Local variables

INTEGER              :: ier, j1, j2, pos1, pos2, best(2), i, iwk(last)
REAL (dp)            :: smax, rssnew, rssmin, save_rss
REAL (dp), PARAMETER :: zero = 0.0_dp

10 best(1) = 0
best(2) = 0
rssmin = rss(nsize)

!     Two loops to place all pairs of variables in positions nsize-1 and nsize.
!     POS1 = destination for variable from position nsize.
!     POS2 = destination for variable from position nsize-1.

DO pos1 = first, nsize
  DO pos2 = pos1, nsize-1
    CALL add2(nsize-1, last, smax, j1, ier)
    j2 = j1 + 1
    IF (j1+j2 > nsize + nsize - 1) THEN
      rssnew = MAX(rss(nsize-2) - smax, zero)
      IF (rssnew < rssmin) THEN
        best(1) = vorder(j1)
        best(2) = vorder(j2)
        iwk(1:nsize-2) = vorder(1:nsize-2)
        rssmin = rssnew
      END IF
    END IF

    CALL vmove(nsize-1, pos2, ier)
  END DO
  CALL vmove(nsize, pos1, ier)
  DO i = pos1, nsize
    CALL report(i, rss(i))
  END DO
END DO

!     If any replacement reduces the RSS, make the best one.

IF (best(1) + best(2) > 0) THEN
  iwk(nsize-1) = best(2)
  iwk(nsize) = best(1)
  save_rss = rss(nsize)
  CALL reordr(iwk, nsize, 1, ier)
  DO i = first, nsize
    CALL report(i, rss(i))
  END DO

!    The calculated value of rssmin above is only a rough approximation to
!    the real residual sum of squares, thiugh usually good enough.
!    The new value of rss(nsize) is more accurate.   It is used below
!    to avoid cycling when several subsets give the same RSS.

  IF (rss(nsize) < save_rss) GO TO 10
END IF

RETURN
END SUBROUTINE replace2



SUBROUTINE seqrep(first, last, ier)

!     Sequential replacement algorithm applied to the variables in
!     positions FIRST, ..., LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: nv, size, start, best, from, i, jmax, count, j
REAL (dp)  :: zero = 0.0_dp, ssred, ss(last), smax

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0 .OR. nbest <= 0) RETURN

nv = MIN(max_size, last-1)

!     Outer loop; SIZE = current size of subset being considered.

DO size = first, nv
  count = 0
  start = first
  10 ssred = zero
  best = 0
  from = 0

!     Find the best variable from those in positions SIZE+1, ..., LAST
!     to replace the one in position SIZE.   Then rotate variables in
!     positions START, ..., SIZE.

  DO i = start, size
    CALL add1(size, last, ss, smax, jmax, ier)
    IF (jmax > size) THEN
      CALL exadd1(size, smax, jmax, ss, last)
      IF (smax > ssred) THEN
        ssred = smax
        best = jmax
        IF (i < size) THEN
          from = size + start - i - 1
        ELSE
          from = size
        END IF
      END IF
    END IF
    IF (i < size) CALL vmove(size, start, ier)
    DO j = start, size-1
      CALL report(j, rss(j))
    END DO
  END DO ! i = start, size

!     If any replacement reduces the RSS, make the best one.
!     Move variable from position FROM to SIZE.
!     Move variable from position BEST to FIRST.

  IF (best > size) THEN
    IF (from < size) CALL vmove(from, size, ier)
    CALL vmove(best, first, ier)
    DO j = first, best-1
      CALL report(j, rss(j))
    END DO
    count = 0
    start = first + 1
  ELSE
    count = count + 1
  END IF

!     Repeat until COUNT = SIZE - START + 1

  IF (count <= size - start) GO TO 10
END DO

RETURN
END SUBROUTINE seqrep


SUBROUTINE xhaust(first, last, ier)

!     Exhaustive search algorithm, using leaps and bounds, applied to
!     the variables in positions FIRST, ..., LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: row, i, jmax, ipt, newpos, iwk(max_size)
REAL (dp)  :: ss(last), smax, temp

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0 .OR. nbest <= 0) RETURN

!     Record subsets contained in the initial ordering, including check
!     for variables which are linearly related to earlier variables.
!     This should be redundant if the user has first called SING and
!     init_subsets.

DO row = first, max_size
  IF (d(row) <= tol(row)) THEN
    ier = -999
    RETURN
  END IF
  CALL report(row, rss(row))
END DO

!     IWK(I) contains the upper limit for the I-th simulated DO-loop for
!     I = FIRST, ..., max_size-1.
!     IPT points to the current DO loop.

iwk(first:max_size) = last

!     Innermost loop.
!     Find best possible variable for position max_size from those in
!     positions max_size, .., IWK(max_size).
30 CALL add1(max_size, iwk(max_size), ss, smax, jmax, ier)
CALL exadd1(max_size, smax, jmax, ss, iwk(max_size))

!     Move to next lower numbered loop which has not been exhausted.

ipt = max_size - 1
40 IF (ipt >= iwk(ipt)) THEN
  ipt = ipt - 1
  IF (ipt >= first) GO TO 40
  RETURN
END IF

!     Lower variable from position IPT to position IWK(IPT).
!     Record any good new subsets found by the move.

newpos = iwk(ipt)
CALL vmove(ipt, newpos, ier)
DO i = ipt, MIN(max_size, newpos-1)
  CALL report(i, rss(i))
END DO

!     Reset all ends of loops for I >= IPT.

iwk(ipt:max_size) = newpos - 1

!     If residual sum of squares for all variables above position NEWPOS
!     is greater than BOUND(I), no better subsets of size I can be found
!     inside the current loop.

temp = rss(newpos-1)
DO i = ipt, max_size
  IF (temp > bound(i)) GO TO 80
END DO
IF (iwk(max_size) > max_size) GO TO 30
ipt = max_size - 1
GO TO 40

80 ipt = i - 1
IF (ipt < first) RETURN
GO TO 40

END SUBROUTINE xhaust



SUBROUTINE random_pick(first, last, npick)
! Pick npick variables at random from those in positions first, ..., last
! and move them to occupy positions starting from first.

INTEGER, INTENT(IN)  :: first, last, npick

! Local variables

INTEGER  :: first2, i, ilist(1:last), j, k, navail
REAL     :: r

navail = last + 1 - first
IF (npick >= navail .OR. npick <= 0) RETURN
DO i = first, last
  ilist(i) = vorder(i)
END DO

first2 = first
DO i = 1, npick
  CALL RANDOM_NUMBER(r)
  k = first2 + r*navail
  IF (k > first2) THEN
    j = ilist(first2)
    ilist(first2) = ilist(k)
    ilist(k) = j
  END IF
  first2 = first2 + 1
  navail = navail - 1
END DO

CALL reordr(ilist(first:), npick, first, i)

RETURN
END SUBROUTINE random_pick

SUBROUTINE legal(vlist,ok,flag)
! Routine to analyse subgroups for data where there are two PCs as the non-backgound basis functions for each mineral
! It assumes that the basis functions are ordered with PC1 first. If VLIST contains a channel number for a PC2 it checks
! to see that the corresponding channel for PC1 is in the list. If not ok=0 else ok=1.
! NB The the number of basis functions (ncol-offset) MUST be even.

INTEGER, INTENT(IN)::vlist(:),flag
INTEGER, INTENT(OUT)::ok
INTEGER ul,n
INTEGER, ALLOCATABLE :: ulist(:)

option=1
ok=1

if(option == 0)RETURN

n=SIZE(vlist)-offset
ALLOCATE (ulist(n))
ulist(1:n)=vlist(offset+1:n+offset)
ul=(ncol-offset)/2

if(option == 1)then
  SELECT CASE(n)
  CASE (1)
    if(ulist(1) > ul) ok=0
    RETURN
  CASE (2)
    if(ABS(ulist(2)-ulist(1)) /= ul) ok=0
    RETURN
  CASE (3)
    RETURN
  CASE (4)
    CALL SHELL(ulist,n)
    if(ulist(3)-ulist(1) /= ul) ok=0
    if(ulist(4)-ulist(2) /= ul) ok=0
    GO TO 100
    !RETURN
  CASE (5)
    RETURN
  CASE (6)
    CALL SHELL(ulist,n)
    if(ulist(4)-ulist(1) /= ul) ok=0
    if(ulist(5)-ulist(2) /= ul) ok=0
    if(ulist(6)-ulist(3) /= ul) ok=0
    RETURN
  CASE (7)
    RETURN
  CASE (8)
    CALL SHELL(ulist,n)
    if(ulist(5)-ulist(1) /= ul) ok=0
    if(vlist(6)-ulist(2) /= ul) ok=0
    if(vlist(7)-ulist(3) /= ul) ok=0
    if(ulist(8)-ulist(4) /= ul) ok=0
    RETURN
  CASE (:0)
    RETURN
  CASE DEFAULT
    !WRITE(*,*),'Bad CASE parameter in LEGAL'
    !stop
  END SELECT
end if

100 if(ulist(1) == 34 .and. ulist(2) == 41 .and. ulist(3) == 94 .and. ulist(4) == 101) then
      ! write(*,*),'#######',flag,' -- ',ok
      ! write(*,*),vlist
    end if

RETURN
END SUBROUTINE legal

SUBROUTINE gxhaust_qr(first, last, ier)
! This exhaustive search code is used to test the group-wise fitting
! This code is revised from xhaust subroutine by putting group information in and rewriting the inner loop.
! Yi Guo @ CMIS 8/1/10

!     Exhaustive search algorithm, using leaps and bounds, applied to
!     the groups in positions FIRST, ..., LAST.
!     If FIRST > 1, groups in positions prior to this are forced in.
!     If LAST < ncol, groups in positions after this are forced out.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: row, i, jmax, ipt, newpos, iwk(max_size)
REAL (dp)  :: ss(last), smax, temp

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0 .OR. nbest <= 0) RETURN

!     Record subsets contained in the initial ordering, including check
!     for variables which are linearly related to earlier variables.
!     This should be redundant if the user has first called SING and
!     init_subsets.

DO row = first, max_size
  IF (d(row) <= tol(row)) THEN
    ier = -999
    RETURN
  END IF
  !CALL report(row, rss(row))   !Yi
  CALL report(row, rss(g_ptr(row+1)-1)) !Yi
END DO

!     IWK(I) contains the upper limit for the I-th simulated DO-loop for
!     I = FIRST, ..., max_size-1.
!     IPT points to the current DO loop.

iwk(first:max_size) = last

!     Innermost loop.
!     Find best possible variable for position max_size from those in
!     positions max_size, .., IWK(max_size).
30 CALL gadd1(max_size, iwk(max_size), ier) !repetance of jmax is only for convenience. Should be removed in final version.
!CALL exadd1(max_size, smax, jmax, ss, iwk(max_size))

!     Move to next lower numbered loop which has not been exhausted.

ipt = max_size - 1
40 if (ipt<1) return
IF (ipt >= iwk(ipt)) THEN
  ipt = ipt - 1
  IF (ipt >= first) GO TO 40
  RETURN
END IF

!     Lower variable from position IPT to position IWK(IPT).
!     Record any good new subsets found by the move.

newpos = iwk(ipt)
! insert debuging code Yi
!write(*,*), "gorder:",gorder
CALL gmove(ipt, newpos, ier)
DO i = ipt, MIN(max_size, newpos-1)
  !CALL report(i, rss(i))   !Yi
  CALL report(i, rss(g_ptr(i+1)-1)) !Yi
END DO

!     Reset all ends of loops for I >= IPT.

iwk(ipt:max_size) = newpos - 1

!     If residual sum of squares for all variables above position NEWPOS
!     is greater than BOUND(I), no better subsets of size I can be found
!     inside the current loop.

temp = rss(g_ptr(newpos)-1)
DO i = ipt, max_size
  IF (temp > bound(i)) GO TO 80
END DO
IF (iwk(max_size) > max_size) GO TO 30
ipt = max_size - 1
GO TO 40

80 ipt = i - 1
IF (ipt < first) RETURN
GO TO 40

END SUBROUTINE gxhaust_qr

SUBROUTINE allsubsets(p, k, count)
! This exhaustive search code is used to test the group-wise fitting
! This code is revised from xhaust subroutine by putting group information in and rewriting the inner loop.
! Yi Guo @ CMIS 8/1/10

!     Exhaustive search algorithm, using leaps and bounds, applied to
!     the groups in positions FIRST, ..., LAST.
!     If FIRST > 1, groups in positions prior to this are forced in.
!     If LAST < ncol, groups in positions after this are forced out.

INTEGER, INTENT(IN)  :: p, k
INTEGER, INTENT(OUT) :: count
integer :: tmp
integer, allocatable, save :: v(:)

!     Local variables

INTEGER    :: row, i, jmax, ipt, newpos,first, iwk(p)

count = 0
first = 1

allocate(v(k))
do i = 1, k
    v(i) = i
end do
iwk(1:p) = k

ipt = p
! write(*,*),"-------------------------------------"
do
   !  write(*,*), v
    count = count + 1
    if (ipt < iwk(ipt) .and. iwk(ipt) > p) then
        tmp = v(ipt)
        v(ipt:iwk(ipt)-1) = v(ipt+1:iwk(ipt))
        v(iwk(ipt)) = tmp
        iwk(ipt) = iwk(ipt) - 1
        iwk(ipt+1:p) = iwk(ipt)
        ipt = p
    else
        ipt = ipt -1
        if (ipt < 1) then
            return
        end if
    end if
end do

END SUBROUTINE allsubsets
! New group exhaustive search code ends here

SUBROUTINE gadd1(first, last, ier)

! Calculate the reduction in residual sum of squares when one variable,
! selected from those in positions FIRST .. LAST, is added in position FIRST,
! given that the variables in positions 1 .. FIRST-1 (if any) are already
! included.

! call qrupdate to get RSS reduction of including each group

INTEGER, INTENT(IN)     :: first, last
INTEGER, INTENT(OUT)    :: ier

!     Local variables

INTEGER    :: j, inc, pos, row, col,i, cnv, nr, m, idx_first, idx_last, ltemp, fposwts, maxval, maxidx, wtsdim, idx(nbest) !weights
REAL (dp)  :: zero = 0.0_dp, diag, dy, ssqx, sxx(ncol), sxy(ncol), wsserr, wy, ss(nbest)
real (dp), allocatable :: wr(:), wd(:),wrhs(:), wx(:), beta(:),ubeta(:), varswts(:) !weights
REAL (dp)               :: weight = 1.0_dp, rssred, rtemp !weights
integer, allocatable :: wrow_ptr(:), subvars(:,:), subnvars(:), subrank(:), subloc(:) !Locators of weights and variables, weights code !weights
real (dp), allocatable :: subwts(:,:) ! weights

!     Check call arguments

ier = 0
IF (first > ncol) ier = 1
IF (last < first) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

ss = -1
idx = 0
! Yi, weights code
call maximum(gnv,ng, maxval, maxidx)
wtsdim = first*maxval+nforce
allocate (subwts(wtsdim,nbest), subvars(wtsdim,nbest),subrank(nbest), subnvars(nbest),subloc(nbest))
subrank = 1
subnvars = 0
! Yi, weights code

!write(*,*),gorder,g_ptr,row_ptr !make things visible

!! The following 3 lines are for debugging. Output data and read those in R to compare results.
!OPEN(10, FILE='data.txt')
!write(10,*),nforce,ncol,r,d,rhs,ng,g_ptr,first,last
!close(10)
!! Debugging code ends here.

do i = first, last
    cnv = gnv(gorder(i))
    nr = cnv*(cnv-1)/2
    if (allocated(wr)) deallocate(wr,wd,wrhs,wx)
    allocate (wr(nr),wd(cnv),wrhs(cnv),wx(cnv))
    idx_first = g_ptr(i)
    idx_last = g_ptr(i+1)-1
    pos = 1
    inc = cnv - 1
    wr = 0
    wd = 0
    wrhs = 0
    wsserr = 0 !wsserr is almost a dummy variable. It is supposed to record the rss after regress out all variables. But in order to speed things up
               !a bit the following loop is from g_ptr(first) to g_ptr(i+1)-1 (skipped the rows which are all zero.
    rssred = 0
    do m=g_ptr(first),g_ptr(i+1)-1 !Change g_ptr(i+1)-1 to g_ptr(last+1)-1 will activate wsserr. But not necessary here.
        wx = 0
        if (m<=idx_last) then
        if (m>=idx_first) then
            if(m>=idx_last) then
                wx(cnv) = 1
            else
                wx(m-idx_first+1) = 1
                wx(m-idx_first+2:cnv) = r(row_ptr(m):row_ptr(m)+cnv-m+idx_first-2)
            end if
        else
            wx = r(row_ptr(m)+idx_first-1-m:row_ptr(m)+idx_first-1-m+cnv-1)
        end if
        end if
        wx = wx*sqrt(d(m))
        wy = rhs(m)*sqrt(d(m))
        call qrgupdate(weight, wx, wy, cnv, wd, wr, wrhs, wsserr) ! weight here may be a potential bug. If it is not 1, it could bring some trouble.
    end do
    ! Get the rss reduction here.
    do m = 1, cnv
        rssred = rssred + wd(m)*wrhs(m)**2
    end do

    if ((rss(g_ptr(first)-1)-rssred)>=bound(first)) then
        deallocate (wr, wd, wrhs, wx)
        cycle
    end if

    ! code for weights
    if (allocated(beta)) deallocate(beta)
    allocate (beta(cnv+g_ptr(first)-1))
    if (allocated(wrow_ptr)) deallocate(wrow_ptr)
    allocate(wrow_ptr(cnv))
    wrow_ptr = 0
    if(cnv>1) wrow_ptr(1) = 1
    DO m = 2, cnv-1
      wrow_ptr(m) = wrow_ptr(m-1) + cnv - m + 1
    END DO

!    call gaddposweights(beta(g_ptr(first)), cnv, corder(g_ptr(i):g_ptr(i+1)-1), wr, wrow_ptr, wrhs, fposwts)
    call gaddwtsconstraint(beta(g_ptr(first)), cnv, &
    corder(g_ptr(i):g_ptr(i+1)-1), &
    wtslb(g_ptr(i):g_ptr(i+1)-1), wr, wrow_ptr, &
    wrhs, fposwts)
    if (fposwts>=-9*cnv) then ! If at least the weights for the vars in group to be added are ok, then get the other weights
       if (allocated(wx)) deallocate(wx)
       allocate (wx(g_ptr(first)-1))
        do m = 1, g_ptr(first)-1
            rtemp = 0
            do j = 1, cnv
                rtemp = rtemp + r(row_ptr(m)+idx_first+j-2-m) * beta(g_ptr(first)-1+j)
            end do
            wx(m) = rhs(m) - rtemp
        end do
!        call gaddposweights(beta, g_ptr(first)-1, corder(1:g_ptr(first)-1), r, row_ptr, wx, fposwts) ! Get the weights of groups before first.
	call gaddwtsconstraint(beta, g_ptr(first)-1, &
	corder(1:g_ptr(first)-1), wtslb(1:g_ptr(first)-1),&
	r, row_ptr, wx, fposwts)

        if (fposwts >= -9 * cnv) then ! Worth going on
            ! if weights are ok, the next is to put things together, say weights, variables in order.
            ! However, in order to make it efficient, it is necessary to check its rank first. But for debugging purpose, test the following first.
            if (allocated(varswts)) deallocate(varswts)
            allocate (varswts(cnv+g_ptr(first)-1))
            varswts =(/ vorder(1:g_ptr(first)-1), vorder(g_ptr(i):g_ptr(i+1)-1) /) ! Assemble the variables.

            ! If the weights are all valid, go ahead and update the data for nbest list
            ! Update nbest largest RSS reduction
            do m=1,nbest
                if (rssred>ss(m)) then
                    if (m < nbest) then
                        ss(m+1:nbest) = ss(m:nbest-1)
                        idx(m+1:nbest) = idx(m:nbest-1)
                        subloc(m+1:nbest) = subloc(m:nbest-1) !weight code
                    end if
                    ss(m) = rssred
                    idx(m) = gorder(i)

                     ! Update data structure
                    call maximum(subrank,nbest,maxval,maxidx)
                    subwts(1:cnv+g_ptr(first)-1,maxidx) = beta
                    subvars(1:cnv+g_ptr(first)-1,maxidx) = varswts
                    where(subrank >= m) subrank = subrank + 1
                    subrank(maxidx) = m
                    subnvars(maxidx) = cnv+g_ptr(first)-1 - nforce
                    subloc(m) = maxidx

                    exit
                end if
            end do
        end if
    end if ! temporary stop if here.
    ! code for weights
    deallocate (wr, wd, wrhs, wx)
 end do

! Update the ress, lop, bound etc.
ltemp = gorder(first)
do i= 1, nbest
    ! It is possible that less than nbest combinations are found, for example first = 5, last = 6, nbest =4.
    if (ss(i)>=0) then
        gorder(first) = idx(i)
        !call report(first,rss(g_ptr(first)-1)-ss(i)) ! update ! Commented for weights update
        call gaddreport(first,rss(g_ptr(first)-1)-ss(i), subwts(:,subloc(i)), subvars(:,subloc(i)), subnvars(subloc(i))) !Weights update
    end if
end do
gorder(first) = ltemp


! The following code is trying to make the calculation more efficient however, it seems some bugs are in it or the theory is wrong.
!do i = first, last
!    cnv = gnv(gorder(i))
!    nr = cnv*(cnv-1)/2
!    allocate (wr(nr),wd(cnv),wrhs(cnv),wx(cnv))
!    idx_first = g_ptr(i)
!    idx_last = g_ptr(i+1)-1
!    pos = 1
!    inc = cnv - 1
!    do m = idx_first, idx_last - 1
!        wr(pos:pos+inc-1) = r(row_ptr(m):row_ptr(m)+inc-1) * sqrt(d(m))
!        wrhs = rhs(m) * sqrt(d(m))
!        pos = pos + inc
!        inc = inc - 1
!    end do
!    wd = d(idx_first:idx_first+cnv-1)
!    !wrhs = rhs(idx_first:idx_first+cnv-1)
!    wsserr = 0
!    rssred = 0
!    do m=idx_first+cnv,ncol
!        wsserr = wsserr + d(m)*rhs(m)**2
!    end do
!    do m=g_ptr(first),g_ptr(i)-1
!        wx = r(row_ptr(m)+idx_first-1-m:row_ptr(m)+idx_first-1-m+cnv-1)
!        wy = rhs(m)
!        call qrgupdate(weight, wx, wy, cnv, wd, wr, wrhs, wsserr)
!    end do
!    do m = 1, cnv
!        rssred = rssred + wd(m)*wrhs(m)**2
!    end do
!    deallocate (wr, wd, wrhs, wx)
!end do

END SUBROUTINE gadd1

SUBROUTINE gxhaust_add2(first, last, ier)
! This exhaustive search code is used to test the group-wise fitting
! This code is revised from xhaust subroutine by putting group information in and rewriting the inner loop.
! Yi Guo @ CMIS 8/1/10

!     Exhaustive search algorithm, using leaps and bounds, applied to
!     the groups in positions FIRST, ..., LAST.
!     If FIRST > 1, groups in positions prior to this are forced in.
!     If LAST < ncol, groups in positions after this are forced out.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: row, i, jmax, ipt, newpos, iwk(max_size)
REAL (dp)  :: ss(last), smax, temp

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0 .OR. nbest <= 0) RETURN

!     Record subsets contained in the initial ordering, including check
!     for variables which are linearly related to earlier variables.
!     This should be redundant if the user has first called SING and
!     init_subsets.

DO row = first, max_size
  IF (d(row) <= tol(row)) THEN
    ier = -999
    RETURN
  END IF
  !CALL report(row, rss(row))   !Yi
  CALL report(row, rss(g_ptr(row+1)-1)) !Yi
END DO

!     IWK(I) contains the upper limit for the I-th simulated DO-loop for
!     I = FIRST, ..., max_size-1.
!     IPT points to the current DO loop.

iwk(first:max_size) = last

!     Innermost loop.
!     Find best possible variable for position max_size from those in
!     positions max_size, .., IWK(max_size).
30 CALL add2(max_size, iwk(max_size), smax, jmax, ier) !repetance of jmax is only for convenience. Removed!
!CALL exadd1(max_size, smax, jmax, ss, iwk(max_size))

!     Move to next lower numbered loop which has not been exhausted.

ipt = max_size - 1
40 if (ipt<1) return
IF (ipt >= iwk(ipt)) THEN
  ipt = ipt - 1
  IF (ipt >= first) GO TO 40
  RETURN
END IF

!     Lower variable from position IPT to position IWK(IPT).
!     Record any good new subsets found by the move.

newpos = iwk(ipt)
! insert debuging code Yi
!write(*,*), "gorder:",gorder
CALL gmove(ipt, newpos, ier)
DO i = ipt, MIN(max_size, newpos-1)
  !CALL report(i, rss(i))   !Yi
  CALL report(i, rss(g_ptr(i+1)-1)) !Yi
END DO

!     Reset all ends of loops for I >= IPT.

iwk(ipt:max_size) = newpos - 1

!     If residual sum of squares for all variables above position NEWPOS
!     is greater than BOUND(I), no better subsets of size I can be found
!     inside the current loop.

temp = rss(g_ptr(newpos)-1)
DO i = ipt, max_size
  IF (temp > bound(i)) GO TO 80
END DO
IF (iwk(max_size) > max_size) GO TO 30
ipt = max_size - 1
GO TO 40

80 ipt = i - 1
IF (ipt < first) RETURN
GO TO 40

END SUBROUTINE gxhaust_add2

SUBROUTINE gxhaust_gmove(first, last, ier)
! This exhaustive search code is used to test the group-wise fitting
! This code is revised from xhaust subroutine by putting group information in and rewriting the inner loop.
! Yi Guo @ CMIS 8/1/10

!     Exhaustive search algorithm, using leaps and bounds, applied to
!     the groups in positions FIRST, ..., LAST.
!     If FIRST > 1, groups in positions prior to this are forced in.
!     If LAST < ncol, groups in positions after this are forced out.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: row, i, jmax, ipt, newpos, iwk(max_size)
REAL (dp)  :: ss(last), smax, temp

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0 .OR. nbest <= 0) RETURN

!     Record subsets contained in the initial ordering, including check
!     for variables which are linearly related to earlier variables.
!     This should be redundant if the user has first called SING and
!     init_subsets.

DO row = first, max_size
  IF (d(row) <= tol(row)) THEN
    ier = -999
    RETURN
  END IF
  !CALL report(row, rss(row))   !Yi
  CALL report(row, rss(g_ptr(row+1)-1)) !Yi
END DO

!     IWK(I) contains the upper limit for the I-th simulated DO-loop for
!     I = FIRST, ..., max_size-1.
!     IPT points to the current DO loop.

iwk(first:max_size) = last

!     Innermost loop.
!     Find best possible variable for position max_size from those in
!     positions max_size, .., IWK(max_size).
!30 CALL add1(max_size, iwk(max_size), ss, smax, jmax, ier)
!CALL exadd1(max_size, smax, jmax, ss, iwk(max_size))

!     Move to next lower numbered loop which has not been exhausted.

30 ipt = max_size
40 if (ipt<1) return
IF (ipt >= iwk(ipt)) THEN
  ipt = ipt - 1
  IF (ipt >= first) GO TO 40
  RETURN
END IF

!     Lower variable from position IPT to position IWK(IPT).
!     Record any good new subsets found by the move.

newpos = iwk(ipt)
! insert debuging code Yi
!write(*,*), "gorder:",gorder
CALL gmove(ipt, newpos, ier)
DO i = ipt, MIN(max_size, newpos-1)
  !CALL report(i, rss(i))   !Yi
  CALL report(i, rss(g_ptr(i+1)-1)) !Yi
END DO

!     Reset all ends of loops for I >= IPT.

iwk(ipt:max_size) = newpos - 1

!     If residual sum of squares for all variables above position NEWPOS
!     is greater than BOUND(I), no better subsets of size I can be found
!     inside the current loop.

temp = rss(g_ptr(newpos)-1)
DO i = ipt, max_size
  IF (temp > bound(i)) GO TO 80
END DO
IF (iwk(max_size) > max_size) GO TO 30
ipt = max_size - 1
GO TO 40

80 ipt = i - 1
IF (ipt < first) RETURN
GO TO 40

END SUBROUTINE gxhaust_gmove

! Add a subroutine to get the maximum of an array.
! Yi 30/01/10
subroutine maximum(a,size,m,idx)
integer, intent (in) :: size
integer, intent (in) :: a(size)
integer, intent (out) :: m, idx
integer :: j

m=a(1)
idx = 1
do j=2,size
    if (m.LT.a(j)) then
        m=a(j)
        idx = j
    end if
end do
END subroutine maximum

SUBROUTINE gaddreport(nv, ssq, subwts, subvars, subnvars)

!     Update record of the best NBEST subsets of NV variables, if
!     necessary, using SSQ.

INTEGER, INTENT(IN)    :: nv, subvars(:), subnvars
REAL (dp), INTENT(IN)  :: ssq, subwts(:)

!     Local variables
INTEGER    :: rank, pos1, j, list(nv), ok, ipos, idx, m, temp(nbest)
REAL (dp)  :: under1 = 0.99999_dp, above1 = 1.00001_dp



!   AAG Mod 21 Sept 2009   ----    Test to see that the combination of basis functions is legal
INTEGER, ALLOCATABLE :: sublist(:)
ALLOCATE (sublist(nv))
!sublist=vorder(1:nv)    !Yi
sublist=gorder(1:nv)    !Yi
!CALL legal(sublist,ok,0) !Commented by Yi Guo
!write(*,*),ok,offset,size(sublist)
!if(ok == 0)RETURN  !Commented by Yi Guo

!     If residual sum of squares (SSQ) for the new subset > the
!     appropriate bound, return.

IF(nv > max_size) RETURN
IF(ssq >= bound(nv)) RETURN
pos1 = (nv*(nv-1))/2 + 1

!! New stuff to incorporate positive weights test
!! Invalid weight? Yes, discard it; No, good and move on.
!allocate (beta(g_ptr(nv+1)-1))
!call validweights(beta, g_ptr(nv+1)-1,ipos)
!if (ipos < -(g_ptr(nv+1)-1)) return
!! Positive weights test ends

!     Find rank of the new subset
DO rank = 1, nbest
  IF(ssq < ress(nv,rank)*above1) THEN
    !list = vorder(1:nv)    !Yi
    list = gorder(1:nv) !Yi
    CALL shell(list, nv)

!     Check list of variables if ssq is almost equal to ress(nv,rank) -
!     to avoid including the same subset twice.

    IF (ssq > ress(nv,rank)*under1) THEN
      IF (same_vars(list, lopt(pos1:,rank), nv)) RETURN
    END IF

!     Record the new subset, and move the others down one place.

    DO j = nbest-1, rank, -1
      ress(nv,j+1) = ress(nv,j)
      lopt(pos1:pos1+nv-1, j+1) = lopt(pos1:pos1+nv-1, j)
    END DO
    ress(nv,rank) = ssq
    lopt(pos1:pos1+nv-1, rank) = list(1:nv)
    bound(nv) = ress(nv,nbest)
    ! Code for weights
    ! Record the weights and correponding variables
    temp = rank_ind(nv,1:nbest)
    m = 0 ! New code for init rank_ind = 0
    idx = 0 ! New code for init rank_ind = 0
    ! New code for init rank_ind = 0
    do j = 1, nbest
        if (temp(j) == 0) then
            idx = j
            m = nbest + 1
            exit
        end if
    end do
    ! New code for init rank_ind = 0
    if (idx==0) then ! New code for init rank_ind = 0
    call maximum(temp,nbest,m,idx)
    end if ! New code for init rank_ind = 0
    wts(wts_ptr(nv):wts_ptr(nv)+subnvars+nforce-1,idx) = subwts(1:subnvars+nforce)
    vwts(wts_ptr(nv):wts_ptr(nv)+subnvars+nforce-1,idx) = subvars(1:subnvars+nforce)
    where(temp >= rank .and. temp /= 0) temp = temp + 1 ! New code for init rank_ind = 0
    temp(idx) = rank
    rank_ind(nv,1:nbest) = temp
    nvwts(nv,idx) = subnvars
    RETURN
  END IF
END DO

RETURN
END SUBROUTINE gaddreport

END MODULE find_subsets
subroutine gss(data,library,consvec,consval,ncase, nvartotal, &
                     nfixed, nvar_max, nbestmix, numberofgroups,&
                     nvaringroups,outgroups,outrss,outwts,&
!                     outwtsvars,& ! To be deleted
!                     outwtsnvars,outwtsrank,& ! To be deleted
		     noutwts, comptime)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'gss_' :: gss

! N.B. First line is necessary for building Windows DLL but no need in UNIX.
!
USE lsq
USE find_subsets

IMPLICIT NONE
INTEGER, INTENT(IN)  :: numberofgroups,nvaringroups(numberofgroups),nbestmix,&
                          ncase, nvartotal,nvar_max,nfixed
integer, intent (in) :: noutwts, consvec(nvartotal+nfixed)
real (dp), intent (in) :: data(ncase),library(nvartotal+nfixed,ncase), &
			consval(nvartotal+nfixed)
real (dp), intent (out) :: outrss(nvar_max*nbestmix),comptime
real (dp), intent (out) :: outwts(noutwts) !weights
INTEGER, INTENT(OUT) :: outgroups((nvar_max+1)*nvar_max/2*nbestmix)!, &
!                          outwtsvars(noutwts), outwtsnvars(nvar_max*nbestmix),&
!                          outwtsrank(nvar_max*nbestmix)

CHARACTER (LEN = 40) :: fname_dat, fname_red, fname_rpt, fname_ctl, fname_y,lstring
CHARACTER (LEN = 1)  :: ans, bel = CHAR(7), yesno
CHARACTER (LEN = 10) :: string
INTEGER              :: first, last, ier, i, nsize, pos, j,  &
                        ypos, in, dimc, nv, i0, best_size, line1, nrepl,    &
                        search_method, criterion, iostatus, rank_deficit, ndf, &
                        iy, fc,npc,nnn,fdisp, nvar, iflag !Yi
INTEGER, ALLOCATABLE :: list(:), seed(:), order_copy(:), force(:)
LOGICAL, ALLOCATABLE :: lindep(:)

REAL (dp)            :: one = 1.0, fin, fout, total_sumsq, r2, msep, press, &
                        y, e, h, eps
REAL (dp), ALLOCATABLE  :: cormat(:), ycorr(:), beta(:), x(:), xcopy(:)
!REAL (dp), ALLOCATABLE  :: library(:,:)	, data(:)
LOGICAL              :: lsel = .false., OK, cvd,Recompute_QR, firstpass
REAL                 :: var, Cp, Cp_last, fmax, f1, f5, f10, zero = 0.0
real (dp) :: t_start, t_finish
real (dp), ALLOCATABLE :: tmpwts(:)
INTEGER, ALLOCATABLE :: tmpvars(:)
INTEGER :: tmpnvars, k
!-------------------------------------------------------------------------------------------------------
!  Settings for batch operation
!
cvd=.FALSE.							! Cross validation
Recompute_QR=.TRUE.					! Switch to use existing QR file. FALSE uses existing file. TRUE creates it.

!eps=1.0E-10

!-------------------------------------------------------------------------------------------------------
nbest = nbestmix
nforce = nfixed
ng = numberofgroups
IF (ALLOCATED(gnv)) deallocate(gnv)
allocate (gnv(ng))
gnv = nvaringroups
if (allocated(corder)) deallocate(corder)
allocate (corder(nvartotal+nfixed))
corder = consvec
if (allocated(wtslb)) deallocate(wtslb)
allocate (wtslb(nvartotal+nfixed))
wtslb = consval

!    ! The following 3 lines are for debugging. Output data and read those in R to compare results.
!    OPEN(dbfid, FILE='trace.txt')
!    write(dbfid,*),wtslb
!    ! debugging code ends

outrss = 0
outgroups = 0
outwts = 0
!outwtsvars = 0
!outwtsnvars = 0
!outwtsrank = 0
! The first nforce(nfixed) variables are always in the front if nforce (nfixed) is larger than 0.
if(nforce .ne. 0)then
  IF (ALLOCATED(force)) deallocate(force)
  ALLOCATE (force(nforce))
  do i=1,nforce
    force(i) = i
  end do
end if


nvar=nvartotal+nforce

firstpass=.TRUE.

CALL start()                 ! Read in data and form QR reduction, if necessary

CALL init_subsets(nvar_max, .false.)
i0 = 0
total_sumsq = rss(1) + d(1)*rhs(1)**2

call cpu_time(t_start)
CALL gxhaust_qr(first, last, ier) !Use local QR only as it is the fastest.
!select case (alg)
!    case (0)
!        CALL gxhaust_qr(first, last, ier)   !Use qr in innermost loop
!    case (1)
!        CALL gxhaust_gmove(first, last, ier)   !Use gmove in all loops
!    case (2)
!        CALL gxhaust_add2(first, last, ier) !Use add2 in innermost loop. N.B. only suitable for 2 variables a group for all groups case.
!    case default
!        write(*,*), "Unknown algorithm in group subset selection!" !This may not be able to print in DLL mode.
!end select
call cpu_time(t_finish) !record the finish time.
comptime = (t_finish - t_start)*1
!j = max_size*(max_size-1)/2 + 1

!! Output best found
!DO i = first, max_size
!    DO j = 1, nbest
!      outrss((i-1)*nbest+j) = ress(i,j)
!    END DO
!END DO
!pos = 1
!i0 = (max_size+1)*max_size/2
!do j = 1, nbest
!    outgroups(pos:pos+i0-1) = lopt(1:i0,j)
!    pos = pos + i0
!end do
!
!pos = size(wts,1)
!do i = 1, nbest
!    outwtsnvars((i-1)*nvar_max+1:i*nvar_max) = nvwts(:,i)
!    outwtsrank((i-1)*nvar_max+1:i*nvar_max) = rank_ind(:,i)
!    outwts((i-1)*pos+1:i*pos) = wts(:,i)
!    outwtsvars((i-1)*pos+1:i*pos) = vwts(:,i)
!end do

! Insert re-out code frmo here
DO i = first, max_size
    DO j = 1, nbest
      outrss((i-1)*nbest+j) = ress(i,j)
    END DO
END DO

pos = 1
i0 = 1
do i = 1, max_size
	i0 = i0 + i - 1
	do j = 1, nbest
		outgroups(pos:pos+i-1) = lopt(i0:i0+i-1,j)
		pos = pos + i
	end do
end do

pos = size(wts,1)
do k = 1, nvar_max
	do i =1, nbest
		iflag = -1
		do j = 1, nbest
			if (rank_ind(k,j) == i) then
				iflag = j
				exit
			end if
		end do
		if (iflag /= -1)  then! Then iflag is the index of the i-th best
			tmpnvars = nvwts(k,iflag) + nforce

			IF (ALLOCATED(tmpvars)) DEALLOCATE(tmpvars)
			ALLOCATE ( tmpvars(tmpnvars) )
			IF (ALLOCATED(tmpwts)) DEALLOCATE(tmpwts)
			ALLOCATE ( tmpwts(tmpnvars) )

			tmpvars = vwts(wts_ptr(k):wts_ptr(k)+tmpnvars-1,iflag)
			tmpwts = wts(wts_ptr(k):wts_ptr(k)+tmpnvars-1,iflag)

			call Acc_Shell_Sort(tmpvars,tmpwts)
			outwts((i-1)*pos+wts_ptr(k):(i-1)*pos+wts_ptr(k)+tmpnvars-1) = tmpwts
		end if
	end do
end do
! Insert re-out code done

DEALLOCATE(bound, ress, lopt, wts, wts_ptr,vwts, rank_ind, nvwts)
!d, rhs, r, tol, rss, wtslb, gorder, g_ptr, gnv)
return
CONTAINS

SUBROUTINE start()

!     This is the starting routine for the SUBSETS package of programs.
!     If a QR reduction has not been formed, then it forms the
!     upper-triangular Banachiewicz factorization of the input data.
!     Free-format input is assumed, i.e. with data fields separated by
!     spaces, CR's, tabs or commas.   N.B. Some Fortran compilers will
!     not accept tabs and/or commas as delimiters.

!     Latest revision - 10 May 2001

CHARACTER (LEN=20)    :: response
REAL (dp)             :: eps
REAL (dp), PARAMETER  :: pt0001 = 0.0001
INTEGER               :: i

!     Does the user already have a QR factorization?

CALL form_QR()

! Set tolerances and test for singularities

eps = EPSILON(pt0001) ** 0.667

CALL tolset(eps)
IF (ALLOCATED(lindep)) DEALLOCATE(lindep)
ALLOCATE ( lindep(1:nvar) )

CALL sing(lindep, ier)
rank_deficit = -ier

IF (nobs > ncol) THEN
  ndf = nobs - ncol + rank_deficit
END IF

!     Set default values for first & last
!     first = first variable which may an be omitted from subsets
!     last  = last variable which may be included in subsets

first = 1
!last = ncol    !Yi
last = ng  !Yi

RETURN
END SUBROUTINE start

SUBROUTINE form_QR()
! Form the QR-factorization from input arguments

REAL (dp)               :: weight = 1.0_dp, y
REAL (dp), ALLOCATABLE  :: x(:)
CHARACTER (LEN = 79)    :: text
INTEGER                 :: ipos, i, iostatus, ier
INTEGER                 :: dex
REAL, ALLOCATABLE       :: xmin(:), xmax(:), xmean(:), xstdev(:)
REAL                    :: ymin, ymax, ymean, ystdev

i0 = 0


!     Allocate memory for arrays X and VNAME.

ALLOCATE ( x(nvar))

! Initialize public variables for QR
CALL startup(nvar, .false.)          ! From LSQ module

!     Allocate arrays for storing descriptive statistics
!ALLOCATE( xmin(nvar), xmax(nvar), xmean(nvar), xstdev(nvar) )

!     Read in data and form the upper-triangular factorization.

dex=1
DO
  IF(firstpass) THEN
    if(dex > ncase) EXIT
    if(.not. allocated(x)) allocate(x(nvar))
    do i=1,nvar
      x(i)=library(i,dex)
    end do
  END IF
! Update descriptive statistics
  y = data(dex)
  CALL includ(weight, x, y)
  dex=dex+1
END DO

!     Test for singularities & form initial array (rss) of sums of squares
!     of residuals

if(.not. allocated(lindep)) ALLOCATE ( lindep(ncol) )
CALL sing(lindep, ier)
CALL ss()

!     Change extension to .red for output file.
!DEALLOCATE ( x, lindep, xmin, xmax, xmean, xstdev )
deallocate (x)
RETURN
END SUBROUTINE form_QR

SUBROUTINE Acc_Shell_Sort(a,b)

  IMPLICIT NONE
  INTEGER :: i, j, increment
  real(dp) :: tempb
  INTEGER :: temp
  INTEGER, INTENT(in out) :: a(:)
  real(dp), INTENT(in out) :: b(:)

  increment = SIZE(a) / 2
  DO WHILE (increment > 0)
      DO i = increment+1, SIZE(a)
         j = i
         temp = a(i)
		 tempb = b(i)
         DO WHILE (j >= increment+1 .AND. a(j-increment) > temp)
            a(j) = a(j-increment)
			b(j) = b(j-increment)
            j = j - increment
         END DO
         a(j) = temp
		 b(j) = tempb
      END DO
      IF (increment == 2) THEN
   	  increment = 1
      ELSE
         increment = increment * 5 / 11
      END IF
  END DO

END SUBROUTINE Acc_Shell_Sort

END subroutine gss
