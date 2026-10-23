!> Interface-mode amplitude: the correct observable for an interfacial instability.  Fortran counterpart of plotting/interface_mode.py.
!>
!> WHY THIS MODULE EXISTS (learned from running the solver)
!! --------------------------------------------------------
!! The obvious NS-MFP observable, |Q'(t)| = |Q_pert - Q_base|, fails for Rayleigh-Plateau in two different ways depending on the field chosen, and both failures only became visible with real solver output:
!!
!! 1. VELOCITY FIELDS give Q'(0) = 0 exactly.  The perturbation lives on the interface, so both twins start from rest.  That breaks the DMD amplitude ranking, which projects each mode onto the FIRST snapshot: all amplitudes come out zero and the ordering is meaningless. Velocities are still a valid state vector, but only after the startup transient is skipped.
!!
!! 2. THE RAW LEVEL SET dilutes the signal.  L0101 is a signed distance function, so Q' ~ eps*sin(k z) across the WHOLE domain while only the near-interface part participates in the instability.  In a real run the raw level-set norm stayed flat to 0.6% while the interface amplitude genuinely grew by 37%.  A flat norm does NOT prove a stable interface - it can mean the observable is wrong.
!!
!! The fix is to measure what the dispersion relation predicts: the amplitude of the k-Fourier component of the interface radius r(z,t).  For r = r0 + a(t) cos(k z), linear theory gives a(t) = a(0) exp(sigma t), so a log-linear fit of a(t) returns sigma directly, immune to both the zero initial condition and the static-offset dilution.
!!
!! Use alongside DMD, not instead of it: DMD gives the spectrum and mode shapes, this gives a robust scalar growth rate to check it against.
!> ORIENTATION
!>
!> This module measures the shape of the interface, which is the quantity the growth rate is fitted to.
!>
!> The solver stores the interface implicitly, as the zero contour of a level-set field that is positive inside the liquid and negative outside (the deck template sets phi = -|r| + 1 + radblob sin(2 pi z/yblob)); only the sign change is used, so the method does not depend on this convention. To find the interface radius at one axial station, walk outwards along that row until the sign changes, then interpolate linearly between the two straddling cell centres. Repeating at every station gives the profile r(z), and the amplitude of its first Fourier harmonic is the observable the analysis uses.
!>
!> That observable was chosen deliberately. A velocity-based measure has no signal at the first instant, because both twin calculations start from rest and their difference is identically zero. The raw level-set difference spreads the discrepancy across the whole domain and dilutes it. The interface amplitude is the quantity the theory predicts.
module interface_mode_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use plotfile_mod
   use dispersion_mod, only: PI
   use util_mod, only: fit_log_linear, cell_centres, argsort_real_asc
   implicit none
   private

   public :: dp
   public :: interface_radius, mode_amplitude, amplitude_series
   public :: fit_growth_rate_amp, local_growth_rate, has_converged

contains

   !> Interface radius r(z) from the zero contour of the level set.
   !>
   !> Linear interpolation of the first sign change in each axial column. Columns with no sign change (interface outside the domain) are set to a negative sentinel rather than a silently wrong value. A cell centre at which the level set is exactly zero is the interface itself; the former test a*b < 0 skipped it and took the next crossing or none. Coordinates are built from prob_lo and prob_hi, not from prob_hi alone.
   subroutine interface_radius(pltdir, ls_field, time, z, radius, nz, ierr)
      character(len=*), intent(in) :: pltdir, ls_field
      real(dp), intent(out) :: time
      real(dp), allocatable, intent(out) :: z(:), radius(:)
      integer, intent(out) :: nz, ierr

      type(plotfile_t) :: pf
      character(len=64) :: want(1)
      integer :: i, j, nr, ihit
      real(dp) :: denom, f
      real(dp), allocatable :: r(:)

      want(1) = ls_field
      call read_plotfile(pltdir, pf, ierr, want)
      if (ierr /= 0) then
         nz = 0
         return
      end if

      nr = pf%nx
      nz = pf%ny
      time = pf%time

      allocate(r(nr), z(nz), radius(nz))
      call cell_centres(pf%prob_lo(1), pf%prob_hi(1), nr, r)
      call cell_centres(pf%prob_lo(2), pf%prob_hi(2), nz, z)

      radius = -1.0_dp                       ! sentinel: not found
      do j = 1, nz
         ihit = 0
         do i = 1, nr - 1
            if (pf%data(i,j,1) == 0.0_dp) then
               radius(j) = r(i)
               exit
            end if
            if (pf%data(i,j,1) * pf%data(i+1,j,1) < 0.0_dp) then
               ihit = i
               exit
            end if
         end do
         if (ihit > 0) then
            denom = pf%data(ihit,j,1) - pf%data(ihit+1,j,1)
            if (denom /= 0.0_dp) then
               f = pf%data(ihit,j,1) / denom
            else
               f = 0.0_dp
            end if
            radius(j) = r(ihit) + f * (r(ihit+1) - r(ihit))
         end if
      end do

      call free_plotfile(pf)
   end subroutine interface_radius

   !> Amplitude of the `harmonic`-th Fourier component of r(z).
   !>
   !> r(z) is fitted by least squares with c0 + a cos(k z) + b sin(k z), k = 2 pi n / lambda, over the columns where an interface was found, and the amplitude is sqrt(a^2 + b^2). This makes the amplitude that of a cos(k z) perturbation, the dispersion-relation convention, independent of phase, and the constant c0 absorbs the mean radius so that a drifting mean (a mass-conservation error, say) cannot leak into the mode.
   !>
   !> When every column is valid and the box holds a whole number of wavelengths, the fit is exactly the former projection a = 2 |<(r - <r>) exp(-i k z)>|. When columns are missing it is not, and the difference matters: the projection treated the remaining columns as if they still formed a complete period, and one missing column out of 72 biased the amplitude by about 2 per cent. The fit has no such bias. n_missing, if requested, returns the number of columns without an interface, so that callers can flag a record built on incomplete profiles. Returns -1 with fewer than 4 valid columns or a singular fit.
   !>
   !> coef, if present, returns the signed coefficients (a, b) of cos(k z) and sin(k z). The magnitude alone cannot follow a crest turning into a trough, so it cannot separate the oscillation frequency of a stable mode (k r0 > 1) from its damping; the signed coefficients can (the old report's recommendation for the stable-band control).
   function mode_amplitude(z, radius, wavelength, harmonic, n_missing, coef) &
         result(a)
      real(dp), intent(in) :: z(:), radius(:), wavelength
      integer, intent(in) :: harmonic
      integer, intent(out), optional :: n_missing
      real(dp), intent(out), optional :: coef(2)
      real(dp) :: a
      integer :: j, n
      real(dp) :: ang, c, sn, m(3,3), rhs(3), sol(3)
      logical :: ok

      a = -1.0_dp
      if (present(coef)) coef = 0.0_dp
      n = count(radius >= 0.0_dp)
      if (present(n_missing)) n_missing = size(radius) - n
      if (n < 4) return                       ! not enough valid columns

      m = 0.0_dp; rhs = 0.0_dp
      do j = 1, size(radius)
         if (radius(j) < 0.0_dp) cycle
         ang = 2.0_dp * PI * real(harmonic,dp) * z(j) / wavelength
         c = cos(ang); sn = sin(ang)
         m(1,1) = m(1,1) + 1.0_dp
         m(1,2) = m(1,2) + c;      m(1,3) = m(1,3) + sn
         m(2,2) = m(2,2) + c*c;    m(2,3) = m(2,3) + c*sn
         m(3,3) = m(3,3) + sn*sn
         rhs(1) = rhs(1) + radius(j)
         rhs(2) = rhs(2) + radius(j) * c
         rhs(3) = rhs(3) + radius(j) * sn
      end do
      m(2,1) = m(1,2); m(3,1) = m(1,3); m(3,2) = m(2,3)
      call solve3(m, rhs, sol, ok)
      if (.not. ok) return
      a = sqrt(sol(2)**2 + sol(3)**2)
      if (present(coef)) coef = sol(2:3)
   end function mode_amplitude

   !> Solve the symmetric 3x3 system m x = b by Gaussian elimination with partial pivoting; ok is false if m is numerically singular.
   pure subroutine solve3(m, b, x, ok)
      real(dp), intent(in) :: m(3,3), b(3)
      real(dp), intent(out) :: x(3)
      logical, intent(out) :: ok
      real(dp) :: a(3,4), tmp(4), f, scale
      integer :: i, k, p

      a(:,1:3) = m
      a(:,4) = b
      scale = maxval(abs(m))
      ok = .false.
      x = 0.0_dp
      if (scale <= 0.0_dp) return
      do k = 1, 3
         p = k - 1 + maxloc(abs(a(k:3, k)), 1)
         if (abs(a(p, k)) <= 1.0e-12_dp * scale) return
         if (p /= k) then
            tmp = a(k,:); a(k,:) = a(p,:); a(p,:) = tmp
         end if
         do i = k + 1, 3
            f = a(i,k) / a(k,k)
            a(i,k:4) = a(i,k:4) - f * a(k,k:4)
         end do
      end do
      do i = 3, 1, -1
         x(i) = (a(i,4) - dot_product(a(i,i+1:3), x(i+1:3))) / a(i,i)
      end do
      ok = .true.
   end subroutine solve3

   !> Interface mode amplitude a(t) over every plotfile in a run.
   !>
   !> `meanr` is returned because a drifting mean radius is the usual signal of a mass-conservation problem, which would invalidate the growth rate.
   !>
   !> nmiss, if present, returns per snapshot the number of axial columns in which no interface was found; any non-zero entry means the amplitude was fitted to an incomplete profile and the record should be inspected.
   !>
   !> coefs, if present, returns the signed cos and sin coefficients per snapshot (2, n), and meanr2 the axial mean of r^2 over the valid columns. The column volume per unit length is pi * mean(r^2), so meanr2 is the geometric volume check of the old report: the mean radius itself must fall as a^2/(4 r0) when the amplitude grows at constant volume, whereas mean(r^2) must stay constant.
   subroutine amplitude_series(run_dir, wavelength, ls_field, harmonic, &
                               times, amps, meanr, n, ierr, prefix, nmiss, &
                               coefs, meanr2)
      character(len=*), intent(in) :: run_dir, ls_field
      real(dp), intent(in) :: wavelength
      integer, intent(in) :: harmonic
      real(dp), allocatable, intent(out) :: times(:), amps(:), meanr(:)
      integer, intent(out) :: n, ierr
      character(len=*), intent(in), optional :: prefix
      integer, allocatable, intent(out), optional :: nmiss(:)
      real(dp), allocatable, intent(out), optional :: coefs(:,:), meanr2(:)

      character(len=MAX_NAME), allocatable :: paths(:)
      real(dp), allocatable :: z(:), radius(:)
      integer, allocatable :: miss(:), order(:)
      real(dp), allocatable :: cf(:,:), m2(:)
      character(len=64) :: pfx
      integer :: i, nz, e, j, cnt
      real(dp) :: t, s

      pfx = 'nddataPLT'
      if (present(prefix)) pfx = prefix

      call find_plotfiles(run_dir, paths, n, trim(pfx))
      if (n == 0) then
         ierr = 1                             ! no plotfiles found
         return
      end if

      allocate(times(n), amps(n), meanr(n), miss(n), cf(2,n), m2(n))
      do i = 1, n
         call interface_radius(trim(paths(i)), ls_field, t, z, radius, nz, e)
         if (e /= 0) then
            ierr = 100 + e
            return
         end if
         times(i) = t
         amps(i) = mode_amplitude(z, radius, wavelength, harmonic, miss(i), &
                                  cf(:,i))
         s = 0.0_dp; cnt = 0
         do j = 1, nz
            if (radius(j) >= 0.0_dp) then
               s = s + radius(j); cnt = cnt + 1
            end if
         end do
         meanr(i) = merge(s/real(max(cnt,1),dp), -1.0_dp, cnt > 0)
         if (cnt > 0) then
            m2(i) = sum(radius**2, mask=radius >= 0.0_dp) / real(cnt, dp)
         else
            m2(i) = -1.0_dp
         end if
         deallocate(z, radius)
      end do

      ! Sort by TIME, not by step.  For twin pairing, step number is the right key (matching steps must be subtracted).  For a time series it is not: checkpoint restarts can emit a plotfile whose step number does not follow the time ordering of the surrounding files, so a step-sorted series can be non-monotonic in t.  Fits that select by time value are unaffected, but anything using "the last element" would be wrong.
      allocate(order(n))
      call argsort_real_asc(times, order)
      times = times(order); amps = amps(order)
      meanr = meanr(order); miss = miss(order)
      cf = cf(:, order); m2 = m2(order)
      if (present(nmiss)) call move_alloc(miss, nmiss)
      if (present(coefs)) call move_alloc(cf, coefs)
      if (present(meanr2)) call move_alloc(m2, meanr2)
      ierr = 0
   end subroutine amplitude_series

   !> Log-linear fit of a(t) over an optional window.
   !>
   !> Thin wrapper over util_mod::fit_log_linear, the single implementation shared with the field-norm fit, requiring at least 3 samples. t_min and t_max may be util_mod's UNSET.
   subroutine fit_growth_rate_amp(times, amps, t_min, t_max, sigma, r2, ierr)
      real(dp), intent(in) :: times(:), amps(:), t_min, t_max
      real(dp), intent(out) :: sigma, r2
      integer, intent(out) :: ierr
      real(dp) :: b

      call fit_log_linear(times, amps, t_min, t_max, 3, sigma, b, r2, ierr)
      if (ierr /= 0) ierr = 1
   end subroutine fit_growth_rate_amp

   !> Instantaneous growth rate d(ln a)/dt over sliding windows.
   !>
   !> This is the diagnostic that says whether the run has reached asymptotic modal growth.  If sigma_local is still CLIMBING at the end of the record, the flow is still in the startup transient and any single fitted growth rate underestimates the true one - which is exactly what the first real runs of this problem showed.
   !>
   !> Each window must hold at least 4 samples, as in the Python implementation (interface_mode.local_growth_rate); the Fortran version previously accepted 3, so the two could compare different windows.
   subroutine local_growth_rate(times, amps, width, centres, sigmas, nout)
      real(dp), intent(in) :: times(:), amps(:), width
      real(dp), allocatable, intent(out) :: centres(:), sigmas(:)
      integer, intent(out) :: nout

      integer, parameter :: NW = 8
      integer :: i, e
      real(dp) :: lo, hi, c, s, r2, b

      allocate(centres(NW), sigmas(NW))
      nout = 0
      lo = minval(times) + 0.5_dp*width
      hi = maxval(times) - 0.5_dp*width
      if (hi <= lo) return

      do i = 1, NW
         c = lo + (hi-lo)*real(i-1,dp)/real(NW-1,dp)
         call fit_log_linear(times, amps, c-0.5_dp*width, c+0.5_dp*width, &
                             4, s, b, r2, e)
         if (e == 0) then
            nout = nout + 1
            centres(nout) = c
            sigmas(nout) = s
         end if
      end do
   end subroutine local_growth_rate

   !> Heuristic: has the local growth rate stopped climbing?
   !>
   !> Compares the last two sliding-window estimates.  True is NECESSARY but not sufficient - also check grid convergence.
   function has_converged(times, amps, width, rel_tol, sigma_last, rel) &
         result(ok)
      real(dp), intent(in) :: times(:), amps(:), width, rel_tol
      real(dp), intent(out) :: sigma_last, rel
      logical :: ok
      real(dp), allocatable :: c(:), s(:)
      integer :: n

      call local_growth_rate(times, amps, width, c, s, n)
      if (n < 2) then
         ok = .false.; sigma_last = -1.0_dp; rel = -1.0_dp
         return
      end if
      sigma_last = s(n)
      rel = abs(s(n) - s(n-1)) / max(abs(s(n)), 1.0e-30_dp)
      ok = (rel <= rel_tol)
   end function has_converged

end module interface_mode_mod
