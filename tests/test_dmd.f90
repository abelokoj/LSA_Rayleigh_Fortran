!> Tests for dmd_mod.
!>
!> Synthetic signals with analytically known eigenvalues, so a failure here means the DMD is wrong rather than the CFD being hard.  The last test is the one that matters most: a Rayleigh-Plateau-like signal (a stationary, exponentially growing interface mode buried in decaying contamination) must return the exact dispersion-relation growth rate.
program test_dmd
   use testing_mod
   use dmd_mod
   use dispersion_mod, only: growth_rate
   implicit none

   integer, parameter :: NS = 200
   real(dp), allocatable :: X(:,:)
   type(dmd_result_t) :: res
   real(dp) :: dt, rate, freq, s, sexact
   integer :: ierr, i, j, nsnap
   real(dp) :: shp(NS), shp2(NS), shp3(NS)
   real(dp) :: t, nan_zero

   call init_shapes()

   call section('single stationary growing mode')
   dt = 0.05_dp; rate = 0.37_dp; nsnap = 60
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      X(:, j) = shp * exp(rate * real(j-1, dp) * dt)
   end do
   call dmd_compute(X, dt, 0, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 0, 'dmd_compute succeeds')
   call check_close_rel(leading_growth_rate(res), rate, 1.0e-8_dp, &
                        'recovers the planted growth rate')
   call check_close(aimag(res%omega(1)), 0.0_dp, 1.0e-9_dp, &
                    'leading mode is stationary')
   call dmd_free(res); deallocate(X)

   call section('decaying mode: sign convention')
   dt = 0.05_dp; rate = -0.42_dp; nsnap = 60
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      X(:, j) = shp * exp(rate * real(j-1, dp) * dt)
   end do
   call check(norm2(X(:, nsnap)) < norm2(X(:, 1)), &
              'planted decay really decays')
   call dmd_compute(X, dt, 0, 1.0e-10_dp, res, ierr)
   call check_close_rel(leading_growth_rate(res), rate, 1.0e-8_dp, &
                        'recovers a negative growth rate')
   call check(leading_growth_rate(res) < 0.0_dp, &
              'negative sigma means decay')
   call dmd_free(res); deallocate(X)

   call section('oscillatory conjugate pair')
   ! A real oscillation is a conjugate pair; built from one complex shape and its conjugate it spans a genuine 2-D subspace (a cos - b sin).  Two independent real shapes would collapse to rank 1 and carry no frequency.
   dt = 0.02_dp; rate = -0.15_dp; freq = 1.3_dp; nsnap = 200
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      t = real(j-1, dp) * dt
      X(:, j) = exp(rate*t) * (shp*cos(2.0_dp*3.141592653589793_dp*freq*t) &
                             - shp2*sin(2.0_dp*3.141592653589793_dp*freq*t))
   end do
   call dmd_compute(X, dt, 0, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 0, 'dmd_compute succeeds')
   call check_close_rel(abs(aimag(res%omega(1)))/(2.0_dp*3.141592653589793_dp), &
                        freq, 1.0e-6_dp, 'recovers the frequency')
   call check_close_rel(real(res%omega(1), dp), rate, 1.0e-6_dp, &
                        'recovers the decay rate of the pair')
   call dmd_free(res); deallocate(X)

   call section('mixed spectrum: growing stationary + decaying oscillatory')
   dt = 0.01_dp; nsnap = 400
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      t = real(j-1, dp) * dt
      X(:, j) = shp * exp(0.42_dp*t) &
              + exp(-0.30_dp*t) * &
                 (shp2*cos(2.0_dp*3.141592653589793_dp*2.0_dp*t) &
                                 - shp3*sin(2.0_dp*3.141592653589793_dp*2.0_dp*t))
   end do
   call dmd_compute(X, dt, 0, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 0, 'dmd_compute succeeds')
   call check_close_rel(leading_growth_rate(res), 0.42_dp, 1.0e-6_dp, &
                        'stationary filter isolates the unstable mode')
   call check(res%rank >= 3, 'rank at least 3 for a 3-mode signal')
   call dmd_free(res); deallocate(X)

   call section('rank truncation')
   dt = 0.05_dp; nsnap = 100
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      t = real(j-1, dp) * dt
      X(:, j) = shp*exp(0.2_dp*t) + shp2*exp(-0.4_dp*t) + shp3*exp(-0.9_dp*t)
   end do
   call dmd_compute(X, dt, 2, 1.0e-10_dp, res, ierr)
   call check_int(res%rank, 2, 'explicit rank is honoured')
   call dmd_free(res)
   call dmd_compute(X, dt, 0, 1.0e-8_dp, res, ierr)
   call check_int(res%rank, 3, 'tolerance finds the true rank 3')
   call dmd_free(res); deallocate(X)

   call section('input validation')
   allocate(X(NS, 1))
   X = 0.0_dp
   call dmd_compute(X, 0.05_dp, 0, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 1, 'rejects fewer than 2 snapshots')
   deallocate(X)
   allocate(X(NS, 5)); X = 0.0_dp
   call dmd_compute(X, -1.0_dp, 0, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 2, 'rejects non-positive dt')
   X = 0.0_dp
   ! Build a NaN at run time; sqrt(-1) as a literal is a compile-time error.
   nan_zero = 0.0_dp
   X(1,1) = nan_zero / nan_zero
   call dmd_compute(X, 0.05_dp, 0, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 4, 'rejects non-finite snapshots')
   deallocate(X)

   call section('rank above the numerical rank (former crash)')
   ! Rank-1 data with an explicit rank of 3 used to divide by a zero singular value; reference LAPACK then stopped the program from inside zgeev.
   dt = 0.05_dp; nsnap = 40
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      X(:, j) = shp * exp(0.25_dp * real(j-1, dp) * dt)
   end do
   call dmd_compute(X, dt, 3, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 0, 'explicit rank 3 on rank-1 data does not fail')
   call check_int(res%rank, 1, 'rank is capped at the numerical rank')
   call check_int(res%rank_requested, 3, 'requested rank is reported')
   call check_close_rel(leading_growth_rate(res), 0.25_dp, 1.0e-8_dp, &
                        'and the rate is still recovered')
   call dmd_free(res)
   X = 0.0_dp
   call dmd_compute(X, dt, 0, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 5, 'all-zero snapshots are rejected, not divided by')
   call dmd_free(res); deallocate(X)

   call section('zero first snapshot (velocity fields with --skip 0)')
   ! Q'(0) = 0 exactly, as for twins that both start from rest, followed by a growing mode. Ranking by the first snapshot gave every amplitude zero, so the amplitude filter kept every mode. (With a single mode the zero column does not disturb the eigenvalue, so the planted rate is recovered exactly and the test isolates the ranking.)
   dt = 0.05_dp; nsnap = 80
   allocate(X(NS, nsnap))
   X(:, 1) = 0.0_dp
   do j = 2, nsnap
      t = real(j-1, dp) * dt
      X(:, j) = shp * exp(0.3_dp*t)
   end do
   call dmd_compute(X, dt, 0, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 0, 'dmd_compute succeeds')
   call check(res%amp_snapshot > 1, 'amplitudes taken from a later snapshot')
   call check(maxval(res%amplitude) > 0.0_dp, 'amplitudes are not all zero')
   call check_close_rel(leading_growth_rate(res), 0.3_dp, 1.0e-6_dp, &
                        'the growing mode is still found')
   call dmd_free(res); deallocate(X)

   call section('nothing survives the filters')
   ! A purely oscillatory pair and no stationary mode: the former code returned -huge, which then entered the tables.
   dt = 0.05_dp; nsnap = 60
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      t = real(j-1, dp) * dt
      X(:, j) = exp(0.1_dp*t) * (shp*cos(3.0_dp*t) + shp2*sin(3.0_dp*t))
   end do
   call dmd_compute(X, dt, 0, 1.0e-10_dp, res, ierr)
   s = leading_growth_rate(res)
   call check(s /= s, 'returns NaN when no stationary mode exists')
   call check_close_rel(leading_growth_rate(res, stationary_only=.false.), &
                        0.1_dp, 1.0e-6_dp, 'the oscillatory rate is available on request')
   call dmd_free(res); deallocate(X)

   call section('stationary test is relative as well as absolute')
   ! A real growing mode whose eigenvalue carries a tiny imaginary part, as noise gives it in practice: |Im| = 5e-5 is above the absolute 1e-6 but far below 1e-3 |Re|.
   dt = 0.05_dp; nsnap = 60
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      t = real(j-1, dp) * dt
      X(:, j) = exp(0.3_dp*t) * (shp*cos(5.0e-5_dp*t) + shp2*sin(5.0e-5_dp*t))
   end do
   call dmd_compute(X, dt, 0, 1.0e-12_dp, res, ierr)
   call check_close_rel(leading_growth_rate(res), 0.3_dp, 1.0e-6_dp, &
                        'a nearly stationary mode is kept')
   call dmd_free(res); deallocate(X)

   call section('dominant versus fastest mode')
   ! A large slow mode and a tiny fast one: leading_growth_rate (fastest that passes the 1e-3 amplitude filter) returns the tiny one, dominant_growth_rate the one that carries the record.
   dt = 0.05_dp; nsnap = 60
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      t = real(j-1, dp) * dt
      X(:, j) = shp*exp(0.3_dp*t) + 2.0e-3_dp*shp2*exp(0.9_dp*t)
   end do
   call dmd_compute(X, dt, 0, 1.0e-12_dp, res, ierr)
   call check_close_rel(dominant_growth_rate(res), 0.3_dp, 1.0e-6_dp, &
                        'dominant mode is the large slow one')
   call check_close_rel(leading_growth_rate(res), 0.9_dp, 1.0e-6_dp, &
                        'fastest-mode rule picks the tiny fast one')
   call dmd_free(res); deallocate(X)

   call section('synthetic Rayleigh-Plateau signal')
   ! Plant the exact inviscid growth rate for k r0 = 0.7 in a stationary sinusoidal interface mode, plus faster-decaying contamination.
   sexact = growth_rate(0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp)
   dt = 0.05_dp; nsnap = 300
   allocate(X(NS, nsnap))
   do j = 1, nsnap
      t = real(j-1, dp) * dt
      do i = 1, NS
         X(i, j) = 1.0e-5_dp * cos(0.7_dp*real(i,dp)*0.1_dp) * exp(sexact*t) &
                 + 5.0e-6_dp * cos(2.1_dp*real(i,dp)*0.1_dp) * exp(-2.0_dp*t)
      end do
   end do
   call dmd_compute(X, dt, 0, 1.0e-10_dp, res, ierr)
   s = leading_growth_rate(res)
   call check_int(ierr, 0, 'dmd_compute succeeds')
   call check_close_rel(s, sexact, 1.0e-6_dp, &
                        'recovers the exact Rayleigh growth rate')
   call check(s > 0.2_dp .and. s < 0.4_dp, &
              'growth rate has the expected magnitude ~0.34')
   call dmd_free(res); deallocate(X)

   call summary('test_dmd')

contains

   subroutine init_shapes()
      integer :: i
      ! Deterministic, mutually independent spatial shapes.
      do i = 1, NS
         shp(i)  = sin(0.11_dp*real(i,dp)) + 0.3_dp*cos(0.023_dp*real(i,dp))
         shp2(i) = cos(0.07_dp*real(i,dp)) - 0.2_dp*sin(0.041_dp*real(i,dp))
         shp3(i) = sin(0.031_dp*real(i,dp)+0.5_dp)
      end do
   end subroutine init_shapes

end program test_dmd
