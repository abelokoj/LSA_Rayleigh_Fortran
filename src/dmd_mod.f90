!> SVD-based Dynamic Mode Decomposition (Schmid 2010).
!>
!> Fortran counterpart of plotting/dmd.py.  Implements the variant Ranjan, Unnikrishnan & Gaitonde (JCP 403, 2020) use in section 2.2 to extract stability modes from an NS-MFP snapshot subspace, with the amplitude ranking of Jovanovic, Schmid & Nichols (2014).
!>
!> Algorithm, for snapshots X (n_space x n_snap) sampled every dt: A = X(:,1:m-1),  B = X(:,2:m) A = U S V^H                                   (LAPACK dgesvd) S_tilde = U_r^H B V_r S_r^{-1}                (r = truncation rank) S_tilde w = mu w                              (LAPACK zgeev) Phi = B V_r S_r^{-1} w                        (exact DMD modes) omega = log(mu)/dt Re(omega) is the growth rate, Im(omega) the circular frequency.
!>
!> For Rayleigh-Plateau the mode of interest is STATIONARY (Im(omega) = 0): the capillary instability grows monotonically rather than oscillating, hence leading_growth_rate's stationary_only default.
!> ORIENTATION
!>
!> Dynamic Mode Decomposition: the step that turns a sequence of flow fields into growth rates.
!>
!> The idea in one sentence: assume each snapshot follows from the previous one by multiplication by a fixed matrix A, then obtain the eigenvalues of A without ever forming it. If an eigenvalue lambda satisfies |lambda| > 1 the corresponding mode grows from snapshot to snapshot, and sigma = log(lambda)/dt converts that per-snapshot factor into a rate per unit time.
!>
!> The obstacle is size. Each snapshot holds thousands of values, so A would hold millions, and forming it is neither possible nor necessary. The projection below reduces the problem to the few directions the snapshots actually span, found by a singular value decomposition, and takes the eigenvalues of the resulting small matrix. Those are the eigenvalues sought; the large matrix never exists. This is the sense in which the method is Jacobian-free.
module dmd_mod
   use table_mod
   use util_mod, only: argsort_real_desc, nan_value
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private

   public :: dp
   public :: dmd_result_t, dmd_compute, dmd_free
   public :: leading_growth_rate, dominant_growth_rate, dmd_summary

   !> Outcome of a DMD, sorted by descending amplitude.
   type :: dmd_result_t
      integer :: n_modes = 0
      integer :: rank = 0
      real(dp) :: dt = 0.0_dp
      complex(dp), allocatable :: ritz(:)       !< discrete eigenvalues mu
      complex(dp), allocatable :: omega(:)      !< log(mu)/dt
      complex(dp), allocatable :: modes(:,:)    !< (n_space, n_modes)
      real(dp),    allocatable :: amplitude(:)  !< |b|
      real(dp),    allocatable :: svals(:)      !< singular values of A
      integer :: rank_requested = 0             !< rank asked for (<= 0: by tolerance)
      integer :: amp_snapshot = 1               !< snapshot the amplitudes were fitted to
      !> Relative residual of the reduced linear map, ||Y - U_r S_tilde U_r^T X||_F / ||Y||_F. Small means the retained subspace evolves linearly from snapshot to snapshot; a large value means the data are not a linear sequence in the retained directions (transient, nonlinearity or too small a rank), which no eigenvalue printed alongside it can reveal.
      real(dp) :: residual = -1.0_dp
   end type dmd_result_t

contains

   subroutine dmd_free(res)
      type(dmd_result_t), intent(inout) :: res
      if (allocated(res%ritz))      deallocate(res%ritz)
      if (allocated(res%omega))     deallocate(res%omega)
      if (allocated(res%modes))     deallocate(res%modes)
      if (allocated(res%amplitude)) deallocate(res%amplitude)
      if (allocated(res%svals))     deallocate(res%svals)
      res%n_modes = 0
      res%rank = 0
   end subroutine dmd_free

   !> Compute the DMD of a snapshot matrix.
   !>
   !> @param X      (n_space, n_snap) snapshots; column j is time (j-1)*dt @param dt     sampling interval @param rank   truncation rank; <= 0 selects by tolerance @param tol    relative singular-value cutoff when rank <= 0 @param res    result @param ierr   0 on success, non-zero on failure (5: the snapshots are identically zero)
   !>
   !> An explicit rank is capped at the numerical rank of the data, the number of singular values above max(n, m) * epsilon * S(1). Asking for more previously divided by a zero singular value, filled S_tilde with NaN and made LAPACK stop the whole program; res%rank reports the rank actually used and res%rank_requested the one asked for.
   !>
   !> Amplitudes are the least-squares coefficients of the modes in one snapshot. Ranjan et al. use the first snapshot, and so does this routine whenever that snapshot carries signal. With velocity fields and skip = 0 it does not: both twins start from rest, Q'(0) = 0 exactly, every amplitude came out zero and the ranking was arbitrary. When the first snapshot is negligible (norm below 1e-8 of the largest), the snapshot of largest norm is used instead, and the amplitudes are the contributions of the modes to that snapshot. They are deliberately not referred back to t = 0 by dividing by mu^(k-1): when Q'(0) = 0 the t = 0 coefficients of the growing and decaying modes cancel one another, and the division would inflate the strongly decaying ones until they led the ranking. res%amp_snapshot records which snapshot was used, and the driver warns.
   subroutine dmd_compute(X, dt, rank, tol, res, ierr)
      real(dp), intent(in) :: X(:,:)
      real(dp), intent(in) :: dt, tol
      integer,  intent(in) :: rank
      type(dmd_result_t), intent(out) :: res
      integer, intent(out) :: ierr

      integer :: n, m, mm, r, i, j, info, lwork, nnz, kref
      real(dp) :: cnorm, cmax
      real(dp), allocatable :: A(:,:), B(:,:), U(:,:), VT(:,:), S(:)
      real(dp), allocatable :: rwork(:), work(:)
      complex(dp), allocatable :: St(:,:), w(:,:), vl(:,:), cwork(:)
      complex(dp), allocatable :: BVS(:,:), rhs(:), lhs(:,:)
      complex(dp), allocatable :: mu(:)
      real(dp), allocatable :: BVSr(:,:)
      integer, allocatable :: order(:)

      ierr = 0
      n = size(X, 1)
      m = size(X, 2)

      if (m < 2) then
         ierr = 1; return            ! need at least two snapshots
      end if
      if (dt <= 0.0_dp) then
         ierr = 2; return
      end if
      if (n < 1) then
         ierr = 3; return
      end if
      if (any(X /= X)) then          ! NaN check
         ierr = 4; return
      end if
      res%rank_requested = rank

      mm = m - 1
      allocate(A(n, mm), B(n, mm))
      A = X(:, 1:mm)
      B = X(:, 2:m)

      ! ---- SVD of A:  A = U S V^T -------------------------------------
      allocate(U(n, min(n, mm)), VT(min(n, mm), mm), S(min(n, mm)))
      lwork = max(1, 10 * max(n, mm))
      allocate(work(lwork))
      call dgesvd('S', 'S', n, mm, A, n, S, U, n, VT, min(n, mm), &
                  work, lwork, info)
      if (info /= 0) then
         ierr = 10 + info; return
      end if
      deallocate(work)
      if (.not. (S(1) > 0.0_dp)) then
         ierr = 5; return            ! all snapshots zero: nothing to decompose
      end if

      ! ---- truncation rank --------------------------------------------
      ! nnz is the numerical rank; no singular value at or below it may be inverted below.
      nnz = count(S > real(max(n, mm), dp) * epsilon(1.0_dp) * S(1))
      nnz = max(nnz, 1)
      if (rank > 0) then
         r = min(rank, size(S), nnz)
      else
         r = 0
         do i = 1, size(S)
            if (S(i) > tol * S(1)) r = r + 1
         end do
         r = min(max(r, 1), nnz)
      end if

      ! ---- BVS = B V_r S_r^{-1}  (real) -------------------------------- V = VT^T, so V(:,i) = VT(i,:).
      allocate(BVSr(n, r))
      BVSr = 0.0_dp
      do i = 1, r
         do j = 1, mm
            BVSr(:, i) = BVSr(:, i) + B(:, j) * VT(i, j)
         end do
         BVSr(:, i) = BVSr(:, i) / S(i)
      end do

      ! ---- S_tilde = U_r^T * BVS  (r x r) ------------------------------
      allocate(St(r, r))
      do i = 1, r
         do j = 1, r
            St(i, j) = cmplx(dot_product(U(:, i), BVSr(:, j)), 0.0_dp, dp)
         end do
      end do

      ! ---- residual of the reduced linear map ---------------------------
      ! U_r^T A = S_r V_r^T, so U_r S_tilde U_r^T A = U_r (S_tilde (S_r V_r^T)).
      block
         real(dp), allocatable :: svt(:,:), str(:,:), pred(:,:)
         real(dp) :: ynorm
         allocate(svt(r, mm), str(r, r))
         do i = 1, r
            svt(i, :) = S(i) * VT(i, :)
         end do
         str = real(St, dp)
         pred = matmul(U(:, 1:r), matmul(str, svt))
         ynorm = sqrt(sum(B**2))
         if (ynorm > 0.0_dp) then
            res%residual = sqrt(sum((B - pred)**2)) / ynorm
         else
            res%residual = 0.0_dp
         end if
      end block

      ! ---- eigen-decomposition of S_tilde ------------------------------
      allocate(mu(r), w(r, r), vl(1, 1))
      lwork = max(1, 4 * r)
      allocate(cwork(lwork), rwork(max(1, 2 * r)))
      call zgeev('N', 'V', r, St, r, mu, vl, 1, w, r, cwork, lwork, rwork, info)
      if (info /= 0) then
         ierr = 20 + info; return
      end if
      deallocate(cwork, rwork)

      ! ---- exact DMD modes:  Phi = BVS * w ------------------------------
      allocate(BVS(n, r))
      BVS = cmplx(BVSr, 0.0_dp, dp)
      allocate(res%modes(n, r))
      res%modes = matmul(BVS, w)

      ! ---- amplitudes: least-squares fit of the modes to one snapshot ---
      ! The first snapshot unless it is negligible; see the note above.
      kref = 1
      cmax = 0.0_dp
      do j = 1, m
         cnorm = sqrt(sum(X(:, j)**2))
         if (cnorm > cmax) then
            cmax = cnorm; kref = j
         end if
      end do
      if (sqrt(sum(X(:, 1)**2)) > 1.0e-8_dp * cmax) kref = 1
      res%amp_snapshot = kref

      allocate(lhs(n, r), rhs(max(n, r)))
      lhs = res%modes
      rhs = cmplx(0.0_dp, 0.0_dp, dp)
      rhs(1:n) = cmplx(X(:, kref), 0.0_dp, dp)
      allocate(res%amplitude(r))
      call solve_lstsq(lhs, rhs, n, r, res%amplitude, info)
      if (info /= 0) then
         ierr = 30 + info; return
      end if

      ! ---- continuous-time eigenvalues ---------------------------------
      allocate(res%ritz(r), res%omega(r))
      res%ritz = mu
      do i = 1, r
         if (abs(mu(i)) > 0.0_dp) then
            res%omega(i) = log(mu(i)) / dt
         else
            ! a mode annihilated in one step: -inf growth
            res%omega(i) = cmplx(-huge(1.0_dp), 0.0_dp, dp)
         end if
      end do

      res%n_modes = r
      res%rank = r
      res%dt = dt
      allocate(res%svals(size(S)))
      res%svals = S

      ! ---- sort by descending amplitude --------------------------------
      allocate(order(r))
      call argsort_real_desc(res%amplitude, order)
      res%ritz      = res%ritz(order)
      res%omega     = res%omega(order)
      res%amplitude = res%amplitude(order)
      res%modes     = res%modes(:, order)
   end subroutine dmd_compute

   !> Least-squares solve of lhs * b = rhs, returning |b|.
   subroutine solve_lstsq(lhs, rhs, n, r, absb, info)
      complex(dp), intent(inout) :: lhs(:,:), rhs(:)
      integer, intent(in) :: n, r
      real(dp), intent(out) :: absb(:)
      integer, intent(out) :: info

      complex(dp), allocatable :: cwork(:)
      real(dp), allocatable :: rwork(:), sv(:)
      integer :: lwork, rnk, i

      allocate(sv(min(n, r)), rwork(max(1, 5 * min(n, r))))
      lwork = max(1, 2 * min(n, r) + max(n, r) + 64 * max(n, r))
      allocate(cwork(lwork))
      call zgelss(n, r, 1, lhs, size(lhs, 1), rhs, size(rhs), sv, &
                  -1.0_dp, rnk, cwork, lwork, rwork, info)
      if (info == 0) then
         do i = 1, r
            absb(i) = abs(rhs(i))
         end do
      end if
   end subroutine solve_lstsq

   !> Largest growth rate among physically retained modes.
   !>
   !> For Rayleigh-Plateau use stationary_only = .true. (the default): the capillary mode is non-oscillatory, and restricting to it avoids picking up weak spurious oscillatory pairs. A mode counts as stationary when |Im omega| <= max(freq_tol, freq_rel * |Re omega|). The relative part matters with real data: noise gives the physical mode a small imaginary part, and with the former purely absolute test of 1e-6 the mode was discarded, nothing survived, and the routine returned -huge, which went into the tables and the most-unstable-k search while the program still exited 0. When no mode survives the filters the result is now NaN, which callers test with util_mod::is_nan and report as a failure.
   function leading_growth_rate(res, stationary_only, freq_tol, amp_frac, &
                                freq_rel) result(s)
      type(dmd_result_t), intent(in) :: res
      logical,  intent(in), optional :: stationary_only
      real(dp), intent(in), optional :: freq_tol, amp_frac, freq_rel
      real(dp) :: s
      logical  :: stat_only, found
      real(dp) :: ftol, frel, afrac, amax, g, f
      integer  :: i

      stat_only = .true.
      if (present(stationary_only)) stat_only = stationary_only
      ftol = 1.0e-6_dp
      if (present(freq_tol)) ftol = freq_tol
      frel = 1.0e-3_dp
      if (present(freq_rel)) frel = freq_rel
      afrac = 1.0e-3_dp
      if (present(amp_frac)) afrac = amp_frac

      s = nan_value()
      found = .false.
      if (res%n_modes < 1) return
      amax = maxval(res%amplitude)

      do i = 1, res%n_modes
         if (res%amplitude(i) < afrac * amax) cycle
         f = aimag(res%omega(i))
         g = real(res%omega(i), dp)
         if (stat_only .and. abs(f) > max(ftol, frel * abs(g))) cycle
         if (.not. found .or. g > s) s = g
         found = .true.
      end do
   end function leading_growth_rate

   !> Growth rate of the DOMINANT growing, near-stationary mode: among the modes with Re(omega) > 0 and |Im(omega)| <= freq_rel |Re(omega)| (default 0.1), the one with the largest amplitude.
   !>
   !> Added 2026-10-08 after the first campaign data. leading_growth_rate returns the LARGEST growth rate among modes that pass a weak amplitude filter; with a rank close to the number of snapshots, which the default tolerance gives, small spurious modes then won (meas08 at rank 49: +194 %). The physical capillary mode is the one that carries the record, so the amplitude, not the growth rate, is the selection criterion here. The relative frequency tolerance is wider than the stationary filter of leading_growth_rate because two close real eigenvalues that DMD cannot separate come out as a complex pair with a small imaginary part; their real part still estimates the rate. NaN if no mode qualifies.
   function dominant_growth_rate(res, freq_rel) result(s)
      type(dmd_result_t), intent(in) :: res
      real(dp), intent(in), optional :: freq_rel
      real(dp) :: s, frel, best, g
      integer :: i
      frel = 0.1_dp
      if (present(freq_rel)) frel = freq_rel
      s = nan_value()
      best = -1.0_dp
      do i = 1, res%n_modes
         g = real(res%omega(i), dp)
         if (.not. (g > 0.0_dp)) cycle
         if (abs(aimag(res%omega(i))) > frel * g) cycle
         if (res%amplitude(i) > best) then
            best = res%amplitude(i); s = g
         end if
      end do
   end function dominant_growth_rate

   !> Print a table of the leading modes. The frequency column is Im(omega)/(2 pi), in cycles per unit time, not the angular frequency.
   subroutine dmd_summary(res, nshow, unit)
      type(dmd_result_t), intent(in) :: res
      integer, intent(in), optional  :: nshow, unit
      integer :: n, u, i
      type(table_t) :: t

      n = min(10, res%n_modes)
      if (present(nshow)) n = min(nshow, res%n_modes)
      u = 6
      if (present(unit)) u = unit

      if (u == 6) then
         call t%init([cel('mode '), cel('growth sigma'), cel('freq Im(w)/2pi'), &
            cel('|lambda| '), cel('amplitude ')], 'rrrrr')
         do i = 1, n
            call t%row([cel(fmt_i(i-1)), cel(fmt_f(real(res%omega(i), dp), 6)), &
               cel(fmt_f(aimag(res%omega(i))/(2.0_dp*3.141592653589793_dp), 6)), &
               cel(fmt_f(abs(res%ritz(i)), 6)), cel(fmt_e(res%amplitude(i)))])
         end do
         call t%emit('Decomposition: rank '//fmt_i(res%rank)//', dt '// &
                     fmt_e(res%dt)//', '//fmt_i(res%n_modes)//' modes')
      else
         write(u, '(a,i0,a,es12.5,a,i0)') 'rank=', res%rank, '  dt=', res%dt, &
            '  n_modes=', res%n_modes
         write(u, '(a4,a16,a16,a12,a14)') '#', 'growth sigma', 'freq Im(w)/2pi', &
            '|lambda|', 'amplitude'
         do i = 1, n
            write(u, '(i4,es16.6,es16.6,f12.6,es14.4)') i - 1, &
               real(res%omega(i), dp), aimag(res%omega(i)) / &
               (2.0_dp * 3.141592653589793_dp), &
               abs(res%ritz(i)), res%amplitude(i)
         end do
      end if
   end subroutine dmd_summary

end module dmd_mod
