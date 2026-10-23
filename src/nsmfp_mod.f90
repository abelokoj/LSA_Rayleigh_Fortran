!> NS-MFP snapshot assembly for the Rayleigh-Plateau instability.
!>
!> Fortran counterpart of plotting/nsmfp.py.  Implements the subspace-generation half of Ranjan, Unnikrishnan & Gaitonde (JCP 403, 2020) on top of the amrex_implicit_interfaces two-phase solver, WITHOUT modifying the solver.
!>
!> The B_f body force, and why the solver requires no modification
!! -----------------------------------------------------
!! Ranjan et al. add a restraining body force
!!
!! B_f = dQbar/dt - F(Qbar)                              (their Eq. 3)
!!
!! whose sole role, visible in their Eq. (6), is to annihilate the bracketed base-state drift term.  We achieve the same cancellation by running an unperturbed TWIN alongside the perturbed case from an identical restart, with identical discretisation, and subtracting snapshot by snapshot:
!!
!! Q'(t_n) = Q_perturbed(t_n) - Q_unperturbed(t_n).
!!
!! Both runs suffer the same drift, so the difference removes it to machine precision and Q' obeys the linearised dynamics to O(eps^2), exactly as in their Eq. (7)-(8).  This is the "subtract the base flow from the total flow" step of their Fig. 1; the drift that B_f cancels by stored forcing is here cancelled by a computed twin.
!!
!! For this problem the substitution is more than a convenience.  A quiescent liquid column is an EXACT equilibrium, so the only base-state drift is discretisation error - principally the parasitic (spurious) currents of the discrete curvature/surface-tension operator.  Those are precisely what the twin measures and removes, and unlike a frozen B_f the twin tracks their evolution.
!!
!! Validity conditions:
!! - identical grid, domain, field set and times (ns.fixed_dt, max_level=0); enforced by plotfile_mod::compatible
!! - same MPI rank count for both runs; NOT checked here, because a plotfile does not record the rank count. It is the user's responsibility: run_all.sh records mpi_ranks in each case's manifest.txt and warns on a mismatch for --twin-of.
module nsmfp_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use plotfile_mod
   use util_mod, only: fit_log_linear, cell_centres
   implicit none
   private

   public :: dp
   public :: build_snapshot_matrix, perturbation_norms
   public :: fit_exponential_growth, uniform_dt, build_krylov_matrix
   public :: base_drift

contains

   !> Assemble the NS-MFP subspace Q'(t_n) = Q_pert(t_n) - Q_unpert(t_n).
   !>
   !> Snapshots are paired strictly by STEP NUMBER, so a plotfile missing on one side cannot silently shift the pairing.
   !>
   !> @param pert_dir, base_dir  directories holding the twin runs' plt* output @param fields   component names forming the state vector @param skip     discard this many leading snapshots (shed the transient) @param stride   use every stride-th snapshot @param maxsnap  cap on retained snapshots (<= 0 for no cap) @param X        (n_space, n_snap) subspace, allocated here @param times    snapshot times @param dt       uniform sampling interval @param ierr     0 on success @param weight   'none' (default) or 'volume': with 'volume' every row is multiplied by sqrt(r) at its cell centre, so that the Euclidean inner product of the snapshot matrix becomes the axisymmetric volume integral (the factor 2 pi dr dz is constant and dropped). The unweighted norm counts a cell near the axis as much as one at the outer wall, although it holds far less fluid; the old report noted that it therefore should not be called an energy, and dissertation section 4.4.3 makes the choice of inner product explicit. Exact DMD eigenvalues are unchanged by an invertible row scaling; what changes is which directions the rank truncation keeps. @param grid  if present, returns nx, ny and ncomp of the plotfiles, so that a mode vector can be reshaped into fields @param domain  if present, returns prob_lo(1:2) and prob_hi(1:2)
   subroutine build_snapshot_matrix(pert_dir, base_dir, fields, skip, stride, &
                                    maxsnap, X, times, dt, ierr, verbose, &
                                    prefix, weight, grid, domain)
      character(len=*), intent(in) :: pert_dir, base_dir
      character(len=*), intent(in) :: fields(:)
      integer, intent(in) :: skip, stride, maxsnap
      real(dp), allocatable, intent(out) :: X(:,:), times(:)
      real(dp), intent(out) :: dt
      integer, intent(out) :: ierr
      logical, intent(in), optional :: verbose
      !> Plotfile directory prefix.  The solver does NOT honour amr.plot_file: with ns.visual_nddata_format=1 it writes nddataPLT00000023 (and MOF_PLT*), not plt00023.  Confirmed by running the solver.  Default therefore is 'nddataPLT'.
      character(len=*), intent(in), optional :: prefix
      character(len=*), intent(in), optional :: weight
      integer, intent(out), optional :: grid(3)
      real(dp), intent(out), optional :: domain(4)

      character(len=MAX_NAME), allocatable :: ppaths(:), bpaths(:)
      real(dp), allocatable :: rc(:), w(:)
      logical :: weighted
      integer :: np, nb, i, j, k, n_pair, n_use, nsp, idx
      integer, allocatable :: pidx(:), bidx(:)
      type(plotfile_t) :: pa, pb
      character(len=256) :: msg
      logical :: verb
      real(dp) :: nrm
      character(len=64) :: pfx

      verb = .false.
      if (present(verbose)) verb = verbose
      weighted = .false.
      if (present(weight)) then
         select case (trim(weight))
         case ('', 'none'); weighted = .false.
         case ('volume');   weighted = .true.
         case default
            ierr = 7; return                 ! unknown weight
         end select
      end if
      pfx = 'nddataPLT'
      if (present(prefix)) pfx = prefix
      ierr = 0
      dt = 0.0_dp
      if (skip < 0 .or. stride < 1) then
         ierr = 6; return                    ! skip must be >= 0, stride >= 1
      end if

      call find_plotfiles(pert_dir, ppaths, np, trim(pfx))
      call find_plotfiles(base_dir, bpaths, nb, trim(pfx))
      if (np == 0) then
         ierr = 1; return                    ! no nddataPLT* under perturbed dir
      end if
      if (nb == 0) then
         ierr = 2; return                    ! no nddataPLT* under base dir
      end if

      ! ---- pair by step number ------------------------------------------
      allocate(pidx(np), bidx(np))
      n_pair = 0
      do i = 1, np
         do j = 1, nb
            if (plotfile_step(ppaths(i)) == plotfile_step(bpaths(j))) then
               n_pair = n_pair + 1
               pidx(n_pair) = i
               bidx(n_pair) = j
               exit
            end if
         end do
      end do
      if (n_pair < 2) then
         ierr = 3; return                    ! need >= 2 matching steps
      end if

      ! ---- apply skip / stride / cap ------------------------------------
      n_use = 0
      do i = skip + 1, n_pair, max(stride, 1)
         n_use = n_use + 1
         if (maxsnap > 0 .and. n_use >= maxsnap) exit
      end do
      if (n_use < 2) then
         ierr = 4; return                    ! too few after skip/stride
      end if

      ! ---- read and subtract ---------------------------------------------
      allocate(times(n_use))
      k = 0
      do i = skip + 1, n_pair, max(stride, 1)
         k = k + 1
         if (k > n_use) exit

         call read_plotfile(trim(ppaths(pidx(i))), pa, ierr, fields)
         if (ierr /= 0) then
            ierr = 1000 + ierr; return
         end if
         call read_plotfile(trim(bpaths(bidx(i))), pb, ierr, fields)
         if (ierr /= 0) then
            ierr = 2000 + ierr; return
         end if
         if (.not. compatible(pa, pb, 1.0e-9_dp, msg)) then
            if (verb) write(*,'(a)') 'incompatible plotfiles: '//trim(msg)
            ierr = 5; return
         end if

         nsp = pa%nx * pa%ny * pa%ncomp
         if (.not. allocated(X)) then
            allocate(X(nsp, n_use))
            if (present(grid)) grid = [pa%nx, pa%ny, pa%ncomp]
            if (present(domain)) domain = [pa%prob_lo, pa%prob_hi]
            if (weighted) then
               ! sqrt(r) per row, r constant along each z column of a component
               allocate(rc(pa%nx), w(pa%nx * pa%ny))
               call cell_centres(pa%prob_lo(1), pa%prob_hi(1), pa%nx, rc)
               do j = 1, pa%ny
                  w((j-1)*pa%nx+1:j*pa%nx) = sqrt(max(rc, 0.0_dp))
               end do
            end if
         end if

         idx = 0
         do j = 1, pa%ncomp
            X(idx+1:idx+pa%nx*pa%ny, k) = &
               reshape(pa%data(:,:,j) - pb%data(:,:,j), [pa%nx*pa%ny])
            if (weighted) X(idx+1:idx+pa%nx*pa%ny, k) = &
               X(idx+1:idx+pa%nx*pa%ny, k) * w
            idx = idx + pa%nx * pa%ny
         end do
         times(k) = pa%time

         if (verb .and. (mod(k-1, 25) == 0 .or. k == n_use)) then
            nrm = sqrt(sum(X(:,k)**2))
            write(*,'(a,i0,a,i0,a,es12.5,a,es12.5)') '  [', k, '/', n_use, &
               '] t=', times(k), '  |Q''|=', nrm
         end if

         call free_plotfile(pa)
         call free_plotfile(pb)
      end do

      call uniform_dt(times, dt, ierr)
      if (ierr /= 0) return

      if (verb) then
         write(*,'(a,i0,a,i0,a,es12.5)') 'subspace: ', size(X,1), &
            ' dof x ', size(X,2), ' snapshots, dt=', dt
      end if
   end subroutine build_snapshot_matrix

   !> Snapshot matrix from the solver's built-in LSA driver (restarting Krylov cycle).
   !>
   !> Member 0 ends in S_alt = N^M(S0); member k >= 1 ends in N^M(S0 + dS_{k-1}). The k-th column is the difference of the two final plotfiles, dS_k = N^M(S0 + dS_{k-1}) - S_alt, the "discretize then differentiate" Jacobian-vector product of the draft article, so the columns form the Krylov sequence dS_k ~ A^k dS_0 directly, with no growing amplitude to limit its length. The solver masks the level-set part of each member to a narrow band before it is added back (fort_LS_evec_max), so consecutive columns are related by A P with a fixed projection P; the nonzero eigenvalues are those of P A. The sampling interval is M dt, the time each member runs; the returned times are k M dt. The final plotfile of each member is the one with the largest step; all members must end at the same step and time. ierr: 1 no LSA plotfiles, 2 fewer than 2 members after member 0, 3 inconsistent final steps or times, 1000+/2000+ read errors, 5 incompatible plotfiles.
   subroutine build_krylov_matrix(dir, fields, first, last, X, times, dt, ierr, &
                                  weight, grid, domain)
      character(len=*), intent(in) :: dir
      character(len=*), intent(in) :: fields(:)
      integer, intent(in) :: first, last
      real(dp), allocatable, intent(out) :: X(:,:), times(:)
      real(dp), intent(out) :: dt
      integer, intent(out) :: ierr
      character(len=*), intent(in), optional :: weight
      integer, intent(out), optional :: grid(3)
      real(dp), intent(out), optional :: domain(4)

      character(len=MAX_NAME), allocatable :: paths(:)
      integer, allocatable :: steps(:), members(:), lastpath(:)
      integer :: n, i, k, kmax, mstep, nm, col, j, idx
      type(plotfile_t) :: p0, pk
      character(len=256) :: msg
      real(dp), allocatable :: rc(:), w(:)
      logical :: weighted

      ierr = 0; dt = 0.0_dp
      weighted = .false.
      if (present(weight)) weighted = (trim(weight) == 'volume')
      call find_lsa_members(dir, paths, steps, members, n)
      if (n == 0) then
         ierr = 1; return
      end if
      kmax = maxval(members)
      mstep = maxval(steps)
      allocate(lastpath(0:kmax))
      lastpath = 0
      do i = 1, n
         if (steps(i) == mstep) lastpath(members(i)) = i
      end do
      if (lastpath(0) == 0) then
         ierr = 3; return
      end if
      ! Members first..last (default 1..kmax) that reached the final step.
      nm = 0
      do k = max(first, 1), merge(kmax, min(last, kmax), last <= 0)
         if (lastpath(k) > 0) nm = nm + 1
      end do
      if (nm < 2) then
         ierr = 2; return
      end if
      call read_plotfile(trim(paths(lastpath(0))), p0, ierr, fields)
      if (ierr /= 0) then
         ierr = 2000 + ierr; return
      end if
      if (present(grid)) grid = [p0%nx, p0%ny, p0%ncomp]
      if (present(domain)) domain = [p0%prob_lo, p0%prob_hi]
      allocate(X(p0%nx*p0%ny*p0%ncomp, nm), times(nm))
      if (weighted) then
         allocate(rc(p0%nx), w(p0%nx*p0%ny))
         call cell_centres(p0%prob_lo(1), p0%prob_hi(1), p0%nx, rc)
         do j = 1, p0%ny
            w((j-1)*p0%nx+1:j*p0%nx) = sqrt(max(rc, 0.0_dp))
         end do
      end if
      dt = p0%time
      if (.not. (dt > 0.0_dp)) then
         ierr = 3; return
      end if
      col = 0
      do k = max(first, 1), merge(kmax, min(last, kmax), last <= 0)
         if (lastpath(k) == 0) cycle
         call read_plotfile(trim(paths(lastpath(k))), pk, ierr, fields)
         if (ierr /= 0) then
            ierr = 1000 + ierr; return
         end if
         if (.not. compatible(pk, p0, 1.0e-9_dp, msg)) then
            write(*,'(a)') 'incompatible Krylov member: '//trim(msg)
            ierr = 5; return
         end if
         col = col + 1
         idx = 0
         do j = 1, p0%ncomp
            X(idx+1:idx+p0%nx*p0%ny, col) = &
               reshape(pk%data(:,:,j) - p0%data(:,:,j), [p0%nx*p0%ny])
            if (weighted) X(idx+1:idx+p0%nx*p0%ny, col) = &
               X(idx+1:idx+p0%nx*p0%ny, col) * w
            idx = idx + p0%nx*p0%ny
         end do
         times(col) = real(k, dp) * dt
         call free_plotfile(pk)
      end do
      write(*,'(a,i0,a,i0,a,i0,a,es12.5)') 'Krylov sequence: members ', &
         max(first,1), '..', merge(kmax, min(last, kmax), last <= 0), &
         ' (', nm, ' columns), sampling interval M dt = ', dt
      call free_plotfile(p0)
   end subroutine build_krylov_matrix

   !> Sampling interval, insisting the snapshots are evenly spaced.
   !>
   !> DMD assumes a constant dt; non-uniform sampling biases every eigenvalue, so this is checked rather than silently averaged.
   subroutine uniform_dt(times, dt, ierr)
      real(dp), intent(in) :: times(:)
      real(dp), intent(out) :: dt
      integer, intent(out) :: ierr
      integer :: i, n
      real(dp) :: d, spread

      ierr = 0
      n = size(times)
      if (n < 2) then
         ierr = 10; dt = 0.0_dp; return
      end if
      do i = 2, n
         if (times(i) <= times(i-1)) then
            ierr = 11; return                ! not strictly increasing
         end if
      end do
      dt = (times(n) - times(1)) / real(n - 1, dp)
      spread = 0.0_dp
      do i = 2, n
         d = times(i) - times(i-1)
         spread = max(spread, abs(d - dt))
      end do
      if (spread > 1.0e-6_dp * max(dt, tiny(1.0_dp))) then
         ierr = 12                           ! non-uniform sampling
      end if
   end subroutine uniform_dt

   !> Column-wise L2 norms of the subspace.
   !>
   !> On a linear, exponentially growing mode this is a straight line on a log plot with slope equal to the growth rate, giving an independent check on the DMD result.
   subroutine perturbation_norms(X, norms)
      real(dp), intent(in) :: X(:,:)
      real(dp), allocatable, intent(out) :: norms(:)
      integer :: j
      allocate(norms(size(X,2)))
      do j = 1, size(X,2)
         norms(j) = sqrt(sum(X(:,j)**2))
      end do
   end subroutine perturbation_norms

   !> Least-squares slope of log|Q'| versus t over an optional window.
   !>
   !> A poor r^2 usually means the window still contains the initial transient, or that nonlinearity has set in at the late end. t_min and t_max may each be util_mod's UNSET (-huge) to leave that side open; a set t_max at or below a set t_min is reported as ierr = 21 rather than silently ignored. The fit itself is util_mod::fit_log_linear, shared with the interface-amplitude fit so that the two cannot drift apart. At least two samples are required.
   subroutine fit_exponential_growth(times, norms, t_min, t_max, &
                                     slope, intercept, r2, ierr)
      real(dp), intent(in) :: times(:), norms(:)
      real(dp), intent(in) :: t_min, t_max
      real(dp), intent(out) :: slope, intercept, r2
      integer, intent(out) :: ierr

      call fit_log_linear(times, norms, t_min, t_max, 2, slope, intercept, &
                          r2, ierr)
   end subroutine fit_exponential_growth

   !> Drift of one run from its first plotfile: for every plotfile k and every requested component c, dmax(c,k) = max |Q_c(t_k) - Q_c(t_1)| and drms(c,k) = sqrt( sum r (Q_c(t_k) - Q_c(t_1))^2 / sum r ), the root-mean-square change weighted by the cell-centre radius (the axisymmetric volume element). q0max(c) = max |Q_c(t_1)| gives the scale. For an unperturbed twin of a column at rest this is the discretisation drift the twin subtraction removes. Plotfiles are taken in step order. ierr: 1 fewer than two plotfiles, 100+e reading the first, 200+e reading a later one, 3 grid changed.
   subroutine base_drift(dir, fields, times, dmax, drms, q0max, grid, ierr, &
                         prefix)
      character(len=*), intent(in) :: dir
      character(len=*), intent(in) :: fields(:)
      real(dp), allocatable, intent(out) :: times(:), dmax(:,:), drms(:,:)
      real(dp), allocatable, intent(out) :: q0max(:)
      integer, intent(out) :: grid(3), ierr
      character(len=*), intent(in), optional :: prefix
      character(len=MAX_NAME), allocatable :: paths(:)
      character(len=64) :: pfx
      type(plotfile_t) :: p0, pk
      real(dp), allocatable :: rc(:), wr(:,:)
      integer :: n, e, k, c, nx, ny, nc

      ierr = 0
      grid = 0
      pfx = 'nddataPLT'
      if (present(prefix)) pfx = prefix
      call find_plotfiles(dir, paths, n, trim(pfx))
      if (n < 2) then
         ierr = 1; return
      end if
      call read_plotfile(trim(paths(1)), p0, e, fields)
      if (e /= 0) then
         ierr = 100 + e; return
      end if
      nx = p0%nx; ny = p0%ny; nc = p0%ncomp
      grid = [nx, ny, nc]
      allocate(rc(nx), times(n), dmax(nc, n), drms(nc, n), q0max(nc), wr(nx, ny))
      call cell_centres(p0%prob_lo(1), p0%prob_hi(1), nx, rc)
      wr = spread(rc, 2, ny)
      times(1) = p0%time
      dmax(:,1) = 0.0_dp
      drms(:,1) = 0.0_dp
      do c = 1, nc
         q0max(c) = maxval(abs(p0%data(:,:,c)))
      end do
      do k = 2, n
         call read_plotfile(trim(paths(k)), pk, e, fields)
         if (e /= 0) then
            ierr = 200 + e; return
         end if
         if (pk%nx /= nx .or. pk%ny /= ny) then
            ierr = 3; return
         end if
         times(k) = pk%time
         do c = 1, nc
            dmax(c,k) = maxval(abs(pk%data(:,:,c) - p0%data(:,:,c)))
            drms(c,k) = sqrt(sum(wr * (pk%data(:,:,c) - p0%data(:,:,c))**2) &
                             / sum(wr))
         end do
         call free_plotfile(pk)
      end do
      call free_plotfile(p0)
   end subroutine base_drift

end module nsmfp_mod
