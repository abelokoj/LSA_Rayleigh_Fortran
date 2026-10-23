!> Tests for the features added for the 2026-10 run campaign from the recommendations of the September 2026 evidence report: the admissible-window analysis of amplitude records (record_mod), the DMD residual and archives (dmd_mod, archive_mod), the volume-weighted inner product (nsmfp_mod), the signed Fourier coefficients and mean(r^2) of the interface record (interface_mode_mod), and the long-box and eps = c dr options of the deck generator (make_inputs_mod).
program test_campaign
   use testing_mod
   use util_mod, only: UNSET, cell_centres
   use record_mod
   use dmd_mod
   use archive_mod
   use nsmfp_mod
   use plotfile_writer_mod
   use plotfile_mod, only: MAX_NAME, find_lsa_members
   use interface_mode_mod
   use make_inputs_mod
   use dispersion_mod, only: PI
   implicit none

   character(len=512) :: sroot

   sroot = scratch_root()
   call test_record()
   call test_residual_and_archive()
   call test_weight()
   call test_signed_coefficients()
   call test_decks()
   call test_krylov()
   call test_drift()
   call execute_command_line("rm -rf '"//trim(sroot)//"'")
   call summary('test_campaign')

contains

   !> a = a0 cosh(s t): the two-mode startup of an interface displaced from rest. The local rate s tanh(s t) climbs to s; the window must start late enough that the fit is within the settling tolerance of s, and must end at the amplitude bound.
   subroutine test_record()
      integer, parameter :: N = 200
      real(dp) :: t(N), a(N)
      real(dp), allocatable :: tr(:), ar(:), c(:), r(:)
      type(record_window_t) :: win
      real(dp), parameter :: S = 0.33_dp, A0 = 0.005_dp
      integer :: i, u, e, nout
      character(len=:), allocatable :: path

      call section('admissible window of an amplitude record')
      do i = 1, N
         t(i) = 0.06_dp * real(i - 1, dp)
         a(i) = A0 * cosh(S * t(i))
      end do
      call analyse_record(t, a, 1.0_dp, 0.1_dp, 0.8_dp, 0.005_dp, UNSET, &
                          UNSET, win, c, r, nout)
      call check_int(win%ierr, 0, 'window fit succeeds')
      call check(win%bounded, 'amplitude bound a/r0 = 0.1 is detected')
      call check(win%t_hi <= acosh(0.1_dp/A0)/S + 1.0e-12_dp, &
                 'window ends before a exceeds the bound')
      call check(win%settled, 'a settled start is found')
      call check(abs(win%rate - S)/S < 0.005_dp, &
                 'rate within the settling tolerance of the true rate')
      call check(win%rate < S, 'and below it, as the cosh transient implies')
      call analyse_record(t, a, 1.0_dp, 0.1_dp, 0.8_dp, 0.005_dp, 1.0_dp, &
                          UNSET, win, c, r, nout)
      call check(win%rate < 0.98_dp*S, &
                 'a window starting at t = 1 is biased low by the transient')

      ! A record file with comments, a -1 missing-value row and 6 columns.
      path = trim(sroot)//'/rec.dat'
      open(newunit=u, file=path, status='replace', action='write')
      write(u,'(a)') '# time  amplitude  mean_radius  a_cos  a_sin  mean_r2'
      do i = 1, N
         if (i == 50) then
            write(u,'(6es24.16)') t(i), -1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp
         else
            write(u,'(6es24.16)') t(i), a(i), 1.0_dp, 0.0_dp, a(i), 1.0_dp
         end if
      end do
      close(u)
      call read_amplitude_record(path, tr, ar, e)
      call check_int(e, 0, 'record file is read')
      call check_int(size(tr), N - 1, 'missing-value row is dropped')
      call check(abs(ar(10) - a(10)) < 1.0e-15_dp, 'values read exactly')
   end subroutine test_record

   subroutine test_residual_and_archive()
      integer, parameter :: NS = 120, M = 40
      real(dp) :: X(NS, M), Xn(NS, M), tt
      real(dp), allocatable :: Xr(:,:), times(:)
      type(dmd_result_t) :: res
      real(dp) :: dt, dt2, times0(M)
      integer :: i, j, e, grid(3), nlines, u, ios
      character(len=256) :: line
      character(len=64) :: f(2)

      call section('DMD residual and archives')
      dt = 0.1_dp
      do j = 1, M
         tt = real(j-1, dp) * dt
         times0(j) = tt
         do i = 1, NS
            X(i, j) = sin(0.1_dp*i)*exp(0.3_dp*tt) + cos(0.05_dp*i)*exp(-0.5_dp*tt)
            ! a quadratic term makes the sequence nonlinear
            Xn(i, j) = X(i, j) + 0.2_dp*(sin(0.1_dp*i)*exp(0.3_dp*tt))**2
         end do
      end do
      call dmd_compute(X, dt, 0, 1.0e-10_dp, res, e)
      call check(res%residual < 1.0e-10_dp, 'a linear sequence has a tiny residual')
      call dmd_free(res)
      call dmd_compute(Xn, dt, 2, 1.0e-10_dp, res, e)
      call check(res%residual > 1.0e-4_dp, &
                 'a nonlinear sequence at rank 2 shows a clear residual')
      call dmd_free(res)

      call dmd_compute(X, dt, 0, 1.0e-10_dp, res, e)
      f(1) = 'x_velocity'; f(2) = 'y_velocity'
      call write_dmd_archive(trim(sroot)//'/arch a', res, times0, 0.3_dp, f, &
                             'none', [10, 6, 2], 0, 1, 5, e)
      call check_int(e, 0, 'DMD archive written (path with a space)')
      nlines = 0
      open(newunit=u, file=trim(sroot)//'/arch a/eigenvalues.txt', &
           status='old', action='read', iostat=ios)
      call check_int(ios, 0, 'eigenvalues.txt exists')
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:1) /= '#') nlines = nlines + 1
      end do
      close(u)
      call check_int(nlines, res%n_modes, 'one eigenvalue line per mode')
      call dmd_free(res)

      call write_perturbation_export(trim(sroot)//'/exp', X, times0, dt, f, &
           [10, 6, 2], [0.0_dp, 0.0_dp, 4.0_dp, 9.0_dp], 'p', 'b', 'none', e)
      call check_int(e, 0, 'perturbation export written')
      call read_perturbation_export(trim(sroot)//'/exp', Xr, times, dt2, grid, e)
      call check_int(e, 0, 'perturbation export read back')
      call check(all(Xr == X) .and. all(times == times0) .and. dt2 == dt, &
                 'export round-trips exactly (lossless)')
      call check(all(grid == [10, 6, 2]), 'grid shape preserved')
   end subroutine test_residual_and_archive

   !> The volume weight multiplies each row by sqrt(r) at its cell centre.
   subroutine test_weight()
      integer, parameter :: NX = 8, NY = 6
      real(dp) :: dp_(NX, NY, 1), db(NX, NY, 1), rc(NX), dtt
      real(dp), allocatable :: X0(:,:), X1(:,:), tt(:)
      character(len=64) :: nm(1)
      integer :: e, i, j, k, grid(3)
      real(dp) :: dom(4), err

      call section('volume-weighted inner product')
      nm(1) = 'x_velocity'
      do j = 1, NY
         do i = 1, NX
            db(i, j, 1) = 0.1_dp * i + j
         end do
      end do
      do k = 0, 2
         dp_ = db + real(k + 1, dp) * 0.01_dp
         call write_plotfile(trim(sroot)//'/w/p/nddataPLT0000'//char(48+k), &
              dp_, nm, 0.1_dp*k, k, [0.0_dp, 0.0_dp], [4.0_dp, 3.0_dp])
         call write_plotfile(trim(sroot)//'/w/b/nddataPLT0000'//char(48+k), &
              db, nm, 0.1_dp*k, k, [0.0_dp, 0.0_dp], [4.0_dp, 3.0_dp])
      end do
      call build_snapshot_matrix(trim(sroot)//'/w/p', trim(sroot)//'/w/b', nm, &
           0, 1, 0, X0, tt, dtt, e, grid=grid, domain=dom)
      call check_int(e, 0, 'unweighted subspace built')
      call check(all(grid == [NX, NY, 1]), 'grid shape returned')
      call check(all(dom == [0.0_dp, 0.0_dp, 4.0_dp, 3.0_dp]), 'domain returned')
      call build_snapshot_matrix(trim(sroot)//'/w/p', trim(sroot)//'/w/b', nm, &
           0, 1, 0, X1, tt, dtt, e, weight='volume')
      call check_int(e, 0, 'volume-weighted subspace built')
      call cell_centres(0.0_dp, 4.0_dp, NX, rc)
      err = 0.0_dp
      do k = 1, 3
         do j = 1, NY
            do i = 1, NX
               err = max(err, abs(X1(i+(j-1)*NX, k) - sqrt(rc(i))*X0(i+(j-1)*NX, k)))
            end do
         end do
      end do
      call check(err < 1.0e-14_dp, 'rows are scaled by sqrt(r)')
      call build_snapshot_matrix(trim(sroot)//'/w/p', trim(sroot)//'/w/b', nm, &
           0, 1, 0, X1, tt, dtt, e, weight='bogus')
      call check_int(e, 7, 'an unknown weight is rejected')
   end subroutine test_weight

   !> R(z) = 1 + a cos(kz) + b sin(kz): the fit returns (a, b) with their signs.
   subroutine test_signed_coefficients()
      integer, parameter :: NR = 48, NZ = 64
      real(dp) :: d(NR, NZ, 1), rc(NR), zc(NZ), lam, cf(2), amp
      real(dp), allocatable :: z(:), radius(:), times(:), amps(:), meanr(:)
      real(dp), allocatable :: coefs(:,:), m2(:)
      character(len=64) :: nm(1)
      integer :: i, j, e, nzr, n
      real(dp) :: t
      real(dp), parameter :: AC = -0.012_dp, BS = 0.007_dp

      call section('signed Fourier coefficients and mean(r^2)')
      nm(1) = 'L0101'
      lam = 2.0_dp*PI/0.7_dp
      call cell_centres(0.0_dp, 4.0_dp, NR, rc)
      call cell_centres(0.0_dp, lam, NZ, zc)
      do j = 1, NZ
         do i = 1, NR
            d(i, j, 1) = rc(i) - (1.0_dp + AC*cos(0.7_dp*zc(j)) + BS*sin(0.7_dp*zc(j)))
         end do
      end do
      call write_plotfile(trim(sroot)//'/sc/nddataPLT00000', d, nm, 0.0_dp, 0, &
                          [0.0_dp, 0.0_dp], [4.0_dp, lam], coord_sys=1)
      call interface_radius(trim(sroot)//'/sc/nddataPLT00000', 'L0101', t, z, &
                            radius, nzr, e)
      amp = mode_amplitude(z, radius, lam, 1, coef=cf)
      call check(abs(cf(1) - AC) < 1.0e-12_dp .and. abs(cf(2) - BS) < 1.0e-12_dp, &
                 'cos and sin coefficients recovered with their signs')
      call check(abs(amp - sqrt(AC**2 + BS**2)) < 1.0e-12_dp, &
                 'amplitude is their magnitude')
      call amplitude_series(trim(sroot)//'/sc', lam, 'L0101', 1, times, amps, &
                            meanr, n, e, coefs=coefs, meanr2=m2)
      call check_int(e, 0, 'series with coefficients succeeds')
      call check(abs(coefs(1,1) - AC) < 1.0e-12_dp, 'series carries the coefficients')
      call check(abs(m2(1) - (1.0_dp + 0.5_dp*(AC**2 + BS**2))) < 1.0e-9_dp, &
                 'mean(r^2) = r0^2 + a^2/2 for a sinusoid')
   end subroutine test_signed_coefficients

   !> base_drift: a run whose level set shifts by D(t) = 1e-3 t and whose velocity grows as 2e-4 t everywhere; max and rms changes from the first plotfile must equal those values.
   subroutine test_drift()
      integer, parameter :: NR = 12, NZ = 10, NT = 4
      real(dp) :: d(NR, NZ, 2), rc(NR), t
      real(dp), allocatable :: tt(:), dmax(:,:), drms(:,:), q0(:)
      character(len=64) :: nm(2)
      character(len=MAX_NAME) :: path
      integer :: i, j, k, e, grid(3)

      call section('base-state drift')
      nm = [character(len=64) :: 'L0101', 'x_velocity']
      call cell_centres(0.0_dp, 4.0_dp, NR, rc)
      do k = 0, NT - 1
         t = 0.5_dp * k
         do j = 1, NZ
            do i = 1, NR
               d(i, j, 1) = 1.0_dp - rc(i) + 1.0e-3_dp * t
               d(i, j, 2) = 2.0e-4_dp * t
            end do
         end do
         write(path, '(a,i5.5)') trim(sroot)//'/drift/nddataPLT', 10 * k
         call write_plotfile(trim(path), d, nm, t, 10 * k, [0.0_dp, 0.0_dp], &
                             [4.0_dp, 3.0_dp], coord_sys=1)
      end do
      call base_drift(trim(sroot)//'/drift', nm, tt, dmax, drms, q0, grid, e)
      call check_int(e, 0, 'base_drift succeeds')
      call check_int(size(tt), NT, 'all plotfiles read')
      call check(abs(dmax(1, NT) - 1.5e-3_dp) < 1.0e-12_dp, 'level-set max change = 1e-3 t')
      call check(abs(drms(2, NT) - 3.0e-4_dp) < 1.0e-12_dp, 'velocity rms change = 2e-4 t')
      call check(dmax(1, 1) == 0.0_dp .and. q0(2) == 0.0_dp, 'first plotfile is the reference')
      call check_int(grid(1), NR, 'grid returned')
   end subroutine test_drift

   !> Members of the solver's built-in LSA driver: nddataPLT<step>LSA<k>, each restarting at step 0. Member k ends in base + lambda^k shape, so the differences from member 0 form an exact Krylov sequence with growth rate log(lambda)/T.
   subroutine test_krylov()
      integer, parameter :: NX = 8, NY = 10, NK = 12
      real(dp), parameter :: T = 0.08_dp, SIG = 0.33_dp
      real(dp) :: d(NX, NY, 1), b(NX, NY, 1), shp(NX, NY, 1), dtk
      real(dp), allocatable :: XX(:,:), tt(:)
      type(dmd_result_t) :: res
      character(len=64) :: nm(1)
      character(len=256) :: kd, p
      character(len=MAX_NAME), allocatable :: paths(:)
      integer, allocatable :: steps(:), members(:)
      integer :: i, j, k, e, n

      call section('built-in LSA driver: Krylov members')
      nm(1) = 'L0101'
      kd = trim(sroot)//'/kry run'
      do j = 1, NY
         do i = 1, NX
            b(i, j, 1) = real(i, dp) * 0.1_dp - 1.0_dp
            shp(i, j, 1) = 1.0e-3_dp * sin(0.3_dp*i) * cos(0.6_dp*j)
         end do
      end do
      do k = 0, NK
         if (k == 0) then
            d = b
         else
            d = b + shp * exp(SIG * T * real(k, dp))
         end if
         write(p, '(a,a,i8.8,a,i8.8)') trim(kd), '/nddataPLT', 50, 'LSA', k
         call write_plotfile(trim(p), d, nm, T, 50, [0.0_dp, 0.0_dp], [4.0_dp, 9.0_dp])
         ! an intermediate plotfile of the same member, which must not be used
         write(p, '(a,a,i8.8,a,i8.8)') trim(kd), '/nddataPLT', 25, 'LSA', k
         call write_plotfile(trim(p), d * 0.5_dp, nm, 0.5_dp*T, 25, [0.0_dp, 0.0_dp], [4.0_dp, 9.0_dp])
      end do
      call find_lsa_members(trim(kd), paths, steps, members, n)
      call check_int(n, 2*(NK+1), 'all member plotfiles found')
      call check(members(1) == 0 .and. members(n) == NK, 'sorted by member')
      call build_krylov_matrix(trim(kd), nm, 1, 0, XX, tt, dtk, e)
      call check_int(e, 0, 'Krylov matrix built')
      call check_int(size(XX, 2), NK, 'one column per member after member 0')
      call check_close_rel(dtk, T, 1.0e-12_dp, 'sampling interval is the member length')
      call dmd_compute(XX, dtk, 0, 1.0e-10_dp, res, e)
      call check_close_rel(dominant_growth_rate(res), SIG, 1.0e-8_dp, &
                           'DMD of the members recovers the growth rate')
      call dmd_free(res)
      call build_krylov_matrix(trim(kd), nm, 3, 6, XX, tt, dtk, e)
      call check_int(size(XX, 2), 4, 'member range 3..6 selects 4 columns')
   end subroutine test_krylov

   subroutine test_decks()
      type(run_plan_t) :: p1, p8, pe
      character(len=:), allocatable :: t1, t8
      character(len=*), parameter :: TMPL = 'inputs/inputs_rayleigh_template.txt'
      integer :: e
      logical :: exists

      call section('long box and eps = c dr')
      p1%kr0 = 0.7_dp; p1%eps = 5.0e-3_dp; p1%cells_per_r0 = 8
      p8 = p1; p8%n_waves = 8; p8%rmax = 9.0_dp
      call plan_run(p1); call plan_run(p8)
      call check(abs(p8%zlen - 8.0_dp*p1%lambda) < 1.0e-12_dp, 'box holds 8 wavelengths')
      call check(abs(p8%lambda - p1%lambda) < 1.0e-14_dp, 'the imposed wavelength is unchanged')
      call check(abs(p8%dz - p8%dr) < 0.05_dp*p8%dr, 'cells stay nearly square')
      call check(p8%fixed_dt == p1%fixed_dt, 'same cell size gives the same time step')
      pe%kr0 = 0.7_dp; pe%eps = 1.0_dp; pe%cells_per_r0 = 16; pe%eps_cells = 0.08_dp
      call plan_run(pe)
      call check(abs(pe%eps - 0.08_dp*pe%dr) < 1.0e-15_dp, 'eps = c dr when eps_cells is set')
      pe%eps = 0.0_dp; call plan_run(pe)
      call check(pe%eps == 0.0_dp, 'the twin keeps eps = 0')

      inquire(file=TMPL, exist=exists)
      if (.not. exists) return
      call render_deck(p1, TMPL, t1, e)
      call render_deck(p8, TMPL, t8, e)
      call check_int(e, 0, 'long-box deck renders')
      call check(index(t8, 'geometry.prob_hi     = 9.000000000 71.80783') > 0, &
                 'long-box deck has prob_hi = 9 x 8 lambda')
      call check(index(t1, 'geometry.prob_hi     = 4.000000000 8.975979') > 0, &
                 'single-wavelength deck unchanged')
   end subroutine test_decks

end program test_campaign
