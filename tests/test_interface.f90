!> Tests for interface_mode_mod and for the robustness checks of plotfile_mod.
!>
!> The interface-mode amplitude is the observable the growth rate is fitted to, yet before 2026-10-07 no test exercised it. Every case here plants an interface r = R(z) with a closed form in a genuine on-disk AMReX plotfile, as a level set that is exactly linear in r, so the recovered radius is exact to round-off and the recovered amplitude and growth rate can be checked to tight tolerances.
!>
!> The plotfiles are written under a scratch directory whose name contains a space, multi-box and big-endian files are mixed into one run, and decoy names are placed beside the real plotfiles, so that the reader is tested under the conditions it meets on a workstation folder and on the cluster rather than under the most convenient ones.
program test_interface
   use testing_mod
   use plotfile_mod
   use plotfile_writer_mod
   use interface_mode_mod
   use util_mod, only: UNSET, cell_centres
   use dispersion_mod, only: PI
   implicit none

   integer, parameter :: NR = 64, NZ = 72
   real(dp), parameter :: KR0 = 0.7_dp, RMAX = 4.0_dp, A0 = 0.02_dp
   real(dp) :: lambda, lo(2), hi(2), sigma_planted
   real(dp), allocatable :: z(:), radius(:), times(:), amps(:), meanr(:)
   real(dp), allocatable :: cen(:), sig(:)
   integer, allocatable :: nmiss(:)
   character(len=512) :: sroot, rdir
   character(len=64) :: names(1)
   character(len=MAX_NAME), allocatable :: paths(:)
   type(plotfile_t) :: pf
   real(dp) :: t, a, r2, rel, slast
   integer :: ierr, nzr, n, i, nm, nout
   logical :: conv

   lambda = 2.0_dp * PI / KR0
   lo = [0.0_dp, 0.0_dp]
   hi = [RMAX, lambda]
   names(1) = 'L0101'
   sroot = scratch_root()

   ! ------------------------------------------------------------------
   call section('interface radius and mode amplitude, full profile')
   call write_interface(trim(sroot)//'/lsa_if/full/nddataPLT00000', &
                        cosine(A0, 0.0_dp), 0.0_dp, 0, lo, hi)
   call interface_radius(trim(sroot)//'/lsa_if/full/nddataPLT00000', 'L0101', &
                         t, z, radius, nzr, ierr)
   call check_int(ierr, 0, 'plotfile in a path containing a space is read')
   call check_int(nzr, NZ, 'one radius per axial column')
   call check(maxval(abs(radius - (1.0_dp + A0*cos(KR0*z)))) < 1.0e-12_dp, &
              'radius matches the planted interface to round-off')
   a = mode_amplitude(z, radius, lambda, 1, nm)
   call check_close_rel(a, A0, 1.0e-10_dp, 'mode amplitude equals the planted one')
   call check_int(nm, 0, 'no missing columns reported')

   call section('mode amplitude is phase independent')
   call write_interface(trim(sroot)//'/lsa_if/phase/nddataPLT00000', &
                        cosine(A0, 1.1_dp), 0.0_dp, 0, lo, hi)
   call interface_radius(trim(sroot)//'/lsa_if/phase/nddataPLT00000', &
                         'L0101', t, z, radius, nzr, ierr)
   call check_close_rel(mode_amplitude(z, radius, lambda, 1), A0, 1.0e-10_dp, &
                        'a phase-shifted mode gives the same amplitude')

   call section('missing columns no longer bias the amplitude')
   ! Column 20 has no interface (the level set does not change sign), as when the interface touches the axis. The former projection treated the 71 remaining columns as a complete period and was off by about 2 per cent; the least-squares fit is exact.
   call interface_radius(trim(sroot)//'/lsa_if/full/nddataPLT00000', 'L0101', &
                         t, z, radius, nzr, ierr)
   radius(20) = -1.0_dp
   radius(41) = -1.0_dp
   a = mode_amplitude(z, radius, lambda, 1, nm)
   call check_close_rel(a, A0, 1.0e-10_dp, &
                        'amplitude exact with two columns missing')
   call check_int(nm, 2, 'the two missing columns are reported')
   call check(abs(old_projection(z, radius, lambda) - A0) > 1.0e-3_dp * A0, &
              'whereas the former projection was biased')
   radius(1:NZ-3) = -1.0_dp
   call check_close(mode_amplitude(z, radius, lambda, 1), -1.0_dp, 0.0_dp, &
                    'fewer than 4 valid columns gives -1')

   call section('a level set exactly zero at a cell centre')
   ! R(z) = r_30, a cell centre: the level set is 0.0 there, and the former test a*b < 0 missed the crossing in every column.
   block
      real(dp) :: rc(NR)
      call cell_centres(0.0_dp, RMAX, NR, rc)
      call write_interface(trim(sroot)//'/lsa_if/zero/nddataPLT00000', &
                           flat(rc(30)), 0.0_dp, 0, lo, hi)
      call interface_radius(trim(sroot)//'/lsa_if/zero/nddataPLT00000', &
                            'L0101', t, z, radius, nzr, ierr)
      call check(all(radius == rc(30)), &
                 'an exact zero is taken as the interface')
   end block

   call section('domain not starting at zero')
   ! prob_lo = (0, 5): the axial coordinates must start at 5 + dz/2, not at dz/2. The planted interface uses the absolute z, so an offset axis would shift the phase and the radius comparison would fail.
   call write_interface(trim(sroot)//'/lsa_if/off/nddataPLT00000', &
                        cosine(A0, 0.0_dp), 0.0_dp, 0, [0.0_dp, 5.0_dp], &
                        [RMAX, 5.0_dp + lambda])
   call interface_radius(trim(sroot)//'/lsa_if/off/nddataPLT00000', 'L0101', &
                         t, z, radius, nzr, ierr)
   call check_close(z(1), 5.0_dp + 0.5_dp*lambda/real(NZ,dp), 1.0e-12_dp, &
                    'axial coordinates start from prob_lo')
   call check(maxval(abs(radius - (1.0_dp + A0*cos(KR0*z)))) < 1.0e-12_dp, &
              'radius correct on the offset domain')

   ! ------------------------------------------------------------------
   call section('growth record: amplitude series and fits')
   ! A run of 41 plotfiles with a(t) = A0 exp(sigma t). Every third file is split into three boxes and every fourth is big-endian, so the series also exercises the multi-box scatter and the byte swap.
   sigma_planted = 0.3289_dp
   rdir = trim(sroot)//'/lsa_if/run'
   do i = 0, 40
      t = 0.1_dp * real(i, dp)
      call write_interface(trim(rdir)//'/nddataPLT'//step_name(10*i), &
                           cosine(A0*exp(sigma_planted*t), 0.0_dp), t, 10*i, &
                           lo, hi, nsplit=merge(3, 1, mod(i,3) == 0), &
                           big_endian=(mod(i,4) == 0))
   end do
   ! Decoys beside the real plotfiles: a renamed plotfile directory and a plain file whose name matches the old glob. Both must be ignored.
   call write_interface(trim(rdir)//'/nddataPLT00000099.old', &
                        cosine(1.0_dp, 0.0_dp), 9.9_dp, 99, lo, hi)
   call execute_command_line("touch '"//trim(rdir)//"/nddataPLT00009999'")

   call find_plotfiles(trim(rdir), paths, n)
   call check_int(n, 41, 'decoy directory and plain file are ignored')

   call amplitude_series(trim(rdir), lambda, 'L0101', 1, times, amps, meanr, &
                         n, ierr, nmiss=nmiss)
   call check_int(ierr, 0, 'amplitude series succeeds')
   call check_int(n, 41, 'one amplitude per plotfile')
   call check(all(times(2:n) > times(1:n-1)), 'series is sorted by time')
   call check(maxval(abs(amps/(A0*exp(sigma_planted*times)) - 1.0_dp)) &
              < 1.0e-9_dp, 'every amplitude matches the planted growth')
   call check(maxval(abs(meanr - 1.0_dp)) < 1.0e-12_dp, &
              'mean radius stays at r0')
   call check_int(sum(nmiss), 0, 'no missing columns in the record')

   call fit_growth_rate_amp(times, amps, UNSET, UNSET, a, r2, ierr)
   call check_int(ierr, 0, 'growth-rate fit succeeds')
   call check_close_rel(a, sigma_planted, 1.0e-9_dp, &
                        'fitted growth rate is the planted one')
   call check(r2 > 0.9999999_dp, 'fit is essentially perfect')
   call fit_growth_rate_amp(times, amps, 2.0_dp, 1.0_dp, a, r2, ierr)
   call check(ierr /= 0, 'an inverted window is an error, not ignored')

   call local_growth_rate(times, amps, 0.8_dp, cen, sig, nout)
   call check(nout >= 2, 'local growth rate has several windows')
   call check(maxval(abs(sig(1:nout) - sigma_planted)) < 1.0e-9_dp, &
              'local growth rate is flat at the planted value')
   conv = has_converged(times, amps, 0.8_dp, 1.0e-3_dp, slast, rel)
   call check(conv, 'a pure exponential is reported as converged')
   call local_growth_rate(times(1:4), amps(1:4), 0.25_dp, cen, sig, nout)
   call check_int(nout, 0, 'windows of fewer than 4 samples are not fitted')

   ! ------------------------------------------------------------------
   call section('plotfile consistency checks')
   call write_interface(trim(sroot)//'/lsa_if/hole/nddataPLT00000', &
                        cosine(A0, 0.0_dp), 0.0_dp, 0, lo, hi, nsplit=3, &
                        drop_box=2)
   call read_plotfile(trim(sroot)//'/lsa_if/hole/nddataPLT00000', pf, ierr)
   call check_int(ierr, 13, 'boxes not covering the domain are rejected')
   call free_plotfile(pf)

   call write_interface(trim(sroot)//'/lsa_if/single/nddataPLT00000', &
                        cosine(A0, 0.0_dp), 0.0_dp, 0, lo, hi)
   call execute_command_line("sed -i 's/),(8, (8 7 6 5 4 3 2 1)/),(4, (4 3 2 1)/' '"// &
      trim(sroot)//"/lsa_if/single/nddataPLT00000/Level_0/Cell_D_00000'")
   call read_plotfile(trim(sroot)//'/lsa_if/single/nddataPLT00000', pf, ierr)
   call check_int(ierr, 612, 'single-precision data is refused, not misread')
   call free_plotfile(pf)

   call check_int(plotfile_step('/a b/nddataPLT00000120'), 120, &
                  'step parsed from a path containing a space')
   call check_int(plotfile_step('/x/nddataPLT'), -1, 'no digits gives -1')

   call execute_command_line("rm -rf '"//trim(sroot)//"'")
   call summary('test_interface')

contains

   !> Interface shapes, encoded as (kind, amplitude, phase/level).
   function cosine(amp, phase) result(spec)
      real(dp), intent(in) :: amp, phase
      real(dp) :: spec(3)
      spec = [1.0_dp, amp, phase]
   end function cosine

   function flat(level) result(spec)
      real(dp), intent(in) :: level
      real(dp) :: spec(3)
      spec = [2.0_dp, 0.0_dp, level]
   end function flat

   !> Write a plotfile whose level set L0101 = r - R(z) is negative in the liquid; R(z) = 1 + amp cos(k z + phase) for a cosine spec, a constant for a flat one.
   subroutine write_interface(path, spec, time, step, plo, phi, nsplit, &
                              drop_box, big_endian)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: spec(3), time, plo(2), phi(2)
      integer, intent(in) :: step
      integer, intent(in), optional :: nsplit, drop_box
      logical, intent(in), optional :: big_endian
      real(dp) :: d(NR, NZ, 1), rc(NR), zc(NZ), rz
      integer :: i, j
      call cell_centres(plo(1), phi(1), NR, rc)
      call cell_centres(plo(2), phi(2), NZ, zc)
      do j = 1, NZ
         if (spec(1) < 1.5_dp) then
            rz = 1.0_dp + spec(2) * cos(KR0 * zc(j) + spec(3))
         else
            rz = spec(3)
         end if
         do i = 1, NR
            d(i, j, 1) = rc(i) - rz
         end do
      end do
      call write_plotfile(path, d, names, time, step, plo, phi, &
                          coord_sys=1, big_endian=big_endian, &
                          nsplit=nsplit, drop_box=drop_box)
   end subroutine write_interface

   function step_name(k) result(s)
      integer, intent(in) :: k
      character(len=8) :: s
      write(s, '(i8.8)') k
   end function step_name

   !> The former estimator: a = 2 |<(r - <r>) exp(-i k z)>| over the valid columns only.
   function old_projection(zz, rr, lam) result(amp)
      real(dp), intent(in) :: zz(:), rr(:), lam
      real(dp) :: amp, rbar, sr, si, ang
      integer :: j, m
      m = count(rr >= 0.0_dp)
      rbar = sum(rr, mask=rr >= 0.0_dp) / real(m, dp)
      sr = 0.0_dp; si = 0.0_dp
      do j = 1, size(rr)
         if (rr(j) < 0.0_dp) cycle
         ang = 2.0_dp * PI * zz(j) / lam
         sr = sr + (rr(j) - rbar) * cos(ang)
         si = si - (rr(j) - rbar) * sin(ang)
      end do
      amp = 2.0_dp * sqrt(sr*sr + si*si) / real(m, dp)
   end function old_projection

end program test_interface
