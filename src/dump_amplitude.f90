!> Write the interface-mode amplitude record of a calculation as plain text.
!>
!> Columns (since 2026-10-07): time, mode amplitude, mean interface radius, signed cos coefficient, signed sin coefficient, mean of r^2. Readers that use only the first three columns (all the plotting scripts) are unaffected; the signed coefficients serve the stable-band analysis (plotting/fit_stable_band.py) and mean(r^2) the volume check.
!>
!> The record is the quantity the dispersion relation predicts, namely the amplitude of the k-Fourier component of the interface radius; see interface_mode_mod for why a field norm is unsuitable for an interfacial instability. The mean interface radius is included because a drifting mean signals a mass-conservation error, which would invalidate the growth rate.
!>
!> Writing plain text rather than plotting directly keeps this implementation free of any Python dependency; the figure scripts in plotting/ read these files, and they may equally be plotted elsewhere.
!>
!> ./bin/dump_amplitude <rundir> <wavelength> <outfile>
program dump_amplitude
   use interface_mode_mod
   use plotfile_mod, only: find_plotfiles, MAX_NAME
   use util_mod, only: read_real, read_int, argsort_real_asc
   implicit none

   character(len=512) :: rundir, lambda_text, outfile
   real(dp), allocatable :: times(:), amps(:), mean_radius(:)
   real(dp) :: wavelength
   integer :: n, ierr, i, unit_out, ios, nbad
   integer, allocatable :: nmiss(:)
   real(dp), allocatable :: coefs(:,:), meanr2(:)
   character(len=64) :: nharm_text
   integer :: nharm

   if (command_argument_count() < 3) then
      write(*,'(a)') 'usage: dump_amplitude <rundir> <wavelength> <outfile> [nharm]'
      write(*,'(a)') '  nharm > 1: also write harmonics 1..nharm of <wavelength> '// &
         '(broadband run, <wavelength> = box length) to '// &
         '<outfile without .dat>_harmonics/hNN.dat'
      stop 1
   end if
   call get_command_argument(1, rundir)
   call get_command_argument(2, lambda_text)
   call get_command_argument(3, outfile)
   call read_real(lambda_text, wavelength, ierr)
   if (ierr /= 0 .or. .not. (wavelength > 0.0_dp)) then
      write(*,'(a)') 'ERROR: <wavelength> must be a positive number'
      stop 1
   end if

   nharm = 1
   if (command_argument_count() >= 4) then
      call get_command_argument(4, nharm_text)
      call read_int(nharm_text, nharm, ierr)
      if (ierr /= 0 .or. nharm < 1) then
         write(*,'(a)') 'ERROR: [nharm] must be a positive integer'
         stop 1
      end if
   end if

   ! The level set of the first material, L0101, carries the interface; the first harmonic is the mode imposed by the deck.
   call amplitude_series(trim(rundir), wavelength, 'L0101', 1, &
                         times, amps, mean_radius, n, ierr, nmiss=nmiss, &
                         coefs=coefs, meanr2=meanr2)
   if (ierr /= 0) then
      write(*,'(a,i0)') 'ERROR reading the run, ierr=', ierr
      if (ierr == 1) write(*,'(a)') &
         '  no plotfiles found; the solver writes nddataPLT* directories'
      stop 1
   end if

   open(newunit=unit_out, file=trim(outfile), status='replace', &
        action='write', iostat=ios)
   if (ios /= 0) then
      write(*,'(a,a)') 'ERROR: cannot open for writing: ', trim(outfile)
      stop 1
   end if
   ! A value of -1 in the amplitude or mean-radius column marks a snapshot in which too few columns held an interface to fit the mode; the plotting scripts treat negative values as missing.
   write(unit_out,'(a)') '# time  amplitude  mean_radius  a_cos  a_sin  mean_r2'
   do i = 1, n
      write(unit_out,'(6es24.16)') times(i), amps(i), mean_radius(i), &
         coefs(1,i), coefs(2,i), meanr2(i)
   end do
   close(unit_out)

   write(*,'(a,i0,a,a)') 'wrote ', n, ' rows to ', trim(outfile)

   if (nharm > 1) call write_harmonics()

   ! Snapshots whose interface profile had gaps. The amplitude is still a valid least-squares estimate, but a gap usually means the interface reached the axis or the outer wall, which is worth knowing before a growth rate is fitted to the record.
   nbad = count(nmiss > 0)
   if (nbad > 0) then
      write(*,'(a,i0,a,i0,a)') 'WARNING: ', nbad, ' of ', n, &
         ' snapshots had axial columns with no interface (max ' // &
         trim(itoa(maxval(nmiss))) // ' columns); inspect the record'
   end if

contains

   !> Broadband run: amplitude of each harmonic n = 1..nharm of the box (wavelength / n) at every snapshot, from one interface profile per plotfile, one record per harmonic in the usual six-column format.
   subroutine write_harmonics()
      character(len=MAX_NAME), allocatable :: paths(:)
      character(len=512) :: base, hdir, hfile
      real(dp), allocatable :: z(:), radius(:), tt(:), amp(:,:), cf(:,:,:), mr(:), m2(:)
      integer, allocatable :: order(:)
      integer :: np, k, h, e, nz, u, cnt
      real(dp) :: t, c2(2)
      base = outfile
      if (len_trim(base) > 4) then
         if (base(len_trim(base)-3:len_trim(base)) == '.dat') base = base(1:len_trim(base)-4)
      end if
      hdir = trim(base)//'_harmonics'
      call execute_command_line("mkdir -p '"//trim(hdir)//"'")
      call find_plotfiles(trim(rundir), paths, np, 'nddataPLT')
      allocate(tt(np), amp(nharm, np), cf(2, nharm, np), mr(np), m2(np), order(np))
      do k = 1, np
         call interface_radius(trim(paths(k)), 'L0101', t, z, radius, nz, e)
         if (e /= 0) then
            write(*,'(a,a)') 'ERROR reading ', trim(paths(k)); stop 1
         end if
         tt(k) = t
         cnt = count(radius >= 0.0_dp)
         mr(k) = -1.0_dp; m2(k) = -1.0_dp
         if (cnt > 0) then
            mr(k) = sum(radius, mask=radius >= 0.0_dp) / cnt
            m2(k) = sum(radius**2, mask=radius >= 0.0_dp) / cnt
         end if
         do h = 1, nharm
            amp(h, k) = mode_amplitude(z, radius, wavelength, h, coef=c2)
            cf(:, h, k) = c2
         end do
      end do
      call argsort_real_asc(tt, order)
      do h = 1, nharm
         write(hfile,'(a,a,i2.2,a)') trim(hdir), '/h', h, '.dat'
         open(newunit=u, file=trim(hfile), status='replace', action='write')
         write(u,'(a,i0,a,es24.16)') '# harmonic ', h, '  wavelength ', wavelength / real(h, dp)
         write(u,'(a)') '# time  amplitude  mean_radius  a_cos  a_sin  mean_r2'
         do k = 1, np
            write(u,'(6es24.16)') tt(order(k)), amp(h, order(k)), mr(order(k)), &
               cf(1, h, order(k)), cf(2, h, order(k)), m2(order(k))
         end do
         close(u)
      end do
      write(*,'(a,i0,a,a)') 'wrote ', nharm, ' harmonic records to ', trim(hdir)
   end subroutine write_harmonics

   function itoa(k) result(t)
      integer, intent(in) :: k
      character(len=16) :: t
      write(t,'(i0)') k
   end function itoa
end program dump_amplitude
