!> Archives of a DMD result and of the matched perturbation sequence.
!>
!> The September 2026 evidence report found that the study had implemented DMD but kept none of its output: no eigenvalue table, no singular values, no retained rank, no modes, no residual, and no uniformly sampled matched perturbation sequence from which any of these could be recomputed once the plotfiles were gone. Its first recommendation (section 11.1) was to save exactly those. This module writes them, in plain formats that need no Fortran to read.
!>
!> DMD archive, a directory:
!>
!> - summary.txt: key = value lines (growth rate, rank, residual, dt, snapshot range, fields, weight, grid)
!> - singular_values.txt: index, sigma_i, sigma_i/sigma_1, cumulative energy fraction
!> - eigenvalues.txt: one line per mode in amplitude order: Re omega, Im omega, |mu|, arg mu, amplitude
!> - times.txt: the snapshot times used
!> - modes.bin + modes.txt: the leading n_keep exact modes, complex128, stored as real and imaginary parts, column-major (n_space, n_keep); modes.txt gives the shape and the field layout
!>
!> Perturbation export, a directory: perturbation.bin (real64, column-major (n_space, n_snap), the twin difference Q_p - Q_b exactly as fed to DMD), perturbation.txt (shape, fields, grid, domain, dt, source directories) and times.txt. In Python: np.fromfile('perturbation.bin', dtype='<f8').reshape((n_space, n_snap), order='F'). run_lsa analyse --from-export DIR runs the decomposition from it again, so every rank, start-time and stride study can be repeated after the plotfiles have been deleted.
module archive_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use dmd_mod, only: dmd_result_t
   implicit none
   private

   public :: dp
   public :: write_dmd_archive, write_perturbation_export
   public :: read_perturbation_export, make_dir

contains

   !> mkdir -p, with the path quoted for the shell.
   subroutine make_dir(path)
      character(len=*), intent(in) :: path
      call execute_command_line("mkdir -p '" // trim(path) // "'")
   end subroutine make_dir

   subroutine write_dmd_archive(dir, res, times, rate, fields, weight, grid, &
                                skip, stride, n_keep, ierr)
      character(len=*), intent(in) :: dir, fields(:), weight
      type(dmd_result_t), intent(in) :: res
      real(dp), intent(in) :: times(:), rate
      integer, intent(in) :: grid(3), skip, stride, n_keep
      integer, intent(out) :: ierr
      integer :: u, i, k, ios
      real(dp) :: etot, ecum
      real(dp), allocatable :: buf(:)

      ierr = 0
      call make_dir(dir)
      open(newunit=u, file=trim(dir)//'/summary.txt', status='replace', &
           action='write', iostat=ios)
      if (ios /= 0) then
         ierr = 1; return
      end if
      write(u,'(a,es24.16)') 'leading_stationary_growth_rate = ', rate
      write(u,'(a,i0)')      'rank = ', res%rank
      write(u,'(a,i0)')      'rank_requested = ', res%rank_requested
      write(u,'(a,es24.16)') 'reduced_map_residual = ', res%residual
      write(u,'(a,es24.16)') 'dt = ', res%dt
      write(u,'(a,i0)')      'n_snapshots = ', size(times)
      write(u,'(a,es24.16)') 't_first = ', times(1)
      write(u,'(a,es24.16)') 't_last = ', times(size(times))
      write(u,'(a,i0)')      'skip = ', skip
      write(u,'(a,i0)')      'stride = ', stride
      write(u,'(a,i0)')      'amplitude_snapshot = ', res%amp_snapshot
      write(u,'(a,a)')       'weight = ', trim(weight)
      write(u,'(a,i0,1x,i0,1x,i0)') 'grid_nx_ny_ncomp = ', grid
      write(u,'(a)', advance='no') 'fields ='
      do i = 1, size(fields)
         write(u,'(1x,a)', advance='no') trim(fields(i))
      end do
      write(u,'(a)') ''
      close(u)

      open(newunit=u, file=trim(dir)//'/singular_values.txt', &
           status='replace', action='write')
      write(u,'(a)') '# i  sigma_i  sigma_i/sigma_1  cumulative_energy_fraction'
      etot = sum(res%svals**2)
      ecum = 0.0_dp
      do i = 1, size(res%svals)
         ecum = ecum + res%svals(i)**2
         write(u,'(i6,3es24.16)') i, res%svals(i), &
            res%svals(i) / max(res%svals(1), tiny(1.0_dp)), &
            ecum / max(etot, tiny(1.0_dp))
      end do
      close(u)

      open(newunit=u, file=trim(dir)//'/eigenvalues.txt', status='replace', &
           action='write')
      write(u,'(a)') '# mode  Re_omega  Im_omega  abs_mu  arg_mu  amplitude'
      do i = 1, res%n_modes
         write(u,'(i6,5es24.16)') i, real(res%omega(i), dp), &
            aimag(res%omega(i)), abs(res%ritz(i)), &
            atan2(aimag(res%ritz(i)), real(res%ritz(i), dp)), res%amplitude(i)
      end do
      close(u)

      open(newunit=u, file=trim(dir)//'/times.txt', status='replace', &
           action='write')
      do i = 1, size(times)
         write(u,'(es24.16)') times(i)
      end do
      close(u)

      k = min(max(n_keep, 0), res%n_modes)
      open(newunit=u, file=trim(dir)//'/modes.txt', status='replace', &
           action='write')
      write(u,'(a,i0)') 'n_space = ', size(res%modes, 1)
      write(u,'(a,i0)') 'n_modes = ', k
      write(u,'(a,i0,1x,i0,1x,i0)') 'grid_nx_ny_ncomp = ', grid
      write(u,'(a)') 'layout = for each mode: n_space real parts, then '// &
         'n_space imaginary parts; within a part, component-major, then '// &
         'z (ny), then r (nx) fastest'
      write(u,'(a)') "python = m = np.fromfile('modes.bin','<f8')"// &
         ".reshape(n_modes, 2, n_space); phi = m[:,0] + 1j*m[:,1]"
      close(u)
      open(newunit=u, file=trim(dir)//'/modes.bin', status='replace', &
           access='stream', form='unformatted', action='write')
      allocate(buf(size(res%modes, 1)))
      do i = 1, k
         buf = real(res%modes(:, i), dp)
         write(u) buf
         buf = aimag(res%modes(:, i))
         write(u) buf
      end do
      close(u)
   end subroutine write_dmd_archive

   subroutine write_perturbation_export(dir, X, times, dt, fields, grid, &
                                        domain, pert_dir, base_dir, weight, &
                                        ierr)
      character(len=*), intent(in) :: dir, fields(:), pert_dir, base_dir, weight
      real(dp), intent(in) :: X(:,:), times(:), dt, domain(4)
      integer, intent(in) :: grid(3)
      integer, intent(out) :: ierr
      integer :: u, i, ios

      ierr = 0
      call make_dir(dir)
      open(newunit=u, file=trim(dir)//'/perturbation.bin', status='replace', &
           access='stream', form='unformatted', action='write', iostat=ios)
      if (ios /= 0) then
         ierr = 1; return
      end if
      write(u) X
      close(u)

      open(newunit=u, file=trim(dir)//'/perturbation.txt', status='replace', &
           action='write')
      write(u,'(a,i0)') 'n_space = ', size(X, 1)
      write(u,'(a,i0)') 'n_snap = ', size(X, 2)
      write(u,'(a,i0,1x,i0,1x,i0)') 'grid_nx_ny_ncomp = ', grid
      write(u,'(a,4es24.16)') 'domain_lo_r_lo_z_hi_r_hi_z = ', domain
      write(u,'(a,es24.16)') 'dt = ', dt
      write(u,'(a,a)') 'weight = ', trim(weight)
      write(u,'(a,a)') 'perturbed = ', trim(pert_dir)
      write(u,'(a,a)') 'base = ', trim(base_dir)
      write(u,'(a)', advance='no') 'fields ='
      do i = 1, size(fields)
         write(u,'(1x,a)', advance='no') trim(fields(i))
      end do
      write(u,'(a)') ''
      write(u,'(a)') "python = X = np.fromfile('perturbation.bin','<f8')"// &
         ".reshape((n_space, n_snap), order='F')"
      close(u)

      open(newunit=u, file=trim(dir)//'/times.txt', status='replace', &
           action='write')
      do i = 1, size(times)
         write(u,'(es24.16)') times(i)
      end do
      close(u)
   end subroutine write_perturbation_export

   !> Read back an export written above. ierr = 1 missing header, 2 missing data, 3 inconsistent sizes.
   subroutine read_perturbation_export(dir, X, times, dt, grid, ierr, domain, &
                                       weight)
      character(len=*), intent(in) :: dir
      real(dp), allocatable, intent(out) :: X(:,:), times(:)
      real(dp), intent(out) :: dt
      integer, intent(out) :: grid(3), ierr
      real(dp), intent(out), optional :: domain(4)
      character(len=*), intent(out), optional :: weight
      real(dp) :: dom(4)
      character(len=64) :: wt
      integer :: u, ios, nsp, nsn, i, p
      character(len=1024) :: line, key, val

      ierr = 0; nsp = 0; nsn = 0; dt = 0.0_dp; grid = 0; dom = 0.0_dp
      wt = 'none'
      open(newunit=u, file=trim(dir)//'/perturbation.txt', status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 1; return
      end if
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) exit
         p = index(line, '=')
         if (p == 0) cycle
         key = adjustl(line(1:p-1)); val = adjustl(line(p+1:))
         select case (trim(key))
         case ('n_space');           read(val, *, iostat=ios) nsp
         case ('n_snap');            read(val, *, iostat=ios) nsn
         case ('dt');                read(val, *, iostat=ios) dt
         case ('grid_nx_ny_ncomp');  read(val, *, iostat=ios) grid
         case ('domain_lo_r_lo_z_hi_r_hi_z'); read(val, *, iostat=ios) dom
         case ('weight');            wt = trim(val)
         end select
      end do
      close(u)
      if (present(domain)) domain = dom
      if (present(weight)) weight = wt
      if (nsp < 1 .or. nsn < 2) then
         ierr = 3; return
      end if
      allocate(X(nsp, nsn), times(nsn))
      open(newunit=u, file=trim(dir)//'/perturbation.bin', status='old', &
           access='stream', form='unformatted', action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 2; return
      end if
      read(u, iostat=ios) X
      close(u)
      if (ios /= 0) then
         ierr = 3; return
      end if
      open(newunit=u, file=trim(dir)//'/times.txt', status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 2; return
      end if
      do i = 1, nsn
         read(u, *, iostat=ios) times(i)
         if (ios /= 0) then
            ierr = 3; close(u); return
         end if
      end do
      close(u)
   end subroutine read_perturbation_export

end module archive_mod
