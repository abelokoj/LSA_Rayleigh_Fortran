!> make_inputs -- generate matched AMReX decks for the NS-MFP study.
!>
!> Fortran counterpart of make_inputs.py (the Python implementation used as an independent check during development; not part of this repository).
!>
!> Usage: make_inputs single    --kr0 0.7 --eps 1e-3 --outdir runs/decks make_inputs sweep-eps --kr0 0.7 --eps 1e-4,1e-3,1e-2 --outdir runs/decks make_inputs sweep-k   --kr0 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
!>                         --eps 1e-3 --outdir runs/decks
!>
!> Options: --rmax --cells-per-r0 --sigma --rho-l --rho-g --mu-l --mu-g
!>          --n-periods --snapshots --safety --template
program make_inputs
   use make_inputs_mod
   use dispersion_mod, only: dp
   use util_mod, only: nice, read_real, read_int, parse_real_list
   implicit none

   character(len=32)  :: cmd
   character(len=512) :: outdir, template
   real(dp) :: kr0_list(64), eps_list(64)
   integer  :: n_kr0, n_eps
   type(run_plan_t) :: base_p
   integer :: i, j, ierr, ndeck, total_steps
   character(len=512) :: path
   character(len=64)  :: tag_k, tag_e

   call parse_args(cmd, outdir, template, kr0_list, n_kr0, eps_list, n_eps, &
                   base_p, ierr)
   if (ierr /= 0) then
      call usage()
      stop 1
   end if

   ! Quoted so that an output directory containing a space is created as one directory.
   call execute_command_line("mkdir -p '" // trim(outdir) // "'")
   ndeck = 0
   total_steps = 0

   select case (trim(cmd))

   case ('single')
      call emit(kr0_list(1), 0.0_dp, 'unperturbed twin (base state)', .true.)
      call emit(kr0_list(1), eps_list(1), 'perturbed', .false.)

   case ('sweep-eps')
      ! One twin serves every amplitude: the base state is eps-independent.
      call emit(kr0_list(1), 0.0_dp, 'unperturbed twin (base state)', .true.)
      do j = 1, n_eps
         call emit(kr0_list(1), eps_list(j), 'perturbed', .false.)
      end do
      write(*,'(a)') ''
      write(*,'(a)') 'All amplitudes share the single base deck above: the'
      write(*,'(a)') 'twin is independent of eps, so it is run once and reused.'

   case ('sweep-k')
      ! Each wavenumber needs its own twin: the box length changes with k, so the base state is a different discrete problem each time.
      do i = 1, n_kr0
         call emit(kr0_list(i), 0.0_dp, 'unperturbed twin (base state)', .true.)
         call emit(kr0_list(i), eps_list(1), 'perturbed', .false.)
      end do
      write(*,'(a)') ''
      write(*,'(a)') 'Each wavenumber needs its own twin: the box length'
      write(*,'(a)') 'changes with k, so the base state differs each time.'

   case default
      call usage()
      stop 1
   end select

   write(*,'(a)') ''
   write(*,'(a,i0,a,i0)') 'wrote ', ndeck, ' deck(s); total time steps ', &
      total_steps

contains

   subroutine emit(kr0, eps, role, is_base)
      real(dp), intent(in) :: kr0, eps
      character(len=*), intent(in) :: role
      logical, intent(in) :: is_base
      type(run_plan_t) :: p

      p = base_p
      p%kr0 = kr0
      p%eps = eps
      p%role = role
      call plan_run(p)

      write(tag_k,'(a)') trim(nice(kr0))
      if (is_base) then
         path = trim(outdir) // '/inputs.k' // trim(tag_k) // '.base'
      else
         ! Named from the amplitude actually used, which differs from the --eps value when --eps-cells is given.
         write(tag_e,'(a)') trim(nice(p%eps))
         path = trim(outdir) // '/inputs.k' // trim(tag_k) // &
                '.eps' // trim(tag_e)
      end if

      call write_deck(p, trim(template), trim(path), ierr)
      if (ierr /= 0) then
         write(*,'(a,a,a,i0)') 'ERROR writing ', trim(path), ', ierr=', ierr
         stop 2
      end if
      write(*,'(a)') trim(path)
      call describe_plan(p)
      write(*,'(a)') ''
      ndeck = ndeck + 1
      total_steps = total_steps + p%max_step
   end subroutine emit

   subroutine usage()
      write(*,'(a)') 'usage: make_inputs <single|sweep-eps|sweep-k> [options]'
      write(*,'(a)') '  --kr0 V[,V...]      wavenumber(s) k*r0'
      write(*,'(a)') '  --eps V[,V...]      perturbation amplitude(s)'
      write(*,'(a)') '  --outdir DIR        output directory'
      write(*,'(a)') '  --template PATH     deck template'
      write(*,'(a)') '  --rmax V            domain radius (default 4 r0)'
      write(*,'(a)') '  --cells-per-r0 N    radial resolution (default 32)'
      write(*,'(a)') '  --sigma V --rho-l V --rho-g V --mu-l V --mu-g V'
      write(*,'(a)') '  --n-periods V       e-folding times (default 5)'
      write(*,'(a)') '  --snapshots N       target plotfile count'
      write(*,'(a)') '  --safety V          fraction of capillary dt limit'
      write(*,'(a)') '  --n-waves N         wavelengths in the axial box (default 1)'
      write(*,'(a)') '  --eps-cells C       perturbed eps = C * dr (overrides --eps)'
      write(*,'(a)') '  --broadband M       sum of M axial harmonics of the box, eps each'
   end subroutine usage

   subroutine parse_args(cmd, outdir, template, kl, nk, el, ne, p, ierr)
      character(len=*), intent(out) :: cmd, outdir, template
      real(dp), intent(out) :: kl(:), el(:)
      integer, intent(out) :: nk, ne, ierr
      type(run_plan_t), intent(out) :: p

      integer :: na, i
      character(len=512) :: a, v

      ierr = 0
      cmd = ''
      outdir = 'runs/decks'
      template = 'inputs/inputs_rayleigh_template.txt'
      kl(1) = 0.7_dp; nk = 1
      el(1) = 1.0e-3_dp; ne = 1

      na = command_argument_count()
      if (na < 1) then
         ierr = 1; return
      end if
      call get_command_argument(1, cmd)

      i = 2
      do while (i <= na)
         call get_command_argument(i, a)
         select case (trim(a))
         case ('--kr0')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call parse_real_list(v, kl, nk, ierr); if (ierr /= 0) return
         case ('--eps')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call parse_real_list(v, el, ne, ierr); if (ierr /= 0) return
         case ('--outdir')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            outdir = v
         case ('--template')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            template = v
         case ('--rmax')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_real(v, p%rmax, ierr); if (ierr /= 0) return
         case ('--cells-per-r0')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_int(v, p%cells_per_r0, ierr); if (ierr /= 0) return
         case ('--sigma')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_real(v, p%sigma, ierr); if (ierr /= 0) return
         case ('--rho-l')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_real(v, p%rho_l, ierr); if (ierr /= 0) return
         case ('--rho-g')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_real(v, p%rho_g, ierr); if (ierr /= 0) return
         case ('--mu-l')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_real(v, p%mu_l, ierr); if (ierr /= 0) return
         case ('--mu-g')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_real(v, p%mu_g, ierr); if (ierr /= 0) return
         case ('--n-periods')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_real(v, p%n_periods, ierr); if (ierr /= 0) return
         case ('--snapshots')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_int(v, p%target_snapshots, ierr); if (ierr /= 0) return
         case ('--safety')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_real(v, p%safety, ierr); if (ierr /= 0) return
         case ('--n-waves')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_int(v, p%n_waves, ierr); if (ierr /= 0) return
         case ('--broadband')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_int(v, p%broadband, ierr); if (ierr /= 0) return
         case ('--eps-cells')
            call next_val(i, na, v, ierr); if (ierr /= 0) return
            call read_real(v, p%eps_cells, ierr); if (ierr /= 0) return
         case default
            write(*,'(a,a)') 'unknown option: ', trim(a)
            ierr = 2; return
         end select
         i = i + 1
      end do
      ! Values that would give a meaningless deck are refused here rather than producing a division by zero or an empty grid further on.
      if (any(kl(1:nk) <= 0.0_dp) .or. any(el(1:ne) < 0.0_dp) .or. &
          p%rmax <= 0.0_dp .or. p%cells_per_r0 < 1 .or. &
          p%sigma <= 0.0_dp .or. p%rho_l <= 0.0_dp .or. &
          p%safety <= 0.0_dp .or. p%n_periods <= 0.0_dp .or. &
          p%target_snapshots < 1 .or. p%n_waves < 1 .or. &
          p%eps_cells < 0.0_dp .or. p%broadband < 0) then
         write(*,'(a)') 'error: kr0, rmax, cells-per-r0, sigma, rho-l, '// &
            'safety, n-periods and snapshots must be positive, eps >= 0'
         ierr = 4
      end if
   end subroutine parse_args

   subroutine next_val(i, na, v, ierr)
      integer, intent(inout) :: i
      integer, intent(in) :: na
      character(len=*), intent(out) :: v
      integer, intent(out) :: ierr
      ierr = 0
      i = i + 1
      if (i > na) then
         ierr = 3; return
      end if
      call get_command_argument(i, v)
   end subroutine next_val

end program make_inputs
