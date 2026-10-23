!> run_lsa -- NS-MFP + DMD stability analysis driver.
!>
!> Fortran counterpart of run_lsa.py (the Python implementation used as an independent check during development; not part of this repository).
!>
!> Subcommands: analyse    one (perturbed, twin) pair -> growth rate vs theory linearity  amplitude sweep -> proportionality test verdict dispersion wavenumber sweep -> measured vs exact dispersion curve converge   subspace-size convergence of the growth rate
!>
!> Examples: run_lsa analyse --perturbed runs/k0.7_eps0.001 --base runs/k0.7_base \
!>                   --kr0 0.7
!> run_lsa linearity --base runs/k0.7_base --kr0 0.7 \
!>           --case 1e-4=runs/a --case 1e-3=runs/b --case 1e-2=runs/c
!> run_lsa dispersion --root runs --eps 1e-3 --kr0 0.3,0.5,0.7,0.9
!>
!> Subcommands added 2026-10-07 from the recommendations of the September 2026 evidence report: record    admissible-window growth rate of an interface-amplitude record (no plotfiles needed) export    write the matched perturbation sequence Q_p - Q_b as a lossless binary archive and stop Options added to analyse/converge: --archive DIR (save singular values, eigenvalues, residual and leading modes), --export DIR, --from-export DIR (analyse an export instead of plotfiles), --weight volume, --strides, --ranks, --csv PATH.
!>
!> Exit status: 0 success; 1 bad command line; 2 the data could not be read or decomposed; 3 the linearity test did not pass; 4 the analysis ran but produced no usable DMD growth rate (no mode survived the filters), or in dispersion and converge some wavenumber or configuration produced none. A failed norm fit alone is reported as a warning with exit 0, since the DMD rate is the headline result. A non-zero status is what lets a Slurm job or a sweep script notice a failure that is otherwise only printed.
!>
!> Theory columns. "viscous exact" is Chandrasekhar's full viscous dispersion relation (dispersion_mod::growth_rate_viscous_exact), the benchmark of the dissertation; "viscous approx" is the low-Oh Weber-type estimate that was the only viscous reference before 2026-10-07; "two-fluid inviscid" includes the ambient gas density (--rho-g). Departures are measured against the exact viscous rate.
program run_lsa
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use dispersion_mod
   use table_mod
   use dmd_mod
   use nsmfp_mod
   use proportionality_mod
   use util_mod, only: nice, read_real, read_int, parse_real_list, &
                       parse_name_list, is_nan, nan_value, UNSET, &
                       fit_log_linear, argsort_real_asc
   use archive_mod
   use record_mod
   implicit none

   character(len=32)  :: cmd
   character(len=512) :: pert_dir, base_dir, root, out_json
   character(len=64)  :: fields(8)
   character(len=64)  :: plt_prefix
   integer :: nfields
   real(dp) :: kr0_list(64), eps_list(64)
   integer  :: n_kr0, n_eps
   character(len=512) :: case_dirs(32)
   real(dp) :: case_eps(32)
   integer :: n_cases
   integer :: skip, stride, maxsnap, rank
   real(dp) :: tol, t_min, t_max
   real(dp) :: r0, sigma, rho_l, mu_l, rho_g, mu_g, r_wall
   integer :: starts(16), n_starts
   integer :: strides(16), n_strides, ranks(16), n_ranks
   character(len=512) :: archive_dir, export_dir, from_export, csv_path
   character(len=512) :: record_path
   character(len=16)  :: weight
   integer :: n_keep
   character(len=16) :: select_mode
   integer :: rscan_lo, rscan_hi
   real(dp) :: band
   character(len=512) :: krylov_dir
   real(dp) :: a_max, width, rel_tol, starts_t(32)
   integer :: n_starts_t
   integer :: ierr

   call defaults()
   call parse_args(ierr)
   if (ierr /= 0) then
      call usage(); stop 1
   end if

   select case (trim(cmd))
   case ('analyse');    call do_analyse()
   case ('linearity');  call do_linearity()
   case ('dispersion'); call do_dispersion()
   case ('converge');   call do_converge()
   case ('record');     call do_record()
   case ('export');     call do_export()
   case ('theory');     call do_theory()
   case ('drift');      call do_drift()
   case default;        call usage(); stop 1
   end select

contains

   subroutine defaults()
      cmd = ''
      pert_dir = ''; base_dir = ''; root = 'runs'; out_json = ''
      fields(1) = 'x_velocity'; fields(2) = 'y_velocity'; nfields = 2
      plt_prefix = 'nddataPLT'
      kr0_list(1) = 0.7_dp; n_kr0 = 1
      eps_list(1) = 1.0e-3_dp; n_eps = 1
      n_cases = 0
      skip = 0; stride = 1; maxsnap = 0; rank = 0
      tol = 1.0e-10_dp; t_min = UNSET; t_max = UNSET
      r0 = 1.0_dp; sigma = 1.0_dp; rho_l = 1.0_dp; mu_l = 0.02_dp
      rho_g = 1.225e-3_dp
      mu_g = 2.6e-4_dp
      r_wall = 0.0_dp
      starts(1) = 0; n_starts = 1
      strides(1) = 1; n_strides = 1
      ranks(1) = -1; n_ranks = 1            ! -1: use --rank
      archive_dir = ''; export_dir = ''; from_export = ''; csv_path = ''
      record_path = ''
      weight = 'none'
      n_keep = 10
      select_mode = 'dominant'
      rscan_lo = 0; rscan_hi = -1
      band = 0.0_dp
      krylov_dir = ''
      a_max = 0.1_dp; width = 0.8_dp; rel_tol = 0.005_dp
      n_starts_t = 0
   end subroutine defaults

   !> The twin-difference subspace, from plotfiles or from an export.
   !>
   !> With --from-export the columns are selected here with the same --skip, --stride and --max-snapshots semantics that build_snapshot_matrix applies to plotfiles, so a study repeated from an export uses the same snapshots as one run from plotfiles.
   subroutine get_subspace(pd, bd, X, times, dt, grid, domain, e)
      character(len=*), intent(in) :: pd, bd
      real(dp), allocatable, intent(out) :: X(:,:), times(:)
      real(dp), intent(out) :: dt, domain(4)
      integer, intent(out) :: grid(3), e
      real(dp), allocatable :: XA(:,:), TA(:), rc(:)
      integer :: n, i, k, nuse, j, c, off
      integer, allocatable :: cols(:)
      character(len=64) :: wt_in

      domain = 0.0_dp
      if (len_trim(from_export) > 0) then
         call read_perturbation_export(trim(from_export), XA, TA, dt, grid, e, &
                                       domain=domain, weight=wt_in)
         if (e /= 0) then
            write(*,'(a,i0,a,a)') 'ERROR reading export (ierr=', e, ') in ', &
               trim(from_export)
            e = 3000 + e; return
         end if
         n = size(TA)
         allocate(cols(n))
         nuse = 0
         do i = skip + 1, n, max(stride, 1)
            nuse = nuse + 1
            cols(nuse) = i
            if (maxsnap > 0 .and. nuse >= maxsnap) exit
         end do
         if (nuse < 2) then
            e = 4; return
         end if
         allocate(X(size(XA,1), nuse), times(nuse))
         do k = 1, nuse
            X(:, k) = XA(:, cols(k))
            times(k) = TA(cols(k))
         end do
         dt = dt * real(max(stride, 1), dp)
         ! Volume weighting of an unweighted export: rows times sqrt(r), as build_snapshot_matrix does for plotfiles.
         if (trim(weight) == 'volume' .and. trim(wt_in) /= 'volume') then
            allocate(rc(grid(1)))
            do i = 1, grid(1)
               rc(i) = domain(1) + (real(i,dp) - 0.5_dp)*(domain(3) - domain(1))/real(grid(1),dp)
            end do
            do c = 1, grid(3)
               do j = 1, grid(2)
                  off = (c-1)*grid(1)*grid(2) + (j-1)*grid(1)
                  do i = 1, grid(1)
                     X(off+i, :) = X(off+i, :) * sqrt(max(rc(i), 0.0_dp))
                  end do
               end do
            end do
         end if
         write(*,'(a,i0,a,i0,a,es12.5,a,a)') 'subspace (from export): ', &
            size(X,1), ' dof x ', size(X,2), ' snapshots, dt=', dt, '  ', &
            trim(from_export)
      else if (len_trim(krylov_dir) > 0) then
         ! The solver's built-in restarting Krylov driver: columns are dS_k for members skip+1 .. skip+max-snapshots.
         call build_krylov_matrix(trim(krylov_dir), fields(1:nfields), skip + 1, &
              merge(skip + maxsnap, 0, maxsnap > 0), X, times, dt, e, &
              weight=trim(weight), grid=grid, domain=domain)
         if (e /= 0) then
            write(*,'(a,i0,a,a)') 'ERROR reading the Krylov members (ierr=', &
               e, ') in ', trim(krylov_dir)
            e = 4000 + e
         end if
      else
         call build_snapshot_matrix(pd, bd, fields(1:nfields), skip, stride, &
                                    maxsnap, X, times, dt, e, &
                                    verbose=.true., prefix=trim(plt_prefix), &
                                    weight=trim(weight), grid=grid, &
                                    domain=domain)
      end if
      if (e == 0 .and. band > 0.0_dp .and. trim(cmd) /= 'export') &
         call apply_band(X, grid, domain)
   end subroutine get_subspace

   !> --band H: keep only the rows within H radial cells of the undisturbed interface r = r0, in every component and at every axial station.
   !>
   !> Added 2026-10-08 from the first campaign data. With the level set L0101 as the observable this is the narrow-band state of dissertation eq:4_norm (|phi_bar| < c dx): the raw level-set difference is dominated, away from the interface, by components the redistancing creates that do not grow, and the decomposition then ranks a neutral mode first; the velocity difference carries a second structure of 1-2 % that makes the rate move by several per cent with the rank. In the band |r - r0| < 2 dr the level-set difference gave rank-insensitive rates on meas08, meas16 and sw080_16 (interquartile range over ranks 2-12 of 0.1-0.5 %).
   subroutine apply_band(X, grid, domain)
      real(dp), allocatable, intent(inout) :: X(:,:)
      integer, intent(in) :: grid(3)
      real(dp), intent(in) :: domain(4)
      real(dp), allocatable :: Y(:,:)
      logical, allocatable :: keep(:)
      real(dp) :: drr, rc
      integer :: i, j, c, row, n
      if (product(grid) /= size(X, 1) .or. domain(3) <= domain(1)) then
         write(*,'(a)') 'WARNING: --band ignored: grid or domain unknown'
         return
      end if
      drr = (domain(3) - domain(1)) / real(grid(1), dp)
      allocate(keep(size(X, 1)))
      row = 0
      do c = 1, grid(3)
         do j = 1, grid(2)
            do i = 1, grid(1)
               row = row + 1
               rc = domain(1) + (real(i, dp) - 0.5_dp) * drr
               keep(row) = abs(rc - r0) < band * drr
            end do
         end do
      end do
      n = count(keep)
      if (n == 0) then
         write(*,'(a)') 'WARNING: --band selects no rows; ignored'
         return
      end if
      allocate(Y(n, size(X, 2)))
      do j = 1, size(X, 2)
         Y(:, j) = pack(X(:, j), keep)
      end do
      call move_alloc(Y, X)
      write(*,'(a,f6.2,a,i0,a)') 'narrow band |r - r0| < ', band, &
         ' dr: ', n, ' rows kept'
   end subroutine apply_band

   !> Build the subspace for one twin pair and extract the growth rate.
   !>
   !> Two independent estimates are produced deliberately: rate_fit - slope of log|Q'(t)|, assumption-light, sees only the
   !>              dominant mode;
   !> rate_dmd - the leading stationary DMD eigenvalue, which also yields
   !>              the rest of the spectrum.
   !> They should agree once the transient has shed; disagreement usually means the fit window still contains the non-modal transient.
   subroutine analyse_pair(pd, bd, rate_dmd, rate_fit, r2, res, ok)
      character(len=*), intent(in) :: pd, bd
      real(dp), intent(out) :: rate_dmd, rate_fit, r2
      type(dmd_result_t), intent(out) :: res
      logical, intent(out) :: ok

      real(dp), allocatable :: X(:,:), times(:), norms(:)
      real(dp) :: dt, intercept, tlo, domain(4)
      integer :: e, n, grid(3)
      character(len=512) :: adir

      ok = .false.
      rate_dmd = nan_value(); rate_fit = nan_value(); r2 = 0

      call get_subspace(pd, bd, X, times, dt, grid, domain, e)
      if (e /= 0) then
         write(*,'(a,i0)') 'ERROR building subspace, ierr=', e
         call explain_nsmfp_error(e)
         return
      end if

      call perturbation_norms(X, norms)
      if (maxval(abs(norms)) <= 0.0_dp) then
         write(*,'(a)') 'ERROR: the perturbation field is identically zero.'
         write(*,'(a)') 'The two runs produced identical output; check that'
         write(*,'(a)') 'the perturbed deck really has ns.radblob /= 0.'
         return
      end if

      ! With velocity fields both twins start from rest and Q'(0) = 0, so with --skip 0 the first snapshot carries no signal. dmd_compute then ranks the modes by their contribution to the snapshot of largest norm instead, but the startup transient is still inside the subspace; the remedy is --skip, and saying so here is cheaper than a wrong rate.
      if (norms(1) <= 1.0e-8_dp * maxval(norms)) then
         write(*,'(a)') 'WARNING: the first snapshot of the subspace is '// &
            'negligible (Q''(0) = 0, as for velocity fields from rest).'
         write(*,'(a,i0,a)') '         --skip ', skip, ' leaves the '// &
            'startup transient in the subspace; use --skip > 0 '// &
            '(scripts/run_nsmfp.sh uses 60).'
      end if

      n = size(times)
      tlo = t_min
      if (.not. (tlo > UNSET/2)) tlo = times(max(1, n/3))
      call fit_exponential_growth(times, norms, tlo, t_max, &
                                  rate_fit, intercept, r2, e)
      if (e /= 0) then
         rate_fit = nan_value()
         if (e == 21) then
            write(*,'(a)') 'WARNING: growth-rate fit not done: --t-max '// &
               'does not exceed the start of the window'
         else
            write(*,'(a)') 'WARNING: growth-rate fit failed (fewer than '// &
               'two positive samples in the window)'
         end if
      end if

      call dmd_compute(X, dt, rank, tol, res, e)
      if (e /= 0) then
         write(*,'(a,i0)') 'ERROR in DMD, ierr=', e
         if (e == 5) write(*,'(a)') '  every snapshot is zero'
         return
      end if
      if (rank > 0 .and. res%rank < rank) write(*,'(a,i0,a,i0,a)') &
         'note: --rank ', rank, ' exceeds the numerical rank of the '// &
         'data; ', res%rank, ' modes were used'
      if (res%amp_snapshot > 1) write(*,'(a,i0,a)') 'note: mode '// &
         'amplitudes are contributions to snapshot ', res%amp_snapshot, &
         ' (the first is negligible)'
      rate_dmd = pick_rate(res)
      if (rscan_hi >= rscan_lo .and. rscan_lo > 0) call rank_scan(X, dt)
      if (is_nan(rate_dmd)) then
         write(*,'(a)') 'WARNING: no DMD mode passed the stationary and '// &
            'amplitude filters; no DMD growth rate'
      end if
      write(*,'(a,es10.3,a)') 'reduced-map residual ', res%residual, &
         '  (||Y - U S~ U^T X|| / ||Y||; small means a linear sequence)'

      ! Archives, one directory per pair when several pairs are analysed (dispersion): the archive name gets the run directory's last component.
      if (len_trim(archive_dir) > 0) then
         adir = trim(archive_dir)
         if (trim(cmd) == 'dispersion') adir = trim(archive_dir)//'/'//trim(basename(pd))
         call write_dmd_archive(trim(adir), res, times, rate_dmd, &
              fields(1:nfields), trim(weight), grid, skip, stride, n_keep, e)
         if (e == 0) then
            write(*,'(a,a)') 'DMD archive written to ', trim(adir)
         else
            write(*,'(a,a)') 'WARNING: could not write DMD archive in ', trim(adir)
         end if
      end if
      if (len_trim(export_dir) > 0 .and. len_trim(from_export) == 0) then
         adir = trim(export_dir)
         if (trim(cmd) == 'dispersion') adir = trim(export_dir)//'/'//trim(basename(pd))
         call write_perturbation_export(trim(adir), X, times, dt, &
              fields(1:nfields), grid, domain, pd, bd, trim(weight), e)
         if (e == 0) write(*,'(a,a)') 'perturbation export written to ', trim(adir)
      end if
      ok = .true.
   end subroutine analyse_pair

   subroutine explain_nsmfp_error(e)
      integer, intent(in) :: e
      select case (e)
      case (1); write(*,'(a)') '  no nddataPLT* directories under the perturbed run'
      case (2); write(*,'(a)') '  no nddataPLT* directories under the base run'
      case (3); write(*,'(a)') &
         '  fewer than 2 matching steps; check amr.plot_int'
      case (4); write(*,'(a)') '  too few snapshots after skip/stride'
      case (5); write(*,'(a)') &
         '  incompatible plotfiles (grid/time/domain)'
      case (6); write(*,'(a)') '  --skip must be >= 0 and --stride >= 1'
      case (1013, 2013); write(*,'(a)') &
         '  the boxes of a plotfile do not cover the domain'
      case (1602, 2602); write(*,'(a)') &
         '  a plotfile box lies outside the domain'
      case (1612, 2612); write(*,'(a)') &
         '  a plotfile holds single-precision data; set fab.precision = DOUBLE'
      case default
         if (e > 2000) then
            write(*,'(a,i0)') '  failed reading a BASE plotfile, code ', e-2000
         else if (e > 1000) then
            write(*,'(a,i0)') &
               '  failed reading a PERTURBED plotfile, code ', e-1000
         end if
      end select
      if (e == 12) write(*,'(a)') '  snapshot spacing is not uniform'
   end subroutine explain_nsmfp_error

   !> Text for a rate, "n/a" when it is NaN (no result).
   function rate_text(x) result(t)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: t
      if (is_nan(x)) then
         t = 'n/a'
      else
         t = trim(fmt_f(x, 6))
      end if
   end function rate_text

   subroutine report(kr0, rate_dmd, rate_fit, r2)
      real(dp), intent(in) :: kr0, rate_dmd, rate_fit, r2
      real(dp) :: inv, vis, visx, two, twov, spread
      type(table_t) :: t, c
      inv = growth_rate(kr0/r0, r0, sigma, rho_l)
      vis = growth_rate_viscous(kr0/r0, r0, sigma, rho_l, mu_l)
      visx = growth_rate_viscous_exact(kr0/r0, r0, sigma, rho_l, mu_l)
      two = growth_rate_two_fluid_inviscid(kr0/r0, r0, sigma, rho_l, rho_g)
      twov = growth_rate_two_fluid_viscous(kr0/r0, r0, sigma, rho_l, mu_l, &
                                           rho_g, mu_g)

      write(*,'(a)') ''
      call t%init([cel('quantity '), cel('growth rate '), cel('note ')], 'lrl')
      call t%row([cel('DMD: leading stationary mode'), cel(rate_text(rate_dmd)), &
         cel(' ')])
      call t%row([cel('log-norm slope fit '), cel(rate_text(rate_fit)), &
         cel('r^2 = '//fmt_f(r2, 5))])
      call t%row([cel('viscous exact (Chandrasekhar)'), cel(rate_text(visx)), &
         cel('theory, benchmark ')])
      call t%row([cel('viscous approx (Weber) '), cel(fmt_f(vis, 6)), &
         cel('theory ')])
      call t%row([cel('exact inviscid (Rayleigh) '), cel(fmt_f(inv, 6)), cel('theory ')])
      call t%row([cel('two-fluid inviscid '), cel(fmt_f(two, 6)), &
         cel('theory, rho_g = '//trim(fmt_e(rho_g)))])
      call t%row([cel('two-fluid viscous (Tomotika) '), cel(rate_text(twov)), &
         cel('theory, rho_g, mu_g = '//trim(fmt_e(rho_g))//', '// &
         trim(fmt_e(mu_g)))])
      call t%emit('Rayleigh-Plateau growth rate at k r0 = '//fmt_f(kr0, 4))

      call c%init([cel('comparison '), cel('departure '), cel('verdict ')], 'lrl')
      if (.not. is_nan(rate_dmd)) then
         if (visx > 0.0_dp) call c%row([cel('measured vs viscous exact '), &
            cel(fmt_pct((rate_dmd-visx)/visx)), cel(' ')])
         if (inv > 0.0_dp) call c%row([cel('measured vs inviscid '), &
            cel(fmt_pct((rate_dmd-inv)/inv)), cel(' ')])
         if (twov > 0.0_dp) call c%row([cel('measured vs two-fluid viscous '), &
            cel(fmt_pct((rate_dmd-twov)/twov)), cel(' ')])
      end if
      if (.not. is_nan(rate_dmd) .and. .not. is_nan(rate_fit) .and. &
          abs(rate_dmd) > 0.0_dp) then
         spread = abs(rate_dmd-rate_fit)/abs(rate_dmd)
         if (spread < 0.05_dp) then
            call c%row([cel('DMD vs norm fit '), cel(fmt_pct(spread)), &
               cel('consistent ')])
         else
            call c%row([cel('DMD vs norm fit '), cel(fmt_pct(spread)), &
               cel('DISAGREE: check fit window ')])
         end if
      end if
      call c%emit('')
   end subroutine report

   subroutine do_analyse()
      real(dp) :: rd, rf, r2
      type(dmd_result_t) :: res
      logical :: ok

      call analyse_pair(trim(pert_dir), trim(base_dir), rd, rf, r2, res, ok)
      if (.not. ok) stop 2
      write(*,'(a)') ''
      call dmd_summary(res, 10)
      call report(kr0_list(1), rd, rf, r2)
      if (len_trim(out_json) > 0) call save_json(out_json, kr0_list(1), rd, &
         rf, r2)
      call dmd_free(res)
      if (is_nan(rd)) stop 4
   end subroutine do_analyse

   !> JSON number, or null for NaN, which JSON cannot represent.
   function jnum(x) result(t)
      real(dp), intent(in) :: x
      character(len=32) :: t
      if (is_nan(x)) then
         t = 'null'
      else
         write(t,'(es20.10)') x
         t = adjustl(t)
      end if
   end function jnum

   subroutine save_json(path, kr0, rd, rf, r2)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: kr0, rd, rf, r2
      integer :: u, ios
      open(newunit=u, file=trim(path), status='replace', action='write', &
           iostat=ios)
      if (ios /= 0) then
         write(*,'(a,a)') 'ERROR: cannot write ', trim(path)
         stop 2
      end if
      write(u,'(a)') '{'
      write(u,'(a,f12.6,a)') '  "kr0": ', kr0, ','
      write(u,'(a,a,a)') '  "growth_rate_dmd": ', trim(jnum(rd)), ','
      write(u,'(a,a,a)') '  "growth_rate_fit": ', trim(jnum(rf)), ','
      write(u,'(a,a,a)') '  "r_squared": ', trim(jnum(r2)), ','
      write(u,'(a,a,a)') '  "growth_rate_exact_inviscid": ', &
         trim(jnum(growth_rate(kr0/r0, r0, sigma, rho_l))), ','
      write(u,'(a,a,a)') '  "growth_rate_viscous_exact": ', &
         trim(jnum(growth_rate_viscous_exact(kr0/r0, r0, sigma, rho_l, mu_l))), ','
      write(u,'(a,a,a)') '  "growth_rate_two_fluid_inviscid": ', &
         trim(jnum(growth_rate_two_fluid_inviscid(kr0/r0, r0, sigma, rho_l, &
         rho_g))), ','
      write(u,'(a,a,a)') '  "growth_rate_two_fluid_viscous": ', &
         trim(jnum(growth_rate_two_fluid_viscous(kr0/r0, r0, sigma, rho_l, &
         mu_l, rho_g, mu_g))), ','
      write(u,'(a,a)') '  "growth_rate_viscous_estimate": ', &
         trim(jnum(growth_rate_viscous(kr0/r0, r0, sigma, rho_l, mu_l)))
      write(u,'(a)') '}'
      close(u)
      write(*,'(a,a)') 'wrote ', trim(path)
   end subroutine save_json

   subroutine do_linearity()
      type(amplitude_run_t), allocatable :: runs(:)
      type(linearity_report_t) :: rep
      real(dp), allocatable :: X(:,:), times(:), norms(:)
      real(dp) :: dt, slope, intercept, r2, tlo
      integer :: i, e, n

      if (n_cases < 2) then
         write(*,'(a)') 'linearity needs at least two --case EPS=DIR options'
         stop 1
      end if
      allocate(runs(n_cases))
      do i = 1, n_cases
         write(*,'(a,es10.3,a,a)') '[eps=', case_eps(i), '] reading ', &
            trim(case_dirs(i))
         call build_snapshot_matrix(trim(case_dirs(i)), trim(base_dir), &
              fields(1:nfields), skip, stride, maxsnap, X, times, dt, e, &
              verbose=.false., prefix=trim(plt_prefix))
         if (e /= 0) then
            write(*,'(a,i0)') 'ERROR building subspace, ierr=', e
            call explain_nsmfp_error(e)
            stop 2
         end if
         call perturbation_norms(X, norms)
         n = size(times)
         tlo = t_min
         if (.not. (tlo > UNSET/2)) tlo = times(1)
         call fit_exponential_growth(times, norms, tlo, t_max, slope, &
                                     intercept, r2, e)
         ! A failed fit is recorded as NaN, which assess_linearity leaves out of the rate statistics; it was previously recorded as a slope of 0, which entered the mean and the spread as if it were a measurement.
         if (e /= 0) then
            write(*,'(a,es10.3,a)') 'WARNING: growth-rate fit failed for '// &
               'eps=', case_eps(i), '; excluded from the rate test'
            slope = nan_value()
         end if
         runs(i)%eps = case_eps(i)
         runs(i)%growth_rate = slope
         runs(i)%r_squared = r2
         ! --t-max also limits the proportionality comparison, not only the fit: beyond the linear window of the largest amplitude the scaled norms are expected to diverge, and comparing them there gives a spurious NONLINEAR verdict.
         if (t_max > UNSET/2.0_dp) n = count(times <= t_max)
         allocate(runs(i)%times(n), runs(i)%norms(n))
         runs(i)%times = times(1:n)
         runs(i)%norms = norms(1:n)
         deallocate(X, times, norms)
      end do

      call assess_linearity(runs, rep, e)
      if (e /= 0) then
         write(*,'(a,i0)') 'ERROR assessing linearity, ierr=', e
         if (e == 5) write(*,'(a)') '  the runs share fewer than two '// &
            'snapshot times; check that they used the same ns.fixed_dt '// &
            'and plot interval'
         stop 2
      end if
      write(*,'(a)') ''
      call print_report(runs, rep)
      write(*,'(a)') ''
      write(*,'(a,f10.6)') 'exact inviscid growth rate : ', &
         growth_rate(kr0_list(1)/r0, r0, sigma, rho_l)
      write(*,'(a,f10.6)') 'viscous exact (Chandrasekhar): ', &
         growth_rate_viscous_exact(kr0_list(1)/r0, r0, sigma, rho_l, mu_l)
      write(*,'(a,f10.6)') 'viscous approx (Weber)     : ', &
         growth_rate_viscous(kr0_list(1)/r0, r0, sigma, rho_l, mu_l)
      write(*,'(a,es12.3)') 'recommended amplitude      : ', &
         recommend_amplitude(runs, rep)

      ! Non-zero exit so a driving script notices the sweep is unusable.
      if (rep%verdict /= VERDICT_LINEAR) stop 3
   end subroutine do_linearity

   subroutine do_dispersion()
      real(dp) :: rd, rf, r2, inv, vis, best_k, best_s, kx, sx
      type(dmd_result_t) :: res
      logical :: ok
      integer :: i, nrow, nfail
      character(len=512) :: pd, bd
      character(len=64) :: tk, te
      real(dp) :: rows_k(64), rows_m(64), rows_f(64)

      nrow = 0; nfail = 0
      do i = 1, n_kr0
         write(tk,'(a)') trim(nice(kr0_list(i)))
         write(te,'(a)') trim(nice(eps_list(1)))
         pd = trim(root)//'/k'//trim(tk)//'_eps'//trim(te)
         bd = trim(root)//'/k'//trim(tk)//'_base'
         write(*,'(a,f6.3,a)') '[k r0=', kr0_list(i), '] reading ...'
         call analyse_pair(trim(pd), trim(bd), rd, rf, r2, res, ok)
         if (.not. ok .or. is_nan(rd)) then
            write(*,'(a,f6.3,a)') '[k r0=', kr0_list(i), '] FAILED, skipping'
            nfail = nfail + 1
            if (ok) call dmd_free(res)
            cycle
         end if
         nrow = nrow + 1
         rows_k(nrow) = kr0_list(i)
         rows_m(nrow) = rd
         rows_f(nrow) = rf
         call dmd_free(res)
      end do

      if (nrow == 0) then
         write(*,'(a)') 'no runs could be analysed'
         stop 2
      end if

      write(*,'(a)') ''
      write(*,'(a8,a14,a14,a14,a14,a14,a12)') 'k r0', 'DMD', 'norm fit', &
         'inviscid', 'visc exact', 'visc approx', 'err vs exact'
      best_s = -huge(1.0_dp); best_k = 0
      do i = 1, nrow
         inv = growth_rate(rows_k(i)/r0, r0, sigma, rho_l)
         vis = growth_rate_viscous_exact(rows_k(i)/r0, r0, sigma, rho_l, mu_l)
         if (vis > 0.0_dp) then
            write(*,'(f8.3,f14.6,f14.6,f14.6,f14.6,f14.6,f11.2,a)') rows_k(i), &
               rows_m(i), rows_f(i), inv, vis, &
               growth_rate_viscous(rows_k(i)/r0, r0, sigma, rho_l, mu_l), &
               100.0_dp*(rows_m(i)-vis)/vis, '%'
         else
            write(*,'(f8.3,f14.6,f14.6,f14.6,f14.6,f14.6,a12)') rows_k(i), &
               rows_m(i), rows_f(i), inv, vis, &
               growth_rate_viscous(rows_k(i)/r0, r0, sigma, rho_l, mu_l), &
               '   n/a'
         end if
         if (rows_m(i) > best_s) then
            best_s = rows_m(i); best_k = rows_k(i)
         end if
      end do

      call most_unstable_wavenumber(r0, sigma, rho_l, kx, sx)
      write(*,'(a)') ''
      write(*,'(a,f7.3,a,f10.6,a)') 'most unstable sampled k r0 : ', best_k, &
         '  (sigma = ', best_s, ')'
      write(*,'(a,f7.3,a,f10.6,a)') 'theoretical peak (inviscid): ', kx*r0, &
         '  (sigma = ', sx, ')'
      if (nfail > 0) then
         write(*,'(a,i0,a)') 'WARNING: ', nfail, ' wavenumber(s) failed '// &
            'and are missing from the table'
         stop 4
      end if
   end subroutine do_dispersion

   !> Sensitivity of the DMD growth rate to the snapshot count, the start of the record, the sampling interval and the rank.
   !>
   !> Sizes 50, 100, 200, 400, 800, 1600 and all, from every offset in --starts (default 0), for every sub-sampling factor in --strides (default 1; stride s doubles nothing but the interval, dt -> s dt, the test of dissertation section 5.7) and every rank in --ranks (default: --rank). This is the set of variations the September 2026 report asked for before a DMD eigenvalue is quoted (its section 11.1). Each line also gives the reduced-map residual. --csv PATH writes the table for plotting.
   subroutine do_converge()
      real(dp), allocatable :: X(:,:), times(:), Xs(:,:)
      real(dp) :: dt, g, vis, lo, hi, mean, domain(4), dts
      type(dmd_result_t) :: res
      integer :: e, i, n, sizes(7), ns, cnt, is, st, js, ir, rk, grid(3), u
      integer :: sd, nsub, k
      real(dp), allocatable :: rates(:)

      call get_subspace(trim(pert_dir), trim(base_dir), X, times, dt, grid, &
                        domain, e)
      if (e /= 0) then
         write(*,'(a,i0)') 'ERROR building subspace, ierr=', e
         call explain_nsmfp_error(e); stop 2
      end if
      vis = growth_rate_viscous_exact(kr0_list(1)/r0, r0, sigma, rho_l, mu_l)
      allocate(rates(size(sizes) * n_starts * n_strides * n_ranks))
      u = 0
      if (len_trim(csv_path) > 0) then
         open(newunit=u, file=trim(csv_path), status='replace', action='write')
         write(u,'(a)') 'stride,dt,start,n_snap,rank_requested,rank,'// &
            'growth,err_vs_exact_pct,residual'
      end if

      write(*,'(a)') ''
      write(*,'(a7,a8,a9,a7,a7,a14,a12,a12)') 'stride', 'start', 'n_snap', &
         'rank', 'r_used', 'growth', 'err exact', 'residual'
      cnt = 0
      do js = 1, n_strides
         sd = max(strides(js), 1)
         nsub = (size(X, 2) - 1) / sd + 1
         allocate(Xs(size(X,1), nsub))
         do k = 1, nsub
            Xs(:, k) = X(:, 1 + (k-1)*sd)
         end do
         dts = dt * real(sd, dp)
         n = nsub
         sizes = [50, 100, 200, 400, 800, 1600, n]
         do is = 1, n_starts
            ! --starts is in snapshots of the full sequence, so that a given start is the same instant for every stride.
            st = starts(is) / sd
            do i = 1, size(sizes)
               ns = sizes(i)
               if (i < size(sizes) .and. ns >= n - st) cycle   ! "all" covers it
               ns = min(ns, n - st)
               if (ns < 3 .or. st < 0) cycle
               do ir = 1, n_ranks
                  rk = ranks(ir)
                  if (rk < 0) rk = rank
                  call dmd_compute(Xs(:, st+1:st+ns), dts, rk, tol, res, e)
                  if (e /= 0) cycle
                  g = pick_rate(res)
                  if (is_nan(g)) then
                     write(*,'(i7,i8,i9,i7,i7,a14,a12,es12.3)') sd, starts(is), ns, &
                        rk, res%rank, 'n/a', '', res%residual
                  else
                     cnt = cnt + 1
                     rates(cnt) = g
                     write(*,'(i7,i8,i9,i7,i7,f14.6,f11.3,a,es12.3)') sd, starts(is), &
                        ns, rk, res%rank, g, 100.0_dp*(g-vis)/max(vis,tiny(1.0_dp)), &
                        '%', res%residual
                  end if
                  if (u /= 0) write(u,'(i0,a,es16.8,a,i0,a,i0,a,i0,a,i0,'// &
                     'a,es20.12,a,es14.6,a,es14.6)') &
                     sd, ',', dts, ',', starts(is), ',', ns, ',', rk, ',', res%rank, ',', &
                     g, ',', 100.0_dp*(g-vis)/max(vis,tiny(1.0_dp)), ',', res%residual
                  call dmd_free(res)
               end do
            end do
         end do
         deallocate(Xs)
      end do
      if (u /= 0) then
         close(u)
         write(*,'(a,a)') 'wrote ', trim(csv_path)
      end if

      if (cnt > 1) then
         lo = minval(rates(1:cnt)); hi = maxval(rates(1:cnt))
         mean = sum(rates(1:cnt))/real(cnt,dp)
         write(*,'(a)') ''
         write(*,'(a,es12.4,a,f7.2,a)') 'spread over configurations: ', &
            hi-lo, '  (', 100.0_dp*(hi-lo)/max(abs(mean),tiny(1.0_dp)), &
            '% of the mean)'
         write(*,'(a)') 'Converged if this spread is small compared with '// &
            'the difference from theory.'
      else if (cnt == 0) then
         write(*,'(a)') 'no configuration produced a growth rate'
         stop 4
      end if
   end subroutine do_converge

   !> The DMD growth rate under the selection rule of --select: 'dominant' (default since 2026-10-08), the largest-amplitude growing near-stationary mode; 'max', the former rule, the largest growth rate among stationary modes passing the amplitude filter.
   function pick_rate(r) result(g)
      type(dmd_result_t), intent(in) :: r
      real(dp) :: g
      if (trim(select_mode) == 'max') then
         g = leading_growth_rate(r)
      else
         g = dominant_growth_rate(r)
      end if
   end function pick_rate

   !> --rank-scan LO:HI: the selected growth rate at every rank from LO to HI on the same snapshots, with the median and the interquartile range. On the first campaign data the rate moved by several per cent from one rank to the next, so a single rank cannot be quoted; the median over a range of ranks with its spread is the reportable value, and a large spread says the decomposition is not converged for that window.
   subroutine rank_scan(X, dt)
      real(dp), intent(in) :: X(:,:), dt
      type(dmd_result_t) :: r
      real(dp), allocatable :: g(:), gs(:)
      integer, allocatable :: idx(:)
      integer :: k, e, n, u
      real(dp) :: med, q1, q3, visx
      allocate(g(rscan_hi - rscan_lo + 1))
      n = 0
      write(*,'(a)') ''
      write(*,'(a,i0,a,i0)') 'rank scan ', rscan_lo, ' .. ', rscan_hi
      do k = rscan_lo, rscan_hi
         call dmd_compute(X, dt, k, tol, r, e)
         if (e /= 0) cycle
         ! Above the numerical rank the decomposition no longer changes: record the capped rank once and stop.
         if (r%rank < k .and. n > 0) then
            call dmd_free(r); exit
         end if
         n = n + 1
         g(n) = pick_rate(r)
         write(*,'(a,i4,a,a,a,es10.3)') '  rank', r%rank, '  rate ', rate_text(g(n)), &
            '  residual ', r%residual
         if (r%rank < k) then
            call dmd_free(r); exit             ! numerical rank reached
         end if
         call dmd_free(r)
      end do
      gs = pack(g(1:n), .not. is_nan(g(1:n)))
      if (size(gs) == 0) then
         write(*,'(a)') '  no rank gave a growing near-stationary mode'; return
      end if
      allocate(idx(size(gs)))
      call argsort_real_asc(gs, idx)
      gs = gs(idx)
      med = quant(gs, 0.5_dp); q1 = quant(gs, 0.25_dp); q3 = quant(gs, 0.75_dp)
      visx = growth_rate_viscous_exact(kr0_list(1)/r0, r0, sigma, rho_l, mu_l)
      write(*,'(a,f10.6,a,es10.3,a,i0,a,i0,a)') '  median rate ', med, &
         '   IQR ', q3 - q1, '   (', size(gs), ' of ', n, ' ranks gave a rate)'
      if (visx > 0.0_dp) write(*,'(a,f8.3,a)') '  median vs viscous exact ', &
         100.0_dp*(med - visx)/visx, ' %'
      if (len_trim(archive_dir) > 0) then
         call make_dir(trim(archive_dir))
         open(newunit=u, file=trim(archive_dir)//'/rank_scan.txt', status='replace', action='write')
         write(u,'(a,i0,a,i0)') '# rank scan ', rscan_lo, ' .. ', rscan_hi
         write(u,'(a,es24.16)') 'median_rate = ', med
         write(u,'(a,es24.16)') 'iqr = ', q3 - q1
         write(u,'(a,i0)') 'ranks_with_rate = ', size(gs)
         write(u,'(a,es24.16)') 'visc_exact = ', visx
         close(u)
      end if
   end subroutine rank_scan

   !> Linear-interpolated quantile p of a sorted vector.
   pure function quant(v, p) result(q)
      real(dp), intent(in) :: v(:), p
      real(dp) :: q, h
      integer :: lo
      h = p * real(size(v) - 1, dp) + 1.0_dp
      lo = min(max(int(h), 1), size(v))
      q = v(lo) + (h - real(lo, dp)) * (v(min(lo + 1, size(v))) - v(lo))
   end function quant

   !> Last component of a path, ignoring a trailing slash.
   function basename(path) result(b)
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: b
      integer :: n, p
      n = len_trim(path)
      if (n > 1) then
         if (path(n:n) == '/') n = n - 1
      end if
      p = index(path(1:n), '/', back=.true.)
      b = path(p+1:n)
   end function basename

   !> Write the matched perturbation sequence and stop; no decomposition.
   subroutine do_export()
      real(dp), allocatable :: X(:,:), times(:)
      real(dp) :: dt, domain(4)
      integer :: e, grid(3)
      if (len_trim(export_dir) == 0) then
         write(*,'(a)') 'export needs --export DIR'; stop 1
      end if
      call get_subspace(trim(pert_dir), trim(base_dir), X, times, dt, grid, &
                        domain, e)
      if (e /= 0) then
         write(*,'(a,i0)') 'ERROR building subspace, ierr=', e
         call explain_nsmfp_error(e); stop 2
      end if
      call write_perturbation_export(trim(export_dir), X, times, dt, &
           fields(1:nfields), grid, domain, trim(pert_dir), trim(base_dir), &
           trim(weight), e)
      if (e /= 0) then
         write(*,'(a)') 'ERROR writing the export'; stop 2
      end if
      write(*,'(a,a,a,i0,a,i0,a)') 'wrote ', trim(export_dir), ' (', &
         size(X,1), ' dof x ', size(X,2), ' snapshots)'
   end subroutine do_export

   !> Admissible-window growth rate of one interface-amplitude record (record_mod). Needs no plotfiles, so it runs on the archived records as well as on new ones.
   subroutine do_record()
      real(dp), allocatable :: t(:), a(:), c(:), r(:)
      type(record_window_t) :: win
      real(dp) :: visx, inv, twov, b, sk, r2k, rmin, rmax_
      integer :: e, nout, i, u, m, nk
      logical :: bnd
      integer :: iend

      if (len_trim(record_path) == 0) then
         write(*,'(a)') 'record needs --record FILE'; stop 1
      end if
      call read_amplitude_record(trim(record_path), t, a, e)
      if (e /= 0) then
         write(*,'(a,i0,a,a)') 'ERROR reading record (ierr=', e, '): ', &
            trim(record_path)
         stop 2
      end if
      call analyse_record(t, a, r0, a_max, width, rel_tol, t_min, t_max, &
                          win, c, r, nout)
      visx = growth_rate_viscous_exact(kr0_list(1)/r0, r0, sigma, rho_l, mu_l)
      inv = growth_rate(kr0_list(1)/r0, r0, sigma, rho_l)
      twov = growth_rate_two_fluid_viscous(kr0_list(1)/r0, r0, sigma, rho_l, &
                                           mu_l, rho_g, mu_g)

      write(*,'(a,a)') 'record: ', trim(record_path)
      write(*,'(a,i0,a,f10.4,a,f10.4)') '  samples ', size(t), '   t = ', &
         t(1), ' .. ', t(size(t))
      write(*,'(a,f6.3,a,a)') '  amplitude bound a/r0 <= ', a_max, &
         merge('   (reached)    ', '   (not reached)', win%bounded)
      if (win%ierr /= 0) then
         write(*,'(a,i0)') 'ERROR: no fit in the window, ierr=', win%ierr
         stop 4
      end if
      write(*,'(a,f10.4,a,f10.4,a,i0,a,f6.2,a)') '  window t = ', win%t_lo, &
         ' .. ', win%t_hi, '  (', win%n, ' samples, ', win%efolds, ' e-folds)'
      if (.not. win%settled) write(*,'(a)') '  WARNING: the local growth '// &
         'rate never settled to within the tolerance; the window is the '// &
         'second half of the linear part and the rate is a finite-window '// &
         'value, not an asymptotic one'
      if (abs(win%efolds) < 0.5_dp) write(*,'(a)') '  WARNING: the window spans '// &
         'less than half an e-fold; the rate is poorly constrained'
      write(*,'(a,f10.6,a,es10.3)') '  growth rate s = ', win%rate, &
         '   1 - r^2 = ', 1.0_dp - win%r2
      if (visx > 0.0_dp) write(*,'(a,f10.6,a,f8.3,a)') &
         '  viscous exact  = ', visx, '   (', &
         100.0_dp*(win%rate-visx)/visx, ' %)'
      if (inv > 0.0_dp) write(*,'(a,f10.6,a,f8.3,a)') &
         '  inviscid       = ', inv, '   (', 100.0_dp*(win%rate-inv)/inv, ' %)'
      if (twov > 0.0_dp) write(*,'(a,f10.6,a,f8.3,a)') &
         '  two-fluid visc = ', twov, '   (', 100.0_dp*(win%rate-twov)/twov, &
         ' %)'

      ! Start-time sensitivity to the same upper end (report Table 4).
      call linear_end(a, a_max*r0, iend, bnd)
      if (t_max > UNSET/2.0_dp) iend = count(t <= t_max)
      if (n_starts_t > 0) then
         write(*,'(a)') ''
         write(*,'(a)') '  start-time sensitivity (upper end fixed):'
         write(*,'(a12,a14,a10,a12)') 't_start', 'rate', 'samples', '1 - r^2'
         rmin = huge(1.0_dp); rmax_ = -huge(1.0_dp); nk = 0
         do i = 1, n_starts_t
            call fit_log_linear(t, a, starts_t(i), t(iend), 3, sk, b, r2k, e, m)
            if (e /= 0) then
               write(*,'(f12.4,a14)') starts_t(i), 'n/a'; cycle
            end if
            write(*,'(f12.4,f14.6,i10,es12.3)') starts_t(i), sk, m, 1.0_dp-r2k
            rmin = min(rmin, sk); rmax_ = max(rmax_, sk); nk = nk + 1
         end do
         if (nk > 1) write(*,'(a,es10.3,a,f7.3,a)') '  spread ', rmax_-rmin, &
            '  (', 100.0_dp*(rmax_-rmin)/abs(win%rate), ' % of the rate)'
      end if

      if (len_trim(csv_path) > 0) then
         open(newunit=u, file=trim(csv_path), status='replace', action='write')
         write(u,'(a)') 'record,kr0,t_lo,t_hi,n,rate,one_minus_r2,'// &
            'efolds,settled,bounded,visc_exact,inviscid,visc_two_fluid'
         write(u,'(a,a,f8.4,a,es16.8,a,es16.8,a,i0,a,es16.8,a,es12.4,'// &
            'a,f8.4,a,l1,a,l1,a,es16.8,a,es16.8,a,es16.8)') &
            trim(record_path), ',', kr0_list(1), ',', win%t_lo, ',', win%t_hi, &
            ',', win%n, ',', win%rate, ',', 1.0_dp-win%r2, ',', win%efolds, ',', &
            win%settled, ',', win%bounded, ',', visx, ',', inv, ',', twov
         close(u)
         write(*,'(a,a)') '  wrote ', trim(csv_path)
      end if
   end subroutine do_record

   !> Theoretical growth rates at each --kr0: inviscid (Rayleigh), Weber's estimate, Chandrasekhar's exact viscous relation, two-fluid inviscid, and the two-fluid viscous relation of Tomotika (1935), with the departure of the two-fluid rates from Chandrasekhar's. Uses --r0 --sigma --rho-l --mu-l --rho-g --mu-g; writes the same table to --csv. With --r-wall R (added 2026-10-22), a second table gives the two-fluid viscous rate with a no-slip wall at r = R (growth_rate_two_fluid_viscous_confined, the problem the solver integrates) and its departure from the unconfined rate, and the CSV gains a column two_fluid_viscous_confined; without --r-wall the output is unchanged.
   subroutine do_theory()
      type(table_t) :: t, tw
      real(dp) :: k, inv, web, chx, two, twov, twc
      real(dp) :: twv_all(64), twc_all(64)
      integer :: i, u
      call t%init([cel('k r0 '), cel('inviscid '), cel('Weber '), &
         cel('Chandrasekhar '), cel('two-fluid inv. '), cel('Tomotika '), &
         cel('Tomotika vs Ch. ')], 'rrrrrrr')
      if (len_trim(csv_path) > 0) then
         open(newunit=u, file=trim(csv_path), status='replace', action='write')
         if (r_wall > 0.0_dp) then
            write(u,'(a)') 'kr0,inviscid,weber,chandrasekhar,two_fluid_inviscid,'// &
               'two_fluid_viscous,two_fluid_viscous_confined'
         else
            write(u,'(a)') 'kr0,inviscid,weber,chandrasekhar,two_fluid_inviscid,'// &
               'two_fluid_viscous'
         end if
      end if
      do i = 1, n_kr0
         k = kr0_list(i) / r0
         inv = growth_rate(k, r0, sigma, rho_l)
         web = growth_rate_viscous(k, r0, sigma, rho_l, mu_l)
         chx = growth_rate_viscous_exact(k, r0, sigma, rho_l, mu_l)
         two = growth_rate_two_fluid_inviscid(k, r0, sigma, rho_l, rho_g)
         twov = growth_rate_two_fluid_viscous(k, r0, sigma, rho_l, mu_l, rho_g, &
                                              mu_g)
         twc = 0.0_dp
         if (r_wall > 0.0_dp) twc = growth_rate_two_fluid_viscous_confined(k, r0, &
            sigma, rho_l, mu_l, rho_g, mu_g, r_wall)
         twv_all(i) = twov
         twc_all(i) = twc
         if (chx > 0.0_dp .and. .not. is_nan(twov)) then
            call t%row([cel(fmt_f(kr0_list(i), 3)), cel(fmt_f(inv, 6)), &
               cel(fmt_f(web, 6)), cel(rate_text(chx)), cel(fmt_f(two, 6)), &
               cel(rate_text(twov)), cel(fmt_pct((twov - chx) / chx, 3))])
         else
            call t%row([cel(fmt_f(kr0_list(i), 3)), cel(fmt_f(inv, 6)), &
               cel(fmt_f(web, 6)), cel(rate_text(chx)), cel(fmt_f(two, 6)), &
               cel(rate_text(twov)), cel('n/a')])
         end if
         if (len_trim(csv_path) > 0) then
            if (r_wall > 0.0_dp) then
               write(u,'(f8.4,6(a,es20.12))') kr0_list(i), ',', inv, ',', web, ',', &
                  chx, ',', two, ',', twov, ',', twc
            else
               write(u,'(f8.4,5(a,es20.12))') kr0_list(i), ',', inv, ',', web, ',', &
                  chx, ',', two, ',', twov
            end if
         end if
      end do
      call t%emit('Growth rates: Oh = '//trim(fmt_f(mu_l / sqrt(rho_l * sigma * r0), &
         4))//', rho_g/rho_l = '//trim(fmt_e(rho_g / rho_l))//', mu_g/mu_l = '// &
         trim(fmt_e(merge(mu_g / mu_l, 0.0_dp, mu_l > 0.0_dp))))
      if (r_wall > 0.0_dp) then
         call tw%init([cel('k r0 '), cel('Tomotika '), cel('with wall '), &
            cel('wall effect ')], 'rrrr')
         do i = 1, n_kr0
            if (twv_all(i) > 0.0_dp .and. .not. is_nan(twc_all(i))) then
               call tw%row([cel(fmt_f(kr0_list(i), 3)), cel(rate_text(twv_all(i))), &
                  cel(rate_text(twc_all(i))), &
                  cel(fmt_pct((twc_all(i) - twv_all(i)) / twv_all(i), 3))])
            else
               call tw%row([cel(fmt_f(kr0_list(i), 3)), cel(rate_text(twv_all(i))), &
                  cel(rate_text(twc_all(i))), cel('n/a')])
            end if
         end do
         call tw%emit('Two-fluid viscous rate with a no-slip wall at r = '// &
            trim(fmt_f(r_wall, 3)))
      end if
      if (len_trim(csv_path) > 0) then
         close(u)
         write(*,'(a,a)') 'wrote ', trim(csv_path)
      end if
   end subroutine do_theory

   !> Drift of an unperturbed run (the twin) from its initial state, computed by nsmfp_mod::base_drift: for every plotfile of --base, the largest and the radius-weighted root-mean-square change of each --fields component since the first plotfile. For a column at rest this measures the discretisation drift (parasitic currents, level-set redistancing) that the twin subtraction removes. Writes every snapshot to --csv; prints a summary.
   subroutine do_drift()
      real(dp), allocatable :: tt(:), dmax(:,:), drms(:,:), q0max(:)
      integer :: e, k, c, u, n, nc
      integer :: grid(3)

      if (len_trim(base_dir) == 0) then
         write(*,'(a)') 'drift needs --base DIR'; stop 1
      end if
      call base_drift(trim(base_dir), fields(1:nfields), tt, dmax, drms, q0max, &
                      grid, e, trim(plt_prefix))
      if (e /= 0) then
         write(*,'(a,i0,a,a)') 'ERROR in base_drift (ierr=', e, ') for ', &
            trim(base_dir)
         stop 2
      end if
      n = size(tt); nc = size(dmax, 1)
      write(*,'(a,a)') 'drift of ', trim(base_dir)
      write(*,'(a,i0,a,f10.4,a,f10.4,a,i0,a,i0)') '  plotfiles ', n, '   t = ', &
         tt(1), ' .. ', tt(n), '   grid ', grid(1), ' x ', grid(2)
      do c = 1, nc
         write(*,'(a,a20,a,es11.3,a,es11.3,a,es11.3)') '  ', trim(fields(c)), &
            '  max |dQ| over run ', maxval(dmax(c,:)), '   final rms ', &
            drms(c,n), '   max |Q0| ', q0max(c)
      end do
      if (len_trim(csv_path) > 0) then
         open(newunit=u, file=trim(csv_path), status='replace', action='write')
         write(u,'(a)', advance='no') 'time'
         do c = 1, nc
            write(u,'(a)', advance='no') ',max_'//trim(fields(c))//',rms_'// &
               trim(fields(c))
         end do
         write(u,'(a)') ''
         do k = 1, n
            write(u,'(es20.12)', advance='no') tt(k)
            do c = 1, nc
               write(u,'(2(a,es14.6))', advance='no') ',', dmax(c,k), ',', &
                  drms(c,k)
            end do
            write(u,'(a)') ''
         end do
         close(u)
         write(*,'(a,a)') '  wrote ', trim(csv_path)
      end if
   end subroutine do_drift

   subroutine usage()
      write(*,'(a)') 'usage: run_lsa <analyse|linearity|dispersion|converge|record|export|'// &
         'theory|drift> [options]'
      write(*,'(a)') '  --perturbed DIR   perturbed run directory'
      write(*,'(a)') '  --base DIR        unperturbed twin directory'
      write(*,'(a)') '  --kr0 V[,V...]    wavenumber(s) k*r0'
      write(*,'(a)') '  --eps V           perturbation amplitude'
      write(*,'(a)') '  --case EPS=DIR    amplitude case (repeatable)'
      write(*,'(a)') '  --root DIR        root for dispersion sweep'
      write(*,'(a)') '  --fields A,B      plotfile components (default velocities)'
      write(*,'(a)') '  --prefix NAME     plotfile prefix (default nddataPLT)'
      write(*,'(a)') '  --skip N --stride N --max-snapshots N'
      write(*,'(a)') '  --rank N --tol V --t-min V --t-max V'
      write(*,'(a)') '  --r0 V --sigma V --rho-l V --mu-l V --rho-g V --mu-g V'
      write(*,'(a)') '  --starts N[,N...] snapshot offsets for converge'
      write(*,'(a)') '  --strides N[,N...] sub-sampling factors for converge'
      write(*,'(a)') '  --ranks N[,N...]  truncation ranks for converge'
      write(*,'(a)') '  --weight none|volume  inner product of the snapshots'
      write(*,'(a)') '  --archive DIR     save singular values, eigenvalues, modes'
      write(*,'(a)') '  --n-keep N        modes kept in the archive (default 10)'
      write(*,'(a)') '  --select dominant|max  DMD mode selection (default dominant)'
      write(*,'(a)') '  --rank-scan LO:HI rate at every rank LO..HI, median and IQR'
      write(*,'(a)') '  --band H          keep rows with |r - r0| < H dr (narrow band)'
      write(*,'(a)') '  --krylov DIR      read the members nddataPLT<step>LSA<k> of the solver''s'
      write(*,'(a)') '                    built-in LSA driver (--skip/--max-snapshots count members)'
      write(*,'(a)') '  --export DIR      save the matched perturbation sequence'
      write(*,'(a)') '  --from-export DIR analyse a saved sequence, not plotfiles'
      write(*,'(a)') '  --csv PATH        machine-readable table (converge, record)'
      write(*,'(a)') ' record: --record FILE [--a-max 0.1] [--width 0.8]'
      write(*,'(a)') '         [--rel-tol 0.005] [--starts-t T,T,...] [--t-min] [--t-max]'
      write(*,'(a)') ' theory: --kr0 V,V,... [--mu-l --rho-g --mu-g] [--r-wall R] [--csv PATH]'
      write(*,'(a)') '         --r-wall R adds the two-fluid viscous rate with a no-slip wall at r = R'
      write(*,'(a)') ' drift:  --base DIR [--fields A,B] [--csv PATH]'
      write(*,'(a)') '  --save PATH       write results as JSON'
   end subroutine usage

   subroutine parse_args(ierr)
      integer, intent(out) :: ierr
      integer :: na, i, p
      character(len=512) :: a, v

      ierr = 0
      na = command_argument_count()
      if (na < 1) then
         ierr = 1; return
      end if
      call get_command_argument(1, cmd)

      i = 2
      do while (i <= na)
         call get_command_argument(i, a)
         select case (trim(a))
         case ('--perturbed'); call nextv(i,na,v,ierr); pert_dir = v
         case ('--base');      call nextv(i,na,v,ierr); base_dir = v
         case ('--root');      call nextv(i,na,v,ierr); root = v
         case ('--save');      call nextv(i,na,v,ierr); out_json = v
         case ('--kr0')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            call parse_real_list(v, kr0_list, n_kr0, ierr)
         case ('--eps')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            call parse_real_list(v, eps_list, n_eps, ierr)
         case ('--case')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            p = index(v, '=')
            if (p == 0) then
               write(*,'(a)') '--case expects EPS=DIR'; ierr = 2; return
            end if
            if (n_cases >= size(case_eps)) then
               write(*,'(a,i0,a)') 'at most ', size(case_eps), &
                  ' --case options are accepted'
               ierr = 2; return
            end if
            n_cases = n_cases + 1
            call read_real(v(1:p-1), case_eps(n_cases), ierr)
            case_dirs(n_cases) = v(p+1:)
         case ('--prefix')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            if (len_trim(v) > len(plt_prefix)) then
               write(*,'(a)') '--prefix is too long'; ierr = 2; return
            end if
            plt_prefix = trim(v)
         case ('--fields')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            call parse_name_list(v, fields, nfields, ierr)
         case ('--starts')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            call parse_int_list(v, starts, n_starts, ierr)
         case ('--skip');   call nextv(i,na,v,ierr); call ri(v, skip, ierr)
         case ('--stride'); call nextv(i,na,v,ierr); call ri(v, stride, ierr)
         case ('--max-snapshots'); call nextv(i,na,v,ierr); call ri(v, maxsnap, ierr)
         case ('--rank');   call nextv(i,na,v,ierr); call ri(v, rank, ierr)
         case ('--tol');    call nextv(i,na,v,ierr); call rr(v, tol, ierr)
         case ('--t-min');  call nextv(i,na,v,ierr); call rr(v, t_min, ierr)
         case ('--t-max');  call nextv(i,na,v,ierr); call rr(v, t_max, ierr)
         case ('--r0');     call nextv(i,na,v,ierr); call rr(v, r0, ierr)
         case ('--sigma');  call nextv(i,na,v,ierr); call rr(v, sigma, ierr)
         case ('--rho-l');  call nextv(i,na,v,ierr); call rr(v, rho_l, ierr)
         case ('--mu-l');   call nextv(i,na,v,ierr); call rr(v, mu_l, ierr)
         case ('--rho-g');  call nextv(i,na,v,ierr); call rr(v, rho_g, ierr)
         case ('--mu-g');   call nextv(i,na,v,ierr); call rr(v, mu_g, ierr)
         case ('--r-wall'); call nextv(i,na,v,ierr); call rr(v, r_wall, ierr)
         case ('--strides')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            call parse_int_list(v, strides, n_strides, ierr)
         case ('--ranks')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            call parse_int_list(v, ranks, n_ranks, ierr)
         case ('--weight')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            weight = trim(v)
            if (weight /= 'none' .and. weight /= 'volume') then
               write(*,'(a)') '--weight must be none or volume'; ierr = 2
            end if
         case ('--archive');  call nextv(i,na,v,ierr); archive_dir = v
         case ('--export');   call nextv(i,na,v,ierr); export_dir = v
         case ('--from-export'); call nextv(i,na,v,ierr); from_export = v
         case ('--csv');      call nextv(i,na,v,ierr); csv_path = v
         case ('--record');   call nextv(i,na,v,ierr); record_path = v
         case ('--n-keep');   call nextv(i,na,v,ierr); call ri(v, n_keep, ierr)
         case ('--select')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            select_mode = trim(v)
            if (select_mode /= 'dominant' .and. select_mode /= 'max') then
               write(*,'(a)') '--select must be dominant or max'; ierr = 2
            end if
         case ('--krylov');   call nextv(i,na,v,ierr); krylov_dir = v
         case ('--band');     call nextv(i,na,v,ierr); call rr(v, band, ierr)
         case ('--rank-scan')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            p = index(v, ':')
            if (p == 0) then
               write(*,'(a)') '--rank-scan expects LO:HI'; ierr = 2; return
            end if
            call read_int(v(1:p-1), rscan_lo, ierr); if (ierr /= 0) return
            call read_int(v(p+1:), rscan_hi, ierr)
         case ('--a-max');    call nextv(i,na,v,ierr); call rr(v, a_max, ierr)
         case ('--width');    call nextv(i,na,v,ierr); call rr(v, width, ierr)
         case ('--rel-tol');  call nextv(i,na,v,ierr); call rr(v, rel_tol, ierr)
         case ('--starts-t')
            call nextv(i,na,v,ierr); if (ierr /= 0) return
            call parse_real_list(v, starts_t, n_starts_t, ierr)
         case default
            write(*,'(a,a)') 'unknown option: ', trim(a)
            ierr = 3; return
         end select
         if (ierr /= 0) return
         i = i + 1
      end do
      if (skip < 0 .or. stride < 1 .or. maxsnap < 0) then
         write(*,'(a)') '--skip and --max-snapshots must be >= 0, --stride >= 1'
         ierr = 5; return
      end if
      if (r0 <= 0.0_dp .or. sigma <= 0.0_dp .or. rho_l <= 0.0_dp .or. &
          mu_l < 0.0_dp .or. rho_g < 0.0_dp .or. mu_g < 0.0_dp) then
         write(*,'(a)') '--r0, --sigma, --rho-l must be positive, '// &
            '--mu-l, --rho-g and --mu-g non-negative'
         ierr = 5; return
      end if
      if (r_wall /= 0.0_dp .and. r_wall <= r0) then
         write(*,'(a)') '--r-wall must be larger than --r0 (the wall lies outside the column)'
         ierr = 5; return
      end if

   end subroutine parse_args

   ! Read the value of the option just consumed, unless nextv already failed (ierr non-zero), in which case nothing is read and the error stands.
   subroutine ri(text, k, ierr)
      character(len=*), intent(in) :: text
      integer, intent(inout) :: k
      integer, intent(inout) :: ierr
      if (ierr /= 0) return
      call read_int(text, k, ierr)
   end subroutine ri

   subroutine rr(text, x, ierr)
      character(len=*), intent(in) :: text
      real(dp), intent(inout) :: x
      integer, intent(inout) :: ierr
      if (ierr /= 0) return
      call read_real(text, x, ierr)
   end subroutine rr

   !> Parse "a,b,c" into integers; a non-integer such as 1.6 is an error rather than being rounded.
   subroutine parse_int_list(s, v, n, ierr)
      character(len=*), intent(in) :: s
      integer, intent(out) :: v(:)
      integer, intent(out) :: n, ierr
      character(len=32) :: tok(size(v))
      integer :: k
      call parse_name_list(s, tok, n, ierr)
      if (ierr /= 0) return
      do k = 1, n
         call read_int(tok(k), v(k), ierr)
         if (ierr /= 0) return
      end do
   end subroutine parse_int_list

   subroutine nextv(i, na, v, ierr)
      integer, intent(inout) :: i
      integer, intent(in) :: na
      character(len=*), intent(out) :: v
      integer, intent(out) :: ierr
      ierr = 0
      i = i + 1
      if (i > na) then
         ierr = 4; return
      end if
      call get_command_argument(i, v)
   end subroutine nextv

end program run_lsa
