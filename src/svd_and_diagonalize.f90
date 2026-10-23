!> svd_and_diagonalize: companion-matrix eigenanalysis of a Krylov (or evenly sampled snapshot) sequence, as in Prof. Sussman's unpublished draft description of the restarting method (Section "Equivalency of the above analysis to the least squares approach and a new generalization").
!>
!> Steps:
!>   1. Read the sequence v_1 .. v_N from a lossless export written by `run_lsa export` (perturbation.bin, perturbation.txt, times.txt); for the solver's restarting Krylov driver use `run_lsa export --krylov DIR`.
!>   2. Solve the article's eq. (67), (B^T W B) a = B^T W v_N with B = [v_1 .. v_(N-1)], for a by the SVD (pseudo-inverse) of B^T W B. The same problem is also solved directly from sqrt(W) B (dgelss) as a check, and the two solutions are compared.
!>   3. Assemble the companion matrix S (ones on the subdiagonal, a in the last column).
!>   4. Diagonalise S = X Gamma X^(-1); eigenvalues sorted by decreasing |mu|. Continuous growth rates s = ln(mu)/tau, with tau the sampling interval of the export. Ritz amplitudes |b_j| ||P_j||_w with P = B X and X b = e_1.
!>
!> Usage:
!>   ./bin/svd_and_diagonalize --export DIR [--skip K] [--max-members N] [--comp-weights w1,w2,...] [--volume] [--band H] [--r0 R] [--rcond EPS] [--kr0 X] [--mu MU] [--nshow N] [--out DIR]
!>
!> Options:
!>   --export DIR        export directory (required)
!>   --skip K            drop the first K columns (start-up transient); default 0
!>   --max-members N     use at most N columns after the skip (N = number of Krylov vectors, so B has N-1 columns); default all
!>   --comp-weights LIST weight of each field block in the export, in its field order (e.g. 10,1,1 makes the level set dominate the velocities, Remark 2 of the article); default 1 for every block
!>   --volume            multiply every weight by the cell-centre radius r (axisymmetric volume element; note w_i weights a squared residual, so this is r, not sqrt(r))
!>   --band H            set the weight to zero outside |r - r0| < H dr (restrict to a band around the undisturbed interface)
!>   --r0 R              undisturbed radius for --band and the theory; default 1
!>   --rcond EPS         relative singular-value cutoff for the pseudo-inverse of B^T W B; default 1e-12
!>   --kr0 X             wavenumber k r0, to compare growth rates with Chandrasekhar's exact viscous rate
!>   --mu MU             liquid viscosity for that comparison (sigma = rho = 1); default 0.02
!>   --nshow N           number of eigenvalues printed; default 12
!>   --out DIR           write a.txt, singular_values.txt, eigenvalues.txt and summary.txt there
!>
!> Exit status: 0 success, 1 bad command line, 2 export unreadable or too short, 3 linear algebra failure.
program svd_and_diagonalize
   use companion_mod
   use archive_mod, only: read_perturbation_export, make_dir
   use dispersion_mod, only: growth_rate_viscous_exact
   use util_mod, only: read_real, read_int, parse_real_list, cell_centres
   implicit none

   character(len=512) :: export_dir, out_dir
   real(dp) :: rcond, r0, band, kr0, mu_l
   real(dp) :: cw(16)
   integer :: ncw, skip, maxm, nshow
   logical :: volume
   real(dp), allocatable :: X(:,:), times(:), w(:), B(:,:), vN(:), a(:), &
                            a_chk(:), sv(:), S(:,:), amp(:), pnorm(:), rc(:)
   complex(dp), allocatable :: mu(:), Xev(:,:)
   real(dp) :: dt, tau, dom(4), visc, res, diff, dr
   integer :: grid(3), ierr, nsp, ncol, first, last, nvec, n, rank, rank_chk
   integer :: nx, ny, nc, c, j, i, k, row
   character(len=64) :: wt

   call defaults()
   call parse_args(ierr)
   if (ierr /= 0) then
      call usage(); stop 1
   end if

   ! ---- 1. the sequence -------------------------------------------------
   call read_perturbation_export(trim(export_dir), X, times, dt, grid, ierr, &
                                 dom, wt)
   if (ierr /= 0) then
      write(*,'(a,i0,a)') 'ERROR reading the export (ierr=', ierr, ') in '// &
         trim(export_dir)
      stop 2
   end if
   nsp = size(X, 1); ncol = size(X, 2)
   first = skip + 1
   last = ncol
   if (maxm > 0) last = min(ncol, skip + maxm)
   nvec = last - first + 1
   if (nvec < 3) then
      write(*,'(a,i0,a)') 'ERROR: only ', nvec, ' vectors after --skip/'// &
         '--max-members; at least 3 are needed'
      stop 2
   end if
   n = nvec - 1
   B = X(:, first:last - 1)
   vN = X(:, last)
   ! The export's columns are equally spaced in time, so tau is the spacing between consecutive vectors.
   tau = dt
   if (size(times) >= 2) tau = times(2) - times(1)

   ! ---- weights ------------------------------------------------------------
   nx = grid(1); ny = grid(2); nc = grid(3)
   allocate(w(nsp))
   w = 1.0_dp
   if (nx * ny * nc == nsp .and. nc >= 1) then
      allocate(rc(nx))
      call cell_centres(dom(1), dom(3), nx, rc)
      dr = (dom(3) - dom(1)) / real(nx, dp)
      do c = 1, nc
         do j = 1, ny
            do i = 1, nx
               row = (c - 1) * nx * ny + (j - 1) * nx + i
               if (c <= ncw) w(row) = cw(c)
               if (volume) w(row) = w(row) * rc(i)
               if (band > 0.0_dp) then
                  if (.not. (abs(rc(i) - r0) < band * dr)) w(row) = 0.0_dp
               end if
            end do
         end do
      end do
   else if (ncw > 0 .or. volume .or. band > 0.0_dp) then
      write(*,'(a)') 'WARNING: the export does not record its grid; '// &
         '--comp-weights, --volume and --band are ignored'
   end if
   if (trim(wt) == 'volume') write(*,'(a)') 'note: the export rows already'// &
      ' carry sqrt(r) weights (export written with --weight volume)'

   ! ---- 2. eq. (67) by SVD -------------------------------------------------
   allocate(a(n), a_chk(n), sv(n))
   call solve_normal_svd(B, vN, w, rcond, a, sv, rank, ierr)
   if (ierr /= 0) then
      write(*,'(a,i0)') 'ERROR in the SVD solve of B^T W B a = B^T W v_N,'// &
         ' ierr=', ierr
      stop 3
   end if
   call solve_weighted_lstsq(B, vN, w, rcond, a_chk, rank_chk, ierr)
   diff = -1.0_dp
   if (ierr == 0) diff = sqrt(sum((a - a_chk)**2)) / &
                         max(sqrt(sum(a_chk**2)), tiny(1.0_dp))
   res = weighted_residual(B, vN, w, a)

   ! ---- 3. companion matrix, 4. its eigen-decomposition -----------------
   allocate(S(n, n), mu(n), Xev(n, n), amp(n), pnorm(n))
   call assemble_companion(a, S)
   call diagonalize_companion(S, mu, Xev, ierr)
   if (ierr /= 0) then
      write(*,'(a,i0)') 'ERROR diagonalising S, ierr=', ierr
      stop 3
   end if
   call ritz_amplitudes(B, Xev, w, amp, pnorm, ierr)
   if (ierr /= 0) write(*,'(a)') 'WARNING: X is singular; amplitudes not computed'

   visc = -1.0_dp
   if (kr0 > 0.0_dp) visc = growth_rate_viscous_exact(kr0 / r0, r0, 1.0_dp, &
                                                      1.0_dp, mu_l)
   call report()
   if (len_trim(out_dir) > 0) call write_outputs()

contains

   subroutine defaults()
      export_dir = ''; out_dir = ''
      rcond = 1.0e-12_dp; r0 = 1.0_dp; band = 0.0_dp; kr0 = -1.0_dp
      mu_l = 0.02_dp; cw = 1.0_dp; ncw = 0; skip = 0; maxm = 0; nshow = 12
      volume = .false.
   end subroutine defaults

   !> Continuous rate s = ln(mu)/tau of eigenvalue mu.
   complex(dp) function srate(m)
      complex(dp), intent(in) :: m
      if (abs(m) > 0.0_dp) then
         srate = log(m) / tau
      else
         srate = cmplx(-huge(1.0_dp), 0.0_dp, dp)
      end if
   end function srate

   subroutine report()
      complex(dp) :: s
      integer :: jd
      write(*,'(a)') 'svd_and_diagonalize: companion-matrix eigenanalysis '// &
         '(draft article eq. 67)'
      write(*,'(a,a)') '  export          ', trim(export_dir)
      write(*,'(a,i0,a,i0,a,i0,a)') '  vectors         v_', first, ' .. v_', &
         last, '  (N = ', nvec, ')'
      write(*,'(a,i0,a,i0)') '  space size M    ', nsp, '   B is M x ', n
      write(*,'(a,es14.6)') '  tau             ', tau
      write(*,'(a,i0,a,i0,a)') '  weights         ', count(w > 0.0_dp), &
         ' of ', nsp, ' rows non-zero'
      write(*,'(a,i0,a,i0,a,es10.3,a)') '  rank of B^T W B ', rank, ' of ', &
         n, '  (rcond ', rcond, ')'
      write(*,'(a,es12.4,a,es12.4)') '  singular values ', sv(1), ' .. ', sv(n)
      write(*,'(a,es12.4)') '  residual ||v_N - B a||_w / ||v_N||_w ', res
      if (diff >= 0.0_dp) write(*,'(a,es12.4,a,i0,a)') &
         '  check: rel. difference from direct sqrt(W)B solve ', diff, &
         '  (rank ', rank_chk, ')'
      write(*,'(a)') ''
      write(*,'(a)') '   j       Re mu        Im mu        |mu|      Re s'// &
         '         Im s       amplitude'
      do jd = 1, min(nshow, n)
         s = srate(mu(jd))
         write(*,'(i4,3f13.7,2f13.7,es14.4)') jd, real(mu(jd), dp), &
            aimag(mu(jd)), abs(mu(jd)), real(s, dp), aimag(s), amp(jd)
      end do
      jd = dominant()
      if (jd > 0) then
         s = srate(mu(jd))
         write(*,'(a)') ''
         write(*,'(a,i0,a,f12.7)') '  dominant growing mode: j = ', jd, &
            '   s = ', real(s, dp)
         if (visc > 0.0_dp) write(*,'(a,f12.7,a,f8.3,a)') &
            '  viscous exact (Chandrasekhar) ', visc, '   (', &
            100.0_dp * (real(s, dp) / visc - 1.0_dp), ' %)'
      end if
   end subroutine report

   !> The growing, near-stationary mode (|Im s| <= 0.1 Re s) of largest amplitude, as dmd_mod::dominant_growth_rate; 0 if none.
   integer function dominant()
      integer :: jd
      real(dp) :: best
      complex(dp) :: s
      dominant = 0; best = -1.0_dp
      do jd = 1, n
         s = srate(mu(jd))
         if (.not. (real(s, dp) > 0.0_dp)) cycle
         if (abs(aimag(s)) > 0.1_dp * real(s, dp)) cycle
         if (amp(jd) > best) then
            best = amp(jd); dominant = jd
         end if
      end do
   end function dominant

   subroutine write_outputs()
      integer :: u, jd
      complex(dp) :: s
      call make_dir(trim(out_dir))
      open(newunit=u, file=trim(out_dir)//'/a.txt', status='replace')
      write(u,'(a)') '# k  a_k   (v_N = sum_k a_k v_k + r, eq. 67)'
      do jd = 1, n
         write(u,'(i5,es24.16)') jd, a(jd)
      end do
      close(u)
      open(newunit=u, file=trim(out_dir)//'/singular_values.txt', &
           status='replace')
      write(u,'(a)') '# i  sigma_i(B^T W B)  sigma_i/sigma_1'
      do jd = 1, n
         write(u,'(i5,2es24.16)') jd, sv(jd), sv(jd) / sv(1)
      end do
      close(u)
      open(newunit=u, file=trim(out_dir)//'/eigenvalues.txt', &
           status='replace')
      write(u,'(a)') '# j  Re_mu  Im_mu  |mu|  Re_s  Im_s  amplitude'// &
         '  (sorted by |mu|; s = ln(mu)/tau)'
      do jd = 1, n
         s = srate(mu(jd))
         write(u,'(i5,6es24.16)') jd, real(mu(jd), dp), aimag(mu(jd)), &
            abs(mu(jd)), real(s, dp), aimag(s), amp(jd)
      end do
      close(u)
      open(newunit=u, file=trim(out_dir)//'/summary.txt', status='replace')
      write(u,'(a,a)') 'export = ', trim(export_dir)
      write(u,'(a,i0)') 'first_vector = ', first
      write(u,'(a,i0)') 'n_vectors = ', nvec
      write(u,'(a,es24.16)') 'tau = ', tau
      write(u,'(a,i0)') 'rank = ', rank
      write(u,'(a,es24.16)') 'rcond = ', rcond
      write(u,'(a,es24.16)') 'residual = ', res
      write(u,'(a,es24.16)') 'check_rel_difference = ', diff
      jd = dominant()
      if (jd > 0) then
         write(u,'(a,es24.16)') 'dominant_growth_rate = ', real(srate(mu(jd)), dp)
      else
         write(u,'(a)') 'dominant_growth_rate = nan'
      end if
      if (visc > 0.0_dp) write(u,'(a,es24.16)') 'visc_exact = ', visc
      close(u)
      write(*,'(a,a)') '  written to ', trim(out_dir)
   end subroutine write_outputs

   subroutine usage()
      write(*,'(a)') 'usage: svd_and_diagonalize --export DIR [--skip K] '// &
         '[--max-members N] [--comp-weights w1,w2,..] [--volume] [--band H]'
      write(*,'(a)') '       [--r0 R] [--rcond EPS] [--kr0 X] [--mu MU] '// &
         '[--nshow N] [--out DIR]'
   end subroutine usage

   subroutine parse_args(ierr)
      integer, intent(out) :: ierr
      integer :: na, ia
      character(len=512) :: opt, v
      ierr = 0
      na = command_argument_count()
      ia = 1
      do while (ia <= na)
         call get_command_argument(ia, opt)
         select case (trim(opt))
         case ('--volume')
            volume = .true.
         case ('--export', '--skip', '--max-members', '--comp-weights', &
               '--band', '--r0', '--rcond', '--kr0', '--mu', '--nshow', '--out')
            if (ia == na) then
               write(*,'(a,a)') 'error: missing value for ', trim(opt)
               ierr = 2; return
            end if
            ia = ia + 1
            call get_command_argument(ia, v)
            select case (trim(opt))
            case ('--export');       export_dir = v
            case ('--out');          out_dir = v
            case ('--skip');         call read_int(v, skip, ierr)
            case ('--max-members');  call read_int(v, maxm, ierr)
            case ('--nshow');        call read_int(v, nshow, ierr)
            case ('--band');         call read_real(v, band, ierr)
            case ('--r0');           call read_real(v, r0, ierr)
            case ('--rcond');        call read_real(v, rcond, ierr)
            case ('--kr0');          call read_real(v, kr0, ierr)
            case ('--mu');           call read_real(v, mu_l, ierr)
            case ('--comp-weights'); call parse_real_list(v, cw, ncw, ierr)
            end select
            if (ierr /= 0) return
         case default
            write(*,'(a,a)') 'error: unknown option ', trim(opt)
            ierr = 3; return
         end select
         ia = ia + 1
      end do
      if (len_trim(export_dir) == 0) then
         write(*,'(a)') 'error: --export DIR is required'
         ierr = 4
      else if (skip < 0 .or. maxm < 0 .or. nshow < 1 .or. rcond < 0.0_dp &
               .or. r0 <= 0.0_dp .or. band < 0.0_dp) then
         write(*,'(a)') 'error: option out of range'
         ierr = 5
      end if
      if (ierr == 0) then
         do k = 1, ncw
            if (cw(k) < 0.0_dp) then
               write(*,'(a)') 'error: --comp-weights must be non-negative'
               ierr = 5
            end if
         end do
      end if
   end subroutine parse_args

end program svd_and_diagonalize
