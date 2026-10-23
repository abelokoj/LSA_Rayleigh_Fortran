!> Small shared utilities: the log-linear growth fit, sorting, cell-centre axes, compact number formatting, NaN handling, and safe parsing of command-line values.
!>
!> These routines previously existed in two or three slightly different copies (the log-linear fit in nsmfp_mod and interface_mode_mod, insertion sorts in three modules, nice() in two programs). Copies that differ in a minimum sample count or a window rule give different answers for the same record, which is the opposite of what a reference implementation should do, so each now exists once, here. This module depends on nothing else in the project, so every other module may use it.
module util_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                            ieee_is_nan
   implicit none
   private

   public :: dp
   public :: fit_log_linear, window_active
   public :: argsort_real_asc, argsort_real_desc, argsort_int_asc
   public :: cell_centres, nice, nan_value, is_nan
   public :: read_real, read_int, parse_real_list, parse_name_list

   !> Sentinel meaning "no bound" for the optional time window of a fit. Callers pass -huge for an unset t_min or t_max.
   real(dp), parameter, public :: UNSET = -huge(1.0_dp)

contains

   !> True when a window bound has been set, i.e. is not the UNSET sentinel.
   pure logical function window_active(t)
      real(dp), intent(in) :: t
      window_active = t > UNSET / 2.0_dp
   end function window_active

   !> Least-squares fit of log(y) = slope * t + intercept over the samples with y > 0 inside [t_min, t_max].
   !>
   !> Either bound may be UNSET. A set t_max that does not exceed a set t_min is an error (ierr = 21), not a silent removal of the upper bound as before: a window the caller asked for and did not get is worse than a failed fit. Non-positive and NaN values of y are skipped, since their logarithm does not exist. The fit is done about the mean of t, which keeps it well conditioned when the window lies far from t = 0. ierr = 20 when fewer than min_pts samples qualify, 22 when the selected times are all equal.
   subroutine fit_log_linear(t, y, t_min, t_max, min_pts, slope, intercept, &
                             r2, ierr, npts)
      real(dp), intent(in) :: t(:), y(:), t_min, t_max
      integer, intent(in) :: min_pts
      real(dp), intent(out) :: slope, intercept, r2
      integer, intent(out) :: ierr
      integer, intent(out), optional :: npts

      integer :: i, m
      real(dp) :: tbar, ybar, stt, sty, syy, ssres, dt, dy

      slope = 0.0_dp; intercept = 0.0_dp; r2 = 0.0_dp; ierr = 0
      if (present(npts)) npts = 0
      if (window_active(t_min) .and. window_active(t_max)) then
         if (t_max <= t_min) then
            ierr = 21; return
         end if
      end if

      m = 0; tbar = 0.0_dp; ybar = 0.0_dp
      do i = 1, size(t)
         if (.not. use_sample(i)) cycle
         m = m + 1
         tbar = tbar + t(i)
         ybar = ybar + log(y(i))
      end do
      if (present(npts)) npts = m
      if (m < max(min_pts, 2)) then
         ierr = 20; return
      end if
      tbar = tbar / real(m, dp)
      ybar = ybar / real(m, dp)

      stt = 0.0_dp; sty = 0.0_dp; syy = 0.0_dp
      do i = 1, size(t)
         if (.not. use_sample(i)) cycle
         dt = t(i) - tbar
         dy = log(y(i)) - ybar
         stt = stt + dt * dt
         sty = sty + dt * dy
         syy = syy + dy * dy
      end do
      if (stt <= 0.0_dp) then
         ierr = 22; return
      end if

      slope = sty / stt
      intercept = ybar - slope * tbar
      ssres = max(syy - slope * sty, 0.0_dp)
      if (syy > 0.0_dp) then
         r2 = 1.0_dp - ssres / syy
      else
         r2 = 0.0_dp
      end if

   contains

      logical function use_sample(i)
         integer, intent(in) :: i
         use_sample = .false.
         if (.not. (y(i) > 0.0_dp)) return          ! also rejects NaN
         if (window_active(t_min) .and. t(i) < t_min) return
         if (window_active(t_max) .and. t(i) > t_max) return
         use_sample = .true.
      end function use_sample

   end subroutine fit_log_linear

   !> Indices that sort v into ascending order (stable insertion sort; the arrays here hold at most a few thousand entries).
   pure subroutine argsort_real_asc(v, idx)
      real(dp), intent(in) :: v(:)
      integer, intent(out) :: idx(:)
      integer :: i, j, key
      do i = 1, size(v)
         idx(i) = i
      end do
      do i = 2, size(v)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (v(idx(j)) <= v(key)) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = key
      end do
   end subroutine argsort_real_asc

   !> Indices that sort v into descending order (stable insertion sort).
   pure subroutine argsort_real_desc(v, idx)
      real(dp), intent(in) :: v(:)
      integer, intent(out) :: idx(:)
      integer :: i, j, key
      do i = 1, size(v)
         idx(i) = i
      end do
      do i = 2, size(v)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (v(idx(j)) >= v(key)) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = key
      end do
   end subroutine argsort_real_desc

   !> Indices that sort integer v into ascending order (stable insertion sort).
   pure subroutine argsort_int_asc(v, idx)
      integer, intent(in) :: v(:)
      integer, intent(out) :: idx(:)
      integer :: i, j, key
      do i = 1, size(v)
         idx(i) = i
      end do
      do i = 2, size(v)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (v(idx(j)) <= v(key)) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = key
      end do
   end subroutine argsort_int_asc

   !> Cell-centre coordinates of n uniform cells spanning [lo, hi].
   !>
   !> The lower edge is used explicitly. The axes were previously built as (i - 1/2) hi / n, which is right only when prob_lo = 0; the decks generated here do set prob_lo = 0, but a deck that does not would have given every radius and every axial position an offset without any warning.
   pure subroutine cell_centres(lo, hi, n, x)
      real(dp), intent(in) :: lo, hi
      integer, intent(in) :: n
      real(dp), intent(out) :: x(:)
      integer :: i
      do i = 1, n
         x(i) = lo + (real(i, dp) - 0.5_dp) * (hi - lo) / real(n, dp)
      end do
   end subroutine cell_centres

   !> Compact decimal text for filenames: 0.7 -> "0.7", 1e-3 -> "0.001". Avoids the exponent notation that list-directed output produces, so that names stay easy to type and match those written by scripts/run_all.sh, which formats with printf '%.10f' and strips trailing zeros in the same way.
   function nice(x) result(s)
      real(dp), intent(in) :: x
      character(len=64) :: s
      character(len=64) :: t
      integer :: i, n
      if (x == 0.0_dp) then
         s = '0'; return
      end if
      write(t, '(f24.10)') x
      t = adjustl(t)
      n = len_trim(t)
      do i = n, 1, -1
         if (t(i:i) /= '0') exit
         t(i:i) = ' '
      end do
      n = len_trim(t)
      if (n > 0) then
         if (t(n:n) == '.') t(n:n) = ' '
      end if
      s = trim(t)
   end function nice

   !> A quiet NaN, used to mean "no result" where a number is expected.
   function nan_value() result(x)
      real(dp) :: x
      x = ieee_value(x, ieee_quiet_nan)
   end function nan_value

   elemental logical function is_nan(x)
      real(dp), intent(in) :: x
      is_nan = ieee_is_nan(x)
   end function is_nan

   !> True when text holds a single token: no comma, slash or embedded blank. List-directed input stops at the first separator and ignores the rest, so "0.7,0.9" would otherwise be read silently as 0.7.
   pure logical function single_token(text)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: t
      t = adjustl(text)
      single_token = len_trim(t) > 0 .and. scan(trim(t), ', /') == 0
   end function single_token

   !> Read a real from text, reporting failure through ierr instead of aborting with a run-time error that names neither the option nor the value.
   subroutine read_real(text, x, ierr)
      character(len=*), intent(in) :: text
      real(dp), intent(out) :: x
      integer, intent(out) :: ierr
      integer :: ios
      x = 0.0_dp
      ios = 1
      if (single_token(text)) read(text, *, iostat=ios) x
      ierr = merge(0, 1, ios == 0)
      if (ierr /= 0) write(*, '(a,a,a)') 'error: "', trim(text), &
         '" is not a number'
   end subroutine read_real

   subroutine read_int(text, k, ierr)
      character(len=*), intent(in) :: text
      integer, intent(out) :: k
      integer, intent(out) :: ierr
      integer :: ios
      k = 0
      ios = 1
      if (single_token(text)) read(text, *, iostat=ios) k
      ierr = merge(0, 1, ios == 0)
      if (ierr /= 0) write(*, '(a,a,a)') 'error: "', trim(text), &
         '" is not an integer'
   end subroutine read_int

   !> Parse "a,b,c" into v(1:n). ierr = 1 on a malformed number, 2 when the list holds more entries than v can store; the fixed arrays of the drivers previously overflowed without a check.
   subroutine parse_real_list(s, v, n, ierr)
      character(len=*), intent(in) :: s
      real(dp), intent(out) :: v(:)
      integer, intent(out) :: n, ierr
      integer :: p
      character(len=len(s)) :: t
      t = s; n = 0; ierr = 0
      do
         p = index(t, ',')
         n = n + 1
         if (n > size(v)) then
            write(*, '(a,i0,a)') 'error: at most ', size(v), &
               ' values may be given in one list'
            ierr = 2; return
         end if
         if (p == 0) then
            call read_real(t, v(n), ierr)
            return
         end if
         call read_real(t(1:p-1), v(n), ierr)
         if (ierr /= 0) return
         t = t(p+1:)
      end do
   end subroutine parse_real_list

   !> Parse "a,b,c" into a list of names, with the same overflow check.
   subroutine parse_name_list(s, v, n, ierr)
      character(len=*), intent(in) :: s
      character(len=*), intent(out) :: v(:)
      integer, intent(out) :: n, ierr
      integer :: p
      character(len=len(s)) :: t
      t = s; n = 0; ierr = 0
      do
         p = index(t, ',')
         n = n + 1
         if (n > size(v)) then
            write(*, '(a,i0,a)') 'error: at most ', size(v), &
               ' names may be given in one list'
            ierr = 2; return
         end if
         if (p == 0) then
            v(n) = trim(adjustl(t)); return
         end if
         v(n) = trim(adjustl(t(1:p-1)))
         t = t(p+1:)
      end do
   end subroutine parse_name_list

end module util_mod
