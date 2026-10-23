!> Growth rate from an interface-amplitude record, with an explicit admissible window.
!>
!> This module carries over the central methodological result of the September 2026 evidence report (NS_fluids_lib/analysis_report.typ in the solver fork at commit dddb827c, sections 5.3-5.4 and 6.3): a growth rate fitted to a(t) is meaningful only over an interval that begins after the startup transient has decayed and ends before the amplitude leaves the small-disturbance range. The interface is displaced while the fluid is at rest, so the record starts as a combination of a growing and a decaying mode (for the inviscid two-mode model a = a0 cosh(s t), whose logarithmic slope s tanh(s t) approaches s only from below); waiting removes that bias, but the amplitude keeps growing and eventually leaves the linear range. A fit outside that window can have r^2 = 0.99999 and still be biased, which is how the archived rates came to disagree by more than a per cent depending only on where the fit started.
!>
!> The window used here is reproducible and is printed with every result:
!>
!> - the upper end is the last sample before a/r0 exceeds a_max (default 0.1, the operational bound of the study);
!> - the lower end is the earliest centre of a sliding local-rate window from which every later local rate, up to the upper end, stays within rel_tol (run_lsa default 0.5 per cent) of the last one: the point from which the record has stopped climbing. On the archived meas32 record this gives t >= 8.64 and s = 0.32794, against 0.32804 for the report's fixed cutoff t >= 8.81; a 1 per cent tolerance starts at 7.77 and gives 0.32745, still on the rising part;
!> - the fit is the least-squares slope of ln a over that interval (util_mod::fit_log_linear).
!>
!> The start-time sensitivity table of the report (its Table 4) is produced by refitting from a list of fixed start times to the same upper end. The spread over that table, not the number of printed digits, is the honest statement of how precisely the rate is known.
module record_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use util_mod, only: fit_log_linear, UNSET, argsort_real_asc
   implicit none
   private

   public :: dp
   public :: record_window_t, read_amplitude_record, linear_end
   public :: sliding_rates, settled_start, analyse_record

   type :: record_window_t
      real(dp) :: t_lo = 0.0_dp, t_hi = 0.0_dp   !< window actually used
      real(dp) :: rate = 0.0_dp, r2 = 0.0_dp
      integer  :: n = 0                          !< samples in the window
      real(dp) :: a_lo = 0.0_dp, a_hi = 0.0_dp   !< amplitude at the window ends
      real(dp) :: efolds = 0.0_dp                !< rate * (t_hi - t_lo)
      logical  :: settled = .false.              !< a settled start was found
      logical  :: bounded = .false.              !< the amplitude bound was reached
      integer  :: ierr = 0
   end type record_window_t

contains

   !> Read time and amplitude (columns 1 and 2) from a record written by dump_amplitude, sorted by time. Comment lines are skipped; rows with a non-positive amplitude (the -1 missing-value marker) are dropped. ierr = 1 if the file cannot be opened, 2 if fewer than 3 usable rows remain.
   subroutine read_amplitude_record(path, t, a, ierr)
      character(len=*), intent(in) :: path
      real(dp), allocatable, intent(out) :: t(:), a(:)
      integer, intent(out) :: ierr
      integer :: u, ios, n, i
      character(len=1024) :: line
      real(dp) :: x(2)
      real(dp), allocatable :: tt(:), aa(:)
      integer, allocatable :: idx(:)

      ierr = 0
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 1; allocate(t(0), a(0)); return
      end if
      n = 0
      allocate(tt(1024), aa(1024))
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         read(line, *, iostat=ios) x
         if (ios /= 0) cycle
         if (.not. (x(2) > 0.0_dp)) cycle
         n = n + 1
         if (n > size(tt)) then
            tt = [tt, tt]; aa = [aa, aa]
         end if
         tt(n) = x(1); aa(n) = x(2)
      end do
      close(u)
      if (n < 3) then
         ierr = 2; allocate(t(0), a(0)); return
      end if
      allocate(idx(n))
      call argsort_real_asc(tt(1:n), idx)
      t = tt(idx); a = aa(idx)
      do i = 2, n
         if (t(i) <= t(i-1)) then
            ierr = 3; return                  ! duplicate times
         end if
      end do
   end subroutine read_amplitude_record

   !> Index of the last sample before a exceeds a_max; size(a) if it never does. bounded is true when the bound was reached.
   subroutine linear_end(a, a_max, iend, bounded)
      real(dp), intent(in) :: a(:), a_max
      integer, intent(out) :: iend
      logical, intent(out) :: bounded
      integer :: i
      iend = size(a); bounded = .false.
      do i = 1, size(a)
         if (a(i) > a_max) then
            iend = max(i - 1, 1); bounded = .true.; return
         end if
      end do
   end subroutine linear_end

   !> Local growth rate d(ln a)/dt over sliding windows of the given width, centred on every sample whose window lies inside [t(1), t(n)] and holds at least 4 samples.
   subroutine sliding_rates(t, a, n, width, centres, rates, nout)
      real(dp), intent(in) :: t(:), a(:), width
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: centres(:), rates(:)
      integer, intent(out) :: nout
      integer :: i, e, m
      real(dp) :: s, b, r2

      allocate(centres(n), rates(n))
      nout = 0
      do i = 1, n
         if (t(i) - 0.5_dp*width < t(1) .or. t(i) + 0.5_dp*width > t(n)) cycle
         call fit_log_linear(t(1:n), a(1:n), t(i) - 0.5_dp*width, &
                             t(i) + 0.5_dp*width, 4, s, b, r2, e, m)
         if (e /= 0) cycle
         nout = nout + 1
         centres(nout) = t(i)
         rates(nout) = s
      end do
   end subroutine sliding_rates

   !> Earliest window centre from which every later local rate stays within rel_tol of the last one. found is false when even the last two disagree, i.e. the record never settles.
   subroutine settled_start(centres, rates, nout, rel_tol, t_start, found)
      real(dp), intent(in) :: centres(:), rates(:), rel_tol
      integer, intent(in) :: nout
      real(dp), intent(out) :: t_start
      logical, intent(out) :: found
      integer :: i
      real(dp) :: ref

      found = .false.
      t_start = UNSET
      if (nout < 2) return
      ref = rates(nout)
      if (.not. (abs(ref) > 0.0_dp)) return
      do i = nout, 1, -1
         if (abs(rates(i) - ref) > rel_tol * abs(ref)) exit
         t_start = centres(i)
      end do
      found = t_start > UNSET / 2.0_dp .and. t_start < centres(nout)
   end subroutine settled_start

   !> The admissible-window analysis of one record. t_min_user / t_max_user, when set (not UNSET), replace the automatic lower and upper ends.
   subroutine analyse_record(t, a, r0, a_max, width, rel_tol, t_min_user, &
                             t_max_user, win, centres, rates, nout)
      real(dp), intent(in) :: t(:), a(:), r0, a_max, width, rel_tol
      real(dp), intent(in) :: t_min_user, t_max_user
      type(record_window_t), intent(out) :: win
      real(dp), allocatable, intent(out) :: centres(:), rates(:)
      integer, intent(out) :: nout
      integer :: iend, e, m, i
      real(dp) :: b, t_start
      logical :: found

      call linear_end(a, a_max * r0, iend, win%bounded)
      if (t_max_user > UNSET / 2.0_dp) then
         iend = 0
         do i = 1, size(t)
            if (t(i) <= t_max_user) iend = i
         end do
         iend = max(iend, 1)
      end if
      call sliding_rates(t, a, iend, width, centres, rates, nout)
      call settled_start(centres, rates, nout, rel_tol, t_start, found)
      win%settled = found
      if (t_min_user > UNSET / 2.0_dp) then
         t_start = t_min_user
      else if (.not. found) then
         ! No settled start: fall back on the second half of the linear part, and say so through settled = .false.
         t_start = t(1) + 0.5_dp * (t(iend) - t(1))
      end if

      call fit_log_linear(t, a, t_start, t(iend), 3, win%rate, b, win%r2, e, m)
      win%ierr = e
      win%n = m
      if (e /= 0) return
      win%t_lo = huge(1.0_dp); win%t_hi = -huge(1.0_dp)
      do i = 1, iend
         if (t(i) >= t_start) then
            if (t(i) < win%t_lo) then
               win%t_lo = t(i); win%a_lo = a(i)
            end if
            win%t_hi = t(i); win%a_hi = a(i)
         end if
      end do
      win%efolds = win%rate * (win%t_hi - win%t_lo)
   end subroutine analyse_record

end module record_mod
