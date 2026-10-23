!> Write genuine AMReX plotfiles, for testing the reader.
!>
!> The production reader consumes plotfiles written by the AMReX solver.  We cannot run the solver inside a unit test, so we synthesise plotfiles in the real on-disk format and read them back through the same code path.  That exercises the parser, the FAB offsets, the endianness handling and the multi-grid scatter, and fails loudly if any format assumption is wrong.
!>
!> Files produced here are byte-compatible with the ones the Python tree writes (tests/synthetic_plotfile.py) and are readable by yt, which is how the format was cross-checked independently of this code.
module plotfile_writer_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64, int64
   implicit none
   private

   public :: write_plotfile

   ! Real descriptors.  yt decides endianness from the second block: descending "8 7 ... 1" is LITTLE endian, ascending "1 2 ... 8" is BIG.
   character(len=*), parameter :: DESC_LE = &
      '((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))'
   character(len=*), parameter :: DESC_BE = &
      '((8, (64 11 52 0 1 12 0 1023)),(8, (1 2 3 4 5 6 7 8)))'

contains

   !> Write one plotfile directory.
   !>
   !> @param path      directory to create, e.g. ".../plt00100" @param data      (nx, ny, ncomp) @param names     component names, size ncomp @param time      simulation time @param step      step number @param prob_lo, prob_hi   domain extent @param big_endian  write big-endian data (to test the swap path) @param nsplit  number of boxes the domain is cut into along y (default 1), all stored one after another in one Cell_D file at their own offsets, as an MPI run writes them; this exercises the multi-grid scatter of the reader, which a single box does not @param drop_box  if > 0, leave that box out of Cell_H and the data, so that the boxes no longer cover the domain (to test the coverage check)
   !>
   !> The path is quoted for the shell, so that a scratch directory containing a space works; the tests use one deliberately.
   subroutine write_plotfile(path, data, names, time, step, prob_lo, prob_hi, &
                             coord_sys, big_endian, nsplit, drop_box)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: data(:,:,:)
      character(len=*), intent(in) :: names(:)
      real(dp), intent(in) :: time
      integer,  intent(in) :: step
      real(dp), intent(in) :: prob_lo(2), prob_hi(2)
      integer,  intent(in), optional :: coord_sys
      logical,  intent(in), optional :: big_endian
      integer,  intent(in), optional :: nsplit, drop_box

      integer :: nx, ny, ncomp, u, i, j, c, cs, nb, ib, drop, nkept
      integer :: jlo(64), jhi(64)
      logical :: be
      character(len=256) :: box, leveldir, gbox(64)
      real(dp) :: dx(2)
      real(dp), allocatable :: buf(:)
      integer(int64) :: offset(64), pos

      nx = size(data, 1); ny = size(data, 2); ncomp = size(data, 3)
      cs = 0
      if (present(coord_sys)) cs = coord_sys
      be = .false.
      if (present(big_endian)) be = big_endian
      nb = 1
      if (present(nsplit)) nb = max(1, min(nsplit, ny, 64))
      drop = 0
      if (present(drop_box)) drop = drop_box

      dx(1) = (prob_hi(1) - prob_lo(1)) / real(nx, dp)
      dx(2) = (prob_hi(2) - prob_lo(2)) / real(ny, dp)

      call execute_command_line("mkdir -p '" // trim(path) // "/Level_0'")
      leveldir = trim(path) // '/Level_0'

      write(box, '(a,i0,a,i0,a)') '((0,0) (', nx-1, ',', ny-1, ') (0,0))'
      nkept = 0
      do ib = 1, nb
         jlo(ib) = ((ib - 1) * ny) / nb
         jhi(ib) = (ib * ny) / nb - 1
         write(gbox(ib), '(a,i0,a,i0,a,i0,a)') '((0,', jlo(ib), ') (', nx-1, &
            ',', jhi(ib), ') (0,0))'
         if (ib /= drop) nkept = nkept + 1
      end do

      ! ---- global Header ------------------------------------------------
      open(newunit=u, file=trim(path)//'/Header', status='replace', &
           action='write')
      write(u,'(a)') 'HyperCLaw-V1.1'
      write(u,'(i0)') ncomp
      do i = 1, ncomp
         write(u,'(a)') trim(names(i))
      end do
      write(u,'(i0)') 2                            ! dimensionality
      write(u,'(es24.16)') time
      write(u,'(i0)') 0                            ! finest level
      write(u,'(es24.16,1x,es24.16)') prob_lo(1), prob_lo(2)
      write(u,'(es24.16,1x,es24.16)') prob_hi(1), prob_hi(2)
      write(u,'(a)') ''                            ! ref ratios (none)
      write(u,'(a)') trim(box)                     ! index space
      write(u,'(i0)') step                         ! level steps
      write(u,'(es24.16,1x,es24.16)') dx(1), dx(2) ! cell sizes
      write(u,'(i0)') cs                           ! coord_sys
      write(u,'(i0)') 0                            ! bwidth
      write(u,'(i0,1x,i0,1x,es24.16)') 0, nkept, time  ! level, ngrids, time
      write(u,'(i0)') step
      do ib = 1, nb
         if (ib == drop) cycle
         write(u,'(es24.16,1x,es24.16)') prob_lo(1), prob_hi(1)
         write(u,'(es24.16,1x,es24.16)') prob_lo(2) + dx(2)*jlo(ib), &
                                         prob_lo(2) + dx(2)*(jhi(ib)+1)
      end do
      write(u,'(a)') 'Level_0/Cell'
      close(u)

      ! ---- FAB data ------------------------------------------------------ One FAB per box, each an ASCII header line then each component in column-major order, all in one file at increasing offsets.
      open(newunit=u, file=trim(leveldir)//'/Cell_D_00000', &
           status='replace', access='stream', form='unformatted', &
           action='write')
      pos = 0_int64
      do ib = 1, nb
         if (ib == drop) cycle
         offset(ib) = pos
         block
            character(len=512) :: fabhdr
            integer :: n, nloc
            if (be) then
               fabhdr = 'FAB ' // DESC_BE // trim(gbox(ib))
            else
               fabhdr = 'FAB ' // DESC_LE // trim(gbox(ib))
            end if
            n = len_trim(fabhdr)
            write(fabhdr(n+1:), '(1x,i0)') ncomp
            n = len_trim(fabhdr)
            write(u) fabhdr(1:n), char(10)
            pos = pos + int(n, int64) + 1_int64
            nloc = nx * (jhi(ib) - jlo(ib) + 1)
            if (allocated(buf)) deallocate(buf)
            allocate(buf(nloc))
            do c = 1, ncomp
               do j = jlo(ib) + 1, jhi(ib) + 1
                  do i = 1, nx
                     buf(i + (j-1-jlo(ib))*nx) = data(i, j, c)
                  end do
               end do
               if (be) call byteswap8_w(buf)
               write(u) buf
               pos = pos + int(nloc, int64) * 8_int64
            end do
         end block
      end do
      close(u)

      ! ---- Cell_H --------------------------------------------------------
      open(newunit=u, file=trim(leveldir)//'/Cell_H', status='replace', &
           action='write')
      write(u,'(i0)') 1                    ! version
      write(u,'(i0)') 1                    ! how
      write(u,'(i0)') ncomp
      write(u,'(i0)') 0                    ! nghost
      write(u,'(a,i0,a)') '(', nkept, ' 0'
      do ib = 1, nb
         if (ib /= drop) write(u,'(a)') trim(gbox(ib))
      end do
      write(u,'(a)') ')'
      write(u,'(i0)') nkept
      do ib = 1, nb
         if (ib /= drop) write(u,'(a,i0)') 'FabOnDisk: Cell_D_00000 ', &
            offset(ib)
      end do
      write(u,'(a)') ''
      write(u,'(i0,a,i0)') nkept, ',', ncomp
      do ib = 1, nkept
         do c = 1, ncomp
            write(u,'(es24.16,a)', advance='no') minval(data(:,:,c)), ','
         end do
         write(u,'(a)') ''
      end do
      write(u,'(a)') ''
      write(u,'(i0,a,i0)') nkept, ',', ncomp
      do ib = 1, nkept
         do c = 1, ncomp
            write(u,'(es24.16,a)', advance='no') maxval(data(:,:,c)), ','
         end do
         write(u,'(a)') ''
      end do
      close(u)
   end subroutine write_plotfile

   subroutine byteswap8_w(a)
      real(dp), intent(inout) :: a(:)
      integer(int64) :: t, v
      integer :: i
      do i = 1, size(a)
         t = transfer(a(i), t)
         v = 0_int64
         v = ior(ishft(iand(t, int(z'00000000000000FF', int64)), 56), v)
         v = ior(ishft(iand(t, int(z'000000000000FF00', int64)), 40), v)
         v = ior(ishft(iand(t, int(z'0000000000FF0000', int64)), 24), v)
         v = ior(ishft(iand(t, int(z'00000000FF000000', int64)),  8), v)
         v = ior(iand(ishft(t, -8),  int(z'00000000FF000000', int64)), v)
         v = ior(iand(ishft(t, -24), int(z'0000000000FF0000', int64)), v)
         v = ior(iand(ishft(t, -40), int(z'000000000000FF00', int64)), v)
         v = ior(iand(ishft(t, -56), int(z'00000000000000FF', int64)), v)
         a(i) = transfer(v, a(i))
      end do
   end subroutine byteswap8_w

end module plotfile_writer_mod
