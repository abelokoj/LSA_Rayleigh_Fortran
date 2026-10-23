! Original and corrected mode selection of fort_process_LSA on one Gram
! matrix (standalone, system LAPACK).
module lsa_fix
  implicit none
  integer, parameter :: dp = kind(1.0d0)
contains
  subroutine process(G, n, rcond, delta, fixed, lam_sel, amp_sel, nkeep)
    integer, intent(in) :: n
    real(dp), intent(in) :: G(n, n), rcond, delta
    logical, intent(in) :: fixed
    complex(dp), intent(out) :: lam_sel
    real(dp), intent(out) :: amp_sel
    integer, intent(out) :: nkeep
    integer :: nm1, i, j, info, max_index, lwork
    real(dp), allocatable :: A(:,:), ACOPY(:,:), U(:,:), VT(:,:), S(:)
    real(dp), allocatable :: B(:), UTB(:), AVEC(:), WORK(:), RWORK(:)
    complex(dp), allocatable :: SC(:,:), BC(:,:), ALPHA(:), BETA(:)
    complex(dp), allocatable :: VL(:,:), VR(:,:), LAM(:), WORKC(:)
    complex(dp), allocatable :: XC(:,:), CB(:), AX(:)
    integer, allocatable :: IPIV(:)
    real(dp), allocatable :: AMP(:)
    real(dp) :: cutoff, pnorm, ln_abs, ph, best
    nm1 = n - 1
    allocate(A(nm1,nm1), ACOPY(nm1,nm1), U(nm1,nm1), VT(nm1,nm1))
    allocate(S(nm1), B(nm1), UTB(nm1), AVEC(nm1), WORK(20*nm1))
    allocate(RWORK(20*nm1), SC(nm1,nm1), BC(nm1,nm1), ALPHA(nm1))
    allocate(BETA(nm1), VL(nm1,nm1), VR(nm1,nm1), LAM(nm1))
    allocate(WORKC(20*nm1), XC(nm1,nm1), CB(nm1), AX(nm1), IPIV(nm1))
    allocate(AMP(nm1))
    A = G(1:nm1, 1:nm1)
    B = G(1:nm1, n)
    ACOPY = A
    lwork = 20*nm1
    call dgesvd('A', 'A', nm1, nm1, A, nm1, S, U, nm1, VT, nm1, WORK, &
                lwork, info)
    ! --- R2: relative cutoff (original: keep every S(i) > 0) ---
    if (fixed) then
      cutoff = rcond*S(1)
    else
      cutoff = 0.0_dp
    end if
    nkeep = 0
    do i = 1, nm1
      UTB(i) = dot_product(U(:,i), B)
      if (S(i) .gt. cutoff) then
        UTB(i) = UTB(i)/S(i)
        nkeep = nkeep + 1
      else
        UTB(i) = 0.0_dp
      end if
    end do
    do i = 1, nm1
      AVEC(i) = dot_product(VT(:,i), UTB)
    end do
    SC = (0.0_dp, 0.0_dp)
    BC = (0.0_dp, 0.0_dp)
    do i = 1, nm1
      if (i .gt. 1) SC(i, i-1) = (1.0_dp, 0.0_dp)
      SC(i, nm1) = cmplx(AVEC(i), 0.0_dp, dp)
      BC(i, i) = (1.0_dp, 0.0_dp)
    end do
    call zggev('V', 'V', nm1, SC, nm1, BC, nm1, ALPHA, BETA, VL, nm1, VR, &
               nm1, WORKC, lwork, RWORK, info)
    LAM = ALPHA/BETA
    if (.not. fixed) then
      max_index = 1
      do i = 1, nm1
        if (real(LAM(i)) .gt. real(LAM(max_index))) max_index = i
      end do
    else
      ! --- R3: Ritz amplitudes, X c = e_1, amp_j = |c_j| sqrt(x_j^H A x_j)
      XC = VR
      CB = (0.0_dp, 0.0_dp)
      CB(1) = (1.0_dp, 0.0_dp)
      call zgesv(nm1, 1, XC, nm1, IPIV, CB, nm1, info)
      do j = 1, nm1
        AX = matmul(cmplx(ACOPY, 0.0_dp, dp), VR(:,j))
        pnorm = sqrt(max(real(dot_product(VR(:,j), AX)), 0.0_dp))
        AMP(j) = abs(CB(j))*pnorm
      end do
      ! growing and near-stationary: ln|lam| > delta, |arg lam| <= 0.1 ln|lam|
      max_index = 0
      best = -1.0_dp
      do j = 1, nm1
        ln_abs = log(abs(LAM(j)))
        ph = abs(atan2(aimag(LAM(j)), real(LAM(j))))
        if (ln_abs .gt. delta .and. ph .le. 0.1_dp*ln_abs) then
          if (AMP(j) .gt. best) then
            best = AMP(j)
            max_index = j
          end if
        end if
      end do
      if (max_index .eq. 0) then
        print *, "no growing near-stationary mode; largest amplitude used"
        max_index = maxloc(AMP, 1)
      end if
    end if
    lam_sel = LAM(max_index)
    if (fixed) then
      amp_sel = AMP(max_index)
    else
      amp_sel = -1.0_dp
    end if
    if (.not. fixed) then
      print '(a)', '  spectrum (|lam|):'
      print '(5f10.5)', abs(LAM)
    end if
  end subroutine process
end module lsa_fix

program test_fix
  use lsa_fix
  implicit none
  integer :: n
  real(dp), allocatable :: G(:,:), g1(:)
  real(dp) :: amp, tau, rc, dl
  character(len=256) :: fname, arg
  complex(dp) :: lam
  integer :: nk
  tau = 0.25_dp
  call get_command_argument(1, fname)
  call get_command_argument(2, arg)
  read(arg, *) n
  allocate(G(n, n), g1(n*n))
  open(10, file=trim(fname))
  read(10, *) g1
  close(10)
  G = reshape(g1, [n, n])
  call get_command_argument(3, arg)
  read(arg, *) rc
  call get_command_argument(4, arg)
  read(arg, *) dl
  call process(G, n, rc, dl, .false., lam, amp, nk)
  print '(a,2f10.6,a,f9.5,a,i3)', 'original : lam = ', lam, &
    '  s = ', log(abs(lam))/tau, '  kept ', nk
  call process(G, n, rc, dl, .true., lam, amp, nk)
  print '(a,2f10.6,a,f9.5,a,i3)', 'corrected: lam = ', lam, &
    '  s = ', log(abs(lam))/tau, '  kept ', nk
end program test_fix
