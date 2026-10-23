!> Tests of companion_mod (svd_and_diagonalize): eq. (67) of Prof. Sussman's unpublished draft description of the restarting method solved by SVD, the companion matrix S, and its eigenvalues, on Krylov sequences built from known eigenvalues.
program test_companion
   use testing_mod
   use companion_mod, only: solve_normal_svd, solve_weighted_lstsq, &
                            assemble_companion, diagonalize_companion, &
                            ritz_amplitudes, weighted_residual
   implicit none

   call test_known_spectrum()
   call test_companion_structure()
   call test_weights()
   call test_noise_and_rank()
   call summary('test_companion')

contains

   !> Five modes: three real eigenvalues (1.08 growing, 0.9, 0.5) and a complex pair 0.8 exp(+-0.6i). With N - 1 = 5 vectors the companion matrix reproduces them exactly.
   subroutine test_known_spectrum()
      integer, parameter :: M = 300, NV = 6
      real(dp) :: B(M, NV - 1), vN(M), w(M), a(NV - 1), a2(NV - 1), sv(NV - 1)
      real(dp) :: S(NV - 1, NV - 1), amp(NV - 1), pn(NV - 1)
      complex(dp) :: mu(NV - 1), X(NV - 1, NV - 1), want(NV - 1)
      real(dp) :: V(M, NV)
      integer :: rank, rank2, ierr, j

      want = [cmplx(1.08_dp, 0.0_dp, dp), cmplx(0.9_dp, 0.0_dp, dp), &
              0.8_dp * exp(cmplx(0.0_dp, 0.6_dp, dp)), &
              0.8_dp * exp(cmplx(0.0_dp, -0.6_dp, dp)), &
              cmplx(0.5_dp, 0.0_dp, dp)]
      call krylov_sequence(want, [1.0_dp, 0.3_dp, 0.2_dp, 0.05_dp], V)
      B = V(:, 1:NV - 1); vN = V(:, NV); w = 1.0_dp

      call section('eq. (67) by SVD on a sequence with a known spectrum')
      call solve_normal_svd(B, vN, w, 1.0e-14_dp, a, sv, rank, ierr)
      call check_int(ierr, 0, 'SVD solve succeeds')
      call check_int(rank, NV - 1, 'B^T W B has full rank')
      call check(weighted_residual(B, vN, w, a) < 1.0e-8_dp, &
                 'v_N is reproduced by B a (residual < 1e-8)')
      call solve_weighted_lstsq(B, vN, w, 1.0e-14_dp, a2, rank2, ierr)
      call check(sqrt(sum((a - a2)**2)) / sqrt(sum(a2**2)) < 1.0e-6_dp, &
                 'normal-equation SVD agrees with the direct sqrt(W)B solve')

      call assemble_companion(a, S)
      call diagonalize_companion(S, mu, X, ierr)
      call check_int(ierr, 0, 'dgeev succeeds')
      call check_close(real(mu(1), dp), 1.08_dp, 1.0e-8_dp, &
                       'largest |mu| is the growing eigenvalue 1.08')
      call check_close(aimag(mu(1)), 0.0_dp, 1.0e-8_dp, '... and is real')
      call check_close(abs(mu(2)), 0.9_dp, 1.0e-8_dp, 'second |mu| = 0.9')
      call check_close(abs(mu(3)), 0.8_dp, 1.0e-8_dp, 'complex pair |mu| = 0.8')
      call check_close(abs(aimag(mu(3))), 0.8_dp * sin(0.6_dp), 1.0e-8_dp, &
                       '... with the right phase')
      call check_close(real(mu(3) - conjg(mu(4)), dp), 0.0_dp, 1.0e-10_dp, &
                       'the pair is complex-conjugate')
      call check_close(abs(mu(5)), 0.5_dp, 1.0e-8_dp, 'smallest |mu| = 0.5')
      do j = 1, NV - 2
         call check(abs(mu(j)) >= abs(mu(j + 1)), 'sorted by decreasing |mu|')
      end do

      call ritz_amplitudes(B, X, w, amp, pn, ierr)
      call check_int(ierr, 0, 'Ritz amplitudes computed')
      call check(maxloc(amp, 1) == 1, &
                 'largest amplitude belongs to the mode planted largest')
   end subroutine test_known_spectrum

   subroutine test_companion_structure()
      real(dp) :: S(4, 4)
      call section('companion matrix layout')
      call assemble_companion([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], S)
      call check(all(S(:, 4) == [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]), &
                 'last column holds a')
      call check(S(2, 1) == 1.0_dp .and. S(3, 2) == 1.0_dp .and. &
                 S(4, 3) == 1.0_dp, 'ones on the subdiagonal')
      call check(S(1, 1) == 0.0_dp .and. S(1, 2) == 0.0_dp .and. &
                 S(2, 2) == 0.0_dp .and. S(3, 1) == 0.0_dp, 'zeros elsewhere')
   end subroutine test_companion_structure

   !> Rows with zero weight must not influence a: corrupt them and the answer is unchanged.
   subroutine test_weights()
      integer, parameter :: M = 200, NV = 4
      real(dp) :: V(M, NV), B(M, NV - 1), vN(M), w(M), a1(NV - 1), a2(NV - 1)
      real(dp) :: sv(NV - 1)
      integer :: rank, ierr
      complex(dp) :: want(NV - 1)
      call section('weights')
      want = [cmplx(1.05_dp, 0.0_dp, dp), cmplx(0.7_dp, 0.0_dp, dp), &
              cmplx(0.3_dp, 0.0_dp, dp)]
      call krylov_sequence(want, [1.0_dp, 0.5_dp], V)
      B = V(:, 1:NV - 1); vN = V(:, NV)
      w = 1.0_dp; w(M/2 + 1:) = 0.0_dp
      call solve_normal_svd(B, vN, w, 1.0e-14_dp, a1, sv, rank, ierr)
      B(M/2 + 1:, :) = B(M/2 + 1:, :) + 3.0_dp
      vN(M/2 + 1:) = -7.0_dp
      call solve_normal_svd(B, vN, w, 1.0e-14_dp, a2, sv, rank, ierr)
      call check(maxval(abs(a1 - a2)) < 1.0e-9_dp, &
                 'zero-weight rows do not affect a')
   end subroutine test_weights

   !> More vectors than modes: B^T W B is rank deficient, and the cutoff keeps the solve finite and the leading eigenvalue right.
   subroutine test_noise_and_rank()
      integer, parameter :: M = 250, NV = 9
      real(dp) :: V(M, NV), B(M, NV - 1), vN(M), w(M), a(NV - 1), sv(NV - 1)
      real(dp) :: S(NV - 1, NV - 1)
      complex(dp) :: mu(NV - 1), X(NV - 1, NV - 1), want(3)
      integer :: rank, ierr
      call section('rank-deficient sequence')
      want = [cmplx(1.1_dp, 0.0_dp, dp), cmplx(0.6_dp, 0.0_dp, dp), &
              cmplx(0.2_dp, 0.0_dp, dp)]
      call krylov_sequence(want, [1.0_dp, 0.4_dp], V)
      B = V(:, 1:NV - 1); vN = V(:, NV); w = 1.0_dp
      call solve_normal_svd(B, vN, w, 1.0e-12_dp, a, sv, rank, ierr)
      call check_int(ierr, 0, 'SVD solve succeeds')
      call check(rank < NV - 1, 'rank is reduced by the cutoff')
      call check(all(abs(a) < huge(1.0_dp)), 'a is finite')
      call assemble_companion(a, S)
      call diagonalize_companion(S, mu, X, ierr)
      call check(minval(abs(mu - cmplx(1.1_dp, 0.0_dp, dp))) < 1.0e-6_dp, &
                 'the growing eigenvalue 1.1 is among the eigenvalues of S')
   end subroutine test_noise_and_rank

   !> v_k = sum_j c_j Re( mu_j^(k-1) phi_j ), with independent pseudo-random real shapes; a complex-conjugate pair of eigenvalues is represented by the real and imaginary parts of one complex mode, so the sequence is real.
   subroutine krylov_sequence(mus, cvals, V)
      complex(dp), intent(in) :: mus(:)
      real(dp), intent(in)    :: cvals(:)
      real(dp), intent(out)   :: V(:,:)
      integer :: m, nv, nm, j, k, i
      complex(dp), allocatable :: phi(:,:)
      real(dp) :: c
      m = size(V, 1); nv = size(V, 2); nm = size(mus)
      allocate(phi(m, nm))
      do j = 1, nm
         do i = 1, m
            phi(i, j) = cmplx(sin(1.3_dp * i * j + 0.7_dp * j), &
                              cos(0.9_dp * i + 1.1_dp * j * j), dp)
         end do
         ! a real eigenvalue needs a real shape; a conjugate pair shares conjugate shapes
         if (aimag(mus(j)) == 0.0_dp) phi(:, j) = cmplx(real(phi(:, j), dp), 0.0_dp, dp)
         if (j > 1) then
            if (abs(mus(j) - conjg(mus(j - 1))) < 1.0e-14_dp .and. &
                aimag(mus(j)) /= 0.0_dp) phi(:, j) = conjg(phi(:, j - 1))
         end if
      end do
      V = 0.0_dp
      do k = 1, nv
         do j = 1, nm
            c = cvals(min(j, size(cvals)))
            V(:, k) = V(:, k) + c * real(mus(j)**(k - 1) * phi(:, j), dp)
         end do
      end do
   end subroutine krylov_sequence

end program test_companion
