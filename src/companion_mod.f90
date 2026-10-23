!> Companion-matrix (Krylov) eigenanalysis of a snapshot sequence, following Section "Equivalency of the above analysis to the least squares approach and a new generalization" of Prof. Sussman's unpublished draft description of the restarting method and Schmid (2010), eqs. (2.5)-(2.8).
!>
!> Given a Krylov sequence v_1, v_2 = A v_1, ..., v_N = A^(N-1) v_1 (each of length M), write B = V_1^(N-1) = [v_1 ... v_(N-1)] (M x (N-1)). The last vector is expanded in the earlier ones, v_N = B a + r, and the coefficients a minimise the weighted residual J(a) = sum_i w_i ((B a)_i - v_(N,i))^2. Setting dJ/da_k = 0 gives the normal equations of the article's eq. (67):
!>
!>    (B^T W B) a = B^T W v_N,    W = diag(w_1, ..., w_M),
!>
!> whose matrix and right-hand side are weighted inner products of the Krylov vectors among themselves, (B^T W B)_ij = (v_i, v_j)_w and (B^T W v_N)_i = (v_i, v_N)_w (eqs. 69-70). The system is solved here by the singular value decomposition of the (N-1) x (N-1) matrix G = B^T W B, i.e. by its Moore-Penrose pseudo-inverse (Remark 1 of the article): a = V_G S_G^+ U_G^T (B^T W v_N), with singular values below rcond * s_1 discarded.
!>
!> The coefficients fill the last column of the companion matrix (the article's eq. for S, Schmid eq. 2.8): ones on the subdiagonal, a in the last column, zeros elsewhere. Then A V_1^(N-1) = V_1^(N-1) S + r e_(N-1)^T, so with r neglected the eigenvalues of S approximate the dominant eigenvalues of A. S is diagonalised, S = X Gamma X^(-1), with eigenvalues sorted by decreasing magnitude, and the Ritz vectors (approximate eigenvectors of A) are the columns of P = V_1^(N-1) X.
!>
!> Note on conditioning. Forming G squares the condition number of B, so singular values of B below about sqrt(eps) * s_1(B) are lost; the pseudo-inverse with rcond then gives the minimum-norm solution in the remaining subspace. This is the formulation the article asks for (it never forms the M x M matrix U and needs only inner products of the Krylov vectors); solve_weighted_lstsq below solves the same problem directly from sqrt(W) B and is used as an independent check.
module companion_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: dp, solve_normal_svd, solve_weighted_lstsq, assemble_companion, &
             diagonalize_companion, ritz_amplitudes, weighted_residual

contains

   !> Solve (B^T W B) a = B^T W vN (eq. 67 of Prof. Sussman's unpublished draft description of the restarting method) by the SVD of G = B^T W B.
   !> B(M, n) holds v_1..v_n with n = N-1; vN(M) is v_N; w(M) are the weights w_i >= 0.
   !> On return: a(n); sv(n) the singular values of G (descending); rank the number kept (sv > rcond * sv(1)); ierr = 0 success, 1 bad sizes, 2 G is zero, 10+info dgesvd failure.
   subroutine solve_normal_svd(B, vN, w, rcond, a, sv, rank, ierr)
      real(dp), intent(in)  :: B(:,:), vN(:), w(:), rcond
      real(dp), intent(out) :: a(:), sv(:)
      integer,  intent(out) :: rank, ierr
      real(dp), allocatable :: G(:,:), WB(:,:), rhs(:), U(:,:), VT(:,:), work(:)
      real(dp) :: coef
      integer :: m, n, i, lwork, info

      m = size(B, 1); n = size(B, 2)
      a = 0.0_dp; sv = 0.0_dp; rank = 0; ierr = 0
      if (n < 1 .or. size(vN) /= m .or. size(w) /= m .or. size(a) /= n &
          .or. size(sv) /= n) then
         ierr = 1; return
      end if

      ! G = B^T W B and rhs = B^T W vN: weighted inner products of the Krylov vectors (eqs. 69-70).
      allocate(WB(m, n), G(n, n), rhs(n))
      do i = 1, n
         WB(:, i) = w * B(:, i)
      end do
      G = matmul(transpose(B), WB)
      rhs = matmul(transpose(WB), vN)

      ! SVD of the small symmetric matrix G = U S V^T.
      allocate(U(n, n), VT(n, n))
      lwork = max(1, 10 * n)
      allocate(work(lwork))
      call dgesvd('A', 'A', n, n, G, n, sv, U, n, VT, n, work, lwork, info)
      if (info /= 0) then
         ierr = 10 + info; return
      end if
      if (.not. (sv(1) > 0.0_dp)) then
         ierr = 2; return
      end if

      ! Pseudo-inverse: a = sum over kept i of v_i (u_i . rhs) / s_i.
      do i = 1, n
         if (sv(i) <= rcond * sv(1)) exit
         rank = i
         coef = dot_product(U(:, i), rhs) / sv(i)
         a = a + coef * VT(i, :)
      end do
   end subroutine solve_normal_svd

   !> Independent check: the same weighted least-squares problem solved directly, min || sqrt(W) (B a - vN) ||_2, by LAPACK dgelss (an SVD of sqrt(W) B), without forming G. rcond as above, but relative to the singular values of sqrt(W) B.
   subroutine solve_weighted_lstsq(B, vN, w, rcond, a, rank, ierr)
      real(dp), intent(in)  :: B(:,:), vN(:), w(:), rcond
      real(dp), intent(out) :: a(:)
      integer,  intent(out) :: rank, ierr
      real(dp), allocatable :: Bs(:,:), rhs(:), s(:), work(:)
      integer :: m, n, i, lwork, info

      m = size(B, 1); n = size(B, 2)
      a = 0.0_dp; rank = 0; ierr = 0
      if (n < 1 .or. size(vN) /= m .or. size(w) /= m .or. size(a) /= n) then
         ierr = 1; return
      end if
      allocate(Bs(m, n), rhs(max(m, n)), s(min(m, n)))
      do i = 1, n
         Bs(:, i) = sqrt(max(w, 0.0_dp)) * B(:, i)
      end do
      rhs = 0.0_dp
      rhs(1:m) = sqrt(max(w, 0.0_dp)) * vN
      lwork = max(1, 3 * min(m, n) + max(2 * min(m, n), max(m, n), 1) + 64 * n)
      allocate(work(lwork))
      call dgelss(m, n, 1, Bs, m, rhs, max(m, n), s, rcond, rank, work, &
                  lwork, info)
      if (info /= 0) then
         ierr = 10 + info; return
      end if
      a = rhs(1:n)
   end subroutine solve_weighted_lstsq

   !> Companion matrix of the article's eq. for S (Schmid 2010 eq. 2.8): S(i+1, i) = 1, S(:, n) = a, zero elsewhere.
   subroutine assemble_companion(a, S)
      real(dp), intent(in)  :: a(:)
      real(dp), intent(out) :: S(:,:)
      integer :: n, i
      n = size(a)
      S = 0.0_dp
      do i = 1, n - 1
         S(i + 1, i) = 1.0_dp
      end do
      S(:, n) = a
   end subroutine assemble_companion

   !> Diagonalise S = X Gamma X^(-1) with LAPACK dgeev. Eigenvalues mu(n) (complex) and right eigenvectors X(n, n) (complex, unit 2-norm columns) are returned sorted by decreasing |mu| (the article's eq. for the spectrum of S). ierr = 10+info on dgeev failure.
   subroutine diagonalize_companion(S, mu, X, ierr)
      real(dp), intent(in)     :: S(:,:)
      complex(dp), intent(out) :: mu(:), X(:,:)
      integer, intent(out)     :: ierr
      real(dp), allocatable :: Sc(:,:), wr(:), wi(:), vr(:,:), work(:), mag(:)
      real(dp) :: vl(1, 1)
      complex(dp), allocatable :: mu0(:), X0(:,:)
      integer, allocatable :: order(:)
      integer :: n, j, lwork, info, i, k, tmp

      n = size(S, 1); ierr = 0
      allocate(Sc(n, n))
      Sc = S
      allocate(wr(n), wi(n), vr(n, n), mu0(n), X0(n, n), mag(n), order(n))
      lwork = max(1, 8 * n)
      allocate(work(lwork))
      call dgeev('N', 'V', n, Sc, n, wr, wi, vl, 1, vr, n, work, lwork, info)
      if (info /= 0) then
         ierr = 10 + info; return
      end if

      ! dgeev stores a complex-conjugate pair as two real columns (re, im) at j, j+1, with wi(j) > 0.
      j = 1
      do while (j <= n)
         if (wi(j) == 0.0_dp .or. j == n) then
            mu0(j) = cmplx(wr(j), 0.0_dp, dp)
            X0(:, j) = cmplx(vr(:, j), 0.0_dp, dp)
            j = j + 1
         else
            mu0(j) = cmplx(wr(j), wi(j), dp)
            mu0(j + 1) = cmplx(wr(j + 1), wi(j + 1), dp)
            X0(:, j) = cmplx(vr(:, j), vr(:, j + 1), dp)
            X0(:, j + 1) = conjg(X0(:, j))
            j = j + 2
         end if
      end do
      do j = 1, n
         X0(:, j) = X0(:, j) / max(sqrt(sum(abs(X0(:, j))**2)), tiny(1.0_dp))
      end do

      ! Stable insertion sort of the indices by decreasing |mu|.
      mag = abs(mu0)
      order = [(i, i = 1, n)]
      do i = 2, n
         k = i
         do while (k > 1)
            if (mag(order(k - 1)) >= mag(order(k))) exit
            tmp = order(k); order(k) = order(k - 1); order(k - 1) = tmp
            k = k - 1
         end do
      end do
      mu = mu0(order)
      X = X0(:, order)
   end subroutine diagonalize_companion

   !> Amplitude of each Ritz vector in the first snapshot. Since B = P X^(-1), v_1 = B e_1 = P b with X b = e_1; solved by LAPACK zgesv. amp(j) = |b_j| * ||P(:, j)||_w, the size of mode j's contribution to v_1. ierr = 20+info if X is singular.
   subroutine ritz_amplitudes(B, X, w, amp, pnorm, ierr)
      real(dp), intent(in)    :: B(:,:), w(:)
      complex(dp), intent(in) :: X(:,:)
      real(dp), intent(out)   :: amp(:), pnorm(:)
      integer, intent(out)    :: ierr
      complex(dp), allocatable :: Xc(:,:), rhs(:,:), P(:)
      integer, allocatable :: ipiv(:)
      integer :: n, j, info

      n = size(X, 1); ierr = 0
      allocate(Xc(n, n), P(size(B, 1)))
      Xc = X
      allocate(rhs(n, 1), ipiv(n))
      rhs = (0.0_dp, 0.0_dp)
      rhs(1, 1) = (1.0_dp, 0.0_dp)
      call zgesv(n, 1, Xc, n, ipiv, rhs, n, info)
      if (info /= 0) then
         ierr = 20 + info; amp = 0.0_dp; pnorm = 0.0_dp; return
      end if
      do j = 1, n
         P = cmplx(matmul(B, real(X(:, j), dp)), &
                   matmul(B, aimag(X(:, j))), dp)
         pnorm(j) = sqrt(sum(w * abs(P)**2))
         amp(j) = abs(rhs(j, 1)) * pnorm(j)
      end do
   end subroutine ritz_amplitudes

   !> Relative weighted residual || vN - B a ||_w / || vN ||_w (the neglected r of the article, measured).
   function weighted_residual(B, vN, w, a) result(rel)
      real(dp), intent(in) :: B(:,:), vN(:), w(:), a(:)
      real(dp) :: rel, den
      real(dp), allocatable :: r(:)
      r = vN - matmul(B, a)
      den = sqrt(sum(w * vN**2))
      if (den > 0.0_dp) then
         rel = sqrt(sum(w * r**2)) / den
      else
         rel = 0.0_dp
      end if
   end function weighted_residual

end module companion_mod
