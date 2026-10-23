!> Exact linear dispersion relation for the Rayleigh-Plateau instability.
!>
!> Fortran counterpart of plotting/dispersion.py, and a generalisation of the solver's own utility NS_fluids_lib/run2d/rayleigh_capillary_growth.F90, which implements the inviscid Rayleigh (1878) result quoted in Popinet (2009) section 6.5.
!>
!> For an axisymmetric liquid column of radius r0, density rho and surface tension sigma, perturbed as r(z,t) = r0 + eps exp(s t) cos(k z),
!>
!> s^2 = (sigma/(rho r0^3)) * x * I1(x)/I0(x) * (1 - x^2),   x = k r0
!>
!> unstable for x < 1, peaking at x ~= 0.697 (lambda ~= 9.02 r0).
!>
!> The modified Bessel functions use the same truncated power series as the original .F90 so this module reproduces it bit for bit; test_dispersion checks that claim against an independent continued-fraction evaluation.
module dispersion_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use util_mod, only: nan_value
   implicit none
   private

   public :: dp
   public :: mod_bessel_i, growth_rate, growth_rate_squared
   public :: growth_rate_viscous, growth_rate_viscous_exact, ohnesorge
   public :: growth_rate_two_fluid_inviscid, bessel_i1_over_i0
   public :: bessel_k1_over_k0, viscous_residual
   public :: most_unstable_wavenumber, most_unstable_wavelength
   public :: box_length_for_mode
   public :: growth_rate_two_fluid_viscous, two_fluid_viscous_det
   public :: bessel_i_scaled, bessel_k_scaled
   public :: growth_rate_two_fluid_inviscid_confined
   public :: growth_rate_two_fluid_viscous_confined, two_fluid_viscous_confined_det

   real(dp), parameter, public :: PI = 3.141592653589793238462643383279502884_dp

contains

   !> Modified Bessel function of the first kind, I_n(x), n >= 0.
   !>
   !> Truncated power series, matching mod_bessel_first_kind in rayleigh_capillary_growth.F90: I_n(x) = sum_{m=0}^{M} (x/2)^{n+2m} / (m! (m+n)!), with at most 100 terms.
   !>
   !> The sum now stops once a term falls below a quarter of the unit roundoff of the partial sum. Such a term cannot change the sum under round-to-nearest, and every later term is smaller still, so the result is bit for bit the one of the full 100-term loop; what changes is that the factorials are no longer carried to i = 98 and beyond, where m! (m+n)! overflowed to Inf. That overflow was harmless (the term became zero) but it raised the IEEE overflow flag and aborted any build compiled with -ffpe-trap=overflow.
   pure function mod_bessel_i(x, n) result(val)
      real(dp), intent(in) :: x
      integer,  intent(in) :: n
      real(dp) :: val
      real(dp) :: m_fact, mpn_fact, xpower, term
      integer  :: i

      m_fact = 1.0_dp
      mpn_fact = 1.0_dp
      do i = 1, n
         mpn_fact = mpn_fact * real(i, dp)
      end do

      xpower = (0.5_dp * x) ** n
      val = 0.0_dp
      do i = 0, 100
         if (i > 0) then
            m_fact = m_fact * real(i, dp)
            mpn_fact = mpn_fact * real(i + n, dp)
            xpower = xpower * x * x / 4.0_dp
         end if
         term = xpower / (m_fact * mpn_fact)
         val = val + term
         if (i > 0 .and. term <= 0.25_dp * epsilon(1.0_dp) * abs(val)) exit
      end do
   end function mod_bessel_i

   !> s^2(k); negative in the stable band x > 1 (oscillatory, not growing).
   pure function growth_rate_squared(k, r0, sigma, rho) result(s2)
      real(dp), intent(in) :: k, r0, sigma, rho
      real(dp) :: s2, x, i0, i1

      x = k * r0
      if (x <= 0.0_dp) then
         s2 = 0.0_dp
         return
      end if
      i0 = mod_bessel_i(x, 0)
      i1 = mod_bessel_i(x, 1)
      s2 = sigma * i1 * x * (1.0_dp - x * x) / (rho * r0**3 * i0)
   end function growth_rate_squared

   !> Inviscid growth rate; zero rather than imaginary in the stable band.
   pure function growth_rate(k, r0, sigma, rho) result(s)
      real(dp), intent(in) :: k, r0, sigma, rho
      real(dp) :: s, s2

      s2 = growth_rate_squared(k, r0, sigma, rho)
      s = merge(sqrt(s2), 0.0_dp, s2 > 0.0_dp)
   end function growth_rate

   !> Ohnesorge number Oh = mu / sqrt(rho sigma r0).
   pure function ohnesorge(mu, rho, sigma, r0) result(oh)
      real(dp), intent(in) :: mu, rho, sigma, r0
      real(dp) :: oh
      oh = mu / sqrt(rho * sigma * r0)
   end function ohnesorge

   !> APPROXIMATE viscous growth rate (Weber-type, low Ohnesorge number).
   !>
   !> Solves s^2 + 3 nu k^2 s - s_inviscid^2 = 0 for its larger root, with the exact inviscid s_inviscid^2 (Rayleigh, I_1/I_0 kept). This is the long-wave viscous correction of Weber, the form of eq:3_longwave in the dissertation, and NOT Chandrasekhar's full viscous relation; that one is growth_rate_viscous_exact below, which the drivers now use as the viscous benchmark. The two agree closely at the Ohnesorge numbers of this study, which is why this cheaper form is kept for planning run lengths in make_inputs_mod.
   !>
   !> In the stable band, x > 1, s_inviscid^2 < 0 and the function returns the real part of the least-damped root, which is NEGATIVE: -a + sqrt(a^2 + s2) while the discriminant is positive, and -a once it turns negative (a damped oscillation). The previous comment said the result was clamped at zero, which it never was.
   pure function growth_rate_viscous(k, r0, sigma, rho, mu) result(s)
      real(dp), intent(in) :: k, r0, sigma, rho, mu
      real(dp) :: s, nu, a, s2, disc

      nu = mu / rho
      s2 = growth_rate_squared(k, r0, sigma, rho)
      a = 1.5_dp * nu * k * k
      disc = max(a * a + s2, 0.0_dp)
      s = -a + sqrt(disc)
   end function growth_rate_viscous

   !> I_1(z)/I_0(z) for z > 0, by backward recurrence of the ratio continued fraction I_{n+1}/I_n = 1 / (2(n+1)/z + I_{n+2}/I_{n+1}).
   !>
   !> The ratio of the two functions is all the viscous relation needs, and computing it directly avoids the power series, which loses accuracy and overflows at the large arguments l = sqrt(x^2 + s/Oh) reached at small Ohnesorge number. Backward recurrence of this ratio is stable because I_n is the minimal solution of the recurrence. The starting depth grows with z; at z = 0.7 the result agrees with the power series to round-off, and test_dispersion checks it across 0 < z < 200.
   pure function bessel_i1_over_i0(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      integer :: n, nstart
      if (z <= 0.0_dp) then
         q = 0.0_dp; return
      end if
      nstart = 60 + 2 * int(min(z, 1.0e6_dp))
      q = 0.0_dp
      do n = nstart, 0, -1
         q = 1.0_dp / (2.0_dp * real(n + 1, dp) / z + q)
      end do
   end function bessel_i1_over_i0

   !> Residual of Chandrasekhar's viscous dispersion relation for a liquid column in a passive ambient, in units of the capillary time sqrt(rho r0^3 / sigma), in which the kinematic viscosity equals Oh.
   !>
   !> Chandrasekhar, Hydrodynamic and Hydromagnetic Stability (1961), ch. XII, eqs. (158), (159) and (89); dissertation eq:3_chandra. Written with q(z) = I_1(z)/I_0(z) only, using I_1'(z)/I_1(z) = 1/q(z) - 1/z and hence I_1'(x)/I_0(x) = 1 - q(x)/x:
   !>
   !>   s^2 + 2 Oh x^2 s [ (1 - q(x)/x) - 2 x l/(x^2 + l^2) q(x) (1/q(l) - 1/l) ] = (l^2 - x^2)/(l^2 + x^2) x (1 - x^2) q(x),   l^2 = x^2 + s/Oh.
   !>
   !> This is the same form as dispersion.growth_rate_viscous_exact in the Python tree. As Oh -> 0 it reduces to the Rayleigh relation.
   pure function viscous_residual(s, x, oh) result(f)
      real(dp), intent(in) :: s, x, oh
      real(dp) :: f, l, qx, ql, bracket
      l = sqrt(x * x + s / oh)
      qx = bessel_i1_over_i0(x)
      ql = bessel_i1_over_i0(l)
      bracket = (1.0_dp - qx / x) &
              - (2.0_dp * x * l / (x * x + l * l)) * qx * (1.0_dp / ql - 1.0_dp / l)
      f = s * s + 2.0_dp * oh * x * x * s * bracket &
        - ((l * l - x * x) / (l * l + x * x)) * x * (1.0_dp - x * x) * qx
   end function viscous_residual

   !> EXACT viscous growth rate (Chandrasekhar), by bisection on viscous_residual.
   !>
   !> The inviscid rate bounds the root from above, since viscosity can only retard the capillary instability, and zero bounds it from below in the unstable band, so a bracket always exists and no starting guess is needed. Bisection is used for its guarantee rather than its speed: 200 halvings take the bracket below round-off, and the routine is called a handful of times per analysis. Returned in the same units as growth_rate. In the stable band (x >= 1) there is no growing root and 0 is returned; NaN is returned if the bracket fails, which would mean the relation was evaluated outside its range of validity.
   function growth_rate_viscous_exact(k, r0, sigma, rho, mu) result(s)
      real(dp), intent(in) :: k, r0, sigma, rho, mu
      real(dp) :: s
      real(dp) :: x, oh, t_c, lo, hi, mid, flo, fmid, upper
      integer :: it

      x = k * r0
      s = 0.0_dp
      if (x <= 0.0_dp .or. x >= 1.0_dp) return
      t_c = sqrt(rho * r0**3 / sigma)
      oh = ohnesorge(mu, rho, sigma, r0)
      upper = sqrt(x * (1.0_dp - x * x) * bessel_i1_over_i0(x))
      if (oh <= 0.0_dp) then
         s = upper / t_c; return
      end if

      ! The lower end of the bracket is raised until the residual there is negative. Very close to s = 0 the bracket term of the residual is a difference of nearly equal quantities, and at large Oh rounding can make the residual positive at a tiny s although the true value is negative; a fixed lower end of 1e-14 * upper then had both ends of the same sign and returned NaN (for x >= 0.97 at Oh = 10, x >= 0.25 at Oh = 100).
      hi = upper
      lo = 1.0e-14_dp * upper
      flo = viscous_residual(lo, x, oh)
      do while (flo >= 0.0_dp .and. lo < 0.1_dp * hi)
         lo = 10.0_dp * lo
         flo = viscous_residual(lo, x, oh)
      end do
      if (flo * viscous_residual(hi, x, oh) > 0.0_dp) then
         s = nan_value(); return
      end if
      do it = 1, 200
         mid = 0.5_dp * (lo + hi)
         fmid = viscous_residual(mid, x, oh)
         if (flo * fmid <= 0.0_dp) then
            hi = mid
         else
            lo = mid; flo = fmid
         end if
         if (hi - lo <= 4.0_dp * epsilon(1.0_dp) * hi) exit
      end do
      s = 0.5_dp * (lo + hi) / t_c
   end function growth_rate_viscous_exact

   !> K_1(x)/K_0(x) for x > 0. For 0 < x <= 2 from the ascending series of Abramowitz and Stegun 9.6.13 (K_0) and 9.6.11 (K_1); for x > 2 from Steed's continued fraction (Numerical Recipes, bessik, CF2), in which the factor exp(-x) cancels from the ratio, so large arguments neither overflow nor underflow. NaN for x <= 0.
   function bessel_k1_over_k0(x) result(rk)
      real(dp), intent(in) :: x
      real(dp) :: rk
      real(dp) :: k0, k1

      if (x <= 0.0_dp) then
         rk = nan_value(); return
      end if
      if (x > 2.0_dp) then
         rk = k1_over_k0_cf2(x); return
      end if
      call k01_series(x, k0, k1)
      rk = k1 / k0
   end function bessel_k1_over_k0

   !> K_0(x) and K_1(x) for 0 < x <= 2 from the ascending series of Abramowitz and Stegun 9.6.13 (K_0) and 9.6.11 (K_1). Moved here unchanged from bessel_k1_over_k0 (2026-10-22) so that bessel_k_scaled can use the same values; the arithmetic is the same, so K_1/K_0 is unchanged bit for bit.
   subroutine k01_series(x, k0, k1)
      real(dp), intent(in) :: x
      real(dp), intent(out) :: k0, k1
      real(dp), parameter :: EULER = 0.57721566490153286060651209_dp
      real(dp) :: q, term0, term1, psi1, psi2, lg
      integer :: m
      q = 0.25_dp * x * x
      lg = log(0.5_dp * x)
      ! K_0 = -(ln(x/2) + gamma) I_0 + sum_m q^m/(m!)^2 H_m,  H_m harmonic number.
      ! K_1 = 1/x + ln(x/2) I_1 - (x/4) sum_m [psi(m+1) + psi(m+2)] q^m/(m!(m+1)!)
      k0 = -(lg + EULER) * mod_bessel_i(x, 0)
      k1 = 1.0_dp / x + lg * mod_bessel_i(x, 1)
      term0 = 1.0_dp                          ! q^m/(m!)^2
      term1 = 1.0_dp                          ! q^m/(m!(m+1)!)
      psi1 = -EULER                           ! psi(m+1)
      psi2 = 1.0_dp - EULER                   ! psi(m+2)
      do m = 0, 60
         if (m > 0) then
            term0 = term0 * q / real(m * m, dp)
            term1 = term1 * q / real(m * (m + 1), dp)
            psi1 = psi1 + 1.0_dp / real(m, dp)
            psi2 = psi2 + 1.0_dp / real(m + 1, dp)
            k0 = k0 + term0 * (psi1 + EULER)  ! H_m = psi(m+1) + gamma
         end if
         k1 = k1 - 0.25_dp * x * (psi1 + psi2) * term1
         if (m > 0 .and. term1 < epsilon(1.0_dp) * 1.0e-3_dp) exit
      end do
   end subroutine k01_series

   !> Inviscid growth rate with an ambient fluid of density rho_outer.
   !>
   !> Potential flow inside and outside the column, the outer disturbance decaying as K_0(kr), gives s^2 t_c^2 = x (1 - x^2) / [ I_0(x)/I_1(x) + (rho_outer/rho) K_0(x)/K_1(x) ]. The decks carry an ambient phase of density 1.225e-3 against a column of density 1, and the single-fluid relation omits it; this quantifies the omission. The ambient adds inertia and so always lowers the rate. Same form as dispersion.growth_rate_two_fluid_inviscid in the Python tree. Zero in the stable band and at k = 0.
   function growth_rate_two_fluid_inviscid(k, r0, sigma, rho, rho_outer) &
         result(s)
      real(dp), intent(in) :: k, r0, sigma, rho, rho_outer
      real(dp) :: s, x, t_c, s2
      x = k * r0
      s = 0.0_dp
      if (x <= 0.0_dp .or. x >= 1.0_dp) return
      t_c = sqrt(rho * r0**3 / sigma)
      s2 = x * (1.0_dp - x * x) / (1.0_dp / bessel_i1_over_i0(x) &
           + (rho_outer / rho) / bessel_k1_over_k0(x))
      s = sqrt(max(s2, 0.0_dp)) / t_c
   end function growth_rate_two_fluid_inviscid

   !> Fastest-growing wavenumber, by dense sampling plus parabolic refinement.
   subroutine most_unstable_wavenumber(r0, sigma, rho, k_max, s_max)
      real(dp), intent(in)  :: r0, sigma, rho
      real(dp), intent(out) :: k_max, s_max
      integer,  parameter   :: n = 200001
      real(dp) :: x, k, s, dk, y0, y1, y2, denom, delta
      real(dp) :: kbest, sbest
      integer  :: i

      kbest = 0.0_dp; sbest = -1.0_dp
      dk = (1.0_dp / r0) / real(n - 1, dp)

      do i = 1, n
         x = real(i - 1, dp) / real(n - 1, dp)
         if (x <= 0.0_dp) x = 1.0e-9_dp
         k = x / r0
         s = growth_rate(k, r0, sigma, rho)
         if (s > sbest) then
            sbest = s
            kbest = k
         end if
      end do

      ! Parabolic vertex through the discrete maximum and its two neighbours, re-evaluated here.
      y0 = growth_rate(max(kbest - dk, 1.0e-12_dp), r0, sigma, rho)
      y1 = sbest
      y2 = growth_rate(min(kbest + dk, 1.0_dp / r0), r0, sigma, rho)
      denom = y0 - 2.0_dp * y1 + y2
      if (abs(denom) > 0.0_dp) then
         delta = 0.5_dp * (y0 - y2) / denom
      else
         delta = 0.0_dp
      end if

      k_max = kbest + delta * dk
      s_max = growth_rate(k_max, r0, sigma, rho)
   end subroutine most_unstable_wavenumber

   !> Wavelength of the fastest-growing mode.
   function most_unstable_wavelength(r0, sigma, rho) result(lam)
      real(dp), intent(in) :: r0, sigma, rho
      real(dp) :: lam, k_max, s_max
      call most_unstable_wavenumber(r0, sigma, rho, k_max, s_max)
      lam = 2.0_dp * PI / k_max
   end function most_unstable_wavelength

   !> Axial box length holding exactly n wavelengths of wavenumber k.
   pure function box_length_for_mode(k, n) result(L)
      real(dp), intent(in) :: k
      integer,  intent(in) :: n
      real(dp) :: L
      L = 2.0_dp * PI * real(n, dp) / k
   end function box_length_for_mode

   !> K_1(x)/K_0(x) for x >= 2 by Steed's continued fraction CF2 (Numerical Recipes 6.7, bessik, with mu = 0): K_1/K_0 = (x + 1/2 - h)/x.
   function k1_over_k0_cf2(x) result(rk)
      real(dp), intent(in) :: x
      real(dp) :: rk, k0s
      call k_cf2(x, rk, k0s)
   end function k1_over_k0_cf2

   !> Steed's continued fraction CF2 for x >= 2 (Numerical Recipes 6.7, bessik, mu = 0): returns K_1/K_0 and exp(x) K_0(x) = sqrt(pi/(2x))/S, where S is the series sum the same recurrence accumulates (rkmu in bessik without its factor exp(-x)). Split out of k1_over_k0_cf2 on 2026-10-22 with the same arithmetic, so the ratio is unchanged bit for bit.
   subroutine k_cf2(x, rk, k0s)
      real(dp), intent(in) :: x
      real(dp), intent(out) :: rk, k0s
      real(dp) :: a, a1, b, c, d, delh, dels, h, q, q1, q2, qnew, ss
      integer :: i
      b = 2.0_dp * (1.0_dp + x)
      d = 1.0_dp / b
      h = d
      delh = d
      q1 = 0.0_dp
      q2 = 1.0_dp
      a1 = 0.25_dp
      q = a1
      c = a1
      a = -a1
      ss = 1.0_dp + q * delh
      do i = 2, 10000
         a = a - real(2 * (i - 1), dp)
         c = -a * c / real(i, dp)
         qnew = (q1 - b * q2) / a
         q1 = q2
         q2 = qnew
         q = q + c * qnew
         b = b + 2.0_dp
         d = 1.0_dp / (b + a * d)
         delh = (b * d - 1.0_dp) * delh
         h = h + delh
         dels = q * delh
         ss = ss + dels
         if (abs(dels / ss) < epsilon(1.0_dp)) exit
      end do
      h = a1 * h
      rk = (x + 0.5_dp - h) / x
      k0s = sqrt(PI / (2.0_dp * x)) / ss
   end subroutine k_cf2

   !> Determinant of the linear system for the axisymmetric capillary instability of a viscous column (density rho_1, viscosity mu_1) in an unbounded viscous ambient (rho_2, mu_2), Tomotika (1935), in capillary units: R0 = rho_1 = sigma = 1, so time is in t_c = sqrt(rho_1 R0^3/sigma), oh1 = mu_1/sqrt(rho_1 sigma R0), rr = rho_2/rho_1, mr = mu_2/mu_1.
   !>
   !> Derivation (written out here because the code is the only record of it). Normal modes exp(s t + i x z). In each phase the velocity is the sum of a potential part, grad of phi with phi = A I_0(x r) inside and A K_0(x r) outside and pressure p = -rho s phi, and a vortical part satisfying rho s u = mu lap u with l^2 = x^2 + s/nu: inside u_r = C I_1(l r), u_z = i (l/x) C I_0(l r); outside u_r = C K_1(l r), u_z = -i (l/x) C K_0(l r) (both divergence free). The interface r = 1 + eta exp(s t + i x z) moves with s eta = u_r. The four conditions at r = 1 are continuity of u_r and u_z, continuity of the tangential stress mu (du_r/dz + du_z/dr), and the jump of the normal stress, (-p_1 + 2 mu_1 du_r1/dr) - (-p_2 + 2 mu_2 du_r2/dr) = -sigma kappa' with kappa' = (x^2 - 1) eta. Each column is divided by its own Bessel function (I_1(x), I_1(l_1), K_1(x), K_1(l_2)), which leaves only the ratios q = I_1/I_0 and p = K_1/K_0 and keeps the matrix finite at any viscosity. The unknowns are (A_1, C_1, A_2, C_2).
   !>
   !> With mr = 0 the ambient exerts no shear: the u_z continuity row and the outer vortical column are dropped and the 3 x 3 system remains; with rr = 0 as well it reduces to Chandrasekhar's relation. det is real for real s.
   function two_fluid_viscous_det(s, x, oh1, rr, mr) result(det)
      real(dp), intent(in) :: s, x, oh1, rr, mr
      real(dp) :: det
      real(dp) :: M(4,4), mu1, mu2, nu1, nu2, l1, l2, qk, q1, pk, p2
      integer :: n

      mu1 = oh1
      mu2 = mr * oh1
      nu1 = mu1
      qk = bessel_i1_over_i0(x)
      pk = bessel_k1_over_k0(x)
      l1 = sqrt(x * x + s / nu1)
      q1 = bessel_i1_over_i0(l1)
      M = 0.0_dp
      ! columns: A1, C1, A2, C2; rows: u_r, u_z, tangential stress, normal stress
      M(1,1) = x
      M(2,1) = x / qk
      M(3,1) = 2.0_dp * mu1 * x * x
      M(4,1) = s / qk + 2.0_dp * mu1 * x * x * (1.0_dp / qk - 1.0_dp / x) &
               + (x * x - 1.0_dp) * x / s
      M(1,2) = 1.0_dp
      M(2,2) = (l1 / x) / q1
      M(3,2) = mu1 * (x + l1 * l1 / x)
      M(4,2) = 2.0_dp * mu1 * l1 * (1.0_dp / q1 - 1.0_dp / l1) &
               + (x * x - 1.0_dp) / s
      M(1,3) = x
      M(2,3) = -x / pk
      M(3,3) = 2.0_dp * mu2 * x * x
      M(4,3) = -rr * s / pk - 2.0_dp * mu2 * x * x * (1.0_dp / pk + 1.0_dp / x)
      if (mr > 0.0_dp .and. rr > 0.0_dp) then
         nu2 = mu2 / rr
         l2 = sqrt(x * x + s / nu2)
         p2 = bessel_k1_over_k0(l2)
         M(1,4) = -1.0_dp
         M(2,4) = (l2 / x) / p2
         M(3,4) = -mu2 * (x + l2 * l2 / x)
         M(4,4) = 2.0_dp * mu2 * l2 * (1.0_dp / p2 + 1.0_dp / l2)
         n = 4
         det = det_small(M, n)
      else
         ! inviscid (or absent) ambient: rows u_r, tangential, normal; columns A1, C1, A2
         M(2,:) = M(3,:)
         M(3,:) = M(4,:)
         n = 3
         det = det_small(M(1:3,1:3), n)
      end if
   end function two_fluid_viscous_det

   !> Determinant of a small real matrix by Gaussian elimination with partial pivoting.
   function det_small(A, n) result(det)
      integer, intent(in) :: n
      real(dp), intent(in) :: A(:,:)
      real(dp) :: det, W(n,n), tmp(n), f
      integer :: i, k, p
      W = A(1:n,1:n)
      det = 1.0_dp
      do k = 1, n
         p = k - 1 + maxloc(abs(W(k:n,k)), 1)
         if (W(p,k) == 0.0_dp) then
            det = 0.0_dp; return
         end if
         if (p /= k) then
            tmp = W(k,:); W(k,:) = W(p,:); W(p,:) = tmp
            det = -det
         end if
         det = det * W(k,k)
         do i = k + 1, n
            f = W(i,k) / W(k,k)
            W(i,k:n) = W(i,k:n) - f * W(k,k:n)
         end do
      end do
   end function det_small

   !> Growth rate of the axisymmetric capillary instability of a viscous column in a viscous ambient, Tomotika (1935), without his Stokes approximation (inertia of both phases kept). Dimensional arguments as growth_rate_viscous_exact, plus the ambient density rho_outer and viscosity mu_outer; returns s in 1/time.
   !>
   !> Limits: mu_outer = rho_outer = 0 gives growth_rate_viscous_exact (Chandrasekhar); mu = mu_outer = 0 gives growth_rate_two_fluid_inviscid; for Oh >> 1 the most unstable wavenumber tends to Tomotika's Stokes values (x = 0.568 for mu/mu_outer = 0.91). The inviscid two-fluid rate bounds the root from above; the determinant is scanned downward from it on a logarithmic grid and the largest sign change is refined by bisection. Zero for x <= 0 or x >= 1; NaN if mu = 0 while mu_outer > 0 (not supported) or if no root is bracketed.
   function growth_rate_two_fluid_viscous(k, r0, sigma, rho, mu, rho_outer, &
                                          mu_outer) result(s)
      real(dp), intent(in) :: k, r0, sigma, rho, mu, rho_outer, mu_outer
      real(dp) :: s, x, t_c, oh1, rr, mr, upper, a, b, fa, fb, fm, mid
      integer, parameter :: NSCAN = 800
      integer :: i, it

      x = k * r0
      s = 0.0_dp
      if (x <= 0.0_dp .or. x >= 1.0_dp) return
      t_c = sqrt(rho * r0**3 / sigma)
      rr = rho_outer / rho
      if (mu <= 0.0_dp) then
         if (mu_outer > 0.0_dp) then
            s = nan_value()
         else
            s = growth_rate_two_fluid_inviscid(k, r0, sigma, rho, rho_outer)
         end if
         return
      end if
      oh1 = mu / sqrt(rho * sigma * r0)
      mr = mu_outer / mu
      upper = growth_rate_two_fluid_inviscid(k, r0, sigma, rho, rho_outer) * t_c
      upper = upper * (1.0_dp + 1.0e-9_dp)
      b = upper
      fb = two_fluid_viscous_det(b, x, oh1, rr, mr)
      s = nan_value()
      do i = 1, NSCAN
         a = upper * 10.0_dp**(-14.0_dp * real(i, dp) / real(NSCAN, dp))
         fa = two_fluid_viscous_det(a, x, oh1, rr, mr)
         if (fa * fb <= 0.0_dp) then
            do it = 1, 200
               mid = 0.5_dp * (a + b)
               fm = two_fluid_viscous_det(mid, x, oh1, rr, mr)
               if (fa * fm <= 0.0_dp) then
                  b = mid; fb = fm
               else
                  a = mid; fa = fm
               end if
               if (b - a <= 4.0_dp * epsilon(1.0_dp) * b) exit
            end do
            s = 0.5_dp * (a + b) / t_c
            return
         end if
         b = a
         fb = fa
      end do
   end function growth_rate_two_fluid_viscous

   !> exp(-x) I_n(x) for n = 0 or 1 and x >= 0, finite for any x. For x <= 30 from the power series mod_bessel_i (all terms positive, so accurate to a few ulp); for x > 30 from the asymptotic expansion of Abramowitz and Stegun 9.7.1, summed until the terms fall below a quarter of the unit roundoff (they decrease until about k = 2x, so this happens long before the expansion diverges). NaN for x < 0 or other n.
   function bessel_i_scaled(x, n) result(val)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: val, term, mu, z8
      integer :: k
      if (x < 0.0_dp .or. (n /= 0 .and. n /= 1)) then
         val = nan_value(); return
      end if
      if (x <= 30.0_dp) then
         val = mod_bessel_i(x, n) * exp(-x)
         return
      end if
      mu = 4.0_dp * real(n * n, dp)
      z8 = 8.0_dp * x
      term = 1.0_dp
      val = 1.0_dp
      do k = 1, 200
         term = -term * (mu - real((2 * k - 1)**2, dp)) / (real(k, dp) * z8)
         val = val + term
         if (abs(term) <= 0.25_dp * epsilon(1.0_dp) * abs(val)) exit
      end do
      val = val / sqrt(2.0_dp * PI * x)
   end function bessel_i_scaled

   !> exp(x) K_n(x) for n = 0 or 1 and x > 0, finite for any x: the series of k01_series for x <= 2, Steed's continued fraction (k_cf2) for x > 2. NaN for x <= 0 or other n.
   function bessel_k_scaled(x, n) result(val)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: val, k0, k1, rk, k0s
      if (x <= 0.0_dp .or. (n /= 0 .and. n /= 1)) then
         val = nan_value(); return
      end if
      if (x <= 2.0_dp) then
         call k01_series(x, k0, k1)
         if (n == 0) then
            val = k0 * exp(x)
         else
            val = k1 * exp(x)
         end if
      else
         call k_cf2(x, rk, k0s)
         if (n == 0) then
            val = k0s
         else
            val = rk * k0s
         end if
      end if
   end function bessel_k_scaled

   !> Inviscid growth rate of the column inside a coaxial wall at r = r_wall, with an ambient of density rho_outer between them (potential flow, u_r = 0 at the wall):
   !>
   !> s^2 t_c^2 = x (1 - x^2) / [ I_0(x)/I_1(x) + (rho_outer/rho) G ],  G = [K_0(x) I_1(x R) + I_0(x) K_1(x R)] / [K_1(x) I_1(x R) - I_1(x) K_1(x R)],
   !>
   !> with x = k r0 and R = r_wall/r0. G tends to K_0(x)/K_1(x) as R grows, recovering growth_rate_two_fluid_inviscid. Evaluated with the scaled Bessel functions, numerator and denominator both divided by exp(x R - x), so E = exp(-2 x (R - 1)) <= 1 is the only exponential left. It bounds the viscous confined root from above. Zero for x <= 0 or x >= 1; NaN if r_wall <= r0.
   function growth_rate_two_fluid_inviscid_confined(k, r0, sigma, rho, rho_outer, r_wall) &
         result(s)
      real(dp), intent(in) :: k, r0, sigma, rho, rho_outer, r_wall
      real(dp) :: s, x, rw, t_c, e, g, s2
      x = k * r0
      s = 0.0_dp
      if (r_wall <= r0) then
         s = nan_value(); return
      end if
      if (x <= 0.0_dp .or. x >= 1.0_dp) return
      rw = r_wall / r0
      t_c = sqrt(rho * r0**3 / sigma)
      e = exp(-2.0_dp * x * (rw - 1.0_dp))
      g = (bessel_k_scaled(x, 0) * bessel_i_scaled(x * rw, 1) &
           + bessel_i_scaled(x, 0) * bessel_k_scaled(x * rw, 1) * e) &
          / (bessel_k_scaled(x, 1) * bessel_i_scaled(x * rw, 1) &
           - bessel_i_scaled(x, 1) * bessel_k_scaled(x * rw, 1) * e)
      s2 = x * (1.0_dp - x * x) / (1.0_dp / bessel_i1_over_i0(x) + (rho_outer / rho) * g)
      s = sqrt(s2) / t_c
   end function growth_rate_two_fluid_inviscid_confined

   !> Determinant of the linear system for the axisymmetric capillary instability of a viscous column in a viscous ambient bounded by a coaxial no-slip wall at r = rw (rw = r_wall/R0 > 1), capillary units as in two_fluid_viscous_det.
   !>
   !> Derivation. As in two_fluid_viscous_det, but the ambient now occupies 1 < r < rw, so neither Bessel family can be dropped there: potential part phi = A_2 K_0(x r) + B_2 I_0(x r), vortical part u_r = C_2 K_1(l_2 r) + D_2 I_1(l_2 r), u_z = i (l_2/x) [-C_2 K_0(l_2 r) + D_2 I_0(l_2 r)]. Unknowns (A_1, C_1, A_2, B_2, C_2, D_2); rows: continuity of u_r, u_z, tangential stress and the normal-stress jump (with s eta = u_r) at r = 1, then u_r = 0 and u_z = 0 at r = rw. For each basis function the columns hold (u_r, u_z/i, d(u_z/i)/dr, du_r/dr, -p). Every basis function is multiplied by a constant that makes it of order one where it is largest (I-type functions of a r by exp(-a R_ref), K-type by exp(a), with R_ref = 1 in the liquid and rw in the ambient), so all entries stay finite at any viscosity; each column is then divided by its largest entry. Constant column factors do not move the roots. Same system as plotting/confined_dispersion.py, which was checked against two_fluid_viscous_det (unconfined limit) and against growth_rate_two_fluid_inviscid_confined (inviscid limit). Requires oh1, rr, mr > 0 and rw > 1; NaN otherwise.
   function two_fluid_viscous_confined_det(s, x, oh1, rr, mr, rw) result(det)
      real(dp), intent(in) :: s, x, oh1, rr, mr, rw
      real(dp) :: det, M(6,6), mu1, mu2, l1, l2, f(5), cmax
      integer :: j

      if (oh1 <= 0.0_dp .or. rr <= 0.0_dp .or. mr <= 0.0_dp .or. rw <= 1.0_dp) then
         det = nan_value(); return
      end if
      mu1 = oh1
      mu2 = mr * oh1
      l1 = sqrt(x * x + s / mu1)
      l2 = sqrt(x * x + s / (mu2 / rr))
      M = 0.0_dp
      do j = 1, 2
         if (j == 1) then
            f = pot_i(x, 1.0_dp, x, s)
         else
            f = vor_i(l1, 1.0_dp, l1, x)
         end if
         M(1,j) = f(1)
         M(2,j) = f(2)
         M(3,j) = mu1 * (f(3) + x * f(1))
         M(4,j) = f(5) + 2.0_dp * mu1 * f(4) - (1.0_dp - x * x) * f(1) / s
      end do
      do j = 3, 6
         f = gas_basis(j, 1.0_dp)
         M(1,j) = -f(1)
         M(2,j) = -f(2)
         M(3,j) = -mu2 * (f(3) + x * f(1))
         M(4,j) = -(f(5) + 2.0_dp * mu2 * f(4))
         f = gas_basis(j, rw)
         M(5,j) = f(1)
         M(6,j) = f(2)
      end do
      do j = 1, 6
         cmax = maxval(abs(M(:,j)))
         if (cmax > 0.0_dp) M(:,j) = M(:,j) / cmax
      end do
      det = det_small(M, 6)

   contains

      function gas_basis(col, r) result(g)
         integer, intent(in) :: col
         real(dp), intent(in) :: r
         real(dp) :: g(5)
         select case (col)
         case (3)
            g = pot_k(x, r, x, rr * s)
         case (4)
            g = pot_i(x, r, x * rw, rr * s)
         case (5)
            g = vor_k(l2, r, l2, x)
         case default
            g = vor_i(l2, r, l2 * rw, x)
         end select
      end function gas_basis

   end function two_fluid_viscous_confined_det

   !> I_n(a r) exp(-c) and K_n(a r) exp(c) from the scaled functions; used with c chosen so that a r - c <= 0 (I) and c - a r <= 0 (K), so the exponential never overflows.
   function i_sc(n, a, r, c) result(v)
      integer, intent(in) :: n
      real(dp), intent(in) :: a, r, c
      real(dp) :: v
      v = bessel_i_scaled(a * r, n) * exp(a * r - c)
   end function i_sc

   function k_sc(n, a, r, c) result(v)
      integer, intent(in) :: n
      real(dp), intent(in) :: a, r, c
      real(dp) :: v
      v = bessel_k_scaled(a * r, n) * exp(c - a * r)
   end function k_sc

   !> Basis functions of the confined system: (u_r, u_z/i, d(u_z/i)/dr, du_r/dr, -p) at r, scaled by exp(-c) (I-type) or exp(c) (K-type). Potential parts phi = I_0(x r) or K_0(x r) with -p = rho s phi (rhos = rho s); vortical parts with u_r = I_1(l r) or K_1(l r).
   function pot_i(x, r, c, rhos) result(f)
      real(dp), intent(in) :: x, r, c, rhos
      real(dp) :: f(5), i0, i1
      i0 = i_sc(0, x, r, c)
      i1 = i_sc(1, x, r, c)
      f = [x * i1, x * i0, x * x * i1, x * x * (i0 - i1 / (x * r)), rhos * i0]
   end function pot_i

   function pot_k(x, r, c, rhos) result(f)
      real(dp), intent(in) :: x, r, c, rhos
      real(dp) :: f(5), k0, k1
      k0 = k_sc(0, x, r, c)
      k1 = k_sc(1, x, r, c)
      f = [-x * k1, x * k0, -x * x * k1, x * x * (k0 + k1 / (x * r)), rhos * k0]
   end function pot_k

   function vor_i(l, r, c, x) result(f)
      real(dp), intent(in) :: l, r, c, x
      real(dp) :: f(5), i0, i1
      i0 = i_sc(0, l, r, c)
      i1 = i_sc(1, l, r, c)
      f = [i1, (l / x) * i0, (l * l / x) * i1, l * (i0 - i1 / (l * r)), 0.0_dp]
   end function vor_i

   function vor_k(l, r, c, x) result(f)
      real(dp), intent(in) :: l, r, c, x
      real(dp) :: f(5), k0, k1
      k0 = k_sc(0, l, r, c)
      k1 = k_sc(1, l, r, c)
      f = [k1, -(l / x) * k0, (l * l / x) * k1, -l * (k0 + k1 / (l * r)), 0.0_dp]
   end function vor_k

   !> Growth rate of the axisymmetric capillary instability of a viscous column (rho, mu) in a viscous ambient (rho_outer, mu_outer) inside a coaxial no-slip wall at r = r_wall: the problem the solver integrates (outer boundary NoSlipWall, hi_bc = 5, zero wall velocity). Dimensional arguments as growth_rate_two_fluid_viscous, plus r_wall; returns s in 1/time. The confined inviscid rate bounds the root from above; the root is bracketed and refined exactly as in growth_rate_two_fluid_viscous. Tends to growth_rate_two_fluid_viscous as r_wall grows. Zero for x <= 0 or x >= 1; NaN unless mu, rho_outer, mu_outer > 0 and r_wall > r0, or if no root is bracketed.
   function growth_rate_two_fluid_viscous_confined(k, r0, sigma, rho, mu, rho_outer, &
                                                   mu_outer, r_wall) result(s)
      real(dp), intent(in) :: k, r0, sigma, rho, mu, rho_outer, mu_outer, r_wall
      real(dp) :: s, x, t_c, oh1, rr, mr, rw, upper, a, b, fa, fb, fm, mid
      integer, parameter :: NSCAN = 800
      integer :: i, it

      x = k * r0
      s = 0.0_dp
      if (x <= 0.0_dp .or. x >= 1.0_dp) return
      if (mu <= 0.0_dp .or. rho_outer <= 0.0_dp .or. mu_outer <= 0.0_dp &
          .or. r_wall <= r0) then
         s = nan_value(); return
      end if
      t_c = sqrt(rho * r0**3 / sigma)
      rr = rho_outer / rho
      oh1 = mu / sqrt(rho * sigma * r0)
      mr = mu_outer / mu
      rw = r_wall / r0
      upper = growth_rate_two_fluid_inviscid_confined(k, r0, sigma, rho, rho_outer, r_wall)
      upper = upper * t_c * (1.0_dp + 1.0e-9_dp)
      b = upper
      fb = two_fluid_viscous_confined_det(b, x, oh1, rr, mr, rw)
      s = nan_value()
      do i = 1, NSCAN
         a = upper * 10.0_dp**(-14.0_dp * real(i, dp) / real(NSCAN, dp))
         fa = two_fluid_viscous_confined_det(a, x, oh1, rr, mr, rw)
         if (fa * fb <= 0.0_dp) then
            do it = 1, 200
               mid = 0.5_dp * (a + b)
               fm = two_fluid_viscous_confined_det(mid, x, oh1, rr, mr, rw)
               if (fa * fm <= 0.0_dp) then
                  b = mid; fb = fm
               else
                  a = mid; fa = fm
               end if
               if (b - a <= 4.0_dp * epsilon(1.0_dp) * b) exit
            end do
            s = 0.5_dp * (a + b) / t_c
            return
         end if
         b = a
         fb = fa
      end do
   end function growth_rate_two_fluid_viscous_confined

end module dispersion_mod
