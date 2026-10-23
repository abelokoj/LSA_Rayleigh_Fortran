!> Tests for dispersion_mod.
!>
!> The module reproduces the truncated power series used by the solver's own run2d/rayleigh_capillary_growth.F90, so testing it against that same series would be circular.  Instead I_n(x) is cross-checked against an INDEPENDENT evaluation (numerical integration of the standard integral representation
!>
!> I_n(x) = (1/pi) * integral_0^pi exp(x cos t) cos(n t) dt
!>
!> by the trapezoid rule, which converges spectrally for this smooth periodic integrand).  Agreement to ~1e-14 pins the series independently.
!>
!> The headline value, sigma = 0.34332680120934783 at the single-precision k = 0.7 the original program actually uses, is checked against the number printed by that program when compiled and run.
program test_dispersion
   use testing_mod
   use dispersion_mod
   use util_mod, only: read_real, read_int, is_nan
   implicit none

   real(dp) :: k, k_max, s_max, lam, s1, s2, s4
   integer :: i, ierr, i_val

   call section('modified Bessel I_n against an independent integral')
   do i = 0, 10
      k = real(i, dp) / 10.0_dp
      ! Mixed tolerance: I_1(0) is exactly zero, so a purely relative test against the quadrature's round-off can never pass there.
      call check_close_mixed(mod_bessel_i(k, 0), bessel_i_quad(k, 0), &
                             1.0e-14_dp, 1.0e-13_dp, 'I_0 matches quadrature')
      call check_close_mixed(mod_bessel_i(k, 1), bessel_i_quad(k, 1), &
                             1.0e-14_dp, 1.0e-13_dp, 'I_1 matches quadrature')
   end do

   call section('growth rate against the original Fortran utility')
   ! rayleigh_capillary_growth.F90 stores k in SINGLE precision, so the value it actually evaluates is 0.69999998807907104, and it prints 0.34332680120934783.  Feeding the same value must reproduce it exactly.
   call check_close(growth_rate(0.69999998807907104_dp, 1.0_dp, 1.0_dp, &
      1.0_dp), &
                    0.34332680120934783_dp, 1.0e-15_dp, &
                    'reproduces the shipped utility to 1e-15')

   call section('stability boundary')
   call check(growth_rate_squared(0.999_dp, 1.0_dp, 1.0_dp, 1.0_dp) > 0.0_dp, &
              'unstable just below k r0 = 1')
   call check_close(growth_rate_squared(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp), &
                    0.0_dp, 1.0e-12_dp, 'marginal at k r0 = 1')
   call check(growth_rate_squared(1.001_dp, 1.0_dp, 1.0_dp, 1.0_dp) < 0.0_dp, &
              'stable just above k r0 = 1')
   call check_close(growth_rate(1.5_dp, 1.0_dp, 1.0_dp, 1.0_dp), 0.0_dp, &
                    0.0_dp, 'growth rate is zero (not NaN) in the stable band')

   call section('most unstable mode')
   call most_unstable_wavenumber(1.0_dp, 1.0_dp, 1.0_dp, k_max, s_max)
   ! The .F90 header states "anecdotedly, k r0 = .7 is the critical value".
   call check_close(k_max, 0.697_dp, 5.0e-3_dp, 'peak at k r0 ~= 0.697')
   lam = most_unstable_wavelength(1.0_dp, 1.0_dp, 1.0_dp)
   ! inputs.growthrate.LSA uses yblob = 9.0 for "wave length = 2 pi r0/.7"
   call check_close(lam, 9.0_dp, 0.1_dp, 'lambda ~= 9.0 r0')
   call check_close_rel(lam, 2.0_dp*PI/0.7_dp, 0.01_dp, &
                        'lambda ~= 2 pi r0 / 0.7')
   call check(s_max >= growth_rate(0.5_dp,1.0_dp,1.0_dp,1.0_dp), &
              'peak growth exceeds k r0 = 0.5')
   call check(s_max >= growth_rate(0.9_dp,1.0_dp,1.0_dp,1.0_dp), &
              'peak growth exceeds k r0 = 0.9')

   call section('scaling with physical parameters')
   s1 = growth_rate(0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp)
   s2 = growth_rate(0.7_dp, 1.0_dp, 4.0_dp, 1.0_dp)
   call check_close_rel(s2, 2.0_dp*s1, 1.0e-12_dp, 'sigma ~ sqrt(surface tension)')
   s2 = growth_rate(0.7_dp, 1.0_dp, 1.0_dp, 4.0_dp)
   call check_close_rel(s2, 0.5_dp*s1, 1.0e-12_dp, 'sigma ~ 1/sqrt(density)')
   ! at fixed k*r0, doubling r0 scales sigma by r0^{-3/2}
   s4 = growth_rate(0.35_dp, 2.0_dp, 1.0_dp, 1.0_dp)
   call check_close_rel(s4, s1*2.0_dp**(-1.5_dp), 1.0e-12_dp, &
                        'sigma ~ r0^(-3/2) at fixed k r0')

   call section('viscous correction')
   call check_close(growth_rate_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp), &
                    s1, 1.0e-14_dp, 'reduces to inviscid as mu -> 0')
   call check(growth_rate_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp) < s1, &
              'viscosity reduces the growth rate')
   call check(growth_rate_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp) > &
      0.0_dp, &
              'still unstable at Oh = 0.02')
   ! stable band: must be real and non-positive, never NaN
   s2 = growth_rate_viscous(1.5_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp)
   call check(s2 == s2, 'viscous rate is not NaN in the stable band')
   call check(s2 <= 0.0_dp, 'viscous rate is non-positive in the stable band')

   call section('exact viscous relation (Chandrasekhar)')
   ! Reference values from scipy (iv, brentq) through plotting/dispersion.py growth_rate_viscous_exact, an independent implementation of the same relation; agreement to 1e-9 shows the two codes solve the same equation.
   call check_close_rel(growth_rate_viscous_exact(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp), &
      0.3293671031521454_dp, 1.0e-9_dp, 'x=0.7, Oh=0.02 matches the scipy solution')
   call check_close_rel(growth_rate_viscous_exact(0.5_dp,1.0_dp,1.0_dp,1.0_dp,0.1_dp), &
      0.26723370477999997_dp, 1.0e-9_dp, 'x=0.5, Oh=0.1 matches the scipy solution')
   call check_close_rel(growth_rate_viscous_exact(0.3_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp), &
      0.10798305831705056_dp, 1.0e-9_dp, 'x=0.3, Oh=1 matches the scipy solution')
   call check_close_rel(growth_rate_viscous_exact(0.9_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp), &
      0.24272977221143652_dp, 1.0e-9_dp, 'x=0.9, Oh=0.02 matches the scipy solution')
   call check_close_rel(growth_rate_viscous_exact(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp), &
      s1, 1.0e-14_dp, 'exact viscous reduces to Rayleigh at Oh = 0')
   call check_close_rel(growth_rate_viscous_exact(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.0e-6_dp), &
      s1, 1.0e-4_dp, 'and tends to Rayleigh as Oh -> 0')
   ! Long-wave limit (dissertation eq:3_longwave): s = sqrt(x^2/2 - x^4/2 + 9/4 Oh^2 x^4) - 3/2 Oh x^2, independent of Bessel functions.
   s2 = 0.05_dp
   s4 = sqrt(0.5_dp*s2**2*(1.0_dp - s2**2) + 2.25_dp*s2**4) - 1.5_dp*s2**2
   call check_close_rel(growth_rate_viscous_exact(s2,1.0_dp,1.0_dp,1.0_dp,1.0_dp), &
      s4, 1.0e-3_dp, 'agrees with the long-wave formula at x = 0.05, Oh = 1')
   ! High Ohnesorge number, where a fixed lower bracket of 1e-14 * upper used to give NaN.
   call check_close_rel(growth_rate_viscous_exact(0.98_dp,1.0_dp,1.0_dp,1.0_dp,10.0_dp), &
      0.0006562396926536359_dp, 1.0e-8_dp, 'x=0.98, Oh=10 is solved (not NaN)')
   call check_close_rel(growth_rate_viscous_exact(0.5_dp,1.0_dp,1.0_dp,1.0_dp,100.0_dp), &
      0.0012494527504896898_dp, 1.0e-8_dp, 'x=0.5, Oh=100 is solved (not NaN)')
   call check_close(growth_rate_viscous_exact(1.2_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp), &
      0.0_dp, 0.0_dp, 'zero in the stable band')
   ! Scaling: dimensional rate = (dimensionless rate)/t_c with t_c = sqrt(rho r0^3/sigma).
   call check_close_rel(growth_rate_viscous_exact(0.35_dp,2.0_dp,1.0_dp,1.0_dp,0.02_dp*sqrt(2.0_dp)), &
      0.3293671031521454_dp*2.0_dp**(-1.5_dp), 1.0e-9_dp, &
      'scales with the capillary time at fixed x and Oh')
   call check(abs(growth_rate_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp) - &
      growth_rate_viscous_exact(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp)) < &
      0.01_dp*growth_rate_viscous_exact(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp), &
      'the Weber approximation is within 1% of exact at Oh = 0.02')

   call section('Bessel ratios')
   call check_close_rel(bessel_i1_over_i0(0.01_dp), 0.004999937501041649_dp, &
      1.0e-13_dp, 'I1/I0 at 0.01')
   call check_close_rel(bessel_i1_over_i0(1.0_dp), 0.4463899658965345_dp, &
      1.0e-13_dp, 'I1/I0 at 1')
   call check_close_rel(bessel_i1_over_i0(7.3_dp), 0.9287379172265551_dp, &
      1.0e-13_dp, 'I1/I0 at 7.3')
   call check_close_rel(bessel_i1_over_i0(50.0_dp), 0.9899489673784977_dp, &
      1.0e-13_dp, 'I1/I0 at 50 (beyond the power series)')
   call check_close_rel(bessel_i1_over_i0(180.0_dp), 0.9972183425758241_dp, &
      1.0e-13_dp, 'I1/I0 at 180')
   call check_close_rel(bessel_i1_over_i0(0.7_dp), &
      mod_bessel_i(0.7_dp,1)/mod_bessel_i(0.7_dp,0), 1.0e-14_dp, &
      'ratio agrees with the power series at 0.7')
   call check_close_rel(bessel_k1_over_k0(0.3_dp), 2.22665279763796_dp, &
      1.0e-12_dp, 'K1/K0 at 0.3')
   call check_close_rel(bessel_k1_over_k0(0.7_dp), 1.5900862321504068_dp, &
      1.0e-12_dp, 'K1/K0 at 0.7')
   call check_close_rel(bessel_k1_over_k0(0.95_dp), 1.4498495858097538_dp, &
      1.0e-12_dp, 'K1/K0 at 0.95')
   ! The series now stops early; the result must be bit for bit that of the full 100-term sum, which is what makes it identical to the solver's rayleigh_capillary_growth.F90.
   call check(mod_bessel_i(0.7_dp,0) == full_series(0.7_dp,0) .and. &
              mod_bessel_i(0.7_dp,1) == full_series(0.7_dp,1) .and. &
              mod_bessel_i(0.3_dp,1) == full_series(0.3_dp,1) .and. &
              mod_bessel_i(5.0_dp,0) == full_series(5.0_dp,0), &
              'early-stopped series is bit-identical to the 100-term sum')

   call section('two-fluid inviscid relation')
   call check_close_rel(growth_rate_two_fluid_inviscid(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp), &
      s1, 1.0e-13_dp, 'reduces to Rayleigh with no ambient density')
   call check_close_rel(growth_rate_two_fluid_inviscid(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.225e-3_dp), &
      0.34328314375322494_dp, 1.0e-10_dp, 'x=0.7, rho_g=1.225e-3 matches scipy')
   call check_close_rel(growth_rate_two_fluid_inviscid(0.3_dp,1.0_dp,1.0_dp,1.0_dp,1.225e-3_dp), &
      0.20122826095072902_dp, 1.0e-10_dp, 'x=0.3, rho_g=1.225e-3 matches scipy')
   call check(growth_rate_two_fluid_inviscid(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.225e-3_dp) < s1, &
      'the ambient gas lowers the rate')

   call section('two-fluid viscous relation (Tomotika)')
   call check_close_rel(growth_rate_two_fluid_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,0.0_dp,0.0_dp), &
      growth_rate_viscous_exact(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp), 1.0e-12_dp, &
      'no ambient: equals Chandrasekhar at x=0.7')
   call check_close_rel(growth_rate_two_fluid_viscous(0.3_dp,1.0_dp,1.0_dp,1.0_dp,0.2_dp,0.0_dp,0.0_dp), &
      growth_rate_viscous_exact(0.3_dp,1.0_dp,1.0_dp,1.0_dp,0.2_dp), 1.0e-12_dp, &
      'no ambient: equals Chandrasekhar at x=0.3, Oh=0.2')
   call check_close_rel(growth_rate_two_fluid_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.0e-10_dp,0.5_dp,1.0e-10_dp), &
      growth_rate_two_fluid_inviscid(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.5_dp), 1.0e-4_dp, &
      'vanishing viscosities: tends to the two-fluid inviscid rate')
   call check(growth_rate_two_fluid_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp) < &
      growth_rate_viscous_exact(0.7_dp,1.0_dp,1.0_dp,1.0_dp, &
            0.02_dp), 'the ambient gas lowers the viscous rate')
   call check(growth_rate_two_fluid_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-5_dp) > &
      growth_rate_two_fluid_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp), &
      'a less viscous ambient raises it')
   call check(growth_rate_two_fluid_viscous(1.1_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp) == 0.0_dp, &
      'zero in the stable band')
   call check(abs(stokes_xmax(0.91_dp) - 0.568_dp) < 0.004_dp, &
      'Stokes limit, mu/mu_outer = 0.91: fastest x = 0.568 (Tomotika 1935, eq. 41)')
   call check(abs(stokes_xmax(1.0_dp) - 2.0_dp*PI/11.16_dp) < 0.003_dp, &
      'Stokes limit, equal viscosities: fastest wavelength 11.16 R0')
   call check_close_rel(bessel_k1_over_k0(5.0_dp), 1.0957750456413309_dp, 1.0e-12_dp, &
      'K1/K0 at x = 5 (continued fraction) matches scipy')
   call check_close_rel(bessel_k1_over_k0(50.0_dp), 1.0099509699867857_dp, 1.0e-12_dp, &
      'K1/K0 at x = 50 matches scipy')
   call check_close_rel(bessel_k1_over_k0(2.0_dp), 1.2280369298189084_dp, 1.0e-12_dp, &
      'K1/K0 at x = 2 (series) matches scipy')

   call section('scaled Bessel functions (used by the confined relation)')
   ! Reference values from scipy.special.ive / kve (exponentially scaled I_n and K_n), independent of the series, continued fraction and asymptotic expansion used here.
   call check_close_rel(bessel_i_scaled(0.5_dp, 0), 0.6450352704491501_dp, 1.0e-13_dp, &
         'exp(-x) I_0 at x = 0.5 (series)')
   call check_close_rel(bessel_i_scaled(5.0_dp, 1), 0.16397226694454237_dp, 1.0e-13_dp, &
         'exp(-x) I_1 at x = 5 (series)')
   call check_close_rel(bessel_i_scaled(30.0_dp, 0), 0.07314594648223727_dp, 1.0e-13_dp, &
         'exp(-x) I_0 at x = 30 (series, last point)')
   call check_close_rel(bessel_i_scaled(35.0_dp, 1), 0.06670443172949145_dp, 1.0e-13_dp, &
         'exp(-x) I_1 at x = 35 (asymptotic)')
   call check_close_rel(bessel_i_scaled(2000.0_dp, 0), 0.008921178276439672_dp, 1.0e-13_dp, &
         'exp(-x) I_0 at x = 2000 (asymptotic)')
   call check_close_rel(bessel_k_scaled(0.1_dp, 1), 10.890182683049698_dp, 1.0e-13_dp, &
         'exp(x) K_1 at x = 0.1 (series)')
   call check_close_rel(bessel_k_scaled(2.0_dp, 0), 0.8415682150707713_dp, 1.0e-13_dp, &
         'exp(x) K_0 at x = 2 (series, last point)')
   call check_close_rel(bessel_k_scaled(2.5_dp, 0), 0.7595486903280996_dp, 1.0e-13_dp, &
         'exp(x) K_0 at x = 2.5 (continued fraction)')
   call check_close_rel(bessel_k_scaled(100.0_dp, 1), 0.12579995047957854_dp, 1.0e-13_dp, &
         'exp(x) K_1 at x = 100 (continued fraction)')
   call check_close_rel(bessel_k_scaled(2000.0_dp, 0), 0.02802320501460432_dp, 1.0e-13_dp, &
         'exp(x) K_0 at x = 2000 (continued fraction)')
   call check_close_rel(bessel_k_scaled(5.0_dp, 1) / bessel_k_scaled(5.0_dp, 0), bessel_k1_over_k0(5.0_dp), &
      1.0e-15_dp, 'scaled K_1/K_0 agrees with bessel_k1_over_k0')

   call section('two-fluid viscous relation with a no-slip outer wall')
   ! Campaign fluid properties: Oh = 0.02, rho_2/rho_1 = 1.225e-3, mu_2 = 2.6e-4. Reference values from the independent Python implementation plotting/confined_dispersion.py (scipy Bessel functions), which reproduces the unconfined Tomotika values of this module to 5e-13.
   call check_close_rel(growth_rate_two_fluid_viscous_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp,1.5_dp), &
      0.31830531087778385_dp, 1.0e-10_dp, 'x=0.7, wall at 1.5 R0 matches the Python implementation')
   call check_close_rel(growth_rate_two_fluid_viscous_confined(0.5_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp,1.5_dp), &
      0.28274258047133055_dp, 1.0e-10_dp, 'x=0.5, wall at 1.5 R0 matches the Python implementation')
   call check_close_rel(growth_rate_two_fluid_viscous_confined(0.8_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp,2.0_dp), &
      0.3064469154131196_dp, 1.0e-10_dp, 'x=0.8, wall at 2 R0 matches the Python implementation')
   call check_close_rel(growth_rate_two_fluid_viscous_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,0.02_dp,4.0_dp), &
      0.2984047482442836_dp, 1.0e-10_dp, 'x=0.7, equal viscosities, wall at 4 R0 matches the Python implementation')
   call check_close_rel(growth_rate_two_fluid_viscous_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp,40.0_dp), &
      growth_rate_two_fluid_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp), 1.0e-10_dp, &
      'wall at 40 R0: equals the unconfined Tomotika rate')
   call check(growth_rate_two_fluid_viscous_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp,1.5_dp) < &
      growth_rate_two_fluid_viscous_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp, &
            2.0_dp), &
      'a closer wall lowers the rate')
   call check_close_rel(growth_rate_two_fluid_inviscid_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.5_dp), &
      0.2564549522576867_dp, 1.0e-12_dp, 'inviscid confined, dense ambient, wall at 1.5 R0 matches scipy')
   call check_close_rel(growth_rate_two_fluid_inviscid_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,40.0_dp), &
      growth_rate_two_fluid_inviscid(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp), 1.0e-12_dp, &
      'inviscid confined, wall at 40 R0: equals the unconfined two-fluid inviscid rate')
   ! Independent check of the wall rows: with a dense ambient the wall lowers the inviscid rate by 18 %, and the viscous confined rate must tend to the confined (not the unconfined) inviscid rate as both viscosities vanish.
   call check_close_rel(growth_rate_two_fluid_viscous_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.0e-10_dp,1.0_dp,0.5e-10_dp,1.5_dp), &
      growth_rate_two_fluid_inviscid_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.5_dp), 1.0e-4_dp, &
      'vanishing viscosities: tends to the confined inviscid rate')
   call check_close_rel(growth_rate_two_fluid_viscous_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,1.0e-10_dp,1.0_dp,0.5e-10_dp,1.5_dp), &
      0.25644922904899836_dp, 1.0e-9_dp, 'vanishing viscosities: matches the Python implementation')
   call check(growth_rate_two_fluid_viscous_confined(1.1_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp,1.5_dp) == 0.0_dp, &
      'zero in the stable band')
   call check(is_nan(growth_rate_two_fluid_viscous_confined(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp,1.225e-3_dp,2.6e-4_dp,0.9_dp)), &
      'a wall inside the column is refused (NaN)')

   call section('command-line number parsing')
   call read_real('0.7,0.9', s2, ierr)
   call check(ierr /= 0, 'a list given where one number is expected is refused')
   call read_real(' 2.5e-3 ', s2, ierr)
   call check(ierr == 0 .and. s2 == 2.5e-3_dp, 'a single number is read')
   call read_int('1.6', i_val, ierr)
   call check(ierr /= 0, 'a non-integer is refused where an integer is expected')

   call section('Ohnesorge number and box length')
   call check_close(ohnesorge(0.02_dp,1.0_dp,1.0_dp,1.0_dp), 0.02_dp, &
                    1.0e-14_dp, 'Oh = mu for unit rho, sigma, r0')
   call check_close(box_length_for_mode(0.7_dp, 1), 2.0_dp*PI/0.7_dp, &
                    1.0e-12_dp, 'box holds one wavelength')
   call check_close(box_length_for_mode(0.7_dp, 8), 8.0_dp*2.0_dp*PI/0.7_dp, &
                    1.0e-12_dp, 'box holds eight wavelengths')

   call summary('test_dispersion')

contains

   !> Most unstable wavenumber of the two-fluid viscous relation at Oh = 100 with equal densities and viscosity ratio lam = mu/mu_outer (the Stokes regime of Tomotika 1935).
   function stokes_xmax(lam) result(xb)
      real(dp), intent(in) :: lam
      real(dp) :: xb, xx, s, best
      best = -1.0_dp
      xb = 0.0_dp
      xx = 0.50_dp
      do while (xx < 0.62_dp)
         s = growth_rate_two_fluid_viscous(xx, 1.0_dp, 1.0_dp, 1.0_dp, 100.0_dp, &
                                           1.0_dp, 100.0_dp / lam)
         if (s > best) then
            best = s; xb = xx
         end if
         xx = xx + 0.0005_dp
      end do
   end function stokes_xmax

   !> I_n(x) by trapezoid quadrature of the integral representation. Independent of the power series under test.
   function bessel_i_quad(x, n) result(val)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: val, t, h
      integer :: i, m
      m = 20000
      h = PI / real(m, dp)
      val = 0.0_dp
      do i = 0, m
         t = real(i, dp) * h
         if (i == 0 .or. i == m) then
            val = val + 0.5_dp * exp(x*cos(t)) * cos(real(n,dp)*t)
         else
            val = val + exp(x*cos(t)) * cos(real(n,dp)*t)
         end if
      end do
      val = val * h / PI
   end function bessel_i_quad

   !> The original 100-term power series, kept here as the reference for the early-stopping check. Terms past the point of overflow are skipped, since they are zero in the original too.
   function full_series(x, n) result(val)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: val, m_fact, mpn_fact, xpower, d
      integer :: i
      m_fact = 1.0_dp; mpn_fact = 1.0_dp
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
         if (m_fact > huge(1.0_dp) / mpn_fact) exit   ! product would overflow
         d = m_fact * mpn_fact
         val = val + xpower / d
      end do
   end function full_series

end program test_dispersion
