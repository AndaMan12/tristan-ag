module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  use m_helpers
  implicit none

  !--- USER PARAMETERS (from <problem>) ------------------------!
  real, private :: B_0, B_amp, mode, spread, shift_gamma, TT, psi
  real, private :: mult1, mult2                 ! upstream / downstream multiplicities
  integer, private :: ramp_width                ! in units of Alfven wavelengths
  real, private :: profile_Lx                   ! "virtual" big-box length used for profiles

  !--- DERIVED / INTERNAL --------------------------------------!
  
  ! Integration step for density profiles (used consistently everywhere)
  real, private, parameter :: dx_int = 0.25d0   ! you can tune this (e.g. 0.1)

  real, private :: jA_ec_max                    ! max |jA/ec| over profile_Lx
  real, private :: weight_factor                ! macroparticle weight
  real, private :: max_density                  ! max(n_plus, n_minus) over profile_Lx
  integer, private :: init_x_boundary, fin_x_boundary
  logical, private :: norms_initialized = .false.

contains
  !=====================================================================
  !  INPUT
  !=====================================================================
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'B_0',          B_0)
    call getInput('problem', 'psi',          psi)
    call getInput('problem', 'B_amplitude',  B_amp)
    call getInput('problem', 'spread',       spread)
    call getInput('problem', 'mode',         mode)
    call getInput('problem', 'shift_gamma',  shift_gamma)
    call getInput('problem', 'temperature',  TT)
    call getInput('problem', 'multiplicity_1', mult1)
    call getInput('problem', 'multiplicity_2', mult2)
    call getInput('problem', 'ramp_width',   ramp_width)

    ! Virtual profile length (big-box length).
    call getInput('problem', 'profile_Lx', profile_Lx, 0.0)

    if (profile_Lx <= 0.0) then
      profile_Lx = real(global_mesh%sx)
    end if

    ! Positions where density transition happens (in physical x)
    init_x_boundary = int(5.0 * M_PI / mode) + 1
    fin_x_boundary  = init_x_boundary + int(ramp_width * 2.0 * M_PI / mode)
  end subroutine userReadInput

  !=====================================================================
  !  BASIC PULSE SHAPE
  !=====================================================================
  function b_wave(x_glob, y_glob, z_glob) result(bw)
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real :: bw

    if (x_glob <= 4.0 * M_PI / mode) then
      bw = (B_amp / 0.869619) * sin(x_glob * mode) * sin(x_glob * mode / 4.0)**2
    else
      bw = 0.0
    end if
  end function b_wave

  function db_wave_dx(x_glob, y_glob, z_glob) result(db)
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real :: db

    if (x_glob <= 4.0 * M_PI / mode) then
      db = - (B_amp / 0.869619) * 0.125 * &
           ( mode * ( cos(0.5*mode*x_glob) - 4.0*cos(mode*x_glob) + 3.0*cos(1.5*mode*x_glob) ) )
    else
      db = 0.0
    end if
  end function db_wave_dx

  !=====================================================================
  !  BISECTION FOR beta±
  !=====================================================================
  subroutine solve_for_beta(k, beta_plus, beta_minus)
    real, intent(in)  :: k
    real, intent(out) :: beta_plus, beta_minus
    real :: x, y

    call bisection_method(k, x, y, tolerance=1.0e-6, max_iterations=9000)
    beta_plus  = (1.0 - x**2) / (1.0 + x**2)
    beta_minus = (1.0 - y**2) / (1.0 + y**2)
  end subroutine solve_for_beta

  subroutine bisection_method(k, x, y, tolerance, max_iterations)
    real, intent(in)    :: k, tolerance
    integer, intent(in) :: max_iterations
    real, intent(out)   :: x, y
    real :: a, b, c, fa, fb, fc
    integer :: i

    a  = 1.0e-6
    b  = 2.0 - 1.0e-6
    fa = f(a, k)
    fb = f(b, k)

    if (fa * fb > 0.0) then
      print *, 'ERROR in bisection_method: function does not change sign'
      stop
    end if

    do i = 1, max_iterations
      c  = 0.5 * (a + b)
      fc = f(c, k)
      if (abs(fc) < tolerance) then
        x = c
        y = 2.0 - x
        return
      end if
      if (fa * fc < 0.0) then
        b  = c
        fb = fc
      else
        a  = c
        fa = fc
      end if
    end do

    print *, 'WARNING: maximum iterations exceeded in bisection_method'
    x = 0.5 * (a + b)
    y = 2.0 - x
  end subroutine bisection_method

  real function f(x, k)
    real, intent(in) :: x, k
    real :: y
    y = 2.0 - x
    f = (1.0 / x**2) - (1.0 / y**2) - 2.0 * k
  end function f

  !=====================================================================
  !  BACKGROUND DENSITY PROFILE n0(x)
  !=====================================================================
  subroutine bg_density_scaled(mult1, mult2, sx_prof, xg, n_0_local)
    real, intent(in)  :: mult1, mult2, sx_prof, xg
    real, intent(out) :: n_0_local
    real :: k

    if (mult1 == mult2) then
      n_0_local = mult1
    else
      ! smooth transition from mult1 → mult2 over [init_x_boundary, fin_x_boundary]
      k = log((mult1 - mult2) / (0.01 * mult2)) / &
          real(fin_x_boundary - init_x_boundary)

      if (xg <= real(init_x_boundary)) then
        n_0_local = mult1
      else
        n_0_local = mult1 + (mult2 - mult1) * &
                    (1.0 - (1.0 + k*(xg - real(init_x_boundary))) * &
                           exp(-k*(xg - real(init_x_boundary))))
      end if
    end if
  end subroutine bg_density_scaled

  !=====================================================================
  !  E_parallel PROFILE (same as old code)
  !=====================================================================
  function e_parr(x, jA_ec_max_in) result(ep)
    real, intent(in) :: x, jA_ec_max_in
    real :: ep
    real :: sx_prof
    real :: jA_ec, n_0_local, q1, q2, k, beta_plus1, beta_plus2, beta_minus

    sx_prof = profile_Lx

    call bg_density_scaled(mult1, mult2, sx_prof, x - 1.0e-1, n_0_local)
    n_0_local = n_0_local * jA_ec_max_in
    jA_ec     = -(db_wave_dx(x - 1.0e-1, 0.5, 0.5) * sin(psi) * B_norm) / unit_ch
    k         = jA_ec / n_0_local
    call solve_for_beta(k, beta_plus1, beta_minus)
    q1 = sqrt( (1.0 - beta_plus1) / (1.0 + beta_plus1) )

    call bg_density_scaled(mult1, mult2, sx_prof, x + 1.0e-1, n_0_local)
    n_0_local = n_0_local * jA_ec_max_in
    jA_ec     = -(db_wave_dx(x + 1.0e-1, 0.5, 0.5) * sin(psi) * B_norm) / unit_ch
    k         = jA_ec / n_0_local
    call solve_for_beta(k, beta_plus2, beta_minus)
    q2 = sqrt( (1.0 - beta_plus2) / (1.0 + beta_plus2) )

    ep = - cos(psi) * c_omp / sqrt(sigma) * (q2 - q1) / (2.0e-1)
  end function e_parr

  !=====================================================================
  !  NORMALIZATION (jA_ec_max, max_density, weight_factor)
  !=====================================================================
  subroutine init_normalization()
    implicit none
    real :: xg, sx_prof
    real :: n_0_local, jA_ec, k
    real :: beta_p, beta_m, n_plus, n_minus
    real :: pos_tot_dens, ele_tot_dens

    if (norms_initialized) return

    sx_prof       = profile_Lx
    jA_ec_max     = 0.0
    pos_tot_dens  = 0.0
    ele_tot_dens  = 0.0
    max_density   = 0.0

    ! 1) Find max |jA/ec| and integrate densities over [0, profile_Lx]
    do xg = 0.0, sx_prof, dx_int
      ! jA/e c
      jA_ec = -(db_wave_dx(xg, 0.5, 0.5) * sin(psi) * B_norm) / unit_ch
      jA_ec_max = max(jA_ec_max, abs(jA_ec))

      ! background multiplicity profile
      call bg_density_scaled(mult1, mult2, sx_prof, xg, n_0_local)
      n_0_local = n_0_local * jA_ec_max
      if (n_0_local <= 0.0) cycle

      k = jA_ec / n_0_local
      call solve_for_beta(k, beta_p, beta_m)

      n_plus  = n_0_local / (1.0 - beta_p)
      n_minus = n_0_local / (1.0 - beta_m)

      ! accumulate integrals ∫ n(x) dx with step dx_int
      pos_tot_dens = pos_tot_dens + n_plus  * dx_int
      ele_tot_dens = ele_tot_dens + n_minus * dx_int

      max_density = max(max_density, n_plus, n_minus)
    end do

    ! 2) Set particle weight so that average total ppc ≈ ppc0
    ! For a uniform pair plasma, this gives weight = 2 n0 / ppc0
    weight_factor = (pos_tot_dens + ele_tot_dens) / (ppc0 * sx_prof)

    if (mpi_rank == 0) then
      print *, 'init_normalization: jA_ec_max      = ', jA_ec_max
      print *, 'init_normalization: max_density    = ', max_density
      print *, 'init_normalization: weight_factor  = ', weight_factor
    end if

    norms_initialized = .true.
  end subroutine init_normalization


  !=====================================================================
  !  PARTICLE INJECTION IN A GIVEN X–SLAB [xmin, xmax]
  !=====================================================================
  subroutine inject_plasma_in_slab(xmin, xmax)
    implicit none
    real, intent(in) :: xmin, xmax
    real :: sx_slab, sy_glob
    real :: xg, yg
    real :: x_int_min, x_int_max, x
    real :: n_0_local, jA_ec, k
    real :: beta_p, beta_m, n_plus, n_minus
    real :: pos_tot_slab, ele_tot_slab
    integer :: npart_pos_target, npart_ele_target
    integer :: accepted
    real :: density_normalized, gamma, u_parr, v_perp, u_, v_, w_
    type(maxwellian) :: mwl

    sx_slab = max(xmax - xmin, 0.0)
    if (sx_slab <= 0.0) return

    sy_glob = real(global_mesh%sy)

    call init_normalization()

  !------------------ 1) integrate densities over the slab ----------------!
  pos_tot_slab = 0.0
  ele_tot_slab = 0.0

  x_int_min = max(0.0,        xmin)
  x_int_max = min(profile_Lx, xmax)

  if (x_int_max <= x_int_min) then
    return
  end if

  do x = x_int_min, x_int_max, dx_int
    call bg_density_scaled(mult1, mult2, profile_Lx, x, n_0_local)
    n_0_local = n_0_local * jA_ec_max
    if (n_0_local <= 0.0) cycle

    jA_ec = -(db_wave_dx(x, 0.5, 0.5) * sin(psi) * B_norm) / unit_ch
    k     = jA_ec / n_0_local
    call solve_for_beta(k, beta_p, beta_m)

    n_plus  = n_0_local / (1.0 - beta_p)
    n_minus = n_0_local / (1.0 - beta_m)

    ! Again, integrals ∫_slab n(x) dx with same dx_int as in init_normalization
    pos_tot_slab = pos_tot_slab + n_plus  * dx_int
    ele_tot_slab = ele_tot_slab + n_minus * dx_int
  end do

    ! If we are in a region where the profile is effectively zero, nothing to do.
    if (pos_tot_slab <= 0.0 .and. ele_tot_slab <= 0.0) return

    !------------------ 2) target macro counts in this slab -----------------!
    npart_pos_target = INT( sy_glob * pos_tot_slab / weight_factor )
    npart_ele_target = INT( sy_glob * ele_tot_slab / weight_factor )

    ! Maxwellian template
    mwl%dimension   = 2
    mwl%generated   = .false.
    mwl%temperature = TT
    mwl%shift_flag  = .true.

    !------------------ 3) positrons (species 1) ----------------------------!
    accepted = 0
    do while (accepted < npart_pos_target)
      xg = xmin + random(dseed) * sx_slab    ! absolute coordinate
      yg = random(dseed) * sy_glob

      call bg_density_scaled(mult1, mult2, profile_Lx, xg, n_0_local)
      n_0_local = n_0_local * jA_ec_max
      if (n_0_local <= 0.0) cycle

      jA_ec = -(db_wave_dx(xg, yg, 0.5) * sin(psi) * B_norm) / unit_ch
      k     = jA_ec / n_0_local
      call solve_for_beta(k, beta_p, beta_m)

      n_plus = n_0_local / (1.0 - beta_p)
      density_normalized = n_plus / max_density
      if (density_normalized <= 0.0) cycle
      if (density_normalized > 1.0) density_normalized = 1.0   ! safety

      if (random(dseed) <= density_normalized) then
        accepted = accepted + 1

        gamma = 1.0 / sqrt(1.0 - beta_p**2)
        if (abs(beta_p) < 1.0e-3) then
          mwl%shift_gamma = 1.0
          mwl%shift_dir   = 1
        else
          mwl%shift_gamma = gamma
          mwl%shift_dir   = INT(sign(1.0, beta_p))
        end if

        call generateFromMaxwellian(mwl, u_parr, v_perp, w_)
        u_ = u_parr * cos(psi) - v_perp * sin(psi)
        v_ = u_parr * sin(psi) + v_perp * cos(psi)

        call injectParticleGlobally(1, xg, yg, 0.5, u_, v_, w_, weight_factor)
      end if
    end do

    !------------------ 4) electrons (species 2) ----------------------------!
    accepted = 0
    do while (accepted < npart_ele_target)
      xg = xmin + random(dseed) * sx_slab
      yg = random(dseed) * sy_glob

      call bg_density_scaled(mult1, mult2, profile_Lx, xg, n_0_local)
      n_0_local = n_0_local * jA_ec_max
      if (n_0_local <= 0.0) cycle

      jA_ec = -(db_wave_dx(xg, yg, 0.5) * sin(psi) * B_norm) / unit_ch
      k     = jA_ec / n_0_local
      call solve_for_beta(k, beta_p, beta_m)

      n_minus = n_0_local / (1.0 - beta_m)
      density_normalized = n_minus / max_density
      if (density_normalized <= 0.0) cycle
      if (density_normalized > 1.0) density_normalized = 1.0

      if (random(dseed) <= density_normalized) then
        accepted = accepted + 1

        gamma = 1.0 / sqrt(1.0 - beta_m**2)
        if (abs(beta_m) < 1.0e-3) then
          mwl%shift_gamma = 1.0
          mwl%shift_dir   = 1
        else
          mwl%shift_gamma = gamma
          mwl%shift_dir   = INT(sign(1.0, beta_m))
        end if

        call generateFromMaxwellian(mwl, u_parr, v_perp, w_)
        u_ = u_parr * cos(psi) - v_perp * sin(psi)
        v_ = u_parr * sin(psi) + v_perp * cos(psi)

        call injectParticleGlobally(2, xg, yg, 0.5, u_, v_, w_, weight_factor)
      end if
    end do

  end subroutine inject_plasma_in_slab

  !=====================================================================
  !  FIELD INITIALISATION IN A SLAB [xmin, xmax]
  !=====================================================================
  subroutine init_fields_in_slab(xmin, xmax)
    implicit none
    real, intent(in) :: xmin, xmax
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    integer :: i_start, i_end
    real :: ex_parr, ey_parr

    ! Only overwrite fields in [xmin, xmax] on this meshblock.
    i_start = max(0, int(ceiling(xmin - real(this_meshblock%ptr%x0))))
    i_end   = min(this_meshblock%ptr%sx - 1, int(floor(xmax - real(this_meshblock%ptr%x0))))
    if (i_end < i_start) return

    do i = i_start, i_end
      i_glob  = i + this_meshblock%ptr%x0
      ex_parr = e_parr(real(i_glob),       jA_ec_max)
      ey_parr = e_parr(real(i_glob) - 0.5, jA_ec_max)

      do j = 0, this_meshblock%ptr%sy - 1
        j_glob = j + this_meshblock%ptr%y0
        do k = 0, this_meshblock%ptr%sz - 1
          k_glob = k + this_meshblock%ptr%z0

          ! Background B only in this slab
          bx(i, j, k) = B_0 * cos(psi)
          by(i, j, k) = B_0 * sin(psi)

          ex(i, j, k) = -sin(psi) * b_wave(real(i_glob),       real(j_glob-0.5), real(k_glob-0.5)) + &
                         cos(psi) * ex_parr
          ey(i, j, k) =  cos(psi) * b_wave(real(i_glob-0.5),   real(j_glob),     real(k_glob-0.5)) + &
                         sin(psi) * ey_parr
          bz(i, j, k) =  b_wave(real(i_glob),                 real(j_glob),     real(k_glob-0.5))
        end do
      end do
    end do
  end subroutine init_fields_in_slab

  !=====================================================================
  !  INITIALISATION AT t = 0
  !=====================================================================
  subroutine userInitParticles()
    implicit none
    real :: xmin, xmax

    call init_normalization()

    xmin = real(global_mesh%x0)
    xmax = xmin + real(global_mesh%sx)

    call inject_plasma_in_slab(xmin, xmax)
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    real :: xmin, xmax

    call init_normalization()

    ! Zero E initially, fill background B over the whole window
    ex(:, :, :) = 0.0
    ey(:, :, :) = 0.0
    ez(:, :, :) = 0.0
    bx(:, :, :) = B_0 * cos(psi)
    by(:, :, :) = B_0 * sin(psi)
    bz(:, :, :) = 0.0

    xmin = real(global_mesh%x0)
    xmax = xmin + real(global_mesh%sx)

    call init_fields_in_slab(xmin, xmax)
  end subroutine userInitFields

  !=====================================================================
  !  CALLED BY MOVING WINDOW: REFILL NEW RIGHT-HAND REGION
  !=====================================================================
  subroutine userFillNewRegion(xmin, xmax)
    implicit none
    real, intent(in) :: xmin, xmax

    call init_normalization()
    call init_fields_in_slab(xmin, xmax)
    call inject_plasma_in_slab(xmin, xmax)
  end subroutine userFillNewRegion

  !=====================================================================
  !  EXTRA HOOKS
  !=====================================================================
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! No custom current deposition for the moving-window run.
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userDriveParticles

  !--- Open particle boundary on the LEFT side of the moving window ----!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    integer :: s, ti, tj, tk, p
    real :: x_l, x_g
    real :: x_left, x_right

    x_left  = real(global_mesh%x0)
    x_right = x_left + real(global_mesh%sx)

    do s = 1, nspec
      do tk = 1, species(s)%tile_nz
        do tj = 1, species(s)%tile_ny
          do ti = 1, species(s)%tile_nx
            p = 1
            do while (p <= species(s)%prtl_tile(ti, tj, tk)%npart_sp)
              x_l = real(species(s)%prtl_tile(ti, tj, tk)%xi(p)) + &
                    species(s)%prtl_tile(ti, tj, tk)%dx(p)
              x_g = x_l + real(this_meshblock%ptr%x0)

              ! Open boundary at left; right side is fed by userFillNewRegion.
              if (x_g < x_left) t
