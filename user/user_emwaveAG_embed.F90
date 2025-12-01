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
! #ifdef USROUTPUT
!   use m_writeusroutput
! #endif
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: B_0, B_amp, mode, spread, shift_gamma, TT, psi, mult1, mult2
  logical, private :: use_moving_window = .false.
  integer, private :: mw_shift_start = 0, mw_shift_interval = 0
  real(kind=dprec), private :: mw_speed = 0.0, mw_gamma_param = 0.0
  real, private :: cached_weight_factor = 1.0, cached_max_density = 0.0
  real, private :: cached_pos_tot_dens = 0.0, cached_ele_tot_dens = 0.0
  real, private :: cached_jA_ec_max = 0.0
  logical, private :: profile_cache_ready = .false.
  integer, private :: init_x_boundary, fin_x_boundary, ramp_width ! ramp width in units of Alfven wavelengths
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  ! private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
    subroutine userReadInput()
      implicit none
      integer :: movwin_flag

      call getInput('problem', 'B_0', B_0)
    call getInput('problem', 'psi', psi)
    call getInput('problem', 'B_amplitude', B_amp)
    call getInput('problem', 'spread', spread)
    call getInput('problem', 'mode', mode)
    call getInput('problem', 'shift_gamma', shift_gamma)
    call getInput('problem', 'temperature', TT)
    call getInput('problem', 'multiplicity_1', mult1)
      call getInput('problem', 'multiplicity_2', mult2)
      call getInput('problem', 'ramp_width', ramp_width)

    call getInput('moving_window', 'movwin', movwin_flag, 0)
    call getInput('moving_window', 'shiftstart', mw_shift_start, 0)
    call getInput('moving_window', 'shiftinterval', mw_shift_interval, 0)
    call getInput('moving_window', 'movwingam', mw_gamma_param, 0.0)
    
    use_moving_window = (movwin_flag == 1)
    mw_speed = compute_window_speed(mw_gamma_param)
    ! call getInput('problem', 'init_x_boundary', init_x_boundary)
    !AG: Overriding given init x boundary from input (exception case...remove later)
    init_x_boundary = int(5 * M_PI / mode) + 1
    fin_x_boundary = init_x_boundary + int(ramp_width * 2 * M_PI / mode)
  end subroutine userReadInput

  real(kind=dprec) function compute_window_speed(gamma_like)
    implicit none
    real(kind=dprec), intent(in) :: gamma_like
    real(kind=dprec) :: gamma_val

    if (.not. use_moving_window) then
      compute_window_speed = 0.0
      return
    end if

    if (gamma_like > 10000.0d0) then
      compute_window_speed = CC
    else if (gamma_like < 1.0d0) then
      compute_window_speed = CC * max(gamma_like, 0.0d0)
    else
      gamma_val = max(gamma_like, 1.0d0)
      compute_window_speed = CC * sqrt(1.0d0 - 1.0d0 / (gamma_val * gamma_val))
    end if
  end function compute_window_speed

  function b_wave(x_glob, y_glob, z_glob) !AG: Waveform pulse through here!!
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real :: b_wave

    ! b_wave = B_amp * sin(x_glob * mode) * exp(-0.5*(x_glob - 4*spread)**2 /spread**2)
    if(x_glob .le. 4 * M_PI / mode) then
      b_wave = (B_amp / 0.869619) * sin(x_glob * mode) * sin(x_glob * mode/4)**2 
    else 
      b_wave = 0
    endif
  end function

  function db_wave_dx(x_glob, y_glob, z_glob) !AG: Waveform pulse gradient through here!!
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real :: db_wave_dx

    ! db_wave_dx = (b_wave(x_glob + 1e-3, y_glob, z_glob) - b_wave(x_glob - 1e-3, y_glob, z_glob))/(2e-3)
    ! db_wave_dx = B_amp * (mode * cos(x_glob * mode) - sin(x_glob * mode) * (x_glob- 4*spread)/spread**2) * exp(-0.5*(x_glob - 4*spread)**2 /spread**2)
    if(x_glob .le. 4 * M_PI / mode) then
      db_wave_dx = - (B_amp / 0.869619) * 0.125*(mode*(Cos((mode*x_glob)/2.) &
      - 4*Cos(mode*x_glob) + 3*Cos((3*mode*x_glob)/2.)))
    else
      db_wave_dx = 0
    endif
  end function

  !--- Bisection Method Solver for beta+ and beta- ---!
  subroutine solve_for_beta(k, beta_plus, beta_minus)
    real, intent(in) :: k
    real, intent(out) :: beta_plus, beta_minus
    real :: x, y

    ! Call the bisection solver to solve x + y = 2 and the second equation
    call bisection_method(k, x, y, tolerance=1.0e-6, max_iterations=9000)

    ! Convert x and y into beta_plus and beta_minus
    beta_plus = (1.0 - x**2) / (1.0 + x**2)
    beta_minus = (1.0 - y**2) / (1.0 + y**2)
  end subroutine solve_for_beta

  !--- Bisection Method for Solving x+y=2 and the second equation ---!
  subroutine bisection_method(k, x, y, tolerance, max_iterations)
    real, intent(in) :: k, tolerance
    integer, intent(in) :: max_iterations
    real, intent(out) :: x, y
    real :: a, b, c, fa, fb, fc
    integer :: i

    a = 1e-6      ! Lower bound
    b = 2 - 1e-6  ! Upper bound

    fa = f(a, k)
    fb = f(b, k)

    if (fa * fb > 0.0) then
      print *, 'Error: Function does not change sign over the interval'
      stop
    end if

    ! Bisection iteration
    do i = 1, max_iterations
      c = (a + b) / 2.0
      fc = f(c, k)

      if (abs(fc) < tolerance) then
        x = c
        y = 2.0 - x
        return
      end if

      if (fa * fc < 0.0) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do

    ! If we reach here, the maximum iterations were exceeded
    print *, 'Error: Maximum iterations exceeded'
    x = (a + b) / 2.0
    y = 2.0 - x
  end subroutine bisection_method

  ! Function to evaluate f(x) = (1/x^2 - 1/(2-x)^2) - 2k
  real function f(x, k)
    real, intent(in) :: x, k
    real :: y
    y = 2.0 - x
    f = (1.0 / x**2) - (1.0 / y**2) - 2.0 * k
  end function f

  subroutine bg_density_scaled (mult1, mult2, sx_glob, xg, n_0_local)
    real, intent(in) :: mult1, mult2, sx_glob, xg
    real, intent(out) :: n_0_local
    real :: t
    real :: k
    if (mult1 .eq. mult2) then
      n_0_local = mult1

    else 
      k = log((mult1-mult2) / (0.01 * mult2))/ (fin_x_boundary - init_x_boundary)
      if (xg <= init_x_boundary) then
        n_0_local = mult1
      else
        n_0_local = mult1 + (mult2 - mult1) * &
              (1.0 - (1.0 + k*(xg - init_x_boundary)) * &
              exp(-k*(xg - init_x_boundary)))
        ! n_0_local = mult1 * exp(- k * (xg - init_x_boundary))
      endif

    endif
    
    ! if (xg <= init_x_boundary) then
    !   n_0_local = mult1
    ! else if (xg >= fin_x_boundary) then
    !   n_0_local = mult2
    ! else 
    !   t = (xg - init_x_boundary)/(fin_x_boundary - init_x_boundary)
    !   n_0_local = mult1 + (mult2 - mult1) * (3.0*t**2 - 2.0*t**3)

    ! endif

  end subroutine bg_density_scaled

  subroutine ensure_profile_cache()
    implicit none
    real :: xg, sx_glob
    real :: jA_ec, n_0_local, k, beta_plus_val, beta_minus_val
    real :: n_plus_val, n_minus_val

    if (profile_cache_ready) return

    sx_glob = REAL(global_mesh%sx)

    cached_jA_ec_max = 0.0
    do xg = 0, sx_glob
      jA_ec = -(db_wave_dx(xg, 0.5, 0.5) * sin(psi) * B_norm)/unit_ch
      cached_jA_ec_max = max(cached_jA_ec_max, abs(jA_ec))
    end do

    cached_pos_tot_dens = 0.0
    cached_ele_tot_dens = 0.0
    cached_max_density = 0.0
    do xg = 0, sx_glob
      call bg_density_scaled(mult1, mult2, sx_glob, xg, n_0_local)
      n_0_local = n_0_local * cached_jA_ec_max

      jA_ec = -(db_wave_dx(xg, 0.5, 0.5) * sin(psi) * B_norm)/unit_ch
      if (n_0_local .ne. 0.0) then
        k = jA_ec / n_0_local
      else
        k = 0.0
      end if
      call solve_for_beta(k, beta_plus_val, beta_minus_val)
      n_plus_val = n_0_local / (1.0 - beta_plus_val)
      n_minus_val = n_0_local / (1.0 - beta_minus_val)
      cached_pos_tot_dens = cached_pos_tot_dens + n_plus_val
      cached_ele_tot_dens = cached_ele_tot_dens + n_minus_val
      cached_max_density = max(cached_max_density, n_plus_val, n_minus_val)
    end do

    cached_weight_factor = (cached_pos_tot_dens + cached_ele_tot_dens) / (ppc0 * sx_glob)
    profile_cache_ready = .true.
  end subroutine ensure_profile_cache


  !--- Particle Initialization ---!
  subroutine userInitParticles()
    implicit none
    integer :: npart, accepted_particles
    real :: xg, yg, u_, v_, w_, gamma, density_value
    real :: u_parr, v_perp
    real :: sx_glob, sy_glob, max_density, pos_tot_dens, ele_tot_dens
    real :: k, beta_plus, beta_minus, n_plus, n_minus
    real :: beta_plus_val, beta_minus_val
    real :: n_plus_val, n_minus_val
    real :: jA_ec, n_0_local, jA_ec_max
    real :: weight_factor
    real :: density_normalized
    real :: ramp_factor  ! Linear ramp factor
    type(maxwellian) :: mwl
    mwl % dimension = 2
    mwl % generated = .false.
    mwl % temperature = TT 
    mwl % shift_flag = .true.
    ! mwl % shift_gamma = 1/(1 - 0.0**2)
    ! mwl % shift_dir = INT(sign(1.0, 0.0))


    ! call generateFromMaxwellian(mwl, u_parr, v_perp, w_)

    ! print *, "u_parr1:", u_parr

    ! mwl % shift_gamma = 1/(1 - 0.5**2)
    ! mwl % shift_dir = INT(sign(1.0, -0.5))

    ! call generateFromMaxwellian(mwl, u_parr, v_perp, w_)

    ! print *, "u_parr2:", u_parr
    

    sx_glob = global_mesh%sx
    sy_glob = global_mesh%sy

    call ensure_profile_cache()
    jA_ec_max = cached_jA_ec_max
    weight_factor = cached_weight_factor
    print *, "Max critical density (jA/ec max):", jA_ec_max

    ! ! Step 2: Compute the max density over the entire grid (after jA_ec_max is known)
    ! max_density = 0.0  ! Initialize the maximum density value

    ! do xg = 0, sx_glob
    !     ! ! Compute the linear ramp for n_0 from init_x_boundary
    !     ! if (xg <= init_x_boundary) then
    !     !     n_0_local = mult1 * jA_ec_max
    !     ! else
    !     !     ! ramp_factor = (xg - init_x_boundary) / (sx_glob - init_x_boundary)
    !     !     n_0_local = mult2 * jA_ec_max !+ ramp_factor * (mult2 * jA_ec_max - mult1 * jA_ec_max)
    !     ! endif  
    !     call bg_density_scaled(mult1, mult2, sx_glob, xg, n_0_local) 
    !     n_0_local = n_0_local * jA_ec_max     

    !     ! Calculate k = jA_ec / n_0
    !     jA_ec = -(db_wave_dx(xg, 0.5, 0.5) * sin(psi) * B_norm)/unit_ch
    !     k = jA_ec / (n_0_local)

    !     ! Solve for beta_plus and beta_minus
    !     call solve_for_beta(k, beta_plus_val, beta_minus_val)

    !     ! Compute n_plus and n_minus based on beta values
    !     n_plus_val = n_0_local / (1.0 - beta_plus_val)
    !     n_minus_val = n_0_local / (1.0 - beta_minus_val)

    !     ! Find the maximum density between positrons and electrons at this point
    !     max_density = max(max_density, n_plus_val, n_minus_val)
    ! end do

    ! ! Step 3: Compute the weight factor for particles (using the max_density)
    ! weight_factor = max_density / (0.50 * ppc0)
    ! print *, "Weight factor:", weight_factor

    ele_tot_dens = 0.0
    pos_tot_dens = 0.0
    pos_tot_dens = cached_pos_tot_dens
    ele_tot_dens = cached_ele_tot_dens
    max_density = cached_max_density
    print *, "Weight factor:", weight_factor
    
    ! Step 4: Inject positrons and electrons using normalized densities and rejection sampling
    npart = INT(sx_glob * sy_glob * ppc0 * pos_tot_dens/(pos_tot_dens + ele_tot_dens))
    accepted_particles = 0

    ! Loop for positrons
    do while (accepted_particles < npart) 
        ! Generate random positions for xg and yg (uniformly distributed)
        xg = random(dseed) * sx_glob
        yg = random(dseed) * sy_glob

        ! Calculate the density at this random position        
        
        call bg_density_scaled(mult1, mult2, sx_glob, xg, n_0_local)
        n_0_local = n_0_local * jA_ec_max

        jA_ec = -(db_wave_dx(xg, yg, 0.5) * sin(psi) * B_norm) / unit_ch
        k = jA_ec / (n_0_local)
        call solve_for_beta(k, beta_plus, beta_minus)
        
        n_plus = n_0_local / (1.0 - beta_plus)

        ! Normalize the density with max_density
        density_normalized = n_plus / max_density

        ! Generate a random acceptance threshold
        if (random(dseed) <= density_normalized) then
            accepted_particles = accepted_particles + 1

            ! Assign velocities for positrons
            gamma = 1.0 / sqrt(1.0 - beta_plus**2)            

            if (abs(beta_plus) .lt. 1e-3) then
              mwl % shift_gamma = 1.0              
              mwl % shift_dir = 1
              call generateFromMaxwellian(mwl, u_parr, v_perp, w_)

            else
              mwl % shift_gamma = gamma
              mwl % shift_dir = INT(sign(1.0, beta_plus))
              call generateFromMaxwellian(mwl, u_parr, v_perp, w_)
            endif            

            u_ = u_parr * cos(psi) - v_perp * sin(psi)
            v_ = u_parr * sin(psi) + v_perp * cos(psi)
            ! u_ = gamma * beta_plus * cos(psi)
            ! v_ = gamma * beta_plus * sin(psi)
            ! w_ = 0.0

            ! Inject the positron into the system
            call injectParticleGlobally(1, xg, yg, 0.5, u_, v_, w_, weight_factor)
        end if
    end do

    ! Loop for electrons
    npart = INT(sx_glob * sy_glob * ppc0 * ele_tot_dens/(pos_tot_dens + ele_tot_dens))
    accepted_particles = 0
    do while (accepted_particles < npart) 
        ! Generate random positions for xg and yg (uniformly distributed)
        xg = random(dseed) * sx_glob
        yg = random(dseed) * sy_glob

        ! Calculate the density at this random position   
        call bg_density_scaled(mult1, mult2, sx_glob, xg, n_0_local) 
        n_0_local = n_0_local * jA_ec_max

        jA_ec = -(db_wave_dx(xg, yg, 0.5) * sin(psi) * B_norm) / unit_ch
        k = jA_ec / (n_0_local)
        call solve_for_beta(k, beta_plus, beta_minus)
        
        n_minus = n_0_local / (1.0 - beta_minus)

        ! Normalize the density with max_density
        density_normalized = n_minus / max_density

        ! Generate a random acceptance threshold
        if (random(dseed) <= density_normalized) then
            accepted_particles = accepted_particles + 1

            ! Assign velocities for electrons with temperature
            gamma = 1.0 / sqrt(1.0 - beta_minus**2)             

            if (abs(beta_minus) .lt. 1e-3) then
              mwl % shift_gamma = 1.0
              mwl % shift_dir = 1
              call generateFromMaxwellian(mwl, u_parr, v_perp, w_)

            else
              mwl % shift_gamma = gamma
              mwl % shift_dir = INT(sign(1.0, beta_minus))
              call generateFromMaxwellian(mwl, u_parr, v_perp, w_)
            endif                                    

            u_ = u_parr * cos(psi) - v_perp * sin(psi)
            v_ = u_parr * sin(psi) + v_perp * cos(psi)

            ! if (accepted_particles .eq. 1) then
            !   print *, "u_1:", u_, "beta_minus:", beta_minus
            !   ! print *, "sign expr:", sign(1.0, beta_minus)
            ! endif
            ! u_ = gamma * beta_minus * cos(psi)
            ! v_ = gamma * beta_minus * sin(psi)
            ! if (accepted_particles .eq. 1) then
            !   print *, "u_2:", u_
            ! endif
            ! w_ = 0.0

            ! Inject the electron into the system
            call injectParticleGlobally(2, xg, yg, 0.5, u_, v_, w_, weight_factor)
        end if
    end do

  end subroutine userInitParticles

  function e_parr(x, jA_ec_max)
    real, intent(in), optional :: x, jA_ec_max
    real :: e_parr
    real :: sx_glob
    real :: jA_ec, n_0_local, q1, q2, k, beta_plus_1, beta_plus_2, beta_minus

    sx_glob = REAL(global_mesh%sx)
    call bg_density_scaled(mult1, mult2, sx_glob, x - 1e-1, n_0_local) 
    n_0_local = n_0_local * jA_ec_max 
    ! Calculate k = jA_ec / n_0
    jA_ec = -(db_wave_dx(x - 1e-1, 0.5, 0.5) * sin(psi) * B_norm)/unit_ch
    k = jA_ec / (n_0_local)
    ! Solve for beta_plus and beta_minus
    call solve_for_beta(k, beta_plus_1, beta_minus)
    q1 = SQRT((1 - beta_plus_1)/(1 + beta_plus_1))

    call bg_density_scaled(mult1, mult2, sx_glob, x + 1e-1, n_0_local) 
    n_0_local = n_0_local * jA_ec_max 
    ! Calculate k = jA_ec / n_0
    jA_ec = -(db_wave_dx(x + 1e-1, 0.5, 0.5) * sin(psi) * B_norm)/unit_ch
    k = jA_ec / (n_0_local)
    ! Solve for beta_plus and beta_minus
    call solve_for_beta(k, beta_plus_2, beta_minus)
    q2 = SQRT((1 - beta_plus_2)/(1 + beta_plus_2))

    e_parr = - cos(psi) * c_omp / sqrt(sigma) * (q2 - q1)/(2e-1) 
  end function


  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    real :: sx_glob, sy_glob
    ! real :: kx
    real :: ex_norm, ey_norm, ez_norm, exyz_norm
    real :: bx_norm, by_norm, bz_norm, bxyz_norm
    real :: xg, jA_ec, jA_ec_max, e_parr1, e_parr2
    
    
    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)
    
    ! Step 1: Compute jA_ec_max over the entire grid
    jA_ec_max = 0.0  ! Initialize jA_ec_max

    do xg = 0, sx_glob
        ! Compute the current density jA_ec at this position
        jA_ec = -(db_wave_dx(xg, 0.5, 0.5) * sin(psi) * B_norm)/unit_ch
        jA_ec_max = max(jA_ec_max, abs(jA_ec))  ! Store the maximum value of jA_ec
    end do
    
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = B_0 * cos(psi); by(:, :, :) = B_0 * sin(psi); bz(:, :, :) = 0


    ! kx = mode
    do i = 0, this_meshblock % ptr % sx - 1
      i_glob = i + this_meshblock % ptr % x0

      e_parr1 = e_parr(REAL(i_glob), jA_ec_max)
      e_parr2 = e_parr(REAL(i_glob - 0.5), jA_ec_max)

      do j = 0, this_meshblock % ptr % sy - 1
        j_glob = j + this_meshblock % ptr % y0
        do k = 0, this_meshblock % ptr % sz - 1
          k_glob = k + this_meshblock % ptr % z0
          ex(i, j, k) = -sin(psi) * b_wave(REAL(i_glob), REAL(j_glob-0.5), REAL(k_glob-0.5)) + cos(psi) * e_parr1
          ey(i, j, k) =  cos(psi) * b_wave(REAL(i_glob-0.5), REAL(j_glob), REAL(k_glob-0.5)) + sin(psi) * e_parr2
          bz(i, j, k) = b_wave(REAL(i_glob), REAL(j_glob), REAL(k_glob-0.5))
        end do
      end do
    end do

    ! ex_p(:, :, :) = ex(:, :, :); ey_p(:, :, :) = ey(:, :, :); ez_p(:, :, :) = ez(:, :, :)
    ! bx_p(:, :, :) = bx(:, :, :); by_p(:, :, :) = by(:, :, :); bz_p(:, :, :) = bz(:, :, :)

  end subroutine userInitFields

  !--- driving ------------------------------------------------!
  ! subroutine userCurrentDeposit(step)
  !   implicit none
  !   integer, optional, intent(in) :: step
  !   ! called after particles move and deposit ...
  !   ! ... and before the currents are added to the electric field
  ! end subroutine userCurrentDeposit

  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
    integer :: s, ti, tj, tk, p
    real :: u_, v_, w_, x_l, y_l, z_l, x_g, inv_gamma, gamma
    real :: tfrac, xcolis, xnew, x0
    
    ! type(maxwellian) :: mwl
    ! mwl % dimension = 2
    ! mwl % generated = .false.
    ! mwl % temperature = TT
    ! mwl % shift_gamma = shift_gamma
    ! mwl % shift_flag = .true.
    ! mwl % shift_dir = -1

    do s = 1, nspec
     do ti = 1, species(s) % tile_nx
       do tj = 1, species(s) % tile_ny
         do tk = 1, species(s) % tile_nz
           do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
             x_l = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p)) + species(s) % prtl_tile(ti, tj, tk) % dx(p)
             x_g = x_l + REAL(this_meshblock % ptr % x0)
             if (x_g .lt. 0) then
               u_ = species(s) % prtl_tile(ti, tj, tk) % u(p)
               v_ = species(s) % prtl_tile(ti, tj, tk) % v(p)
               w_ = species(s) % prtl_tile(ti, tj, tk) % w(p)
               gamma = sqrt(1.0 + u_**2 + v_**2 + w_**2)
               inv_gamma = 1.0 / gamma
               x0 = x_l - u_ * inv_gamma * CC
               tfrac = abs((x0 - 0) / (u_ * inv_gamma * CC))
               if (tfrac .lt. 1) then
                 xcolis = x0 + u_ * inv_gamma * CC * tfrac
               end if
               y_l = REAL(species(s) % prtl_tile(ti, tj, tk) % yi(p)) + species(s) % prtl_tile(ti, tj, tk) % dy(p)
               z_l = REAL(species(s) % prtl_tile(ti, tj, tk) % zi(p)) + species(s) % prtl_tile(ti, tj, tk) % dz(p)
               call depositCurrentsFromSingleParticle(s, species(s) % prtl_tile(ti, tj, tk), p, &
                                                      xcolis, y_l, z_l, x_l, y_l, z_l, -1.0)
               ! spawning new particle from same thermal dist.
              !  call generateFromMaxwellian(mwl, u_, v_, w_)
               !u_ = 0; v_ = 0; w_ = 0
               species(s) % prtl_tile(ti, tj, tk) % u(p) = -u_
               species(s) % prtl_tile(ti, tj, tk) % v(p) = v_
               species(s) % prtl_tile(ti, tj, tk) % w(p) = w_
              !  u_ = species(s) % prtl_tile(ti, tj, tk) % u(p)
               tfrac = min(abs((x_l - xcolis) / max(abs(x_l - x0), 1e-9)), 1.0)
               xnew = xcolis + (- u_) * inv_gamma * CC * tfrac
               species(s) % prtl_tile(ti, tj, tk) % xi(p) = INT(FLOOR(xnew), 2)
               species(s) % prtl_tile(ti, tj, tk) % dx(p) = xnew - REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p))
               call depositCurrentsFromSingleParticle(s, species(s) % prtl_tile(ti, tj, tk), p, &
                                                      xcolis, y_l, z_l, xnew, y_l, z_l, 1.0)
             end if

             if (x_g .gt. global_mesh%sx) then
              u_ = species(s) % prtl_tile(ti, tj, tk) % u(p)
              v_ = species(s) % prtl_tile(ti, tj, tk) % v(p)
              w_ = species(s) % prtl_tile(ti, tj, tk) % w(p)
              gamma = sqrt(1.0 + u_**2 + v_**2 + w_**2)
              inv_gamma = 1.0 / gamma
              x0 = x_l - u_ * inv_gamma * CC
              tfrac = abs((x0 - global_mesh%sx) / (u_ * inv_gamma * CC))
              if (tfrac .lt. 1) then
                xcolis = x0 + u_ * inv_gamma * CC * tfrac
              end if
              y_l = REAL(species(s) % prtl_tile(ti, tj, tk) % yi(p)) + species(s) % prtl_tile(ti, tj, tk) % dy(p)
              z_l = REAL(species(s) % prtl_tile(ti, tj, tk) % zi(p)) + species(s) % prtl_tile(ti, tj, tk) % dz(p)
              call depositCurrentsFromSingleParticle(s, species(s) % prtl_tile(ti, tj, tk), p, &
                                                     xcolis, y_l, z_l, x_l, y_l, z_l, -1.0)
              ! spawning new particle from same thermal dist.
              ! call generateFromMaxwellian(mwl, u_, v_, w_)
              !u_ = 0; v_ = 0; w_ = 0
              species(s) % prtl_tile(ti, tj, tk) % u(p) = -u_
              species(s) % prtl_tile(ti, tj, tk) % v(p) = v_
              species(s) % prtl_tile(ti, tj, tk) % w(p) = w_
             !  u_ = species(s) % prtl_tile(ti, tj, tk) % u(p)
              tfrac = min(abs((x_l - xcolis) / max(abs(x_l - x0), 1e-9)), 1.0)
              xnew = xcolis + (- u_) * inv_gamma * CC * tfrac
              species(s) % prtl_tile(ti, tj, tk) % xi(p) = INT(FLOOR(xnew), 2)
              species(s) % prtl_tile(ti, tj, tk) % dx(p) = xnew - REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p))
              call depositCurrentsFromSingleParticle(s, species(s) % prtl_tile(ti, tj, tk), p, &
                                                     xcolis, y_l, z_l, xnew, y_l, z_l, 1.0)
            end if
           end do
         end do
       end do
     end do
    end do
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
    ! ... dummy loop ...
    ! integer :: s, ti, tj, tk, p
    ! do s = 1, nspec
    !   do ti = 1, species(s)%tile_nx
    !     do tj = 1, species(s)%tile_ny
    !       do tk = 1, species(s)%tile_nz
    !         do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
    !           ...
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do
  end subroutine userDriveParticles
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ! real :: kx, ky, kz

    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical :: updateE_, updateB_

    if (present(updateE)) then
      updateE_ = updateE
    else
      updateE_ = .true.
    end if

    if (present(updateB)) then
      updateB_ = updateB
    else
      updateB_ = .true.
    end if

  end subroutine userFieldBoundaryConditions

  subroutine userFillNewRegion(xmin, xmax)
    implicit none
    real, intent(in) :: xmin, xmax
    real :: jA_ec_max, max_density, pos_tot_dens, ele_tot_dens, weight_factor
    real :: slab_min, slab_max, slab_extent, sy_glob
    real :: xg, yg, density_normalized, kappa, beta_plus, beta_minus
    real :: n_plus, n_minus, n_0_local, jA_ec
    integer :: npart, accepted_particles
    type(maxwellian) :: mwl
    real :: u_parr, v_perp, u_, v_, w_, gamma
    integer :: i, j, k, i_glob, j_glob, k_glob
    real :: e_parr1, e_parr2

    call ensure_profile_cache()

    jA_ec_max = cached_jA_ec_max
    max_density = cached_max_density
    pos_tot_dens = cached_pos_tot_dens
    ele_tot_dens = cached_ele_tot_dens
    weight_factor = cached_weight_factor

    slab_min = max(xmin, REAL(this_meshblock % ptr % x0))
    slab_max = min(xmax, REAL(this_meshblock % ptr % x0 + this_meshblock % ptr % sx - 1))

    if (slab_max < slab_min) return

    do i = 0, this_meshblock % ptr % sx - 1
      i_glob = i + this_meshblock % ptr % x0
      if ((REAL(i_glob) < slab_min) .or. (REAL(i_glob) > slab_max)) cycle

      e_parr1 = e_parr(REAL(i_glob), jA_ec_max)
      e_parr2 = e_parr(REAL(i_glob - 0.5), jA_ec_max)

      do j = 0, this_meshblock % ptr % sy - 1
        j_glob = j + this_meshblock % ptr % y0
        do k = 0, this_meshblock % ptr % sz - 1
          k_glob = k + this_meshblock % ptr % z0
          ex(i, j, k) = -sin(psi) * b_wave(REAL(i_glob), REAL(j_glob - 0.5), REAL(k_glob - 0.5)) + cos(psi) * e_parr1
          ey(i, j, k) =  cos(psi) * b_wave(REAL(i_glob - 0.5), REAL(j_glob), REAL(k_glob - 0.5)) + sin(psi) * e_parr2
          bz(i, j, k) = b_wave(REAL(i_glob), REAL(j_glob), REAL(k_glob - 0.5))
        end do
      end do
    end do

    slab_extent = slab_max - slab_min + 1.0
    if (slab_extent <= 0.0) return

    sy_glob = REAL(global_mesh % sy)
    mwl % dimension = 2
    mwl % generated = .false.
    mwl % temperature = TT
    mwl % shift_flag = .true.

    npart = INT(slab_extent * sy_glob * ppc0 * pos_tot_dens / (pos_tot_dens + ele_tot_dens))
    accepted_particles = 0
    do while (accepted_particles < npart)
      xg = slab_min + random(dseed) * slab_extent
      yg = random(dseed) * sy_glob

      call bg_density_scaled(mult1, mult2, REAL(global_mesh % sx), xg, n_0_local)
      n_0_local = n_0_local * jA_ec_max

      jA_ec = -(db_wave_dx(xg, yg, 0.5) * sin(psi) * B_norm) / unit_ch
      if (n_0_local .ne. 0.0) then
        kappa = jA_ec / n_0_local
      else
        kappa = 0.0
      end if
      call solve_for_beta(kappa, beta_plus, beta_minus)

      n_plus = n_0_local / (1.0 - beta_plus)
      density_normalized = n_plus / max_density

      if (random(dseed) <= density_normalized) then
        accepted_particles = accepted_particles + 1

        gamma = 1.0 / sqrt(1.0 - beta_plus**2)
        if (abs(beta_plus) .lt. 1e-3) then
          mwl % shift_gamma = 1.0
          mwl % shift_dir = 1
        else
          mwl % shift_gamma = gamma
          mwl % shift_dir = INT(sign(1.0, beta_plus))
        end if
        call generateFromMaxwellian(mwl, u_parr, v_perp, w_)

        u_ = u_parr * cos(psi) - v_perp * sin(psi)
        v_ = u_parr * sin(psi) + v_perp * cos(psi)

        call injectParticleGlobally(1, xg, yg, 0.5, u_, v_, w_, weight_factor)
      end if
    end do

    npart = INT(slab_extent * sy_glob * ppc0 * ele_tot_dens / (pos_tot_dens + ele_tot_dens))
    accepted_particles = 0
    do while (accepted_particles < npart)
      xg = slab_min + random(dseed) * slab_extent
      yg = random(dseed) * sy_glob

      call bg_density_scaled(mult1, mult2, REAL(global_mesh % sx), xg, n_0_local)
      n_0_local = n_0_local * jA_ec_max

      jA_ec = -(db_wave_dx(xg, yg, 0.5) * sin(psi) * B_norm) / unit_ch
      if (n_0_local .ne. 0.0) then
        kappa = jA_ec / n_0_local
      else
        kappa = 0.0
      end if
      call solve_for_beta(kappa, beta_plus, beta_minus)

      n_minus = n_0_local / (1.0 - beta_minus)
      density_normalized = n_minus / max_density

      if (random(dseed) <= density_normalized) then
        accepted_particles = accepted_particles + 1

        gamma = 1.0 / sqrt(1.0 - beta_minus**2)
        if (abs(beta_minus) .lt. 1e-3) then
          mwl % shift_gamma = 1.0
          mwl % shift_dir = 1
        else
          mwl % shift_gamma = gamma
          mwl % shift_dir = INT(sign(1.0, beta_minus))
        end if
        call generateFromMaxwellian(mwl, u_parr, v_perp, w_)

        u_ = u_parr * cos(psi) - v_perp * sin(psi)
        v_ = u_parr * sin(psi) + v_perp * cos(psi)

        call injectParticleGlobally(2, xg, yg, 0.5, u_, v_, w_, weight_factor)
      end if
    end do
  end subroutine userFillNewRegion

#include "optional.F"






