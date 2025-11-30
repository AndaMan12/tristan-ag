module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private :: B_0, B_amp, mode, spread, shift_gamma, TT, psi, mult1, mult2
  real, private :: jA_ec_max, weight
  integer, private :: init_x_boundary
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  ! private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'B_0', B_0)
    call getInput('problem', 'psi', psi)
    call getInput('problem', 'B_amplitude', B_amp)
    call getInput('problem', 'spread', spread)
    call getInput('problem', 'mode', mode)
    call getInput('problem', 'shift_gamma', shift_gamma)
    call getInput('problem', 'temperature', TT)
    call getInput('problem', 'multiplicity_1', mult1)
    call getInput('problem', 'multiplicity_2', mult2)    
    call getInput('problem', 'init_x_boundary', init_x_boundary)
  end subroutine userReadInput

  function b_wave(x_glob, y_glob, z_glob) !AG: Waveform pulse through here!!
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real :: b_wave

    ! b_wave = B_amp * sin(x_glob * mode) * exp(-0.5*(x_glob - 4*spread)**2 /spread**2)
    if(x_glob .le. 2 * M_PI / mode) then
      b_wave = B_amp * (8 * SQRT(3.0)/9) * sin(x_glob * mode) * sin(x_glob * mode/2)**2 
    else 
      b_wave = 0
    endif
  end function

  function db_wave_dx(x_glob, y_glob, z_glob) !AG: Waveform pulse gradient through here!!
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real :: db_wave_dx

    ! db_wave_dx = (b_wave(x_glob + 1e-3, y_glob, z_glob) - b_wave(x_glob - 1e-3, y_glob, z_glob))/(2e-3)
    ! db_wave_dx = B_amp * (mode * cos(x_glob * mode) - sin(x_glob * mode) * (x_glob- 4*spread)/spread**2) * exp(-0.5*(x_glob - 4*spread)**2 /spread**2)
    if(x_glob .le. 2 * M_PI / mode) then
      db_wave_dx = B_amp * (8 * SQRT(3.0)/9) * (mode/2) * (cos(mode * x_glob) - cos(2 * mode * x_glob))
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
    call bisection_method(k, x, y, tolerance=1.0e-5, max_iterations=1000)

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

    a = 0.0001  ! Lower bound
    b = 1.9999  ! Upper bound

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

  !--- Particle Initialization ---!
  subroutine userInitParticles()
    implicit none
    integer :: npart, accepted_particles
    real :: xg, yg, u_, v_, w_, gamma, density_value
    real :: sx_glob, sy_glob, max_density
    real :: k, beta_plus, beta_minus, n_plus, n_minus
    real :: beta_plus_val, beta_minus_val
    real :: n_plus_val, n_minus_val
    real :: jA_ec, n_0_local, jA_ec_max
    real :: weight_factor
    real :: density_normalized
    real :: ramp_factor  ! Linear ramp factor

    !AG: Overriding given init x boundary from input (exception case...remove later)
    init_x_boundary = 2 * M_PI / mode
    
    sx_glob = global_mesh%sx
    sy_glob = global_mesh%sy

    ! Step 1: Compute jA_ec_max over the entire grid
    jA_ec_max = 0.0  ! Initialize jA_ec_max

    do xg = 0, sx_glob
        ! Compute the current density jA_ec at this position
        jA_ec = -(db_wave_dx(xg, 0.5, 0.5) * sin(psi) * B_norm)/unit_ch
        jA_ec_max = max(jA_ec_max, abs(jA_ec))  ! Store the maximum value of jA_ec
    end do

    ! Step 2: Compute the max density over the entire grid (after jA_ec_max is known)
    max_density = 0.0  ! Initialize the maximum density value

    do xg = 0, sx_glob
        ! Compute the linear ramp for n_0 from init_x_boundary
        if (xg <= init_x_boundary) then
            n_0_local = mult1 * jA_ec_max
        else
            ! ramp_factor = (xg - init_x_boundary) / (sx_glob - init_x_boundary)
            n_0_local = mult2 * jA_ec_max !mult1 * jA_ec_max + ramp_factor * (mult2 * jA_ec_max - mult1 * jA_ec_max)
        endif        

        ! Calculate k = jA_ec / n_0
        jA_ec = -(db_wave_dx(xg, 0.5, 0.5) * sin(psi) * B_norm)/unit_ch
        k = jA_ec / n_0_local

        ! Solve for beta_plus and beta_minus
        call solve_for_beta(k, beta_plus_val, beta_minus_val)

        ! Compute n_plus and n_minus based on beta values
        n_plus_val = n_0_local / (1.0 - beta_plus_val)
        n_minus_val = n_0_local / (1.0 - beta_minus_val)

        ! Find the maximum density between positrons and electrons at this point
        max_density = max(max_density, n_plus_val, n_minus_val)
    end do

    ! Step 3: Compute the weight factor for particles (using the max_density)
    weight_factor = max_density / (0.50 * ppc0)
    print *, "Weight factor:", weight_factor

    ! Step 4: Inject positrons and electrons using normalized densities and rejection sampling
    npart = INT(sx_glob * sy_glob * 0.50 * ppc0)
    accepted_particles = 0

    ! Loop for positrons
    do while (accepted_particles < npart) 
        ! Generate random positions for xg and yg (uniformly distributed)
        xg = random(dseed) * sx_glob
        yg = random(dseed) * sy_glob

        ! Calculate the density at this random position
        jA_ec = -(db_wave_dx(xg, yg, 0.5) * sin(psi) * B_norm) / unit_ch

        ! Compute the linear ramp for n_0 from init_x_boundary
        if (xg <= init_x_boundary) then
            n_0_local = mult1 * jA_ec_max
        else
            ! ramp_factor = (xg - init_x_boundary) / (sx_glob - init_x_boundary)
            n_0_local = mult2 * jA_ec_max !mult1 * jA_ec_max + ramp_factor * (mult2 * jA_ec_max - mult1 * jA_ec_max)
        endif

        k = jA_ec / n_0_local
        call solve_for_beta(k, beta_plus, beta_minus)
        ! print *, beta_plus, beta_minus
        n_plus = n_0_local / (1.0 - beta_plus)

        ! Normalize the density with max_density
        density_normalized = n_plus / max_density

        ! Generate a random acceptance threshold
        if (random(dseed) <= density_normalized) then
            accepted_particles = accepted_particles + 1

            ! Assign velocities for positrons
            gamma = 1.0 / sqrt(1.0 - beta_plus**2)
            u_ = gamma * beta_plus * cos(psi)
            v_ = gamma * beta_plus * sin(psi)
            w_ = 0.0

            ! Inject the positron into the system
            call injectParticleGlobally(1, xg, yg, 0.5, u_, v_, w_, weight_factor)
        end if
    end do

    ! Loop for electrons
    accepted_particles = 0
    do while (accepted_particles < npart) 
        ! Generate random positions for xg and yg (uniformly distributed)
        xg = random(dseed) * sx_glob
        yg = random(dseed) * sy_glob

        ! Calculate the density at this random position
        jA_ec = -(db_wave_dx(xg, yg, 0.5) * sin(psi) * B_norm) / unit_ch

        ! Compute the linear ramp for n_0 from init_x_boundary
        if (xg <= init_x_boundary) then
            n_0_local = mult1 * jA_ec_max
        else
            ramp_factor = (xg - init_x_boundary) / (sx_glob - init_x_boundary)
            n_0_local = mult1 * jA_ec_max + ramp_factor * (mult2 * jA_ec_max - mult1 * jA_ec_max)
        endif

        k = jA_ec / n_0_local
        call solve_for_beta(k, beta_plus, beta_minus)
        n_minus = n_0_local / (1.0 - beta_minus)

        ! Normalize the density with max_density
        density_normalized = n_minus / max_density

        ! Generate a random acceptance threshold
        if (random(dseed) <= density_normalized) then
            accepted_particles = accepted_particles + 1

            ! Assign velocities for electrons
            gamma = 1.0 / sqrt(1.0 - beta_minus**2)
            u_ = gamma * beta_minus * cos(psi)
            v_ = gamma * beta_minus * sin(psi)
            w_ = 0.0

            ! Inject the electron into the system
            call injectParticleGlobally(2, xg, yg, 0.5, u_, v_, w_, weight_factor)
        end if
    end do

  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    real :: sx_glob, sy_glob
    real :: kx
    real :: ex_norm, ey_norm, ez_norm, exyz_norm
    real :: bx_norm, by_norm, bz_norm, bxyz_norm
    
    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)
    
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = B_0 * cos(psi); by(:, :, :) = B_0 * sin(psi); bz(:, :, :) = 0

    kx = mode
    do i = 0, this_meshblock % ptr % sx - 1
      i_glob = i + this_meshblock % ptr % x0
      do j = 0, this_meshblock % ptr % sy - 1
        j_glob = j + this_meshblock % ptr % y0
        do k = 0, this_meshblock % ptr % sz - 1
          k_glob = k + this_meshblock % ptr % z0
          ex(i, j, k) = -sin(psi) * b_wave(REAL(i_glob), REAL(j_glob-0.5), REAL(k_glob-0.5))
          ey(i, j, k) =  cos(psi) * b_wave(REAL(i_glob-0.5), REAL(j_glob), REAL(k_glob-0.5))
          bz(i, j, k) = b_wave(REAL(i_glob), REAL(j_glob), REAL(k_glob-0.5))
        end do
      end do
    end do

  end subroutine userInitFields

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
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
  
#include "optional.F"
  !............................................................!
end module m_userfile










! module m_userfile
!   use m_globalnamespace
!   use m_aux
!   use m_readinput
!   use m_domain
!   use m_particles
!   use m_fields
!   use m_thermalplasma
!   use m_particlelogistics
!   implicit none

!   !--- PRIVATE variables -----------------------------------------!
!   real, private :: B_0, B_amp, mode, spread, shift_gamma, TT, psi, mult1, mult2, beta
!   real, private :: jA_ec_max, weight
!   integer, private :: init_x_boundary
!   !...............................................................!

!   !--- PRIVATE functions -----------------------------------------!
!   private :: userSpatialDistribution
!   !...............................................................!
! contains
!   !--- initialization -----------------------------------------!
!   subroutine userReadInput()
!     implicit none
!     call getInput('problem', 'B_0', B_0)
!     call getInput('problem', 'psi', psi)
!     call getInput('problem', 'B_amplitude', B_amp)
!     call getInput('problem', 'spread', spread)
!     call getInput('problem', 'mode', mode)
!     call getInput('problem', 'shift_gamma', shift_gamma)
!     call getInput('problem', 'temperature', TT)
!     call getInput('problem', 'multiplicity_1', mult1)
!     call getInput('problem', 'multiplicity_2', mult2)
!     ! call getInput('problem', 'weights', weight)
!     call getInput('problem', 'beta', beta)
!     call getInput('problem', 'init_x_boundary', init_x_boundary)
!   end subroutine userReadInput

!   function b_wave(x_glob, y_glob, z_glob) !AG: Waveform pulse through here!!
!     real, intent(in), optional :: x_glob, y_glob, z_glob
!     real :: b_wave

!     b_wave = B_amp * sin(x_glob * mode) * exp(-0.5*(x_glob - 4*spread)**2 /spread**2)
!   end function

!   function db_wave_dx(x_glob, y_glob, z_glob) !AG: Waveform pulse gradient through here!!
!     real, intent(in), optional :: x_glob, y_glob, z_glob
!     real :: db_wave_dx

!     db_wave_dx = (b_wave(x_glob + 1e-3, y_glob, z_glob) - b_wave(x_glob - 1e-3, y_glob, z_glob))/(2e-3)
!   end function

!   function find_max_abs_value(x_start, x_end, x_step, y_glob, z_glob) result(max_abs_x)
!     real, intent(in) :: x_start, x_end, x_step, y_glob, z_glob
!     real :: max_abs_x
!     real :: x_current, max_abs_value, db_wave_current

!     ! Initialize the search
!     max_abs_value = 0.0
!     max_abs_x = x_start
!     x_current = x_start

!     ! Loop through x values to find the maximum absolute value of db_wave_dx
!     do while (x_current <= x_end)
!         ! Compute the current db_wave_dx and its absolute value
!         db_wave_current = abs(db_wave_dx(x_current, y_glob, z_glob))

!         ! Check if this is the largest absolute value so far
!         if (db_wave_current > max_abs_value) then
!             max_abs_value = db_wave_current
!             max_abs_x = x_current
!         end if

!         ! Move to the next x value
!         x_current = x_current + x_step
!     end do
!   end function

  
!   function userSpatialDistribution(x_glob, y_glob, z_glob, &
!     dummy1, dummy2, dummy3)
!   real :: userSpatialDistribution
!   real, intent(in), optional :: x_glob, y_glob, z_glob
!   real, intent(in), optional :: dummy1, dummy2, dummy3
!   real           :: sx_glob, jA_ec

!   sx_glob = REAL(global_mesh%sx)
!   jA_ec = (-db_wave_dx(x_glob, y_glob, z_glob) * sin(psi) * B_norm)/unit_ch

!   if (x_glob <= init_x_boundary) then
!   if (dummy1 .eq. 1.0) then  ! AG: positron density distribution    
!   userSpatialDistribution = (mult1 * dummy2 + 0.5 * (1/beta + 1) * jA_ec)/(mult1 * dummy2 + 0.5 * (1/beta + 1) * dummy2)
!   else                  ! AG: electron density distribution
!   userSpatialDistribution = (mult1 * dummy2 + 0.5 * (1/beta - 1) * jA_ec)/(mult1 * dummy2 + 0.5 * (1/beta + 1) * dummy2)
!   endif

!   else                    ! AG: density ramp down outside initialisation region
!   userSpatialDistribution = 1/(1 + 0.5 * (1/beta + 1) / mult1) * &
!                             (1 - (x_glob - init_x_boundary)/(sx_glob - init_x_boundary) * (1 - mult2/mult1))
!   endif
!   return
!   end function

!   subroutine userInitParticles()
!     implicit none
!     integer :: npart, accepted_particles
!     real :: xg, yg, u_,  v_,  w_, gamma, density_value
!     real :: sx_glob, sy_glob
!     real :: u_acceptance
!     procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
!     spat_distr_ptr => userSpatialDistribution
    
!     sx_glob = global_mesh%sx
!     sy_glob = global_mesh%sy
    
!     ! AG: Critical number density for the pulse being setup.
!     jA_ec_max = find_max_abs_value(0.0, REAL(init_x_boundary), 0.01, 0.0, 0.0) * sin(psi) * B_norm / unit_ch
!     ! AG: Weight of each macroparticle to normalise the number density.
!     weight = (mult1 * jA_ec_max + 0.5 * (1/beta + 1) * jA_ec_max)

!     npart = INT(sx_glob * sy_glob * 0.50 * ppc0)

!     accepted_particles = 0
!     ! AG: Loop for positrons
!     ! Loop until we have accepted the required number of particles
!     do while (accepted_particles < npart) 
!         ! Generate random positions for xg and yg (uniformly distributed)        
!         xg = random(dseed) * sx_glob        
!         yg = random(dseed) * sy_glob

!         ! Evaluate the PDF at this random position using userSpatialDistribution
!         density_value = spat_distr_ptr(x_glob=xg, y_glob=yg, z_glob=0.5, dummy1 = 1.0, dummy2 = jA_ec_max)

!         ! Generate a random acceptance threshold between 0 and 1
!         u_acceptance = random(dseed)

!         ! If the random number is less than the PDF value, accept the particle
!         if (u_acceptance <= density_value) then
!             accepted_particles = accepted_particles + 1

!             ! Set velocities as functions of the position (not needed in this try)

!             ! Set velocities as functions of the position (not needed in this try)
!           if (xg < init_x_boundary) then
!             gamma = 1/sqrt(1 - beta**2)
!             u_ = gamma * beta * cos(psi)
!             v_ = gamma * beta * sin(psi)
!             w_ = 0
!           else
!             u_ = 0
!             v_ = 0
!             w_ = 0
!           endif
!             ! Inject the particle into the system
!             call injectParticleGlobally(1, xg, yg, 0.5, u_, v_, w_, weight)
!         end if
!     end do

!     accepted_particles = 0
    
!     ! AG: Loop for electrons
!     do while (accepted_particles < npart) 
!       ! Generate random positions for xg and yg (uniformly distributed)        
!       xg = random(dseed) * sx_glob        
!       yg = random(dseed) * sy_glob

!       ! Evaluate the PDF at this random position using userSpatialDistribution
!       density_value = spat_distr_ptr(x_glob=xg, y_glob=yg, z_glob=0.5, dummy1 = 2.0, dummy2 = jA_ec_max)

!       ! Generate a random acceptance threshold between 0 and 1
!       u_acceptance = random(dseed)

!       ! If the random number is less than the PDF value, accept the particle
!       if (u_acceptance <= density_value) then
!           accepted_particles = accepted_particles + 1

!           ! Set velocities as functions of the position (not needed in this try)
!           if (xg < init_x_boundary) then
!             gamma = 1/sqrt(1 - beta**2)
!             u_ = -gamma * beta * cos(psi)
!             v_ = -gamma * beta * sin(psi)
!             w_ = 0
!           else
!             u_ = 0
!             v_ = 0
!             w_ = 0
!           endif
!           ! Inject the particle into the system
!           call injectParticleGlobally(2, xg, yg, 0.5, u_, v_, w_, weight)
!       end if
!   end do

!   end subroutine userInitParticles

  
!   ! subroutine userInitParticles()
!   !   implicit none
!   !   ! //AG: Initialising variables for injecting thermal plasma
!   !   real           :: ne, np           ! no. density of electron and positron
!   !   type(region)   :: fill_region      ! defines a region (on a domain or meshblock) where the plasma will be sprinkled
!   !   real           :: sx_glob, sy_glob, sz_glob
!   !   procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
!   !   spat_distr_ptr => userSpatialDistribution
!   !   ! *************************************************************************************
    
!   !   ne = 0.5 * ppc0
!   !   np = 0.5 * ppc0
!   !   sx_glob = REAL(global_mesh%sx)
!   !   sy_glob = REAL(global_mesh%sy)
!   !   !sz_glob = REAL(global_mesh%sz)
!   !   fill_region%x_min = REAL(0)
!   !   fill_region%x_max = REAL(sx_glob)
!   !   fill_region%y_min = REAL(0)
!   !   fill_region%y_max = REAL(sy_glob)
!   !   !fill_region%z_min = REAL(0)
!   !   !fill_region%z_max = REAL(sz_glob)

!   !   ! //AG: Parameters of "fillRegionWithThermalPlasma()" for reference
!   !   ! fillRegionWithThermalPlasma(fill_region, fill_species, num_species, ndens_sp, &
!   !   ! temperature, shift_gamma, shift_dir, zero_current, &
!   !   ! dimension, weights, spat_distr_ptr, &
!   !   ! dummy1, dummy2, dummy3)

!   !   ! filling the positrons...
!   !   call fillRegionWithThermalPlasma(fill_region, fill_species = (/1/), num_species = 1,&
!   !    ndens_sp = np, temperature = TT, shift_gamma = shift_gamma, shift_dir = +1, dimension = 2, &
!   !    weights = weight, spat_distr_ptr = spat_distr_ptr)
    
!   !   ! filling the electrons...
!   !   call fillRegionWithThermalPlasma(fill_region, fill_species = (/2/), num_species = 1,&
!   !    ndens_sp = ne, temperature = TT, shift_gamma = shift_gamma, shift_dir = -1, dimension = 2, &
!   !    weights = weight, spat_distr_ptr = spat_distr_ptr)

! !     ! //AG: Old code to inject particles one by one with a random number generator....

! ! !     u_ = shift_gamma * sqrt(1.0 - shift_gamma**(-2))
! ! !     v_ = 0.0
! ! !     w_ = 0.0

! ! !     npart = INT(global_mesh % sx * global_mesh % sy * 0.5 * ppc0)
! ! !     do n = 1, npart
! ! !       xg = random(dseed) * global_mesh % sx
! ! !       ! //AG: exciting a particular mode by displacing particles as amp * sin (2*pi*m*xg / L ) 
! ! !       ! xg = xg + amplitude * sin(2.0 * acos(-1.0) * mode * xg / global_mesh % sx)
! ! !       yg = 0.5
! ! ! #ifdef twoD
! ! !       yg = random(dseed) * global_mesh % sy
! ! ! #endif
! ! !       call injectParticleGlobally(1, xg, yg, 0.5, u_, v_, w_)
! ! !       call injectParticleGlobally(2, xg, yg, 0.5, -u_, v_, w_)
! ! !     end do
! !   end subroutine userInitParticles

!   subroutine userInitFields()
!     implicit none
!     integer :: i, j, k
!     integer :: i_glob, j_glob, k_glob
!     real :: sx_glob, sy_glob
!     real :: kx
!     real :: ex_norm, ey_norm, ez_norm, exyz_norm
!     real :: bx_norm, by_norm, bz_norm, bxyz_norm
    
!     sx_glob = REAL(global_mesh%sx)
!     sy_glob = REAL(global_mesh%sy)
    
!     ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
!     bx(:, :, :) = B_0 * cos(psi); by(:, :, :) = B_0 * sin(psi); bz(:, :, :) = 0

!     kx = mode
!     do i = 0, this_meshblock % ptr % sx - 1
!       i_glob = i + this_meshblock % ptr % x0
!       do j = 0, this_meshblock % ptr % sy - 1
!         j_glob = j + this_meshblock % ptr % y0
!         do k = 0, this_meshblock % ptr % sz - 1
!           k_glob = k + this_meshblock % ptr % z0
!           ex(i, j, k) = -sin(psi) * b_wave(REAL(i_glob), REAL(j_glob-0.5), REAL(k_glob-0.5))!B_amp * sin((i_glob) * kx) !* exp(-0.5*((i_glob) - 3*spread)**2 /spread**2)
!           ey(i, j, k) =  cos(psi) * b_wave(REAL(i_glob-0.5), REAL(j_glob), REAL(k_glob-0.5))!B_amp * sin((i_glob - 0.5) * kx) !* exp(-0.5*((i_glob - 0.5) - 3*spread)**2 /spread**2)
!           ! ez(i, j, k) = 0
!           ! bx(i, j, k) = 0
!           ! by(i, j, k) = 0
!           bz(i, j, k) = b_wave(REAL(i_glob), REAL(j_glob), REAL(k_glob-0.5))!B_amp * sin(i_glob * kx) !* exp(-0.5*(i_glob - 3*spread)**2 /spread**2)
!         end do
!       end do
!     end do

!   end subroutine userInitFields
!   !............................................................!

!   !--- driving ------------------------------------------------!
!   subroutine userCurrentDeposit(step)
!     implicit none
!     integer, optional, intent(in) :: step
!     ! called after particles move and deposit ...
!     ! ... and before the currents are added to the electric field
!   end subroutine userCurrentDeposit

!   subroutine userDriveParticles(step)
!     implicit none
!     integer, optional, intent(in) :: step
!     ! ... dummy loop ...
!     ! integer :: s, ti, tj, tk, p
!     ! do s = 1, nspec
!     !   do ti = 1, species(s)%tile_nx
!     !     do tj = 1, species(s)%tile_ny
!     !       do tk = 1, species(s)%tile_nz
!     !         do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
!     !           ...
!     !         end do
!     !       end do
!     !     end do
!     !   end do
!     ! end do
!   end subroutine userDriveParticles
!   !............................................................!

!   !--- boundaries ---------------------------------------------!
!   subroutine userParticleBoundaryConditions(step)
!     implicit none
!     integer, optional, intent(in) :: step
!   end subroutine userParticleBoundaryConditions

!   subroutine userFieldBoundaryConditions(step, updateE, updateB)
!     implicit none
!     integer :: i, j, k
!     integer :: i_glob, j_glob, k_glob
!     ! real :: kx, ky, kz

!     integer, optional, intent(in) :: step
!     logical, optional, intent(in) :: updateE, updateB
!     logical :: updateE_, updateB_

!     if (present(updateE)) then
!       updateE_ = updateE
!     else
!       updateE_ = .true.
!     end if

!     if (present(updateB)) then
!       updateB_ = updateB
!     else
!       updateB_ = .true.
!     end if
        
!   end subroutine userFieldBoundaryConditions
  
! #include "optional.F"
!   !............................................................!
! end module m_userfile
