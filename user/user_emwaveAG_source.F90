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
  real, private :: B_0, B_amp, duration, freq, shift_gamma, TT, psi, weigh
  real, private :: mult1, mult2, mode, init_x_boundary, fin_x_boundary, ramp_width
  integer, private :: wall_x_location
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'B_0', B_0)
    call getInput('problem', 'psi', psi)
    call getInput('problem', 'B_amplitude', B_amp)
    call getInput('problem', 'duration', duration)
    call getInput('problem', 'frequency', freq)

    mode = 2 * M_PI * freq / (cc * Cos(psi))

    call getInput('problem', 'shift_gamma', shift_gamma)
    call getInput('problem', 'temperature', TT)
    ! call getInput('problem', 'density_decay', density_decay)
    ! call getInput('problem', 'density_scale', density_scale)

    call getInput('problem', 'multiplicity_1', mult1)
    call getInput('problem', 'multiplicity_2', mult2)  
    call getInput('problem', 'ramp_width', ramp_width)    
    ! call getInput('problem', 'weights', weigh)
    call getInput('problem', 'wall_x', wall_x_location, 5)

    init_x_boundary = 5 * M_PI / mode + wall_x_location
    fin_x_boundary = init_x_boundary + ramp_width * 2 * M_PI / mode
    weigh = (mode * sin(psi) * B_norm/unit_ch) / (0.5 * ppc0)

  end subroutine userReadInput

  subroutine bg_density_scaled (mult1, mult2, sx_glob, xg, n_0_local)
    real, intent(in) :: mult1, mult2, sx_glob, xg
    real, intent(out) :: n_0_local
    real :: t
    real :: k
    k = log((mult1-mult2) / (0.01 * mult2))/ (fin_x_boundary - init_x_boundary)
    if (xg <= wall_x_location) then
      n_0_local = 0
    else if (xg <= init_x_boundary) then
      n_0_local = mult1
    else
      n_0_local = mult1 + (mult2 - mult1) * &
            (1.0 - (1.0 + k*(xg - init_x_boundary)) * &
            exp(-k*(xg - init_x_boundary)))
      ! n_0_local = mult1 * exp(- k * (xg - init_x_boundary))
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

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    real           :: sx_glob
    sx_glob = REAL(global_mesh%sx)

    ! userSpatialDistribution = 1 - x_glob/sx_glob * (1 - density_decay)
    call bg_density_scaled(mult1, mult2, sx_glob, x_glob, userSpatialDistribution)
    return
  end function

  subroutine userInitParticles()
    implicit none
    ! //AG: Initialising variables for injecting thermal plasma
    real           :: ne, np           ! no. density of electron and positron
    type(region)   :: fill_region      ! defines a region (on a domain or meshblock) where the plasma will be sprinkled
    real           :: sx_glob, sy_glob, sz_glob
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
    ! *************************************************************************************
    
    ne = 0.5 * ppc0
    np = 0.5 * ppc0
    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)
    !sz_glob = REAL(global_mesh%sz)
    fill_region%x_min = REAL(0)
    fill_region%x_max = REAL(sx_glob)
    fill_region%y_min = REAL(0)
    fill_region%y_max = REAL(sy_glob)
    !fill_region%z_min = REAL(0)
    !fill_region%z_max = REAL(sz_glob)

    ! //AG: Parameters of "fillRegionWithThermalPlasma()" for reference
    ! fillRegionWithThermalPlasma(fill_region, fill_species, num_species, ndens_sp, &
    ! temperature, shift_gamma, shift_dir, zero_current, &
    ! dimension, weights, spat_distr_ptr, &
    ! dummy1, dummy2, dummy3)

    ! filling the positrons...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/1/), num_species = 1,&
     ndens_sp = np, temperature = TT, shift_gamma = shift_gamma, shift_dir = +1, dimension = 2, &
     weights = weigh, spat_distr_ptr = spat_distr_ptr)
    
    ! filling the electrons...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/2/), num_species = 1,&
     ndens_sp = ne, temperature = TT, shift_gamma = shift_gamma, shift_dir = -1, dimension = 2, &
     weights = weigh, spat_distr_ptr = spat_distr_ptr)


    ! //AG: Old code to inject particles one by one with a random number generator....

!     u_ = shift_gamma * sqrt(1.0 - shift_gamma**(-2))
!     v_ = 0.0
!     w_ = 0.0

!     npart = INT(global_mesh % sx * global_mesh % sy * 0.5 * ppc0)
!     do n = 1, npart
!       xg = random(dseed) * global_mesh % sx
!       ! //AG: exciting a particular mode by displacing particles as amp * sin (2*pi*m*xg / L ) 
!       ! xg = xg + amplitude * sin(2.0 * acos(-1.0) * mode * xg / global_mesh % sx)
!       yg = 0.5
! #ifdef twoD
!       yg = random(dseed) * global_mesh % sy
! #endif
!       call injectParticleGlobally(1, xg, yg, 0.5, u_, v_, w_)
!       call injectParticleGlobally(2, xg, yg, 0.5, -u_, v_, w_)
!     end do
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    ! integer :: i, j, k
    ! integer :: i_glob, j_glob, k_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = B_0 * cos(psi); by(:, :, :) = B_0 * sin(psi); bz(:, :, :) = 0

  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
    integer :: s, ti, tj, tk, p
    real :: u_, v_, w_, x_l, y_l, z_l, x_g, inv_gamma, gamma
    real :: tfrac, xcolis, xnew, x0
    type(maxwellian) :: mwl
    mwl % dimension = 2
    mwl % generated = .false.
    mwl % temperature = TT
    mwl % shift_gamma = shift_gamma
    mwl % shift_flag = .true.
    mwl % shift_dir = -1

    do s = 1, nspec
     do ti = 1, species(s) % tile_nx
       do tj = 1, species(s) % tile_ny
         do tk = 1, species(s) % tile_nz
           do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
             x_l = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p)) + species(s) % prtl_tile(ti, tj, tk) % dx(p)
             x_g = x_l + REAL(this_meshblock % ptr % x0)
             if (x_g .lt. wall_x_location) then
               u_ = species(s) % prtl_tile(ti, tj, tk) % u(p)
               v_ = species(s) % prtl_tile(ti, tj, tk) % v(p)
               w_ = species(s) % prtl_tile(ti, tj, tk) % w(p)
               gamma = sqrt(1.0 + u_**2 + v_**2 + w_**2)
               inv_gamma = 1.0 / gamma
               x0 = x_l - u_ * inv_gamma * CC
               tfrac = abs((x0 - wall_x_location) / (u_ * inv_gamma * CC))
               if (tfrac .lt. 1) then
                 xcolis = x0 + u_ * inv_gamma * CC * tfrac
               end if
               y_l = REAL(species(s) % prtl_tile(ti, tj, tk) % yi(p)) + species(s) % prtl_tile(ti, tj, tk) % dy(p)
               z_l = REAL(species(s) % prtl_tile(ti, tj, tk) % zi(p)) + species(s) % prtl_tile(ti, tj, tk) % dz(p)
               call depositCurrentsFromSingleParticle(s, species(s) % prtl_tile(ti, tj, tk), p, &
                                                      xcolis, y_l, z_l, x_l, y_l, z_l, -1.0)
               ! reflecting particle
              !  call generateFromMaxwellian(mwl, u_, v_, w_)
               u_ = 0; v_ = 0; w_ = 0
               species(s) % prtl_tile(ti, tj, tk) % u(p) = u_
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
    if (step .le. 3.5*duration) then ! AG: Initial ramp up
      do i = 0, this_meshblock % ptr % sx - 1
        i_glob = i + this_meshblock % ptr % x0
        if (i_glob <= wall_x_location) then
          ey(i, :, :) = B_amp * sin(step * freq) * exp(-0.5*(step -0.5*cc - 3.5*duration)**2 /duration**2)
          bz(i, :, :) = B_amp * sin(step * freq) * exp(-0.5*(step - 3.5*duration)**2 /duration**2) 
          !ez(i, :, :) = B_amp * exp(-0.5*(step -0.5*cc - 3*duration)**2 /duration**2) * cos(step * freq)
          !by(i, :, :) = B_amp * exp(-0.5*(step - 3*duration)**2 /duration**2) * cos(step * freq)

        end if
      end do
    else ! AG: After the initial ramping up is done, a plane wave.
      do i = 0, this_meshblock % ptr % sx - 1
        i_glob = i + this_meshblock % ptr % x0
        if (i_glob <= wall_x_location) then
          ey(i, :, :) = B_amp * sin(step * freq)
          bz(i, :, :) = B_amp * sin(step * freq) !* exp(-0.5*(step - 3.5*duration)**2 /duration**2) 
          !ez(i, :, :) = B_amp * exp(-0.5*(step -0.5*cc - 3*duration)**2 /duration**2) * cos(step * freq)
          !by(i, :, :) = B_amp * exp(-0.5*(step - 3*duration)**2 /duration**2) * cos(step * freq)

        end if
      end do      
    end if
  end subroutine userFieldBoundaryConditions
  
#include "optional.F"
!............................................................!
end module m_userfile
