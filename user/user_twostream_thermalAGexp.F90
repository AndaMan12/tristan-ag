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
  real, private :: shift_gamma
  ! //AG: Temperature of both electron and positron
  real, private :: TT
  ! //AG: Required to excite a particular mode with a given amplitude
  real, private :: amplitude  
  real, private :: mode
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'shift_gamma', shift_gamma)
    ! //AG: Read temperature from input file (units: m_e c**2)
    call getInput('problem', 'temperature', TT)
    ! //AG: Read which mode to excite and the amplitude from the input file.
    call getInput('problem', 'amplitude', amplitude)    
    call getInput('problem', 'mode', mode)
    
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    
    ! dummy1 -> amplitude of perturbation
    ! dummy2 -> spatial mode of perturbation
    
    userSpatialDistribution = 1 + dummy1 * sin(dummy2 * x_glob)
    return
  end function

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    return
  end function

  subroutine userInitParticles()
    implicit none
    ! //AG: Initialising variables for injecting thermal plasma
    real        :: n_ele, n_posi, n_pro, n_apro ! no. density of electron, positron, proton, antiproton
    type(region):: fill_region                  ! defines a region (on a domain or meshblock) where the plasma will be sprinkled
    real        :: sx_glob
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
    ! *************************************************************************************
    
    n_ele             = 0.25 * ppc0
    n_posi            = 0.25 * ppc0
    n_pro             = 0.25 * ppc0
    n_apro            = 0.25 * ppc0
    
    sx_glob           = REAL(global_mesh%sx)
    fill_region%x_min = REAL(0)
    fill_region%x_max = REAL(sx_glob)

    ! //AG: Parameters of "fillRegionWithThermalPlasma()" for reference
    ! fillRegionWithThermalPlasma(fill_region, fill_species, num_species, ndens_sp, &
    ! temperature, shift_gamma, shift_dir, zero_current, &
    ! dimension, weights, spat_distr_ptr, &
    ! dummy1, dummy2, dummy3)

    ! filling the positrons...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/1/), num_species = 1,&
     ndens_sp = n_posi, temperature = TT, shift_gamma = shift_gamma, shift_dir = +1, dimension = 1,&
     spat_distr_ptr = spat_distr_ptr, dummy1 = amplitude, dummy2 = mode)
    
    ! filling the electrons...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/2/), num_species = 1,&
     ndens_sp = n_ele, temperature = TT, shift_gamma = shift_gamma, shift_dir = -1, dimension = 1,&
     spat_distr_ptr = spat_distr_ptr, dummy1 = amplitude, dummy2 = mode)

     ! filling the antiprotons...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/3/), num_species = 1,&
     ndens_sp = n_apro, temperature = TT, shift_gamma = shift_gamma, shift_dir = +1, dimension = 1,&
     spat_distr_ptr = spat_distr_ptr, dummy1 = amplitude, dummy2 = mode)
    
    ! filling the protonss...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/4/), num_species = 1,&
     ndens_sp = n_pro, temperature = TT, shift_gamma = shift_gamma, shift_dir = -1, dimension = 1,&
     spat_distr_ptr = spat_distr_ptr, dummy1 = amplitude, dummy2 = mode)


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
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0
    ! ... dummy loop ...
    ! do i = 0, this_meshblock%ptr%sx - 1
    !   i_glob = i + this_meshblock%ptr%x0
    !   do j = 0, this_meshblock%ptr%sy - 1
    !     j_glob = j + this_meshblock%ptr%y0
    !     do k = 0, this_meshblock%ptr%sz - 1
    !       k_glob = k + this_meshblock%ptr%z0
    !       ...
    !     end do
    !   end do
    ! end do
  end subroutine userInitFields
  !............................................................!

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

  subroutine userExternalFields(xp, yp, zp, &
                                ex_ext, ey_ext, ez_ext, &
                                bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in) :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    ! some functions of xp, yp, zp
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
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
  !............................................................!

#include "optional.F"
end module m_userfile
