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
  ! //AG: Strength of background magnetic field (chosen to be along x-axis)
  real, private :: b0
  ! //AG: Required to excite a particular mode with a given amplitude
  real, private :: amplitude  
  real, private :: mode_x
  real, private :: mode_y
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
    ! //AG: Read background magnetic field fro input
    call getInput('problem', 'B_field', b0)
    ! //AG: Read which mode to excite and the amplitude from the input file.
    call getInput('problem', 'amplitude', amplitude)    
    call getInput('problem', 'mode_x', mode_x)
    call getInput('problem', 'mode_y', mode_y)
    
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    
    ! dummy1 -> amplitude of perturbation
    ! dummy2 -> x mode of perturbation
    ! dummy3 -> y mode of perturbation
    
    userSpatialDistribution = (1 + dummy1 * sin(dummy2 * x_glob + dummy2 * y_glob))/(1 + dummy1)
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
    real           :: ne, np           ! no. density of electron and positron
    type(region)   :: fill_region      ! defines a region (on a domain or meshblock) where the plasma will be sprinkled
    real           :: sx_glob, sy_glob
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
    ! *************************************************************************************
    
    ne = 0.5 * ppc0
    np = 0.5 * ppc0
    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)
    fill_region%x_min = REAL(0)
    fill_region%x_max = REAL(sx_glob)
    fill_region%y_min = REAL(0)
    fill_region%y_max = REAL(sy_glob)

    ! //AG: Parameters of "fillRegionWithThermalPlasma()" for reference
    ! fillRegionWithThermalPlasma(fill_region, fill_species, num_species, ndens_sp, &
    ! temperature, shift_gamma, shift_dir, zero_current, &
    ! dimension, weights, spat_distr_ptr, &
    ! dummy1, dummy2, dummy3)

    ! filling the positrons...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/1/), num_species = 1,&
     ndens_sp = np, temperature = TT, shift_gamma = shift_gamma, shift_dir = +1, dimension = 2,&
     spat_distr_ptr = spat_distr_ptr, dummy1 = amplitude, dummy2 = 2 * mode_x * M_PI / sx_glob,&
     dummy3 = 2 * mode_y * M_PI / sy_glob)
    
    ! filling the electrons...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/2/), num_species = 1,&
     ndens_sp = ne, temperature = TT, shift_gamma = shift_gamma, shift_dir = -1, dimension = 2,&
     spat_distr_ptr = spat_distr_ptr, dummy1 = amplitude, dummy2 = 2 * mode_x * M_PI / sx_glob,&
     dummy3 = 2 * mode_y * M_PI / sy_glob)


  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = b0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0
    ! AG: loop through all cells...
    ! do i = 0, this_meshblock%ptr%sx - 1
    !   i_glob = i + this_meshblock%ptr%x0
    !   do j = 0, this_meshblock%ptr%sy - 1
    !     j_glob = j + this_meshblock%ptr%y0
    !     do k = 0, this_meshblock%ptr%sz - 1
    !       k_glob = k + this_meshblock%ptr%z0
    !       ! Setting the constant magnetic field along x-axis
    !       ! No need for staggering, as field is uniform.
          
    !       bx(i, j, k) = b0
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
