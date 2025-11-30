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
  real, private :: shift_gamma, beta
  ! //AG: Temperature of both electron and positron
  real, private :: TT
  ! //AG: decay rate of number density
  real, private :: decay_rate
  ! //AG: radius limit of uniform distribution
  real, private :: rad_lim

  ! //AG: Required to excite a particular mode with a given amplitude (inactive)
  ! real, private :: amplitude  
  ! real, private :: mode_x
  ! real, private :: mode_y
  ! real, private :: mode_z
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'shift_gamma', shift_gamma)
    beta = SQRT(1 - 1/shift_gamma**2)
    ! //AG: Read temperature from input file (units: m_e c**2)
    call getInput('problem', 'temperature', TT)
    ! //AG: Read radius limit of density profile from input
    call getInput('problem', 'rad_lim', rad_lim)
    ! //AG: Read decay rate of density profile from input
    call getInput('problem', 'decay_rate', decay_rate)

    ! //AG: Read which mode to excite and the amplitude from the input file (inactive).
    ! call getInput('problem', 'amplitude', amplitude)    
    ! call getInput('problem', 'mode_x', mode_x)
    ! call getInput('problem', 'mode_y', mode_y)
    ! call getInput('problem', 'mode_z', mode_z)
    
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real :: dist_by_rad
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    
    ! dummy1 -> radius of central uniform dist
    ! dummy2 -> rate of exponential decay of density
    dist_by_rad = SQRT((x_glob-dummy3)**2 + (y_glob-dummy3)**2)/dummy1
    !print *, dist_by_rad
    if (dist_by_rad <= 1) then
      userSpatialDistribution = 1 
    else
      userSpatialDistribution = EXP(-dummy2 * (dist_by_rad-1))
    endif

    ! number density profile looks like this:
    !^
    !|
    !|-------------
    !|############|#\
    !|############|##\
    !|############|####\
    !|############|#######\
    !|############|###########\
    !|------------|----------------> axial radius
    ! (uniform) (rad) (exp. decay)

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
    real           :: sx_glob, sy_glob, sz_glob
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
    ! *************************************************************************************
    
    ne = 0.5 * ppc0
    np = 0.5 * ppc0

    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)
    sz_glob = REAL(global_mesh%sz)

    fill_region%x_min = REAL(0)
    fill_region%x_max = REAL(sx_glob)
    fill_region%y_min = REAL(0)
    fill_region%y_max = REAL(sy_glob)
    fill_region%z_min = REAL(0)
    fill_region%z_max = REAL(sz_glob)

    ! //AG: Parameters of "fillRegionWithThermalPlasma()" for reference
    ! fillRegionWithThermalPlasma(fill_region, fill_species, num_species, ndens_sp, &
    ! temperature, shift_gamma, shift_dir, zero_current, &
    ! dimension, weights, spat_distr_ptr, &
    ! dummy1, dummy2, dummy3)

    ! filling the positrons...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/1/), num_species = 1,&
     ndens_sp = np, temperature = TT, shift_gamma = shift_gamma, shift_dir = +3, dimension = 3,&
     spat_distr_ptr = spat_distr_ptr, dummy1 = rad_lim, dummy2 = decay_rate, dummy3 = sx_glob/2.0)
    
    ! filling the electrons...
    call fillRegionWithThermalPlasma(fill_region, fill_species = (/2/), num_species = 1,&
     ndens_sp = ne, temperature = TT, shift_gamma = shift_gamma, shift_dir = -3, dimension = 3,&
     spat_distr_ptr = spat_distr_ptr, dummy1 = rad_lim, dummy2 = decay_rate, dummy3 = sx_glob/2.0)


  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    real    :: b_phi_mag, rad_ax1, rad_ax2, bx_stag, by_stag, x_cen, y_cen
    ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
    bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0

    x_cen = REAL(global_mesh%sx)/2
    y_cen = REAL(global_mesh%sy)/2
     
    print *, x_cen
    ! AG: loop through all cells...
    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      do j = 0, this_meshblock%ptr%sy - 1
        j_glob = j + this_meshblock%ptr%y0
        ! AG: To take staggerring into account, we compute B-field at appropriate location where
        ! Bx and By are assigned.

        rad_ax1 = SQRT((i_glob-0.5-x_cen)**2 + (j_glob-y_cen)**2) ! AG: radial loc. of bx assignment, staggerred along x
        rad_ax2 = SQRT((j_glob-0.5-y_cen)**2 + (i_glob-x_cen)**2) ! AG: radial loc. of by assignment, staggerred along y
        
        if (rad_ax1 <= rad_lim) then
          b_phi_mag = (unit_ch / B_norm) * beta * ppc0 * M_PI * rad_ax1
          bx_stag   = -b_phi_mag * (j_glob-y_cen)/rad_ax1
        else
          b_phi_mag = (unit_ch / B_norm) * beta * ppc0 * M_PI * rad_lim**2 * (1/rad_ax1 + &
          (1 - EXP(-decay_rate * (rad_ax1/rad_lim - 1)**2))/(decay_rate * rad_ax1))
          bx_stag   = -b_phi_mag * (j_glob-y_cen)/rad_ax1
        endif

        if (rad_ax2 <= rad_lim) then
          b_phi_mag = (unit_ch / B_norm) * beta * ppc0 * M_PI * rad_ax2
          by_stag   = b_phi_mag * (i_glob-x_cen)/rad_ax2
        else
          b_phi_mag = (unit_ch / B_norm) * beta * ppc0 * M_PI * rad_lim**2 * (1/rad_ax2 + &
          (1 - EXP(-decay_rate * (rad_ax2/rad_lim - 1)**2))/(decay_rate * rad_ax2))
          by_stag   = b_phi_mag * (i_glob-x_cen)/rad_ax2
        endif

        do k = 0, this_meshblock%ptr%sz - 1
          k_glob = k + this_meshblock%ptr%z0
          ! Setting the azimuthal magnetic field
          bx(i, j, k) = bx_stag
          by(i, j, k) = by_stag
        end do
      end do
    end do
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
!   real, private :: shift_gamma, beta
!   ! //AG: Temperature of both electron and positron
!   real, private :: TT
!   ! //AG: decay rate of number density
!   real, private :: decay_rate
!   ! //AG: radius limit of uniform distribution
!   real, private :: rad_lim

!   ! //AG: Required to excite a particular mode with a given amplitude (inactive)
!   ! real, private :: amplitude  
!   ! real, private :: mode_x
!   ! real, private :: mode_y
!   ! real, private :: mode_z
!   !...............................................................!

!   !--- PRIVATE functions -----------------------------------------!
!   private :: userSpatialDistribution
!   !...............................................................!
! contains
!   !--- initialization -----------------------------------------!
!   subroutine userReadInput()
!     implicit none
!     call getInput('problem', 'shift_gamma', shift_gamma)
!     beta = SQRT(1 - 1/shift_gamma**2)
!     ! //AG: Read temperature from input file (units: m_e c**2)
!     call getInput('problem', 'temperature', TT)
!     ! //AG: Read radius limit of density profile from input
!     call getInput('problem', 'rad_lim', rad_lim)
!     ! //AG: Read decay rate of density profile from input
!     call getInput('problem', 'decay_rate', decay_rate)

!     ! //AG: Read which mode to excite and the amplitude from the input file (inactive).
!     ! call getInput('problem', 'amplitude', amplitude)    
!     ! call getInput('problem', 'mode_x', mode_x)
!     ! call getInput('problem', 'mode_y', mode_y)
!     ! call getInput('problem', 'mode_z', mode_z)
    
!   end subroutine userReadInput

!   function userSpatialDistribution(x_glob, y_glob, z_glob, &
!                                    dummy1, dummy2, dummy3)
!     real :: userSpatialDistribution
!     real :: dist_by_rad
!     real, intent(in), optional :: x_glob, y_glob, z_glob
!     real, intent(in), optional :: dummy1, dummy2, dummy3
    
!     ! dummy1 -> radius of central uniform dist
!     ! dummy2 -> rate of exponential decay of density
!     dist_by_rad = SQRT((x_glob-dummy3)**2 + (y_glob-dummy3)**2)/dummy1
        
!     if (dist_by_rad <= 1) then
!       userSpatialDistribution = 1 
!     else
!       userSpatialDistribution = EXP(-dummy2 * dist_by_rad)
!     endif

!     ! number density profile looks like this:
!     !^
!     !|
!     !|-------------
!     !|############|#\
!     !|############|##\
!     !|############|####\
!     !|############|#######\
!     !|############|###########\
!     !|------------|----------------> axial radius
!     ! (uniform) (rad) (exp. decay)

!     return
!   end function

!   function userSLBload(x_glob, y_glob, z_glob, &
!                        dummy1, dummy2, dummy3)
!     real :: userSLBload
!     ! global coordinates
!     real, intent(in), optional :: x_glob, y_glob, z_glob
!     ! global box dimensions
!     real, intent(in), optional :: dummy1, dummy2, dummy3
!     return
!   end function

!   subroutine userInitParticles()
!     implicit none
!     ! //AG: Initialising variables for injecting thermal plasma
!     real           :: ne, np           ! no. density of electron and positron
!     type(region)   :: fill_region      ! defines a region (on a domain or meshblock) where the plasma will be sprinkled
!     real           :: sx_glob, sy_glob, sz_glob
!     procedure(spatialDistribution), pointer :: spat_distr_ptr => null()
!     spat_distr_ptr => userSpatialDistribution
!     ! *************************************************************************************
    
!     ne = 0.5 * ppc0
!     np = 0.5 * ppc0

!     sx_glob = REAL(global_mesh%sx)
!     sy_glob = REAL(global_mesh%sy)
!     sz_glob = REAL(global_mesh%sz)

!     fill_region%x_min = REAL(0)
!     fill_region%x_max = REAL(sx_glob)
!     fill_region%y_min = REAL(0)
!     fill_region%y_max = REAL(sy_glob)
!     fill_region%z_min = REAL(0)
!     fill_region%z_max = REAL(sz_glob)

!     ! //AG: Parameters of "fillRegionWithThermalPlasma()" for reference
!     ! fillRegionWithThermalPlasma(fill_region, fill_species, num_species, ndens_sp, &
!     ! temperature, shift_gamma, shift_dir, zero_current, &
!     ! dimension, weights, spat_distr_ptr, &
!     ! dummy1, dummy2, dummy3)

!     ! filling the positrons...
!     call fillRegionWithThermalPlasma(fill_region, fill_species = (/1/), num_species = 1,&
!      ndens_sp = np, temperature = TT, shift_gamma = shift_gamma, shift_dir = +3, dimension = 3,&
!      spat_distr_ptr = spat_distr_ptr, dummy1 = rad_lim, dummy2 = decay_rate, dummy3 = sx_glob/2)
    
!     ! filling the electrons...
!     call fillRegionWithThermalPlasma(fill_region, fill_species = (/2/), num_species = 1,&
!      ndens_sp = ne, temperature = TT, shift_gamma = shift_gamma, shift_dir = -3, dimension = 3,&
!      spat_distr_ptr = spat_distr_ptr, dummy1 = rad_lim, dummy2 = decay_rate, dummy3 = sx_glob/2)


!   end subroutine userInitParticles

!   subroutine userInitFields()
!     implicit none
!     integer :: i, j, k
!     integer :: i_glob, j_glob, k_glob
!     real    :: b_phi_mag, rad_ax1, rad_ax2, bx_stag, by_stag, x_cen, y_cen
!     ex(:, :, :) = 0; ey(:, :, :) = 0; ez(:, :, :) = 0
!     bx(:, :, :) = 0; by(:, :, :) = 0; bz(:, :, :) = 0
!     jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0

!     x_cen = this_meshblock%ptr%sx / 2
!     y_cen = this_meshblock%ptr%sy / 2
!     ! AG: loop through all cells...
!     do i = 0, this_meshblock%ptr%sx - 1
!       i_glob = i + this_meshblock%ptr%x0
!       do j = 0, this_meshblock%ptr%sy - 1
!         j_glob = j + this_meshblock%ptr%y0
!         ! AG: To take staggerring into account, we compute B-field at appropriate location where
!         ! Bx and By are assigned.

!         rad_ax1 = SQRT((i_glob-0.5-x_cen)**2 + (j_glob-y_cen)**2) ! AG: radial loc. of bx assignment, staggerred along x
!         rad_ax2 = SQRT((j_glob-0.5-y_cen)**2 + (i_glob-x_cen)**2) ! AG: radial loc. of by assignment, staggerred along y
        
!         if (rad_ax1 <= rad_lim) then
!           b_phi_mag = (unit_ch / B_norm) * beta * ppc0 * M_PI * rad_ax1
!           bx_stag   = -b_phi_mag * (j_glob-y_cen)/rad_ax1
!         else
!           b_phi_mag = (unit_ch / B_norm) * beta * ppc0 * M_PI * rad_lim**2 * (1/rad_ax1 + &
!           (EXP(-decay_rate) - EXP(-decay_rate * (rad_ax1/rad_lim)**2))/(decay_rate * rad_ax1))
!           bx_stag   = -b_phi_mag * (j_glob-y_cen)/rad_ax1
!         endif

!         if (rad_ax2 <= rad_lim) then
!           b_phi_mag = (unit_ch / B_norm) * beta * ppc0 * M_PI * rad_ax2
!           by_stag   = b_phi_mag * (i_glob-x_cen)/rad_ax2
!         else
!           b_phi_mag = (unit_ch / B_norm) * beta * ppc0 * M_PI * rad_lim**2 * (1/rad_ax2 + &
!           (EXP(-decay_rate) - EXP(-decay_rate * (rad_ax2/rad_lim)**2))/(decay_rate * rad_ax2))
!           by_stag   = b_phi_mag * (i_glob-x_cen)/rad_ax2
!         endif

!         do k = 0, this_meshblock%ptr%sz - 1
!           k_glob = k + this_meshblock%ptr%z0
!           ! Setting the azimuthal magnetic field
!           bx(i, j, k) = bx_stag
!           by(i, j, k) = by_stag
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

!   subroutine userExternalFields(xp, yp, zp, &
!                                 ex_ext, ey_ext, ez_ext, &
!                                 bx_ext, by_ext, bz_ext)
!     implicit none
!     real, intent(in) :: xp, yp, zp
!     real, intent(out) :: ex_ext, ey_ext, ez_ext
!     real, intent(out) :: bx_ext, by_ext, bz_ext
!     ! some functions of xp, yp, zp
!     ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
!     bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
!   end subroutine userExternalFields
!   !............................................................!

!   !--- boundaries ---------------------------------------------!
!   subroutine userParticleBoundaryConditions(step)
!     implicit none
!     integer, optional, intent(in) :: step
!   end subroutine userParticleBoundaryConditions

!   subroutine userFieldBoundaryConditions(step, updateE, updateB)
!     implicit none
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
!   !............................................................!

! #include "optional.F"
! end module m_userfile