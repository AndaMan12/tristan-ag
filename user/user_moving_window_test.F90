module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_helpers
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
#ifdef USROUTPUT
  use m_writeusroutput
#endif
  implicit none

  real :: background_T = 1e-4
  real :: background_n = 1.00
  real :: drift_gamma = 1.00
  real :: field_amplitude = 0.00
  real :: field_wavelength = 64.00
  logical :: zero_current_load = .true.

  private :: userSpatialDistribution, apply_field_profile
contains
  subroutine userReadInput()
    implicit none
    integer :: zero_curr_flag

    call getInput('problem', 'background_n', background_n, 1.00)
    call getInput('problem', 'background_T', background_T, 1e-4)
    call getInput('problem', 'drift_gamma', drift_gamma, 1.00)
    call getInput('problem', 'field_amplitude', field_amplitude, 0.00)
    call getInput('problem', 'field_wavelength', field_wavelength, 64.00)
    call getInput('problem', 'zero_current', zero_curr_flag, 1)

    zero_current_load = (zero_curr_flag == 1)
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    userSpatialDistribution = 1.0
  end function userSpatialDistribution

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    userSLBload = 1.0
  end function userSLBload

  subroutine userInitParticles()
    implicit none
    type(region) :: back_region
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

    back_region%x_min = 0.0
    back_region%x_max = global_mesh%sx
#ifdef twoD
    back_region%y_min = 0.0
    back_region%y_max = global_mesh%sy
#endif
#ifdef threeD
    back_region%y_min = 0.0
    back_region%y_max = global_mesh%sy
    back_region%z_min = 0.0
    back_region%z_max = global_mesh%sz
#endif

    spat_distr_ptr => userSpatialDistribution

    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, background_n, background_T, &
                                     shift_gamma=drift_gamma, shift_dir=1, zero_current=zero_current_load, &
                                     spat_distr_ptr=spat_distr_ptr)
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    ex(:, :, :) = 0.00; ey(:, :, :) = 0.00; ez(:, :, :) = 0.00
    bx(:, :, :) = 0.00; by(:, :, :) = 0.00; bz(:, :, :) = 0.00
    jx(:, :, :) = 0.00; jy(:, :, :) = 0.00; jz(:, :, :) = 0.00

    call apply_field_profile(real(this_meshblock%ptr%x0), real(this_meshblock%ptr%x0 + this_meshblock%ptr%sx - 1))
  end subroutine userInitFields

  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userDriveParticles

  subroutine userExternalFields(xp, yp, zp, &
                                ex_ext, ey_ext, ez_ext, &
                                bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in) :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext

    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields

  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
  end subroutine userFieldBoundaryConditions

  subroutine writeUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine writeUsrRestart

  subroutine readUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine readUsrRestart

#ifdef USROUTPUT
  subroutine userOutput(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userOutput

  logical function userExcludeParticles(s, ti, tj, tk, p)
    implicit none
    integer, intent(in) :: s, ti, tj, tk, p
    userExcludeParticles = .false.
  end function userExcludeParticles
#endif

  subroutine userFillNewRegion(xmin, xmax)
    implicit none
    real, intent(in) :: xmin, xmax
    type(region) :: refill_region

    refill_region%x_min = xmin
    refill_region%x_max = xmax + 1.0
#ifdef twoD
    refill_region%y_min = 0.0
    refill_region%y_max = global_mesh%sy
#endif
#ifdef threeD
    refill_region%y_min = 0.0
    refill_region%y_max = global_mesh%sy
    refill_region%z_min = 0.0
    refill_region%z_max = global_mesh%sz
#endif

    call fillRegionWithThermalPlasma(refill_region, (/1, 2/), 2, background_n, background_T, &
                                     shift_gamma=drift_gamma, shift_dir=1, zero_current=zero_current_load)

    call apply_field_profile(xmin, xmax)
  end subroutine userFillNewRegion

  subroutine apply_field_profile(xmin_glob, xmax_glob)
    implicit none
    real, intent(in) :: xmin_glob, xmax_glob
    integer :: i_glob, j_glob, k_glob
    integer :: i_local, j_local, k_local
    integer :: i_start, i_end
    real :: phase, wavelength_inv

    if (field_amplitude == 0.00) return

    wavelength_inv = 0.00
    if (field_wavelength .ne. 0.00) wavelength_inv = 1.0 / field_wavelength

    i_start = max(int(ceiling(xmin_glob)), this_meshblock%ptr%x0)
    i_end = min(int(floor(xmax_glob)), this_meshblock%ptr%x0 + this_meshblock%ptr%sx - 1)

    if (i_end < i_start) return

    do i_glob = i_start, i_end
      i_local = this_meshblock%ptr%i1 + (i_glob - this_meshblock%ptr%x0)
#ifdef twoD
      do j_glob = 0, this_meshblock%ptr%sy - 1
        j_local = this_meshblock%ptr%j1 + j_glob
        phase = 2.00 * M_PI * (real(i_glob) + 0.50) * wavelength_inv
        ex(i_local, j_local, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
        ey(i_local, j_local, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = field_amplitude * sin(phase)
        ez(i_local, j_local, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
        bx(i_local, j_local, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
        by(i_local, j_local, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
        bz(i_local, j_local, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = field_amplitude * sin(phase)
        ! jx(i_local, j_local, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
        ! jy(i_local, j_local, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
        ! jz(i_local, j_local, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
      end do
#elif defined(threeD)
      do j_glob = 0, this_meshblock%ptr%sy - 1
        j_local = this_meshblock%ptr%j1 + j_glob
        do k_glob = 0, this_meshblock%ptr%sz - 1
          k_local = this_meshblock%ptr%k1 + k_glob
          phase = 2.00 * M_PI * (real(i_glob) + 0.50) * wavelength_inv
          ex(i_local, j_local, k_local) = 0.00
          ey(i_local, j_local, k_local) = field_amplitude * sin(phase)
          ez(i_local, j_local, k_local) = 0.00
          bx(i_local, j_local, k_local) = 0.00
          by(i_local, j_local, k_local) = 0.00
          bz(i_local, j_local, k_local) = field_amplitude * sin(phase)
          ! jx(i_local, j_local, k_local) = 0.00
          ! jy(i_local, j_local, k_local) = 0.00
          ! jz(i_local, j_local, k_local) = 0.00
        end do
      end do
#else
      phase = 2.00 * M_PI * (real(i_glob) + 0.50) * wavelength_inv
      ex(i_local, this_meshblock%ptr%j1:this_meshblock%ptr%j2, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
      ey(i_local, this_meshblock%ptr%j1:this_meshblock%ptr%j2, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = field_amplitude * sin(phase)
      ez(i_local, this_meshblock%ptr%j1:this_meshblock%ptr%j2, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
      bx(i_local, this_meshblock%ptr%j1:this_meshblock%ptr%j2, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
      by(i_local, this_meshblock%ptr%j1:this_meshblock%ptr%j2, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
      bz(i_local, this_meshblock%ptr%j1:this_meshblock%ptr%j2, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = field_amplitude * sin(phase)
      ! jx(i_local, this_meshblock%ptr%j1:this_meshblock%ptr%j2, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
      ! jy(i_local, this_meshblock%ptr%j1:this_meshblock%ptr%j2, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
      ! jz(i_local, this_meshblock%ptr%j1:this_meshblock%ptr%j2, this_meshblock%ptr%k1:this_meshblock%ptr%k2) = 0.00
#endif
    end do
  end subroutine apply_field_profile
    subroutine userDeallocate()
    implicit none
  end subroutine userDeallocate
end module m_userfile
