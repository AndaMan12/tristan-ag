#ifdef MOVING_WINDOW
module m_moving_window
  use m_globalnamespace
  use m_errors
  use m_domain
  use m_fields
  use m_particles
  use m_particlelogistics, only: shiftParticlesX
  use m_readinput, only: getInput
  use m_userfile, only: userFillNewRegion
  implicit none

  real(kind=dprec), save :: mw_residual_shift = 0.0d0
  logical, save :: mw_initialized = .false.
  logical, save :: use_moving_window = .false.
  integer, save :: mw_shift_start = 0, mw_shift_interval = 0
  real(kind=dprec), save :: mw_speed = 0.0d0, mw_gamma_param = 0.0d0

  private :: shift_fields, shift_particles, refill_new_region!, initialize_tile, append_to_tile, clear_tile
contains
  subroutine moving_window_step(it)
    implicit none
    integer, intent(in) :: it
    real(kind=dprec) :: shift_accum
    integer :: shift_cells

    if (.not. mw_initialized) call initialize_moving_window_config()

    if (.not. use_moving_window) return
    if (it < mw_shift_start) return

    if (mw_speed <= 0.0d0) return
    if (mw_shift_interval <= 0) return

    mw_residual_shift = mw_residual_shift + mw_speed
    if (modulo(it - mw_shift_start, mw_shift_interval) .ne. 0) return

    shift_accum = mw_residual_shift
    shift_cells = int(floor(shift_accum))
    mw_residual_shift = shift_accum - real(shift_cells, kind=dprec)

    if (shift_cells <= 0) return

    call shift_fields(shift_cells)
    call shift_particles(shift_cells)

    this_meshblock % ptr % x0 = this_meshblock % ptr % x0 + shift_cells
    global_mesh % x0 = global_mesh % x0 + shift_cells

    call refill_new_region(shift_cells)
  end subroutine moving_window_step

  subroutine initialize_moving_window_config()
    implicit none
    integer :: movwin_flag

    call getInput('moving_window', 'movwin', movwin_flag, 0)
    call getInput('moving_window', 'shiftstart', mw_shift_start, 0)
    call getInput('moving_window', 'shiftinterval', mw_shift_interval, 0)
    call getInput('moving_window', 'movwingam', mw_gamma_param, 0.0d0)

    use_moving_window = (movwin_flag == 1)
    if (use_moving_window) then
      if (mw_gamma_param > 10000.0d0) then
        mw_speed = CC
      else if (mw_gamma_param < 1.0d0) then
        mw_speed = CC * max(mw_gamma_param, 0.0d0)
      else
        mw_speed = CC * sqrt(1.0d0 - 1.0d0 / (mw_gamma_param * mw_gamma_param))
      end if
    else
      mw_speed = 0.0d0
    end if

    mw_initialized = .true.
  end subroutine initialize_moving_window_config

  subroutine shift_fields(shift)
    implicit none
    integer, intent(in) :: shift
    integer :: i1, i2, j1, j2, k1, k2, new_start
    integer :: sy, sz
#ifdef MPI
    integer :: ierr, left_rank, right_rank, tag_base
#ifdef MPI08
    type(MPI_STATUS) :: istat
#else
    integer :: istat(MPI_STATUS_SIZE)
#endif
#endif
    real(kind=dprec), allocatable :: recv_ex(:, :, :), recv_ey(:, :, :), recv_ez(:, :, :)
    real(kind=dprec), allocatable :: recv_bx(:, :, :), recv_by(:, :, :), recv_bz(:, :, :)
    real(kind=dprec), allocatable :: recv_jx(:, :, :), recv_jy(:, :, :), recv_jz(:, :, :)
    real(kind=dprec), allocatable :: send_ex(:, :, :), send_ey(:, :, :), send_ez(:, :, :)
    real(kind=dprec), allocatable :: send_bx(:, :, :), send_by(:, :, :), send_bz(:, :, :)
    real(kind=dprec), allocatable :: send_jx(:, :, :), send_jy(:, :, :), send_jz(:, :, :)

    i1 = this_meshblock % ptr % i1
    i2 = this_meshblock % ptr % i2
    j1 = this_meshblock % ptr % j1
    j2 = this_meshblock % ptr % j2
    k1 = this_meshblock % ptr % k1
    k2 = this_meshblock % ptr % k2
    sy = j2 - j1 + 1
    sz = k2 - k1 + 1

    if (shift >= this_meshblock % ptr % sx) then
      call throwError('Requested moving-window shift is larger than the local domain in x.')
    end if

    allocate (recv_ex(shift, sy, sz), recv_ey(shift, sy, sz), recv_ez(shift, sy, sz))
    allocate (recv_bx(shift, sy, sz), recv_by(shift, sy, sz), recv_bz(shift, sy, sz))
    allocate (recv_jx(shift, sy, sz), recv_jy(shift, sy, sz), recv_jz(shift, sy, sz))

    recv_ex = 0.0d0; recv_ey = 0.0d0; recv_ez = 0.0d0
    recv_bx = 0.0d0; recv_by = 0.0d0; recv_bz = 0.0d0
    recv_jx = 0.0d0; recv_jy = 0.0d0; recv_jz = 0.0d0

#ifdef MPI
    tag_base = 900
    left_rank = MPI_PROC_NULL
    right_rank = MPI_PROC_NULL
    if (associated(this_meshblock % ptr % neighbor(-1, 0, 0) % ptr)) &
      left_rank = this_meshblock % ptr % neighbor(-1, 0, 0) % ptr % rnk
    if (associated(this_meshblock % ptr % neighbor(1, 0, 0) % ptr)) &
      right_rank = this_meshblock % ptr % neighbor(1, 0, 0) % ptr % rnk

    allocate (send_ex(shift, sy, sz), send_ey(shift, sy, sz), send_ez(shift, sy, sz))
    allocate (send_bx(shift, sy, sz), send_by(shift, sy, sz), send_bz(shift, sy, sz))
    allocate (send_jx(shift, sy, sz), send_jy(shift, sy, sz), send_jz(shift, sy, sz))

    send_ex = ex(i1:i1 + shift - 1, j1:j2, k1:k2)
    send_ey = ey(i1:i1 + shift - 1, j1:j2, k1:k2)
    send_ez = ez(i1:i1 + shift - 1, j1:j2, k1:k2)
    send_bx = bx(i1:i1 + shift - 1, j1:j2, k1:k2)
    send_by = by(i1:i1 + shift - 1, j1:j2, k1:k2)
    send_bz = bz(i1:i1 + shift - 1, j1:j2, k1:k2)
    send_jx = jx(i1:i1 + shift - 1, j1:j2, k1:k2)
    send_jy = jy(i1:i1 + shift - 1, j1:j2, k1:k2)
    send_jz = jz(i1:i1 + shift - 1, j1:j2, k1:k2)

    call MPI_SENDRECV(send_ex, shift * sy * sz, default_mpi_real, left_rank, tag_base, &
                      recv_ex, shift * sy * sz, default_mpi_real, right_rank, tag_base, MPI_COMM_WORLD, istat, ierr)
    call MPI_SENDRECV(send_ey, shift * sy * sz, default_mpi_real, left_rank, tag_base + 1, &
                      recv_ey, shift * sy * sz, default_mpi_real, right_rank, tag_base + 1, MPI_COMM_WORLD, istat, ierr)
    call MPI_SENDRECV(send_ez, shift * sy * sz, default_mpi_real, left_rank, tag_base + 2, &
                      recv_ez, shift * sy * sz, default_mpi_real, right_rank, tag_base + 2, MPI_COMM_WORLD, istat, ierr)
    call MPI_SENDRECV(send_bx, shift * sy * sz, default_mpi_real, left_rank, tag_base + 3, &
                      recv_bx, shift * sy * sz, default_mpi_real, right_rank, tag_base + 3, MPI_COMM_WORLD, istat, ierr)
    call MPI_SENDRECV(send_by, shift * sy * sz, default_mpi_real, left_rank, tag_base + 4, &
                      recv_by, shift * sy * sz, default_mpi_real, right_rank, tag_base + 4, MPI_COMM_WORLD, istat, ierr)
    call MPI_SENDRECV(send_bz, shift * sy * sz, default_mpi_real, left_rank, tag_base + 5, &
                      recv_bz, shift * sy * sz, default_mpi_real, right_rank, tag_base + 5, MPI_COMM_WORLD, istat, ierr)
    call MPI_SENDRECV(send_jx, shift * sy * sz, default_mpi_real, left_rank, tag_base + 6, &
                      recv_jx, shift * sy * sz, default_mpi_real, right_rank, tag_base + 6, MPI_COMM_WORLD, istat, ierr)
    call MPI_SENDRECV(send_jy, shift * sy * sz, default_mpi_real, left_rank, tag_base + 7, &
                      recv_jy, shift * sy * sz, default_mpi_real, right_rank, tag_base + 7, MPI_COMM_WORLD, istat, ierr)
    call MPI_SENDRECV(send_jz, shift * sy * sz, default_mpi_real, left_rank, tag_base + 8, &
                      recv_jz, shift * sy * sz, default_mpi_real, right_rank, tag_base + 8, MPI_COMM_WORLD, istat, ierr)

    deallocate (send_ex, send_ey, send_ez, send_bx, send_by, send_bz, send_jx, send_jy, send_jz)
#endif

    ex(i1:i2 - shift, j1:j2, k1:k2) = ex(i1 + shift:i2, j1:j2, k1:k2)
    ey(i1:i2 - shift, j1:j2, k1:k2) = ey(i1 + shift:i2, j1:j2, k1:k2)
    ez(i1:i2 - shift, j1:j2, k1:k2) = ez(i1 + shift:i2, j1:j2, k1:k2)
    bx(i1:i2 - shift, j1:j2, k1:k2) = bx(i1 + shift:i2, j1:j2, k1:k2)
    by(i1:i2 - shift, j1:j2, k1:k2) = by(i1 + shift:i2, j1:j2, k1:k2)
    bz(i1:i2 - shift, j1:j2, k1:k2) = bz(i1 + shift:i2, j1:j2, k1:k2)
    jx(i1:i2 - shift, j1:j2, k1:k2) = jx(i1 + shift:i2, j1:j2, k1:k2)
    jy(i1:i2 - shift, j1:j2, k1:k2) = jy(i1 + shift:i2, j1:j2, k1:k2)
    jz(i1:i2 - shift, j1:j2, k1:k2) = jz(i1 + shift:i2, j1:j2, k1:k2)

    new_start = i2 - shift + 1
    if (.not. associated(this_meshblock % ptr % neighbor(1, 0, 0) % ptr)) then
      ex(new_start:i2, j1:j2, k1:k2) = 0.0
      ey(new_start:i2, j1:j2, k1:k2) = 0.0
      ez(new_start:i2, j1:j2, k1:k2) = 0.0
      bx(new_start:i2, j1:j2, k1:k2) = 0.0
      by(new_start:i2, j1:j2, k1:k2) = 0.0
      bz(new_start:i2, j1:j2, k1:k2) = 0.0
      jx(new_start:i2, j1:j2, k1:k2) = 0.0
      jy(new_start:i2, j1:j2, k1:k2) = 0.0
      jz(new_start:i2, j1:j2, k1:k2) = 0.0
    else
      ex(new_start:i2, j1:j2, k1:k2) = recv_ex
      ey(new_start:i2, j1:j2, k1:k2) = recv_ey
      ez(new_start:i2, j1:j2, k1:k2) = recv_ez
      bx(new_start:i2, j1:j2, k1:k2) = recv_bx
      by(new_start:i2, j1:j2, k1:k2) = recv_by
      bz(new_start:i2, j1:j2, k1:k2) = recv_bz
      jx(new_start:i2, j1:j2, k1:k2) = recv_jx
      jy(new_start:i2, j1:j2, k1:k2) = recv_jy
      jz(new_start:i2, j1:j2, k1:k2) = recv_jz
    end if

    deallocate (recv_ex, recv_ey, recv_ez, recv_bx, recv_by, recv_bz, recv_jx, recv_jy, recv_jz)
  end subroutine shift_fields

  subroutine shift_particles(shift)    
    implicit none
    integer, intent(in) :: shift

    if (shift <= 0) return

    ! We are moving the *window* +shift in x, which in code coordinates
    ! is equivalent to shifting *all particles* by -shift in xi.
    call shiftParticlesX(-shift)

    ! Do NOT do any MPI here. We rely on the standard
    ! particle-exchange step that is already present in the code
    ! (the same one that runs after the normal particle push).
  end subroutine shift_particles



  subroutine refill_new_region(shift)
    implicit none
    integer, intent(in) :: shift
    real(kind=dprec) :: xmin_d, xmax_d
    real :: xmin, xmax

    if (.not. associated(this_meshblock % ptr % neighbor(1, 0, 0) % ptr)) then
      xmin_d = REAL(this_meshblock % ptr % x0 + this_meshblock % ptr % sx - shift, kind=dprec)
      xmax_d = REAL(this_meshblock % ptr % x0 + this_meshblock % ptr % sx - 1, kind=dprec)

      xmin = REAL(xmin_d)
      xmax = REAL(xmax_d)

      call userFillNewRegion(xmin, xmax)
    end if
  end subroutine refill_new_region

  ! subroutine initialize_tile(tile, s, ti, tj, tk)
  !   implicit none
  !   type(particle_tile), intent(inout) :: tile
  !   integer, intent(in) :: s, ti, tj, tk
  !   integer :: maxptl_on_tile

  !   maxptl_on_tile = INT(maxptl_array(s) / INT(species(s) % tile_nx * species(s) % tile_ny * species(s) % tile_nz, 8), 4)
  !   maxptl_on_tile = ((maxptl_on_tile + VEC_LEN - 1) / VEC_LEN) * VEC_LEN

  !   tile % spec = s
  !   tile % x1 = (ti - 1) * species(s) % tile_sx
  !   tile % x2 = min(ti * species(s) % tile_sx, this_meshblock % ptr % sx)
  !   tile % y1 = (tj - 1) * species(s) % tile_sy
  !   tile % y2 = min(tj * species(s) % tile_sy, this_meshblock % ptr % sy)
  !   tile % z1 = (tk - 1) * species(s) % tile_sz
  !   tile % z2 = min(tk * species(s) % tile_sz, this_meshblock % ptr % sz)

  !   call allocateParticlesOnEmptyTile(tile, maxptl_on_tile)
  ! end subroutine initialize_tile

!   subroutine append_to_tile(tile, xi, yi, zi, dx, dy, dz, u, v, w, ind, proc, weight)
!     implicit none
!     type(particle_tile), intent(inout) :: tile
!     integer, intent(in) :: ind, proc
!     integer(kind=2), intent(in) :: xi, yi, zi
!     real, intent(in) :: dx, dy, dz, u, v, w, weight
!     integer :: p_new

! #ifdef DEBUG
!     if ((xi < tile % x1) .or. (xi >= tile % x2)) then
!       call throwError('Moving window attempted to place particle outside tile bounds in x.')
!     end if
! #endif

!     if (tile % npart_sp >= tile % maxptl_sp) then
!       if (resize_tiles) then
!         call reallocTileSize(tile, .true.)
!       else
!         call throwError('Moving window ran out of particle capacity when shifting tiles.')
!       end if
!     end if

!     p_new = tile % npart_sp + 1
!     tile % xi(p_new) = xi
!     tile % yi(p_new) = yi
!     tile % zi(p_new) = zi
!     tile % dx(p_new) = dx
!     tile % dy(p_new) = dy
!     tile % dz(p_new) = dz
!     tile % u(p_new) = u
!     tile % v(p_new) = v
!     tile % w(p_new) = w
!     tile % ind(p_new) = ind
!     tile % proc(p_new) = proc
!     tile % weight(p_new) = weight
!     tile % npart_sp = p_new
!   end subroutine append_to_tile

!   subroutine clear_tile(tile)
!     implicit none
!     type(particle_tile), intent(inout) :: tile
!     if (allocated(tile % xi)) deallocate (tile % xi)
!     if (allocated(tile % yi)) deallocate (tile % yi)
!     if (allocated(tile % zi)) deallocate (tile % zi)
!     if (allocated(tile % dx)) deallocate (tile % dx)
!     if (allocated(tile % dy)) deallocate (tile % dy)
!     if (allocated(tile % dz)) deallocate (tile % dz)
!     if (allocated(tile % u)) deallocate (tile % u)
!     if (allocated(tile % v)) deallocate (tile % v)
!     if (allocated(tile % w)) deallocate (tile % w)
!     if (allocated(tile % ind)) deallocate (tile % ind)
!     if (allocated(tile % proc)) deallocate (tile % proc)
!     if (allocated(tile % weight)) deallocate (tile % weight)
! #ifdef PRTLPAYLOADS
!     if (allocated(tile % payload1)) deallocate (tile % payload1)
!     if (allocated(tile % payload2)) deallocate (tile % payload2)
!     if (allocated(tile % payload3)) deallocate (tile % payload3)
! #endif
!   end subroutine clear_tile
end module m_moving_window
#else
module m_moving_window
  implicit none
contains
  subroutine moving_window_step(it)
    implicit none
    integer, intent(in) :: it
  end subroutine moving_window_step
end module m_moving_window
#endif
