#ifdef MOVING_WINDOW
module m_moving_window
  use m_globalnamespace
  use m_errors
  use m_domain
  use m_fields
  use m_particles
  use m_particlelogistics
  use m_readinput, only: getInput
  use m_userfile, only: userFillNewRegion
  implicit none

  real(kind=dprec), save :: mw_residual_shift = 0.0d0
  logical, save :: mw_initialized = .false.
  logical, save :: use_moving_window = .false.
  integer, save :: mw_shift_start = 0, mw_shift_interval = 0
  real(kind=dprec), save :: mw_speed = 0.0d0, mw_gamma_param = 0.0d0

  private :: shift_fields, shift_particles, refill_new_region, initialize_tile, append_to_tile, clear_tile
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
    integer :: s, p
    integer :: ti, tj, tk
    integer :: left_rank, right_rank
#ifdef MPI
    integer :: ierr
    integer :: send_cnt_left, send_cnt_right
    integer :: recv_cnt_left, recv_cnt_right
#ifdef MPI08
    type(MPI_STATUS) :: istat
#else
    integer :: istat(MPI_STATUS_SIZE)
#endif
#endif
    integer(kind=8) :: global_min_new, global_max_new
    integer(kind=8) :: x0_old_local, x0_new_local
    integer(kind=8) :: x0_new_left, x0_new_right
    integer(kind=8) :: global_x
    integer(kind=2) :: xi_new
    type(prtl_enroute), allocatable :: send_left(:), send_right(:)
    type(prtl_enroute), allocatable :: recv_left(:), recv_right(:)

    if (shift <= 0) return

    !-----------------------------------------------------------------
    ! Define new global window and new local origins (still using OLD x0)
    !-----------------------------------------------------------------
    global_min_new = INT(global_mesh % x0 + shift, kind=8)
    global_max_new = global_min_new + INT(global_mesh % sx - 1, kind=8)

    x0_old_local = INT(this_meshblock % ptr % x0, kind=8)
    x0_new_local = x0_old_local + shift

    x0_new_left  = 0_8
    x0_new_right = 0_8

#ifdef MPI
    left_rank  = MPI_PROC_NULL
    right_rank = MPI_PROC_NULL
    if (associated(this_meshblock % ptr % neighbor(-1, 0, 0) % ptr)) then
      left_rank  = this_meshblock % ptr % neighbor(-1, 0, 0) % ptr % rnk
      x0_new_left = INT(this_meshblock % ptr % neighbor(-1, 0, 0) % ptr % x0 + shift, kind=8)
    end if
    if (associated(this_meshblock % ptr % neighbor(1, 0, 0) % ptr)) then
      right_rank  = this_meshblock % ptr % neighbor(1, 0, 0) % ptr % rnk
      x0_new_right = INT(this_meshblock % ptr % neighbor(1, 0, 0) % ptr % x0 + shift, kind=8)
    end if
#else
    left_rank  = -1
    right_rank = -1
#endif

    !-----------------------------------------------------------------
    ! 1) Backup all local particles into prtl_backup(:)
    !-----------------------------------------------------------------
    call backupParticles()

    !-----------------------------------------------------------------
    ! 2) Clear local tiles (we will repopulate them from the backup)
    !-----------------------------------------------------------------
    do s = 1, nspec
      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            species(s) % prtl_tile(ti, tj, tk) % npart_sp = 0
          end do
        end do
      end do
    end do

    !-----------------------------------------------------------------
    ! 3) For each species, redistribute backed-up particles:
    !    - drop those outside the new global window
    !    - keep local ones
    !    - send others one rank left/right
    !-----------------------------------------------------------------
    do s = 1, nspec

#ifdef MPI
      if (prtl_backup(s) % cnt > 0) then
        allocate (send_left(prtl_backup(s) % cnt))
        allocate (send_right(prtl_backup(s) % cnt))
      else
        allocate (send_left(0))
        allocate (send_right(0))
      end if
      send_cnt_left  = 0
      send_cnt_right = 0
#endif

      do p = 1, prtl_backup(s) % cnt

        ! global x position BEFORE shift (integer cell index)
        global_x = x0_old_local + INT(prtl_backup(s) % enroute(p) % xi, kind=8)

        ! Particle leaves the global moving window?
        if ((global_x < global_min_new) .or. (global_x > global_max_new)) cycle

#ifdef MPI
        ! Does this particle now belong to the LEFT neighbour?
        if ((global_x < x0_new_local) .and. (left_rank .ne. MPI_PROC_NULL)) then
          send_cnt_left = send_cnt_left + 1
          send_left(send_cnt_left) = prtl_backup(s) % enroute(p)
          ! convert xi to the NEW local coordinate on the left rank
          send_left(send_cnt_left) % xi = INT(global_x - x0_new_left, kind=2)
          cycle
        end if

        ! Or to the RIGHT neighbour?
        if ((global_x >= x0_new_local + INT(this_meshblock % ptr % sx, kind=8)) .and. &
            (right_rank .ne. MPI_PROC_NULL)) then
          send_cnt_right = send_cnt_right + 1
          send_right(send_cnt_right) = prtl_backup(s) % enroute(p)
          send_right(send_cnt_right) % xi = INT(global_x - x0_new_right, kind=2)
          cycle
        end if
#endif

        ! Otherwise, the particle stays on THIS rank.
        if ((global_x >= x0_new_local) .and. &
            (global_x < x0_new_local + INT(this_meshblock % ptr % sx, kind=8))) then

          xi_new = INT(global_x - x0_new_local, kind=2)

          call createParticleFromAttributes( s, xi_new, &
                                             prtl_backup(s) % enroute(p) % yi, &
                                             prtl_backup(s) % enroute(p) % zi, &
                                             prtl_backup(s) % enroute(p) % dx, &
                                             prtl_backup(s) % enroute(p) % dy, &
                                             prtl_backup(s) % enroute(p) % dz, &
                                             prtl_backup(s) % enroute(p) % u,  &
                                             prtl_backup(s) % enroute(p) % v,  &
                                             prtl_backup(s) % enroute(p) % w,  &
#ifdef PRTLPAYLOADS
                                             prtl_backup(s) % enroute(p) % payload1, &
                                             prtl_backup(s) % enroute(p) % payload2, &
                                             prtl_backup(s) % enroute(p) % payload3, &
#endif
#ifdef DEBUG
                                             '`shift_particles backup`', &
#endif
                                             prtl_backup(s) % enroute(p) % ind, &
                                             prtl_backup(s) % enroute(p) % proc, &
                                             prtl_backup(s) % enroute(p) % weight )
        end if

      end do  ! p loop over prtl_backup(s)

#ifdef MPI
      !---------------------------------------------------------------
      ! 4) Exchange particles with left/right neighbours
      !---------------------------------------------------------------
      recv_cnt_right = 0
      recv_cnt_left  = 0

      call MPI_SENDRECV(send_cnt_left,  1, MPI_INTEGER, left_rank,  930, &
                        recv_cnt_right, 1, MPI_INTEGER, right_rank, 930, &
                        MPI_COMM_WORLD, istat, ierr)
      call MPI_SENDRECV(send_cnt_right, 1, MPI_INTEGER, right_rank, 932, &
                        recv_cnt_left,  1, MPI_INTEGER, left_rank,  932, &
                        MPI_COMM_WORLD, istat, ierr)

      if (recv_cnt_right > 0) then
        allocate (recv_right(recv_cnt_right))
      else
        allocate (recv_right(0))
      end if
      if (recv_cnt_left > 0) then
        allocate (recv_left(recv_cnt_left))
      else
        allocate (recv_left(0))
      end if

      call MPI_SENDRECV(send_left,  send_cnt_left,  myMPI_ENROUTE, left_rank,  931, &
                        recv_right, recv_cnt_right, myMPI_ENROUTE, right_rank, 931, &
                        MPI_COMM_WORLD, istat, ierr)
      call MPI_SENDRECV(send_right, send_cnt_right, myMPI_ENROUTE, right_rank, 933, &
                        recv_left,  recv_cnt_left,  myMPI_ENROUTE, left_rank,  933, &
                        MPI_COMM_WORLD, istat, ierr)

      ! Create all particles received from the RIGHT
      do p = 1, recv_cnt_right
        call createParticleFromAttributes( s, recv_right(p) % xi, recv_right(p) % yi, recv_right(p) % zi, &
                                           recv_right(p) % dx, recv_right(p) % dy, recv_right(p) % dz, &
                                           recv_right(p) % u,  recv_right(p) % v,  recv_right(p) % w,  &
#ifdef PRTLPAYLOADS
                                           recv_right(p) % payload1, &
                                           recv_right(p) % payload2, &
                                           recv_right(p) % payload3, &
#endif
#ifdef DEBUG
                                           '`shift_particles recv_right`', &
#endif
                                           recv_right(p) % ind, recv_right(p) % proc, recv_right(p) % weight )
      end do

      ! ...and from the LEFT
      do p = 1, recv_cnt_left
        call createParticleFromAttributes( s, recv_left(p) % xi, recv_left(p) % yi, recv_left(p) % zi, &
                                           recv_left(p) % dx, recv_left(p) % dy, recv_left(p) % dz, &
                                           recv_left(p) % u,  recv_left(p) % v,  recv_left(p) % w,  &
#ifdef PRTLPAYLOADS
                                           recv_left(p) % payload1, &
                                           recv_left(p) % payload2, &
                                           recv_left(p) % payload3, &
#endif
#ifdef DEBUG
                                           '`shift_particles recv_left`', &
#endif
                                           recv_left(p) % ind, recv_left(p) % proc, recv_left(p) % weight )
      end do

      if (allocated(send_left))  deallocate (send_left)
      if (allocated(send_right)) deallocate (send_right)
      if (allocated(recv_left))  deallocate (recv_left)
      if (allocated(recv_right)) deallocate (recv_right)
#endif  ! MPI

    end do   ! species loop

    ! prtl_backup can be freed later by deallocateParticleBackup() if desired

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

  subroutine initialize_tile(tile, s, ti, tj, tk)
    implicit none
    type(particle_tile), intent(inout) :: tile
    integer, intent(in) :: s, ti, tj, tk
    integer :: maxptl_on_tile

    maxptl_on_tile = INT(maxptl_array(s) / INT(species(s) % tile_nx * species(s) % tile_ny * species(s) % tile_nz, 8), 4)
    maxptl_on_tile = ((maxptl_on_tile + VEC_LEN - 1) / VEC_LEN) * VEC_LEN

    tile % spec = s
    tile % x1 = (ti - 1) * species(s) % tile_sx
    tile % x2 = min(ti * species(s) % tile_sx, this_meshblock % ptr % sx)
    tile % y1 = (tj - 1) * species(s) % tile_sy
    tile % y2 = min(tj * species(s) % tile_sy, this_meshblock % ptr % sy)
    tile % z1 = (tk - 1) * species(s) % tile_sz
    tile % z2 = min(tk * species(s) % tile_sz, this_meshblock % ptr % sz)

    call allocateParticlesOnEmptyTile(tile, maxptl_on_tile)
  end subroutine initialize_tile

  subroutine append_to_tile(tile, xi, yi, zi, dx, dy, dz, u, v, w, ind, proc, weight)
    implicit none
    type(particle_tile), intent(inout) :: tile
    integer, intent(in) :: ind, proc
    integer(kind=2), intent(in) :: xi, yi, zi
    real, intent(in) :: dx, dy, dz, u, v, w, weight
    integer :: p_new

#ifdef DEBUG
    if ((xi < tile % x1) .or. (xi >= tile % x2)) then
      call throwError('Moving window attempted to place particle outside tile bounds in x.')
    end if
#endif

    if (tile % npart_sp >= tile % maxptl_sp) then
      if (resize_tiles) then
        call reallocTileSize(tile, .true.)
      else
        call throwError('Moving window ran out of particle capacity when shifting tiles.')
      end if
    end if

    p_new = tile % npart_sp + 1
    tile % xi(p_new) = xi
    tile % yi(p_new) = yi
    tile % zi(p_new) = zi
    tile % dx(p_new) = dx
    tile % dy(p_new) = dy
    tile % dz(p_new) = dz
    tile % u(p_new) = u
    tile % v(p_new) = v
    tile % w(p_new) = w
    tile % ind(p_new) = ind
    tile % proc(p_new) = proc
    tile % weight(p_new) = weight
    tile % npart_sp = p_new
  end subroutine append_to_tile

  subroutine clear_tile(tile)
    implicit none
    type(particle_tile), intent(inout) :: tile
    if (allocated(tile % xi)) deallocate (tile % xi)
    if (allocated(tile % yi)) deallocate (tile % yi)
    if (allocated(tile % zi)) deallocate (tile % zi)
    if (allocated(tile % dx)) deallocate (tile % dx)
    if (allocated(tile % dy)) deallocate (tile % dy)
    if (allocated(tile % dz)) deallocate (tile % dz)
    if (allocated(tile % u)) deallocate (tile % u)
    if (allocated(tile % v)) deallocate (tile % v)
    if (allocated(tile % w)) deallocate (tile % w)
    if (allocated(tile % ind)) deallocate (tile % ind)
    if (allocated(tile % proc)) deallocate (tile % proc)
    if (allocated(tile % weight)) deallocate (tile % weight)
#ifdef PRTLPAYLOADS
    if (allocated(tile % payload1)) deallocate (tile % payload1)
    if (allocated(tile % payload2)) deallocate (tile % payload2)
    if (allocated(tile % payload3)) deallocate (tile % payload3)
#endif
  end subroutine clear_tile
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
