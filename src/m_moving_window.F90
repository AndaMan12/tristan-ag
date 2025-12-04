#ifdef MOVING_WINDOW
module m_moving_window
  use m_globalnamespace
  use m_errors
  use m_domain
  use m_fields
  use m_particles
  use m_particlelogistics, only: reallocTileSize, allocateParticlesOnEmptyTile, createParticleFromAttributes
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
    integer :: s, ti, tj, tk, p
    integer :: ti_new, tj_new, tk_new
    integer(kind=2) :: xi_new, yi_new, zi_new
    type(particle_tile), allocatable :: old_tiles(:, :, :)
    type(particle_tile), allocatable :: new_tiles(:, :, :)
    integer :: total_npart
#ifdef MPI
    integer :: ierr, left_rank, right_rank, left_sx
    integer :: recv_cnt_right, recv_cnt_left, send_cnt_left, send_cnt_right
#ifdef MPI08
    type(MPI_STATUS) :: istat
#else
    integer :: istat(MPI_STATUS_SIZE)
#endif
#endif
    type(prtl_enroute), allocatable :: send_left(:), recv_right(:)
    type(prtl_enroute), allocatable :: send_right(:), recv_left(:)

    if (shift <= 0) return

    do s = 1, nspec
      call move_alloc(species(s) % prtl_tile, old_tiles)

      allocate (new_tiles(species(s) % tile_nx, species(s) % tile_ny, species(s) % tile_nz))
      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            call initialize_tile(new_tiles(ti, tj, tk), s, ti, tj, tk)
          end do
        end do
      end do

      call move_alloc(new_tiles, species(s) % prtl_tile)

      species(s) % cntr_sp = 0
      total_npart = 0

#ifdef MPI
      left_rank = MPI_PROC_NULL
      right_rank = MPI_PROC_NULL
      left_sx = 0
      if (associated(this_meshblock % ptr % neighbor(-1, 0, 0) % ptr)) then
        left_rank = this_meshblock % ptr % neighbor(-1, 0, 0) % ptr % rnk
        left_sx = this_meshblock % ptr % neighbor(-1, 0, 0) % ptr % sx
      end if
      if (associated(this_meshblock % ptr % neighbor(1, 0, 0) % ptr)) then
        right_rank = this_meshblock % ptr % neighbor(1, 0, 0) % ptr % rnk
      end if

      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            total_npart = total_npart + old_tiles(ti, tj, tk) % npart_sp
          end do
        end do
      end do

      if (total_npart > 0) then
        allocate (send_left(total_npart))
        allocate (send_right(total_npart))
      else
        allocate (send_left(0))
        allocate (send_right(0))
      end if
      send_cnt_left = 0
      send_cnt_right = 0
#endif

      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            do p = 1, old_tiles(ti, tj, tk) % npart_sp
              xi_new = old_tiles(ti, tj, tk) % xi(p) - shift
              yi_new = old_tiles(ti, tj, tk) % yi(p)
              zi_new = old_tiles(ti, tj, tk) % zi(p)

#ifdef MPI
              if (xi_new < 0) then
                if (left_rank .ne. MPI_PROC_NULL) then
                  send_cnt_left = send_cnt_left + 1
                  send_left(send_cnt_left) % xi = xi_new + INT(left_sx, kind=2)
                  send_left(send_cnt_left) % yi = yi_new
                  send_left(send_cnt_left) % zi = zi_new
                  send_left(send_cnt_left) % dx = old_tiles(ti, tj, tk) % dx(p)
                  send_left(send_cnt_left) % dy = old_tiles(ti, tj, tk) % dy(p)
                  send_left(send_cnt_left) % dz = old_tiles(ti, tj, tk) % dz(p)
                  send_left(send_cnt_left) % u = old_tiles(ti, tj, tk) % u(p)
                  send_left(send_cnt_left) % v = old_tiles(ti, tj, tk) % v(p)
                  send_left(send_cnt_left) % w = old_tiles(ti, tj, tk) % w(p)
                  send_left(send_cnt_left) % ind = old_tiles(ti, tj, tk) % ind(p)
                  send_left(send_cnt_left) % proc = old_tiles(ti, tj, tk) % proc(p)
                  send_left(send_cnt_left) % weight = old_tiles(ti, tj, tk) % weight(p)
#ifdef PRTLPAYLOADS
                  send_left(send_cnt_left) % payload1 = old_tiles(ti, tj, tk) % payload1(p)
                  send_left(send_cnt_left) % payload2 = old_tiles(ti, tj, tk) % payload2(p)
                  send_left(send_cnt_left) % payload3 = old_tiles(ti, tj, tk) % payload3(p)
#endif
                end if
                cycle
              else if (xi_new >= this_meshblock % ptr % sx) then
                if (right_rank .ne. MPI_PROC_NULL) then
                  send_cnt_right = send_cnt_right + 1
                  send_right(send_cnt_right) % xi = xi_new - INT(this_meshblock % ptr % sx, kind=2)
                  send_right(send_cnt_right) % yi = yi_new
                  send_right(send_cnt_right) % zi = zi_new
                  send_right(send_cnt_right) % dx = old_tiles(ti, tj, tk) % dx(p)
                  send_right(send_cnt_right) % dy = old_tiles(ti, tj, tk) % dy(p)
                  send_right(send_cnt_right) % dz = old_tiles(ti, tj, tk) % dz(p)
                  send_right(send_cnt_right) % u = old_tiles(ti, tj, tk) % u(p)
                  send_right(send_cnt_right) % v = old_tiles(ti, tj, tk) % v(p)
                  send_right(send_cnt_right) % w = old_tiles(ti, tj, tk) % w(p)
                  send_right(send_cnt_right) % ind = old_tiles(ti, tj, tk) % ind(p)
                  send_right(send_cnt_right) % proc = old_tiles(ti, tj, tk) % proc(p)
                  send_right(send_cnt_right) % weight = old_tiles(ti, tj, tk) % weight(p)
#ifdef PRTLPAYLOADS
                  send_right(send_cnt_right) % payload1 = old_tiles(ti, tj, tk) % payload1(p)
                  send_right(send_cnt_right) % payload2 = old_tiles(ti, tj, tk) % payload2(p)
                  send_right(send_cnt_right) % payload3 = old_tiles(ti, tj, tk) % payload3(p)
#endif
                end if
                cycle
              end if
#endif

              if ((xi_new < 0) .or. (xi_new >= this_meshblock % ptr % sx)) cycle

              ti_new = FLOOR(REAL(xi_new) / REAL(species(s) % tile_sx)) + 1
              tj_new = FLOOR(REAL(yi_new) / REAL(species(s) % tile_sy)) + 1
              tk_new = FLOOR(REAL(zi_new) / REAL(species(s) % tile_sz)) + 1

              call createParticleFromAttributes(s, xi_new, yi_new, zi_new, &
                                               old_tiles(ti, tj, tk) % dx(p), old_tiles(ti, tj, tk) % dy(p), &
                                               old_tiles(ti, tj, tk) % dz(p), &
                                               old_tiles(ti, tj, tk) % u(p), old_tiles(ti, tj, tk) % v(p), &
                                               old_tiles(ti, tj, tk) % w(p), &
#ifdef DEBUG
                                               '`shift_particles`', &
#endif
                                               old_tiles(ti, tj, tk) % ind(p), mpi_rank, old_tiles(ti, tj, tk) % weight(p))

              species(s) % cntr_sp = species(s) % cntr_sp + 1
            end do

            call clear_tile(old_tiles(ti, tj, tk))
          end do
        end do
      end do

#ifdef MPI
      recv_cnt_right = 0
      recv_cnt_left = 0

      call MPI_SENDRECV(send_cnt_left, 1, MPI_INTEGER, left_rank, 920, &
                        recv_cnt_right, 1, MPI_INTEGER, right_rank, 920, MPI_COMM_WORLD, istat, ierr)
      call MPI_SENDRECV(send_cnt_right, 1, MPI_INTEGER, right_rank, 922, &
                        recv_cnt_left, 1, MPI_INTEGER, left_rank, 922, MPI_COMM_WORLD, istat, ierr)

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

      call MPI_SENDRECV(send_left, send_cnt_left, myMPI_ENROUTE, left_rank, 921, &
                        recv_right, recv_cnt_right, myMPI_ENROUTE, right_rank, 921, MPI_COMM_WORLD, istat, ierr)
      call MPI_SENDRECV(send_right, send_cnt_right, myMPI_ENROUTE, right_rank, 923, &
                        recv_left, recv_cnt_left, myMPI_ENROUTE, left_rank, 923, MPI_COMM_WORLD, istat, ierr)

      do p = 1, recv_cnt_right
        call createParticleFromAttributes(s, recv_right(p) % xi, recv_right(p) % yi, recv_right(p) % zi, &
                                         recv_right(p) % dx, recv_right(p) % dy, recv_right(p) % dz, &
                                         recv_right(p) % u, recv_right(p) % v, recv_right(p) % w, &
#ifdef DEBUG
                                         '`shift_particles recv_right`', &
#endif
                                         recv_right(p) % ind, recv_right(p) % proc, recv_right(p) % weight)
        species(s) % cntr_sp = species(s) % cntr_sp + 1
      end do

      do p = 1, recv_cnt_left
        xi_new = recv_left(p) % xi + INT(this_meshblock % ptr % sx, kind=2)
        call createParticleFromAttributes(s, xi_new, recv_left(p) % yi, recv_left(p) % zi, &
                                         recv_left(p) % dx, recv_left(p) % dy, recv_left(p) % dz, &
                                         recv_left(p) % u, recv_left(p) % v, recv_left(p) % w, &
#ifdef DEBUG
                                         '`shift_particles recv_left`', &
#endif
                                         recv_left(p) % ind, recv_left(p) % proc, recv_left(p) % weight)
        species(s) % cntr_sp = species(s) % cntr_sp + 1
      end do

      if (allocated(send_left)) deallocate (send_left)
      if (allocated(send_right)) deallocate (send_right)
      if (allocated(recv_right)) deallocate (recv_right)
      if (allocated(recv_left)) deallocate (recv_left)
#endif

      if (allocated(old_tiles)) then
        deallocate (old_tiles)
      end if
    end do
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
