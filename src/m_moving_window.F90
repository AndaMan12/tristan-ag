module m_moving_window
  use m_globalnamespace
  use m_errors
  use m_domain
  use m_fields
  use m_particles
  use m_particlelogistics, only: reallocTileSize, allocateParticlesOnEmptyTile
  use m_userfile, only: use_moving_window, mw_shift_start, mw_shift_interval, mw_ncells, userFillNewRegion
  implicit none

  private :: shift_fields, shift_particles, refill_new_region, initialize_tile, append_to_tile, clear_tile
contains
  subroutine moving_window_step(it)
    implicit none
    integer, intent(in) :: it

    if (.not. use_moving_window) return
    if (mw_ncells <= 0) return
    if (mw_shift_interval <= 0) return
    if (it < mw_shift_start) return
    if (modulo(it - mw_shift_start, mw_shift_interval) .ne. 0) return

    call shift_fields(mw_ncells)
    call shift_particles(mw_ncells)

    this_meshblock % ptr % x0 = this_meshblock % ptr % x0 + mw_ncells
    global_mesh % x0 = global_mesh % x0 + mw_ncells

    call refill_new_region(mw_ncells)
  end subroutine moving_window_step

  subroutine shift_fields(shift)
    implicit none
    integer, intent(in) :: shift
    integer :: i1, i2, j1, j2, k1, k2, new_start

    i1 = this_meshblock % ptr % i1
    i2 = this_meshblock % ptr % i2
    j1 = this_meshblock % ptr % j1
    j2 = this_meshblock % ptr % j2
    k1 = this_meshblock % ptr % k1
    k2 = this_meshblock % ptr % k2

    if (shift >= this_meshblock % ptr % sx) then
      call throwError('Requested moving-window shift is larger than the local domain in x.')
    end if

    ex(i1:i2-shift, j1:j2, k1:k2) = ex(i1+shift:i2, j1:j2, k1:k2)
    ey(i1:i2-shift, j1:j2, k1:k2) = ey(i1+shift:i2, j1:j2, k1:k2)
    ez(i1:i2-shift, j1:j2, k1:k2) = ez(i1+shift:i2, j1:j2, k1:k2)
    bx(i1:i2-shift, j1:j2, k1:k2) = bx(i1+shift:i2, j1:j2, k1:k2)
    by(i1:i2-shift, j1:j2, k1:k2) = by(i1+shift:i2, j1:j2, k1:k2)
    bz(i1:i2-shift, j1:j2, k1:k2) = bz(i1+shift:i2, j1:j2, k1:k2)

    new_start = i2 - shift + 1
    ex(new_start:i2, j1:j2, k1:k2) = 0.0
    ey(new_start:i2, j1:j2, k1:k2) = 0.0
    ez(new_start:i2, j1:j2, k1:k2) = 0.0
    bx(new_start:i2, j1:j2, k1:k2) = 0.0
    by(new_start:i2, j1:j2, k1:k2) = 0.0
    bz(new_start:i2, j1:j2, k1:k2) = 0.0

    jx(:, :, :) = 0.0
    jy(:, :, :) = 0.0
    jz(:, :, :) = 0.0
  end subroutine shift_fields

  subroutine shift_particles(shift)
    implicit none
    integer, intent(in) :: shift
    integer :: s, ti, tj, tk, p
    integer :: ti_new, tj_new, tk_new
    integer(kind=2) :: xi_new, yi_new, zi_new
    type(particle_tile), pointer :: old_tiles(:, :, :)
    type(particle_tile), allocatable :: new_tiles(:, :, :)

    if (shift <= 0) return

    do s = 1, nspec
      old_tiles => species(s) % prtl_tile

      allocate (new_tiles(species(s) % tile_nx, species(s) % tile_ny, species(s) % tile_nz))
      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            call initialize_tile(new_tiles(ti, tj, tk), s, ti, tj, tk)
          end do
        end do
      end do

      species(s) % cntr_sp = 0

      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            do p = 1, old_tiles(ti, tj, tk) % npart_sp
              xi_new = old_tiles(ti, tj, tk) % xi(p) - shift
              if ((xi_new < 0) .or. (xi_new >= this_meshblock % ptr % sx)) cycle
              yi_new = old_tiles(ti, tj, tk) % yi(p)
              zi_new = old_tiles(ti, tj, tk) % zi(p)

              ti_new = FLOOR(REAL(xi_new) / REAL(species(s) % tile_sx)) + 1
              tj_new = FLOOR(REAL(yi_new) / REAL(species(s) % tile_sy)) + 1
              tk_new = FLOOR(REAL(zi_new) / REAL(species(s) % tile_sz)) + 1

              call append_to_tile(new_tiles(ti_new, tj_new, tk_new), xi_new, yi_new, zi_new, &
                                   old_tiles(ti, tj, tk) % dx(p), old_tiles(ti, tj, tk) % dy(p), &
                                   old_tiles(ti, tj, tk) % dz(p), &
                                   old_tiles(ti, tj, tk) % u(p), old_tiles(ti, tj, tk) % v(p), &
                                   old_tiles(ti, tj, tk) % w(p), &
                                   old_tiles(ti, tj, tk) % ind(p), 0, old_tiles(ti, tj, tk) % weight(p))

              species(s) % cntr_sp = species(s) % cntr_sp + 1
            end do

            call clear_tile(old_tiles(ti, tj, tk))
          end do
        end do
      end do

      if (allocated(species(s) % prtl_tile)) deallocate (species(s) % prtl_tile)
      call move_alloc(new_tiles, species(s) % prtl_tile)
    end do
  end subroutine shift_particles

  subroutine refill_new_region(shift)
    implicit none
    integer, intent(in) :: shift
    real :: xmin, xmax

    xmin = REAL(this_meshblock % ptr % x0 + this_meshblock % ptr % sx - shift, kind=dprec)
    xmax = REAL(this_meshblock % ptr % x0 + this_meshblock % ptr % sx - 1, kind=dprec)

    call userFillNewRegion(xmin, xmax)
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
