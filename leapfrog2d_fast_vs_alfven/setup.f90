MODULE setup

  USE shared_data
  USE control

  IMPLICIT NONE

CONTAINS

  SUBROUTINE after_control

    REAL(num) :: length

    ALLOCATE(yb(0:ny))
    ALLOCATE(zb(0:nz))
    ALLOCATE(yc(0:ny+1))
    ALLOCATE(zc(0:nz+1))

    ALLOCATE(     ux1(0:ny+1,0:nz+1))
    ALLOCATE(    ux1m(0:ny+1,0:nz+1))
    ALLOCATE(ux1_pred(0:ny+1,0:nz+1))

    ALLOCATE(     ux2(0:ny,0:nz))
    ALLOCATE(    ux2m(0:ny,0:nz))
    ALLOCATE(ux2_pred(0:ny,0:nz))

    ALLOCATE(     u_perp1(0:ny+1,0:nz))
    ALLOCATE(    u_perp1m(0:ny+1,0:nz))
    ALLOCATE(u_perp1_pred(0:ny+1,0:nz))

    ALLOCATE(     u_perp2(0:ny,0:nz+1))
    ALLOCATE(    u_perp2m(0:ny,0:nz+1))
    ALLOCATE(u_perp2_pred(0:ny,0:nz+1))

    ALLOCATE(     bx1(0:ny+1,0:nz))
    ALLOCATE(    bx1m(0:ny+1,0:nz))
    ALLOCATE(bx1_pred(0:ny+1,0:nz))

    ALLOCATE(     bx2(0:ny,1:nz))
    ALLOCATE(    bx2m(0:ny,1:nz))
    ALLOCATE(bx2_pred(0:ny,1:nz))

    ALLOCATE(     b_perp1(0:ny+1,1:nz))
    ALLOCATE(    b_perp1m(0:ny+1,1:nz))
    ALLOCATE(b_perp1_pred(0:ny+1,1:nz))

    ALLOCATE(     b_perp2(0:ny,0:nz))
    ALLOCATE(    b_perp2m(0:ny,0:nz))
    ALLOCATE(b_perp2_pred(0:ny,0:nz))

    ALLOCATE(     b_par1(0:ny+1,1:nz))
    ALLOCATE(    b_par1m(0:ny+1,1:nz))
    ALLOCATE(b_par1_pred(0:ny+1,1:nz))

    ALLOCATE(     b_par2(0:ny,0:nz))
    ALLOCATE(    b_par2m(0:ny,0:nz))
    ALLOCATE(b_par2_pred(0:ny,0:nz))

    ! Create grid
    dy = (y_max - y_min) / REAL(ny, num)
    dz = (z_max - z_min) / REAL(nz, num)
    DO iy = 0, ny
      yb(iy) = y_min + iy * dy
    END DO
    DO iz = 0, nz
      zb(iz) = z_min + iz * dz
    END DO
    DO iy = 0, ny + 1
      yc(iy) = y_min + (iy - 0.5_num) * dy
    END DO
    DO iz = 0, nz + 1
      zc(iz) = z_min + (iz - 0.5_num) * dz
    END DO

    by0 = SIN(alpha)
    bz0 = COS(alpha)

    ! CFL condition
    length = MIN(dy, dz)
    dt = dt_multiplier * length
    IF (n_iters < 0) n_iters = INT(t_max / (REAL(dt, num))) + 2

  END SUBROUTINE after_control

  SUBROUTINE finalize

    DEALLOCATE(yb)
    DEALLOCATE(zb)
    DEALLOCATE(yc)
    DEALLOCATE(zc)

    DEALLOCATE(     ux1)
    DEALLOCATE(    ux1m)
    DEALLOCATE(ux1_pred)

    DEALLOCATE(     ux2)
    DEALLOCATE(    ux2m)
    DEALLOCATE(ux2_pred)

    DEALLOCATE(     u_perp1)
    DEALLOCATE(    u_perp1m)
    DEALLOCATE(u_perp1_pred)

    DEALLOCATE(     u_perp2)
    DEALLOCATE(    u_perp2m)
    DEALLOCATE(u_perp2_pred)

    DEALLOCATE(     bx1)
    DEALLOCATE(    bx1m)
    DEALLOCATE(bx1_pred)

    DEALLOCATE(     bx2)
    DEALLOCATE(    bx2m)
    DEALLOCATE(bx2_pred)

    DEALLOCATE(     b_perp1)
    DEALLOCATE(    b_perp1m)
    DEALLOCATE(b_perp1_pred)

    DEALLOCATE(     b_perp2)
    DEALLOCATE(    b_perp2m)
    DEALLOCATE(b_perp2_pred)

    DEALLOCATE(     b_par1)
    DEALLOCATE(    b_par1m)
    DEALLOCATE(b_par1_pred)

    DEALLOCATE(     b_par2)
    DEALLOCATE(    b_par2m)
    DEALLOCATE(b_par2_pred)

  END SUBROUTINE finalize

END MODULE setup
