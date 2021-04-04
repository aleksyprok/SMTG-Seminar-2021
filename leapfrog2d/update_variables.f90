MODULE update_variables

  USE shared_data
  USE control

  IMPLICIT NONE

CONTAINS

  SUBROUTINE predictor_corrector

    REAL(num) :: w1, w2, w3, w4, wa

    ! Predictor step
    DO iz = 0, nz
      izp = iz + 1
      izm = iz - 1
      DO iy = 0, ny
        iyp = iy + 1
        iym = iy - 1

        IF (iy .NE. 0 .AND. iz .NE. 0) THEN
          w1 = by0 * (bx2(iy,iz) - bx2(iym,iz )) / dy
          w2 = bz0 * (bx1(iy,iz) - bx1(iy ,izm)) / dz
          w3 = kx * b_par1(iy,iz)
          wa = w1 + w2 + w3
          ux1_pred(iy,iz) = ux1m(iy,iz) + 2.0_num * dt * wa
        END IF

        IF (iz .NE. 0 .AND. iz .NE. nz) THEN
          w1 = by0 * (bx1(iyp,iz ) - bx1(iy,iz)) / dy
          w2 = bz0 * (bx2(iy ,izp) - bx2(iy,iz)) / dz
          w3 = kx * b_par2(iy,iz)
          wa = w1 + w2 + w3
          ux2_pred(iy,iz) = ux2m(iy,iz) + 2.0_num * dt * wa
        END IF

        IF (iy .NE. 0) THEN
          w1 = by0 * (ux2(iy,iz ) - ux2(iym,iz)) / dy
          w2 = bz0 * (ux1(iy,izp) - ux1(iy ,iz)) / dz
          bx1_pred(iy,iz) = bx1m(iy,iz) + 2.0_num * dt * (w1 + w2)
        END IF

        IF (iz .NE. 0) THEN
          w1 = by0 * (ux1(iyp,iz) - ux1(iy,iz )) / dy
          w2 = bz0 * (ux2(iy ,iz) - ux2(iy,izm)) / dz
          bx2_pred(iy,iz) = bx2m(iy,iz) + 2.0_num * dt * (w1 + w2)
        END IF

        IF (iy .NE. 0 .AND. iz .NE. 0 .AND. iz .NE. nz) THEN
          w1 = by0 * (b_perp2(iy,iz ) - b_perp2(iym,iz)) / dy
          w2 = bz0 * (b_perp1(iy,izp) - b_perp1(iy ,iz)) / dz
          w3 = bz0 * (b_par2(iy,iz ) - b_par2(iym,iz)) / dy
          w4 = by0 * (b_par1(iy,izp) - b_par1(iy ,iz)) / dz
          wa = w1 + w2 - w3 + w4
          u_perp1_pred(iy,iz) = u_perp1m(iy,iz) + 2.0_num * dt * wa
        END IF

        IF (iz .NE. 0) THEN
          w1 = by0 * (b_perp1(iyp,iz) - b_perp1(iy,iz )) / dy
          w2 = bz0 * (b_perp2(iy ,iz) - b_perp2(iy,izm)) / dz
          w3 = bz0 * (b_par1(iyp,iz) - b_par1(iy,iz )) / dy
          w4 = by0 * (b_par2(iy ,iz) - b_par2(iy,izm)) / dz
          wa = w1 + w2 - w3 + w4
          u_perp2_pred(iy,iz) = u_perp2m(iy,iz) + 2.0_num * dt * wa
        END IF

        IF (iy .NE. 0 .AND. iz .NE. 0) THEN
          w1 = by0 * (u_perp2(iy,iz) - u_perp2(iym,iz )) / dy
          w2 = bz0 * (u_perp1(iy,iz) - u_perp1(iy ,izm)) / dz
          b_perp1_pred(iy,iz) = b_perp1m(iy,iz) + 2.0_num * dt * (w1 + w2)
        END IF

        w1 = by0 * (u_perp1(iyp,iz ) - u_perp1(iy,iz)) / dy
        w2 = bz0 * (u_perp2(iy ,izp) - u_perp2(iy,iz)) / dz
        b_perp2_pred(iy,iz) = b_perp2m(iy,iz) + 2.0_num * dt * (w1 + w2)

        IF (iy .NE. 0 .AND. iz .NE. 0) THEN
          w1 = kx * ux1(iy,iz)
          w2 = bz0 * (u_perp2(iy,iz) - u_perp2(iym,iz)) / dy
          w3 = by0 * (u_perp1(iy,iz) - u_perp1(iy,izm)) / dz
          b_par1_pred(iy,iz) = b_par1m(iy,iz) - 2.0_num * dt * (w1 + w2 - w3)
        END IF

        w1 = kx * ux2(iy,iz)
        w2 = bz0 * (u_perp1(iyp,iz ) - u_perp1(iy,iz)) / dy
        w3 = by0 * (u_perp2(iy ,izp) - u_perp2(iy,iz)) / dz
        b_par2_pred(iy,iz) = b_par2m(iy,iz) - 2.0_num * dt * (w1 + w2 - w3)

      END DO
    END DO

    CALL set_boundary_conditions(PRED_VARIABLES, time + dt)

    ! Overwrite m variables
    ux1m = ux1
    ux2m = ux2
    u_perp1m = u_perp1
    u_perp2m = u_perp2
    bx1m = bx1
    bx2m = bx2
    b_perp1m = b_perp1
    b_perp2m = b_perp2
    b_par1m = b_par1
    b_par2m = b_par2

    ! Corrector step
    DO iz = 0, nz
      izp = iz + 1
      izm = iz - 1
      DO iy = 0, ny
        iyp = iy + 1
        iym = iy - 1

        IF (iy .NE. 0 .AND. iz .NE. 0) THEN
          w1 = by0 * (bx2m(iy,iz) - bx2m(iym,iz )) / dy
          w2 = bz0 * (bx1m(iy,iz) - bx1m(iy ,izm)) / dz
          w3 = kx * b_par1m(iy,iz)
          wa = w1 + w2 + w3
          w1 = by0 * (bx2_pred(iy,iz) - bx2_pred(iym,iz )) / dy
          w2 = bz0 * (bx1_pred(iy,iz) - bx1_pred(iy ,izm)) / dz
          w3 = kx * b_par1_pred(iy,iz)
          wa = 0.5_num * (wa + w1 + w2 + w3)
          ux1(iy,iz) = ux1(iy,iz) + dt * wa
        END IF

        IF (iz .NE. 0 .AND. iz .NE. nz) THEN
          w1 = by0 * (bx1m(iyp,iz ) - bx1m(iy,iz)) / dy
          w2 = bz0 * (bx2m(iy ,izp) - bx2m(iy,iz)) / dz
          w3 = kx * b_par2m(iy,iz)
          wa = w1 + w2 + w3
          w1 = by0 * (bx1_pred(iyp,iz ) - bx1_pred(iy,iz)) / dy
          w2 = bz0 * (bx2_pred(iy ,izp) - bx2_pred(iy,iz)) / dz
          w3 = kx * b_par2_pred(iy,iz)
          wa = 0.5_num * (wa + w1 + w2 + w3)
          ux2(iy,iz) = ux2(iy,iz) + dt * wa
        END IF

        IF (iy .NE. 0) THEN
          w1 = by0 * (ux2m(iy,iz ) - ux2m(iym,iz)) / dy
          w2 = bz0 * (ux1m(iy,izp) - ux1m(iy ,iz)) / dz
          wa = w1 + w2
          w1 = by0 * (ux2_pred(iy,iz ) - ux2_pred(iym,iz)) / dy
          w2 = bz0 * (ux1_pred(iy,izp) - ux1_pred(iy ,iz)) / dz
          wa = 0.5_num * (wa + w1 + w2)
          bx1(iy,iz) = bx1(iy,iz) + dt * wa
        END IF

        IF (iz .NE. 0) THEN
          w1 = by0 * (ux1m(iyp,iz) - ux1m(iy,iz )) / dy
          w2 = bz0 * (ux2m(iy ,iz) - ux2m(iy,izm)) / dz
          wa = w1 + w2
          w1 = by0 * (ux1_pred(iyp,iz) - ux1_pred(iy,iz )) / dy
          w2 = bz0 * (ux2_pred(iy ,iz) - ux2_pred(iy,izm)) / dz
          wa = 0.5_num * (wa + w1 + w2)
          bx2(iy,iz) = bx2(iy,iz) + dt * wa
        END IF

        IF (iy .NE. 0 .AND. iz .NE. 0 .AND. iz .NE. nz) THEN
          w1 = by0 * (b_perp2m(iy,iz ) - b_perp2m(iym,iz)) / dy
          w2 = bz0 * (b_perp1m(iy,izp) - b_perp1m(iy,iz )) / dz
          w3 = bz0 * (b_par2m(iy,iz ) - b_par2m(iym,iz)) / dy
          w4 = by0 * (b_par1m(iy,izp) - b_par1m(iy ,iz)) / dz
          wa = w1 + w2 - w3 + w4
          w1 = by0 * (b_perp2_pred(iy,iz ) - b_perp2_pred(iym,iz)) / dy
          w2 = bz0 * (b_perp1_pred(iy,izp) - b_perp1_pred(iy ,iz)) / dz
          w3 = bz0 * (b_par2_pred(iy,iz ) - b_par2_pred(iym,iz)) / dy
          w4 = by0 * (b_par1_pred(iy,izp) - b_par1_pred(iy ,iz)) / dz
          wa = 0.5_num * (wa + w1 + w2 - w3 + w4)
          u_perp1(iy,iz) = u_perp1(iy,iz) + dt * wa
        END IF

        IF (iz .NE. 0) THEN
          w1 = by0 * (b_perp1m(iyp,iz) - b_perp1m(iy,iz )) / dy
          w2 = bz0 * (b_perp2m(iy ,iz) - b_perp2m(iy,izm)) / dz
          w3 = bz0 * (b_par1m(iyp,iz) - b_par1m(iy,iz )) / dy
          w4 = by0 * (b_par2m(iy ,iz) - b_par2m(iy,izm)) / dz
          wa = w1 + w2 - w3 + w4
          w1 = by0 * (b_perp1_pred(iyp,iz) - b_perp1_pred(iy,iz )) / dy
          w2 = bz0 * (b_perp2_pred(iy ,iz) - b_perp2_pred(iy,izm)) / dz
          w3 = bz0 * (b_par1_pred(iyp,iz) - b_par1_pred(iy,iz )) / dy
          w4 = by0 * (b_par2_pred(iy ,iz) - b_par2_pred(iy,izm)) / dz
          wa = 0.5_num * (wa + w1 + w2 - w3 + w4)
          u_perp2(iy,iz) = u_perp2(iy,iz) + dt * wa
        END IF

        IF (iy .NE. 0 .AND. iz .NE. 0) THEN
          w1 = by0 * (u_perp2m(iy,iz) - u_perp2m(iym,iz )) / dy
          w2 = bz0 * (u_perp1m(iy,iz) - u_perp1m(iy ,izm)) / dz
          wa = w1 + w2
          w1 = by0 * (u_perp2_pred(iy,iz) - u_perp2_pred(iym,iz )) / dy
          w2 = bz0 * (u_perp1_pred(iy,iz) - u_perp1_pred(iy ,izm)) / dz
          wa = 0.5_num * (wa + w1 + w2)
          b_perp1(iy,iz) = b_perp1(iy,iz) + dt * wa
        END IF

        w1 = by0 * (u_perp1m(iyp,iz ) - u_perp1m(iy,iz)) / dy
        w2 = bz0 * (u_perp2m(iy ,izp) - u_perp2m(iy,iz)) / dz
        wa = w1 + w2
        w1 = by0 * (u_perp1_pred(iyp,iz ) - u_perp1_pred(iy,iz)) / dy
        w2 = bz0 * (u_perp2_pred(iy ,izp) - u_perp2_pred(iy,iz)) / dz
        wa = 0.5_num * (wa + w1 + w2)
        b_perp2(iy,iz) = b_perp2(iy,iz) + dt * wa

        IF (iy .NE. 0 .AND. iz .NE. 0) THEN
          w1 = kx * ux1m(iy,iz)
          w2 = bz0 * (u_perp2m(iy,iz) - u_perp2m(iym,iz )) / dy
          w3 = by0 * (u_perp1m(iy,iz) - u_perp1m(iy ,izm)) / dz
          wa = -(w1 + w2 - w3)
          w1 = kx * ux1_pred(iy,iz)
          w2 = bz0 * (u_perp2_pred(iy,iz) - u_perp2_pred(iym,iz )) / dy
          w3 = by0 * (u_perp1_pred(iy,iz) - u_perp1_pred(iy ,izm)) / dz
          wa = 0.5_num * (wa - (w1 + w2 - w3))
          b_par1(iy,iz) = b_par1(iy,iz) + dt * wa
        END IF

        w1 = kx * ux2m(iy,iz)
        w2 = bz0 * (u_perp1m(iyp,iz ) - u_perp1m(iy,iz)) / dy
        w3 = by0 * (u_perp2m(iy ,izp) - u_perp2m(iy,iz)) / dz
        wa = -(w1 + w2 - w3)
        w1 = kx * ux2_pred(iy,iz)
        w2 = bz0 * (u_perp1_pred(iyp,iz ) - u_perp1_pred(iy,iz)) / dy
        w3 = by0 * (u_perp2_pred(iy ,izp) - u_perp2_pred(iy,iz)) / dz
        wa = 0.5_num * (wa - (w1 + w2 - w3))
        b_par2(iy,iz) = b_par2(iy,iz) + dt * wa

      END DO
    END DO

    CALL set_boundary_conditions(CURR_VARIABLES, time + dt)

  END SUBROUTINE predictor_corrector

END MODULE update_variables
