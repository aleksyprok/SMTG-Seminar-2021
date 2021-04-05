MODULE control

  USE shared_data
  USE constants

  IMPLICIT NONE

CONTAINS

  SUBROUTINE control_variables

    alpha = DATAN(0.25_num)
    k_par = 0.5_num * pi
    kx = 0.0_num * k_par
    ny = 512
    nz = 256
    y_min = -4.0_num
    y_max =  8.0_num
    z_min =  0.0_num
    z_max =  6.0_num
    t_max =  6.0_num
    n_iters = -1

    output_ux1     = .TRUE.
    output_ux2     = .TRUE.
    output_u_perp1 = .TRUE.
    output_u_perp2 = .TRUE.
    output_bx1     = .TRUE.
    output_bx2     = .TRUE.
    output_b_perp1 = .TRUE.
    output_b_perp2 = .TRUE.
    output_b_par1  = .TRUE.
    output_b_par2  = .TRUE.

    dt_snapshots = 0.1_num

    dt_multiplier = 0.4_num

  END SUBROUTINE control_variables

  FUNCTION half_cylinder(z)

    REAL(num), INTENT(IN) :: z
    REAL(num) :: half_cylinder

    IF (DABS(z) .LE. 1) THEN
      half_cylinder = DCOS(k_par * z) ** 2.0_num
    ELSE
      half_cylinder = 0.0_num
    END IF

  END FUNCTION half_cylinder

  FUNCTION dome(y, z)

    REAL(num), INTENT(IN) :: y, z
    REAL(num) :: dome
    REAL(num) :: r


    r = DSQRT(y * y + z * z)
    IF (DABS(r) .LE. 1) THEN
      dome = DCOS(k_par * r) ** 2.0_num
    ELSE
      dome = 0.0_num
    END IF

  END FUNCTION dome

  FUNCTION u_ana(y, z, t)

    REAL(num), INTENT(IN) :: z, y, t
    REAL(num) :: u_ana

    u_ana = half_cylinder(by0 * y + bz0 * (z - 4.0_num) + t)

  END FUNCTION

  FUNCTION b_ana(y, z, t)

    REAL(num), INTENT(IN) :: z, y, t
    REAL(num) :: b_ana

    b_ana = half_cylinder(by0 * y + bz0 * (z - 4.0_num) + t)

  END FUNCTION

  SUBROUTINE set_initial_conditions

    REAL(num) :: Lz, kz

    ux1  = 0.0_num
    ux1m = 0.0_num
    ux2  = 0.0_num
    ux2m = 0.0_num
    u_perp1  = 0.0_num
    u_perp1m = 0.0_num
    u_perp2  = 0.0_num
    u_perp2m = 0.0_num
    bx1  = 0.0_num
    bx1m = 0.0_num
    bx2  = 0.0_num
    bx2m = 0.0_num
    b_perp1  = 0.0_num
    b_perp1m = 0.0_num
    b_perp2  = 0.0_num
    b_perp2m = 0.0_num
    b_par1  = 0.0_num
    b_par1m = 0.0_num
    b_par2  = 0.0_num
    b_par2m = 0.0_num

    Lz = z_max - z_min
    kz = 2.0_num * pi / Lz

    ! DO iz = 0, nz
    !   DO iy = 0, ny + 1
    !     u_perp1( iy,iz) = u_ana(yc(iy), zb(iz), 0.0_num)
    !     u_perp1m(iy,iz) = u_ana(yc(iy), zb(iz), -dt)
    !   END DO
    ! END DO
    !
    ! DO iz = 0, nz + 1
    !   DO iy = 0, ny
    !     u_perp2( iy,iz) = u_ana(yb(iy), zc(iz), 0.0_num)
    !     u_perp2m(iy,iz) = u_ana(yb(iy), zc(iz), -dt)
    !   END DO
    ! END DO
    !
    ! DO iz = 1, nz
    !   DO iy = 0, ny + 1
    !     b_perp1( iy,iz) = b_ana(yc(iy), zc(iz), 0.0_num)
    !     b_perp1m(iy,iz) = b_ana(yc(iy), zc(iz), -dt)
    !   END DO
    ! END DO
    !
    ! DO iz = 0, nz
    !   DO iy = 0, ny
    !     b_perp2( iy,iz) = b_ana(yb(iy), zb(iz), 0.0_num)
    !     b_perp2m(iy,iz) = b_ana(yb(iy), zb(iz), -dt)
    !   END DO
    ! END DO

    DO iz = 0, nz + 1
      DO iy = 0, ny + 1
        ux1( iy,iz) = dome(yc(iy) - 1.0_num, zc(iz) - 4.0_num)
        ux1m(iy,iz) = dome(yc(iy) - 1.0_num, zc(iz) - 4.0_num)
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 0, ny
        ux2( iy,iz) = dome(yb(iy) - 1.0_num, zb(iz) - 4.0_num)
        ux2m(iy,iz) = dome(yb(iy) - 1.0_num, zb(iz) - 4.0_num)
      END DO
    END DO

    DO iz = 0, nz
      DO iy = 0, ny + 1
        bx1( iy,iz) = dome(yc(iy) - 1.0_num, zb(iz) - 4.0_num)
        bx1m(iy,iz) = dome(yc(iy) - 1.0_num, zb(iz) - 4.0_num)
      END DO
    END DO

    DO iz = 1, nz
      DO iy = 0, ny
        bx2( iy,iz) = dome(yb(iy) - 1.0_num, zc(iz) - 4.0_num)
        bx2m(iy,iz) = dome(yb(iy) - 1.0_num, zc(iz) - 4.0_num)
      END DO
    END DO

  END SUBROUTINE set_initial_conditions

  SUBROUTINE set_boundary_conditions(variables, time)

    INTEGER, INTENT(IN) :: variables
    REAL(num), INTENT(IN) :: time

    IF (variables == CURR_VARIABLES) THEN
      CALL boundary_conditions(ux1, ux2, u_perp1, u_perp2, bx1, bx2, &
                               b_perp1, b_perp2, b_par1, b_par2, time)
    ELSE IF (variables == PRED_VARIABLES) THEN
      CALL boundary_conditions(ux1_pred, ux2_pred, u_perp1_pred, u_perp2_pred, &
                               bx1_pred, bx2_pred, b_perp1_pred, b_perp2_pred, &
                               b_par1_pred, b_par2_pred, time)
    ELSE IF (variables == M_VARIABLES) THEN
      CALL boundary_conditions(ux1m, ux2m, u_perp1m, u_perp2m, bx1m, bx2m, &
                               b_perp1m, b_perp2m, b_par1m, b_par2m, time)
    END IF

  END SUBROUTINE set_boundary_conditions

  SUBROUTINE boundary_conditions(ux1, ux2, u_perp1, u_perp2, bx1, bx2, &
                                 b_perp1, b_perp2, b_par1, b_par2, time)

    REAL(num), DIMENSION(0:,0:), INTENT(INOUT) :: ux1, ux2, bx1
    REAL(num), DIMENSION(0:,0:), INTENT(INOUT) :: u_perp1, u_perp2, b_perp2, b_par2
    REAL(num), DIMENSION(0:,1:), INTENT(INOUT) :: bx2, b_perp1, b_par1
    REAL(num), INTENT(IN) :: time

    ! y = y_min
    ux1(0,1:nz) = ux1(1,1:nz) - DTAN(alpha) * dy &
                * (ux2(0,1:nz) - ux2(0,0:nz-1)) / dz
    ux1(0,0) = ux1(1,0)
    ux1(0,nz+1) = ux1(1,nz+1)
    u_perp1(0,:) = u_perp1(1,:) - DTAN(alpha) * dy &
                 * (u_perp2(0,1:nz+1) - u_perp2(0,0:nz)) / dz
    bx1(0,1:nz-1) = bx1(1,1:nz-1) - DTAN(alpha) * dy &
                  * (bx2(0,2:nz) - bx2(0,1:nz-1)) / dz
    bx1(0,0) = bx1(1,0)
    bx1(0,nz) = bx1(1,nz)
    b_perp1(0,:) = b_perp1(1,:) - DTAN(alpha) * dy &
                 * (b_perp2(0,1:nz) - b_perp2(0,0:nz-1)) / dz
    b_par1(0,:) = b_par1(1,:) - DTAN(alpha) * dy &
                * (b_par2(0,1:nz) - b_par2(0,0:nz-1)) / dz

    ! y = y_max
    ux1(ny+1,1:nz) = ux1(ny,1:nz) + DTAN(alpha) * dy &
                   * (ux2(ny,1:nz) - ux2(ny,0:nz-1)) / dz
    ux1(ny+1,0) = ux1(ny,0)
    ux1(ny+1,nz+1) = ux1(ny,nz+1)
    u_perp1(ny+1,:) = u_perp1(ny,:) + DTAN(alpha) * dy &
                    * (u_perp2(ny,1:nz+1) - u_perp2(ny,0:nz)) / dz
    bx1(ny+1,1:nz-1) = bx1(ny,1:nz-1) + DTAN(alpha) * dy &
                     * (bx2(ny,2:nz) - bx2(ny,1:nz-1)) / dz
    bx1(ny+1,0) = bx1(ny,0)
    bx1(ny+1,nz) = bx1(ny,nz)
    b_perp1(ny+1,:) = b_perp1(ny,:) + DTAN(alpha) * dy &
                    * (b_perp2(ny,1:nz) - b_perp2(ny,0:nz-1)) / dz
    b_par1(ny+1,:) = b_par1(ny,:) + DTAN(alpha) * dy &
                   * (b_par2(ny,1:nz) - b_par2(ny,0:nz-1)) / dz

    ! ! y = y_min
    ! ux1(0,:) = ux1(1,:)
    ! u_perp1(0,:) = u_perp1(1,:)
    ! bx1(0,:) = bx1(1,:)
    ! b_perp1(0,:) = b_perp1(1,:)
    ! b_par1(0,:) = b_par1(1,:)
    !
    ! ! y = y_max
    ! ux1(ny+1,:) = ux1(ny,:)
    ! u_perp1(ny+1,:) = u_perp1(ny,:)
    ! bx1(ny+1,:) = bx1(ny,:)
    ! b_perp1(ny+1,:) = b_perp1(ny,:)
    ! b_par1(ny+1,:) = b_par1(ny,:)

    ! z = z_min
    ux1(:,0) = -ux1(:,1)
    ux2(:,0) = 0.0_num
    u_perp1(:,0) = 0.0_num
    u_perp2(:,0) = -u_perp2(:,1)

    ! z = z_max
    ux1(:,nz+1) = -ux1(:,nz)
    ux2(:,nz  ) = 0.0_num
    u_perp1(:,nz  ) = 0.0_num
    u_perp2(:,nz+1) = -u_perp2(:,nz)

  END SUBROUTINE boundary_conditions

END MODULE control
