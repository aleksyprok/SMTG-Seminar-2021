MODULE constants

  IMPLICIT NONE

  INTEGER, PARAMETER :: num = KIND(0.0D0)
  REAL(num), PARAMETER :: pi = 4.0D0 * DATAN(1.0D0)

END MODULE constants

MODULE shared_data

  USE constants

  IMPLICIT NONE

  INTEGER :: ny, nz, step = 0
  INTEGER :: iy, iz, iyp, izp, iym, izm
  INTEGER :: n_iters, file_number = 0
  INTEGER :: tstart, tstop, count_rate
  INTEGER, PARAMETER :: dump_unit = 11
  INTEGER, PARAMETER :: CURR_VARIABLES = 1, PRED_VARIABLES = 2, M_VARIABLES = 3
  REAL(num), DIMENSION(:), ALLOCATABLE :: yb, zb, yc, zc
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: ux1, ux1m, ux1_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: ux2, ux2m, ux2_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: u_perp1, u_perp1m, u_perp1_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: u_perp2, u_perp2m, u_perp2_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: bx1, bx1m, bx1_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: bx2, bx2m, bx2_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: b_perp1, b_perp1m, b_perp1_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: b_perp2, b_perp2m, b_perp2_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: b_par1, b_par1m, b_par1_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: b_par2, b_par2m, b_par2_pred
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: rho_bcc, rho_bbb, rho_ccb, rho_cbc
  REAL(num) :: y_min, y_max, z_min, z_max, t_max
  REAL(num) :: dy, dz, dt, dt_multiplier
  REAL(num) :: alpha, kx, k_par, by0, bz0
  REAL(num) :: time = 0.0_num, dt_snapshots, time_dump
  LOGICAL :: output_ux1, output_ux2, output_u_perp1, output_u_perp2
  LOGICAL :: output_bx1, output_bx2, output_b_perp1, output_b_perp2
  LOGICAL :: output_b_par1, output_b_par2

END MODULE shared_data
