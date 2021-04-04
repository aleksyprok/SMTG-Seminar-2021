MODULE diagnostics

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE output_grid

    OPEN(UNIT = dump_unit, FILE = "Data/grid.dat", FORM = "UNFORMATTED",&
         STATUS = "REPLACE", ACTION = "WRITE")

    WRITE(UNIT = dump_unit) t_max
    WRITE(UNIT = dump_unit) dt
    WRITE(UNIT = dump_unit) ny
    WRITE(UNIT = dump_unit) nz
    WRITE(UNIT = dump_unit) yb
    WRITE(UNIT = dump_unit) zb
    WRITE(UNIT = dump_unit) yc
    WRITE(UNIT = dump_unit) zc
    WRITE(UNIT = dump_unit) alpha
    WRITE(UNIT = dump_unit) k_par
    WRITE(UNIT = dump_unit) kx

    CLOSE(UNIT = dump_unit)

  END SUBROUTINE output_grid

  SUBROUTINE output_dump

    CHARACTER (LEN = 10) :: ux1_file
    CHARACTER (LEN = 10) :: ux2_file
    CHARACTER (LEN = 14) :: u_perp1_file
    CHARACTER (LEN = 14) :: u_perp2_file
    CHARACTER (LEN = 10) :: bx1_file
    CHARACTER (LEN = 10) :: bx2_file
    CHARACTER (LEN = 14) :: b_perp1_file
    CHARACTER (LEN = 14) :: b_perp2_file
    CHARACTER (LEN = 13) :: b_par1_file
    CHARACTER (LEN = 13) :: b_par2_file
    CHARACTER (LEN = 11) :: step_file

    IF (output_ux1) THEN
      WRITE(ux1_file , "(A3, I3.3, A4)") "ux1" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // ux1_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(ux1, 1)
      WRITE(UNIT = dump_unit) SIZE(ux1, 2)
      WRITE(UNIT = dump_unit) ux1
      CLOSE(UNIT = dump_unit)
    END IF

    IF (output_ux2) THEN
      WRITE(ux2_file , "(A3, I3.3, A4)") "ux2" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // ux2_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(ux2, 1)
      WRITE(UNIT = dump_unit) SIZE(ux2, 2)
      WRITE(UNIT = dump_unit) ux2
      CLOSE(UNIT = dump_unit)
    END IF

    IF (output_u_perp1) THEN
      WRITE(u_perp1_file , "(A7, I3.3, A4)") "u_perp1" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // u_perp1_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(u_perp1, 1)
      WRITE(UNIT = dump_unit) SIZE(u_perp1, 2)
      WRITE(UNIT = dump_unit) u_perp1
      CLOSE(UNIT = dump_unit)
    END IF

    IF (output_u_perp2) THEN
      WRITE(u_perp2_file , "(A7, I3.3, A4)") "u_perp2" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // u_perp2_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(u_perp2, 1)
      WRITE(UNIT = dump_unit) SIZE(u_perp2, 2)
      WRITE(UNIT = dump_unit) u_perp2
      CLOSE(UNIT = dump_unit)
    END IF

    IF (output_bx1) THEN
      WRITE(bx1_file , "(A3, I3.3, A4)") "bx1" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // bx1_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(bx1, 1)
      WRITE(UNIT = dump_unit) SIZE(bx1, 2)
      WRITE(UNIT = dump_unit) bx1
      CLOSE(UNIT = dump_unit)
    END IF

    IF (output_bx2) THEN
      WRITE(bx2_file , "(A3, I3.3, A4)") "bx2" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // bx2_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(bx2, 1)
      WRITE(UNIT = dump_unit) SIZE(bx2, 2)
      WRITE(UNIT = dump_unit) bx2
      CLOSE(UNIT = dump_unit)
    END IF

    IF (output_b_perp1) THEN
      WRITE(b_perp1_file , "(A7, I3.3, A4)") "b_perp1" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // b_perp1_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(b_perp1, 1)
      WRITE(UNIT = dump_unit) SIZE(b_perp1, 2)
      WRITE(UNIT = dump_unit) b_perp1
      CLOSE(UNIT = dump_unit)
    END IF

    IF (output_b_perp2) THEN
      WRITE(b_perp2_file , "(A7, I3.3, A4)") "b_perp2" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // b_perp2_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(b_perp2, 1)
      WRITE(UNIT = dump_unit) SIZE(b_perp2, 2)
      WRITE(UNIT = dump_unit) b_perp2
      CLOSE(UNIT = dump_unit)
    END IF

    IF (output_b_par1) THEN
      WRITE(b_par1_file , "(A6, I3.3, A4)") "b_par1" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // b_par1_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(b_par1, 1)
      WRITE(UNIT = dump_unit) SIZE(b_par1, 2)
      WRITE(UNIT = dump_unit) b_par1
      CLOSE(UNIT = dump_unit)
    END IF

    IF (output_b_par2) THEN
      WRITE(b_par2_file , "(A6, I3.3, A4)") "b_par2" , file_number, ".dat"
      OPEN(UNIT = dump_unit, FILE = "Data/" // b_par2_file, FORM = "UNFORMATTED",&
           STATUS = "REPLACE", ACTION = "WRITE")
      WRITE(UNIT = dump_unit) SIZE(b_par2, 1)
      WRITE(UNIT = dump_unit) SIZE(b_par2, 2)
      WRITE(UNIT = dump_unit) b_par2
      CLOSE(UNIT = dump_unit)
    END IF

    WRITE(step_file , "(A4, I3.3, A4)") "step" , file_number, ".dat"
    OPEN(UNIT = dump_unit, FILE = "Data/" // step_file, FORM = "UNFORMATTED",&
         STATUS = "REPLACE", ACTION = "WRITE")
    WRITE(UNIT = dump_unit) step
    CLOSE(UNIT = dump_unit)

    file_number = file_number + 1
    time_dump = time_dump + dt_snapshots

  END SUBROUTINE output_dump

END MODULE diagnostics
