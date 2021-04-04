PROGRAM leapfrog

  USE shared_data
  USE diagnostics
  USE setup
  USE control
  USE update_variables

  IMPLICIT NONE

  CALL system_clock(tstart, count_rate)

  CALL control_variables
  CALL after_control
  CALL set_initial_conditions
  CALL set_boundary_conditions(CURR_VARIABLES, time)
  CALL set_boundary_conditions(M_VARIABLES, time - dt)
  CALL output_grid

  PRINT*, "n_iters =", n_iters
  CALL output_dump
  DO step = 1, n_iters
    CALL predictor_corrector
    IF (time .GT. time_dump) THEN
      CALL output_dump
    ENDIF
    time = time + dt
  END DO

  CALL finalize

  CALL system_clock(tstop, count_rate)

  PRINT*, "Computation time =", REAL(tstop - tstart) / REAL(count_rate)

END PROGRAM leapfrog