dumux_add_test(SOURCES test_newton.cc
               LABELS unit nonlinear)
dumux_add_test(NAME test_newton_linesearch
               TARGET test_newton
               COMMAND test_newton
               CMD_ARGS "-Newton.UseLineSearch" "true"
               LABELS unit nonlinear)
