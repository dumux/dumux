dune_symlink_to_source_files(FILES test_py_gridgeometry.py)
dune_python_add_test(NAME test_py_gridgeometry
                     COMMAND ${PYTHON_EXECUTABLE} test_py_gridgeometry.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python unit)

dune_symlink_to_source_files(FILES test_py_fvproblem.py)
dune_python_add_test(NAME test_py_fvproblem
                     COMMAND ${PYTHON_EXECUTABLE} test_py_fvproblem.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python unit)
