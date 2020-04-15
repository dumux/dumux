dune_python_add_test(NAME test_py_gridgeometry
                     COMMAND ${PYTHON_EXECUTABLE} test_py_gridgeometry.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS python unit)
