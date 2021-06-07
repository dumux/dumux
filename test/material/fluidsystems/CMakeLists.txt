dune_symlink_to_source_files(FILES "pythonfluidsystem.py")

dumux_add_test(SOURCES test_pythonfluidsystem.cc
               LINK_LIBRARIES ${DUNE_LIBS} ${Python3_LIBRARIES}
               LABELS unit material)

dumux_add_test(SOURCES test_fluidsystems.cc
              LABELS unit material)
