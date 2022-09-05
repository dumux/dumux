if (${HAVE_EIGEN3})
    dumux_add_test(
        SOURCES test_poisson_tpfa.cc
        LABELS unit new_assembly
        LINK_LIBRARIES Eigen3::Eigen
    )
else()
    dumux_add_test(
        SOURCES test_poisson_tpfa.cc
        LABELS unit new_assembly
    )
endif ()
