include_guard(GLOBAL)

# set variable for config.h
set(HAVE_KOKKOS ${Kokkos_FOUND})

# perform DUNE-specific setup tasks
if (Kokkos_FOUND)
  dune_register_package_flags(
    COMPILE_DEFINITIONS ENABLE_KOKKOS=1
    LIBRARIES Kokkos::kokkos
  )
endif()
