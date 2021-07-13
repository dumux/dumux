add_subdirectory(assembly)
add_subdirectory(common)
add_subdirectory(discretization)
add_subdirectory(material)
add_subdirectory(io)
add_subdirectory(porousmediumflow)

add_python_targets(dumux
  __init__
)
