add_python_targets(assembly
  __init__
)
dune_add_pybind11_module(NAME _assembly)
set_property(TARGET _assembly PROPERTY LINK_LIBRARIES dunecommon dunegrid APPEND)
install(TARGETS _assembly LIBRARY DESTINATION python/dumux/assembly)
