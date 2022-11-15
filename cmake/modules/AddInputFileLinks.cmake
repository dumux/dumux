# Creates symbolic links to all input files in the source directory
include_guard(GLOBAL)

macro(add_input_file_links)
  FILE(GLOB input_files *.input)
  foreach(VAR ${input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()

macro(add_python_file_links)
  FILE(GLOB input_files *.py)
  foreach(VAR ${input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
