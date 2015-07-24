# Creates symbolic links to all input files in the source directory
macro(add_input_file_links)
  FILE(GLOB input_files *.input)
  foreach(VAR ${input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
