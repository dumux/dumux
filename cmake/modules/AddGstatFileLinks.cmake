# Creates symbolic links to all gstat files in the source directory
macro(add_gstat_file_links)
  FILE(GLOB gstat_files gstat*.txt *.gstat)
  foreach(VAR ${gstat_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
