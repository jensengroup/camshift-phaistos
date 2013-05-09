# If any of the following directories exist, add them as subdirectories
foreach(module_subdir bin test unit_tests)  
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${module_subdir}/CMakeLists.txt)

    message(STATUS "    Adding ${module_subdir}")

    # Exclude test directory from ALL (not built on default make)
    if (${module_subdir} STREQUAL test)
      add_subdirectory(${module_subdir} EXCLUDE_FROM_ALL)
    else (${module_subdir} STREQUAL test)
      add_subdirectory(${module_subdir})
    endif(${module_subdir} STREQUAL test)

  endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${module_subdir}/CMakeLists.txt)
endforeach(module_subdir)  

