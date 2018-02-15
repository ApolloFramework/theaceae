######################################################################
#
# CMakeLists.txt for sphinx docs
#
# $Id: CMakeLists.txt 894 2016-06-29 23:28:32Z kruger $
#
######################################################################

set(iconfiles
  ${CMAKE_SOURCE_DIR}/ptsolve.ico
  ${CMAKE_SOURCE_DIR}/ptsolve.png 

)
foreach (iconfile ${iconfiles})
  add_custom_command(OUTPUT ${iconfile}
    COMMAND $(CMAKE_COMMAND) -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${iconfile}"
      "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${iconfile}"
  )
endforeach ()
add_custom_target(copyiconfiles DEPENDS ${iconfiles})
# Optional installation works without "make userdocs" being done
install(FILES ${iconfiles}
  DESTINATION ${USERDOCS_INSTALLDIR}/_images
  OPTIONAL
)
add_subdirectory(templates)

