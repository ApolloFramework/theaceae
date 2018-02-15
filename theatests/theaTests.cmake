#   Functions for actually invoking the tests
include(CMakeParseArguments)

# Need to generalize at a higher level based on mpi, but punt for now
set (MINUS_N "-n")
set (NIM_TIMEOUT 60)
set (MPIEXEC "mpiexec")
if (ENABLE_PARALLEL)
set (MPIEXECCMD "${MPIEXEC} ${MINUS_N}")
set (NPROCS 2)
else()
set (MPIEXECCMD "")
set (NPROCS "")
endif()

# Convenience
set (NIMTESTDIR ${CMAKE_SOURCE_DIR}/nimtests)
set (BINDIR ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------
# addNimTest
# 
# Wrapper around the add_test to handle it's crappy behavior
# 
# Args:
#   TESTNAME:  target for test
#   TEST:      What to execute.  Usually ${TESTNAME}.sh

function(addNimTest)
  set(opts DEBUG) # no-value args
  set(oneValArgs TESTNAME TEST)
  set(multValArgs ) # e.g., lists
  cmake_parse_arguments(FD "${opts}" "${oneValArgs}" "${multValArgs}" ${ARGN})

  message(STATUS "  Add test ${TESTNAME}")
  #TODO  Add in the permutations with cmake options
  add_test(${FD_TESTNAME} ${FD_TEST})
  set_tests_properties(${FD_TESTNAME} PROPERTIES TIMEOUT ${NIM_TIMEOUT})
endfunction()


#----------------------------------------------------------------
# addNimExecTest
# 
# Add executable AND test.
#   - The testname is determined by the first source file name
#   - See addNimTest for just a pure test
# 
# Args:
#   SOURCEFILES:  List of source files.  
#   LINK_LIBS:    Libraries the test needs

function(addNimUnitTest)
  # Parse out the args
  set(opts DEBUG) # no-value args
  set(oneValArgs)
  set(multValArgs SOURCEFILES LINK_LIBS) # e.g., lists
  cmake_parse_arguments(FD "${opts}" "${oneValArgs}" "${multValArgs}" ${ARGN})

  #------------------------------------------------------
  # EXEC and TESTNAME are caps b/c put into nimtest_template.sh.in
  #
  list(GET FD_SOURCEFILES 0 sfile)  # First file determines names
  get_filename_component(sname ${sfile} NAME_WE)
  set(EXEC ${sname})
  message(STATUS "Adding executable ${EXEC} ${FD_SOURCEFILES}")
  add_executable(${EXEC} ${FD_SOURCEFILES} ${NIMTESTDIR}/nimtest_harness.f90)
  target_link_libraries(${EXEC} ${FD_LINK_LIBS})
  install(TARGETS ${EXEC} 
    DESTINATION share/examples/bin
    PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
    GROUP_READ GROUP_EXECUTE GROUP_WRITE 
    WORLD_READ WORLD_EXECUTE)

  # Handle the shell scripts which drive it
  set(TESTNAME ${sname}_test)
  set(scriptname ${sname}.sh)
  set(EXEC ${BINDIR}/${EXEC}) # Make full path for shell script
  configure_file("${NIMTESTDIR}/nim_unit_test.sh.in" "${scriptname}" @ONLY)
  install(FILES ${BINDIR}/${scriptname} 
    DESTINATION share/examples/bin
    PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE
                GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_READ)

  # Now add the test
  addNimTest(TESTNAME ${TESTNAME} TEST ${scriptname})
endfunction()

#----------------------------------------------------------------
# addNimrodTest
# 
# Add test that invokes NIMROD with namelist file
#   - The testname is determined by the first source file name
#   - See addNimTest for just a pure test
# 
# Args:
#   SOURCEFILES:  List of source files.  
#   LINK_LIBS:    Libraries the test needs

function(addNimrodTest)
  # Parse out the args
  set(opts DEBUG) # no-value args
  set(oneValArgs TESTNAME NL_FILE NL_MODS)
  set(multValArgs ) # e.g., lists
  cmake_parse_arguments(FD "${opts}" "${oneValArgs}" "${multValArgs}" ${ARGN})

  #------------------------------------------------------
  # Variables in capitals are in configure'd _file.
  #
  set(EXEC ${BINDIR}/nimrod) # Make full path for shell script
  set(TESTNAME ${FD_TESTNAME})
  set(scriptname ${FD_TESTNAME}.sh)
  set(NL_FILE ${FD_NL_FILE})
  set(NL_MODS ${FD_NL_MODS})
  configure_file("${NIMTESTDIR}/${TEMPLATE_NAME}" "${scriptname}" @ONLY)
  install(FILES ${BINDIR}/${scriptname} 
    DESTINATION share/examples/bin
    PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE
                GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_READ)

  # Now add the test
  addNimTest(TESTNAME ${TESTNAME} TEST ${scriptname})
endfunction()
