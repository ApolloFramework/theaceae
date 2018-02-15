######################################################################
#
# $Id: CMakeLists.txt 4014 2014-08-18 23:06:22Z jking $
# CMakeLists.txt for theatests
#
# Complete harness is made of:
#    theaTests.cmake, theatest_template.sh.in, and theaharness.sh
#
######################################################################
# Make sure we build locally
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Functions that actually do the tests
include(theaTests.cmake)

include_directories(
  ${CMAKE_BINARY_DIR}
)

# Unit tests are normally in subdirectories.  
# This is a generic unit test that is also good for testing the harness itself
#addTheaUnitTest(SOURCEFILES testSIMD.F90 LINK_LIBS nimblk)

# Useful commands
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
        --force-new-ctest-process --timeout 20 --output-on-failure)


