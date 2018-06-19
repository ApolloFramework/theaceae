cmake_minimum_required (VERSION 2.6)
project (SA2DFlatPlate)
# The version numer
set (SA2DFlatPlate_VERSION_MAJOR 1)
set (SA2DFlatPlate_VERSION_MINOR 0)


if (EXTRA_INCLUDE)
  include_directories(${EXTRA_INCLUDE})
endif()

if (OMP_INCLUDE)
  include_directories(${OMP_INCLUDE})
endif()

if(MPI_INCLUDE)
  include_directories(${MPI_INCLUDE})
endif()
if(MPI_LINK_CC)
  link_libraries(${MPI_LINK_CC})
endif()
if (MPI_LINK_FORTRAN)
  link_libraries(${MPI_LINK_FORTRAN})
endif()

if(MPI_DIR)
  MESSAGE("MPI_DIR = ${MPI_DIR}")
  include_directories(${MPI_DIR}/include)
  link_directories(${MPI_DIR})
  link_libraries(mpi)
  # use with openmpi
  link_libraries(mpi_mpifh)
  # use with mpich / intel mpi
  # link_libraries(mpifort) 
endif()

if(SCALAPACK_LIB)
  link_libraries(${SCALAPACK_LIB})
endif()

if(CAMELLIA_DIR)
  MESSAGE("CAMELLIA_DIR = ${CAMELLIA_DIR}")
  include_directories(${CAMELLIA_DIR}/include)
  link_libraries(${CAMELLIA_LIB})
endif()

if(ADDITIONAL_LIBRARY_DIRS)
  MESSAGE("ADDITIONAL_LIBRARY_DIRS: ${ADDITIONAL_LIBRARY_DIRS}")
  MESSAGE("ADDITIONAL_LIBRARIES: ${ADDITIONAL_LIBRARIES}")
  link_directories(${ADDITIONAL_LIBRARY_DIRS})
  set(Trilinos_TPL_LIBRARIES ${Trilinos_TPL_LIBRARIES}
${ADDITIONAL_LIBRARIES})
endif()

MESSAGE(" Looking for Trilinos in ${TRILINOS_PATH}/lib/cmake/Trilinos
  ${TRILINOS_PATH} $ENV{TRILINOS_DIR}")

# Find Trilinos
find_package(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos
  ${TRILINOS_PATH} $ENV{TRILINOS_DIR})

list(REVERSE Trilinos_TPL_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)

list(REVERSE Trilinos_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
list(REVERSE Trilinos_LIBRARIES)

list(REVERSE Trilinos_TPL_INCLUDE_DIRS)
list(REMOVE_DUPLICATES Trilinos_TPL_INCLUDE_DIRS)
list(REVERSE Trilinos_TPL_INCLUDE_DIRS)

# If FIND_PACKAGE successfully found your Trilinos install, it will
# set the Boolean flag Trilinos_FOUND.  The following IF statement
# fails with a FATAL_ERROR if Trilinos was not found.  If it _was_
# found, it prints out the values of some Trilinos configuration
# details.  You may find them useful for building your application
# that uses Trilinos.
IF(Trilinos_FOUND)
   MESSAGE("\nFound Trilinos!  Here are the details: ")
   MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
   MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
   MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
   MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
   MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
   MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
   MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
   MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
   MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
   MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
   MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS =
     ${Trilinos_Fortran_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
   MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
   MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
   MESSAGE("End of Trilinos details\n")
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos!")
ENDIF()

# If you haven't already set the C++ compiler, use the same compiler
# that was used to build Trilinos on your machine.  
#IF(NOT CMAKE_CXX_COMPILER)
  MESSAGE("Setting CXX compiler to ${Trilinos_CXX_COMPILER}.")
  SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
#ENDIF()

# Testing something:
#SET(CMAKE_CXX_COMPILER "echo")

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

IF(NOT CMAKE_CXX_FLAGS)
  # -fPIC allows shared libraries to be built against the Camellia library
  SET(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})
ENDIF() 

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -g")
MESSAGE("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

MESSAGE (STATUS "CMAKE_C_COMPILE = ${CMAKE_C_COMPILER}")
# If you haven't already set the C compiler, use the same compiler
# that was used to build Trilinos on your machine.  
#IF(NOT CMAKE_C_COMPILER)
  MESSAGE("Setting C compiler.")
  SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
#ENDIF()

#IF (NOT CMAKE_C_FLAGS)
  SET(CMAKE_C_FLAGS ${Trilinos_C_COMPILER_FLAGS})
#ENDIF()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -g")
#MESSAGE("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

# If you haven't already set the Fortran compiler, use the same compiler
# that was used to build Trilinos on your machine.  
IF(NOT CMAKE_Fortran_COMPILER)
  MESSAGE("Setting Fortran compiler.")
  SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
ENDIF()

IF(NOT BUILD_SHARED_LIBS)
  MESSAGE("Setting BUILD_SHARED_LIBS to ${Trilinos_BUILD_SHARED_LIBS}")
  SET(BUILD_SHARED_LIBS ${Trilinos_BUILD_SHARED_LIBS})
ENDIF()

# Find MOAB
#find_package(MOAB PATHS ${MOAB_PATH}/lib)

#if (MOAB_FOUND)
#  ADD_DEFINITIONS(-DHAVE_MOAB)
#  include_directories(
#    ${MOAB_INCLUDE_DIRS}
#  )
#  SET(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${MOAB_LIBRARIES})
#else()
#  MESSAGE("MOAB not found; MOAB_PATH = ${MOAB_PATH}.")
#endif(MOAB_FOUND)

include_directories(
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${Trilinos_INCLUDE_DIRS}
  ${CAMELLIA_SOURCE_DIR}/include
)

link_directories(
  ${Trilinos_LIBRARY_DIRS}
  ${Trilinos_TPL_LIBRARY_DIRS}
)

link_libraries(
 ${Trilinos_LIBRARIES}
 ${Trilinos_TPL_LIBRARIES}
 ${ADDITIONAL_LIBRARIES}
)


add_executable(SAFlatPlate IncompressibleSAFlatPlate.cpp)
