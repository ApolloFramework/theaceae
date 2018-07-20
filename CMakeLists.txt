#---------------------------------------------------------------------
#
# $Id: CMakeLists.txt 5633 2018-02-07 21:56:07Z kruger $
# CMakeLists.txt for Theaceae
#
#---------------------------------------------------------------------
# Project information
project(theaceae)
set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
set(THEACEAE_VERSION "${PROJECT_VERSION}-${PROJECT_REV}")

# Required version
cmake_minimum_required(VERSION 2.8.6)
enable_testing()
cmake_policy(SET CMP0004 OLD)

# Fix behavior of CMAKE_CXX_STANDARD when targeting macOS.
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()

#---------------------------------------------------------------------
# Set default values if not otherwise specified.
#---------------------------------------------------------------------
set(ENABLE_DOCS True CACHE BOOL "Whether to build the documentation")

#-----------------------------------------------------------
# Usual initialization stuff
#-----------------------------------------------------------
set(SRCDIR ${CMAKE_SOURCE_DIR})
set(BLDDIR ${CMAKE_BINARY_DIR})
set (CMAKE_CXX_STANDARD 11)
# Various ways of using scimake
if (EXISTS ${CMAKE_SOURCE_DIR}/scimake)
   set(SCIMAKE_DIR ${CMAKE_SOURCE_DIR}/scimake)
elseif (EXISTS ${CMAKE_SOURCE_DIR}/../scimake)
   set(SCIMAKE_DIR ${CMAKE_SOURCE_DIR}/../scimake)
else ()
   # This is useful for spack
   set(SCIMAKE_DIR $ENV{APOLLO_ROOT}/scimake)
endif()
include(${SCIMAKE_DIR}/SciInit.cmake)

# Set permissions before adding subdirectories
set(SCI_GROUP_WRITE GROUP_WRITE)
set(SCI_WORLD_FILE_PERMS WORLD_READ)
set(SCI_WORLD_PROGRAM_PERMS WORLD_READ WORLD_EXECUTE)

# Parse input
if (NOT DEFINED CMAKE_COLOR_MAKEFILE)
  set(CMAKE_COLOR_MAKEFILE TRUE)
endif ()

#-----------------------------------------------------------
# Gone throught the examples, and lib -- worry about documentation
# Documentation is done in the build directory so have to copy
# over the rst files
#-----------------------------------------------------------
# Sphinx/MathJax packages. 
#-----------------------------------------------------------
if (ENABLE_DOCS)
  find_package(SciSphinx)
  include(${SCIMAKE_DIR}/SciDocLogic.cmake)
  include(${SCIMAKE_DIR}/Modules/SciSphinxFunctions.cmake)

  find_package(SciMathJax REQUIRED)
  set(MATHJAXJS "MathJax/MathJax.js?config=TeX-AMS-MML_SVG-full.js")


  #---------------------
  # This does a lot of file copying and stuff to setup for the sphinx
  # Includes running python script for getting the examples
  # There is no add_subdirectory for the docs because it is not needed
  #include(theadocs/docs.cmake)

  #---------------------
  # Create documentation 
  #---------------------
  set(RST_FILE_BASE "Theaceae")
  set(RST_AUTHOR "ApolloFramework Team")
  set(RST_TITLE  "Theaceae Documentation")
  set(RST_MAIN_FILE ${RST_FILE_BASE}) 

  # Copy stuff over from docs because sphinxbuild doesn't allow full paths
  configure_file(${SRCDIR}/theadocs/conf.py.in ${BLDDIR}/conf.py)  # 1 copy
  configure_file(${SRCDIR}/theadocs/build_index.html.in ${BLDDIR}/build_index.html)  # 1 copy

  # Pull together all ReST files for reference manual
  file(GLOB_RECURSE THEACEAE_DEPS RELATIVE ${SRCDIR} *.rst)

  # This sets the targets for each build (html, pdf, ...)
  SciSphinxTarget(TARGET userdocs
                  RST_FILE_BASE ${RST_FILE_BASE}
                  FILE_DEPS ${THEACEAE_DEPS}
                  INSTALL_SUPERDIR doc
                  SPHINX_DOCTREE_DIR "-d${BLDDIR}/doctrees"
                  SPHINX_ADDL_OPTS ${SPHINX_MATHARG}
  )


  # Add top level targets
  add_custom_target(userdocs ALL DEPENDS userdocs-html)
endif (ENABLE_DOCS)

#----------------------------------------------------------
#   Dependencies
#----------------------------------------------------------
# First get the config file
SciFindPackage(PACKAGE "Camellia"
  INSTALL_DIR par
  FIND_CONFIG_FILE
  USE_CONFIG_FILE
  CONFIG_SUBDIRS lib/Camellia/cmake
)   

include_directories(${Camellia_INCLUDE_DIRS} ${Camellia_TPL_INCLUDE_DIRS})

#---------------------------------------------------------------------
#  Add subdirectories -- 
#---------------------------------------------------------------------
add_subdirectory(theatests)
add_subdirectory(fluid_dynamics)

#---------------------------------------------------------------------
# Finalize (installation, packaging)
#---------------------------------------------------------------------
set(COMMON_INSTALL TRUE)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Theaceae: Collection of test cases to evaluate Petrov-Galerkin and other methods")
set(INSTALL_CONFIG_HEADERS FALSE)
include(${SCIMAKE_DIR}/SciFinalize.cmake)
