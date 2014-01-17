# Author: petter.strandmark@gmail.com (Petter Strandmark)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4)
PROJECT(Submodular C CXX)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Default locations to search for on various platforms.
LIST(APPEND SEARCH_LIBS /usr/lib)
LIST(APPEND SEARCH_LIBS /usr/local/lib)
LIST(APPEND SEARCH_LIBS /usr/local/homebrew/lib) # Mac OS X
LIST(APPEND SEARCH_LIBS /opt/local/lib)

LIST(APPEND SEARCH_HEADERS /usr/include)
LIST(APPEND SEARCH_HEADERS /usr/local/include)
LIST(APPEND SEARCH_HEADERS /usr/local/homebrew/include) # Mac OS X
LIST(APPEND SEARCH_HEADERS /opt/local/include)

ENABLE_TESTING()
SET(CTEST_OUTPUT_ON_FAILURE TRUE)

# Multithreading using OpenMP
OPTION(OPENMP
       "Enable multi-threading (requires OpenMP)"
       ON)

# For Clang-based tools.
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF (${OPENMP})
  FIND_PACKAGE(OpenMP)
  IF(${OPENMP_FOUND})
    MESSAGE("-- Found OpenMP.")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    ADD_DEFINITIONS(-DUSE_OPENMP)
  ELSE(${OPENMP_FOUND})
    MESSAGE("-- Can't find OpenMP. Continuing without it.")
  ENDIF(${OPENMP_FOUND})
ELSE (${OPENMP})
  MESSAGE("-- OpenMP is disabled. Continuing without it.")
ENDIF (${OPENMP})


FILE(GLOB SUBMODULAR_HEADERS ${CMAKE_SOURCE_DIR}/source/library/*.h)
INSTALL(FILES ${SUBMODULAR_HEADERS} DESTINATION include)

# C++11 support.
include(cmake/EnableCPP11.cmake)

#
# Libraries
#
SET(COIN_LIB_DIR "" CACHE PATH "Coin (Clp, ...) library directory")

FIND_LIBRARY(CLP_LIB NAMES Clp libClp PATHS ${COIN_LIB_DIR})
FIND_LIBRARY(COINUTILS_LIB NAMES CoinUtils libCoinUtils PATHS ${COIN_LIB_DIR})

IF (NOT EXISTS ${CLP_LIB})
	MESSAGE(FATAL_ERROR
	        "Can't find Clp. Please specify COIN_LIB_DIR")
ELSE (NOT EXISTS ${CLP_LIB})
	MESSAGE("-- Found Coin library (Clp).")
ENDIF (NOT EXISTS ${CLP_LIB})

# Library dependencies
SET (SUBMODULAR_LIBRARY_DEPENDENCIES)

# OpenMP needs to be included as a library on some platforms.
IF (${OPENMP_FOUND})
  IF (NOT MSVC)
    LIST(APPEND SUBMODULAR_LIBRARY_DEPENDENCIES gomp)
  ENDIF (NOT MSVC)
ENDIF (${OPENMP_FOUND})

LIST(APPEND SUBMODULAR_LIBRARY_DEPENDENCIES ${CLP_LIB} ${COINUTILS_LIB})
LIST(APPEND SUBMODULAR_LIBRARY_DEPENDENCIES QPBO)
LIST(APPEND SUBMODULAR_LIBRARY_DEPENDENCIES maxflow3)
LIST(APPEND SUBMODULAR_LIBRARY_DEPENDENCIES petter)

#
# Include directories
#
INCLUDE_DIRECTORIES(
  source/library
  thirdparty/Petter
  )

FIND_PATH(HOCR_INCLUDE_DIR NAMES HOCR.h PATHS thirdparty/HOCR)
IF (NOT EXISTS ${HOCR_INCLUDE_DIR})
	MESSAGE(FATAL_ERROR
	        "Can't find HOCR library. See thirdparty directory on how to download and install.")
ELSE (NOT EXISTS ${HOCR_INCLUDE_DIR})
	MESSAGE("-- Found HOCR include directory.")
	INCLUDE_DIRECTORIES(${HOCR_INCLUDE_DIR})
ENDIF (NOT EXISTS ${HOCR_INCLUDE_DIR})

FIND_PATH(QPBO_INCLUDE_DIR NAMES QPBO.h PATHS thirdparty/QPBO)
IF (NOT EXISTS ${QPBO_INCLUDE_DIR})
	MESSAGE(FATAL_ERROR
	        "Can't find QPBO library. See thirdparty directory on how to download and install.")
ELSE (NOT EXISTS ${QPBO_INCLUDE_DIR})
	MESSAGE("-- Found QPBO include directory.")
	INCLUDE_DIRECTORIES(${QPBO_INCLUDE_DIR})
ENDIF (NOT EXISTS ${QPBO_INCLUDE_DIR})

FIND_PATH(MAXFLOW_INCLUDE_DIR NAMES maxflow.cpp PATHS thirdparty/maxflow-v3.01.src)
IF (NOT EXISTS ${MAXFLOW_INCLUDE_DIR})
	MESSAGE(FATAL_ERROR
	        "Can't find maxflow-v3.01.src library. See thirdparty directory on how to download and install.")
ELSE (NOT EXISTS ${MAXFLOW_INCLUDE_DIR})
	MESSAGE("-- Found maxflow-v3.01.src include directory.")
	INCLUDE_DIRECTORIES(${MAXFLOW_INCLUDE_DIR})
ENDIF (NOT EXISTS ${MAXFLOW_INCLUDE_DIR})

SET(COIN_INCLUDE_DIR "" CACHE PATH "Coin (Clp, ...) include directory (where coin/ is)")

FIND_PATH(COIN_INCLUDE NAMES coin/ClpSimplex.hpp PATHS ${COIN_INCLUDE_DIR})
IF (NOT EXISTS ${COIN_INCLUDE})
	MESSAGE(FATAL_ERROR
	        "Can't find coin include directory. Please specify COIN_INCLUDE_DIR")
ELSE (NOT EXISTS ${COIN_INCLUDE})
	MESSAGE("-- Found Coin include directory.")
ENDIF (NOT EXISTS ${COIN_INCLUDE})

INCLUDE_DIRECTORIES(${COIN_INCLUDE})

# Change the default build type from Debug to Release, while still
# supporting overriding the build type.
IF (NOT CMAKE_BUILD_TYPE)
  MESSAGE("-- No build type specified; defaulting to CMAKE_BUILD_TYPE=Release.")
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

# C++11 support.
#include(cmake/EnableCPP11.cmake)
#include(cmake/CheckGenericLambdas.cmake)

#
# Clang and GCC settings
#
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR CMAKE_COMPILER_IS_GNUCXX)
  # No warnings for C. The Meschach library contains really old code.
  # -fPIC for building a shared library.
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -w -fpic")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-compare -Wno-unused-parameter -Wno-unused-function -fpic")
  # -Wunused-function triggers incorrect warnings for variadic template
  # recursion base case with Clang 3.2.

  IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
  	MESSAGE("-- Debug mode enabled for Gcc/Clang; adding support for Gcov.")
  	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
  	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  ENDIF (CMAKE_BUILD_TYPE STREQUAL "Debug")

ENDIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR CMAKE_COMPILER_IS_GNUCXX)

#
# MSVC settings
#
IF (MSVC)
  # Disable deprecation warning for standard functions.
  ADD_DEFINITIONS("/wd4996")
  # To get rid of annoying min and max macros if windows.h
  # is included.
  ADD_DEFINITIONS("-DNOMINMAX=1")
  # Macro redefinitions.
  ADD_DEFINITIONS("/wd4005")
ENDIF (MSVC)

#
# Packaging
#
SET(CPACK_GENERATOR ZIP)
SET(CPACK_PACKAGE_FILE_NAME "Submodular")
INCLUDE(CPack)

ADD_SUBDIRECTORY(bin)
ADD_SUBDIRECTORY(bin/tests)
ADD_SUBDIRECTORY(thirdparty)
ADD_SUBDIRECTORY(source/library)
ADD_SUBDIRECTORY(source)

