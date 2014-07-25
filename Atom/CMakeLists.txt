# Copyright (c) 2010-2014, Delft University of Technology
# Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
# All rights reserved.
# See http://bit.ly/12SHPLR for license details.

# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(Atom)

# Define the directory with the source code.
set(SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the project root directory.
set(PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Set path to local CMake modules.
set(CMAKE_MODULE_PATH ${PROJECT_ROOT}/CMakeModules)

# Set testing options based on platform.
enable_testing()

# Set platform-specific compiler flags.
if (WIN32)
  if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
    # Because we are using static boost libraries, with static runtime, we need to force MSVC to
    # also use static runtime: (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace).
    foreach(flag_var
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      # Find all dynamic runtime (MD) references and replace with static (MT)
      if(${flag_var} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
  else()
    set(CMAKE_CXX_FLAGS 
        "${CMAKE_CXX_FLAGS} -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  endif()
else()
  set(CMAKE_CXX_FLAGS 
      "${CMAKE_CXX_FLAGS} -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
endif()

# -------------------------------------------------------------------------------------------------
# DO NOT MANUALLY EDIT ANY SETTINGS AFTER THIS POINT IN THE FILE.
# 
# Pass any additional settings directly to CMake (e.g., via command-line).
# -------------------------------------------------------------------------------------------------

# Set CMake build-type. If it not supplied by the user, the default built type is "Release". 
if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(STATUS "WARNING: building release version!")
    set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WARNING: building debug version!")
endif()

# Set lib and bin directories where static libraries, executable, unit tests are built.
if(NOT LIB_ROOT)
  set(LIB_ROOT "${PROJECT_BINARY_DIR}/../lib")
endif()
if(NOT BIN_ROOT)
  set(BIN_ROOT "${PROJECT_BINARY_DIR}/../bin")
endif()
if(NOT TEST_ROOT)
  set(TEST_ROOT "${PROJECT_BINARY_DIR}/../tests")
endif()

# Set the global macros for setting up targets.
macro(setup_executable_target target_name)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BIN_ROOT}")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BIN_ROOT}")
endmacro(setup_executable_target)

macro(setup_library_target target_name)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIB_ROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIB_ROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_target target_name)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${TEST_ROOT}")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${TEST_ROOT}/${target_name}")
endmacro(setup_unit_test_target)

# Include the top-level directories.
include_directories(AFTER
  "${PROJECT_ROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif( )

# Configure Boost libraries.
if(NOT Boost_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
endif()
if(NOT Boost_USE_MULTITHREADED)
  set(Boost_USE_MULTITHREADED ON)
endif()
if(NOT Boost_USE_STATIC_RUNTIME)
  set(Boost_USE_STATIC_RUNTIME ON)
endif()

# Find Boost libraries on local system.
find_package(Boost 1.52.0 
             COMPONENTS 
             thread date_time system unit_test_framework filesystem regex serialization chrono
             timer 
             REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif( )

# Find GSL library on local system.
find_package(GSL REQUIRED)

# Include GSL directories.
# Set CMake flag to suppress GSL warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${GSL_INCLUDE_DIRS}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${GSL_INCLUDE_DIRS}\"")
endif( )
# Find SGP4DanRW library on local system.
find_package(SGP4DanRW 1 REQUIRED)

# Include SGP4DanRW directories.
# Set CMake flag to suppress SGP4DanRW warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${SGP4DANRW_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SGP4DANRW_INCLUDE_DIR}\"")
endif( )

# If documentation should be built, find Doxygen package and setup config file.
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)

  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it correctly.")
  endif()

  # Configure the Template Doxyfile for our specific project.
  configure_file(${PROJECT_ROOT}/Doxyfile.in 
                 ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

  # Add a custom target to run Doxygen when ever the project is built.
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.
  add_custom_target (docs ALL COMMAND ${DOXYGEN_EXECUTABLE} 
                     ${PROJECT_BINARY_DIR}/Doxyfile 
                     SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
endif()

# Set flags to profile using gprof.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

# Add source files for Atom static library.
set(ATOM_SOURCES
  "${SOURCE_ROOT}/convertCartesianStateToTwoLineElements.cpp"
)

# Add Atom static library.
add_library(atom STATIC ${ATOM_SOURCES})
setup_library_target(atom)

# Add Atom application.
add_executable(atom_main "${SOURCE_ROOT}/main.cpp")
setup_executable_target(atom_main)
target_link_libraries(atom_main
                      ${GSL_LIBRARIES}
                      ${SGP4DANRW_LIBRARY}
                      ${TUDAT_CORE_LIBRARIES} 
                      ${Boost_LIBRARIES})
