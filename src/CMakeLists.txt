# set version number
set (${CMAKE_PROJECT_NAME}_VERSION_MAJOR 0)
set (${CMAKE_PROJECT_NAME}_VERSION_MINOR 1)
set (${CMAKE_PROJECT_NAME}_VERSION_PATCH 0)
set (PROJECT_VERSION "${${CMAKE_PROJECT_NAME}_VERSION_MAJOR}.${${CMAKE_PROJECT_NAME}_VERSION_MINOR}.${${CMAKE_PROJECT_NAME}_VERSION_PATCH}")

# Included scripts do automatic cmake_policy PUSH and POP.
# https://cmake.org/cmake/help/latest/policy/CMP0011.html
cmake_policy(SET CMP0011 NEW)

# output configuration details
message (STATUS "Cmake version ${CMAKE_VERSION}")
message (STATUS "Configuring mspario v${PROJECT_VERSION}")
message (STATUS "Installation directory: ${CMAKE_INSTALL_PREFIX}")
message (STATUS "Source directory: ${CMAKE_SOURCE_DIR}")
message (STATUS "Include path: ${CMAKE_INCLUDE_PATH}")
message (STATUS "Library path: ${CMAKE_LIBRARY_PATH}")
message (STATUS "Prefix path: ${CMAKE_PREFIX_PATH}")
message (STATUS "Module path: ${CMAKE_MODULE_PATH}")
message (STATUS "Binary path: ${PROJECT_BINARY_DIR}")
message (STATUS "Build type is: ${CMAKE_BUILD_TYPE}")
message (STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Require MPI for this project:

# add include directories
#include_directories(${MPI_INCLUDE_PATH})

# set all source files
set(sources
        ms.c
        ms.h
        mspar.c
        mspar.h
        rand1.c
        streec.c)

# create targets for executable and libraries
add_executable (${CMAKE_PROJECT_NAME} ${sources})

# get mpi and link to all targets
find_package(MPI REQUIRED)
target_link_libraries (${CMAKE_PROJECT_NAME} ${MPI_C_LIBRARIES})

set(MPI_COMPILE_FLAGS "-O2 -std=gnu99 -I.")
set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")

# install executable, libraries, and includes
install(TARGETS ${CMAKE_PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
