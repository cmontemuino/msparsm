cmake_minimum_required(VERSION 2.6s)
project(workbench_mpi3)

set(SOURCE_FILES main.c)
add_executable(workbench_mpi3 ${SOURCE_FILES})

# Require MPI for this project:
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(workbench_mpi3 ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(workbench_mpi3 PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(workbench_mpi3 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

# Add a test:
enable_testing()
add_test(SimpleTest ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${CMAKE_CURRENT_BINARY_DIR}/workbench_mpi3)