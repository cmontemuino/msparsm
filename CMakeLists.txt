cmake_minimum_required(VERSION 3.5.1)
project(msparsm)

# Require MPI for this project:
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(MPI_COMPILE_FLAGS "-O3 -std=gnu99 -I.")
set(SOURCE_FILES
        ms.c
        ms.h
        mspar.c
        mspar.h
        rand1.c
        streec.c)

add_executable(msparsm ${SOURCE_FILES})
target_link_libraries(msparsm ${MPI_LIBRARIES} -lm)

if(MPI_COMPILE_FLAGS)
    set_target_properties(msparsm PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(msparsm PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

install(TARGETS msparsm DESTINATION ${CMAKE_INSTALL_PREFIX})
