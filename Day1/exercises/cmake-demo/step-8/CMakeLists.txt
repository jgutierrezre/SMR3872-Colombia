# CMake configuration demo
# version 3.10 is bundled with ubuntu18.04
cmake_minimum_required(VERSION 3.10)

# set the project name, version, and primary programming language (C++)
project(pi VERSION 0.1 LANGUAGES C)

# look for MPI C interface
find_package(MPI REQUIRED)

# add the MPI executable
add_executable(pi_mpi pi_mpi.c)
target_link_libraries(pi_mpi PRIVATE MPI::MPI_C)

# look for OpenMP and OpenMP library
find_package(OpenMP REQUIRED)
include(CheckIncludeFile)
check_include_file(omp.h HAVE_OMP_H)
if(!HAVE_OMP_H)
    message(FATAL_ERROR "Must have OpenMP runtime and omp.h header file")
endif()
add_executable(pi_omp pi_omp.c)
target_link_libraries(pi_omp PRIVATE OpenMP::OpenMP_C)
