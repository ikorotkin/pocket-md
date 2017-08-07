cmake_minimum_required(VERSION 2.8)

project(pocket-md)

enable_language(C)

# Set DIR with GROMACS sources here (for XTC output only)
set(HHMD_GROMACS_SOURCES "/home/gromacs-4.6.7" CACHE STRING "DIR with GROMACS 4.6 sources")

# Set DIR with compiled (but not installed) GROMACS here (for XTC output only)
set(HHMD_GROMACS_RELEASE "/home/gromacs-4.6.7_release" CACHE STRING "DIR with compiled (but not installed) GROMACS 4.6")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Reference." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fopenmp -Wall")
set(INCL_DIR "${HHMD_GROMACS_SOURCES}/include")
set(INCL_DIR_REL "${HHMD_GROMACS_RELEASE}/include")
set(LIB_GMX "${HHMD_GROMACS_RELEASE}/src/gmxlib/libgmx.so.8")

FILE(GLOB SOURCES *.c)

include_directories(${INCL_DIR} ${INCL_DIR_REL})

add_library(GMX SHARED IMPORTED)
set_target_properties(GMX PROPERTIES IMPORTED_LOCATION ${LIB_GMX})

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} GMX)

find_library(MATH_LIB m)
target_link_libraries(${PROJECT_NAME} ${MATH_LIB})

