cmake_minimum_required(VERSION 3.10)

project(hydrofem)

include(FetchContent)
FetchContent_Declare(googletest GIT_REPOSITORY https://github.com/google/googletest.git GIT_TAG v1.14.0)
FetchContent_MakeAvailable(googletest)

set(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
set(CMAKE_CXX_STANDARD 20)
set(LINKER_LANGUAGE CXX)

if (${HYDROFEM_USE_MFEM} STREQUAL "ON")
    if (NOT MFEM_DIR)
        message(FATAL_ERROR "MFEM_DIR is not set. Please set it to the path of your MFEM installation.")
    endif()
    find_package(MFEM REQUIRED)
    if (NOT MFEM_FOUND)
        message(FATAL_ERROR "MFEM not found. Please install MFEM or set MFEM_DIR to the correct path.")
    endif()
endif()

set(CMAKE_CXX_FLAGS "-fopenmp -O3 -Wall -Wextra -Wshadow")

set(MAIN "./main/")
set(UTESTS "./tests/")

set(FE "./disc_fe/")
set(MESH "./mesh/")
set(UTILS "./utils/")
set(LIN_ALG "./linear/")

set(DRIVER "./driver/")
set(DRIVER_AFC_STAB "${DRIVER}afc_stab/")
set(DRIVER_PHYSICS "${DRIVER}physics/")

# TPLs for the project
set(HYDRFEM_TPLS "./tpls")
set(EIGEN "${HYDRFEM_TPLS}/eigen/")
if (MFEM_FOUND)
    include_directories(${MAIN}
                        ${EIGEN} ${FE} ${MESH} ${UTILS} ${DRIVER_PHYSICS}
                        ${LIN_ALG} ${DRIVER} ${DRIVER_AFC_STAB} ${MFEM_INCLUDE_DIRS})
else()
    include_directories(${MAIN}
                        ${EIGEN} ${FE} ${MESH} ${UTILS} ${DRIVER_PHYSICS}
                        ${LIN_ALG} ${DRIVER} ${DRIVER_AFC_STAB})
endif()

file(GLOB SOURCES "${FE}*.cpp"
                  "${DRIVER}*.cpp"
                  "${DRIVER_PHYSICS}*.cpp"
                  "${DRIVER_AFC_STAB}*.cpp"
                  "${LIN_ALG}*.cpp"
                  "${MESH}*.cpp"
                  "${UTILS}*.cpp")

# create a library from the code
add_library(hydrofem STATIC ${SOURCES})

# create the main executable for the code
add_executable(hydrofem.exe ${SOURCES} ${MAIN}/Hydrofem_Main.cpp)
if (NOT MFEM_FOUND)
    target_link_libraries(hydrofem.exe PRIVATE boost_program_options)
else()
    target_link_libraries(hydrofem.exe PRIVATE boost_program_options ${MFEM_LIBRARIES} ${MFEM_EXTRA_LIBRARIES} mpi_cxx mpi)
endif()

# unit test using gtest
add_executable(unit_tests ${UTESTS}/Hydrofem_gtests.cpp)
target_link_libraries(unit_tests hydrofem gtest_main gmock_main)

