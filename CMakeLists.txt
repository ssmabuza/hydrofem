cmake_minimum_required(VERSION 3.10)

project(hydrofem)

include(FetchContent)
FetchContent_Declare(googletest GIT_REPOSITORY https://github.com/google/googletest.git GIT_TAG v1.14.0)
FetchContent_MakeAvailable(googletest)

set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_STANDARD 23)
set(LINKER_LANGUAGE CXX)

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

include_directories(${MAIN}
                    ${EIGEN} ${FE} ${MESH} ${UTILS} ${DRIVER_PHYSICS}
                    ${LIN_ALG} ${DRIVER} ${DRIVER_AFC_STAB})

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
target_link_libraries(hydrofem.exe PRIVATE boost_program_options)

# unit test using gtest
add_executable(unit_tests ${UTESTS}/Hydrofem_gtests.cpp)
target_link_libraries(unit_tests hydrofem gtest_main gmock_main)

