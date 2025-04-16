cmake_minimum_required(VERSION 3.10)

project(hydrofem)

include(CTest)

include(FetchContent)
FetchContent_Declare(googletest
					 GIT_REPOSITORY git@github.com:google/googletest.git
					 GIT_TAG v1.14.0)

FetchContent_GetProperties(googletest)

if (NOT googletest_POPULATED)
	FetchContent_Populate(googletest)
	add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BUILD_DIR})
endif()

set(CMAKE_CXX_STANDARD 23)
set(LINKER_LANGUAGE CXX)

set(CMAKE_CXX_FLAGS "-fopenmp -O3 -Wall -Wextra -Wshadow")
set(HOME    "~/"                      )

set(MAIN    "./main/"                 )
set(UTESTS    "./tests/"                 )

set(FE      "./disc_fe/"              )
set(MESH    "./mesh/"                 )
set(UTILS   "./utils/"                )
set(LIN_ALG "./linear/"       )

set(DRIVER          "./driver/"         )
set(DRIVER_AFC_STAB "./driver/afc_stab/")
set(PHYSICS "./driver/physics/"         )

set(HYDRFEM_TPLS "./tpls")

set(EIGEN "${HYDRFEM_TPLS}/eigen/")

set(CLI11_INCLUDE_DIR "${HYDRFEM_TPLS}/CLI11/include")

include_directories(${MAIN}
                    ${EIGEN} ${FE} ${MESH} ${UTILS} ${PHYSICS}
                    ${LIN_ALG} ${DRIVER} ${DRIVER_AFC_STAB} ${CLI11_INCLUDE_DIR})

file(GLOB SOURCES "${FE}*.cpp"
                  "${FE}*.hpp"
                  "${DRIVER}*.cpp"
                  "${DRIVER}*.hpp"
                  "${PHYSICS}*.cpp"
                  "${PHYSICS}*.hpp"
                  "${DRIVER_AFC_STAB}*.cpp"
                  "${DRIVER_AFC_STAB}*.hpp"
                  "${LIN_ALG}*.cpp"
                  "${LIN_ALG}*.hpp"
                  "${MESH}*.cpp"
                  "${MESH}*.hpp"
                  "${UTILS}*.cpp"
                  "${UTILS}*.hpp")

# create a library from the code
add_library(libhydrofem STATIC ${SOURCES})

# create the main executable for the code
add_executable(hydrofem ${SOURCES} ${MAIN}/Hydrofem_Main.cpp)
target_link_libraries(hydrofem PRIVATE boost_program_options)

# unit test using gtest
add_executable(unit_tests ${UTESTS}/utests.cpp)
target_link_libraries(unit_tests libhydrofem gtest_main gmock_main)


add_executable(hydrofem_default_test ${UTESTS}/Hydrofem_CTest_Default.cpp)
target_link_libraries(hydrofem_default_test PRIVATE libhydrofem)

enable_testing()

# add the test to ctest
add_test(NAME default_test
         COMMAND hydrofem_default_test
         )