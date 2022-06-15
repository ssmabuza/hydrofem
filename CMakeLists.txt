project(hydrofem)

cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(LINKER_LANGUAGE CXX)

set(CMAKE_CXX_FLAGS "-fopenmp -O3 -Wall -Wextra")
set(HOME    "~/"                      )

set(MAIN    "./main/"                 )

set(FE      "./disc_fe/"              )
set(MESH    "./mesh/"                 )
set(UTILS   "./utils/"                )
set(LIN_ALG "./linear_solvers/"       )

set(DRIVER          "./driver/"         )
set(DRIVER_AFC_STAB "./driver/afc_stab/")
set(PHYSICS "./driver/physics/"         )

set(HYDRFEM_TPLS "./tpls")

set(EIGEN "${HYDRFEM_TPLS}/eigen/")

include_directories(${MAIN}
                    ${EIGEN} ${FE} ${MESH} ${UTILS} ${PHYSICS}
                    ${LIN_ALG} ${DRIVER} ${DRIVER_AFC_STAB})

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
                  "${MAIN}*.cpp"
                  "${MAIN}*.hpp"
                  "${MESH}*.cpp"
                  "${MESH}*.hpp"
                  "${UTILS}*.cpp"
                  "${UTILS}*.hpp")

add_executable(hydrofem ${SOURCES})
if (WIN32)
  
else()
  target_link_libraries(hydrofem PRIVATE boost_program_options)
endif(WIN32)
