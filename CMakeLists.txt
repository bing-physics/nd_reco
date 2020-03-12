cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B1)
find_package(EDepSim)
find_package(ROOT REQUIRED COMPONENTS EG)
include("${ROOT_USE_FILE}")

include_directories(${PROJECT_SOURCE_DIR})
add_executable(smearit smear.C)
target_link_libraries(smearit PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(sandsmearit sandsmear.C)
target_link_libraries(sandsmearit PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})
