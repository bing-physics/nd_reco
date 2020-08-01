cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B1)
set(CMAKE_BUILD_TYPE Debug)
find_package(EDepSim)
find_package(ROOT REQUIRED COMPONENTS EG)
include("${ROOT_USE_FILE}")

include_directories(${PROJECT_SOURCE_DIR})
add_executable(smearit smear.C)
target_link_libraries(smearit PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

#add_executable(sandsmearit sandsmear.C)
#target_link_libraries(sandsmearit PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})
				   
add_executable(smearNReso smear_neutron3case2.C)
target_link_libraries(smearNReso PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(smearNResoHqe smear_neutron3case2_hydrogen.C)
target_link_libraries(smearNResoHqe PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})
				    
add_executable(testit test.C)
target_link_libraries(testit PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})				 

add_executable(smearit2 smear2.C)
target_link_libraries(smearit2 PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})


add_executable(smearboth smearECAL_STT.C)
target_link_libraries(smearboth PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(smearLar smearECAL_STT_meniscus.C)
target_link_libraries(smearLar PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(smearTopBot  smearECAL_STT_meniscus_topbot.C)
target_link_libraries(smearTopBot PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(smear3dst smearECAL_3dst.C)
target_link_libraries(smear3dst PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(smear3dst2 smearECAL_3dst2.C)
target_link_libraries(smear3dst2 PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(smearHad smearECAL_STT_meniscus_hadreso.C)
target_link_libraries(smearHad PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(smearLarnoFV smearECAL_STT_meniscus_noFV.C)
target_link_libraries(smearLarnoFV PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(smearEcal_check smearECAL_STT_meniscus_check.C)
target_link_libraries(smearEcal_check PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})

add_executable(smearMuid smearECAL_STT_meniscus_check_muidascii.C)
target_link_libraries(smearMuid PUBLIC EDepSim::edepsim_io ${ROOT_LIBRARIES})
