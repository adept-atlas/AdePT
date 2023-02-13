# SPDX-FileCopyrightText: 2022 CERN
# SPDX-License-Identifier: Apache-2.0

find_package(HepMC3 QUIET)

if(HepMC3_FOUND)
  message(STATUS "HepMC3 found ${HEPMC3_INCLUDE_DIR}") 
  add_definitions(-DHEPMC3_FOUND)
endif()


if(NOT TARGET G4HepEm::g4HepEm)
  message(STATUS "Disabling exampleTileTB (needs G4HepEm)")
  return()
endif()

if(Geant4_FOUND)
  if(NOT Geant4_gdml_FOUND)
    message(STATUS "Disabling exampleTileTB (needs Geant4 with GDML support)")
    return()
  endif()
  if(NOT Geant4_flush_FOUND)
    message(STATUS "Disabling exampleTileTB (needs G4VFastSimulationModel::Flush interface available starting with geant4-11-00-ref-01)")
    return()
  endif()
else()
  message(STATUS "Disabling exampleTileTB (needs Geant4)")
  return()
endif()

set(sources_g4
  src/ActionInitialisation.cc
  src/DetectorConstruction.cc
  src/DetectorMessenger.cc
  src/EMShowerMessenger.cc
  src/EMShowerModel.cc
  src/EventAction.cc
  src/EventActionMessenger.cc
  src/SimpleHit.cc
  src/PrimaryGeneratorAction.cc
  src/PrimaryGeneratorMessenger.cc
  src/RunAction.cc
  src/SensitiveDetector.cc
  src/TileSD.cc
  src/PhysListHepEm.cc
  src/TrackingAction.cc
  src/FTFP_BERT_HepEm.cc
)

if(HepMC3_FOUND)
  set(sources_hepmc3
      src/HepMC3G4AsciiReader.cc
      src/HepMC3G4AsciiReaderMessenger.cc
      src/HepMC3G4Interface.cc
     )
endif()

# exampleTileTB
add_executable(exampleTileTB exampleTileTB.cpp AdeptIntegration.cpp AdeptIntegration.cu BasicScoring.cu TileScoring.cu
               ${sources_g4} ${sources_hepmc3})

target_include_directories(exampleTileTB PRIVATE ${PROJECT_SOURCE_DIR}/examples/ExampleTileTB/include 
                          ${PROJECT_SOURCE_DIR}/examples/ExampleTileTB ${HEPMC3_INCLUDE_DIR})

target_link_libraries(exampleTileTB
  PRIVATE
    AdePT
    CopCore::CopCore
    VecGeom::vecgeom
    VecGeom::vecgeomcuda_static
    VecGeom::vgdml
    ${Geant4_LIBRARIES}
    G4HepEm::g4HepEm
    G4HepEm::g4HepEmData
    G4HepEm::g4HepEmInit
    G4HepEm::g4HepEmRun
    CUDA::cudart
    ${HEPMC3_LIBRARIES} ${HEPMC3_FIO_LIBRARIES})

set_target_properties(exampleTileTB PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# Install macros and geometry file

set(TESTEM3_GDML ${CMAKE_BINARY_DIR}/testEm3.gdml)
set(TILETB_2B1EB_GDML ${CMAKE_BINARY_DIR}/TileTB_2B1EB.gdml)
configure_file("macros/testEm3.gdml" "${CMAKE_BINARY_DIR}/testEm3.gdml")
configure_file("macros/TileTB_2B1EB.gdml" "${CMAKE_BINARY_DIR}/TileTB_2B1EB.gdml")
configure_file("macros/testEm3.mac.in" "${CMAKE_BINARY_DIR}/testEm3.mac")
configure_file("macros/testEm3G4.mac.in" "${CMAKE_BINARY_DIR}/testEm3G4.mac")
configure_file("macros/testTile.mac.in" "${CMAKE_BINARY_DIR}/testTile.mac")
configure_file("macros/testTileG4.mac.in" "${CMAKE_BINARY_DIR}/testTileG4.mac")

# Tests
add_test(NAME exampleTileTB-adept-em3
  COMMAND $<TARGET_FILE:exampleTileTB> -m ${CMAKE_BINARY_DIR}/testEm3.mac)
add_test(NAME exampleTileTB-g4-em3
  COMMAND $<TARGET_FILE:exampleTileTB> -m ${CMAKE_BINARY_DIR}/testEm3G4.mac)
