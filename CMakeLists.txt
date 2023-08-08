## ---------------------------------------------------------
## nexus | CMakeLists.txt
##
## The NEXT Collaboration
## --------------------------------------------------------

cmake_minimum_required(VERSION 3.12)

project(nexus 
        DESCRIPTION "Geant4 simulation framework of the NEXT experiment" 
        LANGUAGES C CXX)

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Find dependencies
find_package(Geant4 REQUIRED ui_all vis_all)
find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED)

option (With_Opticks OFF)
if(With_Opticks)

  message(STATUS "WITH_OPTICKS is set ON")
  list(APPEND CMAKE_MODULE_PATH "$ENV{OPTICKS_HOME}/cmake/Modules")
  #list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")
  include(OpticksBuildOptions)
  message(STATUS "---- OPTICKS ----")
  message(STATUS "CMAKE_MODULE_PATH --> ${CMAKE_MODULE_PATH}")
  message(STATUS "OPTICKS HOME --> $ENV{OPTICKS_HOME}")
  message(STATUS "---- OPTICKS END ----")
  find_dependency(PLog) # For Some Reason it was not finding PLog , doing this compiles
  ## These are needed for opticks
  find_package(G4CX REQUIRED CONFIG)
  find_package(U4 REQUIRED CONFIG)
  find_package(ExtG4 REQUIRED CONFIG)
  find_package(CSG_GGeo REQUIRED CONFIG)

endif ()


### Garfield Implementation
option(With_Garfield OFF)
if(With_GarField)
  ## ROOT
  find_package(ROOT QUIET REQUIRED COMPONENTS RIO Net)
  if(NOT ROOT_FOUND)
    message(STATUS "ROOT package not found.")
    return()
  endif()

  #include(${ROOT_USE_FILE})

  ## Garfield
  if(NOT TARGET Garfield::Garfield)
    find_package(Garfield REQUIRED)
  endif()

  if(NOT Garfield_FOUND)
    message(STATUS "Garfield package not found.")
    return()
  endif()

  ## NEST
  find_package(NEST REQUIRED)
  if(NOT NEST_FOUND)
    message(STATUS "NEST package not found.")
    return()
  endif()
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
endif()

# Define list with names of source folders
set(SOURCE_DIRS actions base generators geometries materials
                persistency physics physics_lists sensdet utils)

# Create temporary folder to contain headers
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})


add_library(lib SHARED)
set_target_properties(lib PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
foreach(DIR ${SOURCE_DIRS})
  # Create lists of *.h and *.cc files
  file(GLOB HDRS ${CMAKE_SOURCE_DIR}/source/${DIR}/*.h)
  file(GLOB SRCS ${CMAKE_SOURCE_DIR}/source/${DIR}/*.cc)
  # Copy header files to the 'include' binary folder
  file(COPY ${HDRS} DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
  ### Include Garfield Related Libraries
  if(With_GarField)
    file(GLOB GRHDRS ${CMAKE_SOURCE_DIR}/source/${DIR}/Garfield/*.h)
    file(GLOB GRSRCS ${CMAKE_SOURCE_DIR}/source/${DIR}/Garfield/*.cc)
    file(COPY ${GRHDRS} DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
  elseif (With_Opticks)
    file(GLOB OpRHDRS ${CMAKE_SOURCE_DIR}/source/${DIR}/Opticks/*.h)
    file(GLOB OpSRCS ${CMAKE_SOURCE_DIR}/source/${DIR}/Opticks/*.cc)
    file(COPY ${OpRHDRS} DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
  elseif (With_Opticks AND With_GarField)
    file(GLOB GRHDRS ${CMAKE_SOURCE_DIR}/source/${DIR}/Garfield/*.h)
    file(GLOB GRSRCS ${CMAKE_SOURCE_DIR}/source/${DIR}/Garfield/*.cc)
    file(COPY ${GRHDRS} DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
    file(GLOB OpRHDRS ${CMAKE_SOURCE_DIR}/source/${DIR}/Opticks/*.h)
    file(GLOB OpSRCS ${CMAKE_SOURCE_DIR}/source/${DIR}/Opticks/*.cc)
    file(COPY ${OpRHDRS} DESTINATION ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME})
  endif()


  target_include_directories(lib PUBLIC
          $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}>
          $<INSTALL_INTERFACE:include/${PROJECT_NAME}>)

  if(With_GarField)
    target_sources(lib PRIVATE ${SRCS} ${GRSRCS})
  elseif (With_Opticks)
    target_sources(lib PRIVATE ${SRCS} ${OpSRCS})
  elseif (With_Opticks AND With_GarField)
    target_sources(lib PRIVATE ${SRCS} ${OpSRCS} ${GRSRCS})
  else ()
    target_sources(lib PRIVATE ${SRCS})
  endif()
endforeach()

if(With_GarField)
  target_include_directories(lib PRIVATE ${Geant4_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${NEST_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
  target_link_libraries(lib PUBLIC ${Geant4_LIBRARIES} PRIVATE ${GSL_LIBRARIES} ${HDF5_LIBRARIES} Garfield::Garfield ROOT::RIO ROOT::Net NEST::NESTCore NEST::NESTG4)
elseif (With_GarField AND With_Opticks)
  target_include_directories(lib PRIVATE ${Geant4_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${NEST_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
  target_link_libraries(lib PUBLIC ${Geant4_LIBRARIES} PRIVATE ${GSL_LIBRARIES} ${HDF5_LIBRARIES} Garfield::Garfield ROOT::RIO ROOT::Net NEST::NESTCore NEST::NESTG4 Opticks::G4CX Opticks::U4 Opticks::ExtG4 Opticks::CSG_GGeo)
elseif (With_Opticks)
  target_include_directories(lib PRIVATE ${Geant4_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${NEST_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
  target_link_libraries(lib PUBLIC ${Geant4_LIBRARIES} PRIVATE ${GSL_LIBRARIES} ${HDF5_LIBRARIES} Garfield::Garfield ROOT::RIO ROOT::Net NEST::NESTCore NEST::NESTG4 Opticks::G4CX Opticks::U4 Opticks::ExtG4 Opticks::CSG_GGeo)
else()
  target_include_directories(lib PRIVATE ${Geant4_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
  target_link_libraries(lib PUBLIC
          ${Geant4_LIBRARIES} PRIVATE
          ${GSL_LIBRARIES} ${HDF5_LIBRARIES} )
endif()


add_executable(exe)
set_target_properties(exe PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_sources(exe PRIVATE ${CMAKE_SOURCE_DIR}/source/nexus.cc)
target_link_libraries(exe PRIVATE lib)

add_executable(test)
set_target_properties(test PROPERTIES OUTPUT_NAME ${PROJECT_NAME}-test)

file(GLOB TESTS ${CMAKE_SOURCE_DIR}/source/tests/*/*.cc)
target_sources(test PRIVATE ${TESTS} ${CMAKE_SOURCE_DIR}/source/nexus-test.cc)
target_include_directories(test PRIVATE ${CMAKE_SOURCE_DIR}/source/tests)
target_link_libraries(test PRIVATE lib)


install(TARGETS lib exe test
        RUNTIME DESTINATION bin  
        LIBRARY DESTINATION lib)

install(DIRECTORY ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME} DESTINATION include
        FILES_MATCHING PATTERN "*.h")
