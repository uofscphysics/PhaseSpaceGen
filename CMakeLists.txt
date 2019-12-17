cmake_minimum_required(VERSION 3.5)
project(ElasticGen C CXX Fortran)
set(VERSION 1.0.0)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "-Ofast ${ROOT_CXX_FLAGS}")
set(CMAKE_C_FLAGS "-Ofast ${ROOT_CXX_FLAGS}")

option(PLOTS "make plots output" OFF)
if(PLOTS)
    add_definitions(-DPLOTS)
endif(PLOTS)


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})


set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake
   ${PROJECT_SOURCE_DIR}/cmake/modules)

include_directories(include)


add_library(physics
    src/physics.cxx
    src/radcorr.cxx
)


set(PROGRAMS 
    ElasticGen
)

foreach(prog ${PROGRAMS})
    get_filename_component(progname ${prog} NAME)
    string( REPLACE ".cxx" "" progname ${progname} )
    add_executable( ${progname} ${prog} )
    target_link_libraries(${progname} PUBLIC physics ${ROOT_LIBRARIES})
    install(TARGETS ${progname} EXPORT ${PROJECT_NAME}Targets RUNTIME DESTINATION bin)
endforeach( prog ${PROGRAMS} )