cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

set(NAME_MAIN_PROGRAM runMain)

project(${NAME_MAIN_PROGRAM}) 

add_compile_options(-O2 -march=native -Wall)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/inc)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/inc/*.hh)

add_executable(${NAME_MAIN_PROGRAM} main.cc ${sources} ${headers})
target_link_libraries(${NAME_MAIN_PROGRAM} ${Geant4_LIBRARIES})

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/..)
install(TARGETS ${NAME_MAIN_PROGRAM} DESTINATION bin)

