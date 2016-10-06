MESSAGE("Hello, world")
IF(DEFINED APPLE)
    MESSAGE("APPLE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF(DEFINED APPLE)

set(CMAKE_OSX_DEPLOYMENT_TARGET "")
INCLUDE(SetCompileFlag.cmake)

cmake_minimum_required(VERSION 2.6)
project(directPsi CXX Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")


## include CERN ROOT
include(FindROOT.cmake)
MESSAGE(${ROOT_LIBRARIES})


## set fortran flags
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                 Fortran "-fno-automatic")



include_directories(~/local/include)


include_directories(src)
include_directories(${ROOT_INCLUDE_DIR})
include_directories(/usr/local/Cellar/boost/1.58.0/include/)

include_directories(./)

# adding sources for calcInts
file(GLOB_RECURSE calcInt_SOURCES "src/*.c++" "src/*.cpp" "src/*.for" )
list(REMOVE_ITEM calcInt_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/export_data.c++")
list(REMOVE_ITEM calcInt_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/double_export.c++")
list(REMOVE_ITEM calcInt_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/partonic.c++")
list(REMOVE_ITEM calcInt_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/hadronic.c++")
list(REMOVE_ITEM calcInt_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/check.c++")
list(REMOVE_ITEM calcInt_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/delta.c++")
list(REMOVE_ITEM calcInt_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/mtr0.c++")
add_executable (calcInts.exe "src/calcInts.c++" ${calcInt_SOURCES})
target_link_libraries(calcInts.exe ${ROOT_LIBRARIES})

#adding sources for export_data
add_executable (export_data.exe "src/export_data.c++" "src/utils.c++")
target_link_libraries(export_data.exe ${ROOT_LIBRARIES})

#adding sources for double_export
add_executable (double_export.exe "src/double_export.c++" "src/utils.c++")
target_link_libraries(double_export.exe ${ROOT_LIBRARIES})

#adding sources for partonic
add_executable (partonic.exe "src/partonic.c++" "src/utils.c++")
target_link_libraries(partonic.exe ${ROOT_LIBRARIES})

#adding sources for hadronic
## search for LHAPDF libs
EXEC_PROGRAM(lhapdf-config
    ARGS "--libs"
    OUTPUT_VARIABLE LHAPDF_LIBS
)
MESSAGE("LHAPDF_LIBS" ${LHAPDF_LIBS})

# adding sources for hadronic
#add_executable (hadronic.exe "src/hadronic.c++" "src/utils.c++")
#target_link_libraries(hadronic.exe ${ROOT_LIBRARIES})
#target_link_libraries(hadronic.exe ${LHAPDF_LIBS})

# adding sources for hadronic
file(GLOB_RECURSE hadronic_SOURCES "src/*.c++" "src/*.cpp" "src/*.for" )
list(REMOVE_ITEM hadronic_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/export_data.c++")
list(REMOVE_ITEM hadronic_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/double_export.c++")
list(REMOVE_ITEM hadronic_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/partonic.c++")
list(REMOVE_ITEM hadronic_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/delta.c++")
list(REMOVE_ITEM hadronic_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/check.c++")
list(REMOVE_ITEM hadronic_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/calcInts.c++")
list(REMOVE_ITEM hadronic_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/mtr.c++")
list(REMOVE_ITEM hadronic_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/mtr0.c++")
add_executable (hadronic.exe "src/hadronic.c++" ${hadronic_SOURCES})
target_link_libraries(hadronic.exe ${ROOT_LIBRARIES})
target_link_libraries(hadronic.exe ${LHAPDF_LIBS})


#adding sources for check
add_executable (check.exe "src/check.c++" "src/utils.c++")
target_link_libraries(check.exe ${ROOT_LIBRARIES})

# adding sources for delta
file(GLOB_RECURSE delta_SOURCES "src/*.c++" "src/*.cpp" "src/*.for" )
list(REMOVE_ITEM delta_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/export_data.c++")
list(REMOVE_ITEM delta_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/double_export.c++")
list(REMOVE_ITEM delta_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/partonic.c++")
list(REMOVE_ITEM delta_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/hadronic.c++")
list(REMOVE_ITEM delta_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/check.c++")
list(REMOVE_ITEM delta_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/calcInts.c++")
list(REMOVE_ITEM delta_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/mtr.c++")
add_executable (delta.exe "src/delta.c++" ${delta_SOURCES})
target_link_libraries(delta.exe ${ROOT_LIBRARIES})
