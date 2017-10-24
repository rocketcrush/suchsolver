project(suchsolver)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)


find_package(OpenMP)
if (OPENMP_FOUND)
 set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
 set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


#Bring the headers
include_directories(include)

#Can manually add the sources using the set command as follows:
set(SOURCES mesh.cpp constants.cpp initialconditions.cpp boundaryconditions.cpp fluxes.cpp calculations.cpp postprocess.cpp)


add_executable(${PROJECT_NAME} ${SRC_LIST})
