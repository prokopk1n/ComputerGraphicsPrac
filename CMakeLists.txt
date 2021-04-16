cmake_minimum_required(VERSION 3.7)
project(witch_ball LANGUAGES CXX)

add_compile_options(-fopenmp -O3 -ffast-math -flto -march=native)

include_directories("lib/" "src/")
file(GLOB_RECURSE SRC RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp" "lib/*.h")

FIND_PACKAGE(OpenMP)
  if(OPENMP_FOUND)
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
 include_directories(SYSTEM ${OpenMP_INCLUDE_PATH})

add_executable(main ${SRC} ${OpenMP_CXX_LIBRARIES})
target_link_libraries(main)