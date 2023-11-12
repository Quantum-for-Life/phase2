cmake_minimum_required(VERSION 3.9.0)
project(phase2 VERSION 0.1.0 LANGUAGES C)

#set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
#set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

#option(BUILD_SHARED_LIBS OFF)
add_subdirectory(lib/QuEST)

add_subdirectory(src)

