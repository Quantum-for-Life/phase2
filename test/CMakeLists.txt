set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")

message(STATUS This: ${PHASE2_INCLUDE_DIRS})

link_libraries(phase2)
include_directories(${PHASE2_INCLUDE_DIRS})

add_executable(test_circ circ/main.c)

add_test(NAME test_circ COMMAND $<TARGET_FILE:test_circ>)