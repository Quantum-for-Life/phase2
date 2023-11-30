set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")

set(TEST_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

link_libraries(phase2)
include_directories(${TEST_INCLUDE_DIRS} ${PHASE2_INCLUDE_DIRS})

add_executable(test_circ test_circ.c)

add_test(NAME test_circ COMMAND $<TARGET_FILE:test_circ>)