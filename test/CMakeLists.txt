set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")

set(TEST_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})


message(STATUS This: ${PHASE2_INCLUDE_DIRS})

link_libraries(phase2)
include_directories(${PHASE2_INCLUDE_DIRS})

add_executable(test_circ circ/main.c)
target_include_directories(test_circ PRIVATE ${TEST_INCLUDE_DIRS})

add_test(NAME test_circ COMMAND $<TARGET_FILE:test_circ>)