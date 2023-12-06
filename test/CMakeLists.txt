set(TEST_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})

link_libraries(phase2)
include_directories(${TEST_INCLUDE_DIRS} ${PHASE2_INCLUDE_DIRS})

add_executable(test_circ test_circ.c)
add_test(NAME test_circ COMMAND $<TARGET_FILE:test_circ>)

add_executable(test_healthcheck test_health.c)
add_test(NAME test_healthcheck COMMAND $<TARGET_FILE:test_healthcheck>)

add_executable(test_rayon test_rayon.c)
target_compile_definitions(test_rayon
        PUBLIC
        DATA_DIR="${PROJECT_SOURCE_DIR}/dat"
)
add_test(NAME test_rayon COMMAND $<TARGET_FILE:test_rayon>)

