set(TEST_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/test)
include_directories(${TEST_INCLUDE_DIRS})

macro(more_test name sources)
    string(CONCAT test_name "test-" ${name})
    string(CONCAT test_source ${test_name} ".c")
    add_executable(${test_name} ${test_source})
    target_sources(${test_name} PUBLIC ${${sources}})
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
endmacro()

more_test(caserand PHASE2_SOURCES)

set(TEST_DATA_SRC
        ${PROJECT_SOURCE_DIR}/test/data
)
set(TEST_DATA_SOURCES
        ${DATA_SOURCES}
        ${PHASE2_SOURCES}
        ${TEST_DATA_SRC}/test-data_open.c
        ${TEST_DATA_SRC}/test-data_multidet.c
        ${TEST_DATA_SRC}/test-data_hamil.c
        ${TEST_DATA_SRC}/test-data_trott_steps.c
)
add_compile_definitions(TEST_DATA_SRC="${TEST_DATA_SRC}")
more_test(data TEST_DATA_SOURCES)
