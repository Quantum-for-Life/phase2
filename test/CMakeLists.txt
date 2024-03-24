set(TEST_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/test)
include_directories(${TEST_INCLUDE_DIRS})

macro(more_test name sources)
    string(CONCAT test_name "test_" ${name})
    string(CONCAT test_source ${test_name} ".c")
    add_executable(${test_name} ${test_source})
    target_sources(${test_name} PUBLIC ${${sources}})
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
endmacro()

more_test(circ PH2_COMMON_SOURCES)
more_test(circ_env PH2_COMMON_SOURCES)
more_test(healthchk PH2_COMMON_SOURCES)
more_test(caserand PH2_COMMON_SOURCES)

set(TEST_DATA_SRC
        ${PROJECT_SOURCE_DIR}/test/data
)
set(TEST_DATA_SOURCES
        ${DATA2_SOURCES}
        ${TEST_DATA_SRC}/test_data_open.c
        ${TEST_DATA_SRC}/test_data_multidet.c
        ${TEST_DATA_SRC}/test_data_hamil.c
        ${TEST_DATA_SRC}/test_data_trotter_steps.c
)
add_compile_definitions(TEST_DATA_SRC="${TEST_DATA_SRC}")
more_test(data TEST_DATA_SOURCES)
