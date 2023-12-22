macro(more_test name sources)
    string(CONCAT test_name "test_" ${name})
    string(CONCAT test_source ${test_name} ".c")
    add_executable(${test_name} ${test_source})
    target_sources(${test_name} PUBLIC ${${sources}})
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
endmacro()

more_test(circ CIRC_SOURCES)
more_test(circ_env CIRC_SOURCES)
more_test(healthchk CIRC_SOURCES)
more_test(rayon PH2_COMMON_SOURCES)

set(TEST_DATA_SOURCES
        ${DATA_SOURCES}
        data/test_data_open.c
)
more_test(data TEST_DATA_SOURCES)