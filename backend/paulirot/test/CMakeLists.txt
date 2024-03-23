set(TEST_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/test)
include_directories(${TEST_INCLUDE_DIR})

macro(more_test name sources)
    string(CONCAT test_name "test_" ${name})
    string(CONCAT test_source ${test_name} ".c")
    add_executable(${test_name} ${test_source})
    target_sources(${test_name} PUBLIC ${${sources}})
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
endmacro()

more_test(qreg QREG_SOURCES)
more_test(paulis QREG_SOURCES)