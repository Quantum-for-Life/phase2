# link QuEST
if (WIN32)
    link_libraries(circ QuEST)
else ()
    link_libraries(circ QuEST m)
endif ()

add_executable(test_circ
        test_circ.c
        mock_circ.c
)

add_test(
        NAME test_circ
        COMMAND test_circ
)