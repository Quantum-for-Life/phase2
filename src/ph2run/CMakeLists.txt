message(STATUS "Configuring ph2run")

add_executable(ph2run
        main.c
)
target_link_libraries(ph2run
        phase2
        log
        QuEST m
)
target_include_directories(ph2run PRIVATE
        ${PHASE2_INCLUDE_DIRS}
        ${LOG_INCLUDE_DIRS}
)