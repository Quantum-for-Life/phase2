add_executable(ph2run
        main.c
        datio.c
        log.c
        opt.c
        run.c
)

target_link_libraries(ph2run
        circ
        data
        QuEST m
)

target_include_directories(ph2run PRIVATE
        ${CIRC_INCLUDE_DIRS}
        ${DATA_INCLUDE_DIRS}
)