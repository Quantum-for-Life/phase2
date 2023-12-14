add_executable(ph2run
        main.c
        datio.c
        log.c
        opt.c
        run.c
)

target_link_libraries(ph2run
        phase2
        QuEST m
)

target_include_directories(ph2run PRIVATE
        ${PHASE2_INCLUDE_DIRS}
)