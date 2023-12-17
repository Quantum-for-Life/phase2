add_executable(ph2run
        main.c
        datio.c
        log.c
        opt.c
        run.c
)

target_link_libraries(ph2run
        data
        QuEST m
        circ
)

target_include_directories(ph2run PRIVATE
        ${DATA_INCLUDE_DIRS}
        ${CIRC_INCLUDE_DIRS}
)

if (DISTRIBUTED)
    target_link_libraries(ph2run
            ${MPI_C_LIBRARIES}
    )
    target_include_directories(ph2run PRIVATE
            ${MPI_C_INCLUDE_DIRS}
    )
endif ()