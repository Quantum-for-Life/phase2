# link QuEST
if (WIN32)
    link_libraries(QuEST)
else ()
    link_libraries(QuEST m)
endif ()

add_library(circ SHARED
        circ.c
        logger.c
)

add_library(linen SHARED
        linen.c
        logger.c
)

add_library(rayon SHARED
        rayon.c
        logger.c
)
target_link_libraries(rayon
        hdf5_serial
)

add_executable(qpe-simul
        simul.c
        logger.c
)
target_link_libraries(qpe-simul
        circ
        linen
        rayon
)