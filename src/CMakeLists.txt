find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})


add_library(circ SHARED
        circ.c
        logger.c
)
target_link_libraries(circ
        QuEST m
)


add_library(linen SHARED
        linen.c
        logger.c
)
target_link_libraries(linen
        QuEST m
)


add_library(rayon SHARED
        rayon.c
        logger.c
)
target_link_libraries(rayon
        QuEST m
        ${HDF5_LIBRARIES}
)


add_executable(phase2
        simul.c
        logger.c
)
target_link_libraries(phase2
        circ
        linen
        rayon
)