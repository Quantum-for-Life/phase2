find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

link_libraries(QuEST m)

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
