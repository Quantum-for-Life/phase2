set(CMAKE_C_FLAGS "-Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror")

link_libraries(${HDF5_LIBRARIES})

add_library(phase2 SHARED
        circ.c
        linen.c
        rayon.c
        sdat.c
        log/log.c
)
target_link_libraries(phase2
        QuEST m
)

add_executable(ph2run
        ph2run.c
)
target_link_libraries(ph2run
        phase2
        QuEST m
)