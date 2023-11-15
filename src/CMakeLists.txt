find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

set(CMAKE_C_FLAGS "-Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror")

add_library(phase2 SHARED
        circ.c
        linen.c
        rayon.c
        sdat.c
        log/log.c
)
target_link_libraries(phase2
        QuEST m
        mpi
        ${HDF5_LIBRARIES}
)

add_executable(ph2run
        ph2run.c
)
target_link_libraries(ph2run
        phase2
        QuEST m
        mpi
)