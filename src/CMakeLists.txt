find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

set(CMAKE_C_FLAGS "-Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror")

add_library(circ SHARED
        circ.c
        log/log.c
        simulh5.c
        simulh5.h
)
target_link_libraries(circ
        QuEST m
)


add_library(linen SHARED
        linen.c
        log/log.c
)
target_link_libraries(linen
        QuEST m
)


add_library(rayon SHARED
        rayon.c
        log/log.c
)
target_link_libraries(rayon
        QuEST m
        ${HDF5_LIBRARIES}
)


add_executable(phase2
        phase2.c
        simulh5.c
        log/log.c
)
target_link_libraries(phase2
        circ
        linen
        rayon
        ${HDF5_LIBRARIES}
)