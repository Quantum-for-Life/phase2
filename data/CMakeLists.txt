set(DATA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(DATA_INCLUDE_DIRS ${DATA_INCLUDE_DIRS} PARENT_SCOPE)


if (DISTRIBUTED)
    set(HDF5_PREFER_PARALLEL TRUE)
endif ()
find_package(HDF5 REQUIRED)
if (HDF5_IS_PARALLEL)
    message(STATUS "Parallel HDF5 found")
else ()
    message(STATUS "Serial HDF5 found")
endif ()


add_library(data STATIC
        data.c
)

target_link_libraries(data
        ${HDF5_LIBRARIES}
)

target_include_directories(data PRIVATE
        ${HDF5_INCLUDE_DIRS}
        ${DATA_INCLUDE_DIRS}
)

if (DISTRIBUTED)
    target_include_directories(data PRIVATE
            ${MPI_C_INCLUDE_DIRS}
    )
endif ()
