set(DATA_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(DATA_INCLUDE_DIRS ${DATA_INCLUDE_DIRS} PARENT_SCOPE)

add_library(data SHARED
        data.c
)

target_link_libraries(data)

target_include_directories(data PRIVATE
        ${DATA_INCLUDE_DIRS}
)
