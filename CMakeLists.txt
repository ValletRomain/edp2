
cmake_minimum_required(VERSION 3.6)

# Godunov
project(godunov)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c+11")

set(GODUNOV_FILES src/godunov.c)
add_executable(godunov ${GODUNOV_FILES})

target_link_libraries(godunov PRIVATE m)

# Godunov error
project(godunov_error)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c+11")

set(GODUNOV_ERROR_FILES src/godunov_error.c)
add_executable(godunov_error ${GODUNOV_ERROR_FILES})

target_link_libraries(godunov_error PRIVATE m)