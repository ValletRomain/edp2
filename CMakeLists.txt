
cmake_minimum_required(VERSION 3.6)

# Godunov
project(solve)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c+11")

set(SOLVE_FILES src/solve.c)
add_executable(solve ${SOLVE_FILES})

target_link_libraries(solve PRIVATE m)

# Godunov error
project(compute_error)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c+11")

set(ERROR_FILES src/compute_error.c)
add_executable(compute_error ${ERROR_FILES})

target_link_libraries(compute_error PRIVATE m)