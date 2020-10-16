
cmake_minimum_required(VERSION 2.6)

project(Calcul_Sci)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

add_executable(
	godunov
    src/godunov.c
)

add_executable(
	godunov_error
    src/godunov_error.c
)
