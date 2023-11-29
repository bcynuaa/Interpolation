cmake_minimum_required(VERSION 3.15)
project(Interpolation LANGUAGES C)

set(SRC_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INCLUDE_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(BIN_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(LIB_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(TEST_PATH_ ${CMAKE_CURRENT_SOURCE_DIR}/test)

add_subdirectory(${TEST_PATH_}/CubicSpline)
add_subdirectory(${TEST_PATH_}/BinaryInterpolation)