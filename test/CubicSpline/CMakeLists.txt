cmake_minimum_required(VERSION 3.15)
project(CubicSplineTest)

aux_source_directory(. SRC)
set(EXECUTABLE_OUTPUT_PATH ${BIN_PATH_}/test/${PROJECT_NAME})

include_directories(${INCLUDE_PATH_})

add_executable(${PROJECT_NAME} ${SRC})