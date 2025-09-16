cmake_minimum_required(VERSION 4.1)

set(CMAKE_CXX_STANDARD 14)

project(MVR_K1)

include_directories(${CMAKE_SOURCE_DIR}/src)

set(INCLUDES src/Source.cpp src/TestTask.cpp src/MainTask.cpp)

add_executable(MVR_K1 ${INCLUDES})