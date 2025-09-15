cmake_minimum_required(VERSION 4.1)

set(CMAKE_CXX_STANDARD 14)

project(MVR_K1)


set(INCLUDES src/Source.cpp src/TestTask.cpp)

add_executable(MVR_K1 ${INCLUDES})