cmake_minimum_required(VERSION 3.12)
project(tiled_parallel_fw)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")


set(CMAKE_CXX_STANDARD 14)

add_executable(tiled_parallel_fw main.cpp)