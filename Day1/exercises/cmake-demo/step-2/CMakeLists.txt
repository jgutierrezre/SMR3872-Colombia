# CMake configuration demo

cmake_minimum_required(VERSION 3.10)

project(canvas-draw VERSION 0.2 LANGUAGES CXX)

add_executable(canvas-draw main.cpp)
add_library(canvas canvas.cpp)
target_link_libraries(canvas-draw PRIVATE canvas)
