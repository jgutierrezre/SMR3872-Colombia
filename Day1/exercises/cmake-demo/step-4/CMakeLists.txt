# CMake configuration demo

cmake_minimum_required(VERSION 3.10)
project(canvas-draw VERSION 0.4 LANGUAGES CXX)

option(BUILD_SHARED_LIBS "Build shared library" OFF)
if(BUILD_SHARED_LIBS)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

add_executable(canvas-draw main.cpp)
configure_file(canvas_config.h.in canvas_config.h)
add_library(canvas canvas.cpp)
target_include_directories(canvas-draw PUBLIC ${CMAKE_BINARY_DIR})
target_link_libraries(canvas-draw PRIVATE canvas)
