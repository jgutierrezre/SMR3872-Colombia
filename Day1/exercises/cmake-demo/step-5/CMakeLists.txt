# CMake configuration demo

cmake_minimum_required(VERSION 3.10)
project(canvas-draw VERSION 0.5 LANGUAGES CXX)
configure_file(canvas_config.h.in canvas_config.h)

option(BUILD_SHARED_LIBS "Build shared library" OFF)
if(BUILD_SHARED_LIBS)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

find_package(JPEG QUIET)
option(USE_JPEG "Enable JPEG support" ${JPEG_FOUND})
add_executable(canvas-draw main.cpp)

add_library(canvas canvas.cpp)
if(USE_JPEG)
  find_package(JPEG REQUIRED)
  target_compile_definitions(canvas PRIVATE -DHAVE_JPEG_LIB)
  target_include_directories(canvas PRIVATE ${JPEG_INCLUDE_DIRS})
  target_link_libraries(canvas PRIVATE ${JPEG_LIBRARIES})
endif()

target_include_directories(canvas-draw PUBLIC ${CMAKE_BINARY_DIR})
target_link_libraries(canvas-draw PRIVATE canvas)
