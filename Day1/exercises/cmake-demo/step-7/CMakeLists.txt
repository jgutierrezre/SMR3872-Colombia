# CMake configuration demo

cmake_minimum_required(VERSION 3.10)
project(canvas-draw VERSION 0.7 LANGUAGES CXX)
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

if(ENABLE_TESTING)
  enable_testing()

  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)
  add_executable(test_color test_color.cpp)
  target_link_libraries(test_color PRIVATE GTest::GTestMain GTest::GTest)
  add_test(NAME Color COMMAND test_color)

  # does the application run?
  add_test(NAME Runs COMMAND canvas-draw)

  # does it create the file write.ppm?
  file(WRITE ${CMAKE_BINARY_DIR}/test_write.sh "rm -f *.ppm && ./canvas-draw && test -f white.ppm")
  add_test(NAME Writes COMMAND bash test_write.sh)
endif()
