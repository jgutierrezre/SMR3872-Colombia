# CMake configuration demo

# version 3.10 is bundled with ubuntu18.04
cmake_minimum_required(VERSION 3.10)

# set the project name, version, and primary programming language (C++)
project(canvas-draw VERSION 0.1 LANGUAGES CXX)

# add the executable
add_executable(canvas-draw main.cpp)
