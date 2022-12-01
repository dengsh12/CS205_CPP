cmake_minimum_required(VERSION 3.20)
project(CS205_CppProject)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(CS205_CppProject
        Exceptions.hpp
        main.cpp
        Structure.hpp)
