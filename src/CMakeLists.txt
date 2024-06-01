cmake_minimum_required(VERSION 3.24)
project(src)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(src
        distance.h
        main.cpp
        pqkmeans.cpp
        pqkmeans.h
        segr.h)
