cmake_minimum_required(VERSION 2.8.9)

project(test)
set(SOURCES test.cpp nIMParton.cpp)

add_library(nIMParton SHARED nIMParton.cpp)
add_executable(test ${SOURCES})
