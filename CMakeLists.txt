cmake_minimum_required(VERSION 2.6)

project(libpoly C)

add_executable(roots main.c)
target_link_libraries(roots poly)
