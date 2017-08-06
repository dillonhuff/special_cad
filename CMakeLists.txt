cmake_minimum_required(VERSION 2.6)

project(libpoly C)

add_executable(roots main.c cad_tree.c polynomial_examples.c utils.c)
target_link_libraries(roots poly gmp)
