project(linear_algebra)

cmake_minimum_required(VERSION 3.8)

include_directories("./iml++" "./")

set(CMAKE_CXX_FLAGS "-std=c++14 -g -fopenmp -O3 -Wall -Wextra")

add_executable(linear_algebra DenseVector.hpp 
                              DenseVector.cpp
                              DenseMatrix.hpp
                              DenseMatrix.cpp
                              LinearSolverInterface.hpp
                              LinearSolverInterface.cpp
                              main.cpp)
