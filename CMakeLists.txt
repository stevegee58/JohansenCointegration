cmake_minimum_required(VERSION 2.8)
project (JohansenTest)
add_executable (JohansenTest JohansenTest.cpp JohansenHelper.cpp)
target_link_libraries (JohansenTest gsl gslcblas)
