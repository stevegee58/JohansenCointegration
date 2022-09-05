cmake_minimum_required(VERSION 2.8)
project (JohansenTest)

add_executable (JohansenTest JohansenTest.cpp JohansenHelper.cpp JohansenHelper.h)

if(WIN32)
	include_directories(../gsl)
	target_link_libraries (JohansenTest ../../gsl/lib/x64/Release/gsl ../../gsl/lib/x64/Release/cblas)
elseif(UNIX)
	target_link_libraries (JohansenTest gsl gslcblas)
endif(WIN32)
