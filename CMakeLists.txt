cmake_minimum_required(VERSION 2.8)
project (JohansenTest)

add_executable (JohansenTest JohansenTest.cpp JohansenHelper.cpp)

if(WIN32)
	include_directories(gsl-1.16)
	target_link_libraries (JohansenTest ../gsl-1.16/dll/Win32/Debug/gsl ../gsl-1.16/dll/Win32/Debug/cblas)
elseif(UNIX)
	target_link_libraries (JohansenTest gsl gslcblas)
endif(WIN32)
