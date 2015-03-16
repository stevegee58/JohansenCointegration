cmake_minimum_required(VERSION 2.8)
project (JohansenTest)

add_executable (JohansenTest JohansenTest.cpp JohansenHelper.cpp)
add_executable (PairsFinder PairsFinder.cpp JohansenHelper.cpp)

if(WIN32)
	include_directories(gsl-1.16 c:/boost_1_53_0)
	target_link_libraries (JohansenTest ../gsl-1.16/dll/Win32/Debug/gsl ../gsl-1.16/dll/Win32/Debug/cblas)
	target_link_libraries (PairsFinder ../gsl-1.16/dll/Win32/Debug/gsl ../gsl-1.16/dll/Win32/Debug/cblas c:/boost_1_53_0/stage/lib/libboost_system-vc110-mt-gd-1_53.lib c:/boost_1_53_0/stage/lib/libboost_filesystem-vc110-mt-gd-1_53.lib)
elseif(UNIX)
	target_link_libraries (JohansenTest gsl gslcblas)
	target_link_libraries (PairsFinder gsl gslcblas boost_system boost_filesystem)
endif(WIN32)
