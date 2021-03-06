cmake_minimum_required (VERSION 2.6)
PROJECT(Consus)

#set compiler options

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -std=c++1y")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS} -g -Og")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -mavx")
#set("${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -mavx -fno-omit-frame-pointer")

set(BOOST_INCLUDEDIR "/usr/include/boost/")
set(BOOST_LIBRARYDIR "/usr/lib/x86_64-linux-gnu/")
FIND_PACKAGE( Boost COMPONENTS iostreams filesystem system regex REQUIRED)
FIND_PACKAGE( OpenMP REQUIRED )

if(OPENMP_FOUND)
  message("OPENMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories( ${Boost_INCLUDE_DIR} )
find_package( BZip2 REQUIRED )
include_directories( ${BZip2_INCLUDE_DIR} )
include_directories(./extern/dlib/)
include_directories(./extern/Catch/single_include/)

add_subdirectory(exe)
