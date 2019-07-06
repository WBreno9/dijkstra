cmake_minimum_required(VERSION 3.7.3)
project(dijkstra)

set(CMAKE_CXX_STANDARD 17)

set(SOURCES
	src/main.cpp
	)

add_executable(dijkstra ${SOURCES})
target_link_libraries(dijkstra PRIVATE m pthread)
