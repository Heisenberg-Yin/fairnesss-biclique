cmake_minimum_required(VERSION 3.6)
project(MBEA)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

set(SOURCE_FILES main.cpp utilities/Defines.h utilities/Timer.h utilities/Utility.h src/BipartiteGraph.cpp src/BipartiteGraph.h src/CommandLineParser.cpp src/CommandLineParser.h)
add_executable(MBEA ${SOURCE_FILES})

include_directories(src)


