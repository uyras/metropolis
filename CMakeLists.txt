cmake_minimum_required(VERSION 3.9)

project(metropolis)

set (CMAKE_CXX_STANDARD 17)

set(CMAKE_PREFIX_PATH ./argumentum ${CMAKE_PREFIX_PATH})

find_package(OpenMP REQUIRED)
find_package(Argumentum CONFIG REQUIRED )

include_directories("../partsEngine")
link_directories("../partsEngine")

add_executable(metropolis main.cpp CorrelationCore.cpp)
target_link_libraries(metropolis partsEngine gmp gmpxx OpenMP::OpenMP_CXX Argumentum::argumentum)

add_executable(distanceAnalyser distanceAnalyser.cpp)
target_link_libraries(distanceAnalyser partsEngine Argumentum::argumentum)