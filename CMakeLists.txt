cmake_minimum_required(VERSION 3.9)

project(metropolis)

set(CMAKE_PREFIX_PATH ./argumentum ${CMAKE_PREFIX_PATH})

find_package(Qt5 COMPONENTS Core)
find_package(OpenMP REQUIRED)
find_package( Argumentum CONFIG REQUIRED )

include_directories("../partsEngine")
link_directories("../partsEngine")

add_executable(metropolis main.cpp)
target_link_libraries(metropolis partsEngine Qt5::Core gmp gmpxx OpenMP::OpenMP_CXX)