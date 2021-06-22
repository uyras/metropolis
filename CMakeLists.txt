cmake_minimum_required(VERSION 3.9)

project(metropolis)

set (CMAKE_CXX_STANDARD 17)

add_subdirectory(partsEngine)
add_subdirectory(argumentum-src)

add_custom_target(fix ALL
	COMMENT "fixing problem with argparse"

	COMMAND ${CMAKE_COMMAND}
	   -E copy
	   "${CMAKE_SOURCE_DIR}/argumentum-src/util/faketarget.cpp"
	   "${CMAKE_SOURCE_DIR}/util/faketarget.cpp"
)

find_package(OpenMP REQUIRED)

include_directories("partsEngine" "argumentum-src/include")

add_executable(metropolis main.cpp CorrelationCore.cpp)
target_link_libraries(metropolis partsEngine gmp gmpxx OpenMP::OpenMP_CXX argumentum)

add_executable(distanceAnalyser distanceAnalyser.cpp)
target_link_libraries(distanceAnalyser partsEngine argumentum)