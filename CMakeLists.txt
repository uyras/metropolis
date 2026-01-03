cmake_minimum_required(VERSION 3.9)

project(metropolis VERSION 0.1.0)

set (CMAKE_CXX_STANDARD 17)

set( ARGUMENTUM_BUILD_STATIC_LIBS ON INTERNAL )
add_subdirectory(argumentum-src)
add_subdirectory(inicpp)


find_package(OpenMP REQUIRED)

include_directories("argumentum-src/include" "inicpp/include")

set(METROPOLIS_SRC
	main.cpp
	ConfigManager.cpp
	CalculationParameter.cpp
	CorrelationCore.cpp 
#	CorrelationPointCore.cpp 
#	MagnetisationCore.cpp
#	MagnetisationLengthCore.cpp
	misc.cpp
	PtBalancerFactory.cpp
	PtBalancerManual.cpp
	Worker.cpp
	MagneticSystem.cpp
)

file(STRINGS examples/example.ini example_string_a)
list(JOIN example_string_a  "\\n"  example_string)

configure_file(defines.h.in defines.h)

add_executable(metropolis ${METROPOLIS_SRC})
target_link_libraries(metropolis gmp gmpxx OpenMP::OpenMP_CXX argumentum inicpp)
target_include_directories(metropolis PUBLIC "${PROJECT_BINARY_DIR}")

add_executable(distanceAnalyser distanceAnalyser.cpp)
target_link_libraries(distanceAnalyser argumentum)

add_executable(groundStater groundStater.cpp)
target_link_libraries(groundStater argumentum)