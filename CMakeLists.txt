cmake_minimum_required(VERSION 3.9)

project(metropolis VERSION 0.1.1)


option( BUILD_TESTS      "Build tests"    OFF )

set (CMAKE_CXX_STANDARD 17)

set( ARGUMENTUM_BUILD_STATIC_LIBS ON INTERNAL )
add_subdirectory(argumentum-src)
add_subdirectory(inicpp)


find_package(OpenMP REQUIRED)


set(METROPOLIS_SRC
	main.cpp
	misc.cpp
	ConfigManager.cpp
	caluclationParameters/CalculationParameter.cpp
	caluclationParameters/CorrelationCore.cpp 
	caluclationParameters/CorrelationPointCore.cpp 
	caluclationParameters/MagnetisationCore.cpp
	caluclationParameters/MagnetisationLengthCore.cpp
	PtBalancerFactory.cpp
	PtBalancerManual.cpp
	Worker.cpp
	MagneticSystem.cpp
	dos2.cpp
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

target_include_directories(metropolis PUBLIC ${CMAKE_SOURCE_DIR}/argumentum-src/include PUBLIC ${CMAKE_SOURCE_DIR}/inicpp/include)
target_include_directories(metropolis PRIVATE ${CMAKE_SOURCE_DIR}/caluclationParameters)
include_directories(".")

if( BUILD_TESTS )
   enable_testing()
   add_subdirectory( tests )
endif()