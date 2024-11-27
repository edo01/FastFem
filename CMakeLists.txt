cmake_minimum_required(VERSION 3.10)

# Set the project name
project(FastFEM)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add include directories
include_directories(include)

# Find all source files in the src directory
file(GLOB SOURCES "src/**/*.cpp")

# Find and enable OpenMP
find_package(OpenMP)

# Add the executable
add_executable(FastFEM ${SOURCES})

# Link OpenMP libraries
if(OpenMP_CXX_FOUND)
    target_link_libraries(FastFEM PUBLIC OpenMP::OpenMP_CXX)
endif()
