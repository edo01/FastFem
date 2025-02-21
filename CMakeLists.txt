cmake_minimum_required(VERSION 3.10)
project(FastFem)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Make OpenMP optional
option(USE_OPENMP "Enable OpenMP support" ON)

# Set compiler flags for different build types
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Allow narrowing conversions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")

# Enable optional fast math optimizations.
option(USE_FAST_MATH "Enable fast math optimizations." OFF)

if(USE_FAST_MATH)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -funsafe-math-optimizations")
    message(STATUS "Enabled fast math optimizations (-Ofast -funsafe-math-optimizations)")
endif()

# Add include directories
include_directories(include)

# Collect all source files
file(GLOB_RECURSE LIB_SOURCES "src/*.cpp")
file(GLOB_RECURSE TEST_SOURCES "test/*.cpp")

# Define the shared library
add_library(FastFem SHARED ${LIB_SOURCES})
target_include_directories(FastFem PUBLIC include)

if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(FastFem PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(FastFem PUBLIC HAVE_OPENMP)
    else()
        message(WARNING "OpenMP requested but not found. Building without OpenMP support.")
    endif()
endif()

# Option to enable/disable building tests
option(BUILD_TESTS "Build tests for FastFem" ON)

# Build tests only if BUILD_TESTS is ON
if(BUILD_TESTS)
    message(STATUS "Building tests for FastFem")
    foreach(TEST_FILE ${TEST_SOURCES})
        get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
        add_executable(${TEST_NAME} ${TEST_FILE})
        target_link_libraries(${TEST_NAME} FastFem)
        
        if(OpenMP_CXX_FOUND AND USE_OPENMP)
            target_link_libraries(${TEST_NAME} OpenMP::OpenMP_CXX)
        endif()
    endforeach()
else()
    message(STATUS "Skipping tests for FastFem")
endif()