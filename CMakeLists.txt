# Set the minimum requirement for building this project
cmake_minimum_required(VERSION 3.22.1)

# Set your desired C++ standard
set(CMAKE_CXX_STANDARD 14)  
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set Project Name
project(OFDM)

# Append the -g flag to the compiler options
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Add the library we want to build
# add_library(ofdm SHARED constellation.cpp ofdm.cpp)

# Use a wildcard pattern to include all .cpp files in the src directory
file(GLOB SOURCES "src/*.cpp")

# Include directories
include_directories(/Applications/MATLAB_R2022b.app/extern/include)
include_directories(${CMAKE_SOURCE_DIR}/include/*)

# Add the executable that we want to build
add_executable(ofdm ${SOURCES})

# Specify that we need to link against our created library
target_link_libraries(ofdm PRIVATE  ${CMAKE_SOURCE_DIR}/lib/libfftw3.dylib)
target_link_libraries(ofdm PRIVATE  /Applications/MATLAB_R2022b.app/bin/maci64/libmat.dylib)
target_link_libraries(ofdm PRIVATE  /Applications/MATLAB_R2022b.app/bin/maci64/libmx.dylib)
# target_link_libraries(main PRIVATE ofdm)



