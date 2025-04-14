cmake_minimum_required(VERSION 3.15)
project(libcdgbs)

set(CMAKE_CXX_STANDARD 17)

# Include public headers
include_directories(include)

# Add Eigen (header-only)
include_directories(external/eigen)

# Add Triangle
add_subdirectory(external/Triangle)

# Add OpenMesh
add_subdirectory(external/OpenMesh)

# Library sources
file(GLOB_RECURSE LIBCDGBS_SRC "src/*.cpp")

# Create the library
add_library(libcdgbs ${LIBCDGBS_SRC})
target_link_libraries(libcdgbs Triangle OpenMeshCore OpenMeshTools)

# Test executable
add_executable(libcdgbs_test tests/main_test.cpp)
target_link_libraries(libcdgbs_test libcdgbs)
