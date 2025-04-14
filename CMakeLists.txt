cmake_minimum_required(VERSION 3.15)
project(libcdgbs)

set(CMAKE_CXX_STANDARD 17)

# ---- External Libraries ----

# Eigen (header-only)
include_directories(external/eigen)

# Triangle (build subdir with known output location)
set(TRIANGLE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/Triangle/src)
set(TRIANGLE_BINARY_DIR ${CMAKE_BINARY_DIR}/triangle_build)
add_subdirectory(${TRIANGLE_SOURCE_DIR} ${TRIANGLE_BINARY_DIR})

# Avoid building tools with Qt and custom macros
set(OPENMESH_BUILD_UNIT_TESTS OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_DECIMATER OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_NORMAL_DEVIATION OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_PYTHON_BINDINGS OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_APPS OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_UTILS OFF CACHE BOOL "" FORCE)

# OpenMesh and libgeom
add_subdirectory(external/OpenMesh)
add_subdirectory(external/libgeom)

# ---- Your Library ----
include_directories(include)

file(GLOB_RECURSE LIBCDGBS_SRC "src/*.cpp")

add_library(libcdgbs ${LIBCDGBS_SRC})

# 🔧 Add include paths for Triangle, including generated headers
target_include_directories(libcdgbs
    PRIVATE
        ${TRIANGLE_SOURCE_DIR}/Triangle
        ${TRIANGLE_BINARY_DIR}/Triangle  # This is where triangle_export.h lives
)

# 🔗 Link libraries
target_link_libraries(libcdgbs
    triangle
    OpenMeshCore
    OpenMeshTools
    geom
)

# ---- Test Executable ----
add_executable(libcdgbs_test tests/main_test.cpp)

# 🔧 Also add include path for generated headers here
target_include_directories(libcdgbs_test
    PRIVATE
        ${TRIANGLE_BINARY_DIR}/Triangle
)

target_link_libraries(libcdgbs_test libcdgbs)
