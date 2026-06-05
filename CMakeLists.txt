cmake_minimum_required(VERSION 3.14)
project(RawHash2 VERSION 1.0 LANGUAGES C CXX)

# Force executables and libraries to build in the project root
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

option(NOPOD5 "Disable POD5 support" ON)
find_package(ZLIB REQUIRED)

# We do not build hdf5 and pod5 dependencies
add_definitions(-DNPOD5RH=1)
add_definitions(-DNHDF5RH=1)

set(PYBIND11_FINDPYTHON ON)

# Add SLOW5 Submodules
# Enable multi-threading and build as static so _core.so has no libslow5.so runtime dependency
# CMAKE_POSITION_INDEPENDENT_CODE=ON ensures slow5 objects are compiled with -fPIC,
# which is required when linking a static lib into a shared object (.so).
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(SLOW5_ENABLE_MT ON CACHE BOOL "Enable multi-threading in slow5lib" FORCE)
set(SLOW5_LINK_STATIC ON CACHE BOOL "Build slow5 as static lib" FORCE)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/slow5lib EXCLUDE_FROM_ALL)
add_subdirectory(extern/pybind11)
include_directories(extern/parlaylib/include)

# This tells CMake to look into the src folder for the actual source files
add_subdirectory(src)