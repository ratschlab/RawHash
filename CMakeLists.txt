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

# Add SLOW5 Submodules
# Enable multi-threading in the slow5lib submodule
set(SLOW5_ENABLE_MT ON CACHE BOOL "Enable multi-threading in slow5lib" FORCE)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/slow5lib EXCLUDE_FROM_ALL)

# This tells CMake to look into the src folder for the actual source files
add_subdirectory(src)