cmake_minimum_required(VERSION 3.14) # FetchContent requires 3.11+
include(FetchContent)

# Project-wide Compiler Flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Handle DEBUG/PROFILE via standard CMake build types
# Usage: cmake -DCMAKE_BUILD_TYPE=Debug ..
set(COMMON_FLAGS "-Wall -march=native -pthread -DHAVE_KALLOC")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-g)
#    add_definitions(-g -fsanitize=address)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS} -O2")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_FLAGS} -O2 -Wc++-compat")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS} -O3")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_FLAGS} -O3 -Wc++-compat")
endif()

set(SLOW5_INC "${CMAKE_SOURCE_DIR}/extern/slow5lib/include")

# Architecture Specifics
if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm")
    add_definitions(-D_FILE_OFFSET_BITS=64 -fsigned-char)
    if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
        add_definitions(-mfpu=neon)
    endif()
endif()

# Define the Executable and Sources
set(RAW_HASH_SOURCES
        kthread.c kalloc.c bseq.cpp roptions.c sequence_until.c
        rutils.cpp rsig.cpp revent.cpp rsketch.c rindex.cpp
        lchain.cpp rseed.cpp rmap.cpp dtw.cpp hit.cpp
)

add_library(rawhash2_lib SHARED ${RAW_HASH_SOURCES})
set_target_properties(rawhash2_lib PROPERTIES OUTPUT_NAME "rawhash2")

# Linking Everything Together
target_include_directories(rawhash2_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${SLOW5_INC}
        ${POD5_INCLUDE_DIR}
)

# System Libraries
target_link_libraries(rawhash2_lib PRIVATE
        slow5
        ZLIB::ZLIB
        dl
        m
        pthread
)

add_executable(rawhash2 main.cpp)
target_link_libraries(rawhash2 PRIVATE rawhash2_lib)