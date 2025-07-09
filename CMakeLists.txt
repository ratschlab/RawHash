cmake_minimum_required(VERSION 3.15)
project(RawHash2Root)

add_subdirectory(src)

add_subdirectory(extern/pybind11)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")

set(PYBIND11_FINDPYTHON ON)
pybind11_add_module(_core src/pybindings.cpp)

#target_compile_options(_core PRIVATE -fpermissive -mavx)
target_link_libraries(_core PUBLIC rawhash2_lib pybind11::headers)
install(TARGETS _core DESTINATION pyrawhash)
