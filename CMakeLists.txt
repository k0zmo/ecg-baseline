cmake_minimum_required(VERSION 2.8.11)

project(ecg-baseline)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Qt5Core REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCE_FILES butter.cpp
                 ecgdata.cpp
                 main.cpp
                 sgolay.cpp
                 wavelet.cpp)

set(HEADER_FILES butter.h
                 sgolay.h
                 wavelet.h)

add_executable(ecg-baseline ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(ecg-baseline Qt5::Core)
