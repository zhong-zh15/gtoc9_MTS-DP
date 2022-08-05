cmake_minimum_required(VERSION 2.8.2)
project(GTOC9_parallel)

set(CMAKE_CXX_STANDARD 14)


#set(CMAKE_C_COMPILER   "${CROSS_TOOLCHAIN_PREFIX}mpiicc")
#set(CMAKE_CXX_COMPILER "${CROSS_TOOLCHAIN_PREFIX}mpiicpc")

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_C_FLAGS_RELEASE "$ENV{CFLAGS} -ipo -O3 -no-prec-div")
#set(CMAKE_C_FLAGS_RELEASE "$ENV{CFLAGS}")
#SET(CMAKE_BUILD_TYPE "Debug")  
#SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")  
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -std=c++14 ")

set(PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}/build)

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include/orbit_base
)

file(GLOB SRC_LIST
        "*.cpp" "*.c" "*.h" "*.hpp"
        "src/*.cpp" "src/*.c" "src/*.h"
        "include/orbit_base/*.cpp" "include/orbit_base/*.c" "include/orbit_base/*.h"
        )


set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_executable(GTOC9_parallel ${SRC_LIST})

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${MY_LIB})
