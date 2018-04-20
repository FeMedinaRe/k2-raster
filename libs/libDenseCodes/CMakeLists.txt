cmake_minimum_required(VERSION 2.8.4)
project(DenseCode)
set(PROJECT_VENDOR "Fernando Silva")
set(CMAKE_PROJECT_NAME "libDenseCode")
set(PROJECT_CONTACT "fernando.silva@udc.es")
set(CMAKE_INSTALL_PREFIX /home/fsilva/software)
set(CMAKE_COMPILER_IS_GNUCXX)


# Set output folder
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(./include/)

set(SOURCES
        src/DC/SCDCV2.cpp)

# Flags
set(CMAKE_CXX_FLAGS " -Wall -m64")
set(CMAKE_CXX_FLAGS_RELEASE " -std=c++0x -Wall -O9")
set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -O9 ")
#set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -g3 -DDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -pg ")

# set default compilation and show details for debug build
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif (NOT CMAKE_BUILD_TYPE)

# Create the library
add_library(denseCode STATIC ${SOURCES})
install(TARGETS denseCode DESTINATION lib)
install(DIRECTORY "includes" DESTINATION includes/denseCode)


##### TESTS ####

# Executables
add_executable(TestSCDC ./test/testSCDC.cpp)
add_executable(TestSCDCV2 ./test/testSCDCV2.cpp)

# Target link libraries to executables
target_link_libraries(TestSCDC LINK_PUBLIC denseCode)
target_link_libraries(TestSCDCV2 LINK_PUBLIC denseCode)