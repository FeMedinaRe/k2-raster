cmake_minimum_required(VERSION 2.8)
project("Malloc Count" C)


include_directories(./include)

set (SOURCES
        include/malloc_count/malloc_count.h
        include/malloc_count/memprofile.h
        include/malloc_count/stack_count.h
        src/malloc_count.cpp
        src/stack_count.cpp
        #main.cpp
        )



# Flags
set(CMAKE_CXX_FLAGS " -Wall -m64")
set(CMAKE_CXX_FLAGS_RELEASE " -std=c++0x -Wall -O9 -DNDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -Wall -O9 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -Wall -g3 -DDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG " -std=c++0x -Wall -O3 -pg -DDEBUG")

# set default compilation and show details for debug build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif(NOT CMAKE_BUILD_TYPE)

add_library(malloc_count STATIC ${SOURCES})
target_link_libraries(malloc_count dl)
set_target_properties(malloc_count PROPERTIES LINKER_LANGUAGE C)
install(TARGETS malloc_count DESTINATION lib)

#add_executable(test-memprofile test-memprofile/test.cc)