cmake_minimum_required(VERSION 2.8.4)
project(k2-raster)
set(PROJECT_VENDOR "Fernando Silva")
set(CMAKE_PROJECT_NAME "k2-raster")
set(PROJECT_CONTACT "fernando.silva@udc.es")
set(CMAKE_INSTALL_PREFIX /home/fsilva/software)
set(CMAKE_COMPILER_IS_GNUCXX)


# Set output folder
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


include_directories(./include/)
include_directories(./libs/RePair/)

set(SOURCES
        src/k2raster/k2-raster.cpp
        src/k2raster/plain/k2-raster-plain.cpp
        src/k2raster/plain/k2-raster-plain-DAC.cpp
        src/k2raster/plain/k2-raster-plain-VByte.cpp
        src/k2raster/compresessLeaves/k2-raster-CompressLeaves.cpp
        src/k2raster/compresessLeaves/k2-raster-CompressLeavesH.cpp
        src/k2raster/k2-raster-entropy.cpp
        src/k2raster/k2-raster-opt.cpp
        src/rtree/GetObjectsQueryBasic.cpp
        src/rtree/SpatialJoin/QueryStrategy/GetRoot.cpp
        src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategy.cpp
        src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyLeaves.cpp
        src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyQueue.cpp
        src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyQueueV2.cpp
        src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTree.cpp
        src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTreeV2.cpp
        src/rtree/SpatialJoin/QueryStrategy/GetObjectsQueryStrategyTreeV3.cpp
        src/rtree/SpatialJoin/Visitor/GetObjectsVisitorCells.cpp
        src/rtree/Statistic/GetObjectsQueryTopK.cpp
        src/rtree/Statistic/GetObjectsQueryMaxMBR.cpp
        src/rtree/Statistic/GetObjectsQueryMaxTopK.cpp
        src/rtree/Statistic/GetObjectsQueryMinTopK.cpp
        src/algebra/algebra.cpp
        src/k2raster/util/QNode.cpp
        src/k2raster/util/InfoMinSubMatrix.cpp
        libs/RePair/Coder/arrayg.cpp
        libs/RePair/Coder/basics.cpp
        libs/RePair/Coder/hash.cpp
        libs/RePair/Coder/heap.cpp
        libs/RePair/Coder/records.cpp
        libs/RePair/Coder/dictionary.cpp
        libs/RePair/Coder/IRePair.cpp
        libs/RePair/Coder/CRePair.cpp
        libs/RePair/RePair.cpp
        libs/RePair/rp-utils/LogSequence.cpp
        libs/RePair/rp-utils/DAC_VLS.cpp)

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
add_library(k2raster STATIC ${SOURCES})
install(TARGETS k2raster DESTINATION lib)
install(DIRECTORY "includes" DESTINATION includes/libk2raster)

#Add libraries
# Add libcds
set(LIBCDS_DIR "./libs/libcds2/")
include_directories(${LIBCDS_DIR})
include_directories(${LIBCDS_DIR}includes/)
add_subdirectory(${LIBCDS_DIR})

#Add libUtils
set(LIBUTILS_DIR "./libs/libUtils/")
include_directories(${LIBUTILS_DIR})
include_directories(${LIBUTILS_DIR}include/)
add_subdirectory(${LIBUTILS_DIR})

# Add libEncoders
#set(LIBENCODERS_DIR "./libs/libEncoders")
#include_directories(${LIBENCODERS_DIR})
#include_directories(${LIBENCODERS_DIR}include/)
#LINK_DIRECTORIES(./libs/libEncoders/lib)

# Add VByte
set(LIBVBYTE_DIR "./libs/libVByte/")
include_directories(${LIBVBYTE_DIR})
include_directories(${LIBVBYTE_DIR}include/)
add_subdirectory(${LIBVBYTE_DIR})


# Add libspatialindex
set(LIBSPTATIALINDEX_DIR "./libs/spatialindex/")
include_directories(${LIBSPTATIALINDEX_DIR})
include_directories(${LIBSPTATIALINDEX_DIR}include/)
include_directories(${LIBSPTATIALINDEX_DIR}include/spatialindex)
add_subdirectory(${LIBSPTATIALINDEX_DIR})


# Executables
add_executable(Encodek2Raster ./src/encodek2Raster.cpp)
add_executable(Infok2Raster ./src/infok2Raster.cpp)
add_executable(GetCellValuek2Raster ./src/testGetCellValue.cpp)
add_executable(GetCellsRangeValuek2Raster ./src/testGetCellsRangeValue.cpp)
add_executable(GetValuesWindowk2Raster ./src/testGetValuesWindow.cpp)
add_executable(CheckValuesk2Raster ./src/testCheckValues.cpp)
add_executable(GetCellValueFilesk2Raster ./src/testGetCellValueFiles.cpp)
add_executable(GetCellsRangeValueFilesk2Raster ./src/testGetCellsRangeValueFiles.cpp)
add_executable(GetValuesWindowFilesk2Raster ./src/testGetValuesWindowFiles.cpp)

# Executables RTree
#add_executable(JoinRTreek2Raster ./src/jointRTree.cpp)
add_executable(JoinTreek2Raster src/exec/SpatialJoin/jointRTree.cpp)

# Executables Algebra
add_executable(AlgebraOperation ./src/exec/algebra/algebraOperation.cpp)
add_executable(AlgebraOperationInMemory ./src/exec/algebra/algebraOperationInMemory.cpp)

# Executables Statistic
add_executable(GetMaxMBR ./src/exec/Statistics/statisticRTreeMaxMBR.cpp)
add_executable(GetMaxTopK ./src/exec/Statistics/statisticRTreeMaxTopK.cpp)
add_executable(GetMinTopK ./src/exec/Statistics/statisticRTreeMinTopK.cpp)

# Tools
add_executable(CreateValueQueries ./tools/createValueQueries.cpp)
add_executable(CreateValueQueriesPercentage ./tools/createValueQueriesPercentage.cpp)
add_executable(TransformRaster ./tools/transformRaster.cpp)

# TEST
add_executable(TestEncodek2Raster ./test/k2raster/testEncode.cpp)
add_executable(TestGetCellValue ./test/k2raster/testGetCellValue.cpp)
add_executable(TestGetCellsRangeValue ./test/k2raster/testGetCellsRangeValue.cpp)
add_executable(TestDecompress ./test/k2raster/testDecompress.cpp)

# PROBAS
add_executable(ProbarRePair ./src/exec/Probas/probar-repair.cpp)
target_link_libraries(ProbarRePair LINK_PUBLIC k2raster cds utils vbyte)


# Target libraries
target_link_libraries(Encodek2Raster LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(Infok2Raster LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(GetCellValuek2Raster LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(GetCellsRangeValuek2Raster LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(GetValuesWindowk2Raster LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(CheckValuesk2Raster LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(GetCellValueFilesk2Raster LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(GetCellsRangeValueFilesk2Raster LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(GetValuesWindowFilesk2Raster LINK_PUBLIC k2raster cds utils vbyte)

# Target libraries RTree
#target_link_libraries(JoinRTreek2Raster LINK_PUBLIC k2raster cds utils vbyte spatialindex)
target_link_libraries(JoinTreek2Raster LINK_PUBLIC k2raster cds utils vbyte spatialindex)

# Target libraries Algebra
target_link_libraries(AlgebraOperation LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(AlgebraOperationInMemory LINK_PUBLIC k2raster cds utils vbyte)

# Target libraries Statistic
target_link_libraries(GetMaxMBR LINK_PUBLIC k2raster cds utils vbyte spatialindex)
target_link_libraries(GetMaxTopK LINK_PUBLIC k2raster cds utils vbyte spatialindex)
target_link_libraries(GetMinTopK LINK_PUBLIC k2raster cds utils vbyte spatialindex)

# Target libraries Tools
target_link_libraries(CreateValueQueries LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(CreateValueQueriesPercentage LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(TransformRaster LINK_PUBLIC k2raster cds utils vbyte)

# Target libraries TEST
target_link_libraries(TestEncodek2Raster LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(TestGetCellValue LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(TestGetCellsRangeValue LINK_PUBLIC k2raster cds utils vbyte)
target_link_libraries(TestDecompress LINK_PUBLIC k2raster cds utils vbyte)