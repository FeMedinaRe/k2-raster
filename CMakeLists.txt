cmake_minimum_required(VERSION 2.8.7)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)

project(k2-raster_sdsl)
set(PROJECT_VENDOR "Fernando Silva-Coira")
set(CMAKE_PROJECT_NAME "k2-raster_sdsl")
set(PROJECT_CONTACT "fernando.silva@udc.es")

set(CMAKE_BUILD_TYPE "Release")

# Set output folder
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# C++11 compiler Check
if(NOT CMAKE_CXX_COMPILER_VERSION) # work around for cmake versions smaller than 2.8.10
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
endif()
if(CMAKE_CXX_COMPILER MATCHES ".*clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_COMPILER_IS_CLANGXX 1)
endif()
if( (CMAKE_COMPILER_IS_GNUCXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.7) OR
(CMAKE_COMPILER_IS_CLANGXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 3.2))
    message(FATAL_ERROR "Your C++ compiler does not support C++11. Please install g++ 4.7 (or greater) or clang 3.2 (or greater)")
else()
    message(STATUS "Compiler is recent enough to support C++11.")
endif()
if( CMAKE_COMPILER_IS_GNUCXX )
    append_cxx_compiler_flags("-std=c++11 -Wall -Wextra " "GCC" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-msse4.2 -O3 -ffast-math -funroll-loops -fno-omit-frame-pointer -g" "GCC" CMAKE_CXX_FLAGS_RELEASE)
else()
    append_cxx_compiler_flags("-std=c++11" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-stdlib=libc++" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-msse4.2 -O3  -ffast-math -funroll-loops -DNDEBUG" "CLANG" CMAKE_CXX_FLAGS_RELEASE)
endif()
include(CheckSSE)
FindSSE ()
if( SSE4_2_FOUND )
    if( CMAKE_COMPILER_IS_GNUCXX )
        append_cxx_compiler_flags("-msse4.2" "GCC" CMAKE_CXX_FLAGS)
    else()
        append_cxx_compiler_flags("-msse4.2" "CLANG" CMAKE_CXX_FLAGS)
    endif()
    message(STATUS "CPU does support SSE4.2.")
else()
    message(STATUS "CPU does NOT support SSE4.2")
endif()

add_subdirectory(external/sdsl-lite)

include_directories("${PROJECT_SOURCE_DIR}/external/sdsl-lite/external/googletest/include"
        ${CMAKE_HOME_DIRECTORY}/include
        ${CMAKE_BINARY_DIR}/external/sdsl-lite/include
        ${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort/include/
        ${CMAKE_HOME_DIRECTORY}//external/sdsl-lite/external/googletest/include
        ${CMAKE_HOME_DIRECTORY}/external/sdsl-lite/external/FastPFor/headers/
        ${STXXL_INCLUDE_DIRS}
        )

# Set source files
set(SOURCES
        include/k2_raster_base.hpp
        include/k2_raster.hpp
        include/k2_raster_plain.hpp
        include/k2_raster_heuristic.hpp
        include/temporal/k2_raster_temporal_base.hpp
        include/temporal/k2_raster_temporal.hpp
        include/temporal/t_k2_raster.hpp
        include/temporal/at_k2_raster.hpp
        include/temporal/k2_raster_temporal_global.hpp
        include/temporal/helpers/t_k2_raster_log.hpp
        include/temporal/helpers/t_k2_raster_plain_log.hpp
        include/temporal/helpers/t_k2_raster_heuristic_log.hpp
        include/temporal/helpers/k2_raster_heuristic_log.hpp
        include/utils/utils_heuristic.hpp
        include/utils/utils_query.hpp
        include/utils/query/query.hpp
        include/utils/dac_vector.hpp
        include/utils/utils_data.hpp
        include/utils/utils_time.hpp
        )


# k2-raster
add_executable(encode_k2r src/k2_raster_encode.cpp ${SOURCES})
target_link_libraries(encode_k2r sdsl divsufsort divsufsort64)

add_executable(check_k2r src/k2_raster_check.cpp ${SOURCES})
target_link_libraries(check_k2r sdsl divsufsort divsufsort64)

add_executable(get_cell_k2r src/k2_raster_get_cell.cpp ${SOURCES})
target_link_libraries(get_cell_k2r sdsl divsufsort divsufsort64)

add_executable(get_cells_by_value_k2r src/k2_raster_get_cells_by_value.cpp ${SOURCES})
target_link_libraries(get_cells_by_value_k2r sdsl divsufsort divsufsort64)

add_executable(get_values_window_k2r src/k2_raster_get_values_window.cpp ${SOURCES})
target_link_libraries(get_values_window_k2r sdsl divsufsort divsufsort64)

add_executable(check_values_window_k2r src/k2_raster_check_values_window.cpp ${SOURCES})
target_link_libraries(check_values_window_k2r sdsl divsufsort divsufsort64)


# Temporal
add_executable(encode_tk2r src/temporal/k2_raster_temporal_encode.cpp ${SOURCES})
target_link_libraries(encode_tk2r sdsl divsufsort divsufsort64)

add_executable(get_cell_tk2r src/temporal/k2_raster_temporal_get_cell.cpp ${SOURCES})
target_link_libraries(get_cell_tk2r sdsl divsufsort divsufsort64)

add_executable(get_cells_by_value_tk2r src/temporal/k2_raster_temporal_get_cells_by_value.cpp ${SOURCES})
target_link_libraries(get_cells_by_value_tk2r sdsl divsufsort divsufsort64)

add_executable(get_values_window_tk2r src/temporal/k2_raster_temporal_get_values_window.cpp ${SOURCES})
target_link_libraries(get_values_window_tk2r sdsl divsufsort divsufsort64)

add_executable(space_tk2r src/temporal/k2_raster_temporal_space.cpp ${SOURCES})
target_link_libraries(space_tk2r sdsl divsufsort divsufsort64)

add_executable(info_tk2r src/temporal/k2_raster_temporal_info.cpp ${SOURCES})
target_link_libraries(info_tk2r sdsl divsufsort divsufsort64)

# Utils
add_executable(create_access_queries_tk2r src/utils/create_access_queries_tk2r.cpp ${SOURCES})
target_link_libraries(create_access_queries_tk2r sdsl divsufsort divsufsort64)

add_executable(create_region_queries_tk2r src/utils/create_region_queries_tk2r.cpp ${SOURCES})
target_link_libraries(create_region_queries_tk2r sdsl divsufsort divsufsort64)

add_executable(create_window_queries_tk2r src/utils/create_window_queries_tk2r.cpp ${SOURCES})
target_link_libraries(create_window_queries_tk2r sdsl divsufsort divsufsort64)

add_executable(create_fixed_region_queries_tk2r src/utils/create_fixed_region_queries_tk2r.cpp ${SOURCES})
target_link_libraries(create_fixed_region_queries_tk2r sdsl divsufsort divsufsort64)


# TEST
add_executable(test_k2_raster test/test_k2_raster.cpp ${SOURCES})
target_link_libraries(test_k2_raster sdsl divsufsort divsufsort64 gtest)

add_executable(test_k2_raster_plain test/test_k2_raster_plain.cpp ${SOURCES})
target_link_libraries(test_k2_raster_plain sdsl divsufsort divsufsort64 gtest)

add_executable(test_k2_raster_heuristic test/test_k2_raster_heuristic.cpp ${SOURCES})
target_link_libraries(test_k2_raster_heuristic sdsl divsufsort divsufsort64 gtest)

# TEST (temporal)
add_executable(test_k2_raster_temporal test/temporal/test_k2_raster_temporal.cpp ${SOURCES})
target_link_libraries(test_k2_raster_temporal sdsl divsufsort divsufsort64 gtest)

add_executable(test_k2_raster_temporal_h test/temporal/test_k2_raster_temporal_h.cpp ${SOURCES})
target_link_libraries(test_k2_raster_temporal_h sdsl divsufsort divsufsort64 gtest)

add_executable(test_t_k2_raster test/temporal/test_t_k2_raster.cpp ${SOURCES})
target_link_libraries(test_t_k2_raster sdsl divsufsort divsufsort64 gtest)

add_executable(test_at_k2_raster test/temporal/test_at_k2_raster.cpp ${SOURCES})
target_link_libraries(test_at_k2_raster sdsl divsufsort divsufsort64 gtest)

add_executable(test_ath_k2_raster test/temporal/test_ath_k2_raster.cpp ${SOURCES})
target_link_libraries(test_ath_k2_raster sdsl divsufsort divsufsort64 gtest)

# TEST (Dac)
add_executable(test_dac test/test_new_DAC.cpp ${SOURCES})
target_link_libraries(test_dac sdsl divsufsort divsufsort64 gtest)
