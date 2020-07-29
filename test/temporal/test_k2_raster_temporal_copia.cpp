/*  
 * Created by Fernando Silva on 9/01/19.
 *
 * Copyright (C) 2019-current-year, Fernando Silva, all rights reserved.
 *
 * 
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * DESCRIPTION
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

// External libraries
#include <gtest/gtest.h>

// Own libraries
#include <temporal/k2_raster_temporal.hpp>
#include <temporal/t_k2_raster.hpp>
#include <temporal/k2_raster_temporal_global.hpp>

#include <utils/utils_data.hpp>

using namespace sdsl;
using namespace k2raster;
using testing::Types;

std::string data_folder ="../../test/temporal/Data/";
std::string output_file_name ="../../test/Output/temporal";
size_t n_datasets = 10;
size_t size_x = 1000;
size_t size_y = 1000;
uint k1 = 4;
uint k2 = 2;
uint levels_k1 = 4;
uint plain_levels = 2;
size_t snap_freq = 4;
size_t scale_factor =0;


// TEST CLASS
template<class T>
class test_k2_raster_temporal : public ::testing::Test { };


typedef Types<
        k2_raster_temporal<>
> Implementations;

TYPED_TEST_CASE(test_k2_raster_temporal, Implementations);

// TEST Encode - Create k2_raster temporal (in-memory)
TYPED_TEST(test_k2_raster_temporal, CreateInMemory){
    std::string input_file_name = data_folder + "datasets.txt";

    // k2-raster Temporal
    {
        k2raster::k2_raster_temporal<> k2raster(input_file_name, data_folder, snap_freq, scale_factor);
        ASSERT_TRUE(k2raster.check(input_file_name, data_folder, scale_factor));
        ASSERT_TRUE(store_to_file(k2raster, output_file_name + ".k2rt"));
    }

    // k2-raster Temporal (with k2-raster heuristic)
    {
        k2raster::k2_raster_temporal<int, k2raster::k2_raster_heuristic<>> k2raster(input_file_name, data_folder,
                                                                                      snap_freq, 0);
        ASSERT_TRUE(k2raster.check(input_file_name, data_folder, scale_factor));
        ASSERT_TRUE(store_to_file(k2raster, output_file_name + ".k2rth"));
    }

    // t-k2-raster
    {
        k2raster::t_k2_raster<> k2raster(input_file_name, data_folder, snap_freq, 0);
        ASSERT_TRUE(k2raster.check(input_file_name, data_folder, scale_factor));
        ASSERT_TRUE(store_to_file(k2raster, output_file_name + ".tk2r"));
    }

    // k2-raster Temporal global
    {
        k2raster::k2_raster_temporal_global<> k2raster(input_file_name, data_folder, snap_freq, 0);
        ASSERT_TRUE(k2raster.check(input_file_name, data_folder, scale_factor));
        ASSERT_TRUE(store_to_file(k2raster, output_file_name + ".k2rg"));
    }
}

// TEST File - Load from a file
TYPED_TEST(test_k2_raster_temporal, LoadFromFile){
    std::string input_file_name = data_folder + "datasets.txt";

    // k2-raster Temporal
    {
        k2raster::k2_raster_temporal<> k2raster;
        std::ifstream output_file(output_file_name + ".k2rt");
        k2raster.load(output_file);
        output_file.close();
        ASSERT_TRUE(k2raster.check(input_file_name, data_folder, scale_factor));
    }

    // k2-raster Temporal
    {
//        k2raster::k2_raster_temporal<int, k2raster::k2_raster_heuristic<>> k2raster;
//        std::ifstream output_file(output_file_name + ".k2rth");
//        k2raster.load(output_file);
//        output_file.close();
//        ASSERT_TRUE(k2raster.check(input_file_name, data_folder, scale_factor));
    }

    // k2-raster Temporal
    {
        k2raster::t_k2_raster<> k2raster;
        std::ifstream output_file(output_file_name + ".tk2r");
        k2raster.load(output_file);
        output_file.close();
        ASSERT_TRUE(k2raster.check(input_file_name, data_folder, scale_factor));
    }

    // k2-raster Temporal
    {
        k2raster::k2_raster_temporal_global<> k2raster;
        std::ifstream output_file(output_file_name + ".k2rg");
        k2raster.load(output_file);
        output_file.close();
        ASSERT_TRUE(k2raster.check(input_file_name, data_folder, scale_factor));
    }


}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    // Create data
    {
        // Main file
        std::ofstream main_file(data_folder + "datasets.txt");

        // Create random data
        for(auto d = 0; d < n_datasets; d++) {
            std::string file_name = "dataset_" + std::to_string(d) + ".bin";
            std::ofstream output(data_folder + file_name);

            std::vector<uint> values = k2raster::create_random_raster<uint>(size_x, size_y);
            sdsl::serialize(values, output);

            // Insert information in main file
            main_file << file_name << " ";                                                              // Filename
            main_file << size_x  << " " << size_y  << " ";                                              // Dataset size (rows and columns)
            main_file << k1  << " " << k2  << " " << levels_k1  << " " << plain_levels << std::endl;    // k2-raster params

            output.close();
        }
        main_file.close();
    }

    return RUN_ALL_TESTS();
}