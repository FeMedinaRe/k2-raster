/*  
 * Created by Fernando Silva on 27/03/20.
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
#include <temporal/at_k2_raster.hpp>
#include <utils/utils_data.hpp>

using namespace sdsl;
using namespace k2raster;
using testing::Types;

std::string data_folder = "../../test/temporal/Data/";
std::string output_filename = "../../test/Output/temporal.k2rt";
std::string inputs_filename = data_folder + "datasets.txt";

// Data params
size_t rows = 1000;
size_t cols = 1000;
size_t times = 10;
int min_value = 0;
int max_value = 5;

// Structure params
uint k1 = 4;
uint k2 = 2;
uint levels_k1 = 4;
uint plain_levels = 2;
size_t snap_freq = 4;
size_t scale_factor =0;
uint n_random_queries = 100;

// Random number generator
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<> dis_values(min_value, max_value);
std::uniform_int_distribution<> dis_row(0, rows-1);
std::uniform_int_distribution<> dis_col(0, cols-1);
std::uniform_int_distribution<> dis_time(0, times-1);

size_t n_changes_time = 100;


// TEST CLASS
template<class T>
class test_ath_k2_raster : public ::testing::Test { };


typedef Types<
        k2raster::athk2r_type
//        k2raster::atpk2r_type
> Implementations;

TYPED_TEST_CASE(test_ath_k2_raster, Implementations);

// TEST Encode - Create k2_raster temporal (in-memory)
TYPED_TEST(test_ath_k2_raster, CreateInMemory){

    // Build structure
    TypeParam k2raster(inputs_filename, data_folder, snap_freq, scale_factor);

    // Check values
    ASSERT_TRUE(k2raster.check(inputs_filename, data_folder, scale_factor));

    // Store to file
    ASSERT_TRUE(store_to_file(k2raster, output_filename));
}

// TEST File - Load the k2-raster from a file
TYPED_TEST(test_ath_k2_raster, LoadFromFile){
    // Load structure
    TypeParam k2raster;
    std::ifstream output_file(output_filename);
    k2raster.load(output_file);
    output_file.close();

    // Check values
    ASSERT_TRUE(k2raster.check(inputs_filename, data_folder, scale_factor));
}

// TEST Query - Random Get_Cell
TYPED_TEST(test_ath_k2_raster, QueryGet){
    TypeParam k2raster;
    std::ifstream output_file(output_filename);
    k2raster.load(output_file);
    output_file.close();

    // Get 'n_random_queries' positions
    for (auto n = 0; n < n_random_queries; n++) {
        auto row = dis_row(gen);
        auto col = dis_col(gen);
        auto t = dis_time(gen);

        // Original data
        std::vector<int> values;
        std::string file_name = data_folder + "dataset_" + std::to_string(t) + ".bin";
        read_input_data(file_name, rows, cols, values, scale_factor);

        // Compare values
        ASSERT_EQ(values[row * cols + col], k2raster.get_cell(row, col, t));
    }

    // Bad position
    auto row = rows + 10;
    auto col = cols + 10;
    auto t = dis_time(gen);
    ASSERT_EQ(0, k2raster.get_cell(row, col, t));
}

// TEST Query - Random Get_Cells_by_Value
TYPED_TEST(test_ath_k2_raster, QueryGetCellsByValue){
    TypeParam k2raster;
    std::ifstream output_file(output_filename);
    k2raster.load(output_file);
    output_file.close();

    // Run 'n_random_queries' queries
    for (auto n = 0; n < n_random_queries; n++) {
        auto xini = dis_row(gen);
        auto yini = dis_col(gen);
        auto tini = dis_time(gen);
        auto valini = dis_values(gen);
        std::uniform_int_distribution<> dis_max_row(xini, rows-1);
        std::uniform_int_distribution<> dis_max_col(yini, cols-1);
        std::uniform_int_distribution<> dis_max_time(tini, times-1);
        std::uniform_int_distribution<> dis_max_val(valini, max_value);
        auto xend = dis_max_row(gen);
        auto yend = dis_max_col(gen);
        auto tend = dis_max_time(gen);
        auto valend = dis_max_val(gen);

//        xini = 310;
//        yini = 998;
//        tini = 1;
//        valini = 5;
//        xend = 310;
//        yend = 998;
//        tend = 1;
//        valend = 5;

        // Run query
//        std::cout << "Query: (" << xini << ", " << yini << ")-(" << xend << ", " << yend << ")";
//        std::cout << "-- ["<<valini<<"-"<<valend<<"] a times [" << tini <<"-"<<tend<<"]"<<std::endl;
        std::vector<std::vector<std::pair<size_t, size_t>>> result(tend - tini + 1);
        size_t n_cells = k2raster.get_cells_by_value(xini, xend, yini, yend, valini, valend, tini, tend, result);

        // Compare n_cells with the number of cells in result array.
        {
            size_t n_result = 0;
            for (auto r : result) {
                n_result += r.size();
            }
            ASSERT_EQ(n_cells, n_result);
        }

        // Check values  (cell by cell)
        int value;
        for (auto t = tini; t <= tend; t++) {
            for (const auto &cell : result[t - tini]) {
                value = k2raster.get_cell(cell.first, cell.second, t);
//                std::cout << "Cell (" << cell.first << ", " << cell.second << ") at time " << t << " gets value " << value << std::endl;
                ASSERT_LE(valini, value);
                ASSERT_GE(valend, value);
            } // END FOR cell
        } // END FOR time

//        std::cout << "---------------------------" << std::endl;
        // Count number of cells of the region with a value within the region
        size_t cells_region = 0;
        for (auto t = tini; t <= tend; t++) {
            for (auto x = xini; x <= xend; x++) {
                for (auto y = yini; y <= yend; y++) {
                    value = k2raster.get_cell(x, y, t);
                    if (value >= valini && value <= valend) {
//                        std::cout << "Cell (" << x << ", " << y << ") at time " << t << " gets value " << value << std::endl;
                        cells_region++;
                    }
                } // END FOR y
            } // END FOR x
        } // END FOR t
        ASSERT_EQ(cells_region, n_cells);
    } // END FOR queries
}

// TEST Query - Random Get_Values_Window
TYPED_TEST(test_ath_k2_raster, QueryGetValuesWindow){
    TypeParam k2raster;
    std::ifstream output_file(output_filename);
    k2raster.load(output_file);
    output_file.close();

    // Run 'n_random_queries' queries
    for (auto n = 0; n < n_random_queries; n++) {
        auto xini = dis_row(gen);
        auto yini = dis_col(gen);
        auto tini = dis_time(gen);
        std::uniform_int_distribution<> dis_max_row(xini, rows-1);
        std::uniform_int_distribution<> dis_max_col(yini, cols-1);
        std::uniform_int_distribution<> dis_max_time(tini, times-1);
        auto xend = dis_max_row(gen);
        auto yend = dis_max_col(gen);
        auto tend = dis_max_time(gen);

        // Run query
        std::vector<std::vector<int>> result(tend - tini + 1);
        size_t n_cells = k2raster.get_values_window(xini, xend, yini, yend, tini, tend, result);
        ASSERT_EQ((xend - xini + 1)*(yend - yini + 1)*(tend - tini + 1), n_cells);

        // Check result (cell by cell)
        for (auto t = tini; t <= tend; t++) {
            size_t cells_region = 0;
            for (auto x = xini; x <= xend; x++) {
                for (auto y = yini; y <= yend; y++) {
                    ASSERT_EQ(k2raster.get_cell(x, y, t), result[t - tini][cells_region++]);
                } // END FOR y
            } // END FOR x
        } // END FOR t
    } // END FOR n queries
}
//
// TEST Query - Random Check_Values_Window (strong)
TYPED_TEST(test_ath_k2_raster, QueryCheckValuesWindowStrong){
    TypeParam k2raster;
    std::ifstream output_file(output_filename);
    k2raster.load(output_file);
    output_file.close();

    // Run 'n_random_queries' queries
    for (auto n = 0; n < n_random_queries; n++) {
        auto xini = dis_row(gen);
        auto yini = dis_col(gen);
        auto tini = dis_time(gen);
        auto valini = dis_values(gen);
        std::uniform_int_distribution<> dis_max_row(xini, rows-1);
        std::uniform_int_distribution<> dis_max_col(yini, cols-1);
        std::uniform_int_distribution<> dis_max_time(tini, times-1);
        std::uniform_int_distribution<> dis_max_val(valini, max_value);
        auto xend = dis_max_row(gen);
        auto yend = dis_max_col(gen);
        auto tend = dis_max_time(gen);
        auto valend = dis_max_val(gen);

        // Run query
        bool result = k2raster.check_values_window(xini, xend, yini, yend, valini, valend, tini, tend, true);

        // Count number of cells of the region with a value within the region
        bool found_no_valid_cell = false;
        int val;
        for (auto t = tini; t <= tend; t++) {
            for (auto x = xini; x <= xend; x++) {
                for (auto y = yini; y <= yend; y++) {
                    val = k2raster.get_cell(x, y, t);
                    if (!(val >= valini && val <= valend)) {
                        found_no_valid_cell = true;
                        break;
                    }
                } // END FOR y
                if (found_no_valid_cell) break;
            } // END FOR x
        } // END FOR t

        ASSERT_TRUE((!found_no_valid_cell && result) || (found_no_valid_cell && !result));
    } // END FOR n queries
}


// TEST Query - Random Check_Values_Window (weak)
TYPED_TEST(test_ath_k2_raster, QueryCheckValuesWindowWeak){
    TypeParam k2raster;
    std::ifstream output_file(output_filename);
    k2raster.load(output_file);
    output_file.close();

    // Run 'n_random_queries' queries
    for (auto n = 0; n < n_random_queries; n++) {
        auto xini = dis_row(gen);
        auto yini = dis_col(gen);
        auto tini = dis_time(gen);
        auto valini = dis_values(gen);
        std::uniform_int_distribution<> dis_max_row(xini, rows-1);
        std::uniform_int_distribution<> dis_max_col(yini, cols-1);
        std::uniform_int_distribution<> dis_max_time(tini, times-1);
        std::uniform_int_distribution<> dis_max_val(valini, max_value);
        auto xend = dis_max_row(gen);
        auto yend = dis_max_col(gen);
        auto tend = dis_max_time(gen);
        auto valend = dis_max_val(gen);

//        xini = 869;
//        yini = 27;
//        tini = 0;
//        valini = 2;
//        xend = 869;
//        yend = 32;
//        tend = 0;
//        valend = 2;

        // Run query (weak)
//        std::cout << "Query: (" << xini << ", " << yini << ")-(" << xend << ", " << yend << ")";
//        std::cout << "-- ["<<valini<<"-"<<valend<<"] a times [" << tini <<"-"<<tend<<"]"<<std::endl;
        bool result = k2raster.check_values_window(xini, xend, yini, yend, valini, valend, tini, tend, false);

        // Count number of cells of the region with a value within the region
        bool found_valid_cell = false;
        int val;
        for (auto t = tini; t <= tend; t++) {
            for (auto x = xini; x <= xend; x++) {
                for (auto y = yini; y <= yend; y++) {
                    val = k2raster.get_cell(x, y, t);
                    if (val >= valini && val <= valend) {
                        found_valid_cell = true;
                        break;
                    }
                } // END FOR y
                if (found_valid_cell) break;
            } // END FOR x
        } // END FOR t

        ASSERT_TRUE((!found_valid_cell && !result) || (found_valid_cell && result));
    } // END FOR n queries
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    // Create data
    {
        // Main file
        std::ofstream main_file(data_folder + "datasets.txt");

        // Write common params
        main_file << rows  << " " << cols  << " ";                                              // Dataset size (rows and columns)
        main_file << k1  << " " << k2  << " " << levels_k1  << " " << plain_levels << std::endl;    // k2-raster params

        // Create random data
        std::vector<int> values = k2raster::create_random_raster<int>(rows, cols, min_value, max_value);
        for(size_t d = 0; d < times; d++) {
            std::string file_name = "dataset_" + std::to_string(d) + ".bin";
            std::ofstream output(data_folder + file_name);

            for (size_t c=0; c < n_changes_time; c++) {
                auto xini = dis_row(gen);
                auto yini = dis_col(gen);
                auto diff = dis_values(gen);
                std::uniform_int_distribution<> dis_max_row(xini, rows - 1);
                std::uniform_int_distribution<> dis_max_col(yini, cols - 1);
                auto xend = dis_max_row(gen);
                auto yend = dis_max_col(gen);

                size_t pos;
                for (auto x = xini; x <= xend; x++) {
                    for (auto y = yini; y <= yend; y++) {
                        pos = (x * cols + y);
                        values[pos] = values[pos] + diff;
                    } // END FOR
                } // END FOR x
            } // END FOR n_changes

            sdsl::serialize_vector(values, output);

            // Insert information in main file
            main_file << file_name << std::endl;                                                              // Filename
            output.close();
        }
        main_file.close();
    }

    return RUN_ALL_TESTS();
}