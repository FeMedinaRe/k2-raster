/*  
 * Created by Fernando Silva on 25/10/18.
 *
 * Copyright (C) 2018-current-year, Fernando Silva, all rights reserved.
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
#include <k2_raster_base.hpp>
#include <k2_raster_plain.hpp>


using namespace sdsl;
using namespace k2raster;
using testing::Types;

// Global Params
size_t rows = 800;
size_t cols = 800;
int min_value = -0;
int max_value = 10;
size_t n_random_queries = 100;
std::vector<int> values(rows*cols);

std::vector<ushort> ks = {4, 2};
std::vector<ushort> levels_k1 = {1, 2, 3};
std::vector<ushort> levels_plain = {0, 1, 2};

ushort k2 = 2;
std::string output_filename = "../../test/Output/testing.k2raster";


// Random number generator
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<> dis_values(min_value, max_value);
std::uniform_int_distribution<> dis_row(0, rows-1);
std::uniform_int_distribution<> dis_col(0, cols-1);

// TEST CLASS
template<class T>
class test_k2_raster_plain : public ::testing::Test { };



typedef Types<
        k2_raster_plain<>
> Implementations;

TYPED_TEST_CASE(test_k2_raster_plain, Implementations);

// TEST Encode - Create k2_raster_plain (in-memory)
TYPED_TEST(test_k2_raster_plain, CreateStructure){

    for(auto const &k1 : ks) {
        for(auto const &level_k1 : levels_k1) {
            for(auto const &level_plain : levels_plain) {
                // Run k2-raster
                TypeParam k2raster(values,rows, cols, k1, k2, level_k1, level_plain);

                // Check values
                ASSERT_TRUE(k2raster.check(values, rows, cols));

                // Store to file
                std::string filename = output_filename + std::to_string(k1) + "_" + std::to_string(k2) + "_" + std::to_string(level_k1) + "_" + std::to_string(level_plain);
                ASSERT_TRUE(store_to_file(k2raster, filename));
            } // END FOR level_plain
        } // END FOR level_k1
    } // END FOR k1
}

// TEST File - Load the k2-raster from a file
TYPED_TEST(test_k2_raster_plain, LoadFromFile){

    for(auto const &k1 : ks) {
        for(auto const &level_k1 : levels_k1) {
            for(auto const &level_plain : levels_plain) {

                TypeParam k2raster;
                std::string filename = output_filename + std::to_string(k1) + "_" + std::to_string(k2) + "_" + std::to_string(level_k1) + "_" + std::to_string(level_plain);
                std::ifstream output_file(filename);
                k2raster.load(output_file);
                output_file.close();

                ASSERT_TRUE(k2raster.check(values, rows, cols));
            } // END FOR level_plain
        } // END FOR level_k1
    } // END FOR k1
}

// TEST Query - Random Get_Cell
TYPED_TEST(test_k2_raster_plain, QueryGet){
    for(auto const &k1 : ks) {
        for(auto const &level_k1 : levels_k1) {
            for(auto const &level_plain : levels_plain) {

                TypeParam k2raster;
                std::string filename = output_filename + std::to_string(k1) + "_" + std::to_string(k2) + "_" + std::to_string(level_k1) + "_" + std::to_string(level_plain);
                std::ifstream output_file(filename);
                k2raster.load(output_file);
                output_file.close();

                // Get 'n_random_queries' positions
                for (size_t  n = 0; n < n_random_queries; n++) {
                    auto row = dis_row(gen);
                    auto col = dis_col(gen);
                    ASSERT_EQ(values[row * cols + col], k2raster.get_cell(row, col));
                }

                // Limits
                ASSERT_EQ(values[0], k2raster.get_cell(0, 0));
                ASSERT_EQ(values[0 + (cols-1)], k2raster.get_cell(0, (cols-1)));
                ASSERT_EQ(values[(rows-1) * cols], k2raster.get_cell((rows-1), 0));
                ASSERT_EQ(values[(rows-1) * cols + (cols-1)], k2raster.get_cell((rows-1), (cols-1)));

                // Bad position
                auto row = rows + 10;
                auto col = cols + 10;
                ASSERT_EQ(0, k2raster.get_cell(row, col));
            } // END FOR level_plain
        } // END FOR level_k1
    } // END FOR k1
}

TYPED_TEST(test_k2_raster_plain, QueryGetCellsByValue){

    for(auto const &k1 : ks) {
        for(auto const &level_k1 : levels_k1) {
            for(auto const &level_plain : levels_plain) {

                TypeParam k2raster;
                std::string filename = output_filename + std::to_string(k1) + "_" + std::to_string(k2) + "_" + std::to_string(level_k1) + "_" + std::to_string(level_plain);
                std::ifstream output_file(filename);
                k2raster.load(output_file);
                output_file.close();

                // Run 'n_random_queries' queries
                for (size_t  n = 0; n < n_random_queries; n++) {
                    auto xini = dis_row(gen);
                    auto yini = dis_col(gen);
                    auto valini = dis_values(gen);
                    std::uniform_int_distribution<> dis_max_row(xini, rows-1);
                    std::uniform_int_distribution<> dis_max_col(yini, cols-1);
                    std::uniform_int_distribution<> dis_max_val(valini, max_value);
                    auto xend = dis_max_row(gen);
                    auto yend = dis_max_col(gen);
                    auto valend = dis_max_val(gen);

                    // Run query
                    std::vector<std::pair<size_t, size_t>> result;
                    size_t n_cells = k2raster.get_cells_by_value(xini, xend, yini, yend, valini, valend, result);

                    // Check result (cell by cell)
                    int val;
                    for (auto cell : result) {
                        val = values[cell.first * cols + cell.second];
                        ASSERT_LE(valini, val);
                        ASSERT_GE(valend, val);
                    }

                    // Count number of cells of the region with a value within the region
                    size_t n_cells_real = 0;
                    for (auto x = xini; x <= xend; x++) {
                        for (auto y = yini; y <= yend; y++) {
                            val = values[x * cols + y];
                            if (val >= valini && val <= valend) {
                                n_cells_real++;
                            }
                        }
                    }
                    ASSERT_EQ(n_cells_real, n_cells);
                }
            } // END FOR level_plain
        } // END FOR level_k1
    } // END FOR k1
}

TYPED_TEST(test_k2_raster_plain, QueryGetValuesWindow){
    for(auto const &k1 : ks) {
        for(auto const &level_k1 : levels_k1) {
            for(auto const &level_plain : levels_plain) {

                TypeParam k2raster;
                std::string filename = output_filename + std::to_string(k1) + "_" + std::to_string(k2) + "_" + std::to_string(level_k1) + "_" + std::to_string(level_plain);
                std::ifstream output_file(filename);
                k2raster.load(output_file);
                output_file.close();

                // Run 'n_random_queries' queries
                for (size_t  n = 0; n < n_random_queries; n++) {
                    auto xini = dis_row(gen);
                    auto yini = dis_col(gen);
                    std::uniform_int_distribution<> dis_max_row(xini, rows-1);
                    std::uniform_int_distribution<> dis_max_col(yini, cols-1);
                    auto xend = dis_max_row(gen);
                    auto yend = dis_max_col(gen);

                    // Run query
                    std::vector<int> result;
                    size_t n_cells = k2raster.get_values_window(xini, xend, yini, yend, result);
                    ASSERT_EQ(n_cells, result.size());

                    // Check result (cell by cell)
                    size_t pos = 0;
                    for (auto x = xini; x <= xend; x++) {
                        for (auto y = yini; y <= yend; y++) {
                            ASSERT_EQ(values[x * cols + y], result[pos++]);
                        } // END FOR y
                    } // END FOR x
                } // END FOR n queries
            } // END FOR level_plain
        } // END FOR level_k1
    } // END FOR k1
}

TYPED_TEST(test_k2_raster_plain, QueryCheckValuesWindowStrong){
    for(auto const &k1 : ks) {
        for(auto const &level_k1 : levels_k1) {
            for(auto const &level_plain : levels_plain) {

                TypeParam k2raster;
                std::string filename = output_filename + std::to_string(k1) + "_" + std::to_string(k2) + "_" + std::to_string(level_k1) + "_" + std::to_string(level_plain);
                std::ifstream output_file(filename);
                k2raster.load(output_file);
                output_file.close();

                // Run 'n_random_queries' queries
                for (size_t  n = 0; n < n_random_queries; n++) {
                    auto xini = dis_row(gen);
                    auto yini = dis_col(gen);
                    auto valini = dis_values(gen);
                    std::uniform_int_distribution<> dis_max_row(xini, rows-1);
                    std::uniform_int_distribution<> dis_max_col(yini, cols-1);
                    std::uniform_int_distribution<> dis_max_val(valini, max_value);
                    auto xend = dis_max_row(gen);
                    auto yend = dis_max_col(gen);
                    auto valend = dis_max_val(gen);

                    // Run query
                    bool result = k2raster.check_values_window(xini, xend, yini, yend, valini, valend, true);

                    // Count number of cells of the region with a value within the region
                    bool found_no_valid_cell = false;
                    int val;
                    for (auto x = xini; x <= xend; x++) {
                        for (auto y = yini; y <= yend; y++) {
                            val = values[x * cols + y];
                            if (!(val >= valini && val <= valend)) {
                                found_no_valid_cell = true;
                                break;
                            }
                        } // END FOR y
                        if (found_no_valid_cell) break;
                    } // END FOR x

                    ASSERT_TRUE((!found_no_valid_cell && result) || (found_no_valid_cell && !result));
                } // END FOR n queries
            } // END FOR level_plain
        } // END FOR level_k1
    } // END FOR k1
}

TYPED_TEST(test_k2_raster_plain, QueryCheckValuesWindowWeak){
    for(auto const &k1 : ks) {
        for(auto const &level_k1 : levels_k1) {
            for(auto const &level_plain : levels_plain) {

                TypeParam k2raster;
                std::string filename = output_filename + std::to_string(k1) + "_" + std::to_string(k2) + "_" + std::to_string(level_k1) + "_" + std::to_string(level_plain);
                std::ifstream output_file(filename);
                k2raster.load(output_file);
                output_file.close();

                // Run 'n_random_queries' queries
                for (size_t  n = 0; n < n_random_queries; n++) {
                    auto xini = dis_row(gen);
                    auto yini = dis_col(gen);
                    auto valini = dis_values(gen);
                    std::uniform_int_distribution<> dis_max_row(xini, rows-1);
                    std::uniform_int_distribution<> dis_max_col(yini, cols-1);
                    std::uniform_int_distribution<> dis_max_val(valini, max_value);
                    auto xend = dis_max_row(gen);
                    auto yend = dis_max_col(gen);
                    auto valend = dis_max_val(gen);

                    // Run query (weak)
                    bool result = k2raster.check_values_window(xini, xend, yini, yend, valini, valend, false);

                    // Count number of cells of the region with a value within the region
                    bool found_valid_cell = false;
                    int val;
                    for (auto x = xini; x <= xend; x++) {
                        for (auto y = yini; y <= yend; y++) {
                            val = values[x * cols + y];
                            if (val >= valini && val <= valend) {
                                found_valid_cell = true;
                                break;
                            }
                        } // END FOR y
                        if (found_valid_cell) break;
                    } // END FOR x

                    ASSERT_TRUE((!found_valid_cell && !result) || (found_valid_cell && result));
                } // END FOR n queries
            } // END FOR level_plain
        } // END FOR level_k1
    } // END FOR k1
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    for (uint r = 0; r < rows; r++) {
        for (uint c = 0; c < cols; c++) {
            values[r * cols + c] = dis_values(gen);
        }
    }
    return RUN_ALL_TESTS();
}