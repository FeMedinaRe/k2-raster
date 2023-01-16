/*  
 * Created by Fernando Silva on 05/01/22.
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
#include <k2_raster.hpp>


using namespace sdsl;
using namespace k2raster;
using testing::Types;

// Global Params
size_t min_rows = 2100;
size_t min_cols = 2100;
size_t max_rows = 4000;
size_t max_cols = 4000;
int min_value = -1000;
int max_value = 1000;

ushort k1 = 4;
ushort k2 = 2;
ushort level_k1 = 2;
std::string output_filename = "../../test/Output/testing.k2raster";

// Operations
static const k2raster::OperationRaster operations[] = {k2raster::OperationRaster::OPERATION_SUM,
                                                       k2raster::OperationRaster::OPERATION_SUBT,
                                                       k2raster::OperationRaster::OPERATION_MULT};

static const k2raster::OperationZonalRaster operationsZonal[] = {k2raster::OperationZonalRaster::OPERATION_ZONAL_SUM};

static const uint zones[] = {10, 100, 200, 500};


// Random number generator
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<> dis_values(min_value, max_value);
std::uniform_int_distribution<size_t> dis_rows(min_rows, max_rows);
std::uniform_int_distribution<size_t> dis_cols(min_cols, max_cols);

// TEST CLASS
template<class T>
class test_algebra_k2_raster : public ::testing::Test { };


typedef Types<
        k2_raster<>
> Implementations;

TYPED_TEST_CASE(test_algebra_k2_raster, Implementations);

// TEST algebra with same size
TYPED_TEST(test_algebra_k2_raster, AlgebraSameSize){

    // ************************* //
    // *** Create raster 1   *** //
    // ************************* //
    size_t rows = dis_rows(gen);
    size_t cols = dis_cols(gen);

    std::vector<int> values1(rows*cols);
    for (size_t r = 0; r < rows; r++) {
        for (size_t c = 0; c < cols; c++) {
            values1[r * cols + c] = dis_values(gen);
        }
    }
    TypeParam k2raster1(values1, rows, cols, k1, k2, level_k1); // create k2_raster
    ASSERT_TRUE(k2raster1.check(values1, rows, cols));          // Check values

    // ************************* //
    // *** Create raster 2   *** //
    // ************************* //
    // With same real size

    std::vector<int> values2(rows*cols);
    for (size_t r = 0; r < rows; r++) {
        for (size_t c = 0; c < cols; c++) {
            values2[r * cols + c] = dis_values(gen);
        }
    }
    TypeParam k2raster2(values2, rows, cols, k1, k2, level_k1); // create k2_raster
    ASSERT_TRUE(k2raster2.check(values2, rows, cols));          // Check values

    // ************************* //
    // *** Algebra same size *** //
    // ************************* //
    for (auto op : operations) {
        TypeParam result_sum(k2raster1, k2raster2, op);         // Run Algebra

        // Check result
        size_t n_rows = result_sum.get_n_rows();
        size_t n_cols = result_sum.get_n_cols();

        ASSERT_EQ(n_rows, rows);
        ASSERT_EQ(n_cols, cols);

        for (size_t r = 0; r < n_rows; r++) {
            for (size_t c = 0; c < n_cols; c++) {
                int result1 = 0 ;
                switch (op) {
                    case k2raster::OperationRaster::OPERATION_SUM:
                        result1 = values1[r * n_cols + c] + values2[r * n_cols + c];
                        break;
                    case k2raster::OperationRaster::OPERATION_SUBT:
                        result1 = values1[r * n_cols + c] - values2[r * n_cols + c];
                        break;
                    case k2raster::OperationRaster::OPERATION_MULT:
                        result1 = values1[r * n_cols + c] * values2[r * n_cols + c];
                        break;
                }
                ASSERT_EQ(result_sum.get_cell(r, c), result1);
            } // END n_rows
        } // END n_cols
    }
}


// TEST algebra with different size
TYPED_TEST(test_algebra_k2_raster, AlgebraDifferentSize){

    // ************************* //
    // *** Create raster 1   *** //
    // ************************* //
    size_t rows = dis_rows(gen);
    size_t cols = dis_cols(gen);

    std::vector<int> values1(rows*cols);
    for (size_t r = 0; r < rows; r++) {
        for (size_t c = 0; c < cols; c++) {
            values1[r * cols + c] = dis_values(gen);
        }
    }
    TypeParam k2raster1(values1, rows, cols, k1, k2, level_k1); // create k2_raster
    ASSERT_TRUE(k2raster1.check(values1, rows, cols));          // Check values

    // ************************* //
    // *** Create raster 2   *** //
    // ************************* //
    // Different size
    size_t rows2 = dis_rows(gen);
    size_t cols2 = dis_cols(gen);

    std::vector<int> values2(rows2*cols2);
    for (size_t r = 0; r < rows2; r++) {
        for (size_t c = 0; c < cols2; c++) {
            values2[r * cols2 + c] = dis_values(gen);
        }
    }
    TypeParam k2raster2(values2, rows2, cols2, k1, k2, level_k1); // create k2_raster
    ASSERT_TRUE(k2raster2.check(values2, rows2, cols2));          // Check values


    // ****************************** //
    // *** Algebra different size *** //
    // ****************************** //
    for (auto op : operations) {
        TypeParam result_sum(k2raster1, k2raster2, op);         // Run Algebra

        // Check result
        size_t n_rows = result_sum.get_n_rows();
        size_t n_cols = result_sum.get_n_cols();

        ASSERT_EQ(n_rows, std::min(rows, rows2));
        ASSERT_EQ(n_cols, std::min(cols, cols2));

        for (size_t r = 0; r < n_rows; r++) {
            for (size_t c = 0; c < n_cols; c++) {
                int result1 = 0 ;
                switch (op) {
                    case k2raster::OperationRaster::OPERATION_SUM:
                        result1 = values1[r * cols + c] + values2[r * cols2 + c];
                        break;
                    case k2raster::OperationRaster::OPERATION_SUBT:
                        result1 = values1[r * cols + c] - values2[r * cols2 + c];
                        break;
                    case k2raster::OperationRaster::OPERATION_MULT:
                        result1 = values1[r * cols + c] * values2[r * cols2 + c];
                        break;
                }
                ASSERT_EQ(result_sum.get_cell(r, c), result1);
            } // END FOR n_rows
        } // END FOR n_cols
    } // END FOR operation
}

// TEST algebra - scalar operation
TYPED_TEST(test_algebra_k2_raster, AlgebraScalar){

    // ************************* //
    // *** Create raster 1   *** //
    // ************************* //
    size_t rows = dis_rows(gen);
    size_t cols = dis_cols(gen);

    std::vector<int> values1(rows*cols);
    for (size_t r = 0; r < rows; r++) {
        for (size_t c = 0; c < cols; c++) {
            values1[r * cols + c] = dis_values(gen);
        }
    }
    TypeParam k2raster1(values1, rows, cols, k1, k2, level_k1); // create k2_raster
    ASSERT_TRUE(k2raster1.check(values1, rows, cols));          // Check values


    // ****************************** //
    // *** Algebra Scalar         *** //
    // ****************************** //
    for (auto op : operations) {
        int scalar_value = dis_values(gen);
        TypeParam result_sum(k2raster1, scalar_value, op);         // Run Algebra

        // Check result
        size_t n_rows = result_sum.get_n_rows();
        size_t n_cols = result_sum.get_n_cols();

        ASSERT_EQ(n_rows, rows);
        ASSERT_EQ(n_cols, cols);

        for (size_t r = 0; r < n_rows; r++) {
            for (size_t c = 0; c < n_cols; c++) {
                int result1 = 0 ;
                switch (op) {
                    case k2raster::OperationRaster::OPERATION_SUM:
                        result1 = values1[r * cols + c] + scalar_value;
                        break;
                    case k2raster::OperationRaster::OPERATION_SUBT:
                        result1 = values1[r * cols + c] - scalar_value;
                        break;
                    case k2raster::OperationRaster::OPERATION_MULT:
                        result1 = values1[r * cols + c] * scalar_value;
                        break;
                }
                ASSERT_EQ(result_sum.get_cell(r, c), result1);
            } // END FOR n_rows
        } // END FOR n_cols
    } // END FOR operation
}

// TEST algebra - thresholding operation
TYPED_TEST(test_algebra_k2_raster, AlgebraThresholding){

    // ************************* //
    // *** Create raster 1   *** //
    // ************************* //
    size_t rows = dis_rows(gen);
    size_t cols = dis_cols(gen);

    std::vector<int> values1(rows*cols);
    for (size_t r = 0; r < rows; r++) {
        for (size_t c = 0; c < cols; c++) {
            values1[r * cols + c] = dis_values(gen);
        }
    }
    TypeParam k2raster1(values1, rows, cols, k1, k2, level_k1); // create k2_raster
    ASSERT_TRUE(k2raster1.check(values1, rows, cols));          // Check values


    // ****************************** //
    // *** Algebra Thresholding   *** //
    // ****************************** //
    int threshold_value = dis_values(gen);
    TypeParam result_sum(k2raster1, threshold_value);           // Run Algebra

    // Check result
    size_t n_rows = result_sum.get_n_rows();
    size_t n_cols = result_sum.get_n_cols();

    ASSERT_EQ(n_rows, rows);
    ASSERT_EQ(n_cols, cols);

    for (size_t r = 0; r < n_rows; r++) {
        for (size_t c = 0; c < n_cols; c++) {
            ASSERT_EQ(result_sum.get_cell(r, c), values1[r * cols + c] >= threshold_value);
        } // END FOR n_rows
    } // END FOR n_cols
}

// TEST algebra - Zonal
TYPED_TEST(test_algebra_k2_raster, AlgebraZonal){

    // ************************* //
    // *** Create raster 1   *** //
    // ************************* //
    size_t rows = dis_rows(gen);
    size_t cols = dis_cols(gen);

    std::vector<int> values1(rows * cols);
    for (size_t r = 0; r < rows; r++) {
        for (size_t c = 0; c < cols; c++) {
            values1[r * cols + c] = dis_values(gen);
        }
    }
    TypeParam k2raster1(values1, rows, cols, k1, k2, level_k1); // create k2_raster
    ASSERT_TRUE(k2raster1.check(values1, rows, cols));          // Check values



    std::vector<int> values_zone(rows * cols);
    for (const auto z : zones) {
        // ************************* //
        // *** Create Zone       *** //
        // ************************* //
        std::uniform_int_distribution<int> dis_zones(1, z);
        for (size_t r = 0; r < rows; r++) {
            for (size_t c = 0; c < cols; c++) {
                values_zone[r * cols + c] = dis_zones(gen);
            } // END FOR c
        } // END FOR r

        TypeParam k2raster_zone(values_zone, rows, cols, k1, k2, level_k1); // create k2_raster
        ASSERT_TRUE(k2raster_zone.check(values_zone, rows, cols));          // Check values

        // ****************************** //
        // *** Algebra - Zonal Sum *** //
        // ****************************** //
        for (auto op: operationsZonal) {
            TypeParam result_sum(k2raster1, k2raster_zone, op);             // Run Algebra

            // Check result
            size_t n_rows = result_sum.get_n_rows();
            size_t n_cols = result_sum.get_n_cols();

            ASSERT_EQ(n_rows, rows);
            ASSERT_EQ(n_cols, cols);

            // Calculate the value of each zone
            std::map<int, int> zonal_sum;
            for (size_t p = 0; p < values1.size(); p++) {
                if (zonal_sum.find(values_zone[p]) == zonal_sum.end()) {
                    zonal_sum[values_zone[p]] = values1[p];
                } else {
                    switch (op) {
                        case k2raster::OperationZonalRaster::OPERATION_ZONAL_SUM:
                            zonal_sum[values_zone[p]] += values1[p];
                            break;
                        default:
                            exit(-1);
                    } // END switch operation
                } // END IF zone
            } // END FOR p

            for (size_t r = 0; r < n_rows; r++) {
                for (size_t c = 0; c < n_cols; c++) {
                    int zone_value = zonal_sum[values_zone[r * n_cols + c]];
                    ASSERT_EQ(zone_value, result_sum.get_cell(r, c));
                } // END FOR c
            } // END FOR r
        } // END FOR zone
    } // END FOR operation
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}




