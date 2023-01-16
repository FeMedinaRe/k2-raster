/*
 * Created by Fernando Silva on 27/07/21.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Program to execute map algebra between two k2-raster files
 *
 * k2-raster is a compact data structure to represent raster data that
 * uses compressed space and offers indexing capabilities.
 * It uses min/max values for indexing and improving query performance.
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

// Own libraries
#include <k2_raster.hpp>
#include <utils/utils_time.hpp>
#include <utils/args/utils_args_algebra.hpp>


//**********************************************************************//
//************************** ALGEBRA ***********************************//
//**********************************************************************//
template<typename k2_raster_type,
        typename k2_raster_in_type = k2_raster_type,
        typename value_type=int>
void algebra(k2_raster_in_type &raster1, k2_raster_in_type &raster2, args_algebra &args) {

#ifndef NDEBUG
    // Check if there is a operation with that code
    std::cout <<"Operation: ";
    switch (args.operation) {
        case k2raster::OperationRaster::OPERATION_SUM:
            std::cout << "'Sum'" << std::endl;
            break;
        case k2raster::OperationRaster::OPERATION_SUBT:
            std::cout << "'Subtraction'" << std::endl;
            break;
        case k2raster::OperationRaster::OPERATION_MULT:
            std::cout << "'Multiplication'" << std::endl;
            break;
        default:
            std::cout << "No valid operation " << args.operation << std::endl;
            return;
    }
#endif

    /*********************/
    /* Run operation     */
    /*********************/
    auto t1 = util::time::user::now(); // Start time
    std::vector<value_type> data1;
    std::vector<typename k2_raster_in_type::value_type> data2;
    size_t n_rows1, n_cols1, n_rows2, n_cols2;
    raster1.decompress(data1, n_rows1, n_cols1);
    raster2.decompress(data2, n_rows2, n_cols2);

    size_t n_rows = std::min(n_rows1, n_rows2);
    size_t n_cols = std::min(n_cols1, n_cols2);
    //std::vector<value_type> result(n_rows * n_cols);
    size_t p = 0;

    switch (args.operation) {
        case k2raster::OperationRaster::OPERATION_SUM:
            for (size_t r = 0; r < n_rows; r++) {
                for (size_t c = 0; c < n_cols; c++) {
                    data1[p++] = (value_type)data1[r * n_cols1 + c] + (value_type)data2[r * n_cols2 + c];
                }
            }
            break;
        case k2raster::OperationRaster::OPERATION_SUBT:
            for (size_t r = 0; r < n_rows; r++) {
                for (size_t c = 0; c < n_cols; c++) {
                    data1[p++] = (value_type)data1[r * n_cols1 + c] - (value_type)data2[r * n_cols2 + c];
                }
            }
            break;
        case k2raster::OperationRaster::OPERATION_MULT:
            for (size_t r = 0; r < n_rows; r++) {
                for (size_t c = 0; c < n_cols; c++) {
                    data1[p++] = (value_type)data1[r * n_cols1 + c] * (value_type)data2[r * n_cols2 + c];
                }
            }
            break;
        default:
            std::cout << "No valid operation " << args.operation << std::endl;
            return;
    }

    k2_raster_type k2raster(data1, n_rows, n_cols, raster1.k1, raster1.k2, raster1.level_k1, 0);
    auto t2 = util::time::user::now(); // End time

    // Print time/space
    auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
    std::cout << "k2-raster build time: " << time;
    std::cout << " milliseconds." << std::endl;

    size_t k2_raster_size = sdsl::size_in_bytes(k2raster);
    double ratio = ((double)k2_raster_size * 100.) / (k2raster.get_n_rows() * k2raster.get_n_cols() * sizeof(int));
    std::cout << "k2-rater space:" << k2_raster_size << " bytes (" << ratio << "%)" << std::endl;

    /*********************/
    /* Save structure    */
    /*********************/
    if (!args.output_data.empty()) {
#ifndef NDEBUG
        std::cout << std::endl << "Storing k2-raster structure in file: " << args.output_data << std::endl;
#endif
        sdsl::store_to_file(k2raster, args.output_data);
#ifndef NDEBUG
        std::string file_name = std::string(args.output_data) + ".html";
        sdsl::write_structure<sdsl::format_type::HTML_FORMAT>(k2raster, file_name);
#endif
    }

    /*********************/
    /* Check             */
    /*********************/
    if (args.set_check){
        std::cout << "Checking k2-raster........" << std::endl;

        for (uint x = 0; x < n_rows; x++) {
            for (uint y = 0; y < n_cols; y++) {
                value_type result1 = 0 ;
                switch (args.operation) {
                    case k2raster::OperationRaster::OPERATION_SUM:
                        result1 = (value_type)raster1.get_cell(x, y) + (value_type)raster2.get_cell(x, y);
                        break;
                    case k2raster::OperationRaster::OPERATION_SUBT:
                        result1 = (value_type)raster1.get_cell(x, y) - (value_type)raster2.get_cell(x, y);
                        break;
                    case k2raster::OperationRaster::OPERATION_MULT:
                        result1 = (value_type)raster1.get_cell(x, y) * (value_type)raster2.get_cell(x, y);
                        break;
                }
                value_type result2 = k2raster.get_cell(x, y);
                if (result1 != result2) {
                    std::cout << "Found error at position (" << x << ", " << y  << "), ";
                    std::cout << "expected " << result1 << " and get " << k2raster.get_cell(x, y) << std::endl;
                    exit(-1);
                } // END IF result1
            } // END FOR y
        } // END FOR x
        std::cout << "ALL OK!!!!" << std::endl;
    } // END IF check
}


int main(int argc, char **argv) {

    /*********************/
    /* Reads params      */
    /*********************/
    args_algebra args;
    parse_args_algebra(argc, argv, args);

    /*********************/
    /* Run algebra       */
    /*********************/
    switch (args.k2_raster_type) {
        case k2raster::K2_RASTER_TYPE: {
            k2raster::k2_raster<> raster1, raster2;

            // Load raster 1
            std::ifstream input_file_1(args.raster1);
            assert(input_file_1.is_open() && input_file_1.good());
            raster1.load(input_file_1);

            // Load raster 2
            std::ifstream input_file_2(args.raster2);
            assert(input_file_2.is_open() && input_file_2.good());
            raster2.load(input_file_2);

            // Select if values are stored as integers or long integers
            bool more_than_int = false;
            long max_result;
            switch (args.operation) {
                case k2raster::OperationRaster::OPERATION_SUM:
                    max_result = (long)raster1.max_value + (long)raster2.max_value;
                    break;
                case k2raster::OperationRaster::OPERATION_SUBT:
                    max_result = (long)raster1.max_value - (long)raster2.max_value;
                    break;
                case k2raster::OperationRaster::OPERATION_MULT:
                    max_result = (long)raster1.max_value * (long)raster2.max_value;
                    break;
            }
            if (max_result >= std::numeric_limits<int>::max() || max_result <= std::numeric_limits<int>::min()) {
#ifndef NDEBUG
                std::cout << "Values are stored as long integers." << std::endl;
#endif
                more_than_int = true;
            }
#ifndef NDEBUG
            else {
                std::cout << "Values are stored as integers." << std::endl;
            }
#endif
            if (more_than_int) {
                // Store values as Long Integers
                algebra<k2raster::k2_raster<long>, k2raster::k2_raster<int>, long>(raster1, raster2, args);
            } else {
                // Store values as Integers
                algebra<k2raster::k2_raster<int>, k2raster::k2_raster<int>, int>(raster1, raster2, args);
            }
            break;
        }
        default:
            std::cout << "Invalid type " << args.k2_raster_type << ": " << std::endl;
            print_usage_algebra(argv);
            exit(-1);
    }

    return 0;
}