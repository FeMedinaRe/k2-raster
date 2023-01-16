/*
 * Created by Fernando Silva on 21/07/21.
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
#ifdef SHOW_CONSUMPTION_MEMORY
#include <utils/utils_memory.hpp>
#endif

//**********************************************************************//
//************************** ALGEBRA ***********************************//
//**********************************************************************//
template<typename k2_raster_type,
        typename k2_raster_in_type = k2_raster_type,
        typename value_type=int>
void algebra(k2_raster_in_type &raster1, k2_raster_in_type &raster2, args_algebra &args) {
#ifndef NDEBUG
    // Check if there is an operation with that code
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
            exit(-1);
    }
#endif

    /*********************/
    /* Run operation     */
    /*********************/
    auto t1 = util::time::user::now(); // Start time
    k2_raster_type k2raster(raster1, raster2, args.operation, args.set_memory_opt);
    auto t2 = util::time::user::now(); // End time

    /*********************/
    /* Print info        */
    /*********************/
    // Print time
    auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
    std::cout << "k2-raster build time: " << time;
    std::cout << " milliseconds." << std::endl;

    // Print space
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

#ifdef MIN_CONSUMPTION_MEMORY
        if (args.set_memory_opt) {
            // Memory opt destroy k2raster1 and k2raster2 during execution. We need load them again.

            // Load raster 1
            std::ifstream input_file_1(args.raster1);
            assert(input_file_1.is_open() && input_file_1.good());
            raster1.load(input_file_1);

            // Load raster 2
            std::ifstream input_file_2(args.raster2);
            assert(input_file_2.is_open() && input_file_2.good());
            raster2.load(input_file_2);
        }
#endif

        // Load raster size
        uint n_rows = k2raster.get_n_rows();
        uint n_cols = k2raster.get_n_cols();

        // Check cell by cell
        for (uint x = 0; x < n_rows; x++) {
            for (uint y = 0; y < n_cols; y++) {
                value_type result1 = 0 ;
                value_type result2 = 0 ;
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
                result2 = k2raster.get_cell(x, y);
                if (result1 != result2) {
                    std::cout << "Found error at position (" << x << ", " << y  << "), ";
                    std::cout << "expected " << (value_type)result1 << " and get " << (value_type)result2 << std::endl;
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
            long min_result, max_result;
            switch (args.operation) {
                case k2raster::OperationRaster::OPERATION_SUM:
                    min_result = (long)raster1.min_value + (long)raster2.min_value;
                    max_result = (long)raster1.max_value + (long)raster2.max_value;
                    break;
                case k2raster::OperationRaster::OPERATION_SUBT:
                    min_result = (long)raster1.min_value - (long)raster2.max_value;
                    max_result = (long)raster1.max_value - (long)raster2.min_value;
                    break;
                case k2raster::OperationRaster::OPERATION_MULT:
                    min_result = (long)raster1.min_value * (long)raster2.min_value;
                    max_result = (long)raster1.max_value * (long)raster2.max_value;
                    break;
            }
            if (max_result >= std::numeric_limits<int>::max() ||
                min_result <= std::numeric_limits<int>::min() ||
                max_result - min_result >= std::numeric_limits<int>::max()) {
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
        }
            break;
        default:
            std::cout << "Invalid type " << args.k2_raster_type << ": " << std::endl;
            print_usage_algebra(argv);
            exit(-1);
    }

#ifdef SHOW_CONSUMPTION_MEMORY
    util::print_memory_consumption();
#endif
    return 0;
}