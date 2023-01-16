/*
 * Created by Fernando Silva on 08/10/21.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Program to execute map algebra. Thresholding operation
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
//#include <utils/utils_memory.hpp>
#include <utils/args/utils_args_algebra.hpp>

//**********************************************************************//
//************************** ALGEBRA ***********************************//
//**********************************************************************//
template<typename k2_raster_type,
        typename k2_raster_in_type = k2_raster_type,
        typename value_type=int>
void algebra_thresholding(k2_raster_in_type &raster1, int threshold_value, args_algebra_thr &args) {
    // Check if there is an operation with that code

    /*********************/
    /* Run operation     */
    /*********************/
    auto t1 = util::time::user::now(); // Start time
    k2_raster_type k2raster(raster1, threshold_value, args.set_memory_opt);
    auto t2 = util::time::user::now(); // End time

    /*********************/
    /* Print info        */
    /*********************/
    // Print time
    auto time = util::time::duration_cast<util::time::milliseconds>(t2 - t1);
    std::cout << "k2-raster - Thresholding time: " << time;
    std::cout << " milliseconds." << std::endl;

    // Print space
    size_t k2_raster_size = sdsl::size_in_bytes(k2raster);
    double ratio = ((double) k2_raster_size * 100.) / (k2raster.get_n_rows() * k2raster.get_n_cols() * sizeof(int));
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
        if (args.set_memory_opt) {
            // Memory opt destroy k2raster1 during execution. We need load them again.

            // Load raster 1
            std::ifstream input_file_1(args.raster1);
            assert(input_file_1.is_open() && input_file_1.good());
            raster1.load(input_file_1);
        }

        // Get size n_rows x n_cols
        uint n_rows = k2raster.get_n_rows();
        uint n_cols = k2raster.get_n_cols();

        size_t n_ones = 0;
        for (uint x = 0; x < n_rows; x++) {
            for (uint y = 0; y < n_cols; y++) {
                int result1 = raster1.get_cell(x, y);
                value_type val = k2raster.get_cell(x, y);
                if ((result1 >= threshold_value) != val) {
                    std::cout << std::endl;
                    std:: cout << "Found error at position (" << x << ", " << y <<"), expected " << result1 << " and get " << val << std::endl;
                    exit(-1);
                } // ENF IF result1
                if (val) n_ones++;
            } // END FOR y
        } // END FOR x
        std::cout << (n_ones * 100) / (n_rows * n_cols) << "% of ones" << std::endl;
        std::cout <<"ALL OK!!!!" << std::endl;
    } // END IF check
}


int main(int argc, char **argv) {

    /*********************/
    /* Reads params      */
    /*********************/
    args_algebra_thr args;
    parse_args_algebra_thr(argc, argv, args);

    /*********************/
    /* Run algebra       */
    /*********************/

    auto t1 = util::time::user::now(); // Start time
    for (auto r=0u; r < args.n_reps; r++) {
        // Load raster 1
        k2raster::k2_raster<> raster1;
        std::ifstream input_file_1(args.raster1);
        assert(input_file_1.is_open() && input_file_1.good());
        raster1.load(input_file_1);
        input_file_1.close();

        switch (args.k2_raster_type) {
            case k2raster::K2_RASTER_TYPE: {
                algebra_thresholding<k2raster::k2_raster<ushort>, k2raster::k2_raster<int>, ushort>(raster1, args.thr_value, args);
                break;
            }
            default:
                std::cout << "Invalid type " << args.k2_raster_type << ": " << std::endl;
                print_usage_algebra_thr(argv);
                exit(-1);
        }
    } // END FOR n_reps
    auto t2 = util::time::user::now(); // End time

    // Print total time
    auto time = util::time::duration_cast<util::time::milliseconds>(t2 - t1);
    std::cout << "[TOTAL] k2-raster - Thresholding time: " << time;
    std::cout << " milliseconds";
    std::cout << " in " << args.n_reps << " executions " << "(avg: " << time / args.n_reps << " milliseconds)" << std::endl;

    //util::print_memory_consumption();
    return 0;
}
