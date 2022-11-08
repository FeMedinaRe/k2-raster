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

#include <k2_raster.hpp>
#include <utils/utils_time.hpp>
#include <utils/args/utils_args_algebra.hpp>


//**********************************************************************//
//************************** ALGEBRA ***********************************//
//**********************************************************************//
template<typename k2_raster_type>
void algebra_thresholding(k2_raster_type &raster1, int threshold_value, std::string &output_data,
                          bool set_check, uint n_reps=1, ushort option=1) {


    auto t1 = util::time::user::now(); // Start time
    size_t n_rows = raster1.get_n_rows();
    size_t n_cols = raster1.get_n_cols();
    std::vector<int> result(n_rows * n_cols, 0);
    /*********************/
    /* Run operation     */
    /*********************/
    k2_raster_type k2raster;
    for ( auto r=0u; r < n_reps; r++) {
        switch (option) {
            case 0: {
                std::vector<int> data1;
                size_t n_rows1, n_cols1;
                raster1.decompress(data1, n_rows1, n_cols1);

                for (size_t p = 0; p < n_rows * n_cols; p++) {
                    if (data1[p] >= threshold_value) {
                        result[p] = 1;
                    }
                    //result[p++] = data1[r * n_cols1 + c] >= threshold_value ? 1 : 0;
                }
                break;
            }
            case 1: {
                std::vector<std::pair<size_t, size_t>> result_points;
                raster1.get_cells_by_value(0, n_rows - 1, 0, n_cols - 1, threshold_value, raster1.max_value,
                                           result_points);

                for (auto point: result_points) {
                    result[point.first * n_cols + point.second] = 1;
                }
                break;
            }
            default:
                std::cout << "Invalid option " << option << ": " << std::endl;
                exit(-1);
        }

        k2raster = k2_raster_type(result, n_rows, n_cols, raster1.k1, raster1.k2, raster1.level_k1, 0);
    } // END FOR r (executions)
    auto t2 = util::time::user::now(); // End time

    // Print time/space
    auto time = util::time::duration_cast<util::time::milliseconds>(t2 - t1);
    std::cout << "k2-raster - Thresholding (de) time: " << time;
    std::cout << " milliseconds." << std::endl;

    size_t k2_raster_size = sdsl::size_in_bytes(k2raster);
    double ratio = ((double)k2_raster_size * 100.) / (k2raster.get_n_rows() * k2raster.get_n_cols() * sizeof(int));
    std::cout << "k2-rater space:" << k2_raster_size << " bytes (" << ratio << "%)" << std::endl;

    /*********************/
    /* Save structure    */
    /*********************/
    if (!output_data.empty()) {
#ifndef NDEBUG
        std::cout << std::endl << "Storing k2-raster structure in file: " << output_data << std::endl;
#endif
        sdsl::store_to_file(k2raster, output_data);
#ifndef NDEBUG
        std::string file_name = std::string(output_data) + ".html";
        sdsl::write_structure<sdsl::format_type::HTML_FORMAT>(k2raster, file_name);
#endif
    }

    /*********************/
    /* Check             */
    /*********************/
    if (set_check){
        std::cout << "Checking k2-raster........";
        size_t n_ones = 0;
        for (uint x = 0; x < n_rows; x++) {
            for (uint y = 0; y < n_cols; y++) {
                int result1 = raster1.get_cell(x, y);
                int val = k2raster.get_cell(x, y);
                if ((result1 >= threshold_value) != val) {
                    std::cout << std::endl;
                    std:: cout << "Found error at position (" << x << ", " << y <<", expected " << result1 << " and get " << val << std::endl;
                    exit(-1);
                }
                if (val) n_ones++;
            }
        }

        std::cout << (n_ones * 100) / (n_rows * n_cols) << "% of ones" << std::endl;
        std::cout <<"ALL OK!!!!" << std::endl;
    }
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
        switch (args.k2_raster_type) {
            case k2raster::K2_RASTER_TYPE: {
                k2raster::k2_raster<> raster1;

                // Load raster 1
                std::ifstream input_file_1(args.raster1);
                assert(input_file_1.is_open() && input_file_1.good());
                raster1.load(input_file_1);

                algebra_thresholding<k2raster::k2_raster<>>(raster1, args.thr_value, args.output_data, args.set_check,
                                                            args.n_reps, 1);
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
    std::cout << "[TOTAL] k2-raster - Thresholding (de) time: " << time;
    std::cout << " milliseconds";
    std::cout << " in " << args.n_reps << " executions " << "(avg: " << time / args.n_reps << " milliseconds)" << std::endl;
    return 0;
}