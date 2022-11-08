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
void algebra_scalar(k2_raster_type &raster1, const int scalar_value, std::string &output_data, k2raster::OperationRaster operation ,
                    bool set_check) {
#ifndef NDEBUG
    // Check if there is a operation with that code
    std::cout <<"Operation: ";
    switch (operation) {
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
            std::cout << "No valid operation " << operation << std::endl;
            return;
    }
#endif



    /*********************/
    /* Run operation     */
    /*********************/
    k2_raster_type k2raster;
    std::vector<int> data;
    size_t n_rows, n_cols;

    auto t1 = util::time::user::now(); // Start time
    raster1.decompress(data, n_rows, n_cols);

    for (size_t p = 0; p < n_rows * n_cols; p++) {
        switch (operation) {
            case k2raster::OperationRaster::OPERATION_SUM:
                data[p] += scalar_value;
                break;
            case k2raster::OperationRaster::OPERATION_SUBT:
                data[p] -= scalar_value;
                break;
            case k2raster::OperationRaster::OPERATION_MULT:
                data[p] *= scalar_value;
                break;
        } // END switch operation
    } // END FOR position

    k2raster = k2_raster_type(data, n_rows, n_cols, raster1.k1, raster1.k2, raster1.level_k1, 0);
    auto t2 = util::time::user::now(); // End time

    // Print time/space
    auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
    std::cout << "k2-raster - Scalar (de) time: " << time;
    std::cout << " milliseconds." << std::endl;

    size_t k2_raster_size = sdsl::size_in_bytes(k2raster);
    double ratio = (k2_raster_size * 100.) / (k2raster.get_n_rows() * k2raster.get_n_cols() * sizeof(int));
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
        std::cout << "Checking k2-raster........" << std::endl;
        for (size_t x = 0; x < n_rows; x++) {
            for (size_t y = 0; y < n_cols; y++) {
                int result1 = raster1.get_cell(x, y);
                switch (operation) {
                    case k2raster::OperationRaster::OPERATION_SUM:
                        result1 += scalar_value;
                        break;
                    case k2raster::OperationRaster::OPERATION_SUBT:
                        result1 -= scalar_value;
                        break;
                    case k2raster::OperationRaster::OPERATION_MULT:
                        result1 *= scalar_value;
                        break;
                }

                int val = k2raster.get_cell(x, y);
                if (result1 != val) {
                    std::cout << "Found error at position (" << x << ", " << y  << "), ";
                    std::cout << "expected " << result1 << " and get " << val << std::endl;
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
    args_algebra_scalar args;
    parse_args_algebra_scalar(argc, argv, args);


    /*********************/
    /* Run algebra       */
    /*********************/
    auto t1 = util::time::user::now(); // Start time
    for (auto r=0u; r < args.n_reps; r++) {
        switch (args.k2_raster_type) {
            case k2raster::K2_RASTER_TYPE: {
                k2raster::k2_raster<> raster1, raster2;

                // Load raster 1
                std::ifstream input_file_1(args.raster1);
                assert(input_file_1.is_open() && input_file_1.good());
                raster1.load(input_file_1);

                algebra_scalar<k2raster::k2_raster<>>(raster1, args.scalar_value, args.output_data,
                                                      static_cast<k2raster::OperationRaster>(args.operation),
                                                      args.set_check);
            }
                break;
            default:
                std::cout << "Invalid type " << args.k2_raster_type << ": " << std::endl;
                print_usage_algebra_scalar(argv);
                exit(-1);
        }
    } // END FOR n_reps
    auto t2 = util::time::user::now(); // End time

    // Print total time
    auto time = util::time::duration_cast<util::time::milliseconds>(t2 - t1);
    std::cout << "[TOTAL] k2-raster - Scalar (de) time: " << time;
    std::cout << " milliseconds";
    std::cout << " in " << args.n_reps << " executions " << "(avg: " << time / args.n_reps << " milliseconds)" << std::endl;

    return 0;
}