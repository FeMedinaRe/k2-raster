/*
 * Created by Fernando Silva on 28/10/21.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Program to execute zonal operation between two k2-raster files
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
        typename value_type=long,
        typename value_in_type=int>
void algebra_zonal(k2_raster_in_type &raster1, k2_raster_in_type &raster_zonal, const args_algebra_zonal &args) {
#ifndef NDEBUG
    // Check if there is an operation with that code
    std::cout << "Operation: ";
    switch (args.operation) {
        case k2raster::OperationZonalRaster::OPERATION_ZONAL_SUM:
            std::cout << "'Zonal Sum'" << std::endl;
            break;
        default:
            std::cout << "No valid operation << " << args.operation << std::endl;
            return;
    }
#endif

    /*********************/
    /* Run operation     */
    /*********************/
    auto t1 = util::time::user::now(); // Start time
    k2_raster_type k2raster(raster1, raster_zonal, args.operation, args.set_version_opt/*, args.set_memory_opt*/);
    auto t2 = util::time::user::now(); // End time

    /*********************/
    /* Print info        */
    /*********************/
    // Print time consumption
    auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
    std::cout << "k2-raster build time: " << time;
    std::cout << " milliseconds." << std::endl;

    // Print space
    size_t k2_raster_size = sdsl::size_in_bytes(k2raster);
    double ratio =
            ((double) k2_raster_size * 100.) / (k2raster.get_n_rows() * k2raster.get_n_cols() * sizeof(value_in_type));
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
            // Memory opt destroy k2raster1 and k2raster2 during execution. We need load them again.

            // Load raster 1
            std::ifstream input_file_1(args.raster1);
            assert(input_file_1.is_open() && input_file_1.good());
            raster1.load(input_file_1);

            // Load raster 2
            std::ifstream input_file_2(args.raster_zonal);
            assert(input_file_2.is_open() && input_file_2.good());
            raster_zonal.load(input_file_2);
        }

        size_t n_rows, n_cols, nz_rows, nz_cols;
        std::vector<value_in_type> zonal_values;
        raster_zonal.decompress(zonal_values, nz_rows, nz_cols);
        std::vector<value_in_type> raster_values;
        raster1.decompress(raster_values, n_rows, n_cols);

        // Calculate the value of each zone
        std::map<value_in_type, value_type> zonal_sum;
        size_t p, pz;
        for (size_t r = 0; r < n_rows; r++) {
            for(size_t c = 0; c < n_cols; c++) {
                p = r * n_cols + c;
                pz = r * nz_cols + c;
                if (zonal_sum.find(zonal_values[pz]) == zonal_sum.end()) {
                    zonal_sum[zonal_values[pz]] = raster_values[p];
                } else {
                    zonal_sum[zonal_values[pz]] += raster_values[p];
                }
            }
        }

        for (size_t x = 0; x < n_rows; x++) {
            for (size_t y = 0; y < n_cols; y++) {
                value_type zone_value = zonal_sum[raster_zonal.get_cell(x, y)];
                if (zone_value != k2raster.get_cell(x, y)) {
                    std::cout << "Found error at position (" << x << ", " << y  << "), ";
                    std::cout << "expected " << zone_value << " and get " << k2raster.get_cell(x, y) << std::endl;
                    exit(-1);
                } // END IF result1
            } // END FOR y
        } // END FOR x
        std::cout << "ALL OK!!!!" << std::endl;
    }
}


int main(int argc, char **argv) {

    /*********************/
    /* Reads params      */
    /*********************/
    args_algebra_zonal args;
    parse_args_algebra_zonal(argc, argv, args);

    /*********************/
    /* Run algebra       */
    /*********************/
    switch (args.k2_raster_type) {
        case k2raster::K2_RASTER_TYPE: {
            k2raster::k2_raster<int> raster1, rasterZ;

            // Load raster 1
            std::ifstream input_file_1(args.raster1);
            assert(input_file_1.is_open() && input_file_1.good());
            raster1.load(input_file_1);

            // Load raster 2
            std::ifstream input_file_zonal(args.raster_zonal);
            assert(input_file_zonal.is_open() && input_file_zonal.good());
            rasterZ.load(input_file_zonal);

            algebra_zonal<k2raster::k2_raster<long>, k2raster::k2_raster<int>, long, int>(raster1, rasterZ, args);
            break;
        }
        default:
            std::cout << "Invalid type " << args.k2_raster_type << ": " << std::endl;
            print_usage_algebra_zonal(argv);
            exit(-1);
    }

#ifdef SHOW_CONSUMPTION_MEMORY
    util::print_memory_consumption();
#endif
    return 0;
}