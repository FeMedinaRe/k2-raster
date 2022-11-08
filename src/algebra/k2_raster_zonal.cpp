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

//**********************************************************************//
//************************** ALGEBRA ***********************************//
//**********************************************************************//
template<typename k2_raster_type, typename t_value=int>
void algebra_zonal(k2_raster_type &raster1, k2_raster_type &raster_zonal, const std::string &output_data,
                   const k2raster::OperationZonalRaster &operation, bool set_check) {
#ifndef NDEBUG
    // Check if there is an operation with that code
    std::cout << "Operation: ";
    switch (operation) {
        case k2raster::OperationZonalRaster::OPERATION_ZONAL_SUM:
            std::cout << "'Zonal Sum'" << std::endl;
            break;
        default:
            std::cout << "No valid operation << " << operation << std::endl;
            return;
    }
#endif

    /*********************/
    /* Run operation     */
    /*********************/
    auto t1 = util::time::user::now(); // Start time
    k2_raster_type k2raster(raster1, raster_zonal, operation, false);
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
            ((double) k2_raster_size * 100.) / (k2raster.get_n_rows() * k2raster.get_n_cols() * sizeof(t_value));
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
        size_t n_rows, n_cols;
        std::vector<t_value> zonal_values;
        raster_zonal.decompress(zonal_values, n_rows, n_cols);
        std::vector<t_value> raster_values;
        raster1.decompress(raster_values, n_rows, n_cols);

        // Calculate the value of each zone
        std::map<t_value, t_value> zonal_sum;
        for (size_t p = 0; p < raster_values.size(); p++) {
            if (zonal_sum.find(zonal_values[p]) == zonal_sum.end()) {
                zonal_sum[zonal_values[p]] = raster_values[p];
            } else {
                zonal_sum[zonal_values[p]] += raster_values[p];
            }
        } // END FOR p

        for (size_t x = 0; x < n_rows; x++) {
            for (size_t y = 0; y < n_cols; y++) {
                t_value zone_value = zonal_sum[raster_zonal.get_cell(x, y)];
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
            k2raster::k2_raster<> raster1, rasterZ;

            // Load raster 1
            std::ifstream input_file_1(args.raster1);
            assert(input_file_1.is_open() && input_file_1.good());
            raster1.load(input_file_1);

            // Load raster 2
            std::ifstream input_file_zonal(args.raster_zonal);
            assert(input_file_zonal.is_open() && input_file_zonal.good());
            rasterZ.load(input_file_zonal);

            algebra_zonal<k2raster::k2_raster<>>(raster1, rasterZ, args.output_data,
                                                 static_cast<k2raster::OperationZonalRaster>(args.operation),
                                                 args.set_check);
            break;
        }
        default:
            std::cout << "Invalid type " << args.k2_raster_type << ": " << std::endl;
            print_usage_algebra_zonal(argv);
            exit(-1);
    }
    return 0;
}