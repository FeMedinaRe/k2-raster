/*  
 * Created by Fernando Silva on 9/01/19.
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


#include <temporal/k2_raster_temporal.hpp>
#include <temporal/t_k2_raster.hpp>
#include <temporal/k2_raster_temporal_global.hpp>
#include <temporal/at_k2_raster.hpp>
#include <utils/utils_time.hpp>

void print_help(char * argv0) {
    printf("Usage: %s <inputs_file> <input_path> <scale_factor>  <output_path> <type> <snap_freq> <set_check>\n",
           argv0);
}

template<typename k2_raster_temporal_type>
void run_encode(std::string inputs_filename, std::string input_path, ushort scale_factor, std::string output_filename, size_t  snap_freq, bool set_check) {

    /*********************/
    /* Encodes data      */
    /*********************/
    auto t1 = util::time::user::now(); // Start time
    k2_raster_temporal_type k2raster(inputs_filename, input_path, snap_freq, scale_factor);


    // Print time/space
    auto t2 = util::time::user::now(); // End time
    auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
    std::cout << "k2-raster temporal build time: " << time;
    std::cout << " milliseconds." << std::endl;
    double k2_raster_size = sdsl::size_in_mega_bytes(k2raster);
    std::cout << "k2-raster space: " << std::setprecision(2) << std::fixed << k2_raster_size << " Mbs" << std::endl;

#ifndef NDEBUG
    std::cout << "Space by time: " << std::endl;
    k2raster.print_space_by_time();
#endif

    /*********************/
    /* Save structure    */
    /*********************/
#ifndef NDEBUG
    std::cout << std::endl << "Storing k2-raster temporal structure in file: " << output_filename  << std::endl;
#endif
    sdsl::store_to_file(k2raster, output_filename);
#ifndef NDEBUG
    std::string file_name = std::string(output_filename) + ".html";
    sdsl::write_structure<sdsl::format_type::HTML_FORMAT>(k2raster, file_name);
#endif



    //************************//
    // TEST structure         //
    //************************//
    if (set_check) {
        std::cout << std::endl << "Checking k2-raster structure....... (in memory)" << std::endl;

        // Load values
        if (k2raster.check(inputs_filename, input_path, scale_factor)) {
            std::cout << "Test Values: OK!!" << std::endl;
        } else {
            std::cout << "Test Values: FAILED!!" << std::endl;
        }

        std::cout << std::endl << "Checking k2-raster structure....... (loaded from disk)" << std::endl;

        // Load structure
        std::ifstream input_file(output_filename);
        assert(input_file.is_open() && input_file.good());
        k2_raster_temporal_type k2_raster_test;
        k2_raster_test.load(input_file);
        input_file.close();

        // Load values
        if (k2_raster_test.check(inputs_filename, input_path, scale_factor)) {
            std::cout << "Test Values: OK!!" << std::endl;
        } else {
            std::cout << "Test Values: FAILED!!" << std::endl;
        }
    }
}

int main(int argc, char **argv) {

    if (argc != 8) {
        print_help(argv[0]);
        exit(-1);
    }

    /*********************/
    /* Reads params      */
    /*********************/
    std::string inputs_filename = argv[1];
    std::string input_path = argv[2];
    ushort scale_factor = atoi(argv[3]);
    std::string output_filename = argv[4];
    ushort k2_raster_type = atoi(argv[5]);
    size_t snap_freq = atol(argv[6]);
    bool set_check = atoi(argv[7]);

    /*********************/
    /* Encodes data      */
    /*********************/
    switch (k2_raster_type) {
        case k2raster::K2_RASTER_TEMPORAL_TYPE:
            run_encode<k2raster::k2_raster_temporal<>>(inputs_filename, input_path, scale_factor, output_filename, snap_freq, set_check);
            break;
        case k2raster::K2_RASTER_TEMPORAL_TYPE_H:
            run_encode<k2raster::k2rth_type>(inputs_filename, input_path, scale_factor, output_filename, snap_freq, set_check);
            break;
        case k2raster::T_K2_RASTER_TYPE:
            run_encode<k2raster::t_k2_raster<>>(inputs_filename, input_path, scale_factor, output_filename, snap_freq, set_check);
            break;
        case k2raster::AT_K2_RASTER_TYPE:
            run_encode<k2raster::atk2r_type>(inputs_filename, input_path, scale_factor, output_filename, snap_freq, set_check);
            break;
        case k2raster::ATH_K2_RASTER_TYPE:
            run_encode<k2raster::athk2r_type>(inputs_filename, input_path, scale_factor, output_filename, snap_freq, set_check);
            break;
        case k2raster::K2_RASTER_TEMPORAL_GLOBAL_TYPE:
            run_encode<k2raster::k2_raster_temporal_global<>>(inputs_filename, input_path, scale_factor, output_filename, snap_freq, set_check);
            break;
        default:
            print_help(argv[0]);
            std::cout << "Invalid type " << k2_raster_type << ": " << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TEMPORAL_TYPE << ": k2-raster temporal (set of k2-raster)." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TEMPORAL_TYPE_H << ": k2-raster temporal (set of k2-raster_heuristic)." << std::endl;
            std::cout << "\t Type " << k2raster::T_K2_RASTER_TYPE << ": t_k2-raster: Temporal k2-raster." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TEMPORAL_GLOBAL_TYPE << ": k2-rasters with an unique dictionary." << std::endl;
            exit(-1);
    }
}
