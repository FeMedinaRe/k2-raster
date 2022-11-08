/*  
 * Created by Fernando Silva on 13/06/22.
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


// System libraries

// Own libraries
#include <k2_raster_heuristic.hpp>
#include <utils/args/utils_args_utils.hpp>

// Third libraries


template<typename k2_raster_type>
void extract_values(const std::string& k2raster_filename, const std::string& output_filename) {


    /*********************/
    /* Load structure    */
    /*********************/
    k2_raster_type k2_raster;
    sdsl::load_from_file(k2_raster, k2raster_filename);

    /**********************/
    /* Create output file */
    /**********************/
    std::ofstream last_level_file(output_filename);
    assert(last_level_file.is_open() && last_level_file.good());

    auto n_values = k2_raster.store_last_level(last_level_file);
    last_level_file.close();

    std::cout << "Extracted " << n_values << " and stored in " << output_filename << std::endl;
}

int main(int argc, char **argv) {

    /*********************/
    /* Reads params      */
    /*********************/
    args_extract_last args;
    parse_args_extract_last(argc, argv, args);

    /*********************/
    /*k2-raster Type     */
    /*********************/
    // Load structure
    std::ifstream k2raster_file(args.raster);
    assert(k2raster_file.is_open() && k2raster_file.good());

    ushort k2_raster_type;
    sdsl::read_member(k2_raster_type, k2raster_file);
    k2raster_file.close();

    /*********************/
    /* Extract values    */
    /*********************/
    switch (k2_raster_type) {
        case k2raster::K2_RASTER_TYPE:
            extract_values<k2raster::k2_raster<>>(args.raster, args.output_data);
            break;
        case k2raster::K2_RASTER_TYPE_PLAIN:
            extract_values<k2raster::k2_raster_plain<>>(args.raster, args.output_data);
            break;
        case k2raster::K2_RASTER_TYPE_HEURISTIC:
            extract_values<k2raster::k2_raster_heuristic<>>(args.raster, args.output_data);
            break;
        default:
            print_usage_extract_last(argv);
            std::cout << "Invalid k2-raster type " << k2_raster_type << ": " << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE << ": hybrid k2-raster." << std::endl;
            exit(-1);
    }
}

