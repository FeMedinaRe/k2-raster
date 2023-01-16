/*  
 * Created by Fernando Silva on 21/12/22.
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


// Own libraries
#include <k2_raster.hpp>
#include <k2_raster_heuristic.hpp>
#include <utils/query/query.hpp>
#include <utils/utils_time.hpp>
#include <utils/args/utils_args_raster.hpp>

template<typename k2_raster_type>
void print_info(std::string k2raster_filename) {


    /*********************/
    /* Load structure    */
    /*********************/
    k2_raster_type k2_raster;
    sdsl::load_from_file(k2_raster, k2raster_filename);

    /*********************/
    /* Print info        */
    /*********************/
    k2_raster.print_info_raster();
}

int main(int argc, char **argv) {

    /*********************/
    /* Reads params      */
    /*********************/
    args_show_info args;
    parse_args_show_info (argc, argv, args);

    /*********************/
    /*k2-raster Type     */
    /*********************/
    ushort k2_raster_type = k2raster::get_type(args.input_file);

    /*********************/
    /* Run queries       */
    /*********************/
    switch (k2_raster_type) {
        case k2raster::K2_RASTER_TYPE:
            print_info<k2raster::k2_raster<>>(args.input_file);
            break;
        case k2raster::K2_RASTER_TYPE_PLAIN:
            print_info<k2raster::k2_raster_plain<>>(args.input_file);
            break;
        case k2raster::K2_RASTER_TYPE_HEURISTIC:
            print_info<k2raster::k2_raster_heuristic<>>(args.input_file);
            break;
        default:
            print_usage_get_cell(argv);
            std::cout << "Invalid k2-raster type " << k2_raster_type << ": " << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE << ": hybrid k2-raster." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE_PLAIN << ": plain k2-raster" << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE_HEURISTIC << ": heuristic k2-raster ." << std::endl;
            exit(-1);
    }
}

