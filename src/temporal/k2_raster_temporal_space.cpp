/*  
 * Created by Fernando Silva on 21/01/19.
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

void print_help(char * argv0) {
    printf("Usage: %s <k2raster_filename>\n",
           argv0);
}

template<typename k2_raster_temporal_type>
void run_space(std::string k2raster_filename) {


    /*********************/
    /* Load structure    */
    /*********************/
    k2_raster_temporal_type k2_raster;
    sdsl::load_from_file(k2_raster, k2raster_filename);

    /*********************/
    /* Print space       */
    /*********************/
    double k2_raster_size = sdsl::size_in_mega_bytes(k2_raster);
    std::cout << "k2-raster space: " << std::setprecision(2) << std::fixed << k2_raster_size << " Mbs" << std::endl;
    std::cout << "Space by time: " << std::endl;
    k2_raster.print_space_by_time();
}

int main(int argc, char **argv) {

    if (argc != 2) {
        print_help(argv[0]);
        exit(-1);
    }

    /*********************/
    /* Reads params      */
    /*********************/
    std::string k2raster_filename = argv[1];

    /*********************/
    /*k2-raster Type     */
    /*********************/
    // Load structure
    std::ifstream k2raster_file(k2raster_filename);
    assert(k2raster_file.is_open() && k2raster_file.good());

    ushort k2_raster_type;
    sdsl::read_member(k2_raster_type, k2raster_file);
    k2raster_file.close();


    /*********************/
    /* Run queries       */
    /*********************/
    switch (k2_raster_type) {
        case k2raster::K2_RASTER_TEMPORAL_TYPE:
            run_space<k2raster::k2_raster_temporal<>>(k2raster_filename);
            break;
        case k2raster::T_K2_RASTER_TYPE:
            run_space<k2raster::t_k2_raster<>>(k2raster_filename);
            break;
        default:
            print_help(argv[0]);
            std::cout << "Invalid k2-raster type " << k2_raster_type << ": " << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TEMPORAL_TYPE << ": k2-raster temporal (set of k2-raster)." << std::endl;
            std::cout << "\t Type " << k2raster::T_K2_RASTER_TYPE << ": t_k2-raster: Temporal k2-raster." << std::endl;
            exit(-1);
    }

    /*********************/
    /* Run space         */
    /*********************/
    
}
