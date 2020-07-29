/*  
 * Created by Fernando Silva on 11/06/19.
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

#include <utils/utils_query.hpp>
#include <temporal/k2_raster_temporal.hpp>
#include <temporal/t_k2_raster.hpp>
#include <temporal/k2_raster_temporal_global.hpp>

void print_help(char * argv0) {
    printf("Usage: %s <n_queries> <k2_raster_file> <query_filename>\n",
           argv0);
}

template<typename k2_raster_type>
void run_create_queries(std::string inputs_filename, std::string query_filename, size_t n_queries) {

    /*********************/
    /* Load structure    */
    /*********************/
    k2_raster_type k2_raster;
    sdsl::load_from_file(k2_raster, inputs_filename);

    /*********************/
    /* Get Parameters    */
    /*********************/
    size_t n_rows = k2_raster.get_n_rows();
    size_t n_cols = k2_raster.get_n_cols();
    int valini = k2_raster.get_min_value();
    int valend = k2_raster.get_max_value();
    size_t tend = k2_raster.get_max_t();

    /*********************/
    /* Open query file   */
    /*********************/
    std::ofstream query_file(query_filename);
    assert(query_file.is_open() && query_file.good());

    /*********************/
    /* Create queries    */
    /*********************/
    std::cout << "Creating " << n_queries << " in file " << query_filename << std::endl;
    std::cout << "Nº rows: " << n_rows << " | Nº cols: " << n_cols << " | Times: " << tend << std::endl;
    std::cout << "Values: " << valini << " - " << valend << std::endl;
    k2raster::create_region_queries<size_t, int, true, true>(n_queries, query_file, 0, n_rows -1, 0, n_cols -1, 0, tend-1, valini, valend);
}

int main(int argc, char **argv) {

    if (argc != 4) {
        print_help(argv[0]);
        exit(-1);
    }

    size_t n_queries = atol(argv[1]);
    std::string k2raster_filename = argv[2];
    std::string query_filename = argv[3];

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
    /* Create queries    */
    /*********************/
    switch (k2_raster_type) {
        case k2raster::K2_RASTER_TEMPORAL_TYPE:
            run_create_queries<k2raster::k2_raster_temporal<>>(k2raster_filename, query_filename, n_queries);
            break;
        case k2raster::K2_RASTER_TEMPORAL_TYPE_H:
            run_create_queries<k2raster::k2_raster_temporal<int, k2raster::k2_raster_heuristic<>>>(k2raster_filename, query_filename, n_queries);
            break;
        case k2raster::T_K2_RASTER_TYPE:
            run_create_queries<k2raster::t_k2_raster<>>(k2raster_filename, query_filename, n_queries);
            break;
        case k2raster::K2_RASTER_TEMPORAL_GLOBAL_TYPE:
            run_create_queries<k2raster::k2_raster_temporal_global<>>(k2raster_filename, query_filename, n_queries);
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
