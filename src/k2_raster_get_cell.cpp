/*  
 * Created by Fernando Silva on 4/04/19.
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

#include <k2_raster.hpp>
#include <utils/query/query.hpp>
#include <utils/utils_time.hpp>
#include <k2_raster_heuristic.hpp>

#define NREPS 1

void print_help(char * argv0) {
    printf("Usage: %s <k2_raster_file> <query_file> <check>\n", argv0);
}

template<typename k2_raster_type>
void run_queries(std::string k2raster_filename, std::string query_filename) {

    /*********************/
    /* Read queries      */
    /*********************/
    std::vector<k2raster::query<false>> queries;

    std::ifstream query_file(query_filename);
    assert(query_file.is_open() && query_file.good());

    while (query_file.good()) {

        k2raster::query<false> query(query_file, true);
        if (query_file.good()) {
            queries.push_back(query);
        }
    }

    /*********************/
    /* Load structure    */
    /*********************/
    k2_raster_type k2_raster;
    sdsl::load_from_file(k2_raster, k2raster_filename);

    /*********************/
    /* Run queries       */
    /*********************/
    long total_value;
    auto t1 = util::time::user::now(); // Start time
    for (uint r = 0; r < NREPS; r++) {
        total_value = 0;
#ifndef NDEBUG
        size_t q = 0;
#endif
        for (const auto &query : queries) {
            int value = k2_raster.get_cell(query.xini, query.yini);
            total_value += value;
#ifndef NDEBUG
            std::cout << "Query " << q << " gets " << value << " ";
            query.print();
            q++;
#endif
        } // END FOR queries
    } // END FOR nreps
    auto t2 = util::time::user::now(); // End time

    { // Print Info
        auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
        std::cout << "Time: " << time << " milliseconds. ";
        std::cout << "Queries: " << NREPS * queries.size() << " (" << NREPS << "x" << queries.size() << ") ";
        std::cout << "Total value = " << total_value << ", us/query = " << ((time * 1000.0)/(NREPS*queries.size())) << std::endl;
    }
}

int main(int argc, char **argv) {

    if (argc != 3) {
        print_help(argv[0]);
        exit(-1);
    }

    /*********************/
    /* Reads params      */
    /*********************/
    std::string k2raster_filename = argv[1];
    std::string query_filename = argv[2];


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
        case k2raster::K2_RASTER_TYPE:
            run_queries<k2raster::k2_raster<>>(k2raster_filename, query_filename);
            break;
        case k2raster::K2_RASTER_TYPE_PLAIN:
            run_queries<k2raster::k2_raster_plain<>>(k2raster_filename, query_filename);
            break;
        case k2raster::K2_RASTER_TYPE_HEURISTIC:
            run_queries<k2raster::k2_raster_heuristic<>>(k2raster_filename, query_filename);
            break;
        default:
            print_help(argv[0]);
            std::cout << "Invalid k2-raster type " << k2_raster_type << ": " << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE << ": hybrid k2-raster." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE_PLAIN << ": k2-raster with plain values." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE_HEURISTIC << ": heuristic k2-raster ." << std::endl;
            exit(-1);
    }
}

