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


// Own libraries
#include <k2_raster.hpp>
#include <k2_raster_heuristic.hpp>
#include <utils/query/query.hpp>
#include <utils/utils_time.hpp>
#include <utils/args/utils_args_raster.hpp>

template<typename k2_raster_type>
void run_queries(std::string k2raster_filename, std::string &query_filename, uint n_reps=1) {

    /*********************/
    /* Read queries      */
    /*********************/
    std::vector<k2raster::query<false>> queries;
    std::ifstream query_file(query_filename);
    assert(query_file.is_open() && query_file.good());

    // Read queries
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
    for (uint r = 0; r < n_reps; r++) {
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
    } // END FOR n_reps
    auto t2 = util::time::user::now(); // End time

    { // Print Info
        auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
        std::cout << "Time: " << time << " milliseconds. ";
        std::cout << "Queries: " << n_reps * queries.size() << " (" << n_reps << "x" << queries.size() << ") ";
        std::cout << "Total value = " << total_value << ", us/query = " << ((time * 1000.0)/(double)(n_reps*queries.size())) << std::endl;
    }
}

int main(int argc, char **argv) {

    /*********************/
    /* Reads params      */
    /*********************/
    args_get_cell args;
    parse_args_get_cell(argc, argv, args);

    /*********************/
    /*k2-raster Type     */
    /*********************/
    ushort k2_raster_type = k2raster::get_type(args.input_file);

    /*********************/
    /* Run queries       */
    /*********************/
    switch (k2_raster_type) {
        case k2raster::K2_RASTER_TYPE:
            run_queries<k2raster::k2_raster<>>(args.input_file, args.query_file, args.n_reps);
            break;
        case k2raster::K2_RASTER_TYPE_PLAIN:
            run_queries<k2raster::k2_raster_plain<>>(args.input_file, args.query_file, args.n_reps);
            break;
        case k2raster::K2_RASTER_TYPE_HEURISTIC:
            run_queries<k2raster::k2_raster_heuristic<>>(args.input_file, args.query_file, args.n_reps);
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

