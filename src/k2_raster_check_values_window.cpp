/*  
 * Created by Fernando Silva on 5/04/19.
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
#include <k2_raster_heuristic.hpp>
#include <utils/query/query.hpp>
#include <utils/utils_time.hpp>
#include <utils/args/utils_args_raster.hpp>

template<typename k2_raster_type>
void run_queries(std::string k2raster_filename, const std::string &query_filename, bool strong_check, bool set_check, uint n_reps=1) {

    /*********************/
    /* Read queries      */
    /*********************/
    std::vector<k2raster::query<false>> queries;

    std::ifstream query_file(query_filename);
    assert(query_file.is_open() && query_file.good());

    while (query_file.good()) {
        queries.emplace_back(query_file, false);
    }

    /*********************/
    /* Load structure    */
    /*********************/
    k2_raster_type k2_raster;
    sdsl::load_from_file(k2_raster, k2raster_filename);

    /*********************/
    /* Run queries       */
    /*********************/

    auto t1 = util::time::user::now(); // Start time
    size_t n_queries_true = 0;

    for (uint r = 0; r < n_reps; r++) {
#ifndef NDEBUG
        size_t q = 0;
#endif
        for (const auto &query : queries) {
            bool result = k2_raster.check_values_window(query.xini, query.xend, query.yini, query.yend,
                    query.valini, query.valend, strong_check);
            if (result) {
                n_queries_true++;
            }
#ifndef NDEBUG
            std::cout << "Query " << q << " gets " << result << " cells: ";
            query.print();
            q++;
#endif
            if (set_check) {
                int value;
                // Check region
                size_t query_OK = strong_check;
                for (auto x = query.xini; x <= query.xend; x++) {
                    for (auto y = query.yini; y <= query.yend; y++) {
                        value = k2_raster.get_cell(x, y);
                        if (value >= query.valini &&  value <= query.valend) {
                            if (!strong_check && !result) {
                                std::cout << "[ERROR] Cell (" << x << ", " << y << ") ";
                                std::cout << " expected " << value << std::endl;
                                exit(-1);
                            }
                            if (!strong_check) {
                                query_OK = true;
                                break;
                            }
                        } else {
                            if (strong_check) {
                                query_OK = false;
                                break;
                            }
                        }
                    } // END FOR y
                    if ((strong_check && !query_OK) || (!strong_check && query_OK)) {
                        break;
                    }
                } // END FOR x
                if ((result && !query_OK) || (!result && query_OK)) {
                    std::cout << "[ERROR] Query returns " << result << ", expected " << query_OK << std::endl;
                    exit(-1);
                }
            } // END IF set_check
        } // END FOR queries
    } // END FOR NREPS
    auto t2 = util::time::user::now(); // End time

    { // Print Info
        auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
        std::cout << "Time: " << time << " milliseconds. ";
        std::cout << "Queries: " << n_reps * queries.size() << " (" << n_reps << "x" << queries.size() << ") ";
        std::cout << "Result: true = " << n_queries_true << ", false = " << n_reps * queries.size() - n_queries_true << ", ";
        std::cout << "us/query = " << ((time * 1000.0)/(double)(n_reps*queries.size())) << std::endl;
    }
}

int main(int argc, char **argv) {

    /*********************/
    /* Reads params      */
    /*********************/
    args_check_values args;
    parse_args_check_values(argc, argv, args);

    /*********************/
    /*k2-raster Type     */
    /*********************/
    ushort k2_raster_type = k2raster::get_type(args.input_file);

    /*********************/
    /* Run queries       */
    /*********************/
    switch (k2_raster_type) {
        case k2raster::K2_RASTER_TYPE:
            run_queries<k2raster::k2_raster<>>(args.input_file, args.query_file, args.strong_check, args.set_check, args.n_reps);
            break;
        case k2raster::K2_RASTER_TYPE_PLAIN:
            run_queries<k2raster::k2_raster_plain<>>(args.input_file, args.query_file, args.strong_check, args.set_check, args.n_reps);
            break;
        case k2raster::K2_RASTER_TYPE_HEURISTIC:
            run_queries<k2raster::k2_raster_heuristic<>>(args.input_file, args.query_file, args.strong_check, args.set_check, args.n_reps);
            break;
        default:
            print_usage_check_values(argv);
            std::cout << "Invalid k2-raster type " << k2_raster_type << ": " << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE << ": hybrid k2-raster." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE_PLAIN << ": k2-raster with plain values." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE_HEURISTIC << ": heuristic k2-raster ." << std::endl;
            exit(-1);
    }
}

