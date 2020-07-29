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
#include <utils/query/query.hpp>
#include <utils/utils_time.hpp>
#include <k2_raster_heuristic.hpp>

#define NREPS 1

void print_help(char * argv0) {
    printf("Usage: %s <k2_raster_file> <query_file> <strong_check> <check>\n", argv0);
}


template<typename k2_raster_type>
void run_queries(std::string k2raster_filename, std::string query_filename, bool strong_check, bool set_check) {

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

    for (uint r = 0; r < NREPS; r++) {
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
                            if (!strong_check) {
                                query_OK = true;
                                break;
                            }

                            if (!strong_check && !result) {
                                std::cout << "[ERROR] Cell (" << x << ", " << y << ") ";
                                std::cout << " expected " << value << std::endl;
                                exit(-1);
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
        std::cout << "Queries: " << NREPS * queries.size() << " (" << NREPS << "x" << queries.size() << ") ";
        std::cout << "Result: true = " << n_queries_true << ", false = " << NREPS * queries.size() - n_queries_true << ", ";
        std::cout << "us/query = " << ((time * 1000.0)/(NREPS*queries.size())) << std::endl;
    }
}

int main(int argc, char **argv) {

    if (argc != 5) {
        print_help(argv[0]);
        exit(-1);
    }

    /*********************/
    /* Reads params      */
    /*********************/
    std::string k2raster_filename = argv[1];
    std::string query_filename = argv[2];
    bool strong_check = atoi(argv[3]);
    bool set_check = atoi(argv[3]);

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
            run_queries<k2raster::k2_raster<>>(k2raster_filename, query_filename, strong_check, set_check);
            break;
        case k2raster::K2_RASTER_TYPE_PLAIN:
            run_queries<k2raster::k2_raster_plain<>>(k2raster_filename, query_filename, strong_check, set_check);
            break;
        case k2raster::K2_RASTER_TYPE_HEURISTIC:
            run_queries<k2raster::k2_raster_heuristic<>>(k2raster_filename, query_filename, strong_check, set_check);
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

