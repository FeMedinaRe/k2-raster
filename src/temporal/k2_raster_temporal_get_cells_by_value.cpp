/*  
 * Created by Fernando Silva on 16/01/19.
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
#include <temporal/k2_raster_temporal_global.hpp>
#include <temporal/t_k2_raster.hpp>
#include <temporal/at_k2_raster.hpp>
#include <utils/query/query.hpp>
#include <utils/utils_time.hpp>

#define N_REPS 1

void print_help(char * argv0) {
    printf("Usage: %s <k2_raster_file> <query_file> <check> [nreps]\n", argv0);
}


template<typename k2_raster_temporal_type, bool n_queries_beginning=true>
void run_queries(std::string k2raster_filename, const std::string &query_filename, bool set_check, uint nreps) {

    /*********************/
    /* Read queries      */
    /*********************/
    std::vector<k2raster::query<>> queries;

    std::ifstream query_file(query_filename);
    assert(query_file.is_open() && query_file.good());

    size_t n_queries = 0;
    if (n_queries_beginning) {
        query_file >> n_queries;
    }

    while (query_file.good() && (!n_queries_beginning || queries.size() < n_queries)) {
        queries.emplace_back(query_file, false);
    }
    query_file.close();

    /*********************/
    /* Load structure    */
    /*********************/
    k2_raster_temporal_type k2_raster;
    sdsl::load_from_file(k2_raster, k2raster_filename);

    /*********************/
    /* Run queries       */
    /*********************/

    auto t1 = util::time::user::now(); // Start time
    size_t total_num_cells;
    for (uint r = 0; r < nreps; r++) {
        total_num_cells = 0;
#ifndef NDEBUG
        size_t q = 0;
#endif
        for (const auto &query : queries) {

            std::vector<std::vector<std::pair<size_t, size_t>>> result(query.tend - query.tini + 1);
            size_t n_cells = k2_raster.get_cells_by_value(query.xini, query.xend, query.yini, query.yend,
                                                          query.valini, query.valend, query.tini, query.tend, result);
            total_num_cells += n_cells;
#ifndef NDEBUG
            std::cout << "Query " << q << " gets " << n_cells << " cells";
            query.print();
            q++;
#endif

            if (set_check) {
                int value;
                size_t time = 0;

                // Check result
                for (auto t = query.tini; t < query.tend; t++) {
                    for (const auto &cell : result[t - query.tini]) {
                        value = k2_raster.get_cell(cell.first, cell.second, t);
                        if (value < query.valini || value > query.valend) {
                            std::cout << "[ERROR] Cell (" << cell.first << ", " << cell.second << ")";
                            std::cout << " is not within range [ " << query.valini << " - " << query.valend << "]" << std::endl;
                            exit(-1);
                        }
                    } // END FOR cell
                    time++;
                } // END FOR time

                // Check region
                size_t cells_region = 0;
                for (auto x = query.xini; x <= query.xend; x++) {
                    for (auto y = query.yini; y <= query.yend; y++) {
                        for (auto t = query.tini; t <= query.tend; t++) {
                            value = k2_raster.get_cell(x, y, t);
                            if (value >= query.valini && value <= query.valend) {
                                cells_region++;
#ifndef NDEBUG
                                // Search cell in the result
                                bool found = false;
                                for (const auto &cell : result[t - query.tini]) {
                                    if (cell.first == x && cell.second == y) {
                                        found = true;
                                        break;
                                    }
                                } // END FOR cell

                                if (!found) {
                                    std::cout << "[ERROR] Cell (" << x << ", " << y << ") at time " << t << " not found" << std::endl;
                                    exit(-1);
                                }
#endif
                            }
                        } // END FOR t
                    } // END FOR y
                } // END FOR x

                if (cells_region != n_cells) {
                    std::cout << "[ERROR] Query returns " << n_cells << " cells and expected " << cells_region << " cells" << std::endl;
                    exit(-1);
                }

            } // END IF set_check
        }
    }
    auto t2 = util::time::user::now(); // End time

    { // Print Info
        auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
        std::cout << "Time: " << time << " milliseconds. ";
        std::cout << "Queries: " << nreps * queries.size() << " (" << nreps << "x" << queries.size() << ") ";
        std::cout << "Nº cells = " << total_num_cells << ", us/query = " << ((time * 1000.0)/(nreps*queries.size()));
        std::cout << ", us/cell = " << ((time * 1000.0))/total_num_cells << std::endl;
    }
}

int main(int argc, char **argv) {

    if (argc != 4 && argc != 5) {
        print_help(argv[0]);
        exit(-1);
    }

    /*********************/
    /* Reads params      */
    /*********************/
    std::string k2raster_filename = argv[1];
    std::string query_filename = argv[2];
    bool set_check = atoi(argv[3]);
    uint nreps = argc == 5 ? atoi(argv[4]) : N_REPS;

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
            run_queries<k2raster::k2_raster_temporal<>>(k2raster_filename, query_filename, set_check, nreps);
            break;
        case k2raster::K2_RASTER_TEMPORAL_TYPE_H:
            run_queries<k2raster::k2_raster_temporal<int, k2raster::k2_raster_heuristic<>>>(k2raster_filename, query_filename, set_check, nreps);
            break;
        case k2raster::T_K2_RASTER_TYPE:
            run_queries<k2raster::t_k2_raster<>>(k2raster_filename, query_filename, set_check, nreps);
            break;
        case k2raster::AT_K2_RASTER_TYPE:
            run_queries<k2raster::atk2r_type>(k2raster_filename, query_filename, set_check, nreps);
            break;
        case k2raster::ATH_K2_RASTER_TYPE:
            run_queries<k2raster::athk2r_type>(k2raster_filename, query_filename, set_check, nreps);
            break;
        case k2raster::K2_RASTER_TEMPORAL_GLOBAL_TYPE:
            run_queries<k2raster::k2_raster_temporal_global<>>(k2raster_filename, query_filename, set_check, nreps);
            break;
        default:
            print_help(argv[0]);
            std::cout << "Invalid type " << k2_raster_type << ": " << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TEMPORAL_TYPE << ": k2-raster temporal (set of k2-raster)." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TEMPORAL_TYPE_H << ": k2-raster temporal (set of k2-raster_heuristic)." << std::endl;
            std::cout << "\t Type " << k2raster::T_K2_RASTER_TYPE << ": t_k2-raster: Temporal k2-raster." << std::endl;
            std::cout << "\t Type " << k2raster::AT_K2_RASTER_TYPE << ": at_k2-raster: adaptive Temporal k2-raster." << std::endl;
            std::cout << "\t Type " << k2raster::ATH_K2_RASTER_TYPE << ": ath_k2-raster: adaptive Temporal Heuristic k2-raster." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TEMPORAL_GLOBAL_TYPE << ": k2-rasters with an unique dictionary." << std::endl;
            exit(-1);
    }
}

