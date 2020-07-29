/*  
 * Created by Fernando Silva on 5/02/19.
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
#include <utils/utils_time.hpp>

void print_help(char * argv0) {
    printf("Usage: %s <input_data> <rows> <cols> <k2_raster_file>\n",
           argv0);
}

template<typename k2_raster_type>
void run_check(std::vector<int> values, size_t n_rows, size_t n_cols, std::string k2_raster_file) {

    /*********************/
    /* Load structure    */
    /*********************/
    k2_raster_type k2_raster;
    sdsl::load_from_file(k2_raster, k2_raster_file);


    //************************//
    // TEST structure         //
    //************************//
    auto t1 = util::time::user::now(); // Start time
    if (k2_raster.check(values, n_rows, n_cols)) {
        std::cout << "Test Values: OK!!" << std::endl;
    } else {
        std::cout << "Test Values: FAILED!!" << std::endl;
    }
    auto t2 = util::time::user::now(); // End time

    { // Print Info
        auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
        std::cout << "Time: " << time << " milliseconds. ";
    }
}

int main(int argc, char **argv) {

    if (argc != 5) {
        print_help(argv[0]);
        exit(-1);
    }

    std::string values_filename = argv[1];
    size_t n_rows = atol(argv[2]);
    size_t n_cols = atol(argv[3]);
    std::string k2raster_filename = argv[4];


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
    /* Reads input data  */
    /*********************/
    std::ifstream values_file(values_filename); // Open file
    assert(values_file.is_open() && values_file.good());


    std::vector<int> values(n_rows * n_cols);
    size_t n = 0;
    for (size_t r = 0; r < n_rows; r++) {
        for (size_t c = 0; c < n_cols; c++) {
            sdsl::read_member(values[n], values_file);
            n++;
        }
    }

    /*********************/
    /* Check structure   */
    /*********************/
    switch (k2_raster_type) {
        case k2raster::K2_RASTER_TYPE:
            run_check<k2raster::k2_raster<>>(values, n_rows, n_cols, k2raster_filename);
            break;
        case k2raster::K2_RASTER_TYPE_PLAIN:
            run_check<k2raster::k2_raster_plain<>>(values, n_rows, n_cols, k2raster_filename);
            break;
        case k2raster::K2_RASTER_TYPE_HEURISTIC:
            run_check<k2raster::k2_raster_heuristic<>>(values, n_rows, n_cols, k2raster_filename);
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