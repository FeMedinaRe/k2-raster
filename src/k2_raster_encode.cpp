/*
 * Created by Fernando Silva on 29/10/18.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Program to execute k2-raster over a dataset (sequence of integers of 32bits)
 *
 * k2-raster is a compact data structure to represent raster data that
 * uses compressed space and offers indexing capabilities.
 * It uses min/max values for indexing and improving query performance.
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
#include <utils/utils_time.hpp>

void print_help(char * argv0) {
    printf("Usage: %s <input_data> <rows> <cols> <output_data> <set_check> <type> <k1> [k2 level_k1 [plain_levels]] \n",
           argv0);
}

template<typename k2_raster_type>
void run_encode(std::vector<int> values, size_t n_rows, size_t n_cols,
                ushort k1, ushort k2, ushort level_k1, ushort plains_levels,
                std::string output_data, bool set_check) {

    /*********************/
    /* Encodes data      */
    /*********************/
#ifndef NDEBUG
    std::cout << std::endl << "Creating k2-raster structure for matrix " << n_rows << "x" << n_cols << std::endl;
#endif
    auto t1 = util::time::user::now(); // Start time
    k2_raster_type k2raster(values, n_rows, n_cols, k1, k2, level_k1, plains_levels);


    // Print time/space
    auto t2 = util::time::user::now(); // End time
    auto time = util::time::duration_cast<util::time::milliseconds>(t2-t1);
    std::cout << "k2-raster build time: " << time;
    std::cout << " milliseconds." << std::endl;

    size_t k2_raster_size = sdsl::size_in_bytes(k2raster);
    double ratio = (k2_raster_size * 100.) / (n_rows * n_cols * sizeof(int));
    std::cout << "k2-rater space:" << k2_raster_size << " bytes (" << ratio << "%)" << std::endl;

    /*********************/
    /* Save structure    */
    /*********************/
#ifndef NDEBUG
    std::cout << std::endl << "Storing k2-raster structure in file: " << output_data  << std::endl;
#endif
    sdsl::store_to_file(k2raster, output_data);
#ifndef NDEBUG
    std::string file_name = std::string(output_data) + ".html";
    sdsl::write_structure<sdsl::format_type::HTML_FORMAT>(k2raster, file_name);
#endif

    //************************//
    // TEST structure         //
    //************************//
    if (set_check) {
        std::cout << std::endl << "Checking k2-raster structure......." << std::endl;

        // Load structure
        std::ifstream input_file(output_data);
        assert(input_file.is_open() && input_file.good());
        k2_raster_type k2_raster_test;
        k2_raster_test.load(input_file);
        input_file.close();

        // Check values
        if (k2_raster_test.check(values, n_rows, n_cols)) {
            std::cout << "Test Values: OK!!" << std::endl;
        } else {
            std::cout << "Test Values: FAILED!!" << std::endl;
        }
    }
}



int main(int argc, char **argv) {

    if (argc != 8 && argc != 11) {
        print_help(argv[0]);
        exit(-1);
    }

    /*********************/
    /* Reads params      */
    /*********************/
    std::string values_filename = argv[1];
    size_t rows = atoi(argv[2]);
    size_t cols = atoi(argv[3]);
    std::string output_data = argv[4];
    bool set_check = atoi(argv[5]);
    ushort k2_raster_type = atoi(argv[6]);
    ushort k1 = atoi(argv[7]);
    ushort k2 = (argc >= 10 ? atoi(argv[8]) : k1);
    ushort level_k1 = (argc >= 10 ? atoi(argv[9]) : 0);
    ushort plain_levels = (argc == 11 ? atoi(argv[10]) : 0);


    /*********************/
    /* Reads input data  */
    /*********************/
    std::ifstream values_file(values_filename); // Open file
    assert(values_file.is_open() && values_file.good());


    std::vector<int> values(rows * cols);
    size_t n = 0;
    for (size_t r = 0; r < rows; r++) {
        for (size_t c = 0; c < cols; c++) {
            sdsl::read_member(values[n], values_file);
            n++;
        }
    }

    /*********************/
    /* Encodes data      */
    /*********************/
    switch (k2_raster_type) {
        case k2raster::K2_RASTER_TYPE:
            run_encode<k2raster::k2_raster<>>(values, rows, cols, k1, k2, level_k1, plain_levels, output_data, set_check);
            break;
        case k2raster::K2_RASTER_TYPE_PLAIN:
            run_encode<k2raster::k2_raster_plain<>>(values, rows, cols, k1, k2, level_k1, plain_levels, output_data, set_check);
            break;
        case k2raster::K2_RASTER_TYPE_HEURISTIC:
            run_encode<k2raster::k2_raster_heuristic<>>(values, rows, cols, k1, k2, level_k1, plain_levels, output_data, set_check);
            break;
        default:
            print_help(argv[0]);
            std::cout << "Invalid type " << k2_raster_type << ": " << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE_PLAIN << ": hybrid k2-raster." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE_PLAIN << ": k2-raster with plain values." << std::endl;
            std::cout << "\t Type " << k2raster::K2_RASTER_TYPE_HEURISTIC << ": heuristic k2-raster ." << std::endl;
            exit(-1);
    }

}