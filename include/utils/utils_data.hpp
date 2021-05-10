/*  
 * Created by Fernando Silva on 18/07/19.
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

#ifndef INCLUDE_UTILS_UTILS_DATA
#define INCLUDE_UTILS_UTILS_DATA

//! Namespace for k2-raster library
#include <vector>
#include <cstdio>
#include <random>

namespace k2raster {

    //*******************************************************//
    //***************** CREATE VALUES ***********************//
    //*******************************************************//
    template <typename value_type=int>
    std::vector<value_type> create_random_raster(size_t size_x, size_t size_y, value_type min_value=0, value_type max_value = 100) {
        // Random number generator
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<value_type> dis_values(min_value, max_value);

        std::vector<value_type> values;
        values.reserve(size_x * size_y);

        for (size_t r = 0; r < size_x; r++) {
            for (size_t c = 0; c < size_y; c++) {
                values.push_back(dis_values(gen));
            }
        }
        return values;
    }

    //*******************************************************//
    //****************** READ VALUES ************************//
    //*******************************************************//
    bool read_params_from_file(std::ifstream &inputs_file, size_t &n_rows, size_t &n_cols,
                               size_t &k1, size_t &k2, size_t &level_k1, size_t &plain_levels) {

        if (inputs_file.eof()){
            return false;
        }

        inputs_file >> n_rows;
        inputs_file >> n_cols;
        inputs_file >> k1;
        inputs_file >> k2;
        inputs_file >> level_k1;
        inputs_file >> plain_levels;

#ifndef NDEBUG
        std::cout << "Input params-> " << " Size: " << n_rows << "x" << n_cols << " ";
        std::cout << "k1: " << k1 << " k2: " << k2 << " Level_k1: " << level_k1 << " Plain_levels: "
                  << plain_levels << std::endl;
#endif
        return true;
    }

    bool next_file_path(std::ifstream &inputs_file, std::string &input_path_folder, std::string &file_path_input) {
        if (inputs_file.eof()){
            return false;
        }
        std::string filename;
        inputs_file >> filename;
        if (filename.empty()) {
            // Empty line, read next line
            return next_file_path(inputs_file, input_path_folder, file_path_input);
        }
        file_path_input = input_path_folder + "/" + filename;
#ifndef NDEBUG
        std::cout << "File Path -> " << file_path_input << std::endl;
#endif
        return true;
    }

    template <typename value_type=int, typename size_type=size_t>
    size_type read_input_data_BIN(std::string &file_path_input, const size_t n_rows, const size_t n_cols, std::vector<value_type> &values) {
        /*********************/
        /* Reads input data  */
        /*********************/
        std::ifstream input(file_path_input);
        assert(input.is_open() && input.good());
        size_type n = 0;
        values.resize(n_rows * n_cols);
        for (size_t r = 0; r < n_rows; r++) {
            for (size_t c = 0; c < n_cols; c++) {
                sdsl::read_member(values[n], input);
                n++;
            }
        }
        input.close();
        return n;
    }

    template <typename value_type=int, typename size_type=size_t>
    size_type read_values_ASC(std::string &file_path_input, size_t &n_rows, size_t &n_cols, std::vector<value_type> &values, ushort scale_factor) {

        std::ifstream input(file_path_input);
        assert(input.is_open() && input.good());

        /*********************/
        /* READ HEADER       */
        /*********************/
        std::string key;
        double val;
        input >> key >> n_cols;
        input >> key >> n_rows;
        input >> key >> val;
        input >> key >> val;
        input >> key >> val;
        input >> key >> val;

        /*********************/
        /* READ VALUES       */
        /*********************/
        size_type n_values = 0;
        uint factor = std::pow(10, scale_factor);

        // Read values
        while (!input.eof() && input.good()) {
            input >> val;
            values.push_back(val * factor);
            n_values++;
        } // END WHILE read values

        input.close();
        return n_values;
    }

    template <typename value_type=int, typename size_type=size_t>
    size_type read_input_data(std::string &file_path_input, size_t &n_rows, size_t &n_cols, std::vector<value_type> &values, ushort scale_factor) {
        if (file_path_input.length() > 4 && file_path_input.substr(file_path_input.length() - 4) == ".asc") {
            return read_values_ASC<value_type, size_type>(file_path_input, n_rows, n_cols, values, scale_factor);
        } else {
            values.resize(n_rows * n_cols);
            return read_input_data_BIN<value_type, size_type>(file_path_input, n_rows, n_cols, values);
        }
    }
}

#endif //  INCLUDE_UTILS_UTILS_DATA
