/*  
 * Created by Fernando Silva on 9/01/19.
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

#ifndef INCLUDED_K2_RASTER_TEMPORAL
#define INCLUDED_K2_RASTER_TEMPORAL

#include <k2_raster.hpp>
#include <k2_raster_heuristic.hpp>
#include <temporal/k2_raster_temporal_base.hpp>

//! Namespace for k2-raster library
namespace k2raster {

template <typename t_value=int,
        typename k2_raster_t=k2_raster<t_value>>
class k2_raster_temporal : public k2_raster_temporal_base<t_value> {

public:
    typedef k2_raster_t                            k2_raster_type; // k2_raster type
    typedef t_value                                value_type;
    typedef k2_raster_temporal_base<>::size_type   size_type;

protected:
    std::vector<k2_raster_type> k_rasters;

public:
    //*******************************************************//
    //******************* CONSTRUCTOR ***********************//
    //*******************************************************//
    k2_raster_temporal() = default;

    k2_raster_temporal(const k2_raster_temporal &tr) {
        *this = tr;
    }

    k2_raster_temporal(k2_raster_temporal &&tr) {
        *this = std::move(tr);
    }

    // Only interface (snapshots_freq unused)
    k2_raster_temporal(std::string &inputs_filename, std::string &input_path_folder, size_type snapshots_freq __attribute__((unused)), const ushort scale_factor=0)
    : k2_raster_temporal(inputs_filename, input_path_folder, scale_factor){};

    k2_raster_temporal(std::string &inputs_filename, std::string &input_path_folder, const ushort scale_factor=0) {
        this->k_k2_raster_type = K2_RASTER_TEMPORAL_TYPE;
        this->k_n_times = build(inputs_filename, input_path_folder, scale_factor);

        if (typeid(k2_raster_t) == typeid(k2rh_type)) {
            this->k_k2_raster_type = K2_RASTER_TEMPORAL_TYPE_H;
        }
    }


    //*******************************************************//
    //*************** BASIC OPERATIONS **********************//
    //*******************************************************//

    //! Move assignment operator
    k2_raster_temporal &operator=(k2_raster_temporal &&tr) {
        if (this != &tr) {
            k2_raster_temporal_base<t_value>::operator=(tr);
            k_rasters = std::move(tr.k_rasters);
        }
        return *this;
    }

    //! Assignment operator
    k2_raster_temporal &operator=(const k2_raster_temporal &tr) {
        if (this != &tr) {
            k2_raster_temporal_base<t_value>::operator=(tr);
            k_rasters = tr.k_rasters;
        }
        return *this;
    }

    //! Swap operator
    void swap(k2_raster_temporal &tr) {
        if (this != &tr) {
            k2_raster_temporal_base<t_value>::swap(tr);
            std::swap(k_rasters, tr.k_rasters);
        }
    }

    //! Equal operator
    bool operator==(const k2_raster_temporal &tr) const {
        if (!k2_raster_temporal_base<t_value>::operator==(tr)) {
            return false;
        }
        if (k_rasters.size() != tr.k_rasters.size()) {
            return false;
        }
        size_type c = 0;
        for (auto raster : k_rasters) {
            if (raster != tr.k_rasters[c++]) {
                return false;
            }
        }
        return true;
    }


    //*******************************************************//
    //******************** QUERIES **************************//
    //*******************************************************//
    value_type get_cell(size_type row, size_type col, size_type time) const {
        assert(time <= (k_rasters.size()-1));
        return k_rasters[time].get_cell(row, col);
    }

    size_type get_cells_by_value(size_type xini, size_type xend, size_type yini, size_type yend,
            value_type valini, value_type valend, size_type tmin, size_type tmax,
            std::vector<std::vector<std::pair<size_type, size_type>>> &result) {

        result.resize(tmax - tmin + 1);
        size_type c = 0;
        size_type count_cells = 0;
        for (auto t = tmin; t <= tmax; t++) {
            count_cells += k_rasters[t].get_cells_by_value(xini, xend, yini, yend, valini, valend, result[c++]);
        }
        return count_cells;
    }

    size_type get_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
            size_type tmin, size_type tmax,
            std::vector<std::vector<value_type>> &result) {
        result.resize(tmax - tmin + 1);
        size_type c = 0;
        size_type count_cells = 0;
        for (auto t = tmin; t <= tmax; t++) {
            count_cells += k_rasters[t].get_values_window(xini, xend, yini, yend, result[c++]);
        }
        return count_cells;
    }

    bool check_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                     value_type valini, value_type valend,
                                     size_type tmin, size_type tmax, bool strong_check) {
        bool result;
        for (auto t = tmin; t <= tmax; t++) {
            result = k_rasters[t].check_values_window(xini, xend, yini, yend, valini, valend, strong_check);
            if (strong_check){
                if (!result) return false;
            } else {
                if (result) return true;
            }
        }
        return strong_check;
    }

    //*******************************************************//
    //********************** FILE ***************************//
    //*******************************************************//
    virtual size_type
    serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {

        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        /****** k2-raster temporal base ******/
        written_bytes += k2_raster_temporal_base<t_value>::serialize(out, child, name);

        /****** Vector of rasters ******/
        written_bytes += write_member(k_rasters.size(), out, child, "num_t");
        for (auto l = 0; l < k_rasters.size(); l++) {
            std::string name = "raster_t_" + std::to_string(l);
            written_bytes += k_rasters[l].serialize(out, child, name);
        }

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    virtual void load(std::istream &in) {
        /****** k2-raster base ******/
        k2_raster_temporal_base<t_value>::load(in);

        /****** Vector of rasters ******/
        ulong t;
        sdsl::read_member(t, in);
        k_rasters.resize(t);
        k_rasters.shrink_to_fit();
        for (auto l = 0; l < t; l++) {
            k_rasters[l].load(in);
        }
    }

    void print_space_by_time() const {
        size_type times = k_rasters.size();
        double size;
        double total_size = 0;
        for (auto t = 0; t < times; t++) {
            size = sdsl::size_in_mega_bytes(k_rasters[t]);
            std::cout << "Time " << t << ": " << std::setprecision(2) << std::fixed << size << " Mbs" << std::endl;
            total_size += size;
        }
        std::cout << "Total size times: " << std::setprecision(2) << std::fixed << total_size << " Mbs" << std::endl;
    };

protected:
    //*******************************************************//
    //******************* AUXILIARY *************************//
    //*******************************************************//

    //*****************************//
    //***** BUILD FUNCTIONS   *****//
    //*****************************//
    size_type build(std::string &inputs_filename, std::string &input_path_folder, const ushort scale_factor=0) {

        /**************************/
        /* Reads params from file */
        /**************************/
        std::ifstream inputs_file(inputs_filename);
        assert(inputs_file.is_open() && inputs_file.good());
        size_type time_count = 0;

        this->k_max_size_x = 0;
        this->k_max_size_y = 0;
        this->k_min_value = std::numeric_limits<value_type>::max();
        this->k_max_value = std::numeric_limits<value_type>::min();

        /**************************/
        /* Encode each k2-raster  */
        /**************************/
        {
            size_t n_rows, n_cols, k1, k2, level_k1, plain_levels;
            std::string file_path_input;

            /*********************/
            /* Read parameters   */
            /*********************/
            read_params_from_file(inputs_file, n_rows, n_cols, k1, k2, level_k1, plain_levels);

            /*********************/
            /* Process rasters   */
            /*********************/
            while (next_file_path(inputs_file, input_path_folder, file_path_input)) {

                /*********************/
                /* Encode data       */
                /*********************/
                std::vector<value_type> values;
                read_input_data(file_path_input, n_rows, n_cols, values, scale_factor);

                k_rasters.emplace_back(values, n_rows, n_cols, k1, k2, level_k1, plain_levels);

                this->k_max_size_x = std::max(this->k_max_size_x, k_rasters[time_count].get_n_rows());
                this->k_max_size_y = std::max(this->k_max_size_y, k_rasters[time_count].get_n_cols());
                this->k_min_value = std::min(this->k_min_value, k_rasters[time_count].min_value);
                this->k_max_value = std::max(this->k_max_value, k_rasters[time_count].max_value);
                time_count++;
            } // END WHiLE read_params
            inputs_file.close();
        } // END BLOCK encode k2-rasters
        return time_count;
    }

}; // END CLASS k2_raster_temporal


// Types
typedef k2raster::k2_raster_temporal<int, k2rh_type> k2rth_type;

}; // END NAMESPACE k2raster

#endif // INCLUDED_K2_RASTER_TEMPORAL
