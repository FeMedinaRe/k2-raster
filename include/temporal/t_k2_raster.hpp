/*  
 * Created by Fernando Silva on 10/01/19.
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

#ifndef INCLUDED_T_K2_RASTER
#define INCLUDED_T_K2_RASTER

#include <k2_raster.hpp>
#include <temporal/k2_raster_temporal_base.hpp>
#include <temporal/helpers/t_k2_raster_log.hpp>

//! Namespace for k2-raster library
namespace k2raster {

template <typename t_value=int,
        typename t_bv=sdsl::bit_vector,
        typename t_rank=sdsl::rank_support_v5<1,1>,
        typename t_values_vec=sdsl::dac_vector_dp_opt<t_bv, t_rank, 3>>
class t_k2_raster : public k2_raster_temporal_base<t_value> {

public:
    typedef k2_raster<t_value,t_bv, t_rank, t_values_vec>           k2_raster_type; // k2_raster type
    typedef t_k2_raster_log<t_value, t_bv, t_rank, t_values_vec>    k2_raster_log_type; // k2_raster type
    typedef k2_raster_temporal_base<>::size_type                    size_type;
    typedef t_value                                                 value_type;

protected:
    size_type                       k_snapshots_freq;
    std::vector<k2_raster_type>     k_snapshots;
    std::vector<k2_raster_log_type> k_logs;

public:
    //*******************************************************//
    //******************* CONSTRUCTOR ***********************//
    //*******************************************************//
    t_k2_raster() = default;

    t_k2_raster(const t_k2_raster &tr) {
        *this = tr;
    }

    t_k2_raster(t_k2_raster &&tr) {
        *this = std::move(tr);
    }

    t_k2_raster(std::string &inputs_filename, std::string &input_path_folder, size_type snapshots_freq, const ushort scale_factor=0) {
        k_snapshots_freq = snapshots_freq;
        this->k_k2_raster_type = T_K2_RASTER_TYPE;

        this->k_n_times = build(inputs_filename, input_path_folder, scale_factor);
    }

    //*******************************************************//
    //*************** BASIC OPERATIONS **********************//
    //*******************************************************//

    //! Move assignment operator
    t_k2_raster &operator=(t_k2_raster &&tr) {
        if (this != &tr) {
            k2_raster_temporal_base<t_value>::operator=(tr);
            k_snapshots_freq = tr.k_snapshots_freq;
            k_snapshots = std::move(tr.k_snapshots);
            k_logs = std::move(tr.k_logs);
        }
        return *this;
    }

    //! Assignment operator
    t_k2_raster &operator=(const t_k2_raster &tr) {
        if (this != &tr) {
            k2_raster_temporal_base<t_value>::operator=(tr);
            k_snapshots_freq = tr.k_snapshots_freq;
            k_snapshots = tr.k_snapshots;
            k_logs = tr.k_logs;
        }
        return *this;
    }

    //! Swap operator
    void swap(t_k2_raster &tr) {
        if (this != &tr) {
            k2_raster_temporal_base<t_value>::swap(tr);
            std::swap(k_snapshots_freq, tr.k_snapshots_freq);
            std::swap(k_snapshots, tr.k_snapshots);
            std::swap(k_logs, tr.k_logs);
        }
    }

    //! Equal operator
    bool operator==(const t_k2_raster &tr) const {
        if (!k2_raster_temporal_base<t_value>::operator==(tr)) {
            return false;
        }

        if (k_snapshots_freq != tr.k_snapshots_freq) {
            return false;
        }

        if ((k_snapshots.size() != tr.k_snapshots.size()) ||
                (k_logs.size() != tr.k_logs.size()) ) {
            return false;
        }

        // Snapshots
//        size_type c = 0;
//        for (auto raster : k_snapshots) {
//            if (raster != tr.k_snapshots[c++]) {
//                return false;
//            }
//        }
//
//        // Logs
//        c = 0;
//        for (auto raster : k_logs) {
//            if (raster != tr.k_logs[c++]) {
//                return false;
//            }
//        }
        return true;
    }

    //*******************************************************//
    //******************** QUERIES **************************//
    //*******************************************************//
    value_type get_cell(size_type row, size_type col, size_type time) const {
        assert(time <= (this->k_n_times));
        size_type last_snap = (time / k_snapshots_freq);
        if (time % k_snapshots_freq) {
            // It is a log
            return k_logs[time - last_snap - 1].get_cell(k_snapshots[last_snap], row, col);
        } else {
            // It is a snapshot
            return k_snapshots[last_snap].get_cell(row, col);
        }

    }

    size_type get_cells_by_value(size_type xini, size_type xend, size_type yini, size_type yend,
                                 value_type valini, value_type valend, size_type tmin, size_type tmax,
                                 std::vector<std::vector<std::pair<size_type, size_type>>> &result) {
        assert(tmin >= 0);
        assert(tmax <= (k_snapshots.size() + k_logs.size() -1));
        assert(tmin <= tmax);

        result.resize(tmax - tmin + 1);;
        size_type count_cells = 0;
        for (auto t = tmin; t <= tmax; t++) {
            size_type last_snap = (t / k_snapshots_freq);
            if (t % k_snapshots_freq) {
                // It is a log
                count_cells += k_logs[t - last_snap - 1].get_cells_by_value(k_snapshots[last_snap],
                        xini, xend, yini, yend, valini, valend, result[t - tmin]);
            } else {
                // It is a snapshot
                count_cells += k_snapshots[last_snap].get_cells_by_value(xini, xend, yini, yend, valini, valend, result[t - tmin]);
            }

        }
        return count_cells;
    }

    size_type get_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                size_type tmin, size_type tmax,
                                std::vector<std::vector<value_type>> &result) {
        assert(tmin >= 0);
        assert(tmax <= (k_snapshots.size() + k_logs.size() -1));
        assert(tmin <= tmax);

        result.resize(tmax - tmin + 1);
        size_type count_cells = 0;
        for (auto t = tmin; t <= tmax; t++) {
            size_type last_snap = (t / k_snapshots_freq);
            if (t % k_snapshots_freq) {
                // It is a log
                count_cells += k_logs[t - last_snap - 1].get_values_window(k_snapshots[last_snap],
                                                                            xini, xend, yini, yend, result[t - tmin]);
            } else {
                // It is a snapshot
                count_cells += k_snapshots[last_snap].get_values_window(xini, xend, yini, yend, result[t - tmin]);
            }

        }
        return count_cells;
    }

    bool check_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                             value_type valini, value_type valend,
                             size_type tmin, size_type tmax, bool strong_check) {
        assert(tmin >= 0);
        assert(tmax <= (k_snapshots.size() + k_logs.size() -1));
        assert(tmin <= tmax);

        bool result;
        for (auto t = tmin; t <= tmax; t++) {
            size_type last_snap = (t / k_snapshots_freq);
            if (t % k_snapshots_freq) {
                // It is a log
               result = k_logs[t - last_snap - 1].check_values_window(k_snapshots[last_snap], xini, xend,
                       yini, yend, valini, valend, strong_check);
            } else {
                // It is a snapshot
                result = k_snapshots[last_snap].check_values_window(xini, xend, yini, yend, valini, valend, strong_check);
            }
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


        /****** t k2-raster ******/
        written_bytes += write_member(k_snapshots_freq, out, child, "snapshots_freq");

        /****** Vector of snapshots ******/
        written_bytes += write_member(k_snapshots.size(), out, child, "num_snapshots");
        written_bytes +=sdsl::serialize_vector(k_snapshots, out, child, "snapshots");

        /****** Vector of logs ******/
        written_bytes += write_member(k_logs.size(), out, child, "num_logs");
        written_bytes +=sdsl::serialize_vector(k_logs, out, child, "logs");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    virtual void load(std::istream &in) {
        /****** k2-raster base ******/
        k2_raster_temporal_base<t_value>::load(in);

        /****** t k2-raster ******/
        sdsl::read_member(k_snapshots_freq, in);

        /****** Vector of snapshots ******/
        size_t n;
        sdsl::read_member(n, in);
        k_snapshots.resize(n);
        sdsl::load_vector(k_snapshots, in);

        /****** Vector of logs ******/
        sdsl::read_member(n, in);
        k_logs.resize(n);
        sdsl::load_vector(k_logs, in);
    }

    void print_space_by_time() const {
        size_type times = k_snapshots.size() + k_logs.size();
        double size;
        double total_size = 0;
        std::string type;

        for (auto t = 0; t < times; t++) {
            size_type last_snap = (t / k_snapshots_freq);
            if (t % k_snapshots_freq) {
                // It is a log
                size = sdsl::size_in_mega_bytes(k_logs[t - last_snap - 1]);
                type = "log";
            } else {
                // It is a snapshot
                size = sdsl::size_in_mega_bytes(k_snapshots[last_snap]);
                type = "snapshot";
            }
            std::cout << "Time " << t << "(" << type << ") : " << std::setprecision(2) << std::fixed << size << " Mbs" << std::endl;
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

        this->k_logs.clear();
        this->k_snapshots.clear();

        /**************************/
        /* Encode each k2-raster  */
        /**************************/
        {
            size_t n_rows, n_cols, k1, k2, level_k1, plain_levels;
            std::string file_path_input;
            k2_raster_type last_snapshot;

            /*********************/
            /* Read parameters   */
            /*********************/
            read_params_from_file(inputs_file, n_rows, n_cols, k1, k2, level_k1, plain_levels);

            /*********************/
            /* Process rasters   */
            /*********************/
            while (next_file_path(inputs_file, input_path_folder, file_path_input)) {

                std::vector<value_type> values;
                read_input_data(file_path_input, n_rows, n_cols, values, scale_factor);

                /*********************/
                /* Encode data       */
                /*********************/
                if (time_count % k_snapshots_freq) {
                    // It is a log
                    k_logs.emplace_back(values, last_snapshot);
                    // Global values
                    this->k_max_size_x = std::max(this->k_max_size_x, k_logs[k_logs.size()-1].get_n_rows());
                    this->k_max_size_y = std::max(this->k_max_size_y, k_logs[k_logs.size()-1].get_n_cols());
                    this->k_min_value = std::min(this->k_min_value, k_logs[k_logs.size()-1].min_value);
                    this->k_max_value = std::max(this->k_max_value, k_logs[k_logs.size()-1].max_value);
                } else {
                    // It is an snapshot
                    k_snapshots.emplace_back(values, n_rows, n_cols, k1, k2, level_k1, plain_levels);
                    // Global values
                    this->k_max_size_x = std::max(this->k_max_size_x, k_snapshots[k_snapshots.size()-1].get_n_rows());
                    this->k_max_size_y = std::max(this->k_max_size_y, k_snapshots[k_snapshots.size()-1].get_n_cols());
                    this->k_min_value = std::min(this->k_min_value, k_snapshots[k_snapshots.size()-1].min_value);
                    this->k_max_value = std::max(this->k_max_value, k_snapshots[k_snapshots.size()-1].max_value);
                    last_snapshot = k_snapshots[k_snapshots.size() - 1];
                }
                time_count++;
            } // END WHiLE read_params
            inputs_file.close();
        } // END BLOCK encode k2-rasters
        return time_count;
    }

}; // END CLASS t_k2_raster
} // END NAMESPACE k2raster

#endif // INCLUDED_T_K2_RASTER
