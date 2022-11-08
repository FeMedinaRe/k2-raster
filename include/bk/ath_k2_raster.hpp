/*  
 * Created by Fernando Silva on 19/9/19.
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

#ifndef INCLUDED_ATH_K2_RASTER
#define INCLUDED_ATH_K2_RASTER

#include <k2_raster.hpp>
#include <temporal/k2_raster_temporal_base.hpp>
#include <temporal/helpers/t_k2_raster_log.hpp>
#include <temporal/helpers/t_k2_raster_plain_log.hpp>
#include <temporal/helpers/t_k2_raster_heuristic_log.hpp>

//! Namespace for k2-raster library
namespace k2raster {

template <typename t_value=int,
        typename t_bv=sdsl::bit_vector,
        typename t_rank_1=sdsl::rank_support_v5<1,1>,
        typename t_values_vec=sdsl::dac_vector_dp_opt<t_bv, t_rank_1, 3>>
class ath_k2_raster : public k2_raster_temporal_base<t_value> {

public:
    typedef k2_raster_heuristic<t_value,t_bv, t_rank_1, t_values_vec>           k2_raster_type; // k2_raster type
    typedef t_k2_raster_heuristic_log<t_value, t_bv, t_rank_1, t_values_vec>    k2_raster_log_type; // k2_raster type
//    typedef k2_raster_plain<t_value,t_bv, t_rank_1, t_values_vec>           k2_raster_type; // k2_raster type
//    typedef t_k2_raster_plain_log<t_value, t_bv, t_rank_1, t_values_vec>    k2_raster_log_type; // k2_raster type
    typedef k2_raster_temporal_base<>::size_type                            size_type;
    typedef t_value                                                         value_type;
    typedef t_bv                                                            bit_vector_type;
    typedef t_rank_1                                                        rank_1_type;

protected:
    size_type                       k_snapshots_freq;
    std::vector<k2_raster_type>     k_snapshots;
    std::vector<k2_raster_log_type> k_logs;

    /****** Structures DIFF ******/
    bit_vector_type                 k_is_snapshot;            // If a raster is a snapshot (1) or log (0)
    rank_1_type                     k_is_snapshot_rank1;      // rank support for 1-bits on m_ta

public:
    //*******************************************************//
    //******************* CONSTRUCTOR ***********************//
    //*******************************************************//
    ath_k2_raster() = default;

    ath_k2_raster(const ath_k2_raster &tr) {
        *this = tr;
    }

    ath_k2_raster(ath_k2_raster &&tr) {
        *this = std::move(tr);
    }

    ath_k2_raster(std::string &inputs_filename, std::string &input_path_folder, size_type snapshots_freq, const ushort scale_factor=0) {
        k_snapshots_freq = snapshots_freq;
        this->k_k2_raster_type = ATH_K2_RASTER_TYPE;

        this->k_n_times = build(inputs_filename, input_path_folder, scale_factor);
    }

    //*******************************************************//
    //*************** BASIC OPERATIONS **********************//
    //*******************************************************//

    //! Move assignment operator
    ath_k2_raster &operator=(ath_k2_raster &&tr) {
        if (this != &tr) {
            k2_raster_temporal_base<t_value>::operator=(tr);
            k_snapshots_freq = tr.k_snapshots_freq;
            k_snapshots = std::move(tr.k_snapshots);
            k_logs = std::move(tr.k_logs);
            k_is_snapshot = std::move(tr.k_is_snapshot);
            k_is_snapshot_rank1 = std::move(tr.k_is_snapshot_rank1);
            k_is_snapshot_rank1.set_vector(&k_is_snapshot);
        }
        return *this;
    }

    //! Assignment operator
    ath_k2_raster &operator=(const ath_k2_raster &tr) {
        if (this != &tr) {
            k2_raster_temporal_base<t_value>::operator=(tr);
            k_snapshots_freq = tr.k_snapshots_freq;
            k_snapshots = tr.k_snapshots;
            k_logs = tr.k_logs;
            k_is_snapshot = tr.k_is_snapshot;
            k_is_snapshot_rank1 = tr.k_is_snapshot_rank1;
            k_is_snapshot_rank1.set_vector(&k_is_snapshot);
        }
        return *this;
    }

    //! Swap operator
    void swap(ath_k2_raster &tr) {
        if (this != &tr) {
            k2_raster_temporal_base<t_value>::swap(tr);
            std::swap(k_snapshots_freq, tr.k_snapshots_freq);
            std::swap(k_snapshots, tr.k_snapshots);
            std::swap(k_logs, tr.k_logs);
            std::swap(k_is_snapshot, tr.k_is_snapshot);
            sdsl::util::swap_support(k_is_snapshot_rank1, tr.k_is_snapshot_rank1, &k_is_snapshot, &(tr.k_is_snapshot));
        }
    }

    //! Equal operator
    bool operator==(const ath_k2_raster &tr) const {
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

        if (k_is_snapshot.size() != tr.k_is_snapshot.size()) {
            return false;
        }
        return true;
    }

    //*******************************************************//
    //******************** QUERIES **************************//
    //*******************************************************//

    //*****************************//
    //********** GET CELL *********//
    //*****************************//
    value_type get_cell(size_type row, size_type col, size_type time) const {
        assert(time <= (this->k_n_times));

        size_type pos_snap = time == 0 ? 0 : k_is_snapshot_rank1(time);

        if ( k_is_snapshot[time] == 1) {
            // It is a snapshot
            return k_snapshots[pos_snap].get_cell(row, col);
        } else {
            // It is a log
            return k_logs[time - pos_snap].get_cell(k_snapshots[pos_snap - 1], row, col);
        }
    }

    //*****************************//
    //****** GET CELLS VALUES *****//
    //*****************************//
    size_type get_cells_by_value(size_type xini, size_type xend, size_type yini, size_type yend,
                                 value_type valini, value_type valend, size_type tmin, size_type tmax,
                                 std::vector<std::vector<std::pair<size_type, size_type>>> &result) {
        assert(tmin >= 0);
        assert(tmax <= (k_snapshots.size() + k_logs.size() -1));
        assert(tmin <= tmax);

        result.resize(tmax - tmin + 1);;
        size_type count_cells = 0;
        size_type pos_snap = tmin == 0 ? 0 : k_is_snapshot_rank1(tmin);
        for (auto t = tmin; t <= tmax; t++) {
            if ( k_is_snapshot[t] == 1) {
                // It is a snapshot
                count_cells += k_snapshots[pos_snap].get_cells_by_value(xini, xend, yini, yend, valini, valend, result[t - tmin]);
                pos_snap++;
            } else {
                // It is a log
                count_cells += k_logs[t - pos_snap].get_cells_by_value(k_snapshots[pos_snap-1],
                                                                            xini, xend, yini, yend, valini, valend, result[t - tmin]);
            }
        }
        return count_cells;
    }

    //*****************************//
    //***** GET VALUES WINDOW *****//
    //*****************************//
    size_type get_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                size_type tmin, size_type tmax,
                                std::vector<std::vector<value_type>> &result) {
        assert(tmin >= 0);
        assert(tmax <= (k_snapshots.size() + k_logs.size() -1));
        assert(tmin <= tmax);

        result.resize(tmax - tmin + 1);;
        size_type count_cells = 0;
        size_type pos_snap = tmin == 0 ? 0 : k_is_snapshot_rank1(tmin);
        for (auto t = tmin; t <= tmax; t++) {
            if ( k_is_snapshot[t] == 1) {
                // It is a snapshot
                count_cells += k_snapshots[pos_snap].get_values_window(xini, xend, yini, yend, result[t - tmin]);
                pos_snap++;
            } else {
                // It is a log
                count_cells += k_logs[t - pos_snap].get_values_window(k_snapshots[pos_snap-1],
                                                                           xini, xend, yini, yend, result[t - tmin]);
            }
        }
        return count_cells;
    }

    //*****************************//
    //**** CHECK VALUES WINDOW ****//
    //*****************************//
    bool check_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                     value_type valini, value_type valend,
                                     size_type tmin, size_type tmax, bool strong_check) {
        assert(tmin >= 0);
        assert(tmax <= (k_snapshots.size() + k_logs.size() -1));
        assert(tmin <= tmax);

        bool result;
        size_type pos_snap = tmin == 0 ? 0 : k_is_snapshot_rank1(tmin);
        for (auto t = tmin; t <= tmax; t++) {
            if ( k_is_snapshot[t] == 1) {
                // It is a snapshot
                result = k_snapshots[pos_snap].check_values_window(xini, xend, yini, yend, valini, valend, strong_check);
                pos_snap++;
            } else {
                // It is a log
                result = k_logs[t - pos_snap].check_values_window(k_snapshots[pos_snap-1],
                                                                      xini, xend, yini, yend, valini, valend, strong_check);
            }
            if (strong_check){
                if (!result) return false;
            } else {
                if (result) return true;
            }
        } // END FOR t
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


        /****** at k2-raster ******/
        written_bytes += write_member(k_snapshots_freq, out, child, "snapshots_freq");
        written_bytes += k_is_snapshot.serialize(out, child, "is_snapshot");
        written_bytes += k_is_snapshot_rank1.serialize(out, child, "is_snapshot_rank1");

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
        k_is_snapshot.load(in);
        k_is_snapshot_rank1.load(in);
        k_is_snapshot_rank1.set_vector(&k_is_snapshot);

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

            size_type pos_snap = t == 0 ? 0 : k_is_snapshot_rank1(t);
            if ( k_is_snapshot[t] == 0) {
                // It is a log
                size = sdsl::size_in_mega_bytes(k_logs[t - pos_snap]);
                type = "log";
            } else {
                // It is a snapshot
                size = sdsl::size_in_mega_bytes(k_snapshots[pos_snap]);
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

        /**************************/
        /* Count rasters          */
        /**************************/
        size_type n_rasters;
        {
            std::ifstream file(inputs_filename);
            n_rasters = std::count(std::istreambuf_iterator<char>(file),
                       std::istreambuf_iterator<char>(), '\n') + 1;
            file.close();
        }

        /**************************/
        /* Encode each k2-raster  */
        /**************************/
        {
            size_t n_rows, n_cols, k1, k2, level_k1, plain_levels;
            std::string file_path_input;
            k2_raster_log_type raster_log, raster_log2, raster_log_prev;
            k2_raster_type raster_snap, raster_snap_prev, raster_snap_last; // Current, previous and last stored
            size_type size_snap, size_log=std::numeric_limits<size_type>::max(), size_log2;
            size_type size_snap_prev, size_log_prev=std::numeric_limits<size_type>::max();
            size_type min_size;

            bool has_prev_raster = false;

            // Set structures
            k_is_snapshot.resize(n_rasters);
            k_logs.clear();
            k_snapshots.clear();

#ifndef NDEBUG
            std::cout << "Reading " <<  n_rasters << " rasters..." << std::endl;
#endif
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
                raster_snap = k2_raster_type(values, n_rows, n_cols, k1, k2, level_k1, plain_levels);                    // Regular k2-raster
                size_snap = sdsl::size_in_bytes(raster_snap);
                if (time_count != 0) {
                    raster_log = k2_raster_log_type(values, raster_snap_last); // Log k2-raster
                    size_log = sdsl::size_in_bytes(raster_log);
                }

                /*********************/
                /* Global values     */
                /*********************/
                this->k_max_size_x = std::max(this->k_max_size_x, raster_snap.get_n_rows());
                this->k_max_size_y = std::max(this->k_max_size_y, raster_snap.get_n_cols());
                this->k_min_value = std::min(this->k_min_value, raster_snap.get_min_value());
                this->k_max_value = std::max(this->k_max_value, raster_snap.get_max_value());

                /*********************/
                /* Check best option */
                /*********************/
                if (!has_prev_raster) {
                    /****************************************/
                    /* A -- Previous raster is a snapshot   */
                    /****************************************/
                    if (size_log < size_snap) {
                        // Store current log for the next iteration
                        store_raster_information(raster_snap, raster_log,raster_snap_prev,raster_log_prev,has_prev_raster);
                    } else {
                        // It is an snapshot
                        k_is_snapshot[time_count] = 1;
                        k_snapshots.push_back(raster_snap);
                        has_prev_raster = false;
                        raster_snap_last = raster_snap;
                        time_count++;
                    }
                } else {
                    /****************************************/
                    /* B -- Previous raster is a log        */
                    /****************************************/

                    raster_log2 = k2_raster_log_type(values,raster_snap_prev); // Log k2-raster
                    size_log2 = sdsl::size_in_bytes(raster_log2);
                    size_snap_prev = sdsl::size_in_bytes(raster_snap_prev);
                    size_log_prev = sdsl::size_in_bytes(raster_log_prev);

                    // Check which is the minimum size
                    min_size = std::min(size_snap + size_log_prev,
                                    std::min(size_log + size_log_prev, size_log2 + size_snap_prev));

                    /****************************************/
                    /* B Option 1 -- log_prev + snapshot    */
                    /****************************************/
                    /* Previous raster is a log and the current raster a snapshot */

                    if (min_size == size_log_prev + size_snap) {
                        // Previous raster is a log
                        k_is_snapshot[time_count] = 0;
                        k_logs.push_back(raster_log_prev);
                        time_count++;

                        // Current raster is a snapshot
                        k_is_snapshot[time_count] = 1;
                        k_snapshots.push_back(raster_snap);
                        has_prev_raster = false;
                        raster_snap_last = raster_snap;
                        time_count++;
                        continue;
                    }

                    /****************************************/
                    /* B Option 2 -- log_prev + log         */
                    /****************************************/
                    /* Both raster are logs */

                    if (min_size == size_log_prev + size_log){
                        // Previous raster is a log
                        k_is_snapshot[time_count] = 0;
                        k_logs.push_back(raster_log_prev);
                        time_count++;

                        // Current raster is a log too
                        store_raster_information(raster_snap, raster_log, raster_snap_prev, raster_log_prev, has_prev_raster);
                        continue;
                    }

                    /****************************************/
                    /* B Option 3 -- snap_prev + log2       */
                    /****************************************/
                    /* Previous raster is an snapshot and the current raster a log */

                    // Previous raster is a snapshot
                    k_is_snapshot[time_count] = 1;
                    k_snapshots.push_back(raster_snap_prev);
                    raster_snap_last = raster_snap_prev;
                    time_count++;

                    // Current raster is a log (log of raster_snap_prev)
                    store_raster_information(raster_snap, raster_log2, raster_snap_prev, raster_log_prev, has_prev_raster);
                }
            } // END WHiLE read_params

            // Store last raster (If it is a log)
            if (has_prev_raster){
                // Last raster of the series
                k_is_snapshot[time_count] = 0;
                k_logs.push_back(raster_log_prev);
                time_count++;
            }

#ifndef NDEBUG
            std::cout << "Times: " << time_count << " where " << k_snapshots.size() << " are snapshots and " << k_logs.size() << " are logs" << std::endl;
#endif

            sdsl::util::init_support(this->k_is_snapshot_rank1, &this->k_is_snapshot);
            inputs_file.close();
        } // END BLOCK encode k2-rasters
        return time_count;
    }

    inline void store_raster_information(k2_raster_type &raster_snap, k2_raster_log_type &raster_log,
            k2_raster_type &raster_snap_prev, k2_raster_log_type &raster_log_prev, bool &has_prev_raster) {
        raster_log_prev = raster_log;
        raster_snap_prev = raster_snap;
        has_prev_raster = true;
    }
}; // END class t_k2_raster

// Types
    typedef k2raster::ath_k2_raster<> athk2r_type;


} // END namespacek2raster

#endif // INCLUDED_ATH_K2_RASTER
