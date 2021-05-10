/*
 * Created by Fernando Silva on 4/12/18.
 *
 * Copyright (C) 2018-current-year, Fernando Silva, all rights reserved.
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

#ifndef INCLUDED_K2_RASTER_
#define INCLUDED_K2_RASTER_

// Own libraries
#include <utils/dac_vector.hpp>
#include <k2_raster_base.hpp>

// Third libraries
#include <sdsl/vectors.hpp>

//! Namespace for the k2-raster library
namespace k2raster {

    template <typename t_value=int,
            typename t_bv=sdsl::bit_vector,
            typename t_rank=sdsl::rank_support_v5<1,1>,
            typename t_values_vec=sdsl::dac_vector_dp_opt<t_bv, t_rank, 3>>
    class k2_raster : public k2_raster_base<t_value, t_bv, t_rank, t_values_vec> {

    public:
        typedef k2_raster_base<t_value, t_bv, t_rank, t_values_vec>     k2_raster_p; // k2_raster parent
        typedef t_value                                                 value_type;
        typedef k2_raster_base<>::size_type                             size_type;
        typedef t_bv                                                    bit_vector_type;


    public:
        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//
        k2_raster() = default;

        k2_raster(const k2_raster& tr) : k2_raster_p()
        {
            *this = tr;
        }

        k2_raster(k2_raster&& tr)
        {
            *this = std::move(tr);
        }

        // Only interface (plain_levels unused)
        template<class Container>
        k2_raster(Container &&c_values, size_type n_rows, size_type n_cols,
                ushort k1, ushort k2, ushort level_k1, ushort plain_levels __attribute__((unused)))
        : k2_raster(c_values, n_rows, n_cols, k1, k2, level_k1) {}


        template<class Container>
        k2_raster(Container &&c_values, size_type n_rows, size_type n_cols, ushort k1, ushort k2, ushort level_k1)
        : k2_raster(n_rows, n_cols, k1, k2, level_k1, K2_RASTER_TYPE) {

            this->init_levels();

            // Build k2-raster
            build(c_values);
        }

        //*******************************************************//
        //*************** BASIC OPERATIONS **********************//
        //*******************************************************//

        //! Move assignment operator
        k2_raster& operator=(k2_raster&& tr)
        {
            if (this != &tr) {
                k2_raster_p::operator=(tr);
            }
            return *this;
        }

        //! Assignment operator
        k2_raster& operator=(const k2_raster& tr)
        {
            if (this != &tr) {
                k2_raster_p::operator=(tr);
            }
            return *this;
        }

        //! Swap operator
        void swap(k2_raster& tr)
        {
            if (this != &tr) {
                k2_raster_p::swap(tr);
            }
        }

        //! Equal operator
        bool operator==(const k2_raster& tr) const
        {
            return k2_raster_p::operator==(tr);
        }

        //*******************************************************//
        //******************** QUERIES **************************//
        //*******************************************************//

        //*****************************//
        //***** GET CELL          *****//
        //*****************************//
        value_type get_cell(size_type row, size_type col) const {

            if (row > this->k_real_size_x || col > this->k_real_size_y) {
                return 0;
            }

            // Check first level
            if (this->k_t.empty()) {
                // All cell of the matrix are equal
                return this->k_max_value;
            }

            // ************************************//
            // Searching from level 1 to level l-1 //
            // ************************************//
            size_type size = this->k_size;                      // Current matrix size
            size_type child_pos = 0;                            // Children position on bitmap k_t of the current node
            value_type max_value = this->k_max_value;           // Current max value
            {
                ushort k = this->get_k(0);
                size_type child;

                for (auto l = 1; l < this->k_height; l++) {
                    size /= k;                                  // Submatrix size
                    child = (row / size) * k + (col / size);    // Number of child (following a z-order)
                    child_pos += child;                         // Position in bitmap k_t
                    max_value -= this->get_max_value(l, child_pos);
//                            this->k_list_max[l - 1][child_pos - this->k_accum_max_values[l-1]];

                    // Check if all the cells in the subarray are equal (uniform matrix)
                    if (!this->k_t[child_pos]) {
                        return max_value;
                    }

                    // Go down one level on the tree.
                    row = row % size;                           // Update local row
                    col = col % size;                           // Update local column

                    child_pos = this->get_children_position(child_pos, l);
                    k = this->get_k(l);                         // Get current 'k'
                }
            } // END BLOCK searching from level 1 to level l - 1


            // ************************************//
            // Searching last level                //
            // ************************************//
            {
                // Search last level
                ushort k = this->get_k(this->k_height-1);
                child_pos += row * k + col;
                max_value -= this->get_max_value(this->k_height, child_pos);
                        //this->k_list_max[this->k_height - 1][child_pos - this->k_accum_max_values[this->k_height - 1]];
                return max_value;
            } // END BLOCK last levels
        }

        //*****************************//
        //***** GET CELL BY VALUE *****//
        //*****************************//
        size_type get_cells_by_value(size_type xini, size_type xend, size_type yini, size_type yend,
                                     value_type valini, value_type valend, std::vector<std::pair<size_type, size_type>> &result) {

            // Check values of root node
            if (this->k_min_value > valend ||
                this->k_max_value < valini) {
                // There is no valid value in the matrix
                return 0;
            }

            // Whole matrix if uniform
            if (this->k_t.empty()) {
                // All cells are valid
                for (auto x = xini; x <= xend; x++){
                    for (auto y= yini; y <= yend; y++){
                        result.emplace_back(x, y);
                    }
                }
                return result.size();
            }

            ushort k = this->get_k(0);
            return get_cells_by_value_helper(xini, xend, yini, yend, valini, valend,
                                             0, 0, this->k_min_value, this->k_max_value, 0, this->k_size / k, 1,
                                             result);
        }

        //*****************************//
        //***** GET VALUES WINDOW *****//
        //*****************************//
        size_type get_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                            std::vector<value_type> &result) {

            size_type count_cells = 0; // Number of cells
            result.resize((xend - xini + 1) * (yend - yini + 1));

            // Whole matrix if uniform
            if (this->k_t.empty()) {
                // All cells are valid
                value_type val = this->k_max_value;
                for (auto x = xini; x <= xend; x++){
                    for (auto y= yend; y <= yend; y++){
                        result[count_cells++] = val;
                    }
                }
                return count_cells;
            }

            ushort k = this->get_k(0);
            return get_values_window_helper(xini, xend, yini, yend, 0, 0,
                    this->k_max_value, 0, this->k_size / k, 1,
                    result, xini, yini, (yend - yini + 1));

        }

        //*****************************//
        //**** CHECK VALUES WINDOW ****//
        //*****************************//
        bool check_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                 value_type valini, value_type valend, bool strong_check) {

            // Check values of root node
            if (this->k_min_value > valend ||
                this->k_max_value < valini) {
                // There is no valid value in the matrix
                return false;
            }

            // Whole matrix if uniform
            if (this->k_t.empty()) {
                // All cells are valid
                return true;
            }

            // Size of window is equal to size of matrix
            if (xini == 0 && yini == 0 && xend == (this->k_real_size_x - 1) && yend == (this->k_real_size_y - 1)){
                return (this->k_min_value >= valini && this->k_max_value <= valend) || !strong_check;
            }

            ushort k = this->get_k(0);
            return check_values_window_helper(xini, xend, yini, yend, valini, valend,
                                             0, 0, this->k_min_value, this->k_max_value,
                                             0, this->k_size / k, 1, strong_check);
        }


    protected:
        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//

        k2_raster(size_type n_rows, size_type n_cols, ushort k1, ushort k2, ushort level_k1, ushort k2_raster_type)
                : k2_raster_p(n_rows, n_cols, k1, k2, level_k1, k2_raster_type) {}


        //*******************************************************//
        //********************** BUILD **************************//
        //*******************************************************//
        template<class Container>
        void build(Container &&c_values) {
            std::vector<std::vector<value_type>> min_values_(this->k_height);
            std::vector<std::vector<value_type>> max_values_(this->k_height);
            std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);

            build(c_values, min_values_, max_values_, tmp_t_, this->k_min_value, this->k_max_value, 0, 0, 0, this->k_size);

            // Encode max values;
            size_type n_nodes = this->build_max_values(max_values_, tmp_t_);
            max_values_.clear();
            max_values_.shrink_to_fit();

            // Encode min values
            this->build_min_values(min_values_);
            min_values_.clear();
            min_values_.shrink_to_fit();

            // Encode bitmap T (copy temporal bitmap T)
            this->build_t(tmp_t_, n_nodes);
        }

        template<class Container>
        void build(Container &&c_values, std::vector<std::vector<value_type>> &min_values_, std::vector<std::vector<value_type>> &max_values_, std::vector<sdsl::int_vector<1>> &tmp_t_,
                   value_type &min_value, value_type &max_value, size_type base_row, size_type base_col,
                   ushort level, size_type sub_size) {

            size_type pos_value, child_base_row, child_base_col;
            min_value = std::numeric_limits<value_type>::max();
            max_value = std::numeric_limits<value_type>::min();

            ushort k = this->get_k(level);
            sub_size = sub_size / k;

            /********************/
            /* NORMAL PARTITION */
            /********************/
            if (level == this->k_height) {
                /**************/
                /* LAST LEVEL */
                /**************/
                pos_value = base_row * this->k_real_size_y + base_col;
                min_value = c_values[pos_value];
                max_value = min_value;
            } else {
                /*****************/
                /* INTERNAL NODE */
                /*****************/
                std::vector<value_type> max_values_children(k * k);
                std::vector<value_type> min_values_children(k * k);

                for (uint x = 0; x < k; x++) {
                    for (uint y = 0; y < k; y++) {
                        child_base_row = base_row + x * sub_size;
                        child_base_col = base_col + y * sub_size;
                        if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {

                            // Recursive search
                            build(c_values, min_values_, max_values_, tmp_t_,
                                  min_values_children[x * k + y], max_values_children[x * k + y],
                                  child_base_row, child_base_col, level + 1, sub_size);

                            if (min_values_children[x * k + y] < min_value) {
                                min_value = min_values_children[x * k + y];
                            }
                            if (max_values_children[x * k + y] > max_value) {
                                max_value = max_values_children[x * k + y];
                            }
                        } else {
                            // Positions out of real matrix
                            min_values_children[x * k + y] = std::numeric_limits<value_type>::max();
                            max_values_children[x * k + y] = std::numeric_limits<value_type>::min();
                        }
                    } // END for y
                } // END for x

                // Check if matrix is not empty or uniform
                if (min_value < max_value) {
                    // --------------------------------------------------------------------- //
                    // Apply delta to min and max values of children.
                    // All values ​​are always equal to or greater than 0
                    // delta_max_value = max_value - child_max_value (max_value >= child_max_value)
                    // delta_min_value = min_value - child_min_value (min_value <= child_min_value)
                    // --------------------------------------------------------------------- //
                    if (level != (this->k_height-1)) {                                                 // Except last level
                        tmp_t_[level].resize(tmp_t_[level].size() + (k * k)); // TODO improve this
                    }
                    for (size_type c = 0; c < (k*k); c++) {
                        if (min_values_children[c] > max_values_children[c]) {
                            // It is an empty submatrix, push a 0
                            if ((level != (this->k_height-1))) {
                                tmp_t_[level][max_values_[level].size()] = 0;
                            }
                            max_values_[level].push_back(0);
                        } else {
                            max_values_[level].push_back(max_value - max_values_children[c]);
                            if ((level != (this->k_height-1))) {                                        // Except last level
                                if (max_values_children[c] != min_values_children[c]) {
                                    min_values_[level].push_back(min_values_children[c] - min_value);
                                    tmp_t_[level][max_values_[level].size() - 1] = 1;
                                } else {
                                    tmp_t_[level][max_values_[level].size() - 1] = 0;
                                }
                            }
                        } // END IF min_values_children[c] > max_values_children[c]
                    } // END FOR children

                } // END IF (min_value < max_value)
            } // END IF (level == (this->k_height)
        }

        //*******************************************************//
        //**************** QUERIES HELPERS **********************//
        //*******************************************************//

        //*****************************//
        //***** GET CELL BY VALUE *****//
        //*****************************//
        size_type get_cells_by_value_helper(size_type xini, size_type xend, size_type yini, size_type yend,
                                            value_type valini, value_type valend,
                                            size_type base_x, size_type base_y,
                                            value_type f_min_value, value_type f_max_value,
                                            size_type children_pos, size_type children_size, ushort level,
                                            std::vector<std::pair<size_type, size_type>> &result) {
            size_type pos;
            size_type c_base_x, c_base_y;
            value_type max_value, min_value;
            size_type ones=0;

            ushort k = this->get_k(level-1);                        // Current value of "k"
            bool is_uniform, is_leaf = level == this->k_height;     // True -> They are leaves
            size_type count_cells = 0;

            // Set limits (children with cell that overlap with interesting region)
            size_type limit_i_x = (xini - base_x) / children_size;
            size_type limit_e_x = (xend - base_x) / children_size;
            size_type limit_i_y = (yini - base_y) / children_size;
            size_type limit_e_y = (yend - base_y) / children_size;

            // Check children
            for (auto x = limit_i_x; x <= limit_e_x; x++) {
                c_base_x = base_x + x * children_size;              // Calculate base position of the current child (row)
                for (auto y = limit_i_y; y <= limit_e_y; y++) {
                    pos = children_pos + x * k + y;                 // Get position at Tree (T)
                    c_base_y = base_y + y * children_size;          // Calculate base position of the current child (column)

                    // Get max value
                    max_value = f_max_value - this->get_max_value(level, pos);
//                                this->k_list_max[level-1][pos - this->k_accum_max_values[level-1]];

                    if (valini > max_value) {
                        // Values out of range
                        continue;
                    }

                    is_uniform = is_leaf || !this->k_t[pos];
                    if (!is_uniform) {
                        // If it is not uniform or a leaf, its minimum value is obtained
                        ones = this->k_t_rank1(pos)+1;
                        min_value = f_min_value + this->get_min_value(level, ones-1);
                                //this->k_list_min[level-1].access(ones - this->k_accum_min_values[level - 1] - 1);
                    } else {
                        min_value = max_value;
                    }

                    if (valend < min_value) {
                        continue;                                   // Values out of range
                    }

                    // Current child may contain valid cells
                    // Calculate positions that overlap with interesting region
                    auto c_xini = std::max(c_base_x, xini);
                    auto c_xend = std::min(c_base_x + children_size - 1, xend);
                    auto c_yini = std::max(c_base_y, yini);
                    auto c_yend = std::min(c_base_y + children_size - 1, yend);

                    if (valini <= min_value && valend >= max_value) {
                        // Add all valid positions to final result
                        for (auto x_l = c_xini; x_l <= c_xend; x_l++) {
                            for (auto y_l = c_yini; y_l <= c_yend; y_l++) {
                                result.emplace_back(x_l, y_l);
                            }
                        }
                        count_cells += (c_xend-c_xini+1) * (c_yend-c_yini+1);
                    } else {
                        if (!is_uniform) {
                            // Go down one level on the tree
                            if (this->is_plain_level(level-1)) {
                                // Child is represented with plain values
                                ones -= (this->k_accum_min_values[level-1]+1);
                                count_cells += this->get_cells_by_value_plain(c_xini, c_xend, c_yini, c_yend, valini, valend,
                                                                        children_size, ones * children_size * children_size, max_value, result);
                            } else {
                                size_type new_children_pos = this->get_children_position_ones(ones, level);
                                count_cells += get_cells_by_value_helper(c_xini, c_xend, c_yini, c_yend,
                                                                         valini, valend,
                                                                         c_base_x, c_base_y, min_value, max_value,
                                                                         new_children_pos,
                                                                         children_size / this->get_k(level),
                                                                         level + 1, result);
                            }
                        } // END IF is not uniform
                    } // END IF values min-max
                } // END FOR y
            } // END FOR x

            return count_cells;
        }

        //*****************************//
        //***** GET VALUES WINDOW *****//
        //*****************************//
        size_type get_values_window_helper(size_type xini, size_type xend, size_type yini, size_type yend,
                                            size_type base_x, size_type base_y,
                                            value_type f_max_value,
                                            size_type children_pos, size_type children_size, ushort level,
                                            std::vector<value_type> &result,
                                            size_type or_x, size_type or_y, size_type window_size) {
            size_type pos;
            size_type c_base_x, c_base_y;
            value_type max_value;

            ushort k = this->get_k(level-1);                        // Current value of "k"
            bool is_uniform, is_leaf = level == this->k_height;     // True -> They are leaves
            size_type count_cells = 0;

            // Set limits (children with cell that overlap with the interesting region)
            size_type limit_i_x = (xini - base_x) / children_size;
            size_type limit_e_x = (xend - base_x) / children_size;
            size_type limit_i_y = (yini - base_y) / children_size;
            size_type limit_e_y = (yend - base_y) / children_size;

            // Check children
            for (auto x = limit_i_x; x <= limit_e_x; x++) {
                c_base_x = base_x + x * children_size;                                                  // Calculate base position of the current child (row)
                for (auto y = limit_i_y; y <= limit_e_y; y++) {
                    pos = children_pos + x * k + y;                                                     // Get position at Tree (T)
                    c_base_y = base_y + y * children_size;                                              // Calculate base position of the current child (column)

                    max_value = f_max_value - this->get_max_value(level, pos);
//                            this->k_list_max[level-1].access(pos - this->k_accum_max_values[level-1]);  // Get max value
                                                              // Check if it is a uniform submatrix

                    // Calculate positions with overlap with the interesting region
                    auto c_xini = std::max(c_base_x, xini);
                    auto c_xend = std::min(c_base_x + children_size - 1, xend);
                    auto c_yini = std::max(c_base_y, yini);
                    auto c_yend = std::min(c_base_y + children_size - 1, yend);

                    is_uniform = is_leaf || !this->k_t[pos];
                    if (is_uniform) {
                        size_type cell_pos;
                        for (auto x = c_xini; x <= c_xend; x++) {
                            cell_pos = (x - or_x) * window_size + (c_yini - or_y);
                            for (auto y = c_yini; y <= c_yend; y++) {
                                result[cell_pos++] = max_value;
                            } // END FOR y
                        } // END FOR x
                        count_cells += (c_xend-c_xini+1) * (c_yend-c_yini+1);
                        continue; // Go to next node
                    } else {
                        if (this->is_plain_level(level-1)) {
                            // Child is represented with plain values
                            size_type ones = (this->k_t_rank1(pos) + 1) - (this->k_accum_min_values[level-1]+1);
                            count_cells += this->get_values_window_plain(c_xini, c_xend, c_yini, c_yend,
                                                                   children_size, ones * children_size* children_size,
                                                                   max_value, result, or_x, or_y, window_size);
                        } else {
                            // Go down one level on the tree
                            size_type new_children_pos = this->get_children_position(pos, level);
                            count_cells += get_values_window_helper(c_xini, c_xend, c_yini, c_yend,
                                                                    c_base_x, c_base_y, max_value,
                                                                    new_children_pos,
                                                                    children_size / this->get_k(level),
                                                                    level + 1,
                                                                    result, or_x, or_y, window_size);
                        }
                    } // END IF is_uniform
                } // END FOR y
            } // END FOR x
            return count_cells;
        }

        //*****************************//
        //**** CHECK VALUES WINDOW ****//
        //*****************************//
        virtual bool check_values_window_helper(size_type xini, size_type xend, size_type yini, size_type yend,
                                            value_type valini, value_type valend,
                                            size_type base_x, size_type base_y,
                                            value_type f_min_value, value_type f_max_value,
                                            size_type children_pos, size_type children_size, ushort level,
                                            bool strong_check) {
            size_type pos;
            size_type c_base_x, c_base_y;
            value_type max_value, min_value;
            size_type ones;

            ushort k = this->get_k(level-1);                        // Current value of "k"
            bool is_uniform, is_leaf = level == this->k_height;     // True -> They are leaves

            // Set limits (children with cell that overlap with interesting region)
            size_type limit_i_x = (xini - base_x) / children_size;
            size_type limit_e_x = (xend - base_x) / children_size;
            size_type limit_i_y = (yini - base_y) / children_size;
            size_type limit_e_y = (yend - base_y) / children_size;

            // Check children
            for (auto x = limit_i_x; x <= limit_e_x; x++) {
                c_base_x = base_x + x * children_size;                                                  // Calculate base position of the current child (row)
                for (auto y = limit_i_y; y <= limit_e_y; y++) {
                    pos = children_pos + x * k + y;                                                     // Get position at Tree (T)
                    c_base_y = base_y + y * children_size;                                              // Calculate base position of the current child (column)

                    max_value = f_max_value - this->get_max_value(level, pos);
//                            this->k_list_max[level-1].access(pos - this->k_accum_max_values[level-1]);  // Get max value

                    if (valini > max_value) {
                        // Values out of range
                        if (strong_check) return false; // (strong) At least one cell is within the range of values
                        continue;
                    }

                    is_uniform = is_leaf || !this->k_t[pos];
                    if (!is_uniform) {
                        // If it is not uniform or a leaf, its minimum value is obtained
                        ones = this->k_t_rank1(pos) + 1;
                        min_value = f_min_value + this->get_min_value(level, ones-1);
//                                this->k_list_min[level-1].access(ones - this->k_accum_min_values[level-1] - 1);
                    } else {
                        min_value = max_value;
                    }

                    if (valend < min_value) {
                        // Values out of range
                        if (strong_check) return false; // (strong) At least one cell is not within the range of values
                        continue;
                    }

                    if (valini <= min_value && valend >= max_value) {
                        // All cells of the submatrix is within the range of values
                        if (!strong_check) return true; // (weak) At least one cell is within the range of values
                        continue;
                    } else {

                        // Calculate positions with overlap with interesting region
                        auto c_xini = std::max(c_base_x, xini);
                        auto c_xend = std::min(c_base_x + children_size - 1, xend);
                        auto c_yini = std::max(c_base_y, yini);
                        auto c_yend = std::min(c_base_y + children_size - 1, yend);

                        if (c_base_x >= c_xini && (c_base_x + children_size)  <= c_xend
                            && c_base_y >= c_yini && (c_base_y + children_size) <= c_yend) {
                            // The current submatrix is completely inside the search window
                            // and contains some cell outside the range of values

                            // (strong) At least one cell is not within the range of values
                            // (weak) At least one cell is within the range of values
                            return !strong_check;
                        }

                        // Search their children
                        // If it is uniform, the children has not children
                        if (!is_uniform) {
                            bool result;
                            if (this->is_plain_level(level-1)) {
                                // Child is represented with plain values
                                ones -= (this->k_accum_min_values[level-1]+1);
                                result = this->check_values_window_plain(c_xini, c_xend, c_yini, c_yend, valini, valend,
                                                                   children_size, ones * children_size * children_size, max_value, strong_check);
                            } else {
                                // Go down one level on the tree
                                size_type new_children_pos = this->get_children_position_ones(ones, level);
                                result = check_values_window_helper(c_xini, c_xend, c_yini, c_yend,
                                                                         valini, valend,
                                                                         c_base_x, c_base_y, min_value, max_value,
                                                                         new_children_pos,
                                                                         children_size / this->get_k(level),
                                                                         level + 1, strong_check);
                            }
                            if (strong_check && !result) return false;
                            if (!strong_check && result) return true;
                        } // END IF is not uniform
                    } // END IF check children
                } // END FOR y
            } // END FOR x

            // (strong) No values outside the range of values found
            // (weak) No values found within the range of values
            return strong_check;
        }
    };
} // END NAMESPACE k2raster

#endif // INCLUDED_K2_RASTER_BASE
