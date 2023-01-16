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
#ifdef SHOW_CONSUMPTION_MEMORY
#include <utils/utils_memory.hpp>
#endif

// Third libraries
#include <sdsl/vectors.hpp>
#include <queue>

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

        //*****************************//
        //********** ALGEBRA **********//
        //*****************************//

        // Algebra operation between two k2-rasters
        template<class raster_type>
        k2_raster(raster_type &k2_raster_1, raster_type &k2_raster_2, const OperationRaster operation, const bool memory_opt=false) {
            this->copy_initial_parameters(k2_raster_1, k2_raster_2);

            if (memory_opt) {
                // Less memory consumption but more time
                this->build_operation_opt_mem<typename raster_type::value_type>(k2_raster_1, k2_raster_2, operation);
            } else {
                // More memory consumption but less time
                this->build_operation<typename raster_type::value_type>(k2_raster_1, k2_raster_2, operation);
            }
        }

        // Set a 1 if the cell value is bigger then a thresholding
        template<class raster_type>
        k2_raster(raster_type &k2_raster_1, long thresholding, const bool memory_opt=false) {
            this->copy_initial_parameters(k2_raster_1);

            if (memory_opt) {
                // Less memory consumption but more time
                this->build_operation_thresholding_opt_men<typename raster_type::value_type>(k2_raster_1, thresholding);
            } else {
                // More memory consumption but less time
                this->build_operation_thresholding<typename raster_type::value_type>(k2_raster_1, thresholding);
            }
        }

        // Apply a scalar operation to all cells of the raster
        template<class raster_type>
        k2_raster(raster_type &k2_raster_1, const value_type scalar_value,
                  const OperationRaster &operation, const bool memory_opt=false) {
            this->copy_initial_parameters(k2_raster_1);
            this->build_operation_scalar(k2_raster_1, scalar_value, operation, memory_opt);
        }

        // Algebra zonal operation between two k2-rasters
        template<class raster_type>
        k2_raster(raster_type &k2_raster_1, raster_type &k2_raster_z, const OperationZonalRaster &operation,
                  const bool opt=true/*, const bool memory_opt=false*/) {

            this->copy_initial_parameters(k2_raster_1, k2_raster_z);
            if (opt) {
                this->build_operation_zonal_opt<typename raster_type::value_type>(k2_raster_1, k2_raster_z,
                                                                                     operation);
            } else {
                this->build_operation_zonal<typename raster_type::value_type>(k2_raster_1, k2_raster_z, operation);
            }
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

        //*******************************************************//
        //******************* DECOMPRESS ************************//
        //*******************************************************//
        template<class Container>
        void decompress(Container&& values, size_type &n_rows, size_type &n_cols) {
            size_type count_cells = 0; // Number of cells
            n_rows = this->k_real_size_x;
            n_cols = this->k_real_size_y;
            values.resize(n_rows * n_cols);

            // Whole matrix if uniform
            if (this->k_t.empty()) {
                // All cells are valid
                value_type val = this->k_max_value;
                for (size_type x = 0; x < n_rows; x++){
                    for (size_type y= 0; y < n_cols; y++){
                        values[count_cells++] = val;
                    }
                }
                return;
            }
            ushort k = this->get_k(0);
            decompress_helper(0, 0, this->k_max_value, 0, this->k_size / k, 1, values);
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
            // Initialize the temporal structures
            std::vector<std::vector<value_type>> min_values_(this->k_height);
            std::vector<std::vector<value_type>> max_values_(this->k_height);
            std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);

            // Build conceptual tree
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
                        tmp_t_[level].resize(tmp_t_[level].size() + (k * k));
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


        //*****************************//
        //********** HELPERS **********//
        //*****************************//
        inline void insert_node_conceptual_tree(std::vector<std::vector<value_type>> &min_values_, std::vector<std::vector<value_type>> &max_values_,
                                                std::vector<sdsl::int_vector<1>> &tmp_t_,
                                                std::vector<value_type> &min_values_children, std::vector<value_type> &max_values_children,
                                                value_type &min_value, value_type &max_value,
                                                ushort level, ushort k) {
            // --------------------------------------------------------------------- //
            // Apply delta to min and max values of children.
            // All values are always equal to or greater than 0
            // delta_max_value = max_value - child_max_value (max_value >= child_max_value)
            // delta_min_value = min_value - child_min_value (min_value <= child_min_value)
            // --------------------------------------------------------------------- //
            if (level != (this->k_height-1)) {                                                 // Except last level
                tmp_t_[level].resize(tmp_t_[level].size() + (k * k));
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
        }

        template <class vector_values_type>
        inline void insert_node_conceptual_tree_opt_mem(std::vector<vector_values_type> &min_values_, std::vector<vector_values_type> &max_values_,
                                                std::vector<sdsl::int_vector<1>> &tmp_t_,
                                                size_type &n_min_values, size_type &n_max_values,
                                                std::vector<value_type> &min_values_children, std::vector<value_type> &max_values_children,
                                                value_type &min_value, value_type &max_value,
                                                ushort level, ushort k) {
            // --------------------------------------------------------------------- //
            // Apply delta to min and max values of children.
            // All values are always equal to or greater than 0
            // delta_max_value = max_value - child_max_value (max_value >= child_max_value)
            // delta_min_value = min_value - child_min_value (min_value <= child_min_value)
            // --------------------------------------------------------------------- //
            /*if (level != (this->k_height-1)) {                                                 // Except last level
                tmp_t_[level].resize(tmp_t_[level].size() + (k * k));
            }*/

            if (min_values_[level].size() < (n_min_values + ( k * k ))) {
                min_values_[level].resize(min_values_[level].size() + (k * k));
            }
            if (max_values_[level].size() < (n_max_values + ( k * k ))) {
                max_values_[level].resize(max_values_[level].size() + (k * k));
                if (level != (this->k_height-1)) {
                    tmp_t_[level].resize(tmp_t_[level].size() + (k * k));
                }
            }

            for (size_type c = 0; c < (k*k); c++) {
                if (min_values_children[c] > max_values_children[c]) {
                    // It is an empty submatrix, push a 0
                    if ((level != (this->k_height-1))) {
                        tmp_t_[level][n_max_values] = 0;
                    }
                    max_values_[level][n_max_values] = 0;
                } else {
                    max_values_[level][n_max_values] = max_value - max_values_children[c];
                    if ((level != (this->k_height-1))) {                                        // Except last level
                        if (max_values_children[c] != min_values_children[c]) {
                            min_values_[level][n_min_values] = min_values_children[c] - min_value;
                            n_min_values++;
                            tmp_t_[level][n_max_values] = 1;
                        } else {
                            tmp_t_[level][n_max_values] = 0;
                        }
                    }
                } // END IF min_values_children[c] > max_values_children[c]
                n_max_values++;
            } // END FOR children
        }

        /*inline void encode_conceptual_tree(std::vector<std::vector<value_type>> &max_values_,
                                           std::vector<std::vector<value_type>> &min_values_,
                                           std::vector<sdsl::int_vector<1>> &tmp_t_) {
            // ------------------------------------------------------------------- //
            // Encode the conceptual tree
            // ------------------------------------------------------------------- //

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
        }*/

        template <class vector_values>
        inline void encode_conceptual_tree(std::vector<vector_values> &max_values_,
                                           std::vector<vector_values> &min_values_,
                                           std::vector<sdsl::int_vector<1>> &tmp_t_) {
            // ------------------------------------------------------------------- //
            // Encode the conceptual tree
            // ------------------------------------------------------------------- //

            // Encode min values
            this->build_min_values(min_values_);
            min_values_.clear();
            min_values_.shrink_to_fit();

            // Encode max values;
            size_type n_nodes = this->build_max_values(max_values_, tmp_t_);
            max_values_.clear();
            max_values_.shrink_to_fit();

            // Encode bitmap T (copy temporal bitmap T)
            this->build_t(tmp_t_, n_nodes);
        }


        //************************************************************************************************************//
        //*********************************************** ALGEBRA ****************************************************//
        //************************************************************************************************************//
        //*******************************************************//
        //***************** Point Wise **************************//
        //*******************************************************//
        template< typename value_in_type=value_type, class raster_type>
        void build_operation(raster_type &k2_raster_1, raster_type &k2_raster_2, const OperationRaster operation){

            // ------------------------------------------------------------------- //
            // Initialize temporal structure where store max a min values and
            // the topology of the tree
            // ------------------------------------------------------------------- /
            std::vector<std::vector<value_type>> min_values_(this->k_height);
            std::vector<std::vector<value_type>> max_values_(this->k_height);
            std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);

            // ------------------------------------------------------------------- //
            // Run the algorithm of construction over the matrix
            // Return the conceptual tree that represent the k2-raster
            // ------------------------------------------------------------------- //
            k2raster_submatrix<value_in_type> node1 {k2_raster_1.max_value == k2_raster_1.min_value, k2_raster_1.max_value, 0, 0};
            k2raster_submatrix<value_in_type> node2 {k2_raster_2.max_value == k2_raster_2.min_value, k2_raster_2.max_value, 0, 0};

//            std::vector<size_type> last_access1(k2_raster_1.k_accum_min_values.size(), 0);
//            std::vector<size_type> last_access2(k2_raster_2.k_accum_min_values.size(), 0);

            std::vector<size_type> last_access1(k2_raster_1.m_accum_min_values);
            std::vector<size_type> last_access2(k2_raster_2.m_accum_min_values);

            check_submatrix_operation(min_values_, max_values_, tmp_t_,
                                      this->k_min_value, this->k_max_value,
                                      0, 0, 0, this->k_size,
                                      k2_raster_1, k2_raster_2, node1, node2,
                                      last_access1, last_access2,
                                      operation);

            // ------------------------------------------------------------------- //
            // Encode the conceptual tree
            // ------------------------------------------------------------------- //
            encode_conceptual_tree(max_values_, min_values_, tmp_t_);
        }

        template<typename value_in_type, class raster_type>
        void build_operation_opt_mem(raster_type &k2_raster_1, raster_type &k2_raster_2, const OperationRaster operation){

            // ------------------------------------------------------------------- //
            // Initialize temporal structure where store max a min values and
            // the topology of the tree
            // ------------------------------------------------------------------- /

            ushort width;
            // Calculate max width for a node value
            value_type max_value_op, min_value_op;
            switch (operation) {
                case OPERATION_SUM:
                    max_value_op = (value_type)k2_raster_1.max_value + (value_type)k2_raster_2.max_value;
                    min_value_op = (value_type)k2_raster_1.min_value + (value_type)k2_raster_2.min_value;
                    width = sdsl::bits::hi(max_value_op - min_value_op) + 1;
                    break;
                case OPERATION_SUBT:
                    max_value_op = (value_type)k2_raster_1.max_value - (value_type)k2_raster_2.min_value;
                    min_value_op = (value_type)k2_raster_1.min_value - (value_type)k2_raster_2.max_value;
                    width = sdsl::bits::hi(max_value_op - min_value_op) + 1;
                    break;
                case OPERATION_MULT:
                    max_value_op = (value_type)k2_raster_1.max_value * (value_type)k2_raster_2.max_value;
                    min_value_op = (value_type)k2_raster_1.min_value * (value_type)k2_raster_2.min_value;
                    width = sdsl::bits::hi(max_value_op - min_value_op) + 1;
                    break;
                default:
                    // No valid operation. The process never reaches this point
                    exit(-1);
            }

#ifdef MIN_CONSUMPTION_MEMORY
            std::vector<sdsl::int_vector<>> min_values_(this->k_height);
            std::vector<sdsl::int_vector<>> max_values_(this->k_height);
#else
            std::vector<std::vector<value_type>> min_values_(this->k_height);
            std::vector<std::vector<value_type>> max_values_(this->k_height);
#endif
            std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);
            std::vector<size_type> n_min_values (this->k_height, 0);
            std::vector<size_type> n_max_values (this->k_height, 0);

            // Set width and vector size
            {
                // Choose the biggest k2-raster
                raster_type *k2_raster_template;
                if (k2_raster_1.height >= k2_raster_2.height) {
                    if (k2_raster_1.height == k2_raster_2.height) {
                        if (k2_raster_1.m_accum_max_values[k2_raster_1.height - 1] >=
                            k2_raster_2.m_accum_max_values[k2_raster_2.height - 1]) {
                            k2_raster_template = &k2_raster_2;
                        } else {
                            k2_raster_template = &k2_raster_1;
                        }
                    } else {
                        k2_raster_template = &k2_raster_1;
                    }
                } else {
                    k2_raster_template = &k2_raster_2;
                }


                for (size_type l = 1; l <= max_values_.size(); l++) {
#ifdef MIN_CONSUMPTION_MEMORY
                    max_values_[l-1].width(width);
#endif
                    max_values_[l-1].resize(k2_raster_template->m_list_max[l].size());
                    if (l < this->k_height) {
                        tmp_t_[l-1].resize(k2_raster_template->m_list_max[l].size());
                    } // END IF reserve bitmap
                    if ((l-1) < min_values_.size()) {
#ifdef MIN_CONSUMPTION_MEMORY
                        min_values_[l-1].width(width);
#endif
                        min_values_[l-1].resize(k2_raster_template->m_list_min[l].size());
                    } // END IF reserve min
                } // END FOR reserve for max/min values
            }

            // ------------------------------------------------------------------- //
            // Run the algorithm of construction over the matrix
            // Return the conceptual tree that represent the k2-raster
            // ------------------------------------------------------------------- //
            k2raster_submatrix<value_in_type> node1 {k2_raster_1.max_value == k2_raster_1.min_value, k2_raster_1.max_value, 0, 0};
            k2raster_submatrix<value_in_type> node2 {k2_raster_2.max_value == k2_raster_2.min_value, k2_raster_2.max_value, 0, 0};

            std::vector<size_type> last_access1(k2_raster_1.m_accum_min_values);
            std::vector<size_type> last_access2(k2_raster_2.m_accum_min_values);

            check_submatrix_operation_opt_mem<value_in_type>(min_values_, max_values_, tmp_t_,
                                              n_min_values, n_max_values,
                                              this->k_min_value, this->k_max_value,
                                              0, 0, 0, this->k_size,
                                              k2_raster_1, k2_raster_2, node1, node2,
                                              last_access1, last_access2,
                                              operation);
#ifdef MIN_CONSUMPTION_MEMORY
            // Free input raster (release memory)
            k2_raster_1 = raster_type();
            k2_raster_2 = raster_type();
#endif


            // Resize values
            for (size_type l = 0; l < max_values_.size(); l++) {
                max_values_[l].resize(n_max_values[l]);
#ifdef MIN_CONSUMPTION_MEMORY
                sdsl::util::bit_compress(max_values_[l]);
#endif
                if (l != max_values_.size() -1) { // TODO check this
                    tmp_t_[l].resize(n_max_values[l]);
                }
                if (l < min_values_.size()) {
                    min_values_[l].resize(n_min_values[l]);
#ifdef MIN_CONSUMPTION_MEMORY
                    sdsl::util::bit_compress(min_values_[l]);
#endif
                }

            }

            // ------------------------------------------------------------------- //
            // Encode the conceptual tree
            // ------------------------------------------------------------------- //
            encode_conceptual_tree(max_values_, min_values_, tmp_t_);
        }

        template <class raster_type, typename value_in_type=value_type>
        void check_submatrix_operation(std::vector<std::vector<value_type>> &min_values_, std::vector<std::vector<value_type>> &max_values_,
                                       std::vector<sdsl::int_vector<1>> &tmp_t_,
                                       value_type &min_value, value_type &max_value, size_type base_row, size_type base_col,
                                       ushort level, size_type sub_size,
                                       raster_type &k2_raster_1, raster_type &k2_raster_2,
                                       k2raster_submatrix<value_in_type> &node1, k2raster_submatrix<value_in_type> &node2,
                                       std::vector<size_type> &last_access1, std::vector<size_type> &last_access2,
                                       const OperationRaster &operation) {


            min_value = std::numeric_limits<value_type>::max();
            max_value = std::numeric_limits<value_type>::min();

            // Calculate parameters for the current level.
            ushort k = this->get_k(level);
            sub_size = sub_size / k;

            // ------------------------------------------------------------------- //
            // If both nodes are uniforms, the new node is also uniform
            // ------------------------------------------------------------------- //
            if (node1.is_uniform && node2.is_uniform) {

                // Run operation
                switch (operation) {
                    case OPERATION_SUM:
                        max_value = (value_type)node1.max_value + (value_type)node2.max_value;
                        break;
                    case OPERATION_SUBT:
                        max_value = (value_type)node1.max_value - (value_type)node2.max_value;
                        break;
                    case OPERATION_MULT:
                        max_value = (value_type)node1.max_value * (value_type)node2.max_value;
                        break;
                    default:
                        // No valid operation. The process never reaches this point
                        exit(-1);
                }

                min_value = max_value; // Uniform node
            } else {
                /*****************/
                /* INTERNAL NODE */
                /*****************/
                std::vector<value_type> min_values_children(k * k);
                std::vector<value_type> max_values_children(k * k);
                k2raster_submatrix<value_in_type> child1, child2;
                child1.level = node1.level + 1;
                child2.level = node2.level + 1;

                size_type child_base_row, child_base_col;

                for (uint x = 0; x < k; x++) {
                    child_base_row = base_row + x * sub_size;
                    for (uint y = 0; y < k; y++) {
                        child_base_col = base_col + y * sub_size;

                        if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {

                            k2_raster_1.get_next_child(node1, child1, last_access1);
                            k2_raster_2.get_next_child(node2, child2, last_access2);

                            check_submatrix_operation(min_values_, max_values_, tmp_t_,
                                                      min_values_children[x * k + y], max_values_children[x * k + y],
                                                      child_base_row, child_base_col, level + 1, sub_size,
                                                      k2_raster_1, k2_raster_2, child1, child2 ,
                                                      last_access1, last_access2, operation);

                            // Calculate min/max value of the parent node (current node)
                            if (min_values_children[x * k + y] < min_value) {
                                min_value = min_values_children[x * k + y];
                            }
                            if (max_values_children[x * k + y] > max_value) {
                                max_value = max_values_children[x * k + y];
                            }
                        } else {
                            // More children in raster1? Reset last_access1
                            if (child_base_row < k2_raster_1.m_real_size_x && child_base_col < k2_raster_1.m_real_size_y) {
                                for (size_t l = child1.level-1; l < last_access1.size(); l++) {
                                    last_access1[l] = 0;
                                }
                            }

                           // More children in raster2? Reset last_access2
                           if (child_base_row < k2_raster_2.m_real_size_x && child_base_col < k2_raster_2.m_real_size_y) {
                               for (size_t l = child2.level-1; l < last_access2.size(); l++) {
                                   last_access2[l] = 0;
                               }
                            }

                            // Positions out of real matrix
                            min_values_children[x * k + y] = std::numeric_limits<value_type>::max();
                            max_values_children[x * k + y] = std::numeric_limits<value_type>::min();
                        }// ENF IF within real matrix

                        node1.children_pos++;
                        node2.children_pos++;
                    } // END for y
                } // END for x

                // Check if matrix is not empty or uniform
                if (min_value < max_value) {
                    insert_node_conceptual_tree(min_values_, max_values_, tmp_t_,
                                                min_values_children, max_values_children,
                                                min_value, max_value, level, k);
                } // END IF (min_value < max_value)
            } // END IF leaf-internal node
        }

        template <typename value_in_type=value_type, class vector_values_type, class raster_type>
        void check_submatrix_operation_opt_mem(std::vector<vector_values_type> &min_values_, std::vector<vector_values_type> &max_values_,
                                               std::vector<sdsl::int_vector<1>> &tmp_t_,
                                               std::vector<size_type> &n_min_values, std::vector<size_type> &n_max_values,
                                       value_type &min_value, value_type &max_value,
                                       size_type base_row, size_type base_col,
                                       ushort level, size_type sub_size,
                                       raster_type &k2_raster_1, raster_type &k2_raster_2,
                                       k2raster_submatrix<value_in_type> &node1, k2raster_submatrix<value_in_type> &node2,
                                       std::vector<size_type> &last_access1, std::vector<size_type> &last_access2,
                                       const OperationRaster &operation) {


            min_value = std::numeric_limits<value_type>::max();
            max_value = std::numeric_limits<value_type>::min();

            // Calculate parameters for the current level.
            ushort k = this->get_k(level);
            sub_size = sub_size / k;

            // ------------------------------------------------------------------- //
            // If both nodes are uniforms, the new node is also uniform
            // ------------------------------------------------------------------- //
            if (node1.is_uniform && node2.is_uniform) {

                // Run operation
                switch (operation) {
                    case OPERATION_SUM:
                        max_value = (value_type)node1.max_value + (value_type)node2.max_value;
                        break;
                    case OPERATION_SUBT:
                        max_value = (value_type)node1.max_value - (value_type)node2.max_value;
                        break;
                    case OPERATION_MULT:
                        max_value = (value_type)node1.max_value * (value_type)node2.max_value;
                        break;
                    default:
                        // No valid operation. The process never reaches this point
                        exit(-1);
                }
                min_value = max_value; // Uniform node
            } else {
                /*****************/
                /* INTERNAL NODE */
                /*****************/
                std::vector<value_type> min_values_children(k * k);
                std::vector<value_type> max_values_children(k * k);
                k2raster_submatrix<value_in_type> child1, child2;
                child1.level = node1.level + 1;
                child2.level = node2.level + 1;

                size_type child_base_row, child_base_col;

                for (uint x = 0; x < k; x++) {
                    child_base_row = base_row + x * sub_size;
                    for (uint y = 0; y < k; y++) {
                        child_base_col = base_col + y * sub_size;

                        if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {

                            k2_raster_1.get_next_child(node1, child1, last_access1);
                            k2_raster_2.get_next_child(node2, child2, last_access2);

                            check_submatrix_operation_opt_mem<value_in_type>(min_values_, max_values_, tmp_t_,
                                                              n_min_values, n_max_values,
                                                      min_values_children[x * k + y], max_values_children[x * k + y],
                                                      child_base_row, child_base_col, level + 1, sub_size,
                                                      k2_raster_1, k2_raster_2, child1, child2 ,
                                                      last_access1, last_access2, operation);

                            // Calculate min/max value of the parent node (current node)
                            if (min_values_children[x * k + y] < min_value) {
                                min_value = min_values_children[x * k + y];
                            }
                            if (max_values_children[x * k + y] > max_value) {
                                max_value = max_values_children[x * k + y];
                            }
                        } else {
                            // More children in raster1? Reset last_access1
                            if (child_base_row < k2_raster_1.m_real_size_x && child_base_col < k2_raster_1.m_real_size_y) {
                                for (size_t l = child1.level-1; l < last_access1.size(); l++) {
                                    last_access1[l] = 0;
                                }
                            }

                            // More children in raster2? Reset last_access2
                            if (child_base_row < k2_raster_2.m_real_size_x && child_base_col < k2_raster_2.m_real_size_y) {
                                for (size_t l = child2.level-1; l < last_access2.size(); l++) {
                                    last_access2[l] = 0;
                                }
                            }

                            // Positions out of real matrix
                            min_values_children[x * k + y] = std::numeric_limits<value_type>::max();
                            max_values_children[x * k + y] = std::numeric_limits<value_type>::min();
                        }// ENF IF within real matrix

                        node1.children_pos++;
                        node2.children_pos++;
                    } // END for y
                } // END for x

                // Check if matrix is not empty or uniform
                if (min_value < max_value) {
                    insert_node_conceptual_tree_opt_mem(min_values_, max_values_, tmp_t_,
                                                        n_min_values[level], n_max_values[level],
                                                min_values_children, max_values_children,
                                                min_value, max_value, level, k);
                } // END IF (min_value < max_value)
            } // END IF leaf-internal node
        }

        //*******************************************************//
        //**************** Thresholding *************************//
        //*******************************************************//
        template <typename value_in_type=value_type, class raster_type>
        void build_operation_thresholding(raster_type &k2_raster_1, const value_in_type threshold){

            // ------------------------------------------------------------------- //
            // Initialize temporal structure where store max a min values and
            // the topology of the tree
            // ------------------------------------------------------------------- /
            std::vector<std::vector<value_type>> min_values_(this->k_height);
            std::vector<std::vector<value_type>> max_values_(this->k_height);
            std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);

            // ------------------------------------------------------------------- //
            // Run the algorithm of construction over the matrix
            // Return the conceptual tree that represent the k2-raster
            // ------------------------------------------------------------------- //
            k2raster_node<value_in_type> node1 = k2_raster_1.get_root();

            std::vector<size_type> last_access1(k2_raster_1.m_accum_min_values);

            check_submatrix_operation<value_in_type>(threshold, min_values_, max_values_, tmp_t_,
                                      this->k_min_value, this->k_max_value,
                                      0, 0, 0, this->k_size,
                                      k2_raster_1, node1, last_access1);


            // ------------------------------------------------------------------- //
            // Encode the conceptual tree
            // ------------------------------------------------------------------- //
            encode_conceptual_tree(max_values_, min_values_, tmp_t_);
        }

        template <typename value_in_type=value_type, class raster_type>
        void build_operation_thresholding_opt_men(raster_type &k2_raster_1, const value_in_type threshold){

            // ------------------------------------------------------------------- //
            // Initialize temporal structure where store max a min values and
            // the topology of the tree
            // ------------------------------------------------------------------- /
            std::vector<sdsl::int_vector<1>> min_values_(this->k_height);
            std::vector<sdsl::int_vector<1>> max_values_(this->k_height);
            std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);
            std::vector<size_type> n_min_values (this->k_height, 0);
            std::vector<size_type> n_max_values (this->k_height, 0);

            /*
            ushort l=1;
            for (auto &vector : min_values_) {
                vector.resize(k2_raster_1.m_accum_min_values[l]);
                l++;
            }

            l=1;
            for (auto &vector : max_values_) {
                vector.resize(k2_raster_1.m_accum_max_values[l]);
                if (l < this->k_height) {
                    tmp_t_[l - 1].resize(k2_raster_1.m_accum_max_values[l]);
                }
                l++;
            }
             */

            // ------------------------------------------------------------------- //
            // Run the algorithm of construction over the matrix
            // Return the conceptual tree that represent the k2-raster
            // ------------------------------------------------------------------- //
            k2raster_node<value_in_type> node1 = k2_raster_1.get_root();
            std::vector<size_type> last_access1(k2_raster_1.m_accum_min_values);
            check_submatrix_operation_opt_men<value_in_type>(threshold, min_values_, max_values_, tmp_t_,
                                              n_min_values, n_max_values,
                                              this->k_min_value, this->k_max_value,
                                              0, 0, 0, this->k_size,
                                              k2_raster_1, node1, last_access1);


            // Free input raster (release memory)
            k2_raster_1 = raster_type();

            // Resize values
            for (size_type l = 0; l < min_values_.size(); l++) {
                min_values_[l].resize(n_min_values[l]);
                sdsl::util::bit_compress(min_values_[l]);
            }
            for (size_type l = 0; l < max_values_.size(); l++) {
                max_values_[l].resize(n_max_values[l]);
                sdsl::util::bit_compress(max_values_[l]);
                if (l != max_values_.size() -1) { // TODO check this
                    tmp_t_[l].resize(n_max_values[l]);
                }
            }

            // ------------------------------------------------------------------- //
            // Encode the conceptual tree
            // ------------------------------------------------------------------- //
            encode_conceptual_tree(max_values_, min_values_, tmp_t_);
        }

        template <typename value_in_type=value_type, class vector_values_type, class raster_type>
        void check_submatrix_operation(const value_in_type threshold, std::vector<vector_values_type> &min_values_, std::vector<vector_values_type> &max_values_,
                                       std::vector<sdsl::int_vector<1>> &tmp_t_,
                                       value_type &min_value, value_type &max_value, size_type base_row, size_type base_col,
                                       ushort level, size_type sub_size,
                                       raster_type &k2_raster_1, k2raster_node<value_in_type> &node1, std::vector<size_type> &last_access1) {

            min_value = std::numeric_limits<value_type>::max();
            max_value = std::numeric_limits<value_type>::min();

            // Calculate parameters for the current level.
            ushort k = this->get_k(level);
            sub_size = sub_size / k;

            // ------------------------------------------------------------------- //
            // If both nodes are uniforms, the new node is also uniform
            // ------------------------------------------------------------------- //
            if (node1.min_value == node1.max_value ||
                    node1.max_value < threshold || node1.min_value >= threshold) {

                // Run operation
                max_value = node1.max_value >= threshold ? 1 : 0;
                min_value = max_value; // Uniform node
            } else {
                /*****************/
                /* INTERNAL NODE */
                /*****************/
                std::vector<value_type> min_values_children(k * k);
                std::vector<value_type> max_values_children(k * k);
                k2raster_node<value_in_type> child1;
                child1.level = node1.level + 1;

                size_type child_base_row, child_base_col;

                for (uint x = 0; x < k; x++) {
                    child_base_row = base_row + x * sub_size;
                    for (uint y = 0; y < k; y++) {
                        child_base_col = base_col + y * sub_size;

                        if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {

                            child1 = k2_raster_1.get_child(node1, x * k + y);

                            check_submatrix_operation<value_in_type>(threshold, min_values_, max_values_, tmp_t_,
                                                      min_values_children[x * k + y], max_values_children[x * k + y],
                                                      child_base_row, child_base_col, level + 1, sub_size,
                                                      k2_raster_1, child1, last_access1);

                            // Calculate min/max value of the parent node (current node)
                            if (min_values_children[x * k + y] < min_value) {
                                min_value = min_values_children[x * k + y];
                            }
                            if (max_values_children[x * k + y] > max_value) {
                                max_value = max_values_children[x * k + y];
                            }
                        } else {
                            // More children in raster1? Reset last_access1
                            //if (child_base_row < k2_raster_1.k_real_size_x && child_base_col < k2_raster_1.k_real_size_y) {
                            //    for (size_t l = child1.level-1; l < last_access1.size(); l++) {
                            //        last_access1[l] = 0;
                            //    }
                            //}
                            // Positions out of real matrix
                            min_values_children[x * k + y] = std::numeric_limits<value_type>::max();
                            max_values_children[x * k + y] = std::numeric_limits<value_type>::min();
                        }// ENF IF within real matrix

                        //node1.children_pos++;
                    } // END for y
                } // END for x

                // Check if matrix is not empty or uniform
                if (min_value < max_value) {
                    insert_node_conceptual_tree(min_values_, max_values_, tmp_t_,
                                                min_values_children, max_values_children,
                                                min_value, max_value, level, k);

                } // END IF (min_value < max_value)
            } // END IF leaf-internal node
        }

        template <typename value_in_type=value_type, class vector_values_type, class raster_type>
        void check_submatrix_operation_opt_men(const value_in_type threshold, std::vector<vector_values_type> &min_values_, std::vector<vector_values_type> &max_values_,
                                               std::vector<sdsl::int_vector<1>> &tmp_t_,
                                               std::vector<size_type> &n_min_values, std::vector<size_type> &n_max_values,
                                       value_type &min_value, value_type &max_value, size_type base_row, size_type base_col,
                                       ushort level, size_type sub_size,
                                       raster_type &k2_raster_1, k2raster_node<value_in_type> &node1, std::vector<size_type> &last_access1) {

            min_value = std::numeric_limits<value_type>::max();
            max_value = std::numeric_limits<value_type>::min();

            // Calculate parameters for the current level.
            ushort k = this->get_k(level);
            sub_size = sub_size / k;

            // ------------------------------------------------------------------- //
            // If both nodes are uniforms, the new node is also uniform
            // ------------------------------------------------------------------- //
            if (node1.min_value == node1.max_value ||
                node1.max_value < threshold || node1.min_value >= threshold) {

                // Run operation
                max_value = node1.max_value >= threshold ? 1 : 0;
                min_value = max_value; // Uniform node
            } else {
                /*****************/
                /* INTERNAL NODE */
                /*****************/
                std::vector<value_type> min_values_children(k * k);
                std::vector<value_type> max_values_children(k * k);
                k2raster_node<value_in_type> child1;
                child1.level = node1.level + 1;

                size_type child_base_row, child_base_col;

                for (uint x = 0; x < k; x++) {
                    child_base_row = base_row + x * sub_size;
                    for (uint y = 0; y < k; y++) {
                        child_base_col = base_col + y * sub_size;

                        if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {

                            child1 = k2_raster_1.get_child(node1, x * k + y);

                            check_submatrix_operation_opt_men<value_in_type>(threshold, min_values_, max_values_, tmp_t_,
                                                      n_min_values, n_max_values,
                                                      min_values_children[x * k + y], max_values_children[x * k + y],
                                                      child_base_row, child_base_col, level + 1, sub_size,
                                                      k2_raster_1, child1, last_access1);

                            // Calculate min/max value of the parent node (current node)
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
                        }// ENF IF within real matrix

                        //node1.children_pos++;
                    } // END for y
                } // END for x

                // Check if matrix is not empty or uniform
                if (min_value < max_value) {
                    insert_node_conceptual_tree_opt_mem(min_values_, max_values_, tmp_t_,
                                                        n_min_values[level], n_max_values[level],
                                                        min_values_children, max_values_children,
                                                        min_value, max_value, level, k);
                } // END IF (min_value < max_value)
            } // END IF leaf-internal node
        }

        //*****************************//
        //*********** Scalar **********//
        //*****************************//

        template <class raster_type, typename value_in_type=value_type>
        void build_operation_scalar(raster_type &k2_raster_1, const value_type scalar,
                                    const OperationRaster &operation, const bool memory_opt=false){

            value_type min_value = k2_raster_1.min_value;
            value_type max_value = k2_raster_1.max_value;
            // Run operation
            switch (operation) {
                case OPERATION_SUM:
                    this->copy(k2_raster_1);
                    this->k_min_value = min_value + scalar;
                    this->k_max_value = max_value + scalar;
                    break;
                case OPERATION_SUBT:
                    this->copy(k2_raster_1);
                    this->k_min_value = min_value - scalar;
                    this->k_max_value = max_value - scalar;
                    break;
                case OPERATION_MULT:
                    this->copy_topology(k2_raster_1);

                    // ------------------------------------------------------------------- //
                    // Encode min/max values
                    // ------------------------------------------------------------------- //

                    // global Min-Max values
                    this->k_min_value = min_value * scalar;
                    this->k_max_value = max_value * scalar;

                    if (memory_opt) {
                        this->build_operation_scalar_mult_opt_men(k2_raster_1, scalar);
                    } else {
                        this->build_operation_scalar_mult(k2_raster_1, scalar);
                    }
                    break;
                default:
                    // No valid operation. The process never reaches this point
                    exit(-1);
            }
        }

        template <class vector_values_type=std::vector<value_type>, class raster_type>
        void build_operation_scalar_mult(raster_type &k2_raster_1, const value_type scalar) {

            ushort level = 0;
            // ---------------- //
            // Min values       //
            // ---------------- //
            for (auto &values_level : k2_raster_1.m_list_min) {
                vector_values_type values_min(values_level.size());
                values_min[0] = values_level.access(0) * scalar;
                for (size_t p = 1; p < values_level.size(); p++) {
                    values_min[p] = values_level.next() * scalar;
                }
                this->k_list_min[level] = t_values_vec(values_min);
                level++;
            }

            // Encode min values
            //this->build_min_values(values_);
            level = 0;

            // ---------------- //
            // Max values       //
            // ---------------- //
            for (auto &values_level : k2_raster_1.m_list_max) {
                vector_values_type values_max(values_level.size());
                values_max[0] = values_level.access(0) * scalar;
                for (size_t p = 1; p < values_level.size(); p++) {
                    values_max[p] = values_level.next() * scalar;
                }
                this->k_list_max[level] = t_values_vec(values_max);
                level++;
            }

            // Encode min values
            //this->build_max_values(values_);
            //values_.clear();
            //values_.shrink_to_fit();
        }

        template <class vector_values_type=sdsl::int_vector<>, class raster_type>
        void build_operation_scalar_mult_opt_men(raster_type &k2_raster_1, const value_type scalar) {

            ushort level = 0;

            // Calculate max width
            size_type width = sdsl::bits::hi(this->max_value - this->min_value) + 1;

            // ---------------- //
            // Min values       //
            // ---------------- //
            // Run scalar operations (min values)
            for (auto &values_level : k2_raster_1.m_list_min) {
                vector_values_type values_min(values_level.size(), 0, width);
                values_min[0] = values_level.access(0) * scalar;
                for (size_t p = 1; p < values_level.size(); p++) {
                    values_min[p] = values_level.next() * scalar;
                }
                values_level = t_values_vec();
                this->k_list_min[level] = t_values_vec(values_min);
                level++;
            }

            // ---------------- //
            // Max values       //
            // ---------------- //
            // Run scalar operations (min values)
            level = 0;
            for (auto &values_level : k2_raster_1.m_list_max) {
                vector_values_type values_max(values_level.size(), 0, width);
                values_max[0] = values_level.access(0) * scalar;
                for (size_t p = 1; p < values_level.size(); p++) {
                    values_max[p] = values_level.next() * scalar;
                }
                values_level = t_values_vec();
                this->k_list_max[level] = t_values_vec(values_max);
                level++;
            }

            // Free k2-raster
            k2_raster_1 = raster_type();
       }

        //*****************************//
        //****** Zonal Operation ******//
        //*****************************//
        template <typename value_in_type=value_type, class raster_type>
        void build_operation_zonal(raster_type &k2_raster_1, raster_type &raster_zonal, const OperationZonalRaster &operation){

            // Hash to store the sum of value for each zone
            std::map<value_in_type, value_type> zonal_sum;

            // ------------------------------------------------------------------- //
            // Run the algorithm of construction over the matrix
            // Return the conceptual tree that represent the k2-raster
            // ------------------------------------------------------------------- //
            {
                k2raster_submatrix<value_in_type> node1{k2_raster_1.max_value == k2_raster_1.min_value,
                                                     k2_raster_1.max_value, 0, 0};
                k2raster_submatrix<value_in_type> node2{raster_zonal.max_value == raster_zonal.min_value,
                                                     raster_zonal.max_value, 0, 0};

                std::vector<size_type> last_access1(k2_raster_1.m_accum_min_values);
                std::vector<size_type> last_access2(raster_zonal.m_accum_min_values);

                calculate_operation_zonal<value_in_type>(0, 0, 0, this->k_size,
                                          k2_raster_1, raster_zonal, node1, node2,
                                          last_access1, last_access2,
                                          operation, zonal_sum);
            } // END BLOCK map calculation

            // ------------------------------------------------------------------- //
            // Initialize temporal structure where store max a min values and
            // the topology of the tree
            // ------------------------------------------------------------------- /
            std::vector<std::vector<value_type>> min_values_(this->k_height);
            std::vector<std::vector<value_type>> max_values_(this->k_height);
            std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);

            // ------------------------------------------------------------------- //
            // Run the algorithm of construction over the matrix
            // Return the conceptual tree that represent the k2-raster
            // ------------------------------------------------------------------- //
            k2raster_submatrix<value_in_type> node_z {raster_zonal.max_value == raster_zonal.min_value,
                                                      raster_zonal.max_value, 0, 0};
            std::vector<size_type> last_access_z(raster_zonal.m_accum_min_values);

            create_zonal_raster<value_in_type>(min_values_, max_values_, tmp_t_,
                                this->k_min_value, this->k_max_value,
                                0, 0, 0, this->k_size,
                                raster_zonal, node_z, last_access_z, zonal_sum);

            // ------------------------------------------------------------------- //
            // Encode the conceptual tree
            // ------------------------------------------------------------------- //
            encode_conceptual_tree(max_values_, min_values_, tmp_t_);
        }

        template <typename value_in_type=value_type, class raster_type>
        void build_operation_zonal_opt(raster_type &k2_raster_1, raster_type &raster_zonal, const OperationZonalRaster &operation){

            // Hash to store the sum of value for each zone
            std::map<value_in_type, value_type> zonal_sum;
            std::vector<std::vector<std::pair<uint8_t, value_type>>> queue_nodes(raster_zonal.m_list_max.size());
            // Type format <node_type, max_value>
            // Types:
            //  - 0 -> No uniform node
            //  - 1 -> Uniform node
            //  - 2 -> Node out of real limits.
            std::vector<std::vector<size_type>> parent_nodes(raster_zonal.m_list_max.size()-1);
            // Reserve memory
            for (size_type l = 0; l < raster_zonal.m_list_max.size(); l++) {
                queue_nodes[l].reserve(raster_zonal.m_list_max[l].size());      // Possible number of 'nodes' per level.
                if (l != raster_zonal.m_list_max.size()-1) {
                    parent_nodes[l].reserve(raster_zonal.m_list_min[l].size()); // Possible number of 'parents' per level.
                }
            }
            // ------------------------------------------------------------------- //
            // Run the algorithm of construction over the matrix
            // Return the conceptual tree that represent the k2-raster
            // ------------------------------------------------------------------- //
            {
                k2raster_submatrix<value_in_type> node1{k2_raster_1.max_value == k2_raster_1.min_value,
                                                     k2_raster_1.max_value, 0, 0};
                k2raster_submatrix<value_in_type> node2{raster_zonal.max_value == raster_zonal.min_value,
                                                     raster_zonal.max_value, 0, 0};

                std::vector<size_type> last_access1(k2_raster_1.m_accum_min_values);
                std::vector<size_type> last_access2(raster_zonal.m_accum_min_values);

                calculate_operation_zonal_op<value_in_type>(0, 0, 0, this->k_size,
                                             k2_raster_1, raster_zonal, node1, node2,
                                             last_access1, last_access2,
                                             operation, zonal_sum, queue_nodes, parent_nodes);
                for (size_type l = 0; l < raster_zonal.m_list_max.size(); l++) {
                    queue_nodes[l].shrink_to_fit();
                    if (l != raster_zonal.m_list_max.size()-1) {
                        parent_nodes[l].shrink_to_fit();
                    }
                }
            } // END BLOCK map calculation

            // ------------------------------------------------------------------- //
            // Initialize temporal structure where store max a min values and
            // the topology of the tree
            // ------------------------------------------------------------------- /
            std::vector<std::vector<value_type>> min_values_(this->k_height);
            std::vector<std::vector<value_type>> max_values_(this->k_height);
            std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);

            // Reserve memory
            for (size_type l = 0; l < raster_zonal.m_list_max.size(); l++) {
                max_values_[l].reserve(queue_nodes[l].size());
                if (l != raster_zonal.m_list_max.size()-1) {
                    min_values_[l].reserve(parent_nodes[l].size());
                }
            }

            // ------------------------------------------------------------------- //
            // Run the algorithm of construction over the matrix
            // Return the conceptual tree that represent the k2-raster
            // ------------------------------------------------------------------- //
            create_zonal_raster_op<value_in_type>(min_values_, max_values_, tmp_t_,
                                   zonal_sum, queue_nodes, parent_nodes);


            // Free memory
            {
                zonal_sum.clear();
                queue_nodes.clear();
                queue_nodes.shrink_to_fit();
                parent_nodes.clear();
                parent_nodes.shrink_to_fit();

                for (size_type l = 0; l < min_values_.size(); l++) {
                    min_values_[l].shrink_to_fit();
                }

                for (size_type l = 0; l < max_values_.size(); l++) {
                    max_values_[l].shrink_to_fit();
                }
            } // END BLOCK free memory

            // ------------------------------------------------------------------- //
            // Encode the conceptual tree
            // ------------------------------------------------------------------- //
            encode_conceptual_tree(max_values_, min_values_, tmp_t_);
        }

        template <typename value_in_type=value_type, class raster_type>
        void calculate_operation_zonal(size_type base_row, size_type base_col,
                                       ushort level, size_type sub_size,
                                       raster_type &k2_raster_1, raster_type &raster_zonal,
                                       k2raster_submatrix<value_in_type> &node1, k2raster_submatrix<value_in_type> &nodeZ,
                                       std::vector<size_type> &last_access1, std::vector<size_type> &last_access2,
                                       const OperationZonalRaster &operation, std::map<value_in_type, value_type> &zonal_sum) {

            // ------------------------------------------------------------------- //
            // If both nodes are uniforms, the new node is also uniform
            // ------------------------------------------------------------------- //
            if (node1.is_uniform && nodeZ.is_uniform) {
                size_type n_nodes = 1;
                if (sub_size != 1) {
                    size_t c_x = std::min(sub_size, this->k_real_size_x - base_row);
                    size_t c_y = std::min(sub_size, this->k_real_size_y - base_col);
                    n_nodes = c_x * c_y;
                }
                // Run operation
                switch (operation) {
                    case OPERATION_ZONAL_SUM:
                        if (zonal_sum.find(nodeZ.max_value) == zonal_sum.end()) {
                            zonal_sum[nodeZ.max_value] = node1.max_value * n_nodes;
                        } else {
                            zonal_sum[nodeZ.max_value] += node1.max_value * n_nodes;
                        }
                        break;
                    default:
                        // No valid operation. The process never reaches this point
                        exit(-1);
                }
            } else {

                // Calculate parameters for the current level.
                ushort k = this->get_k(level);
                sub_size = sub_size / k;

                /*****************/
                /* INTERNAL NODE */
                /*****************/
                k2raster_submatrix<value_in_type> child1, childZ;
                child1.level = node1.level + 1;
                childZ.level = nodeZ.level + 1;
                size_type child_base_row, child_base_col;

                for (uint x = 0; x < k; x++) {
                    child_base_row = base_row + x * sub_size;
                    for (uint y = 0; y < k; y++) {
                        child_base_col = base_col + y * sub_size;

                        if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {
                            k2_raster_1.get_next_child(node1, child1, last_access1);
                            raster_zonal.get_next_child(nodeZ, childZ, last_access2);

                            calculate_operation_zonal<value_in_type>(child_base_row, child_base_col, level + 1, sub_size,
                                                      k2_raster_1, raster_zonal, child1, childZ,
                                                      last_access1, last_access2, operation, zonal_sum);
                        } else {
                            // More children in raster1? Reset last_access1
                            if (child_base_row < k2_raster_1.m_real_size_x && child_base_col < k2_raster_1.m_real_size_y) {
                                for (size_t l = child1.level-1; l < last_access1.size(); l++) {
                                    last_access1[l] = 0;
                                }
                            }

                            // More children in raster2? Reset last_access2
                            if (child_base_row < raster_zonal.m_real_size_x && child_base_col < raster_zonal.m_real_size_y) {
                                for (size_t l = childZ.level-1; l < last_access2.size(); l++) {
                                    last_access2[l] = 0;
                                }
                            }
                        }// ENF IF within real matrix

                        node1.children_pos++;
                        nodeZ.children_pos++;
                    } // END for y
                } // END for x
            } // END IF leaf-internal node
        }

        template <typename value_in_type=value_type, class raster_type>
        void calculate_operation_zonal_op(size_type base_row, size_type base_col,
                                       ushort level, size_type sub_size,
                                       raster_type &k2_raster_1, raster_type &raster_zonal,
                                       k2raster_submatrix<value_in_type> &node1, k2raster_submatrix<value_in_type> &nodeZ,
                                       std::vector<size_type> &last_access1, std::vector<size_type> &last_access2,
                                       const OperationZonalRaster &operation, std::map<value_in_type, value_type> &zonal_sum,
                                       std::vector<std::vector<std::pair<uint8_t, value_type>>> &queue_nodes,
                                       std::vector<std::vector<size_type>> &parent_nodes) {

            // ------------------------------------------------------------------- //
            // If both nodes are uniforms, the new node is also uniform
            // ------------------------------------------------------------------- //
            if (node1.is_uniform && nodeZ.is_uniform) {
                size_type n_cells = 1;
                if (sub_size != 1) {
                    // Calculate the number of cells in the node
                    size_t c_x = std::min(sub_size, this->k_real_size_x - base_row);
                    size_t c_y = std::min(sub_size, this->k_real_size_y - base_col);
                    n_cells = c_x * c_y;
                }
                // Run operation
                switch (operation) {
                    case OPERATION_ZONAL_SUM:
                        if (zonal_sum.find(nodeZ.max_value) == zonal_sum.end()) {
                            zonal_sum[nodeZ.max_value] = node1.max_value * n_cells;
                        } else {
                            zonal_sum[nodeZ.max_value] += node1.max_value * n_cells;
                        }
                        break;
                    default:
                        // No valid operation. The process never reaches this point
                        exit(-1);
                }
            } else {

                // Calculate parameters for the current level.
                ushort k = this->get_k(level);
                sub_size = sub_size / k;

                /*****************/
                /* INTERNAL NODE */
                /*****************/
                k2raster_submatrix<value_in_type> child1, childZ;
                child1.level = node1.level + 1;
                childZ.level = nodeZ.level + 1;
                size_type child_base_row, child_base_col;

                for (uint x = 0; x < k; x++) {
                    child_base_row = base_row + x * sub_size;
                    for (uint y = 0; y < k; y++) {
                        child_base_col = base_col + y * sub_size;

                        if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {
                            k2_raster_1.get_next_child(node1, child1, last_access1);

                            if (!nodeZ.is_uniform) {
                                raster_zonal.get_next_child(nodeZ, childZ, last_access2);
                                queue_nodes[nodeZ.level].push_back({childZ.is_uniform, childZ.max_value});
                                if (!childZ.is_uniform){
                                    size_type pos = nodeZ.children_pos - raster_zonal.m_accum_max_values[nodeZ.level];
                                    parent_nodes[nodeZ.level].push_back(pos);
                                }
                            } else {
                                childZ = nodeZ;
                            }

                            calculate_operation_zonal_op<value_in_type>(child_base_row, child_base_col, level + 1, sub_size,
                                                      k2_raster_1, raster_zonal, child1, childZ,
                                                      last_access1, last_access2, operation,
                                                      zonal_sum, queue_nodes, parent_nodes);
                        } else {
                            if (!nodeZ.is_uniform) {
                                // Add a fake node (out of real raster)
                                queue_nodes[nodeZ.level].push_back({2, std::numeric_limits<value_type>::min()});
                            }
                            // More children in raster1? Reset last_access1 (especial case)
                            if (child_base_row < k2_raster_1.m_real_size_x && child_base_col < k2_raster_1.m_real_size_y) {
                                for (size_t l = child1.level-1; l < last_access1.size(); l++) {
                                    last_access1[l] = 0;
                                }
                            }

                            // More children in raster2? Reset last_access2 (especial case)
                            if (child_base_row < raster_zonal.m_real_size_x && child_base_col < raster_zonal.m_real_size_y) {
                                for (size_t l = childZ.level-1; l < last_access2.size(); l++) {
                                    last_access2[l] = 0;
                                }
                            }
                        }// ENF IF within real matrix

                        node1.children_pos++;
                        nodeZ.children_pos++;
                    } // END for y
                } // END for x
            } // END IF leaf-internal node
        }

//        void calculate_operation_zonal_op2(size_type base_row, size_type base_col,
//                                          ushort level, size_type sub_size,
//                                          k2_raster &k2_raster_1, k2_raster &raster_zonal,
//                                          k2raster_submatrix<value_type> &node1, k2raster_submatrix<value_type> &nodeZ,
//                                          std::vector<size_type> &last_access1, std::vector<size_type> &last_access2,
//                                          const OperationZonalRaster &operation, std::map<value_type, value_type> &zonal_sum,
//                                          std::vector<std::vector<std::tuple<bool, value_type, value_type>>> &queue_nodes,
//                                          std::vector<std::vector<size_type>> &parent_nodes) {
//
//            // ------------------------------------------------------------------- //
//            // If both nodes are uniforms, the new node is also uniform
//            // ------------------------------------------------------------------- //
//            if (node1.is_uniform && nodeZ.is_uniform) {
//                size_type n_nodes = 1;
//                if (sub_size != 1) {
//                    size_t c_x = std::min(sub_size, this->k_real_size_x - base_row);
//                    size_t c_y = std::min(sub_size, this->k_real_size_y - base_col);
//                    n_nodes = c_x * c_y;
//                }
//                // Run operation
//                switch (operation) {
//                    case OPERATION_ZONAL_SUM:
//                        if (zonal_sum.find(nodeZ.max_value) == zonal_sum.end()) {
//                            zonal_sum[nodeZ.max_value] = node1.max_value * n_nodes;
//                        } else {
//                            zonal_sum[nodeZ.max_value] += node1.max_value * n_nodes;
//                        }
//                        break;
//                    default:
//                        // No valid operation. The process never reaches this point
//                        exit(-1);
//                }
//            } else {
//
//                // Calculate parameters for the current level.
//                ushort k = this->get_k(level);
//                sub_size = sub_size / k;
//
//                /*****************/
//                /* INTERNAL NODE */
//                /*****************/
//                k2raster_submatrix<value_type> child1, childZ;
//                child1.level = node1.level + 1;
//                childZ.level = nodeZ.level + 1;
//                size_type child_base_row, child_base_col;
//
//                for (uint x = 0; x < k; x++) {
//                    child_base_row = base_row + x * sub_size;
//                    for (uint y = 0; y < k; y++) {
//                        child_base_col = base_col + y * sub_size;
//
//                        if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {
//
//                            k2_raster_1.get_next_child(node1, child1, last_access1);
//
//                            if (!nodeZ.is_uniform) {
//                                raster_zonal.get_next_child(nodeZ, childZ, last_access2);
//                                queue_nodes[nodeZ.level].push_back({childZ.is_uniform, childZ.max_value, childZ.max_value});
//                                if (!childZ.is_uniform){
//                                    size_type pos = nodeZ.children_pos - raster_zonal.k_accum_max_values[nodeZ.level];
//                                    parent_nodes[nodeZ.level].push_back(pos); // TODO change this
//                                }
//                            } else {
//                                childZ = nodeZ;
//                            }
//
//                            calculate_operation_zonal_op2(child_base_row, child_base_col, level + 1, sub_size,
//                                                         k2_raster_1, raster_zonal, child1, childZ,
//                                                         last_access1, last_access2, operation,
//                                                         zonal_sum, queue_nodes, parent_nodes);
//                        } else {
//                            if (!nodeZ.is_uniform) {
//                                // Add a fake node (out of real raster)
//                                queue_nodes[nodeZ.level].push_back({true,
//                                                                    std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()});
//                            }
//                            // More children in raster1? Reset last_access1
//                            if (child_base_row < k2_raster_1.k_real_size_x && child_base_col < k2_raster_1.k_real_size_y) {
//                                for (size_t l = child1.level-1; l < last_access1.size(); l++) {
//                                    last_access1[l] = 0;
//                                }
//                            }
//
//                            // More children in raster2? Reset last_access2
//                            if (child_base_row < raster_zonal.k_real_size_x && child_base_col < raster_zonal.k_real_size_y) {
//                                for (size_t l = childZ.level-1; l < last_access2.size(); l++) {
//                                    last_access2[l] = 0;
//                                }
//                            }
//                        }// ENF IF within real matrix
//
//                        node1.children_pos++;
//                        nodeZ.children_pos++;
//                    } // END for y
//                } // END for x
//            } // END IF leaf-internal node
//        }

        template <typename value_in_type=value_type, class raster_type>
        void create_zonal_raster(std::vector<std::vector<value_type>> &min_values_, std::vector<std::vector<value_type>> &max_values_, std::vector<sdsl::int_vector<1>> &tmp_t_,
                                       value_type &min_value, value_type &max_value, size_type base_row, size_type base_col,
                                       ushort level, size_type sub_size,
                                       raster_type &k2_raster_z, k2raster_submatrix<value_in_type> &node_z,
                                       std::vector<size_type> &last_access_z, std::map<value_in_type, value_type> &zonal_sum) {

            min_value = std::numeric_limits<value_type>::max();
            max_value = std::numeric_limits<value_type>::min();

            // Calculate parameters for the current level.
            ushort k = this->get_k(level);
            sub_size = sub_size / k;

            // ------------------------------------------------------------------- //
            // If both nodes are uniforms, the new node is also uniform
            // ------------------------------------------------------------------- //
            if (node_z.is_uniform) {
                //std::cout << "zone: " << node_z.max_value << " = " << zonal_sum[node_z.max_value] << std::endl;
                //exit(-1);
                max_value = zonal_sum[node_z.max_value];
                min_value = max_value; // Uniform node
            } else {
                /*****************/
                /* INTERNAL NODE */
                /*****************/
                std::vector<value_type> min_values_children(k * k);
                std::vector<value_type> max_values_children(k * k);
                k2raster_submatrix<value_in_type> child_z;
                child_z.level = node_z.level + 1;

                size_type child_base_row, child_base_col;

                for (uint x = 0; x < k; x++) {
                    child_base_row = base_row + x * sub_size;
                    for (uint y = 0; y < k; y++) {
                        child_base_col = base_col + y * sub_size;

                        if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {

                            k2_raster_z.get_next_child(node_z, child_z, last_access_z);
                            create_zonal_raster<value_in_type>(min_values_, max_values_, tmp_t_,
                                                min_values_children[x * k + y], max_values_children[x * k + y],
                                                child_base_row, child_base_col, level + 1, sub_size,
                                                k2_raster_z, child_z, last_access_z, zonal_sum);

                            // Calculate min/max value of the parent node (current node)
                            if (min_values_children[x * k + y] < min_value) {
                                min_value = min_values_children[x * k + y];
                            }
                            if (max_values_children[x * k + y] > max_value) {
                                max_value = max_values_children[x * k + y];
                            }
                        } else {
                            // TODO remove this?
                            // More children in raster1? Reset last_access1
                            if (child_base_row < k2_raster_z.m_real_size_x && child_base_col < k2_raster_z.m_real_size_y) {
                                for (size_t l = child_z.level-1; l < last_access_z.size(); l++) {
                                    last_access_z[l] = 0;
                                }
                            }

                            // Positions out of real matrix
                            min_values_children[x * k + y] = std::numeric_limits<value_type>::max();
                            max_values_children[x * k + y] = std::numeric_limits<value_type>::min();
                        }// ENF IF within real matrix

                        node_z.children_pos++;
                    } // END for y
                } // END for x

                // Check if matrix is not empty or uniform
                if (min_value < max_value) {
                    insert_node_conceptual_tree(min_values_, max_values_, tmp_t_,
                                                min_values_children, max_values_children,
                                                min_value, max_value, level, k);

                } // END IF (min_value < max_value)
            } // END IF leaf-internal node
        }

        template <typename value_in_type=value_type>
        void create_zonal_raster_op(std::vector<std::vector<value_type>> &min_values_, std::vector<std::vector<value_type>> &max_values_,
                                    std::vector<sdsl::int_vector<1>> &tmp_t_,
                                    std::map<value_in_type, value_type> &zonal_sum,
                                    std::vector<std::vector<std::pair<uint8_t, value_type>>> &queue_nodes,
                                    std::vector<std::vector<size_type>> &parent_nodes) {

            value_type min_value, max_value;
            value_type max_value_node;
            std::vector<value_type> min_values_children;
            std::vector<value_type> max_values_children;
            std::queue<value_type> min_values_parents;

            for (short l = this->k_height-1; l >= 0; l--) {
                ushort k = this->get_k(l);
                min_values_children.resize(k * k);
                max_values_children.resize(k * k);
                size_type n_node = 0;
                size_type n_parent = 0;

                for (size_type n = 0; n < queue_nodes[l].size(); n += (k * k)) {
                    min_value = std::numeric_limits<value_type>::max();
                    max_value = std::numeric_limits<value_type>::min();

                    for (ushort c = 0; c < k * k; c++) {
                        //Read the next k nodes in the queue (they are children of the same parent)
                        auto node = queue_nodes[l][n_node++]; // Get the next child of the current node
                        max_value_node = node.second;

                        switch (node.first) {
                            case 0:
                                // Interval node
                                min_values_children[c] = min_values_parents.front();
                                min_values_parents.pop();
                                max_values_children[c] = max_value_node;
                                break;
                            case 1:
                                // Leaf (uniform) node
                                min_values_children[c] = zonal_sum[max_value_node];
                                max_values_children[c] = zonal_sum[max_value_node];
                                break;
                            case 2:
                                // Node out of real bounds (set min > max value)
                                min_values_children[c] = 1;
                                max_values_children[c] = 0;
                                continue;
                        }

                        // Calculate min/max value of the parent node (current node)
                        min_value = std::min(min_value, min_values_children[c]);
                        max_value = std::max(max_value, max_values_children[c]);
                    } // END FOR c

                    if (min_value < max_value) {
                        insert_node_conceptual_tree(min_values_, max_values_, tmp_t_,
                                                    min_values_children, max_values_children,
                                                    min_value, max_value, l, k);

                    } // END IF (min_value < max_value)

                    // Update father
                    if (l != 0) {
                        auto &parent = queue_nodes[l - 1][parent_nodes[l - 1][n_parent++]];
                        parent.second = max_value;
                        min_values_parents.push(min_value);
                    } else {
                        // Update root node
                        this->k_min_value = min_value;
                        this->k_max_value = max_value;
                    }
                } // END WHILE queue

                // Free memory
                if (l != 0) {
                    parent_nodes[l - 1].clear();
                    parent_nodes[l - 1].shrink_to_fit();
                }
                queue_nodes[l].clear();
                queue_nodes[l].shrink_to_fit();
            } // END FOR level l
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
                            if (this->is_plain_level(level)) {
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
                        for (auto dx = c_xini; dx <= c_xend; dx++) {
                            cell_pos = (dx - or_x) * window_size + (c_yini - or_y);
                            for (auto dy = c_yini; dy <= c_yend; dy++) {
                                result[cell_pos++] = max_value;
                            } // END FOR y
                        } // END FOR x
                        count_cells += (c_xend-c_xini+1) * (c_yend-c_yini+1);
                        continue; // Go to next node
                    } else {
                        if (this->is_plain_level(level)) {
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
                            if (this->is_plain_level(level)) {
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

        //*****************************//
        //***** DECOMPRESS HELPER *****//
        //*****************************//
        template <class Container>
        size_type decompress_helper(size_type base_x, size_type base_y, value_type f_max_value,
                                           size_type children_pos, size_type children_size, ushort level,
                                            Container&& result) {
            size_type pos;
            size_type c_base_x, c_base_y;
            value_type max_value;

            ushort k = this->get_k(level-1);                        // Current value of "k"
            bool is_uniform, is_leaf = level == this->k_height;     // True -> They are leaves
            size_type count_cells = 0;

            // Check children
            for (auto x = 0; x < k; x++) {
                c_base_x = base_x + x * children_size;                                                  // Calculate base position of the current child (row)
                if (c_base_x >= this->k_real_size_x) continue;
                for (auto y = 0; y < k; y++) {
                    pos = children_pos + x * k + y;                                                     // Get position at Tree (T)
                    c_base_y = base_y + y * children_size;                                              // Calculate base position of the current child (column)
                    if (c_base_y >= this->k_real_size_y) continue;

                    max_value = f_max_value - this->get_max_value_op(level, pos);
//                            this->k_list_max[level-1].access(pos - this->k_accum_max_values[level-1]);  // Get max value
                    // Check if it is a uniform submatrix

                    is_uniform = is_leaf || !this->k_t[pos];
                    if (is_uniform) {

                        // Calculate positions with overlap with the interesting region
                        auto c_xend = std::min(c_base_x + children_size - 1, this->k_real_size_x-1);
                        auto c_yend = std::min(c_base_y + children_size - 1, this->k_real_size_y-1);

                        size_type cell_pos;
                        for (auto dx = c_base_x; dx <= c_xend; dx++) {
                            cell_pos = dx * this->k_real_size_y + c_base_y;
                            for (auto dy = c_base_y; dy <= c_yend; dy++) {
                                result[cell_pos++] = max_value;
                            } // END FOR y
                        } // END FOR x
                        count_cells += (c_xend-c_base_x+1) * (c_yend-c_base_y+1);
                        continue; // Go to next node
                    } else {
                        if (this->is_plain_level(level)) {
                            // Child is represented with plain values
                            size_type ones = (this->k_t_rank1(pos) + 1) - (this->k_accum_min_values[level-1]+1);
                            //count_cells += this->decompress_plain(children_size, ones * children_size* children_size,
                            //                                             max_value, result);
                        } else {
                            // Go down one level on the tree
                            size_type new_children_pos = this->get_children_position(pos, level);
                            count_cells += decompress_helper(c_base_x, c_base_y, max_value, new_children_pos,
                                                             children_size / this->get_k(level), level + 1, result);
                        }
                    } // END IF is_uniform
                } // END FOR y
            } // END FOR x
            return count_cells;
        }




    };
} // END NAMESPACE k2raster

#endif // INCLUDED_K2_RASTER_
