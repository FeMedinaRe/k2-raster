/*
 * Created by Fernando Silva on 3/07/18.
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

#ifndef INCLUDED_K2_RASTER_PLAIN
#define INCLUDED_K2_RASTER_PLAIN

#include <sdsl/vectors.hpp>
#include <k2_raster.hpp>


//! Namespace for the k2-raster library
namespace k2raster {

    template <typename t_value=int,
            typename t_bv=sdsl::bit_vector,
            typename t_rank=sdsl::rank_support_v5<1,1>,
            typename t_values_vec=sdsl::dac_vector_dp_opt<t_bv, t_rank, 3>,
            typename t_values_last=sdsl::dac_vector_dp_opt<t_bv, t_rank, 3>>
    class k2_raster_plain : public k2_raster<t_value, t_bv, t_rank, t_values_vec>
    {

    public:
        typedef k2_raster<t_value, t_bv, t_rank, t_values_vec>  k2_raster_p; // k2_raster_parent
        typedef t_value                                         value_type;
        typedef k2_raster<>::size_type                          size_type;
        typedef t_bv                                            bit_vector_type;

    protected:
        ushort k_level_k2=0;
        ushort k_level_plain=0;

        ushort k_size_leaves=0;
        t_values_last k_plain_values;

    public:
        const ushort &m_level_plain = k_level_plain;
        const ushort &m_size_leaves = k_size_leaves;
        const t_values_last &m_plain_values = k_plain_values;

    public:

        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//
        k2_raster_plain() = default;

        k2_raster_plain(const k2_raster_plain& tr) : k2_raster_p()
        {
            *this = tr;
        }

        k2_raster_plain(k2_raster_plain&& tr)
        {
            *this = std::move(tr);
        }

        template<class Container>
        k2_raster_plain(Container &&c_values, size_type n_rows, size_type n_cols, ushort k1, ushort k2, ushort level_k1, ushort plains_levels)
                : k2_raster_p(n_rows, n_cols, k1, k2, level_k1, K2_RASTER_TYPE_PLAIN) {

            k_level_plain = plains_levels;
            this->init_levels();

            std::vector<value_type> plain_values_;
            build(c_values, plain_values_);
            encoded_plain_values(plain_values_);
        }

        //*******************************************************//
        //*************** BASIC OPERATIONS **********************//
        //*******************************************************//

        //! Move assignment operator
        k2_raster_plain& operator=(k2_raster_plain&& tr)
        {
            if (this != &tr) {
                k2_raster_p::operator=(tr);

                k_level_k2 = std::move(tr.k_level_k2);
                k_level_plain = std::move(tr.k_level_plain);
                k_size_leaves = std::move(tr.k_size_leaves);
                k_plain_values = std::move(tr.k_plain_values);
            }
            return *this;
        }

        //! Assignment operator
        k2_raster_plain& operator=(const k2_raster_plain& tr)
        {
            if (this != &tr) {
                k2_raster_p::operator=(tr);

                k_level_k2 = tr.k_level_k2;
                k_level_plain = tr.k_level_plain;
                k_size_leaves = tr.k_size_leaves;
                k_plain_values = tr.k_plain_values;
            }
            return *this;
        }

        //! Swap operator
        void swap(k2_raster_plain& tr)
        {
            if (this != &tr) {
                k2_raster_p::swap(tr);

                std::swap(k_level_k2, tr.k_level_k2);
                std::swap(k_level_plain, tr.k_level_plain);
                std::swap(k_size_leaves, tr.k_size_leaves);
                k_plain_values.swap(tr.k_plain_values);
            }
        }

        //! Equal operator
        bool operator==(const k2_raster_plain& tr) const
        {
            if (!k2_raster_p::operator==(tr)) {
                return false;
            }

            if (k_level_k2 != tr.k_level_k2 || k_level_plain != tr.k_level_plain ||
                    k_size_leaves != tr.k_size_leaves) {
                return false;
            }

            if (k_plain_values.size() != tr.k_plain_values.size()) {
                return false;
            }
            return true;
        }

        //*******************************************************//
        //******************** QUERIES **************************//
        //*******************************************************//

        value_type get_cell(size_type  row, size_type col) const {

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
                ushort levels = std::min(this->k_level_k1 + this->k_level_k2+1, this->k_height-1);
                ushort k = this->get_k(0);
                size_type child;

                for (auto l = 1; l <= levels; l++) {
                    size /= k;                                  // Submatriz size
                    child = (row / size) * k + (col / size);    // Number of child (following a z-order)
                    child_pos += child;                         // Position in bitmap k_t
                    max_value -=  this->get_max_value(l, child_pos);
                            //this->k_list_max[l - 1][child_pos - this->k_accum_max_values[l-1]];

                    // Check if all the cells in the subarray are equal (uniform matrix)
                    if (!this->k_t[child_pos]) {
                        return max_value;
                    }

                    // Go down one level on the tree.
                    row = row % size;                           // Update local row
                    col = col % size;                           // Update local column

                    if (this->k_level_plain != 0 && l == levels) {
                        break;
                    }

                    child_pos = this->get_children_position(child_pos, l);
                    k = this->get_k(l);
                }
            } // END BLOCK searching from level 1 to level l - 1


            // ************************************//
            // Searching last level                //
            // ************************************//
            {
                if (this->k_level_plain == 0) {
                    // No plain values, search last level
                    ushort k = this->get_k(this->k_height-1);
                    child_pos += row * k + col;
                    max_value -= this->get_max_value(this->k_height, child_pos);
                            //this->k_list_max[this->k_height - 1][child_pos - this->k_accum_max_values[this->k_height - 1]];
                    return max_value;
                } else {
                    // Value is in plain form
                    size_type ones = this->k_t_rank1(child_pos) + 1;      // Number of non-empty nodes until position 'child_pos'
                    ones = ones - (this->k_accum_min_values[this->k_height - 2] + 1);
                    child_pos = ones * size * size;
                    return get_cell_plain(row, col, size, child_pos, max_value);
                }
            } // END BLOCK last levels
        }

        //*****************************//
        //***** GET CELL BY VALUE *****//
        //*****************************//

        virtual size_type get_cells_by_value_plain(size_type xini, size_type xend, size_type yini, size_type yend,
                                                   value_type valini, value_type valend,
                                                   size_type size, size_type children_pos, value_type father_value,
                                                   std::vector<std::pair<size_type, size_type>> &result) {
            value_type value;
            size_type count_cells = 0;

            for (auto x = xini; x <= xend; x++) {
                size_type pos = children_pos + ((x % size) * size + (yini % size));
                for (auto y = yini; y <= yend; y++) {
                    value = father_value - k_plain_values[pos];
                    if (value >= valini && value <= valend) {
                        result.emplace_back(x, y);
                        count_cells++;
                    }
                    pos++; // Position of the next value
                } // END FOR y
            } // END FOR x
            return count_cells;
        }

        //*****************************//
        //***** GET VALUES WINDOW *****//
        //*****************************//
        virtual size_type get_values_window_plain(size_type xini, size_type xend, size_type yini, size_type yend,
                                                  size_type size, size_type children_pos, value_type father_value,
                                                  std::vector<value_type> &result, size_type or_x, size_type or_y, size_type window_size) {
            size_type count_cells = 0, cell_pos;
            for (auto x = xini; x <= xend; x++) {
                size_type pos = children_pos + ((x % size) * size + (yini % size));
                cell_pos = (x - or_x) * window_size + (yini - or_y);
                for (auto y = yini; y <= yend; y++) {
                    result[cell_pos++] = father_value - k_plain_values[pos];
                    count_cells++;
                    pos++; // Position of the next value
                } // END FOR y
            } // END FOR x
            return count_cells;
        }

        //*****************************//
        //**** CHECK VALUES WINDOW ****//
        //*****************************//
        virtual bool check_values_window_plain(size_type xini, size_type xend, size_type yini, size_type yend,
                                               value_type valini, value_type valend,
                                               size_type size, size_type children_pos, value_type father_value,
                                               bool strong_check) {
            value_type value;
            for (auto x = xini; x <= xend; x++) {
                size_type pos = children_pos + ((x % size) * size + (yini % size));
                for (auto y = yini; y <= yend; y++) {
                    value = father_value - k_plain_values[pos];
                    if (value >= valini && value <= valend) {
                        if (!strong_check) return true;
                    } else {
                        if (strong_check) return false;
                    }
                    pos++; // Position of the next value
                } // END FOR y
            } // END FOR x
            return strong_check;
        }

        //*******************************************************//
        //******************** AUXILIARY ************************//
        //*******************************************************//
        virtual inline bool is_plain_level(ushort level) const {
            return this->k_level_plain > 0 && (level == (this->k_level_k1 + this->k_level_k2));
        }

        virtual short get_cell_n_levels() const {
            return std::min(this->k_level_k1 + this->k_level_k2+1, this->k_height-1);
        }

        //*******************************************************//
        //********************** FILE ***************************//
        //*******************************************************//
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {

            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;

            written_bytes += serialize_base(out, child, name);
            sdsl::serialize(k_plain_values, out, child, "plain_values");

            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            load_base(in);
            sdsl::load(k_plain_values, in);;
        }

        void copy(ushort &level_k2, ushort &level_plain, ushort &size_leaves) const {
            level_k2 = k_level_k2;
            level_plain = k_level_plain;
            size_leaves = k_size_leaves;
        }


        //*******************************************************//
        //********************** INFO ***************************//
        //*******************************************************//
        virtual void print_info() {
            print_info_base();
            std::cout << "Plain values:";
            for (auto v = 0; v < k_plain_values.size(); v++) {
                value_type value = k_plain_values[v];
                std::cout << value << " ";
            }
            std::cout << "\t(" << k_plain_values.size() << " values)" << std::endl;
        }

        //*******************************************************//
        //******************** NAVIGATION ***********************//
        //*******************************************************//
        virtual k2raster_node<value_type> get_child(const k2raster_node<value_type> parent, uint child) const {
            if (parent.level > this->k_level_k1 + k_level_k2) {
                return parent;
            }
            k2raster_node<value_type> node = k2_raster_p::get_child(parent, child);

            if (node.level == this->k_level_k1 + k_level_k2+1) {
                node.children_pos = parent.children_pos + child;
            }
            return node;
        }

        virtual inline void set_children_position(k2raster_node<value_type> &child_node, size_type ones ) const {
            if (child_node.level <= this->k_level_k1 + k_level_k2) {
                child_node.children_pos = this->get_children_position_ones(ones, child_node.level);
            }
        }

        inline size_type get_plain_value(size_type row, size_type col, value_type max_value,
                                         size_type child_pos, size_type size) const {
            size_type ones = this->k_t_rank1(child_pos);      // Number of non-empty nodes until position 'child_pos'
            ones = ones - this->k_accum_min_values[this->k_list_min.size()-1];
            child_pos = ones * k_size_leaves;
            return get_cell_plain(row, col, size, child_pos, max_value);
        }

        virtual std::vector<value_type> get_plain_values(const k2raster_node<value_type> &parent) const {
            if (parent.min_value == parent.max_value) {
                return std::vector<value_type>(k_size_leaves, parent.max_value);
            }

            size_type ones = this->k_t_rank1(parent.children_pos) + 1;      // Number of non-empty nodes until position 'child_pos'
            ones = ones - (this->k_accum_min_values[parent.level-1] + 1);   // TODO check parent.level == 0
            size_type child_pos = ones * k_size_leaves;                     // Start position

            std::vector<value_type> values (k_size_leaves);
            for (auto c = 0; c < k_size_leaves; c++) {
                values[c] = parent.max_value - k_plain_values[child_pos + c];
            }
            return values;
        }

        //*******************************************************//
        //********************** HELPERS ************************//
        //*******************************************************//

        void copy_parameters(ushort &level_k1, ushort &level_k2, ushort &level_plain, ushort &size_leaves) {
            level_k1 = this->k_level_k1;
            level_k2 = k_level_k2;
            level_plain = k_level_plain;
            size_leaves = k_size_leaves;
        }

    protected:

        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//
        template<class Container>
        k2_raster_plain(Container &&c_values, size_type n_rows, size_type n_cols, ushort k1, ushort k2, ushort level_k1)
                : k2_raster_p(c_values, n_rows, n_cols, k1, k2, level_k1) {
        }


        //*******************************************************//
        //******************** AUXILIARY ************************//
        //*******************************************************//
        void init_levels(){
            k2_raster_p::init_levels();


            // Recalculate k_level_k1 and k_level_k2
            if (this->k_height - k_level_plain >= this->k_level_k1) {
                k_level_k2 = this->k_height - this->k_level_k1 - k_level_plain;
            } else {
                k_level_k2 = 0;
                this->k_level_k1 = this->k_height - k_level_plain;
            }
            if (this->k_level_k1 == 0 && k_level_k2 == 0) {
                std::cout << "Error: k1 = " << this->k_k1 << ", k2 = " << this->k_k2 << " and  Level plain: " << k_level_plain << " are not valid params";
                std::cout << " with size " << this->k_real_size_x << "x" << this->k_real_size_x << " (" << this->k_size << "x" << this->k_size << ")" << std::endl;
                exit(-1);
            }

            k_size_leaves = this->k_size;
            for (auto l = 0; l <= this->k_height - k_level_plain; l++) {
                k_size_leaves /= this->get_k(l);
            }
            k_size_leaves *= k_size_leaves;
        }

        //*******************************************************//
        //********************** BUILD **************************//
        //*******************************************************//

        template<class Container>
        void build(Container &&c_values, std::vector<value_type> &plain_values_) {
            std::vector<std::vector<value_type>> min_values_(this->k_height);
            std::vector<std::vector<value_type>> max_values_(this->k_height);
            std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);

            build(c_values, min_values_, max_values_, tmp_t_, plain_values_, this->k_min_value, this->k_max_value, 0, 0, 0, this->k_size);

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
                        std::vector<value_type> &plain_values_, value_type &min_value, value_type &max_value, size_type base_row, size_type base_col,
                        ushort level, size_type sub_size) {

            size_type pos_value, child_base_row, child_base_col;
            value_type value;
            min_value = std::numeric_limits<value_type>::max();
            max_value = std::numeric_limits<value_type>::min();


            if (level == this->k_level_k1 + k_level_k2 + 1) {
                /****************/
                /* PLAIN VALUES */
                /****************/
                // Seek the max and min value
                for (size_type r = base_row; r < base_row + sub_size; r++) {
                    for (size_type c = base_col; c < base_col + sub_size; c++) {
                        if (r < this->k_real_size_x && c < this->k_real_size_y) {
                            pos_value = r * this->k_real_size_y + c;
                            value = c_values[pos_value];
                            if (value < min_value) {
                                min_value = value;
                            }
                            if (value > max_value) {
                                max_value = value;
                            }
                        }
                    } // END for c
                } // END for r

                if (min_value != max_value){
                    // Add values to plain_values array
                    for (size_type r = base_row; r < base_row + sub_size; r++) {
                        for (size_type c = base_col; c < base_col + sub_size; c++) {
                            if (r < this->k_real_size_x && c < this->k_real_size_y) {
                                pos_value = r * this->k_real_size_y + c;
                                plain_values_.push_back(max_value - c_values[pos_value]);
                            } else {
                                // Virtual cell
                                plain_values_.push_back(0);
                            }// END IF
                        } // END FOR col
                    } // END FOR row
                }
            } else {
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
                                build(c_values, min_values_, max_values_, tmp_t_, plain_values_,
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
                                if (level != (this->k_height-1)) {
                                    tmp_t_[level][max_values_[level].size()] = 0;
                                }
                                max_values_[level].push_back(0);
                            } else {
                                max_values_[level].push_back(max_value - max_values_children[c]);
                                if (level != (this->k_height-1)) {                                          // Except last level
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
                } // END if (level == (this->k_height-1))
            } // END IF plain/normal partition
        }

        virtual void encoded_plain_values(std::vector<value_type> plain_values_) {
            k_plain_values = t_values_last(plain_values_);
        }

        //*******************************************************//
        //**************** QUERIES HELPERS **********************//
        //*******************************************************//

        //********************//
        //***** GET CELL *****//
        //********************//
        virtual value_type get_cell_plain(size_type row, size_type col, size_type size,
                                          size_type children_pos, value_type father_value) const {
            row %= size;
            col %= size;
            size_type pos = children_pos + (row * size + col);
            return father_value - k_plain_values[pos];
        }


        //*******************************************************//
        //********************** INFO ***************************//
        //*******************************************************//
        virtual void print_info_base() {
            std::cout << "Size " << this->k_size << "x" << this->k_size << " (" << this->k_real_size_x << "x" << this->k_real_size_y << ")" << std::endl;
            std::cout << "Levels: " << this->k_height << " || " << this->k_level_k1 << " (k1) + " << this->k_level_k2 << " (k2) + " << k_level_plain << " (plain)" << std::endl;

            std::cout << "Tree (T): ";
            for (auto b = 0; b < this->k_t.size(); b++) {
                std::cout << this->k_t[b];
            }
            std::cout << "(\t" << this->k_t.size() << " bits)" << std::endl;

            std::cout << std::endl;
            std::cout << "rMax: " << this->k_max_value << " || rMin: " << this->k_min_value << std::endl;

            std::cout << std::endl;
            std::cout << "LMax:" << std::endl;
            for (auto l = 0; l < this->k_list_max.size(); l++) {
                std::cout << "L" << l+1 << ": ";
                for (auto v = 0; v < this->k_list_max[l].size(); v++) {
                    value_type value = this->k_list_max[l][v];
                    std::cout << value << " ";
                }
                std::cout << "\t(" << this->k_list_max[l].size() << " values)" << std::endl;
            }

            std::cout << std::endl;
            std::cout << "LMin:" << std::endl;
            for (auto l = 0; l < this->k_list_min.size(); l++) {
                std::cout << "L" << l+1 << ": ";
                for (auto v = 0; v < this->k_list_min[l].size(); v++) {
                    value_type value = this->k_list_min[l][v];
                    std::cout << value << " ";
                }
                std::cout << "\t(" << this->k_list_min[l].size() << " values)" << std::endl;
            }

            std::cout << std::endl;
        }

        //*******************************************************//
        //********************** FILE ***************************//
        //*******************************************************//
        size_type serialize_base(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {

            size_type written_bytes = 0;

            // k2-raster base
            written_bytes += k2_raster_p::serialize(out, v, name);

            // Plain parameters
            written_bytes += write_member(k_level_k2, out, v, "level_k2");
            written_bytes += write_member(k_level_plain, out, v, "level_plain");
            written_bytes += write_member(k_size_leaves, out, v, "k_size_leaves");
            return written_bytes;
        }

        void load_base(std::istream& in) {

            // k2-raster base
            k2_raster_p::load(in);

            // Plain parameters
            sdsl::read_member(k_level_k2, in);
            sdsl::read_member(k_level_plain, in);
            sdsl::read_member(k_size_leaves, in);
        }

    }; // END class k2_raster_plain
} // END namespace sdsl

#endif // INCLUDED_K2_RASTER_PLAIN
