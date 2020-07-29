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

#ifndef INCLUDED_T_K2_RASTER_PLAIN_LOG
#define INCLUDED_T_K2_RASTER_PLAIN_LOG

#include <k2_raster_plain.hpp>
#include <sdsl/dac_vector.hpp>


//! Namespace for k2-raster library
namespace k2raster {

    //*******************************************************//
    //******************* NODE HELPER ***********************//
    //*******************************************************//

    template<typename t_value=int,
            typename t_bv=sdsl::bit_vector,
            typename t_rank=sdsl::rank_support_v5<1,1>,
            typename t_values_vec=sdsl::dac_vector_dp<>,
            typename t_values_last=sdsl::dac_vector_dp_opt<t_bv, t_rank, 3>>
    class t_k2_raster_plain_log : public t_k2_raster_log<t_value, t_bv, t_rank, t_values_vec> {

    public:
        typedef t_k2_raster_log<t_value, t_bv, t_rank, t_values_vec>  k2_raster_p; // k2_raster parent
        typedef t_k2_raster_log<>::value_type                       value_type;
        typedef t_k2_raster_log<>::size_type                        size_type;
        typedef t_k2_raster_log<>::bit_vector_type                  bit_vector_type;
        typedef t_k2_raster_log<>::encode_type                      encode_type;
        typedef t_k2_raster_log<>::decode_type                      decode_type;

    protected:
        ushort k_level_k2=0;
        ushort k_level_plain=0;

        ushort k_size_leaves=0;
        t_values_last k_plain_values;

    public:

    //*******************************************************//
    //******************* CONSTRUCTOR ***********************//
    //*******************************************************//
    t_k2_raster_plain_log() = default;

    t_k2_raster_plain_log(const t_k2_raster_plain_log &tr) : k2_raster_p() {
        *this = tr;
    }

    t_k2_raster_plain_log(t_k2_raster_plain_log &&tr) {
        *this = std::move(tr);
    }

    template<class Snapshot_type>
    t_k2_raster_plain_log(Snapshot_type &&raster_snapshot, ushort k2_raster_type)
            : k2_raster_p(raster_snapshot, k2_raster_type) {
        raster_snapshot.copy_parameters(this->k_level_k1, k_level_k2, k_level_plain, k_size_leaves);
    }

    template<class Container, class Snapshot_type>
    t_k2_raster_plain_log(Container &&c_values, Snapshot_type &&raster_snapshot)
            : k2_raster_p(raster_snapshot, K2_RASTER_TYPE) {
        raster_snapshot.copy_parameters(this->k_level_k1, k_level_k2, k_level_plain, k_size_leaves);;

        // Build k2-raster
        std::vector<value_type> plain_values_;
        build(c_values, raster_snapshot, plain_values_);
        encoded_plain_values(plain_values_);
    }

    //*******************************************************//
    //*************** BASIC OPERATIONS **********************//
    //*******************************************************//

    //! Move assignment operator
    t_k2_raster_plain_log& operator=(t_k2_raster_plain_log&& tr)
    {
        if (this != &tr) {
            k2_raster_p::operator=(tr);

            /****** Structures Plain ******/
            k_level_k2 = std::move(tr.k_level_k2);
            k_level_plain = std::move(tr.k_level_plain);
            k_size_leaves = std::move(tr.k_size_leaves);
            k_plain_values = std::move(tr.k_plain_values);
        }
        return *this;
    }

    //! Assignment operator
    t_k2_raster_plain_log& operator=(const t_k2_raster_plain_log& tr)
    {
        if (this != &tr) {
            k2_raster_p::operator=(tr);

            /****** Structures Plain ******/
            k_level_k2 = tr.k_level_k2;
            k_level_plain = tr.k_level_plain;
            k_size_leaves = tr.k_size_leaves;
            k_plain_values = tr.k_plain_values;
        }
        return *this;
    }

    //! Swap operator
    void swap(t_k2_raster_plain_log& tr)
    {
        if (this != &tr) {
            k2_raster_p::swap(tr);

            /****** Structures Plain ******/
            std::swap(k_level_k2, tr.k_level_k2);
            std::swap(k_level_plain, tr.k_level_plain);
            std::swap(k_size_leaves, tr.k_size_leaves);
            k_plain_values.swap(tr.k_plain_values);
        }
    }

    //! Equal operator
    bool operator==(const t_k2_raster_plain_log& tr) const
    {
        if (!k2_raster_p::operator==(tr)) {
            return false;
        }
        if (k_level_k2 != tr.k_level_k2 || k_level_plain != tr.k_level_plain
        || k_size_leaves != tr.k_size_leaves) {
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

    //*****************************//
    //***** GET CELL          *****//
    //*****************************//
    template<class k2_raster_snap>
    value_type get_cell(const k2_raster_snap &raster_snap, size_type row, size_type col) const {
        // Cell outside the matrix
        if (row > this->k_real_size_x || col > this->k_real_size_y) {
            return 0;
        }

        bool has_children = true;                           // The current nodes has not children
        bool is_uniform = false;                            // The current node is uniform

        // Check first level
        if (this->k_t.empty()) {
            if (this->k_min_value == this->k_max_value) {
                // All cell of the matrix are equal
                return this->k_max_value;
            }
            has_children = false;
        }

        // ************************************//
        // Searching from level 1 to level l-1 //
        // ************************************//
        size_type size = this->k_size;                                      // Current matrix size
        size_type pos_log = 0;                                              // Children position on bitmap k_t of the current node
        decode_type diff_max = this->k_max_value - raster_snap.max_value;   // Difference between the current raster and the snapshot

        // ************************************//
        // Snapshot attributes                 //
        // ************************************//
        size_type pos_snap = 0;                                             // Children position on bitmap k_t of the current node (snapshot)
        value_type max_snap = raster_snap.max_value;                        // Current max value (snapshot)
        bool is_uniform_snap = raster_snap.max_value == raster_snap.min_value; // The current node of the snapshot is uniform

        ushort levels = get_cell_n_levels();
        {
            ushort k = this->get_k(0);
            size_type child;

            for (auto l = 1; l <= levels; l++) {
                size /= k;                                      // Submatrix size
                child = (row / size) * k + (col / size);        // Number of child (following a z-order)
                pos_log += has_children ? child : 0;          // Position in bitmap k_t
                pos_snap += !is_uniform_snap ? child : 0; // Position in bitmap k_t (snapshot)

                // Check snapshot node
                if (!is_uniform_snap) {
                    is_uniform_snap = !raster_snap.t[pos_snap];
                    max_snap -= raster_snap.get_max_value(l, pos_snap);
                            //raster_snap.list_max[l - 1][pos_snap - raster_snap.m_accum_max_values[l-1]];
                }

                // Check current raster node
                if (has_children) {
                    has_children = this->k_t[pos_log];
                    if (!has_children) {
                        is_uniform = !this->k_eqB[pos_log - this->k_t_rank1(pos_log)];
                        diff_max = this->__decode(this->get_max_value(l, pos_log));
//                        diff_max = this->__decode(this->k_list_max[l-1][pos_log - this->k_accum_max_values[l-1]]);
                    }
                }

                // Return value
                if ((is_uniform_snap && !has_children) || is_uniform) {
                    return max_snap + diff_max;
                }

                // Go down one level on the tree.
                row = row % size;                           // Update local row
                col = col % size;                           // Update local column

                if (l == levels) {
                    break;
                }

                // ************************************//
                // Update attributes (log raster)      //
                // ************************************//
                if (has_children) {
                    pos_log = this->get_children_position(pos_log, l);
                }

                // ************************************//
                // Update attributes (snapshot)        //
                // ************************************//
                if (!is_uniform_snap) {
                    pos_snap = raster_snap.get_children_position(pos_snap, l);
                }

                k = this->get_k(l);                         // Update k
            }
        } // END BLOCK searching from level 1 to level l - 1


        // ************************************//
        // Searching last level                //
        // ************************************//
        return get_cell_last_level(raster_snap, row, col, pos_log, has_children,
                                   pos_snap, is_uniform_snap, max_snap, diff_max, size);

    }

    //*******************************************************//
    //********************** FILE ***************************//
    //*******************************************************//
    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {

        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        // Plain parameters
        written_bytes += serialize_base(out, child, name);
        sdsl::serialize(k_plain_values, out, child, "plain_values");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }


    void load(std::istream& in) {
        // Plain parameters
        load_base(in);
        sdsl::load(k_plain_values, in);;
    }
protected:

    //*******************************************************//
    //********************** BUILD **************************//
    //*******************************************************//
    template<class Container, class Snapshot_type>
    void build(Container &&c_values, Snapshot_type &&raster_snapshot, std::vector<value_type> &plain_values_) {


        std::vector<std::vector<value_type>> min_values_(this->k_height);
        std::vector<std::vector<value_type>> max_values_(this->k_height);
        std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);
        std::vector<sdsl::int_vector<1>> tmp_eqB_(this->k_height-1);

        // Conceptual nodes
        k2raster_node<value_type> root_snap = raster_snapshot.get_root();
        node<value_type> root_node = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};
        node<value_type> root_diff = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};

        build(c_values, raster_snapshot,root_snap, root_node, root_diff,
                min_values_, max_values_, plain_values_, tmp_t_, tmp_eqB_,
                0, 0, 0, this->k_size);

        // Encode root node
        this->k_min_value = root_node.min_value;
        this->k_max_value = root_node.max_value;

        size_type n_nodes = 0; // Count the total number of nodes

        // Encode max values;
        n_nodes = this->build_max_values(max_values_, tmp_t_);
        max_values_.clear();
        max_values_.shrink_to_fit();

        // Encode min values
        this->build_min_values(min_values_);
        min_values_.clear();
        min_values_.shrink_to_fit();

        // Encode bitmap T (copy temporal bitmap T)
        this->build_t(tmp_t_, n_nodes);
        tmp_t_.clear();
        tmp_t_.shrink_to_fit();

        // Encode bitmap eqB (copy temporal bitmap eqB)
        this->build_eqB(tmp_eqB_, n_nodes);
        tmp_eqB_.clear();
        tmp_eqB_.shrink_to_fit();
    }

        template<class Container, class Snapshot_type>
        void build(Container &&c_values, Snapshot_type &&raster_snapshot,
                k2raster_node<value_type> &snap_node, node<value_type> &child_node, node<value_type> &diff_node,
                std::vector<std::vector<value_type>> &min_values_, std::vector<std::vector<value_type>> &max_values_,
                std::vector<value_type> &plain_values_, std::vector<sdsl::int_vector<1>> &tmp_t_, std::vector<sdsl::int_vector<1>> &tmp_eqB_,
                size_type base_row, size_type base_col, ushort level, size_type sub_size){

            size_type pos_value, child_base_row, child_base_col;

            // Initialize nodes
            child_node = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};
            diff_node = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};

            if (level == this->k_level_k1 + this->k_level_k2 + 1) {
                value_type value;
                /****************/
                /* PLAIN VALUES */
                /****************/
                // Seek the max and min value
                std::vector<value_type> values_snap = raster_snapshot.get_plain_values(snap_node);

                for (size_type r = base_row; r < base_row + sub_size; r++) {
                    for (size_type c = base_col; c < base_col + sub_size; c++) {
                        if (r < this->k_real_size_x && c < this->k_real_size_y) {
                            pos_value = r * this->k_real_size_y + c;
                            value = c_values[pos_value];
                            if (value < child_node.min_value) {
                                child_node.min_value = value;
                            }
                            if (value > child_node.max_value) {
                                child_node.max_value = value;
                            }

                            // Diff values
                            int diff = value -  values_snap[(r-base_row) * sub_size + (c-base_col)];
                            if (diff < diff_node.min_value) {
                                diff_node.min_value = diff;
                            }
                            if (diff > diff_node.max_value) {
                                diff_node.max_value = diff;
                            }
                        }
                    } // END for c
                } // END for r

                if (child_node.min_value != child_node.max_value && diff_node.min_value != diff_node.max_value){
                    // Add values to plain_values array
                    for (size_type r = base_row; r < base_row + sub_size; r++) {
                        for (size_type c = base_col; c < base_col + sub_size; c++) {
                            if (r < this->k_real_size_x && c < this->k_real_size_y) {
                                pos_value = r * this->k_real_size_y + c;
                                plain_values_.push_back(child_node.max_value - c_values[pos_value]);
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
                    child_node.min_value = c_values[pos_value];
                    child_node.max_value = child_node.min_value;

                    diff_node.min_value = child_node.min_value - snap_node.min_value;
                    diff_node.max_value = diff_node.min_value;
                } else {
                    /*****************/
                    /* INTERNAL NODE */
                    /*****************/
                    std::vector<k2raster_node<value_type>> children_snap(k * k);
                    std::vector<node<value_type>> children_node(k * k);
                    std::vector<node<value_type>> children_diff(k * k);
                    size_type pos = 0;

                    for (uint x = 0; x < k; x++) {
                        for (uint y = 0; y < k; y++) {
                            child_base_row = base_row + x * sub_size;
                            child_base_col = base_col + y * sub_size;
                            if (child_base_row < this->k_real_size_x && child_base_col < this->k_real_size_y) {
                                children_snap[pos] = raster_snapshot.get_child(snap_node, x*k + y);

                                // Recursive search
                                build(c_values, raster_snapshot,
                                        children_snap[pos], children_node[pos], children_diff[pos],
                                        min_values_, max_values_, plain_values_, tmp_t_, tmp_eqB_,
                                        child_base_row, child_base_col, level + 1, sub_size);

                                // Log raster
                                if (children_node[pos].min_value < child_node.min_value) {
                                    child_node.min_value = children_node[pos].min_value;
                                }
                                if (children_node[pos].max_value > child_node.max_value) {
                                    child_node.max_value = children_node[pos].max_value;
                                }

                                // diff raster
                                if (children_diff[pos].min_value < diff_node.min_value) {
                                    diff_node.min_value = children_diff[pos].min_value;
                                }
                                if (children_diff[pos].max_value > diff_node.max_value) {
                                    diff_node.max_value = children_diff[pos].max_value;
                                }
                            } else {
                                // Positions out of real matrix
                                children_node[pos] = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};
                                children_diff[pos] = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};
                            }
                            pos++;
                        } // END for y
                    } // END for x

                    // Check if matrix is not empty or uniform
                    if (child_node.min_value < child_node.max_value && diff_node.min_value != diff_node.max_value) {
                        // --------------------------------------------------------------------- //
                        // - If max t and min t are identical (or if we reach a 1×1 submatrix). Being 'zt' the position in the final bitmap T:
                        //      * We set T[zt] ← 0
                        //      * eqB[rank0[T, zt]] ← 0
                        //      * Lmax[zt] ← (maxt − maxs)
                        // -  If all the values in Ms and Md differ only in a unique value α (or if they are identical, hence α = 0),
                        //      * we set T[zt ] ← 0
                        //      * eqB[rank0[T, z t ]] ← 1
                        //      * Lmax[zt] ← (maxt − maxs).
                        // - Otherwise:
                        //      * We set T[zt] ← 1
                        //      * Lmax[zt] ← (maxt − maxs)
                        //      * Lmin[rank1(zt)] ← (mint − mins).
                        // --------------------------------------------------------------------- //
                        if (level != (this->k_height-1)) {                                                 // Except last level
                            tmp_t_[level].resize(tmp_t_[level].size() + (k * k)); // TODO improve this
                        }
                        for (size_type c = 0; c < (k*k); c++) {
                            if (children_node[c].min_value > children_node[c].max_value) {
                                max_values_[level].emplace_back(0); // Submatrix outside of the original matrix
                            } else {
                                max_values_[level].emplace_back(this->__encode(children_node[c].max_value - children_snap[c].max_value));
                            }

                            if ((level == (this->k_height-1))) {
                                continue; // Do not fill t and eqB at the last level
                            }

                            if (children_node[c].min_value > children_node[c].max_value) {
                                // Submatrix outside of the original matrix, push a 0
                                tmp_t_[level][max_values_[level].size() - 1] = 0;
                                tmp_eqB_[level].resize(tmp_eqB_[level].size() + 1); // TODO improve this
                                tmp_eqB_[level][tmp_eqB_[level].size() - 1] = 0;
                            } else {
                                if (children_node[c].min_value != children_node[c].max_value
                                    && children_diff[c].min_value != children_diff[c].max_value) {
                                    min_values_[level].push_back(this->__encode(children_node[c].min_value - children_snap[c].min_value));
                                    tmp_t_[level][max_values_[level].size() - 1] = 1;
                                } else {
                                    tmp_t_[level][max_values_[level].size() - 1] = 0;
                                    tmp_eqB_[level].resize(tmp_eqB_[level].size() + 1); // TODO improve this
                                    tmp_eqB_[level][tmp_eqB_[level].size() - 1] =
                                            (children_node[c].min_value != children_node[c].max_value);
                                } // END IF set t and eqB
                            } // END IF min_values_children[c] > max_values_children[c]
                        } // END FOR children
                    } // END IF (min_value < max_value)
                } // END IF NORMAL PARTITION
            } // END IF PLAIN LEVELS
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
    template<class k2_raster_snap>
    size_type get_cell_last_level(const k2_raster_snap &raster_snap,
                                         size_type row, size_type col, size_type pos_log, bool has_children,
                                         size_type pos_snap, bool is_uniform_snap,
                                         value_type max_snap, value_type diff_log, size_type size) const {
        if (!is_uniform_snap) {
            if (this->k_level_plain == 0) {
                // TODO Check this special case
                ushort k = raster_snap.get_k(this->k_height - 1);
                pos_snap += row * k + col;
                max_snap -=  raster_snap.get_max_value(this->k_height, pos_snap);
                        //raster_snap.list_max[this->k_height - 1][child_pos_snap - raster_snap.m_accum_max_values[this->k_height-2]];
            } else{
                if (!has_children) {
                    // Value is in plain form
                    return raster_snap.get_plain_value(row, col, max_snap, pos_snap, size) + diff_log;
                }
            }
        }

        if (has_children) {
            if (this->k_level_plain == 0) {
                // TODO Check this special case
                ushort k = this->get_k(this->k_height - 1);
                pos_log += row * k + col;
                diff_log = this->__decode(this->get_max_value(this->k_height, pos_log));
//                diff_log = this->__decode(this->k_list_max[this->k_height - 1][pos_log - this->k_accum_max_values[this->k_height-2]]);
            } else {
                // Value is in plain form
                ushort levels = get_cell_n_levels();
                diff_log = this->__decode(this->get_max_value(levels, pos_log));
//                diff_log = this->__decode(this->k_list_max[levels-1][pos_log - this->k_accum_max_values[levels-1]]);
                return get_plain_value(row, col, max_snap + diff_log, pos_log, size);
            }
        }
        return max_snap + diff_log;
    }


    virtual value_type get_cell_plain(size_type row, size_type col, size_type size,
                                      size_type children_pos, value_type father_value) const {
        row %= size;
        col %= size;
        size_type pos = children_pos + (row * size + col);
        return father_value - k_plain_values[pos];
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
        value_type value;
        size_type count_cells = 0, cell_pos;
        for (auto x = xini; x <= xend; x++) {
            size_type pos = children_pos + ((x % size) * size + (yini % size));
            cell_pos = (x - or_x) * window_size + (yini - or_y);
            for (auto y = yini; y <= yend; y++) {
                value = father_value - k_plain_values[pos];
                result[cell_pos++] = value;
                count_cells++;
                pos++; // Position of the next value
            } // END FOR y
        } // END FOR x
        return count_cells;
    }

    //*****************************//
    //***** CHECK VALUES WINDOW ***//
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
    //******************** NAVIGATION ***********************//
    //*******************************************************//

    inline size_type get_plain_value(size_type row, size_type col, value_type max_value,
                                     size_type child_pos, size_type size) const {
        size_type ones = this->k_t_rank1(child_pos) + 1;      // Number of non-empty nodes until position 'child_pos'
        ones = ones - (this->k_accum_min_values[this->k_list_min.size()-1] + 1);
        child_pos = ones * size * size;
        return get_cell_plain(row, col, size, child_pos, max_value);
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

    //*******************************************************//
    //******************** AUXILIARY ************************//
    //*******************************************************//
    virtual inline bool is_plain_level(ushort level) const {
        return this->k_level_plain > 0 && (level == (this->k_level_k1 + this->k_level_k2));
    }

    virtual short get_cell_n_levels() const {
        return std::min(this->k_level_k1 + this->k_level_k2+1, this->k_height-1);
    }

}; // END CLASS t_k2_raster_heuristic_log
} // END NAMESPACE k2raster

#endif // INCLUDED_T_K2_RASTER_PLAIN_LOG
