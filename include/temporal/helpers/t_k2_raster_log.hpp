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

#ifndef INCLUDED_T_K2_RASTER_LOG
#define INCLUDED_T_K2_RASTER_LOG

#include <k2_raster_base.hpp>


//! Namespace for k2-raster library
namespace k2raster {

    template<typename t_value=int,
            typename t_bv=sdsl::bit_vector,
            typename t_rank=sdsl::rank_support_v5<1,1>,
            typename t_values_vec=sdsl::dac_vector_dp_opt<t_bv, t_rank, 3>>
    class t_k2_raster_log : public k2_raster_base<t_value, t_bv, t_rank, t_values_vec> {

    public:
        typedef k2_raster_base<t_value, t_bv, t_rank, t_values_vec>  k2_raster_p; // k2_raster parent
        typedef t_value                                     value_type;
        typedef k2_raster_base<>::size_type                 size_type;
        typedef t_bv                                        bit_vector_type;
        typedef ulong                                       encode_type;
        typedef long                                        decode_type;

    protected:

        /****** Structures DIFF ******/
        bit_vector_type     k_eqB;            // Bit array to eqB

    public:

    //*******************************************************//
    //******************* CONSTRUCTOR ***********************//
    //*******************************************************//
    t_k2_raster_log() = default;

    t_k2_raster_log(const t_k2_raster_log &tr) : k2_raster_p() {
        *this = tr;
    }

    t_k2_raster_log(t_k2_raster_log &&tr) {
        *this = std::move(tr);
    }

    // Only interface (plain_levels unused)
    template<class Container, class Snapshot_type>
    t_k2_raster_log(Container &&c_values, Snapshot_type &&raster_snapshot, ushort plain_levels __attribute__((unused)))
    : t_k2_raster_log(c_values, raster_snapshot)  { };


    template<class Snapshot_type>
    t_k2_raster_log(Snapshot_type &&raster_snapshot, ushort type)
            : k2_raster_p(raster_snapshot.get_n_rows(), raster_snapshot.get_n_cols(),
                          raster_snapshot.k1, raster_snapshot.k2, raster_snapshot.level_k1, type) {
    }

    template<class Container, class Snapshot_type>
    t_k2_raster_log(Container &&c_values, Snapshot_type &&raster_snapshot)
    : t_k2_raster_log(raster_snapshot, K2_RASTER_TYPE) {

        this->init_levels();

        // Build k2-raster
        build(c_values, raster_snapshot);
    }

    //*******************************************************//
    //*************** BASIC OPERATIONS **********************//
    //*******************************************************//

    //! Move assignment operator
    t_k2_raster_log& operator=(t_k2_raster_log&& tr)
    {
        if (this != &tr) {
            k2_raster_p::operator=(tr);
            /****** Structures DIFF ******/
            k_eqB = std::move(tr.k_eqB);
        }
        return *this;
    }

    //! Assignment operator
    t_k2_raster_log& operator=(const t_k2_raster_log& tr)
    {
        if (this != &tr) {
            k2_raster_p::operator=(tr);
            /****** Structures DIFF ******/
            k_eqB = tr.k_eqB;
        }
        return *this;
    }

    //! Swap operator
    void swap(t_k2_raster_log& tr)
    {
        if (this != &tr) {
            k2_raster_p::swap(tr);
            /****** Structures DIFF ******/
            std::swap(k_eqB, tr.k_eqB);
        }
    }

    //! Equal operator
    bool operator==(const t_k2_raster_log& tr) const
    {
        if (!k2_raster_p::operator==(tr)) {
            return false;
        }
        if (k_eqB.size() != tr.eqB.size()) {
            return false;
        }
        for (unsigned i = 0; i < k_eqB.size(); i++)
            if (k_eqB[i] != tr.k_eqB[i])
                return false;
        return true;
    }

    //*******************************************************//
    //******************** QUERIES **************************//
    //*******************************************************//
    value_type get_cell(size_type  row __attribute__((unused)), size_type col __attribute__((unused))) const { return 0;}
    size_type get_cells_by_value(size_type xini __attribute__((unused)), size_type xend __attribute__((unused)),
            size_type yini __attribute__((unused)), size_type yend __attribute__((unused)),
            value_type valini __attribute__((unused)), value_type valend __attribute__((unused)),
            std::vector<std::pair<size_type, size_type>> &result __attribute__((unused))) {return 0;}
    size_type get_values_window(size_type xini __attribute__((unused)), size_type xend __attribute__((unused)),
            size_type yini __attribute__((unused)), size_type yend __attribute__((unused)),
            std::vector<value_type> &result __attribute__((unused))) {return 0; }
    bool check_values_window(size_type xini __attribute__((unused)), size_type xend __attribute__((unused)),
            size_type yini __attribute__((unused)), size_type yend __attribute__((unused)),
            value_type valini __attribute__((unused)), value_type valend __attribute__((unused)),
            bool strong_check __attribute__((unused))) {return false;};

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
        decode_type diff = this->k_max_value - raster_snap.max_value;       // Difference between the current raster and the snapshot

        // ************************************//
        // Snapshot attributes                 //
        // ************************************//
        size_type pos_snap = 0;                                             // Children position on bitmap k_t of the current node (snapshot)
        value_type max_value_snap = raster_snap.max_value;                  // Current max value (snapshot)
        bool is_uniform_snap = raster_snap.max_value == raster_snap.min_value; // The current node of the snapshot is uniform

        {
            ushort k = this->get_k(0);
            size_type child;

            for (auto l = 1; l < this->k_height; l++) {
                size /= k;                                  // Submatrix size
                child = (row / size) * k + (col / size);    // Number of child (following a z-order)
                pos_log += has_children ? child : 0;        // Position in bitmap k_t
                pos_snap += !is_uniform_snap ? child : 0;   // Position in bitmap k_t (snapshot)

                // Check snapshot node
                if (!is_uniform_snap) {
                    is_uniform_snap = !raster_snap.t[pos_snap];
                    max_value_snap -= raster_snap.get_max_value(l, pos_snap);
//                            raster_snap.list_max[l - 1][pos_snap - raster_snap.m_accum_max_values[l-1]];
                }

                // Check current raster node
                if (has_children) {
                    has_children = this->k_t[pos_log];
                    if (!has_children) {
                        is_uniform = !k_eqB[pos_log - this->k_t_rank1(pos_log)];
                        diff = this->__decode(this->get_max_value(l, pos_log));
//                        diff = this->__decode(this->k_list_max[l-1][pos_log - this->k_accum_max_values[l-1]]);
                    }
                }

                // Return value
                if ((is_uniform_snap && !has_children) || is_uniform) {
                    return max_value_snap + diff;
                }

                // Go down one level on the tree.
                row = row % size;                           // Update local row
                col = col % size;                           // Update local column

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
        {
            if (!is_uniform_snap) {
                ushort k = raster_snap.get_k(this->k_height-1);
                pos_snap += row * k + col;
                max_value_snap -= raster_snap.get_max_value(this->k_height, pos_snap);
//                        raster_snap.list_max[this->k_height-1][pos_snap - raster_snap.m_accum_max_values[this->k_height-1]];
            }

            if (has_children) {
                ushort k = this->get_k(this->k_height-1);
                pos_log += row * k + col;
                diff = __decode(this->get_max_value(this->k_height, pos_log));
//                diff = __decode(this->k_list_max[this->k_height-1][pos_log - this->k_accum_max_values[this->k_height-1]]);
            }

            return max_value_snap + diff;
        } // END BLOCK last levels
    }

    //*****************************//
    //***** GET CELL BY VALUE *****//
    //*****************************//
    template<class k2_raster_snap>
    size_type get_cells_by_value(k2_raster_snap &raster_snap, size_type xini, size_type xend, size_type yini, size_type yend,
                                         value_type valini, value_type valend, std::vector<std::pair<size_type, size_type>> &result) {

        // Check values of root node
        if (this->k_min_value > valend ||
            this->k_max_value < valini) {
            // There is no valid value in the matrix
            return 0;
        }

        bool has_children = true;
        // Whole matrix if uniform or/and all values are valid
        if (this->k_t.empty()) {
            if (this->k_min_value >= valini && this->k_max_value <= valend) {
                for (auto x = xini; x <= xend; x++) {
                    for (auto y = yini; y <= yend; y++) {
                        result.emplace_back(x, y);
                    }
                }
                return (xend - xini + 1) * (yend - yini + 1);
            }
            has_children = false;
        }

        ushort k = this->get_k(0);
        return get_cells_by_value_helper(raster_snap, {xini, xend, yini, yend, valini, valend},
                                         {false, 0}, {raster_snap.min_value == raster_snap.max_value, 0}, has_children,
                                         0, 0, raster_snap.min_value, raster_snap.max_value,
                                         this->k_size / k, 1, result);
    }

    //*****************************//
    //***** GET VALUES WINDOW *****//
    //*****************************//
    template<class k2_raster_snap>
    size_type get_values_window(k2_raster_snap &raster_snap, size_type xini, size_type xend, size_type yini, size_type yend,
                                std::vector<value_type> &result) {

        size_type count_cells = 0; // Number of cells
        result.resize((xend - xini + 1) * (yend - yini + 1));

        // Check first level
        value_type max_value = raster_snap.max_value;
        if (this->k_t.empty()) {
            if (this->k_min_value == this->k_max_value) {
                // All cell of the matrix are equal
                for (auto x = xini; x <= xend; x++){
                    for (auto y= yend; y <= yend; y++){
                        result[count_cells++] = this->k_max_value;;
                    }
                }
                return count_cells;
            }
            max_value = this->k_max_value;
        }

        ushort k = this->get_k(0);
        return get_values_window_helper(raster_snap, {xini, xend, yini, yend},
                                        {false, 0}, {raster_snap.max_value == raster_snap.min_value, 0}, !this->k_t.empty(), 0, 0, max_value,
                                        this->k_size / k, 1, result, xini, yini, (yend - yini + 1));
    }

    //*****************************//
    //**** CHECK VALUES WINDOW ****//
    //*****************************//
    template<class k2_raster_snap>
    bool check_values_window(k2_raster_snap &raster_snap, size_type xini, size_type xend, size_type yini, size_type yend,
                                     value_type valini, value_type valend, bool strong_check) {

        // Check values of root node
        if (this->k_min_value > valend ||
            this->k_max_value < valini) {
            // There is no valid value in the matrix
            return false;
        }

        bool has_children = true;

        // Check first level
        if (this->k_t.empty()) {
            if (this->k_min_value == this->k_max_value) {
                // All cells are valid
                return true;
            }
            has_children = false;
        }

        // Size of window is equal to size of matrix
        if (xini == 0 && yini == 0 && xend == (this->k_real_size_x - 1) && yend == (this->k_real_size_y - 1)){
            return (this->k_min_value >= valini && this->k_max_value <= valend) || !strong_check;
        }

        ushort k = this->get_k(0);
        return check_values_window_helper(raster_snap, {xini, xend, yini, yend, valini, valend},
                                          {false, 0}, {raster_snap.min_value == raster_snap.max_value, 0}, has_children,
                                          0, 0, raster_snap.min_value, raster_snap.max_value,
                                          this->k_size / k, 1, strong_check);
    }

    //*******************************************************//
    //********************** FILE ***************************//
    //*******************************************************//
    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {

        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += k2_raster_p::serialize(out, child, name);
        written_bytes += k_eqB.serialize(out, child, "eqB");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }


    void load(std::istream& in) {
        k2_raster_p::load(in);
        k_eqB.load(in);
    }

protected:

    //*******************************************************//
    //********************** BUILD **************************//
    //*******************************************************//
    template<class Container, class Snapshot_type>
    void build(Container &&c_values, Snapshot_type &&raster_snapshot) {
        std::vector<std::vector<value_type>> min_values_(this->k_height);
        std::vector<std::vector<value_type>> max_values_(this->k_height);
        std::vector<sdsl::int_vector<1>> tmp_t_(this->k_height-1);
        std::vector<sdsl::int_vector<1>> tmp_eqB_(this->k_height-1);

        // Conceptual nodes
        k2raster_node<value_type> root_snap = raster_snapshot.get_root();
        node<value_type> root_node = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};
        node<value_type> root_diff = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};

        build(c_values, raster_snapshot, root_snap, root_node, root_diff,
              min_values_, max_values_, tmp_t_, tmp_eqB_,
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
        build_eqB(tmp_eqB_, n_nodes);
        tmp_eqB_.clear();
        tmp_eqB_.shrink_to_fit();
    }

    void build_eqB(std::vector<sdsl::int_vector<1>> &tmp_eqB_, size_type n_nodes) {
        sdsl::int_vector<1> eqB_(n_nodes);
        size_type n = 0;
        for(auto l = 0; l < tmp_eqB_.size(); l++) {
            for (auto c = 0; c < tmp_eqB_[l].size(); c++) {
                eqB_[n++] = tmp_eqB_[l][c];
            }
        }
        eqB_.resize(n);
        this->k_eqB = bit_vector_type(eqB_);
    }

    template<class Container, class Snapshot_type>
    void build(Container &&c_values, Snapshot_type &&raster_snapshot,
            k2raster_node<value_type> &snap_node, node<value_type> &child_node, node<value_type> &diff_node,
            std::vector<std::vector<value_type>> &min_values_, std::vector<std::vector<value_type>> &max_values_,
            std::vector<sdsl::int_vector<1>> &tmp_t_, std::vector<sdsl::int_vector<1>> &tmp_eqB_,
            size_type base_row, size_type base_col, ushort level, size_type sub_size) {

        size_type pos_value, child_base_row, child_base_col;

        // Initialize nodes
        child_node = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};
        diff_node = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min()};

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
                        build(c_values, raster_snapshot, children_snap[pos], children_node[pos], children_diff[pos],
                                min_values_, max_values_, tmp_t_, tmp_eqB_,
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
                        children_snap[pos] = {std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::min(), 0, 0};
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
                        max_values_[level].emplace_back(__encode(children_node[c].max_value - children_snap[c].max_value));
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
                            min_values_[level].push_back(__encode(children_node[c].min_value - children_snap[c].min_value));
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
        } // END IF (level == (this->k_height)
    }


    //*******************************************************//
    //**************** QUERIES HELPERS **********************//
    //*******************************************************//

    //*****************************//
    //***** GET CELL BY VALUE *****//
    //*****************************//
    template<class k2_raster_snap>
    size_type get_cells_by_value_helper(
            k2_raster_snap &raster_snap, query_info_v<size_type, value_type> &&query,
            node_info<size_type> &&node_log, node_info<size_type> &&node_snap, bool has_children,
            size_type base_x, size_type base_y, value_type f_min, value_type f_max,
            size_type children_size, ushort l,
            std::vector<std::pair<size_type, size_type>> &result) {

        size_type c_base_x, c_base_y;
        size_type pos_log=0, pos_snap=0, ones_snap=0, ones_log=0;
        size_type count_cells = 0;
        value_type max_snap, min_snap;
        decode_type diff_max, diff_min;
        bool c_is_uniform, c_has_children, c_is_uniform_snap;

        ushort k = this->get_k(l-1);
        bool is_leaf = l == this->k_height; // True -> They are leaves

        // Set limits (children with cell that overlap with the interesting region)
        size_type limit_i_x = (query.xini - base_x) / children_size;
        size_type limit_e_x = (query.xend - base_x) / children_size;
        size_type limit_i_y = (query.yini - base_y) / children_size;
        size_type limit_e_y = (query.yend - base_y) / children_size;

        // Check children
        for (auto x = limit_i_x; x <= limit_e_x; x++) {
            // Calculate base position of the current child (row)
            c_base_x = base_x + x * children_size;

            for (auto y = limit_i_y; y <= limit_e_y; y++) {

                // Calculate base position of the current child (column)
                c_base_y = base_y + y * children_size;

                // Check snapshot and get max value (snapshot)
                if (!node_snap.is_uniform) {
                    pos_snap = node_snap.children_pos + x * k + y;
                    c_is_uniform_snap = is_leaf || !raster_snap.t[pos_snap];
                    max_snap = f_max - raster_snap.get_max_value(l, pos_snap);
                            //raster_snap.list_max[l-1][pos_snap - raster_snap.m_accum_max_values[l-1]];
                } else {
                    // Parent is uniform
                    c_is_uniform_snap = true;
                    max_snap = f_max; // Get parent's value
                }

                // Check log and get max value (log)
                if (has_children) {
                    pos_log = node_log.children_pos + x * k + y;
                    diff_max = this->__decode(this->get_max_value(l, pos_log));
//                    diff_max = this->__decode(this->k_list_max[l-1][pos_log - this->k_accum_max_values[l-1]]);

                    if (!is_leaf) {
                        c_has_children = this->k_t[pos_log];
                        c_is_uniform = !c_has_children ? !this->k_eqB[pos_log - this->k_t_rank1(pos_log)] : false;
                    } else {
                        c_has_children = false;
                        c_is_uniform = true;
                    }
                } else {
                    c_has_children = false;
                    diff_max = 0;
                    c_is_uniform = node_log.is_uniform;
                }

                if (query.valini > (max_snap + diff_max)) {
                    // Values out of range
                    continue;
                }

                // Get min value (snapshot)
                if (!c_is_uniform_snap && !c_is_uniform) {
//                    ones_snap = raster_snap.t_rank1(pos_snap)+1;
//                    min_snap = f_min + raster_snap.get_min_value(l, ones_snap-1);
                            //raster_snap.list_min[l-1][ones_snap - raster_snap.m_accum_min_values[l-1] - 1];

                    min_snap = f_min + raster_snap.get_min_value(l, pos_snap, ones_snap);
                } else {
                    min_snap = max_snap;
                }

                if (has_children) {
                    // Get min value (log)
                    if (c_has_children) {
                        // If it is not uniform or a leaf, its minimum value is obtained
                        ones_log = this->k_t_rank1(pos_log) + 1;
                        diff_min = this->__decode(this->get_min_value(l, ones_log-1));
//                                this->k_list_min[l - 1][ones_log - this->k_accum_min_values[l - 1] - 1]);
                    } else {
                        // It is uniform (or a leaf)
                        diff_min = c_is_uniform ? (max_snap + diff_max) - min_snap : diff_max;
                    }
                } else {
                    diff_min = 0;
                }

                if (query.valend < (min_snap + diff_min)) {
                    // Values out of range
                    continue;
                }

                // Current child may contain valid cells
                // Calculate positions with overlap with interesting region
                auto c_xini = std::max(c_base_x, query.xini);
                auto c_xend = std::min(c_base_x + children_size - 1, query.xend);
                auto c_yini = std::max(c_base_y, query.yini);
                auto c_yend = std::min(c_base_y + children_size - 1, query.yend);

                if (query.valini <= (min_snap + diff_min) && query.valend >= (max_snap + diff_max)) {
                    // Add all valid positions to final result
                    for (auto x = c_xini; x <= c_xend; x++) {
                        for (auto y = c_yini; y <= c_yend; y++) {
                            result.emplace_back(x, y);
                        }
                    }
                    count_cells += (c_xend - c_xini + 1) * (c_yend - c_yini + 1);
                } else {
                    // Search their children
                    if (c_has_children || !c_is_uniform_snap) {
                        // Go down one level on the tree

                        if (has_children && (!c_has_children || this->is_plain_level(l-1))) {
                            // We reach a leaf in 'log'. Apply the differences permanently
                            max_snap += diff_max;
                            min_snap += diff_min;
                        }

                        if (this->is_plain_level(l-1)) {
                            // Get value from plain values
                            if (c_has_children) {
                                // Check plain values at log
                                ones_log -= (this->k_accum_min_values[l-1]+1);
                                count_cells += this->get_cells_by_value_plain(c_xini, c_xend, c_yini, c_yend, query.valini, query.valend,
                                                                        children_size, ones_log * children_size * children_size, max_snap, result);
                            } else {
                                // Check plain values at snap
                                ones_snap -= raster_snap.m_accum_min_values[l-1];
                                count_cells += raster_snap.get_cells_by_value_plain(c_xini, c_xend, c_yini, c_yend, query.valini, query.valend,
                                                                                    children_size, ones_snap * children_size * children_size, max_snap, result);
                            } // END IF plain_values
                        } else {
                            // Go down one level
                            size_type new_children_pos = c_has_children ? this->get_children_position_ones(ones_log, l) : 0;
                            size_type new_children_pos_snap = !c_is_uniform_snap ? raster_snap.get_children_position_ones( ones_snap+1, l) : 0;
                            count_cells += get_cells_by_value_helper(raster_snap, {c_xini, c_xend, c_yini, c_yend, query.valini, query.valend},
                                                                     {c_is_uniform, new_children_pos}, {c_is_uniform_snap, new_children_pos_snap},
                                                                     c_has_children, c_base_x, c_base_y, min_snap, max_snap, children_size / this->get_k(l), l+1, result);
                        } // END IF plain_values
                    } // END IF is not uniform
                } // END IF check children
            } // END FOR y
        } // END FOR x
        return count_cells;
    }


    //*****************************//
    //***** GET VALUES WINDOW *****//
    //*****************************//
    template<class k2_raster_snap>
    size_type get_values_window_helper(
            k2_raster_snap &raster_snap, const query_info<size_type> &&query,
            node_info<size_type> &&node_log, node_info<size_type> &&node_snap, bool has_children,
            size_type base_x, size_type base_y, value_type f_max,
            size_type children_size, ushort l, std::vector<value_type> &result,
            size_type or_x, size_type or_y, size_type window_size) {


        size_type c_base_x, c_base_y;
        size_type pos_log=0, pos_snap=0;
        size_type count_cells = 0;
        value_type max_snap;
        decode_type diff_max;
        bool c_is_uniform=true, c_has_children, c_is_uniform_snap=true;

        ushort k = this->get_k(l-1);
        bool is_leaf = l == this->k_height; // True -> They are leaves

        // Set limits (children with cell that overlap with the interesting region)
        size_type limit_i_x = (query.xini - base_x) / children_size;
        size_type limit_e_x = (query.xend - base_x) / children_size;
        size_type limit_i_y = (query.yini - base_y) / children_size;
        size_type limit_e_y = (query.yend - base_y) / children_size;

        // Check children
        for (auto x = limit_i_x; x <= limit_e_x; x++) {
            // Calculate base position of the current child (row)
            c_base_x = base_x + x * children_size;

            for (auto y = limit_i_y; y <= limit_e_y; y++) {

                // Calculate base position of the current child (column)
                c_base_y = base_y + y * children_size;

                // Check snapshot and get max value (snapshot)
                if (!node_snap.is_uniform) {
                    pos_snap = node_snap.children_pos + x * k + y;
                    c_is_uniform_snap = is_leaf || !raster_snap.t[pos_snap];
                    max_snap = f_max - raster_snap.get_max_value(l, pos_snap);
                            //raster_snap.list_max[l-1][pos_snap - raster_snap.m_accum_max_values[l-1]];
                } else {
                    // Parent is uniform
                    c_is_uniform_snap = true;
                    max_snap = f_max; // Get parent's value
                }

                // Check log and get max value (log)
                if (has_children) {
                    pos_log = node_log.children_pos + x * k + y;
                    diff_max = this->__decode(this->get_max_value(l, pos_log));
//                    diff_max = this->__decode(this->k_list_max[l-1][pos_log - this->k_accum_max_values[l-1]]);

                    if (!is_leaf) {
                        c_has_children = this->k_t[pos_log];
                        c_is_uniform = !c_has_children ? !k_eqB[pos_log - this->k_t_rank1(pos_log)] : false;
                    } else {
                        c_has_children = false;
                        c_is_uniform = true;
                    }
                } else {
                    c_has_children = false;
                    diff_max = 0;
                    c_is_uniform = node_log.is_uniform;
                }

                // Current child may contain valid cells
                // Calculate positions with overlap with interesting region
                auto c_xini = std::max(c_base_x, query.xini);
                auto c_xend = std::min(c_base_x + children_size - 1, query.xend);
                auto c_yini = std::max(c_base_y, query.yini);
                auto c_yend = std::min(c_base_y + children_size - 1, query.yend);

                if (c_is_uniform || (c_is_uniform_snap && !c_has_children)) {
                    size_type cell_pos;
                    // Add all valid positions to final result
                    for (auto c_x = c_xini; c_x <= c_xend; c_x++) {
                        cell_pos = (c_x - or_x) * window_size + (c_yini - or_y);
                        for (auto c_y = c_yini; c_y <= c_yend; c_y++) {
                            result[cell_pos++] = max_snap + diff_max;
                            count_cells++;
                        }
                    }
                } else {
                    // Go down one level on the tree
                    if ((has_children && !c_has_children) || this->is_plain_level(l-1)) {
                        // We reach a leaf in 'log'. Apply the differences permanently
                        max_snap += diff_max;
                    }

                    if (this->is_plain_level(l-1)) {
                        // Get value from plain values
                        size_type ones;
                        if (c_has_children) {
                            // Check plain values at log
                            ones = this->k_t_rank1(pos_log) - this->k_accum_min_values[l-1];
                            count_cells += this->get_values_window_plain(c_xini, c_xend, c_yini, c_yend,
                                                                   children_size, ones * children_size * children_size, max_snap,
                                                                   result, or_x, or_y, window_size);
                        } else {
                            // Check plain values at snap
                            ones = raster_snap.t_rank1(pos_snap) - raster_snap.m_accum_min_values[l-1];
                            count_cells += raster_snap.get_values_window_plain(c_xini, c_xend, c_yini, c_yend,
                                                                               children_size, ones * children_size * children_size, max_snap,
                                                                               result, or_x, or_y, window_size);
                        } // END IF plain_values
                    } else {
                        // Go down one level
                        size_type new_children_pos = c_has_children ? this->get_children_position(pos_log, l) : 0;
                        size_type new_children_pos_snap = !c_is_uniform_snap ? raster_snap.get_children_position(pos_snap, l) : 0;
                        count_cells += get_values_window_helper(raster_snap, {c_xini, c_xend, c_yini, c_yend},
                                                                {c_is_uniform, new_children_pos},
                                                                {c_is_uniform_snap, new_children_pos_snap},
                                                                c_has_children, c_base_x, c_base_y, max_snap,
                                                                children_size / this->get_k(l), l + 1, result, or_x,
                                                                or_y, window_size);
                    } // END IF is_plain_level
                } // END IF go down one level
            } // END FOR y
        } // END FOR x
        return count_cells;
    }

    //*****************************//
    //**** CHECK VALUES WINDOW ****//
    //*****************************//
    template<class k2_raster_snap>
    bool check_values_window_helper(
            k2_raster_snap &raster_snap, const query_info_v<size_type, value_type> &&query,
            node_info<size_type> &&node_log, node_info<size_type> &&node_snap, bool has_children,
            size_type base_x, size_type base_y, value_type f_min, value_type f_max,
            size_type children_size, ushort l,
            bool strong_check) {

        size_type c_base_x, c_base_y;
        size_type pos_log=0, pos_snap=0, ones_snap=0, ones_log=0;
        value_type max_snap, min_snap;
        decode_type diff_max, diff_min;
        bool c_is_uniform, c_has_children, c_is_uniform_snap;

        ushort k = this->get_k(l-1);
        bool is_leaf = l == this->k_height; // True -> They are leaves

        // Set limits (children with cell that overlap with the interesting region)
        size_type limit_i_x = (query.xini - base_x) / children_size;
        size_type limit_e_x = (query.xend - base_x) / children_size;
        size_type limit_i_y = (query.yini - base_y) / children_size;
        size_type limit_e_y = (query.yend - base_y) / children_size;

        // Check children
        for (auto x = limit_i_x; x <= limit_e_x; x++) {
            // Calculate base position of the current child (row)
            c_base_x = base_x + x * children_size;

            for (auto y = limit_i_y; y <= limit_e_y; y++) {

                // Calculate base position of the current child (column)
                c_base_y = base_y + y * children_size;

                // Check snapshot and get max value (snapshot)
                if (!node_snap.is_uniform) {
                    pos_snap = node_snap.children_pos + x * k + y;
                    c_is_uniform_snap = is_leaf || !raster_snap.t[pos_snap];
                    max_snap = f_max - raster_snap.get_max_value(l, pos_snap);
                            //raster_snap.list_max[l-1][pos_snap -  raster_snap.m_accum_max_values[l-1]];
                } else {
                    // Parent is uniform
                    c_is_uniform_snap = true;
                    max_snap = f_max; // Get parent's value
                }

                // Check log and get max value (log)
                if (has_children) {
                    pos_log = node_log.children_pos + x * k + y;
                    diff_max = this->__decode(this->get_max_value(l, pos_log));
//                    diff_max = this->__decode(this->k_list_max[l-1][pos_log - this->k_accum_max_values[l-1]]);

                    if (!is_leaf) {
                        c_has_children = this->k_t[pos_log];
                        c_is_uniform = !c_has_children ? !this->k_eqB[pos_log - this->k_t_rank1(pos_log)] : false;
                    } else {
                        c_has_children = false;
                        c_is_uniform = true;
                    }
                } else {
                    c_has_children = false;
                    diff_max = 0;
                    c_is_uniform = node_log.is_uniform;
                }

                // Check max value
                if (query.valini > (max_snap + diff_max)) {
                    // Values out of range
                    if (strong_check) return false; // (strong) At least one cell is within the range of values
                    continue;
                }

                // Get min value (snapshot)
                if (!c_is_uniform_snap && !c_is_uniform) {
                    ones_snap = raster_snap.t_rank1(pos_snap)+1;
                    min_snap = f_min + raster_snap.get_min_value(l, ones_snap-1);
                            //raster_snap.list_min[l-1][ones_snap -  raster_snap.m_accum_min_values[l-1] - 1];
                } else {
                    min_snap = max_snap;
                }

                // Get min value (log)
                if (has_children) {
                    if (c_has_children) {
                        // If it is not uniform or a leaf, its minimum value is obtained
                        ones_log = this->k_t_rank1(pos_log)+1;
                        diff_min = this->__decode(this->get_min_value(l, ones_log-1));
//                        diff_min = this->__decode(this->k_list_min[l-1][ones_log - this->k_accum_min_values[l-1] - 1]);
                    } else {
                        // It is uniform (or a leaf)
                        diff_min = c_is_uniform ? (max_snap + diff_max) - min_snap : diff_max;
                    }
                } else {
                    diff_min = 0;
                }

                // Check min value
                if (query.valend < (min_snap + diff_min)) {
                    // Values out of range
                    if (strong_check) return false; // (strong) At least one cell is not within the range of values
                    continue;
                }

                if (query.valini <= (min_snap + diff_min)  && query.valend >= (max_snap + diff_max)) {
                    // All cells of the submatrix is within the range of values
                    if (!strong_check) return true; // (weak) At least one cell is within the range of values
                    continue;
                } else {
                    // Current child may contain valid cells
                    // Calculate positions with overlap with interesting region
                    auto c_xini = std::max(c_base_x, query.xini);
                    auto c_xend = std::min(c_base_x + children_size - 1, query.xend);
                    auto c_yini = std::max(c_base_y, query.yini);
                    auto c_yend = std::min(c_base_y + children_size - 1, query.yend);

                    if (query.valini <= (min_snap + diff_min) && query.valend >= (max_snap + diff_max)) {
                        // The current submatrix is completely inside the search window
                        // and contains some cell outside the range of values

                        // (strong) At least one cell is not within the range of values
                        // (weak) At least one cell is within the range of values
                        return !strong_check;
                    } else {
                        // Go down one level on the tree
                        if (has_children && (!c_has_children || this->is_plain_level(l-1))) {
                            // We reach a leaf in 'log'. Apply the differences permanently
                            max_snap += diff_max;
                            min_snap += diff_min;
                        }

                        bool result;
                        if (this->is_plain_level(l-1)) {
                            // Get value from plain values
//                            size_type ones;
                            if (c_has_children) {
                                // Check plain values at log
//                                ones = this->k_t_rank1(pos_log) + 1;
                                ones_log -= (this->k_accum_min_values[l-1]+1);
                                result = this->check_values_window_plain(c_xini, c_xend, c_yini, c_yend, query.valini, query.valend,
                                                                   children_size, ones_log * children_size * children_size, max_snap, strong_check);
                            } else {
                                // Check plain values at snap
//                                ones = raster_snap.t_rank1(pos_snap) + 1;
                                ones_snap -= (raster_snap.m_accum_min_values[l-1]+1);
                                result = raster_snap.check_values_window_plain(c_xini, c_xend, c_yini, c_yend, query.valini, query.valend,
                                                                               children_size, ones_snap * children_size * children_size, max_snap, strong_check);
                            } // END IF plain_values
                        } else {
                            // Go down one level
                            size_type new_children_pos = c_has_children ? this->get_children_position_ones(ones_log, l) : 0;
                            size_type new_children_pos_snap = !c_is_uniform_snap ? raster_snap.get_children_position_ones( ones_snap, l) : 0;
                            result = check_values_window_helper(raster_snap,
                                                                {c_xini, c_xend, c_yini, c_yend, query.valini,
                                                                 query.valend},
                                                                {c_is_uniform, new_children_pos},
                                                                {c_is_uniform_snap, new_children_pos_snap},
                                                                c_has_children, c_base_x, c_base_y, min_snap,
                                                                max_snap,
                                                                children_size / this->get_k(l), l + 1,
                                                                strong_check);
                        } // END IF plain_values
                        if (strong_check && !result) return false;
                        if (!strong_check && result) return true;

                    } // END IF check children
                } // END IF check values min-max
            } // END FOR y
        } // END FOR x

        // (strong) No values outside the range of values found
        // (weak) No values found within the range of values
        return strong_check;
    }

    //*******************************************************//
    //********************* HELPERS *************************//
    //*******************************************************//
    decode_type __decode(encode_type value) const {
        return (value %2) ? (value+1)/2 : -(value/2);
    }

    encode_type __encode(decode_type value) const {
        return (value <= 0) ? -value * 2 : 2 * value - 1;
    }
}; // END CLASS t_k2_raster_log
} // END NAMESPACE k2raster

#endif // INCLUDED_T_K2_RASTER_LOG
