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
            typename t_rank=typename t_bv::rank_1_type,
            typename t_values_vec=sdsl::dac_vector_dp<>,
            typename t_values_last=sdsl::dac_vector_dp_opt<t_bv, t_rank, 3>>
    class t_k2_raster_plain_log : public k2_raster_base<t_value, t_bv, t_rank, t_values_vec> {

    public:
        typedef k2_raster_base<t_value, t_bv, t_rank, t_values_vec>  k2_raster_p; // k2_raster parent
        typedef t_value                                     value_type;
        typedef k2_raster_base<>::size_type                 size_type;
        typedef t_bv                                        bit_vector_type;
        typedef ulong                                       encode_type;
        typedef long                                        decode_type;

    protected:
        ushort k_level_k2;
        ushort k_level_plain;

        size_type k_size_leaves;
        t_values_last k_plain_values;

    public:
        const size_type &m_size_leaves = k_size_leaves;
        const t_values_last &m_plain_values = k_plain_values;

    protected:

        /****** Structures DIFF ******/
        bit_vector_type     k_eqB;            // Bit array to eqB

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
            : k2_raster_p(nullptr, raster_snapshot.get_n_rows(), raster_snapshot.get_n_cols(),
                          raster_snapshot.k1, raster_snapshot.k2, raster_snapshot.level_k1) {
        this->k_k2_raster_type = k2_raster_type;
        raster_snapshot.copy_parameters(this->k_level_k1, k_level_k2, k_level_plain, k_size_leaves);
    }

    template<class Container, class Snapshot_type>
    t_k2_raster_plain_log(Container &&c_values, Snapshot_type &&raster_snapshot)
            : k2_raster_p(c_values, raster_snapshot.get_n_rows(), raster_snapshot.get_n_cols(),
                    raster_snapshot.k1, raster_snapshot.k2, raster_snapshot.level_k1) {

        this->k_k2_raster_type = K2_RASTER_TYPE;
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

            /****** Structures DIFF ******/
            k_eqB = std::move(tr.k_eqB);
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

            /****** Structures DIFF ******/
            k_eqB = tr.k_eqB;
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

            /****** Structures DIFF ******/
            std::swap(k_eqB, tr.k_eqB);
        }
    }

    //! Equal operator
    bool operator==(const t_k2_raster_plain_log& tr) const
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
        size_type child_pos = 0;                                            // Children position on bitmap k_t of the current node
        decode_type diff_log = this->k_max_value - raster_snap.max_value;       // Difference between the current raster and the snapshot

        // ************************************//
        // Snapshot attributes                 //
        // ************************************//
        size_type child_pos_snap = 0;                                       // Children position on bitmap k_t of the current node (snapshot)
        value_type max_snap = raster_snap.max_value;                  // Current max value (snapshot)
        bool is_uniform_snap = raster_snap.max_value == raster_snap.min_value; // The current node of the snapshot is uniform

        ushort levels = std::min(this->k_level_k1 + this->k_level_k2+1, this->k_height-1);
        {
            ushort k = this->get_k(0);
            size_type child;

            for (auto l = 1; l <= levels; l++) {
                size /= k;                                      // Submatrix size
                child = (row / size) * k + (col / size);        // Number of child (following a z-order)
                child_pos += has_children ? child : 0;          // Position in bitmap k_t
                child_pos_snap += !is_uniform_snap ? child : 0; // Position in bitmap k_t (snapshot)

                // Check snapshot node
                if (!is_uniform_snap) {
                    is_uniform_snap = !raster_snap.t[child_pos_snap];
                    max_snap -= raster_snap.list_max[l - 1][child_pos_snap - raster_snap.m_accum_max_values[l-1]];
                }

                // Check current raster node
                if (has_children) {
                    has_children = this->k_t[child_pos];
                    if (!has_children) {
                        is_uniform = !k_eqB[child_pos - this->k_t_rank1(child_pos)];          // eqB[rank0(child_pos)]
                        diff_log = this->__decode(this->k_list_max[l - 1][child_pos - this->k_accum_max_values[l-1]]);
                    }
                }

//                std::cout << "Level: " << l << " || Pos_snap: " << child_pos_snap << "(" << max_snap << ")";
//                std::cout << " || Pos_log: " << child_pos << "(" << max_snap + diff_log << ")" << std::endl;

                // Return value
                if ((is_uniform_snap && !has_children) || is_uniform) {
                    return max_snap + diff_log;
                }

                // Go down one level on the tree.
                row = row % size;                           // Update local row
                col = col % size;                           // Update local column

                if (this->k_level_plain != 0 && l == levels) {
                    break;
                }

                // ************************************//
                // Update attributes (log raster)      //
                // ************************************//
                if (has_children) {
                    child_pos = this->get_children_position(child_pos, l);
                }

                // ************************************//
                // Update attributes (snapshot)        //
                // ************************************//
                if (!is_uniform_snap) {
                    child_pos_snap = raster_snap.get_children_position(child_pos_snap, l);
                }

                k = this->get_k(l);                         // Update k
            }
        } // END BLOCK searching from level 1 to level l - 1


        // ************************************//
        // Searching last level                //
        // ************************************//
        {
            if (!is_uniform_snap) {
                if (this->k_level_plain == 0) {
                    // TODO Check this special case
                    ushort k = raster_snap.get_k(this->k_height - 1);
                    child_pos_snap += row * k + col;
                    max_snap -= raster_snap.list_max[this->k_height - 1][child_pos_snap - raster_snap.m_accum_max_values[this->k_height-2]];
                } else{
                    if (!has_children) {
                        // Value is in plain form
                        return raster_snap.get_plain_value(row, col, max_snap, child_pos_snap, size) + diff_log;
                    }
                }
            }

            if (has_children) {
                if (this->k_level_plain == 0) {
                    // TODO Check this special case
                    ushort k = this->get_k(this->k_height - 1);
                    child_pos += row * k + col;
                    diff_log = __decode(this->k_list_max[this->k_height - 1][child_pos - this->k_accum_max_values[this->k_height-2]]);
                } else {
                    // Value is in plain form
                    diff_log = __decode(this->k_list_max[levels-1][child_pos - this->k_accum_max_values[levels-1]]);
//                    std::cout << "Level: " << levels << " || Pos_snap: " << child_pos_snap << "(" << max_snap << ")";
//                    std::cout << " || Pos_log: " << child_pos << "(" << max_snap + diff_log << ")" << std::endl;

                    return get_plain_value(row, col, max_snap + diff_log, child_pos, size);
                }
            }
            return max_snap + diff_log;
        } // END BLOCK last levels
    }

    //*****************************//
    //***** GET CELL BY VALUE *****//
    //*****************************//
    template<class k2_raster_snap>
    size_type get_cells_by_value(k2_raster_snap &raster_snap, size_type xini, size_type xend,
            size_type yini, size_type yend, value_type valini, value_type valend,
            std::vector<std::pair<size_type, size_type>> &result) {

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
                    for (auto y = yend; y <= yend; y++) {
                        result.push_back(std::pair<size_type, size_type>(x, y));
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
                                                {0, 0}, {0, 0}, !this->k_t.empty(), 0, 0, max_value,
                                                this->k_size / k, 1, result, xini, yini, (yend - yini + 1));
    }

    //*****************************//
    //**** CHECK VALUES WINDOW ****//
    //*****************************//
    template<class k2_raster_snap>
    bool check_values_window(const k2_raster_snap &raster_snap, size_type xini, size_type xend, size_type yini, size_type yend,
                             value_type valini, value_type valend, bool strong_check) const {

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
        // Encode min and max values
        {
            ushort l=0;
            this->k_accum_max_values[0] = 0;
            for (l = 0; l < max_values_.size(); l++) {
                if (!max_values_[l].empty()) {
                    this->k_list_max[l] = t_values_vec(max_values_[l]);
                    if (l != this->k_height - 1) {
                        n_nodes += max_values_[l].size();
                    }
                    this->k_accum_max_values[l+1] = this->k_accum_max_values[l] + this->k_list_max[l].size();
                } else {
                    break; // Empty levels
                }
            }
            this->k_list_max.resize(l);
            this->k_accum_max_values.resize(l+1);
            tmp_t_.resize(l);


            this->k_count_1s_k1 = 1; // 1 -> root node
            this->k_accum_min_values[0] = 0;
            for (l = 0; l < min_values_.size(); l++) {
                if (!min_values_[l].empty()) {
                    this->k_list_min[l] = t_values_vec(min_values_[l]);
                    if (l + 1 < this->k_level_k1) {
                        this->k_count_1s_k1 += min_values_[l].size();
                    }
                    this->k_accum_min_values[l+1] = this->k_accum_min_values[l] + this->k_list_min[l].size();
                } else {
                    break; // Empty levels
                }
            }
            this->k_list_min.resize(l);
            this->k_accum_min_values.resize(l+1);
        }

        // Encode bitmap T (copy temporal bitmap T)
        {
            sdsl::int_vector<1> t_(n_nodes);
            size_type n = 0;
            for (auto l = 0; l < tmp_t_.size(); l++) {
                for (auto c = 0; c < tmp_t_[l].size(); c++) {
                    t_[n++] = tmp_t_[l][c];
                }
            }
            t_.resize(n);
            this->k_t = bit_vector_type(t_);
            sdsl::util::init_support(this->k_t_rank1, &this->k_t);
        }

        // Encode bitmap eqB (copy temporal bitmap eqB)
        {
            sdsl::int_vector<1> eqB_(n_nodes);
            size_type n = 0;
            for(auto l = 0; l < tmp_eqB_.size(); l++) {
                for (auto c = 0; c < tmp_eqB_[l].size(); c++) {
                    eqB_[n++] = tmp_eqB_[l][c];
                }
            }
            eqB_.resize(n);
            this->k_eqB = bit_vector_type(eqB_);
//            sdsl::util::init_support(k_eqB_rank1, &k_eqB);
        }
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

    template<class k2_raster_snap>
    size_type get_cells_by_value_helper(
            k2_raster_snap &raster_snap, const query_info_v<size_type, value_type> &&query,
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
                    max_snap = f_max - raster_snap.list_max[l-1][pos_snap - raster_snap.m_accum_max_values[l-1]];
                } else {
                    // Parent is uniform
                    c_is_uniform_snap = true;
                    max_snap = f_max; // Get parent's value
                }

                // Check log and get max value (log)
                if (has_children) {
                    pos_log = node_log.children_pos + x * k + y;
                    diff_max = this->__decode(this->k_list_max[l-1][pos_log - this->k_accum_max_values[l-1]]);

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

//                std::cout << "Level: " << l << " || Pos_snap: " << pos_snap << "(" << max_snap << ")";
//                std::cout << " || Pos_log: " << pos_log << "(" << max_snap + diff_max << ")" << std::endl;

                if (query.valini > (max_snap + diff_max)) {
                    // Values out of range
                    continue;
                }

                // Get min value (snapshot)
                if (!c_is_uniform_snap && !c_is_uniform) {
                    ones_snap = raster_snap.t_rank1(pos_snap)+1;
                    min_snap = f_min + raster_snap.list_min[l-1][ones_snap - raster_snap.m_accum_min_values[l-1] - 1];
                } else {
                    min_snap = max_snap;
                }

                if (has_children) {
                    // Get min value (log)
                    if (c_has_children) {
                        // If it is not uniform or a leaf, its minimum value is obtained
                        ones_log = this->k_t_rank1(pos_log) + 1;
                        diff_min = this->__decode(
                                this->k_list_min[l - 1][ones_log - this->k_accum_min_values[l - 1] - 1]);
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
                    // If it is uniform, the children has not children
                    if (c_has_children || !c_is_uniform_snap) {
                        // Go down one level on the tree

                        if (has_children && (!c_has_children || is_plain_level(l-1))) {
                            // We reach a leaf in 'log'. Apply the differences permanently
                            max_snap += diff_max;
                            min_snap += diff_min;
                        }

                        if (is_plain_level(l-1)) {
                            // Get value from plain values
//                            std::cout << "Level: " << l << " || Pos_snap: " << pos_snap << "(" << max_snap << ")";
//                            std::cout << " || Pos_log: " << pos_log << "(" << max_snap << ")" << std::endl;

                            size_type ones;
                            if (c_has_children) {
                                // Check plain values at log
                                ones = (this->k_t_rank1(pos_log) + 1) - (this->k_accum_min_values[l-1]+1);
                                count_cells += get_cells_by_value_plain(c_xini, c_xend, c_yini, c_yend, query.valini, query.valend,
                                                                   children_size, ones * children_size * children_size, max_snap, result);
                            } else {
                                // Check plain values at snap
                                ones = (raster_snap.t_rank1(pos_snap) + 1) - (raster_snap.m_accum_min_values[l-1]+1);
                                count_cells += raster_snap.get_cells_by_value_plain(c_xini, c_xend, c_yini, c_yend, query.valini, query.valend,
                                                                               children_size, ones * children_size * children_size, max_snap, result);
                            } // END IF plain_values
                        } else {
                            // Go down one level
                            size_type new_children_pos = c_has_children ? this->get_children_position_ones(ones_log, l) : 0;
                            size_type new_children_pos_snap = !c_is_uniform_snap ? raster_snap.get_children_position_ones( ones_snap, l) : 0;
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
                    max_snap = f_max - raster_snap.list_max[l-1][pos_snap - raster_snap.m_accum_max_values[l-1]];
                } else {
                    // Parent is uniform
                    c_is_uniform_snap = true;
                    max_snap = f_max; // Get parent's value
                }

                // Check log and get max value (log)
                if (has_children) {
                    pos_log = node_log.children_pos + x * k + y;
                    diff_max = this->__decode(this->k_list_max[l-1][pos_log - this->k_accum_max_values[l-1]]);

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
                    if ((has_children && !c_has_children) || is_plain_level(l-1)) {
                        // We reach a leaf in 'log'. Apply the differences permanently
                        max_snap += diff_max;
                    }

                    if (is_plain_level(l-1)) {
                        // Get value from plain values
                        size_type ones;
                        if (c_has_children) {
                            // Check plain values at log
                            ones = this->k_t_rank1(pos_log) + 1;
                            ones -= (this->k_accum_min_values[l-1]+1);
                            count_cells += get_values_window_plain(c_xini, c_xend, c_yini, c_yend,
                                                                   children_size, ones * children_size * children_size, max_snap,
                                                                   result, or_x, or_y, window_size);
                        } else {
                            // Check plain values at snap
                            ones = raster_snap.t_rank1(pos_snap) + 1;
                            ones -= (raster_snap.m_accum_min_values[l-1]+1);
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
    //***** CHECK VALUES WINDOW ***//
    //*****************************//
    virtual bool check_values_window_plain(size_type xini, size_type xend, size_type yini, size_type yend,
                                           value_type valini, value_type valend,
                                           size_type size, size_type children_pos, value_type father_value,
                                           bool strong_check) const {

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

    template<class k2_raster_snap>
    bool check_values_window_helper(
            const k2_raster_snap &raster_snap, const query_info_v<size_type, value_type> &&query,
            node_info<size_type> &&node_log, node_info<size_type> &&node_snap, bool has_children,
            size_type base_x, size_type base_y, value_type f_min, value_type f_max,
            size_type children_size, ushort l,
            bool strong_check) const {

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
                    max_snap = f_max - raster_snap.list_max[l-1][pos_snap -  raster_snap.m_accum_max_values[l-1]];
                } else {
                    // Parent is uniform
                    c_is_uniform_snap = true;
                    max_snap = f_max; // Get parent's value
                }

                // Check log and get max value (log)
                if (has_children) {
                    pos_log = node_log.children_pos + x * k + y;
                    diff_max = this->__decode(this->k_list_max[l-1][pos_log - this->k_accum_max_values[l-1]]);

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

                // Check max value
                if (query.valini > (max_snap + diff_max)) {
                    // Values out of range
                    if (strong_check) return false; // (strong) At least one cell is within the range of values
                    continue;
                }

                // Get min value (snapshot)
                if (!c_is_uniform_snap && !c_is_uniform) {
                    ones_snap = raster_snap.t_rank1(pos_snap)+1;
                    min_snap = f_min + raster_snap.list_min[l-1][ones_snap -  raster_snap.m_accum_min_values[l-1] - 1];
                } else {
                    min_snap = max_snap;
                }

                // Get min value (log)
                if (has_children) {
                    if (c_has_children) {
                        // If it is not uniform or a leaf, its minimum value is obtained
                        ones_log = this->k_t_rank1(pos_log)+1;
                        diff_min = this->__decode(this->k_list_min[l-1][ones_log - this->k_accum_min_values[l-1] - 1]);
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
                        if (has_children && (!c_has_children || is_plain_level(l-1))) {
                            // We reach a leaf in 'log'. Apply the differences permanently
                            max_snap += diff_max;
                            min_snap += diff_min;
                        }

                        bool result;
                        if (is_plain_level(l-1)) {
                            // Get value from plain values
                            size_type ones;
                            if (c_has_children) {
                                // Check plain values at log
                                ones = this->k_t_rank1(pos_log) + 1;
                                ones -= (this->k_accum_min_values[l-1]+1);
                                result = check_values_window_plain(c_xini, c_xend, c_yini, c_yend, query.valini, query.valend,
                                                                       children_size, ones * children_size * children_size, max_snap, strong_check);
                            } else {
                                // Check plain values at snap
                                ones = raster_snap.t_rank1(pos_snap) + 1;
                                ones -= (raster_snap.m_accum_min_values[l-1]+1);
                                result = raster_snap.check_values_window_plain(c_xini, c_xend, c_yini, c_yend, query.valini, query.valend,
                                        children_size, ones * children_size * children_size, max_snap, strong_check);
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
    //********************* HELPERS *************************//
    //*******************************************************//
    decode_type __decode(encode_type value) const {
        return (value %2) ? (value+1)/2 : -(value/2);
    }

    encode_type __encode(decode_type value) const {
        return (value <= 0) ? -value * 2 : 2 * value - 1;
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

        // Diff parameters
        written_bytes += k_eqB.serialize(out, v, "eqB");
        return written_bytes;
    }

    void load_base(std::istream& in) {

        // k2-raster base
        k2_raster_p::load(in);

        // Plain parameters
        sdsl::read_member(k_level_k2, in);
        sdsl::read_member(k_level_plain, in);
        sdsl::read_member(k_size_leaves, in);

        // Diff parameters
        k_eqB.load(in);
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

    inline bool is_plain_level(ushort level) const {
        return this->k_level_plain > 0 && (level == (this->k_level_k1 + this->k_level_k2));
    }


}; // END CLASS t_k2_raster_heuristic_log
} // END NAMESPACE k2raster

#endif // INCLUDED_T_K2_RASTER_PLAIN_LOG
