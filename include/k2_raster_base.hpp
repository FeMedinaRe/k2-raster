/*  
 * Created by Fernando Silva on 2/07/18.
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

#ifndef INCLUDED_K2_RASTER_BASE
#define INCLUDED_K2_RASTER_BASE

#include <utils/dac_vector.hpp>


//! Namespace for k2-raster library
namespace k2raster {
    // Types
    static const ushort K2_RASTER_TYPE = 1;
    static const ushort K2_RASTER_TYPE_PLAIN = 10;
    static const ushort K2_RASTER_TYPE_HEURISTIC = 11;

    //*******************************************************//
    //******************* NODE HELPER ***********************//
    //*******************************************************//

    template <typename t_value=int>
    struct k2raster_node {
        t_value min_value;
        t_value max_value;
        ushort level;
        sdsl::int_vector<>::size_type children_pos;
    };

    template <typename value_type>
    struct node{
        value_type min_value;
        value_type max_value;
    };

    template <typename size_type>
    struct query_info {
        size_type xini;
        size_type xend;
        size_type yini;
        size_type yend;
    };

    template <typename size_type, typename value_type>
    struct query_info_v {
        size_type xini;
        size_type xend;
        size_type yini;
        size_type yend;
        value_type valini;
        value_type valend;
    };

    template <typename size_type>
    struct node_info {
        bool is_uniform;
        size_type children_pos;
    };

template < typename t_value=int,
        typename t_bv=sdsl::bit_vector,
        typename t_rank_1=typename sdsl::rank_support_v5<1,1>,
        typename t_values_vec=sdsl::dac_vector_dp_opt<t_bv, t_rank_1, 3>>
class k2_raster_base
{

    public:
        typedef t_value                         value_type;
        typedef sdsl::int_vector<>::size_type   size_type;
        typedef t_rank_1                        rank_1_type;
        typedef t_bv                            bit_vector_type;

    protected:
        /****** k2-raster Type ******/
        ushort                      k_k2_raster_type; // k2-raster Type

        /****** Size params ******/
        size_type                   k_real_size_x;  // Real row size of matrix
        size_type                   k_real_size_y;  // Real column size of matrix
        size_type                   k_size;         // Virtual size of the matrix (power of K)

        /****** Basic params ******/
        ushort                      k_k1;           // Division for first levels
        ushort                      k_k2;           // Division for last levels
        ushort                      k_level_k1;     // Number of level that apply K1 subdivision
        ushort                      k_height;       // Height of the tree

        /****** Structures ******/
        //! Bit array to store all the bits of the tree
        bit_vector_type             k_t;
        rank_1_type                 k_t_rank1;      // rank support for 1-bits on m_ta
        size_type                   k_count_1s_k1;  // Number of 1s in k_t on levels that use k1

        /****** Values ******/
        value_type                  k_max_value;        // Max value of whole matrix
        value_type                  k_min_value;        // Min value of whole matrix
        std::vector<t_values_vec>   k_list_max;         // List of max values by level
        std::vector<t_values_vec>   k_list_min;         // List of min values by level
        std::vector<size_type>      k_accum_min_values; // Accumulate of min values by level
        std::vector<size_type>      k_accum_max_values; // Accumulate of max values by level

    public:
        const value_type& max_value =                       k_max_value;
        const value_type& min_value =                       k_min_value;
        const bit_vector_type& t =                          k_t;
        const rank_1_type& t_rank1 =                        k_t_rank1;
        const std::vector<size_type>& m_accum_min_values =  k_accum_min_values;
        const std::vector<size_type>& m_accum_max_values =  k_accum_max_values;

        const ushort& k1 =                                  k_k1;
        const ushort& k2 =                                  k_k2;
        const ushort& level_k1 =                            k_level_k1;


    public:
        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//
        k2_raster_base() = default;

        k2_raster_base(const k2_raster_base& tr)
        {
            *this = tr;
        }

        k2_raster_base(k2_raster_base&& tr)
        {
            *this = std::move(tr);
        }


        //*******************************************************//
        //*************** BASIC OPERATIONS **********************//
        //*******************************************************//

        //! Move assignment operator
        k2_raster_base& operator=(k2_raster_base&& tr)
        {
            if (this != &tr) {
                k_k2_raster_type = std::move(tr.k_k2_raster_type);
                /****** Size params ******/
                k_real_size_x = std::move(tr.k_real_size_x);
                k_real_size_y = std::move(tr.k_real_size_y);
                k_size = std::move(tr.k_size);
                /****** Basic params ******/
                k_k1 = std::move(tr.k_k1);
                k_k2 = std::move(tr.k_k2);
                k_level_k1 = std::move(tr.k_level_k1);
                k_height = std::move(tr.k_height);
                /****** Structures ******/
                k_t = std::move(tr.k_t);
                k_t_rank1 = std::move(tr.k_t_rank1);
                k_t_rank1.set_vector(&k_t);
                k_count_1s_k1 = std::move(tr.k_count_1s_k1);
                /****** Values ******/
                k_max_value = std::move(tr.k_max_value);
                k_min_value = std::move(tr.k_min_value);
                k_list_max = std::move(tr.k_list_max);
                k_list_min = std::move(tr.k_list_min);
                k_accum_min_values = std::move(tr.k_acum_min_values);
                k_accum_max_values = std::move(tr.k_acum_max_values);
            }
            return *this;
        }

        //! Assignment operator
        k2_raster_base& operator=(const k2_raster_base& tr)
        {
            if (this != &tr) {
                k_k2_raster_type = tr.k_k2_raster_type;
                /****** Size params ******/
                k_real_size_x = tr.k_real_size_x;
                k_real_size_y = tr.k_real_size_y;
                k_size = tr.k_size;
                /****** Basic params ******/
                k_k1 = tr.k_k1;
                k_k2 = tr.k_k2;
                k_level_k1 = tr.k_level_k1;
                k_height = tr.k_height;
                /****** Structures ******/
                k_t = tr.k_t;
                k_t_rank1 = tr.k_t_rank1;
                k_t_rank1.set_vector(&k_t);
                k_count_1s_k1 = tr.k_count_1s_k1;
                /****** Values ******/
                k_max_value = tr.k_max_value;
                k_min_value = tr.k_min_value;
                k_list_max = tr.k_list_max;
                k_list_min = tr.k_list_min;
                k_accum_min_values = tr.k_accum_min_values;
                k_accum_max_values = tr.k_accum_max_values;
            }
            return *this;
        }

        //! Swap operator
        void swap(k2_raster_base& tr)
        {
            if (this != &tr) {
                std::swap(k_k2_raster_type, tr.k_k2_raster_type);
                /****** Size params ******/
                std::swap(k_real_size_x, tr.k_real_size_x);
                std::swap(k_real_size_y, tr.k_real_size_y);
                std::swap(k_size, tr.k_size);
                /****** Basic params ******/
                std::swap(k_k1, tr.k_k1);
                std::swap(k_k2, tr.k_k2);
                std::swap(k_level_k1, tr.k_level_k1);
                std::swap(k_height, tr.k_height);
                /****** Structures ******/
                std::swap(k_t, tr.k_t);
                sdsl::util::swap_support(k_t_rank1, tr.k_t_rank1, &k_t, &(tr.k_t));
                std::swap(k_count_1s_k1, tr.k_count_1s_k1);
                /****** Values ******/
                std::swap(k_max_value, tr.k_max_value);
                std::swap(k_min_value, tr.k_min_value);
                std::swap(k_list_max, tr.k_list_max);
                std::swap(k_list_min, tr.k_list_min);
                std::swap(k_accum_min_values, tr.k_acum_min_values);
                std::swap(k_accum_max_values, tr.k_acum_max_values);
            }
        }

        //! Equal operator
        bool operator==(const k2_raster_base& tr) const
        {
            if (k_k2_raster_type != tr.k_k2_raster_type) {
                return false;
            }
            if (k_real_size_x != tr.k_real_size_x || k_real_size_y != tr.k_real_size_y || k_size != tr.k_size) {
                return false;
            }
            if (k_k1 != tr.k_k1 || k_k2 != tr.k_k2 || k_level_k1 != tr.k_level_k1 || k_height != tr.k_height) {
                return false;
            }
            for (unsigned i = 0; i < k_t.size(); i++)
                if (k_t[i] != tr.k_t[i])
                    return false;
            if (k_count_1s_k1 != tr.k_count_1s_k1)
                return false;
            if (k_max_value != tr.k_max_value || k_min_value != tr.k_min_value) {
                return false;
            }
            for (unsigned i = 0; i < k_list_max.size(); i++)
                if (k_list_max[i].size() != tr.k_list_max[i].size())
                    return false;
            for (unsigned i = 0; i < k_list_min.size(); i++)
                if (k_list_min[i].size() != tr.k_list_min[i].size())
                    return false;
            return true;
        }


        //*******************************************************//
        //********************* GETTERS *************************//
        //*******************************************************//
        inline size_type get_n_rows() const { return k_real_size_x;}
        inline size_type get_n_cols() const { return k_real_size_y;}

        //*******************************************************//
        //****************** GETTERS - HELPERS ******************//
        //*******************************************************//

        inline ushort get_k(ushort level) const {
            return (level < k_level_k1 ? k_k1 : k_k2);
        }

        virtual inline bool is_plain_level(ushort level __attribute__((unused))) const {
            return false;
        }

        inline size_type get_children_position_ones(size_type n_ones, ushort l) const {
            if (l < this->k_level_k1) {
                return n_ones * k_k1 * k_k1;
            } else {
                return (k_count_1s_k1 * k_k1 * k_k1) +
                                        ((n_ones - k_count_1s_k1) * k_k2 * k_k2);
            }
        }

        inline size_type get_children_position(size_type child_pos, ushort l) const {
            size_type ones = k_t_rank1(child_pos) + 1;      // Number of non-empty nodes until position 'child_pos'
            return get_children_position_ones(ones, l);
        }


        //*******************************************************//
        //******************** QUERIES **************************//
        //*******************************************************//
        virtual value_type get_cell(size_type  row, size_type col) const=0;
        virtual size_type get_cells_by_value(size_type xini, size_type xend, size_type yini, size_type yend,
                value_type valini, value_type valend, std::vector<std::pair<size_type, size_type>> &result)=0;
        virtual size_type get_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                std::vector<value_type> &result)=0;
        virtual bool check_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                value_type valini, value_type valend, bool strong_check)=0;

        //*****************************//
        //******* PLAIN QUERIES *******//
        //*****************************//
        virtual size_type get_cells_by_value_plain(size_type xini __attribute__((unused)), size_type xend  __attribute__((unused)),
                size_type yini  __attribute__((unused)), size_type yend __attribute__((unused)),
                value_type valini  __attribute__((unused)), value_type valend  __attribute__((unused)),
                size_type size  __attribute__((unused)), size_type children_pos  __attribute__((unused)),
                value_type father_value  __attribute__((unused)), std::vector<std::pair<size_type, size_type>> &result  __attribute__((unused))) {
            return 0;
        }
        virtual size_type get_values_window_plain(size_type xini __attribute__((unused)), size_type xend __attribute__((unused)),
                size_type yini __attribute__((unused)), size_type yend __attribute__((unused)),
                size_type size __attribute__((unused)), size_type children_pos __attribute__((unused)),
                value_type father_value __attribute__((unused)), std::vector<value_type> &result __attribute__((unused)),
                size_type or_x __attribute__((unused)), size_type or_y __attribute__((unused)),
                size_type window_size __attribute__((unused))) {
            return 0;
        }
        virtual bool check_values_window_plain(size_type xini __attribute__((unused)), size_type xend __attribute__((unused)),
                size_type yini __attribute__((unused)), size_type yend __attribute__((unused)),
                value_type valini __attribute__((unused)), value_type valend __attribute__((unused)),
                size_type size __attribute__((unused)), size_type children_pos __attribute__((unused)),
                value_type father_value __attribute__((unused)), bool strong_check __attribute__((unused))) {
            return false;
        }

        //*******************************************************//
        //********************** FILE ***************************//
        //*******************************************************//
        virtual size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {

            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;

            /****** k2-raster Type ******/
            written_bytes += write_member(k_k2_raster_type, out, child, "k2-raster_type");

            /****** Size params ******/
            written_bytes += write_member(k_real_size_x, out, child, "real_size_x");
            written_bytes += write_member(k_real_size_y, out, child, "real_size_y");
            written_bytes += write_member(k_size, out, child, "size");

            /****** Basic params ******/
            written_bytes += write_member(k_k1, out, child, "k1");
            written_bytes += write_member(k_k2, out, child, "k2");
            written_bytes += write_member(k_level_k1, out, child, "level_k1");
            written_bytes += write_member(k_height, out, child, "height");

            /****** Structures ******/
            written_bytes += k_t.serialize(out, child, "t");
            written_bytes += k_t_rank1.serialize(out, child, "t_rank1");
            written_bytes += write_member(k_count_1s_k1, out, child, "count_1s_k1");


            /****** Values ******/
            written_bytes += write_member(k_max_value, out, child, "max_value");
            written_bytes += write_member(k_min_value, out, child, "min_value");

            // Values
            written_bytes += sdsl::serialize(k_list_max, out, child, "LMax");
            written_bytes += sdsl::serialize(k_list_min, out, child, "LMin");
            written_bytes += sdsl::serialize(k_accum_min_values, out, child, "acum_min_values");
            written_bytes += sdsl::serialize(k_accum_max_values, out, child, "acum_max_values");

            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        virtual void load(std::istream& in) {
            /****** k2-raster Type ******/
            ushort k2_raster_type;
            sdsl::read_member(k2_raster_type, in);

            /****** Size params ******/
            sdsl::read_member(k_real_size_x, in);
            sdsl::read_member(k_real_size_y, in);
            sdsl::read_member(k_size, in);

            /****** Basic params ******/
            sdsl::read_member(k_k1, in);
            sdsl::read_member(k_k2, in);
            sdsl::read_member(k_level_k1, in);
            sdsl::read_member(k_height, in);

            /****** Structures ******/
            k_t.load(in);
            k_t_rank1.load(in);
            k_t_rank1.set_vector(&k_t);
            sdsl::read_member(k_count_1s_k1, in);

            /****** Values ******/
            sdsl::read_member(k_max_value, in);
            sdsl::read_member(k_min_value, in);

            // Values
            sdsl::load(k_list_max, in);
            sdsl::load(k_list_min, in);
            sdsl::load(k_accum_min_values, in);
            sdsl::load(k_accum_max_values, in);
        }


        //*******************************************************//
        //********************** TEST ***************************//
        //*******************************************************//
        template<class Container>
        bool check(Container&& values, size_type n_rows, size_type n_cols) {
            for (auto r = 0; r < n_rows; r++) {
                for (auto c = 0; c < n_cols; c++) {
                    size_type pos = r * n_cols + c;
                    value_type v1 = get_cell(r, c);
                    value_type v2 = values[pos];
                    if ( v1 != v2) {
                        std::cout << "Error position (" << r << ", " << c << "): get " << v1 << " and expected " << v2 << std::endl;
                        return false;
                    }
                }
            }
            return true;
        }

    //*******************************************************//
    //******************** NAVIGATION ***********************//
    //*******************************************************//
    k2raster_node<value_type> get_root() const {
        return {k_min_value, k_max_value, 0, 0};
    }

    virtual k2raster_node<value_type> get_child(const k2raster_node<value_type> parent, uint child) const {
        if (this->k_t.empty() || parent.min_value == parent.max_value|| parent.level == this->k_height) {
            // No children
            return parent;
        }

        k2raster_node<value_type> child_node = parent;
        child_node.level += 1;
        size_type child_pos = parent.children_pos + child;

        // Get max value
        child_node.max_value -= get_max_value(child_node.level, child_pos);

        if (child_node.level == this->k_height) {
            child_node.min_value = child_node.max_value;
            return child_node;
        }

        // Check if all the cells in the subarray are equal (uniform matrix)
        if (!this->k_t[child_pos]) {
            child_node.min_value = child_node.max_value;
        } else {
            size_type ones = k_t_rank1(child_pos) + 1;      // Number of non-empty nodes until position 'child_pos'
            child_node.min_value +=  get_min_value(child_node.level, ones-1);
            set_children_position(child_node, ones);
        }

        return child_node;
    }

    virtual inline void set_children_position(k2raster_node<value_type> &child_node, size_type ones ) const{
        child_node.children_pos = this->get_children_position_ones(ones, child_node.level);
    }

    inline value_type get_max_value(ushort level, size_type pos) const {
        return k_list_max[level-1][pos - this->k_accum_max_values[level-1]];
    }

    inline value_type get_min_value(ushort level, size_type ones) const {
        return k_list_min[level-1][ones - this->k_accum_min_values[level-1]];
    }


    inline value_type get_min_value(ushort level, size_type pos, size_type &ones) const {
            ones = this->k_t_rank1(pos);
            return k_list_min[level-1][ones - this->k_accum_min_values[level-1]];
    }

    protected:
        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//

        k2_raster_base(size_type n_rows, size_type n_cols,
                       ushort k1, ushort k2, ushort level_k1, ushort k2_raster_type) {

            // Set type
            this->k_k2_raster_type = k2_raster_type;

            // Set K values and store real sizes
            k_real_size_x = n_rows;
            k_real_size_y = n_cols;
            k_k1 = k1;
            k_k2 = k2;
            k_level_k1 = level_k1;

            // Initialize levels with new values of K
            init_levels();

            // Resize arrays of values
            k_list_max.resize(k_height);
            k_list_min.resize(k_height);
            k_accum_max_values.resize(k_height+1); //TODO change to k_height
            k_accum_min_values.resize(k_height+1); //TODO change to k_height
        }


        //*******************************************************//
        //******************** HELPERS   ************************//
        //*******************************************************//


        void init_levels(){
            // Calculate size, number of levels and other important params
            size_type size = std::max(k_real_size_x, k_real_size_y);
            k_height = std::ceil(std::log(size)/std::log(std::min(k_k1, k_k2))); // TODO fix this for k1 and k2
            k_size = pow(std::min(k_k1, k_k2), k_height);

            size_type sub_size = k_size;
            ushort l = 0;
            while (sub_size != 1) {
                sub_size /= get_k(l);
                l++;
                if (sub_size == 0) {
                    std::cout << "Error: k1 = " << k_k1 << " and k2 = " << k_k2 << " are not valid params";
                    std::cout << " with size " << k_real_size_x << "x" << k_real_size_x << " (" << k_size << "x" << k_size << ")" << std::endl;
                    exit(-1);
                }
            }

            // Update height of the tree
            k_height = l;
            k_level_k1 = k_level_k1 > l ? l : k_level_k1;
        }

        //*******************************************************//
        //*********************** BUILD *************************//
        //*******************************************************//

        size_type build_max_values(const std::vector<std::vector<value_type>> &max_values_,
                std::vector<sdsl::int_vector<1>> &tmp_t_) {
            ushort l;
            size_type n_nodes = 0;
            this->k_list_max.resize(max_values_.size());
            this->k_accum_max_values.resize(max_values_.size()+1);

            for (l = 0; l < max_values_.size(); l++) {
                if (!max_values_[l].empty()) {
                    this->k_list_max[l] = t_values_vec(max_values_[l]);
                    if (l != this->k_height - 1) {
                        n_nodes += max_values_[l].size();
                    }
                    this->k_accum_max_values[l+1] = this->k_accum_max_values[l] + this->k_list_max[l].size();
                } else{
                    break; // Empty levels
                }
            }
            this->k_list_max.resize(l);
            this->k_accum_max_values.resize(l+1);
            tmp_t_.resize(l);

            return n_nodes;
        }

        void build_min_values(const std::vector<std::vector<value_type>> &min_values_){
            ushort l;
            this->k_list_min.resize(min_values_.size());
            this->k_accum_min_values.resize(min_values_.size()+1);

            this->k_count_1s_k1 = 1; // 1 -> root node
            this->k_accum_min_values[0] = 0;
            for (l = 0; l < min_values_.size(); l++) {
                if (!min_values_[l].empty()) {
                    this->k_list_min[l] = t_values_vec(min_values_[l]);
                    if (l+1 < this->k_level_k1) {
                        this->k_count_1s_k1 += min_values_[l].size();
                    }
                    this->k_accum_min_values[l+1] = this->k_accum_min_values[l] + this->k_list_min[l].size();
                } else{
                    break; // Empty levels
                }
            }
            this->k_list_min.resize(l);
            this->k_accum_min_values.resize(l+1);
        }

        void build_t(std::vector<sdsl::int_vector<1>> tmp_t_, size_type n_nodes){
            sdsl::int_vector<1> t_(n_nodes);
            size_type n = 0;
            ushort l;
            for (l = 0; l < tmp_t_.size(); l++) {
                for(auto c = 0; c < tmp_t_[l].size(); c++) {
                    t_[n++] = tmp_t_[l][c];
                }
            }
            t_.resize(n);
            this->k_t = bit_vector_type(t_);
            sdsl::util::init_support(this->k_t_rank1, &this->k_t);
        }

};
} // END NAMESPACE k2raster


#endif // INCLUDED_K2_RASTER_BASE
