/*  
 * Created by Fernando Silva on 5/06/19.
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

#ifndef INCLUDED_K2_RASTER_HEURISTIC_LOG
#define INCLUDED_K2_RASTER_HEURISTIC_LOG

#include <k2_raster_heuristic.hpp>

//! Namespace for the k2-raster library
namespace k2raster {

    template<typename t_value=int,
            typename t_bv=sdsl::bit_vector,
            typename t_rank_1=sdsl::rank_support_v5<1,1>,
            typename t_values_vec=sdsl::dac_vector_dp_opt<t_bv, t_rank_1, 3>,
            typename t_values_last=sdsl::dac_vector_dp_opt<t_bv, t_rank_1, 3>>
    class k2_raster_heuristic_log : public k2_raster_heuristic<t_value, t_bv, t_rank_1, t_values_vec, t_values_last> {

    public:
        typedef k2_raster_base <t_value, t_bv, t_rank_1, t_values_vec>  k2_raster_p_p; // k2_raster parent
        typedef k2_raster_heuristic <t_value, t_bv, t_rank_1, t_values_vec, t_values_last>  k2_raster_p; // k2_raster parent
        typedef k2_raster_plain <t_value, t_bv, t_rank_1, t_values_vec, std::vector<t_value>>  k2_raster_p_t; // k2_raster parent
        typedef t_value                                         value_type;
        typedef k2_raster_base<>::size_type                     size_type;
        typedef t_rank_1                                        rank_1_type;
        typedef uint                                            codeword_type;
        typedef t_bv                                            bit_vector_type;

    protected:
        std::vector<char> *k_voc_p;              // Array that stores vocabulary words

    public:

        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//
        k2_raster_heuristic_log() = default;

        k2_raster_heuristic_log(const k2_raster_heuristic_log &tr) : k2_raster_p()  {
            *this = tr;
        }

        k2_raster_heuristic_log(k2_raster_heuristic_log &&tr) {
            *this = std::move(tr);
        }

        // k2_raster_plain to k2_raster_heuristic
        k2_raster_heuristic_log(const k2_raster_p_t &k2_raster_t, std::unordered_map<std::string, size_type> &hash_words,
                size_type size_word, size_type n_words, std::vector<char> &k_voc) {

            k2_raster_p_p::operator=(k2_raster_t);
            k2_raster_t.copy(this->k_level_k2, this->k_level_plain, this->k_size_leaves);
            this->k_size_word = size_word;
            this->k_n_words = n_words;
            k_voc_p = &k_voc;

            /******************************************************************************/
            /* SIXTH STEP - Compact submatrices and encode structures                     */
            /******************************************************************************/
            {
                auto sub_size = this->m_size_leaves;         // Size of each (leaf) submatrix
                auto n_submatrices = k2_raster_t.m_plain_values.size() / sub_size;
                char * plain_words = (char *)k2_raster_t.m_plain_values.data();

                // Temporal structures
                std::vector<codeword_type> encoded_values_;
                std::vector<value_type> no_compacted_values_;
                bit_vector_type in_voc_(n_submatrices, 0);

                // Check each submatrix
                codeword_type codeword;
                size_type n_c = 0, n_p = 0;
                for (size_type m = 0; m < n_submatrices; m++) {
                    // Get word
                    std::string word(&(plain_words[m * this->k_size_word]), this->k_size_word);
                    codeword = hash_words[word];
                    if (codeword != 0) {
                        // Use dictionary
                        encoded_values_.push_back(codeword-1);
                        n_c++;
                        in_voc_[m] = 1;
                    } else {
                        // Use plain form
                        for (auto c = 0; c < sub_size; c++) {
                            no_compacted_values_.push_back(k2_raster_t.m_plain_values[m * sub_size + c]); // copy values
                            n_p++;
                        }
                        in_voc_[m] = 0;
                    }
                }

                // Encode structures
                this->k_is_in_voc = bit_vector_type(in_voc_);
                sdsl::util::init_support(this->k_is_in_voc_rank1, &this->k_is_in_voc);
                this->k_encoded_values = t_values_last(encoded_values_);
                this->k_plain_values = t_values_last(no_compacted_values_);
            } // END SIXTH STEP
        }

        //*******************************************************//
        //*************** BASIC OPERATIONS **********************//
        //*******************************************************//

        //! Move assignment operator
        k2_raster_heuristic_log &operator=(k2_raster_heuristic_log &&tr) {
            if (this != &tr) {
                k2_raster_p::operator=(tr);
                k_voc_p = std::move(tr.k_voc_p);
            }
            return *this;
        }

        //! Assignment operator
        k2_raster_heuristic_log &operator=(const k2_raster_heuristic_log &tr) {
            if (this != &tr) {
                k2_raster_p::operator=(tr);
                k_voc_p = tr.k_voc_p;
            }
            return *this;
        }

        //! Swap operator
        void swap(k2_raster_heuristic_log &tr) {
            if (this != &tr) {
                k2_raster_p::swap(tr);
                std::swap(k_voc_p, tr.k_voc_p);
            }
        }

        //! Equal operator
        bool operator==(const k2_raster_heuristic_log &tr) const {
            if (!k2_raster_p::operator==(tr)) {
                return false;
            }

            return true;
        }

        //*******************************************************//
        //********************** FILE ***************************//
        //*******************************************************//
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;

            // Seralize k2-raster
            written_bytes += this->serialize_without_dict(out, child, name);

            // Vocabulary
            written_bytes += write_member(this->k_size_word, out, child, "size_word");
            written_bytes += write_member(this->k_n_words, out, child, "n_words");

            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            this->load_base(in);

            // Submatrices encoded
            this->k_is_in_voc.load(in);
            this->k_is_in_voc_rank1.load(in);
            this->k_is_in_voc_rank1.set_vector(&this->k_is_in_voc);
            this->k_encoded_values.load(in);
            this->k_plain_values.load(in);

            // Vocabulary
            sdsl::read_member(this->k_size_word, in);
            sdsl::read_member(this->k_n_words, in);
        }

        void load(std::istream& in, std::vector<char> &k_voc) {
            load(in);
            k_voc_p = &k_voc;
        }


    protected:
        //*******************************************************//
        //**************** QUERIES HELPERS **********************//
        //*******************************************************//

        //*****************************//
        //***** GET CELL          *****//
        //*****************************//
        virtual value_type get_cell_plain(size_type row, size_type col, size_type size,
                                          size_type children_pos, value_type father_value) const {
            row = row % size;
            col = col % size;
            size_type array_pos = row * size + col;
            children_pos = children_pos / size / size;

            // Check in bitmap  submatrix is in vocabulary or not
            if (this->k_is_in_voc[children_pos]) {
                size_type pos = this->k_is_in_voc_rank1(children_pos);
                value_type codeword = this->k_encoded_values[pos];
                const char * word = &(k_voc_p->data()[codeword * this->k_size_word]);
                value_type *values = (value_type *)word;
                return father_value - values[array_pos];
            } else {
                size_type pos = (children_pos - this->k_is_in_voc_rank1(children_pos)) * size * size;
                return father_value - this->k_plain_values[pos + array_pos];
            }
        }


    }; // END class k2_raster_heuristic
} // END namespace sdsl


#endif // INCLUDED_K2_RASTER_HEURISTIC
