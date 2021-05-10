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

#ifndef INCLUDED_K2_RASTER_HEURISTIC
#define INCLUDED_K2_RASTER_HEURISTIC

#include <k2_raster_plain.hpp>
#include <utils/utils_heuristic.hpp>

//! Namespace for the k2-raster library
namespace k2raster {

    template<typename t_value=int,
            typename t_bv=sdsl::bit_vector,
            typename t_rank_1=sdsl::rank_support_v5<1,1>,
            typename t_values_vec=sdsl::dac_vector_dp_opt<t_bv, t_rank_1, 3>,
            typename t_values_last=sdsl::dac_vector_dp_opt<t_bv, t_rank_1, 3>>
    class k2_raster_heuristic : public k2_raster_plain<t_value, t_bv, t_rank_1, t_values_vec, t_values_last> {

    public:
        typedef k2_raster_plain <t_value, t_bv, t_rank_1, t_values_vec, t_values_last>  k2_raster_p; // k2_raster parent
        typedef t_value                                         value_type;
        typedef k2_raster_base<>::size_type                     size_type;
        typedef t_rank_1                                        rank_1_type;
        typedef uint                                            codeword_type;
        typedef t_bv                                            bit_vector_type;

    protected:

        // Vocabulary
        ushort                      k_size_word=0;      // Size (in bytes) of one word of the vocabulary
        size_type                   k_n_words=0;        // Number of words in the vocabulary
        std::vector<char>           k_voc;              // Array that stores vocabulary words

        // Submatrices encoded
        bit_vector_type             k_is_in_voc;        // Submatrix encode as a entry of the vocabulary or as plain form
        rank_1_type                 k_is_in_voc_rank1;  // rank support for 1-bits on k_is_in_voc

        t_values_last               k_encoded_values;   // Encoded vocabulary of leaves

    public:

        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//
        k2_raster_heuristic() = default;

        k2_raster_heuristic(const k2_raster_heuristic &tr) : k2_raster_p() {
            *this = tr;
        }

        k2_raster_heuristic(k2_raster_heuristic &&tr) {
            *this = std::move(tr);
        }

        template<class Container>
        k2_raster_heuristic(Container &&c_values, size_type n_rows, size_type n_cols, ushort k1, ushort k2, ushort level_k1, ushort plains_levels)
                : k2_raster_p(c_values, n_rows, n_cols, k1, k2, level_k1) {

            this->k_k2_raster_type = K2_RASTER_TYPE_HEURISTIC;
            this->k_level_plain = plains_levels;
            this->init_levels();

            std::vector<value_type> plain_values_;
            this->build(c_values, plain_values_);
            encoded_plain_values(plain_values_);
        }

        //*******************************************************//
        //*************** BASIC OPERATIONS **********************//
        //*******************************************************//

        //! Move assignment operator
        k2_raster_heuristic &operator=(k2_raster_heuristic &&tr) {
            if (this != &tr) {
                k2_raster_p::operator=(tr);

                // Vocabulary
                k_size_word = std::move(tr.k_size_word);
                k_n_words = std::move(tr.k_n_words);
                k_voc = std::move(tr.k_voc);

                // Submatrices encoded
                k_is_in_voc = std::move(tr.k_is_in_voc);
                k_is_in_voc_rank1 = std::move(tr.k_is_in_voc_rank1);
                k_is_in_voc_rank1.set_vector(&k_is_in_voc);
                k_encoded_values = std::move(tr.k_encoded_values);
            }
            return *this;
        }

        //! Assignment operator
        k2_raster_heuristic &operator=(const k2_raster_heuristic &tr) {
            if (this != &tr) {
                k2_raster_p::operator=(tr);

                // Vocabulary
                k_size_word = tr.k_size_word;
                k_n_words = tr.k_n_words;
                k_voc = tr.k_voc;

                // Submatrices encoded
                k_is_in_voc = tr.k_is_in_voc;
                k_is_in_voc_rank1 = tr.k_is_in_voc_rank1;
                k_is_in_voc_rank1.set_vector(&k_is_in_voc);
                k_encoded_values = tr.k_encoded_values;
            }
            return *this;
        }

        //! Swap operator
        void swap(k2_raster_heuristic &tr) {
            if (this != &tr) {
                k2_raster_p::swap(tr);

                // Vocabulary
                std::swap(k_size_word, tr.k_size_word);
                std::swap(k_n_words, tr.k_n_words);
                std::swap(k_voc, tr.k_voc);

                // Submatrices encoded
                std::swap(k_is_in_voc, tr.k_is_in_voc);
                sdsl::util::swap_support(k_is_in_voc_rank1, tr.k_is_in_voc_rank1, &k_is_in_voc, &(tr.k_is_in_voc));
                k_encoded_values.swap(tr.k_encoded_values);
            }
        }

        //! Equal operator
        bool operator==(const k2_raster_heuristic &tr) const {
            if (!k2_raster_p::operator==(tr)) {
                return false;
            }

            if (!k2_raster_p::operator==(tr)) {
                return false;
            }

            // Vocabulary
            if (k_size_word != tr.k_size_word || k_n_words != tr.k_n_words || k_voc.size() != tr.k_voc.size()) {
                return false;
            }

            // Submatrices encoded
            return !(k_is_in_voc.size() != tr.k_is_in_voc.size() || k_encoded_values.size() != tr.k_encoded_values.size());
        }

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
            children_pos = children_pos / this->k_size_leaves;

            // Check in bitmap  submatrix is in vocabulary or not
            if (k_is_in_voc[children_pos]) {
                size_type pos = k_is_in_voc_rank1(children_pos);
                value_type codeword = k_encoded_values[pos];
                const char * word = &(k_voc.data()[codeword * k_size_word]);
                auto *values = (value_type *)word;
                return father_value - values[array_pos];
            } else {
                size_type pos = (children_pos - k_is_in_voc_rank1(children_pos)) * this->k_size_leaves;
                return father_value - this->k_plain_values[pos + array_pos];
            }
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
            children_pos = children_pos / this->k_size_leaves;

            // Check in the submatrix is in vocabulary or not
            if (k_is_in_voc[children_pos]) {
                // Values in VOC
                size_type pos = k_is_in_voc_rank1(children_pos);
                value_type codeword = k_encoded_values[pos];
                const char * word = &(k_voc.data()[codeword * k_size_word]);
                auto *values = (value_type *)word;
                for (auto x = xini; x <= xend; x++) {
                    size_type array_pos = ((x % size) * size + (yini % size));
                    for (auto y = yini; y <= yend; y++) {
                        value = father_value - values[array_pos];
                        if (value >= valini && value <= valend) {
                            result.emplace_back(x, y);
                            count_cells++;
                        }
                        array_pos++; // Position of the next value
                    } // END FOR y
                } // END FOR x
            } else {
                // Values in plain
                size_type pos = (children_pos - k_is_in_voc_rank1(children_pos)) * this->k_size_leaves;
                for (auto x = xini; x <= xend; x++) {
                    size_type array_pos = pos + ((x % size) * size + (yini % size));
                    for (auto y = yini; y <= yend; y++) {
                        value = father_value - this->k_plain_values.access(array_pos);
                        if (value >= valini && value <= valend) {
                            result.emplace_back(x, y);
                            count_cells++;
                        }
                        array_pos++; // Position of the next value
                    } // END FOR y
                } // END FOR x
            } // END IF is_in_voc
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
            children_pos = children_pos / this->k_size_leaves;

            // Check in the submatrix is in vocabulary or not
            if (k_is_in_voc[children_pos]) {
                // Values in VOC
                size_type pos = k_is_in_voc_rank1(children_pos);
                value_type codeword = k_encoded_values[pos];
                const char * word = &(k_voc.data()[codeword * k_size_word]);
                auto *values = (value_type *)word;
                for (auto x = xini; x <= xend; x++) {
                    size_type array_pos = ((x % size) * size + (yini % size));
                    cell_pos = (x - or_x) * window_size + (yini - or_y);
                    for (auto y = yini; y <= yend; y++) {
                        value = father_value - values[array_pos];
                        result[cell_pos++] = value;
                        count_cells++;
                        array_pos++; // Position of the next value
                    } // END FOR y
                } // END FOR x
            } else {
                // Values in plain
                size_type pos = (children_pos - k_is_in_voc_rank1(children_pos)) * this->k_size_leaves;
                for (auto x = xini; x <= xend; x++) {
                    size_type array_pos = pos + ((x % size) * size + (yini % size));
                    cell_pos = (x - or_x) * window_size + (yini - or_y);
                    for (auto y = yini; y <= yend; y++) {
                        value = father_value - this->k_plain_values.access(array_pos);
                        result[cell_pos++] = value;
                        count_cells++;
                        array_pos++; // Position of the next value
                    } // END FOR y
                } // END FOR x
            } // END IF is_in_voc
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
            children_pos = children_pos / this->k_size_leaves;

            // Check in the submatrix is in vocabulary or not
            if (k_is_in_voc[children_pos]) {
                // Values in VOC
                size_type pos = k_is_in_voc_rank1(children_pos);
                value_type codeword = k_encoded_values[pos];
                const char * word = &(k_voc.data()[codeword * k_size_word]);
                auto *values = (value_type *)word;
                for (auto x = xini; x <= xend; x++) {
                    size_type array_pos = ((x % size) * size + (yini % size));
                    for (auto y = yini; y <= yend; y++) {
                        value = father_value - values[array_pos];
                        if (value >= valini && value <= valend) {
                            if (!strong_check) return true;
                        } else {
                            if (strong_check) return false;
                        }
                        array_pos++; // Position of the next value
                    } // END FOR y
                } // END FOR x
            } else {
                // Values in plain
                size_type pos = (children_pos - k_is_in_voc_rank1(children_pos)) * size * size;
                for (auto x = xini; x <= xend; x++) {
                    size_type array_pos = pos + ((x % size) * size + (yini % size));
                    for (auto y = yini; y <= yend; y++) {
                        value = father_value - this->k_plain_values.access(array_pos);
                        if (value >= valini && value <= valend) {
                            if (!strong_check) return true;
                        } else {
                            if (strong_check) return false;
                        }
                        array_pos++; // Position of the next value
                    } // END FOR y
                } // END FOR x
            } // END IF is_in_voc
            return strong_check;
        }

     //*******************************************************//
        //********************** FILE ***************************//
        //*******************************************************//
        size_type serialize_without_dict(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {

            size_type written_bytes = 0;
            written_bytes += this->serialize_base(out, v, name);

            // Submatrices encoded
            written_bytes += k_is_in_voc.serialize(out, v, "is_in_voc");
            written_bytes += k_is_in_voc_rank1.serialize(out, v, "is_in_voc_rank1");
            written_bytes += k_encoded_values.serialize(out, v, "encoded_values");
            written_bytes += this->k_plain_values.serialize(out, v, "plain_values");

            return written_bytes;
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;

            // Seralize k2-raster
            written_bytes += serialize_without_dict(out, child, name);

            // Vocabulary
            written_bytes += write_member(k_size_word, out, child, "size_word");
            written_bytes += write_member(k_n_words, out, child, "n_words");
            written_bytes += serialize_vector(k_voc, out, child, "voc");

            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            this->load_base(in);

            // Submatrices encoded
            k_is_in_voc.load(in);
            k_is_in_voc_rank1.load(in);
            k_is_in_voc_rank1.set_vector(&k_is_in_voc);
            k_encoded_values.load(in);
            this->k_plain_values.load(in);

            // Vocabulary
            sdsl::read_member(k_size_word, in);
            sdsl::read_member(k_n_words, in);
            k_voc.resize(k_n_words * k_size_word);
            sdsl::load_vector(k_voc, in);
        }

        //*******************************************************//
        //********************** INFO ***************************//
        //*******************************************************//
        virtual void print_info() {
            this->print_info_base();

            // k_is_in_voc
            std::cout << "Bitmap (k_is_in_voc): ";
            for (size_t b = 0; b < this->k_is_in_voc.size(); b++) {
                std::cout << this->k_is_in_voc[b];
            }
            std::cout << "(\t" << this->k_is_in_voc.size() << " bits)" << std::endl;
            std::cout << std::endl;

            // Encodes values
            std::cout << "Encoded values: ";
            for (size_t p = 0; p < k_encoded_values.size(); p++) {
                value_type value = k_encoded_values[p];
                std::cout << value << " ";
            }
            std::cout << "\t(" << k_encoded_values.size() << " values)" << std::endl;
            std::cout << std::endl;

            // Encodes values
            std::cout << "Plain values:";
            for (size_t v = 0; v < this->k_plain_values.size(); v++) {
                value_type value = this->k_plain_values.access(v);
                std::cout << value << " ";
            }
            std::cout << "\t(" << this->k_plain_values.size() << " values)" << std::endl;
            std::cout << std::endl;

            // Vocabulary
            std::cout << "Vocabulary:" << std::endl;
            std::cout << "Number of entries: " << k_n_words << " || Length of word: " << k_size_word << std::endl;
            auto* tmp =  (value_type*)k_voc.data();
            size_type tmp_size =k_size_word / sizeof(value_type);
            for (size_type w = 0; w < k_n_words; w++) {
                std::cout << "Entry " << w << ": ";
                for (size_type p = 0; p < tmp_size; p++){
                    std::cout << tmp[w * tmp_size + p] << " ";
                }
                std::cout << std::endl;
            }
        }

        //*******************************************************//
        //******************** NAVIGATION ***********************//
        //*******************************************************//
        virtual std::vector<value_type> get_plain_values(const k2raster_node<value_type> &parent) const {
            if (parent.min_value == parent.max_value) {
                return std::vector<value_type>(this->k_size_leaves, parent.max_value);
            }

            // Calculate children position
            size_type children_pos = this->k_t_rank1(parent.children_pos) + 1;      // Number of non-empty nodes until position 'child_pos'
            children_pos = children_pos - (this->k_accum_min_values[parent.level-1] + 1);

            std::vector<value_type> result(this->k_size_leaves);
            // Check in the submatrix is in vocabulary or not
            if (k_is_in_voc[children_pos]) {
                // Values in VOC
                size_type pos = k_is_in_voc_rank1(children_pos);
                value_type codeword = k_encoded_values[pos];
                const char *word = &(k_voc.data()[codeword * k_size_word]);
                auto *values = (value_type *) word;
                for (auto c = 0; c < this->k_size_leaves; c++) {
                    result[c] = parent.max_value - values[c];
                }
            } else {
                // Values in plain
                size_type pos = (children_pos - k_is_in_voc_rank1(children_pos)) * this->k_size_leaves;
                for (auto c = 0; c < this->k_size_leaves; c++) {
                    result[c] = parent.max_value - this->k_plain_values[pos++];
                }
            } // END IF k_is_in_voc
            return result;
        }

    protected:

        //*******************************************************//
        //********************** BUILD **************************//
        //*******************************************************//
        virtual void encoded_plain_values(std::vector<value_type> &plain_values_) {
            if (!plain_values_.empty()) {

                /******************************************************************************/
                /* ZERO STEP - Set some params                                                */
                /******************************************************************************/
                // This steps translate numbers into char values.
                // This is because we need to manage values ​​as words (strings), making it easier to use the hash table.
                // A word consists of the fusion of all values ​​of the submatrix
                char * plain_words = (char *)plain_values_.data();
                auto sub_size = this->m_size_leaves;                  // Size of each (leaf) submatrix
                k_size_word = sub_size * sizeof(value_type);            // Size, in bytes, of each submatrix
                size_t n_submatrices = plain_values_.size() / sub_size;

                /******************************************************************************/
                /* FIRST STEP - Calculate frequency of each "word" and each different "value" */
                /******************************************************************************/
                std::unordered_map<std::string, size_type> hash_words;   // Stores the frequency of each word
                std::unordered_map<value_type , size_type> hash_values;  // Stores the frequency of each value
                {
                    k2raster_heuristic::create_hashes<value_type, size_type>(hash_words, hash_values, n_submatrices, sub_size,
                                                                         k_size_word, plain_values_.data());
#ifndef NDEBUG
                    std::cout << "Found " << hash_words.size() << " words in leaves values and " << hash_values.size() << " values";
                    std::cout << " in " << n_submatrices << " submatrices" << std::endl;
#endif
                } // END FIRST STEP

                /******************************************************************************/
                /* SECOND STEP - Calculate entropy (vocabulary and values)                    */
                /******************************************************************************/
                double entropy_values = 0, entropy_words = 0;
                {
                    k2raster_heuristic::calculate_entropy<value_type, size_type>(hash_words, n_submatrices,
                                                                             hash_values, n_submatrices * sub_size,
                                                                             entropy_words, entropy_values);
                    hash_values.clear();
#ifndef NDEBUG
                    std::cout << "Entropy Words: " << entropy_words << " and Entropy Values " << entropy_values << std::endl;
#endif
                } // END SECOND STEP

                /******************************************************************************/
                /* THIRD STEP - Sort vocabulary by frequency                                  */
                /******************************************************************************/
                std::vector<std::string> sort_words = k2raster_heuristic::sort_vocabulary<size_type>(hash_words);

                /******************************************************************************/
                /* FOURTH STEP - Assign codewords to vocabulary                               */
                /******************************************************************************/
                size_type n_compacted = 0, n_no_compacted = 0;
                {
                    k_n_words = k2raster_heuristic::assign_codewords<value_type, size_type, codeword_type>(sort_words, hash_words, sub_size,
                                                                                                       entropy_words, entropy_values,
                                                                                                       n_compacted, n_no_compacted);
#ifndef NDEBUG
                    std::cout << "Created dictionary with: " << k_n_words << " words " << std::endl;
#endif
                } // END FOURTH STEP

                /******************************************************************************/
                /* FIFTH STEP - Copy words                                                    */
                /******************************************************************************/
                {
                    k2raster_heuristic::copy_words<size_type, codeword_type>(k_voc, sort_words, hash_words, k_n_words, k_size_word);
                    sort_words.clear();
                    sort_words.shrink_to_fit();
                } // END STEP FIFTH

                /******************************************************************************/
                /* SIXTH STEP - Compact submatrices and encode structures                     */
                /******************************************************************************/
                {
                    // Temporal structures
                    sdsl::int_vector<sizeof(codeword_type) * 8> encoded_values_(n_compacted);
                    sdsl::int_vector<sizeof(value_type) * 8> no_compacted_values_(n_no_compacted * sub_size);
                    bit_vector_type in_voc_(n_submatrices, 0);

                    // Check each submatrix
                    codeword_type codeword;
                    size_type n_c = 0, n_p = 0;
                    for (size_t m = 0; m < n_submatrices; m++) {
                        // Get word
                        std::string word(&(plain_words[m * k_size_word]), k_size_word);
                        codeword = hash_words[word];
                        if (codeword != 0) {
                            // Use dictionary
                            encoded_values_[n_c++] = codeword-1;
                            in_voc_[m] = 1;
                        } else {
                            // Use plain form
                            for (auto c = 0; c < sub_size; c++) {
                                no_compacted_values_[n_p++] = plain_values_[m * sub_size + c]; // copy values
                            }
                            in_voc_[m] = 0;
                        }
                    }

                    // Encode structures
                    k_is_in_voc = bit_vector_type(in_voc_);
                    sdsl::util::init_support(this->k_is_in_voc_rank1, &this->k_is_in_voc);
                    k_encoded_values = t_values_last(encoded_values_);
                    this->k_plain_values = t_values_last(no_compacted_values_);
                } // END SIXTH STEP
            } else {
                k_size_word = 0;
                k_n_words = 0;
            } // END IF plain_values > 0
        }
    }; // END class k2_raster_heuristic

    // Types
    typedef k2raster::k2_raster_heuristic<> k2rh_type;

} // END namespace sdsl


#endif // INCLUDED_K2_RASTER_HEURISTIC
