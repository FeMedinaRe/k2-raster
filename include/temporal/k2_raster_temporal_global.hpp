/*  
 * Created by Fernando Silva on 29/10/18.
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

#ifndef INCLUDED_K2_RASTER_TEMPORAL_GLOBAL
#define INCLUDED_K2_RASTER_TEMPORAL_GLOBAL

// Own libraries
#include <k2_raster_base.hpp>
#include <utils/utils_heuristic.hpp>
#include <temporal/helpers/k2_raster_heuristic_log.hpp>
#include <temporal/k2_raster_temporal_base.hpp>

//! Namespace for k2-raster library
namespace k2raster {

    template<typename t_value=int,
            typename t_bv=sdsl::bit_vector,
            typename t_rank_1=sdsl::rank_support_v5<1,1>,
            typename t_values_vec=sdsl::dac_vector_dp_opt<t_bv, t_rank_1, 3>,
            typename t_values_last=sdsl::dac_vector_dp_opt<t_bv, t_rank_1, 3>>
    class k2_raster_temporal_global : public k2_raster_temporal_base<t_value> {

    public:
        typedef t_value                                 value_type;
        typedef k2_raster_temporal_base<>::size_type    size_type;
        typedef uint                                    codeword_type;
        typedef k2_raster_plain<t_value, t_bv, t_rank_1, t_values_vec, std::vector<t_value>> k2_raster_t;

    protected:
        std::vector<k2_raster_heuristic_log<t_value, t_bv, t_rank_1, t_values_vec, t_values_last>> k_rasters;

        // Vocabulary
        size_type                   k_size_word;        // Size (in bytes) of one word of the vocabulary
        size_type                   k_n_words;          // Number of words in the vocabulary
        std::vector<char>           k_voc;              // Array that stores vocabulary words

    public:
        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//
        k2_raster_temporal_global() = default;

        k2_raster_temporal_global(const k2_raster_temporal_global &tr) {
            *this = tr;
        }

        k2_raster_temporal_global(k2_raster_temporal_global &&tr) {
            *this = std::move(tr);
        }

        // Only interface (snapshots_freq unused)
        k2_raster_temporal_global(std::string &inputs_filename, std::string &input_path_folder, const size_type snapshots_freq __attribute__((unused)),
                const ushort scale_factor=0 ){
            this->k_k2_raster_type = K2_RASTER_TEMPORAL_GLOBAL_TYPE;
            this->k_n_times = build(inputs_filename, input_path_folder, scale_factor);
        }


        //*******************************************************//
        //*************** BASIC OPERATIONS **********************//
        //*******************************************************//

        //! Move assignment operator
        k2_raster_temporal_global &operator=(k2_raster_temporal_global &&tr) {
            if (this != &tr) {
                k2_raster_temporal_base<t_value>::operator=(tr);
                k_rasters = std::move(tr.k_rasters);
                k_size_word = std::move(tr.k_size_word);
                k_n_words = std::move(tr.k_n_words);
                k_voc = std::move(tr.k_voc);
            }
            return *this;
        }

        //! Assignment operator
        k2_raster_temporal_global &operator=(const k2_raster_temporal_global &tr) {
            if (this != &tr) {
                k2_raster_temporal_base<t_value>::operator=(tr);
                k_rasters = tr.k_rasters;
                k_size_word = tr.k_size_word;
                k_n_words = tr.k_n_words;
                k_voc = tr.k_voc;
            }
            return *this;
        }

        //! Swap operator
        void swap(k2_raster_temporal_global &tr) {
            if (this != &tr) {
                k2_raster_temporal_base<t_value>::swap(tr);
                std::swap(k_rasters, tr.k_rasters);
                std::swap(k_size_word, tr.k_size_word);
                std::swap(k_n_words, tr.k_n_words);
                std::swap(k_voc, tr.k_voc);
            }
        }

        //! Equal operator
        bool operator==(const k2_raster_temporal_global &tr) const {
            if (!k2_raster_temporal_base<t_value>::operator==(tr)) {
                return false;
            }
            if (k_rasters.size() != tr.k_rasters.size()) {
                return false;
            }
            size_type c = 0;
            for (auto raster : k_rasters) {
                if (raster != tr.k_rasters[c++]) {
                    return false;
                }
            }
            return true;
        }


        //*******************************************************//
        //******************** QUERIES **************************//
        //*******************************************************//
        value_type get_cell(size_type row, size_type col, size_type time) const {
            assert(time <= (k_rasters.size()-1));
            return k_rasters[time].get_cell(row, col);
        }

        size_type get_cells_by_value(size_type xini, size_type xend, size_type yini, size_type yend,
                                     value_type valini, value_type valend, size_type tmin, size_type tmax,
                                     std::vector<std::vector<std::pair<size_type, size_type>>> &result) {

            result.resize(tmax - tmin + 1);
            size_type c = 0;
            size_type count_cells = 0;
            for (auto t = tmin; t <= tmax; t++) {
                count_cells += k_rasters[t].get_cells_by_value(xini, xend, yini, yend, valini, valend, result[c++]);
            }
            return count_cells;
        }

        size_type get_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                    size_type tmin, size_type tmax,
                                    std::vector<std::vector<value_type>> &result) {
            result.resize(tmax - tmin + 1);
            size_type c = 0;
            size_type count_cells = 0;
            for (auto t = tmin; t <= tmax; t++) {
                count_cells += k_rasters[t].get_values_window(xini, xend, yini, yend, result[c++]);
            }
            return count_cells;
        }

        bool check_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                 value_type valini, value_type valend,
                                 size_type tmin, size_type tmax, bool strong_check) {
            bool result;
            for (auto t = tmin; t <= tmax; t++) {
                result = k_rasters[t].check_values_window(xini, xend, yini, yend, valini, valend, strong_check);
                if (strong_check){
                    if (!result) return false;
                } else {
                    if (result) return true;
                }
            }
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

            /****** Vocabulary ******/
            written_bytes += write_member(k_size_word, out, child, "size_word");
            written_bytes += write_member(k_n_words, out, child, "n_words");
            written_bytes += serialize_vector(k_voc, out, child, "voc");


            /****** Vector of rasters ******/
            written_bytes += write_member(k_rasters.size(), out, child, "num_t");
            for (size_t l = 0; l < k_rasters.size(); l++) {
                std::string name_r = "raster_t_" + std::to_string(l);
                written_bytes += k_rasters[l].serialize(out, child, name_r);
            }

            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        virtual void load(std::istream &in) {
            /****** k2-raster base ******/
            k2_raster_temporal_base<t_value>::load(in);

            //****** Vocabulary ******/
            sdsl::read_member(k_size_word, in);
            sdsl::read_member(k_n_words, in);
            k_voc.resize(k_n_words * k_size_word);
            sdsl::load_vector(k_voc, in);

            /****** Vector of rasters ******/
            ulong t;
            sdsl::read_member(t, in);
            k_rasters.resize(t);
            k_rasters.shrink_to_fit();
            for (ulong l = 0; l < t; l++) {
                k_rasters[l].load(in, k_voc);
            }
        }

        void print_space_by_time() const {
            size_type times = k_rasters.size();
            double size;
            double total_size = 0;
            for (size_type t = 0; t < times; t++) {
                size = sdsl::size_in_mega_bytes(k_rasters[t]);
                std::cout << "Time " << t << ": " << std::setprecision(2) << std::fixed << size << " Mbs" << std::endl;
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

            std::vector<k2_raster_t> k2_rasters_plain;
            std::unordered_map<std::string, size_type> hash_words;   // Stores the frequency of each word

            this->k_max_size_x = 0;
            this->k_max_size_y = 0;
            this->k_min_value = std::numeric_limits<value_type>::max();
            this->k_max_value = std::numeric_limits<value_type>::min();

            /**************************/
            /* Encode each k2-raster  */
            /**************************/
            {
                size_t n_rows, n_cols, k1, k2, level_k1, plain_levels;
                std::string file_path_input;

                /*********************/
                /* Read parameters   */
                /*********************/
                read_params_from_file(inputs_file, n_rows, n_cols, k1, k2, level_k1, plain_levels);

                /*********************/
                /* Process rasters   */
                /*********************/
                while (next_file_path(inputs_file, input_path_folder, file_path_input)) {

                    /*********************/
                    /* Encode data       */
                    /*********************/
                    std::vector<value_type> values;
                    read_input_data(file_path_input, n_rows, n_cols, values, scale_factor);

                    k2_rasters_plain.emplace_back(values, n_rows, n_cols, k1, k2, level_k1, plain_levels);

                    this->k_max_size_x = std::max(this->k_max_size_x, k2_rasters_plain[time_count].get_n_rows());
                    this->k_max_size_y = std::max(this->k_max_size_y, k2_rasters_plain[time_count].get_n_cols());
                    this->k_min_value = std::min(this->k_min_value, k2_rasters_plain[time_count].min_value);
                    this->k_max_value = std::max(this->k_max_value, k2_rasters_plain[time_count].max_value);
                    time_count++;
                } // END WHiLE read_params
                inputs_file.close();
            } // END BLOCK encode k2-rasters

            /**************************/
            /* Create dictionary      */
            /**************************/
            {
                create_dictionary(k2_rasters_plain, hash_words);
            } // END BLOCK dictionary

            /**************************/
            /* Encode values          */
            /**************************/
            {
                for(auto const& k2_raster: k2_rasters_plain) {
                    k_rasters.emplace_back(k2_raster, hash_words, k_size_word, k_n_words, k_voc);
                }
            }

            return time_count;
        }


        //*******************************************************//
        //******************* DICTIONARY ************************//
        //*******************************************************//

        void create_dictionary(std::vector<k2_raster_t> &k2_rasters_plain, std::unordered_map<std::string, size_type> &hash_words) {
            /******************************************************************************/
            /* ZERO STEP - Set some params                                                */
            /******************************************************************************/
            // This steps translate numbers into char values.
            // This is because we need to manage values ​​as words (strings), making it easier to use the hash table.
            // A word consists of the fusion of all values ​​of the submatrix
            auto sub_size = k2_rasters_plain[0].m_size_leaves;         // Size of each (leaf) submatrix
            k_size_word = sub_size * sizeof(value_type);            // Size, in bytes, of each submatrix
            size_type total_submatrices=0;

            /******************************************************************************/
            /* FIRST STEP - Calculate frequency of each "word" and each different "value" */
            /******************************************************************************/
            std::unordered_map<value_type , size_type> hash_values;  // Stores the frequency of each value
            {
                for(auto const& k2_raster: k2_rasters_plain) {
                    auto n_submatrices = k2_raster.m_plain_values.size() / sub_size;
                    total_submatrices+=n_submatrices;

                    k2raster_heuristic::create_hashes<value_type, size_type>(hash_words, hash_values, n_submatrices,
                                                                             sub_size,
                                                                             k_size_word, k2_raster.m_plain_values.data());
                }

#ifndef NDEBUG
                std::cout << "Found " << hash_words.size() << " words in leaves values and " << hash_values.size() << " values";
                    std::cout << " in " << total_submatrices << " submatrices" << std::endl;
#endif
            } // END FIRST STEP


            /******************************************************************************/
            /* SECOND STEP - Calculate entropy (vocabulary and values)                    */
            /******************************************************************************/
            double entropy_values = 0, entropy_words = 0;
            {
                k2raster_heuristic::calculate_entropy<value_type, size_type>(hash_words, total_submatrices,
                                                                             hash_values, total_submatrices * sub_size,
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
        }
    }; // END CLASS k2_raster_temporal_global

} // END NAMESPACE k2raster

#endif // INCLUDED_K2_RASTER_HEURISTIC_GLOBAL
