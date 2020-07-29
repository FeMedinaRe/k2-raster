/*  
 * Created by Fernando Silva on 31/10/18.
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

#ifndef SDSL_K2_RASTER_UTILS_HEURISTIC
#define SDSL_K2_RASTER_UTILS_HEURISTIC

// System libraries
#include <cstdio>
#include <string>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <algorithm>

// 3th libraries
#include <sdsl/io.hpp>

//! Namespace for the succint data structure library
namespace k2raster_heuristic {

    template <typename value_type = int, typename size_type = size_t>
    void create_hashes(std::unordered_map<std::string, size_type> &hash_words,
                       std::unordered_map<value_type , size_type> &hash_values,
                       size_type n_submatrices, size_type sub_size, size_type size_word,
                       const void* plain_values) {

        typename std::unordered_map<std::string, size_type>::iterator it_word;
        typename std::unordered_map<value_type, size_type>::iterator it_value;

        char *t_word;

        // Adds submatrices and their values into the hash tables
        for (auto m = 0; m < n_submatrices; m++) {
            // Get word
            t_word = &(((char *)plain_values)[m * size_word]); // TODO check this
            std::string word(t_word, size_word);

            // Insert word
            it_word = hash_words.find(word);
            if (it_word == hash_words.end()) {
                // New word
                hash_words[word] = 1;
            } else {
                // Increase word counter
                it_word->second++;
            }

            // Check all values of the submatrix
            value_type *values = (value_type*)(t_word);
            for (auto v = 0; v < sub_size; v++) {
                // Get value
                value_type value = values[v];

                // Insert value
                it_value = hash_values.find(value);
                if (it_value == hash_values.end()) {
                    // New word
                    hash_values[value] = 1;
                } else {
                    // Increase value counter
                    it_value->second++;
                }
            } // END FOR add values to hash
        } // END FOR add words to hash
    }

    template <typename value_type = int, typename size_type = size_t>
    void calculate_entropy(std::unordered_map<std::string, size_type> hash_words, size_type n_submatrices,
                           std::unordered_map<value_type , size_type> hash_values, size_type n_values,
                           double &entropy_words, double &entropy_values) {
        // Entropy of values
        double freq;
        for (auto it = hash_values.begin(); it != hash_values.end(); ++it ) {
            freq = it->second / (double)n_values;
            entropy_values += (freq * log2(freq));
        }
        entropy_values = -entropy_values;

        // Entropy of words
        for (auto it = hash_words.begin(); it != hash_words.end(); ++it ) {
            freq = it->second / (double)n_submatrices;
            entropy_words += (freq * log2(freq));
        }
        entropy_words = -entropy_words;
    }

    template <typename size_type = size_t>
    std::vector<std::string> sort_vocabulary(std::unordered_map<std::string, size_type> hash_words) {
        std::vector<std::string> sort_words(hash_words.size());
        // Copy words
        uint c = 0;
        for (std::pair<std::string, size_type> word : hash_words) {
            sort_words[c++] = word.first;
        }

        auto Temporal_comparison = [&hash_words] (std::string const & s1, std::string const & s2) -> bool
        {
            auto freq1 = hash_words[s1];
            auto freq2 = hash_words[s2];
            return freq1 > freq2;
        };
        std::sort(sort_words.begin(), sort_words.end(), Temporal_comparison);

        return sort_words;
    }



    template <typename value_type = int, typename size_type = size_t, typename codeword_type = size_t>
    codeword_type assign_codewords(std::vector<std::string> sort_words, std::unordered_map<std::string, size_type> &hash_words,
                               size_type sub_size, double entropy_words, double entropy_values,
                               size_type &n_compacted, size_type &n_no_compacted) {
        codeword_type codeword = 1; // 0 == no codeword
        size_type freq;
        size_type size_plain, size_dict;

        // For each word in the vocabulary
        for (auto word: sort_words) {
            freq = hash_words[word];

            // Calculates the size of the "freq" submatrices of "word" if its values are encoded in plain
            // freq             -> Number of submatrices "word"
            // sub_size         -> Size of each submatriz
            // entropy_values   -> Bits per value
            size_plain = (freq * sub_size) * entropy_values;

            // Calculates the size of the "freq" submatrices of "word" using a dictionary (1 codeword + "freq" entries)
            // (sub_size * sizeof(size_type) * 8) -> Size of the entry in the vocabulary
            // freq             -> Number of submatrices "word"
            // entropy_words    -> Bits per codeword
            size_dict = (sub_size * sizeof(value_type) * 8) + (freq * entropy_words);

            if (size_dict <= size_plain &&
                codeword != std::numeric_limits<codeword_type>::max()) {
                // Encode values as an entry in the dictionary
                n_compacted += freq;
                hash_words[word] = codeword;
                codeword++;
            } else {
                // Encode values as plain values
                hash_words[word] = 0;
                n_no_compacted += freq;
            }
        } // END FOR word : sort_words

        return codeword - 1;
    }

    template <typename size_type = size_t, typename codeword_type = size_t>
    void copy_words(std::vector<char> &voc, std::vector<std::string> sort_words, std::unordered_map<std::string, size_type> hash_words,
                    size_type n_words, size_type size_word) {

        voc.resize(n_words * size_word); // Initialize array k_words to store dictionary words
        codeword_type codeword;
        for (auto word : sort_words) {
            codeword = hash_words[word];
            if (codeword > 0) {
                // Copy word
                for (auto p = 0; p < size_word; p++) {
                    voc[((codeword-1) * size_word) + p] = word.data()[p];
                }
            } else {
                break; // No more valid words
            } // END FOR copy 'word' to 'voc'
        } // END FOR word : sort_words
    }

    template <typename value_type = int, typename size_type = size_t, typename codeword_type = size_t>
    void compact_submatrices() {
        // Temporal structures
//        int_vector<sizeof(codeword_type) * 8> encoded_values_(n_compacted);
//        int_vector<sizeof(value_type) * 8> no_compacted_values_(n_no_compacted * sub_size);
//        bit_vector in_voc_(n_submatrices, 0);
//
//        // Check each submatrix
//        codeword_type codeword;
//        size_type n_c = 0, n_p = 0;
//        for (auto m = 0; m < n_submatrices; m++) {
//            // Get word
//            std::string word(&(plain_words[m * k_size_word]), k_size_word);
//            codeword = hash_words[word];
//            if (codeword != 0) {
//                // Use dictionary
//                encoded_values_[n_c++] = codeword - 1;
//                in_voc_[m] = 1;
//            } else {
//                // Use plain form
//                for (auto c = 0; c < sub_size; c++) {
//                    no_compacted_values_[n_p++] = plain_values_[m * sub_size + c]; // copy values
//                }
//                in_voc_[m] = 0;
//            }
//        } // END FOR n:submatrices
    }

}; // END NAMESPACE sdsl_heuristic

#endif // SDSL_K2_RASTER_UTILS_HEURISTIC
