/*  
 * Created by Fernando Silva on 22/02/19.
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

// External libraries
#include <gtest/gtest.h>
#include <utils/dac_vector.hpp>
#include <sdsl/dac_vector.hpp>

using namespace sdsl;
using testing::Types;


// Global Params
int min_value = 0;
int max_value = 1000000;
size_t max_size = 100000000;
size_t n_queries = 1000000;

// Random number generator
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<uint> dis_values(min_value, max_value);
std::uniform_int_distribution<uint> dis_size(0, max_size-1);


// TEST CLASS
template<class T>
class test_new_DAC : public ::testing::Test { };

typedef Types<
        dac_vector_dp_opt<>
> Implementations;

TYPED_TEST_CASE(test_new_DAC, Implementations);

TYPED_TEST(test_new_DAC, build){

    std::vector<uint> values(max_size);

    for (auto p = 0; p < max_size; p++) {
        values[p] = dis_values(gen);
    }

    // Run New DAC
    TypeParam dac(values);
    for (auto p = 0; p < max_size; p++) {
        ASSERT_TRUE(values[p] == dac[p]);
    }

    // Run old DAC
    dac_vector_dp<> dac_old(values);

    size_t total_size = max_size * sizeof(uint);
    std::cout << "DAC of size " << max_size << std::endl;
    std::cout << "Space:" << std::endl;
    std::cout << "New: " << size_in_mega_bytes(dac) << "Mbs " << "(" <<  (size_in_bytes(dac) * 100.) / total_size << "%)";
    std::cout << "vs Old: " << size_in_mega_bytes(dac_old) << "Mbs " << "(" <<  (size_in_bytes(dac_old) * 100.) / total_size << "%) ";
    std::cout << "Diff: " << size_in_mega_bytes(dac_old) - size_in_mega_bytes(dac) << "Mbs" << std::endl;


    { // Time decompress
        uint val_old=0, val1=0, val2=0, val3=0, val4=0, val5=0, val6=0;


        auto start_time_dac_old = std::chrono::high_resolution_clock::now();
        for (auto p = 0; p < max_size; p++) {
            val_old += dac_old[p];
        }
        auto stop_time_dac_old = std::chrono::high_resolution_clock::now();


        auto start_time_dac = std::chrono::high_resolution_clock::now();
        for (auto p = 0; p < max_size; p++) {
            val1 += dac[p];
        }
        auto stop_time_dac = std::chrono::high_resolution_clock::now();

        auto start_time_dac_iter = std::chrono::high_resolution_clock::now();
        for (auto val : dac) {
            val2 += val;
        }
        auto stop_time_dac_iter = std::chrono::high_resolution_clock::now();


        auto start_time_dac_for = std::chrono::high_resolution_clock::now();
        for (auto p = 0; p < max_size; p++) {
            val3 += dac.accessFor(p);
        }
        auto stop_time_dac_for = std::chrono::high_resolution_clock::now();

        val4 = dac.accessFor(0);
        auto start_time_dac_next = std::chrono::high_resolution_clock::now();
        for (auto p = 1; p < max_size; p++) {
            val4 += dac.next();
        }
        auto stop_time_dac_next = std::chrono::high_resolution_clock::now();

        auto start_time_dac_access = std::chrono::high_resolution_clock::now();
        for (auto p = 0; p < max_size; p++) {
            val5 += dac.access(p);
        }
        auto stop_time_dac_access = std::chrono::high_resolution_clock::now();

        std::cout << "Time:" << std::endl;
        std::cout << "Val: " << val_old << " | " << val1 << " | " << val2 << " | " << val3 << " | " << val4 << " | " << val5 << " | " << val6 << std::endl;

        std::cout << "Dac (old): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                stop_time_dac_old - start_time_dac_old).count();
        std::cout << " milliseconds." << std::endl;

        std::cout << "Dac: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(stop_time_dac - start_time_dac).count();
        std::cout << " milliseconds." << std::endl;

        std::cout << "Dac (iter): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(stop_time_dac_iter - start_time_dac_iter).count();
        std::cout << " milliseconds." << std::endl;

        std::cout << "Dac (For): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(stop_time_dac_for - start_time_dac_for).count();
        std::cout << " milliseconds." << std::endl;

        std::cout << "Dac (next): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                stop_time_dac_next - start_time_dac_next).count();
        std::cout << " milliseconds." << std::endl;

        std::cout << "Dac (access): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                stop_time_dac_access - start_time_dac_access).count();
        std::cout << " milliseconds." << std::endl;
    }

    { // Time random
        uint val1=0, val2=0, val3=0, val4=0;
        std::vector<uint> queries(n_queries);
        for (auto p = 0; p < n_queries; p++) {
            queries[p] = dis_size(gen);
        }

        auto start_time_dac = std::chrono::high_resolution_clock::now();
        for (auto p = 0; p < n_queries; p++) {
            val1 += dac[queries[p]];
        }
        auto stop_time_dac = std::chrono::high_resolution_clock::now();

        auto start_time_dac_for = std::chrono::high_resolution_clock::now();
        for (auto p = 0; p < n_queries; p++) {
            val3 += dac.accessFor(queries[p]);
        }
        auto stop_time_dac_for = std::chrono::high_resolution_clock::now();

        auto start_time_dac_old = std::chrono::high_resolution_clock::now();
        for (auto p = 0; p < n_queries; p++) {
            val2 += dac_old[queries[p]];
        }
        auto stop_time_dac_old = std::chrono::high_resolution_clock::now();

        auto start_time_dac_access = std::chrono::high_resolution_clock::now();
        for (auto p = 0; p < n_queries; p++) {
            val4 += dac.access(queries[p]);
        }
        auto stop_time_dac_access = std::chrono::high_resolution_clock::now();

        std::cout << "Time (random):" << std::endl;
        std::cout << "Val: " << val1 << " | " << val3 << " | " << val2 << " | " << val4 << std::endl;
        std::cout << "Dac: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(stop_time_dac - start_time_dac).count();
        std::cout << " milliseconds." << std::endl;
        std::cout << "Dac (For): "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(stop_time_dac_for - start_time_dac_for).count();
        std::cout << " milliseconds." << std::endl;
        std::cout << "Dac: " << std::chrono::duration_cast<std::chrono::milliseconds>(
                stop_time_dac_old - start_time_dac_old).count();
        std::cout << " milliseconds." << std::endl;
        std::cout << "Dac (access): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                stop_time_dac_access - start_time_dac_access).count();
        std::cout << " milliseconds." << std::endl;
    }





}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
