/*  
 * Created by Fernando Silva on 16/01/19.
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

#ifndef INCLUDE_UTILS_UTILS_QUERY
#define INCLUDE_UTILS_UTILS_QUERY

#include <fstream>
#include <random>

//! Namespace for k2-raster library
namespace k2raster {

    //*******************************************************//
    //************* CREATE ACCESS QUERIES *******************//
    //*******************************************************//
    template <typename size_type=size_t, bool add_time=true, bool add_n_queries_beginning=true>
    void create_access_queries(size_type n_queries, std::ofstream &outFile, size_type xini, size_type xend, size_type yini, size_type yend,
                               size_type tini=0, size_type tend=0) {
        // Seed with a real random value, if available
        std::random_device r;
        std::mt19937 gen(r());

        std::uniform_int_distribution<size_type>    dis_x(xini, xend);
        std::uniform_int_distribution<size_type>    dis_y(yini, yend);
        std::uniform_int_distribution<size_type>    dis_t(tini, tend);

        if (add_n_queries_beginning) {
            outFile << n_queries << std::endl;
        }

        // Create random queries
        size_type x, y, t;
        for (size_type q = 0; q < n_queries; q++) {
            // Add x (rows)
            x = dis_x(gen);
            outFile << x << " ";

            // Add y (cols)
            y = dis_y(gen);
            outFile << y;

            // Add z (time)
            if (add_time) {
                t = dis_t(gen);
                outFile << " " << t;
            }

            if (q != (n_queries - 1)) {
                outFile << std::endl;
            }
        } // END FOR n_queries
    }



    //*******************************************************//
    //************* CREATE REGION QUERIES *******************//
    //*******************************************************//
    template <typename size_type=size_t, typename value_type=int, bool add_time=true, bool add_values=true, bool add_n_queries_beginning=true>
    void create_region_queries(size_type n_queries, std::ofstream &outFile, size_type xini, size_type xend, size_type yini, size_type yend,
            size_type tini=0, size_type tend=0, value_type valini=0, value_type valend=0) {

        // Seed with a real random value, if available
        std::random_device r;
        std::mt19937 gen(r());

        std::uniform_int_distribution<size_type>    dis_x(xini, xend);
        std::uniform_int_distribution<size_type>    dis_y(yini, yend);
        std::uniform_int_distribution<value_type>   dis_val(valini, valend);
        std::uniform_int_distribution<size_type>    dis_t(tini, tend);

        if (add_n_queries_beginning) {
            outFile << n_queries << std::endl;
        }

        size_type x1, x2, y1, y2, t1, t2;
        value_type v1, v2;
        for (size_type q = 0; q < n_queries; q++) {
            // Add x (rows)
            x1 = dis_x(gen);
            x2 = dis_x(gen);
            outFile << std::min(x1, x2) << " " << std::max(x1,x2) << " ";

            // Add y (cols)
            y1 = dis_y(gen);
            y2 = dis_y(gen);
            outFile << std::min(y1, y2) << " " << std::max(y1, y2);


            // Add z (time)
            if (add_time) {
                t1 = dis_t(gen);
                t2 = dis_t(gen);
                outFile << " " << std::min(t1, t2) << " " << std::max(t1, t2);
            }

            // Add values
            if (add_values) {
                v1 = dis_val(gen);
                v2 = dis_val(gen);
                outFile << " " << std::min(v1, v2) << " " << std::max(v1, v2);
            }

            if (q != (n_queries -1)) {
                outFile << std::endl;
            }
        } // END FOR n_queries
    }

    //*******************************************************//
    //*********** CREATE FIXED REGION QUERIES ***************//
    //*******************************************************//
    template <typename size_type=size_t, typename value_type=int, bool add_time=true, bool add_values=true, bool add_n_queries_beginning=true>
    void create_fixed_region_queries(size_type n_queries, std::ofstream &outFile,
                size_type xini, size_type xend, size_type row_size,
                size_type yini, size_type yend, size_type col_size,
                size_type tini=0, size_type tend=0, size_type time_size=0,
                value_type valini=0, value_type valend=0, size_type value_size=0) {

        // Seed with a real random value, if available
        std::random_device r;
        std::mt19937 gen(r());



        // Check max size of the window
        row_size = ((xend - xini + 1) < row_size) ? (xend - xini) : row_size - 1;
        col_size = ((yend - yini + 1) < col_size) ? (yend - yini) : col_size - 1;
        time_size = ((tend - tini + 1) < time_size) ? (tend - tini) : time_size - 1;
        value_size = ((valend - valini + 1) < (value_type)value_size) ? (valend - valini) : (value_type)(value_size - 1);

        std::uniform_int_distribution<size_type>    dis_x(xini, xend - row_size);
        std::uniform_int_distribution<size_type>    dis_y(yini, yend - col_size);
        std::uniform_int_distribution<value_type>   dis_val(valini, valend - value_size);
        std::uniform_int_distribution<size_type>    dis_t(tini, tend - time_size);

        if (add_n_queries_beginning) {
            outFile << n_queries << std::endl;
        }

        size_type x1, y1, t1;
        value_type v1;
        for (size_type q = 0; q < n_queries; q++) {
            // Add x (rows)
            x1 = dis_x(gen);
            outFile << x1 << " " << x1 + row_size << " ";

            // Add y (cols)
            y1 = dis_y(gen);
            outFile << y1 << " " << y1 + col_size;

            // Add z (time)
            if (add_time) {
                t1 = dis_t(gen);
                outFile << " " << t1 << " " << t1 + time_size;
            }

            // Add values
            if (add_values) {
                v1 = dis_val(gen);
                outFile << " " << v1 << " " << v1 + value_size;
            }

            if (q != (n_queries -1)) {
                outFile << std::endl;
            }
        } // END FOR n_queries
    }
} // END NAMESPACE k2raster

#endif // INCLUDE_UTILS_UTILS_QUERY
