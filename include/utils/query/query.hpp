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

#ifndef INCLUDED_UTILS_QUERY
#define INCLUDED_UTILS_QUERY

#include <fstream>
#include <iostream>

//! Namespace for k2-raster library
namespace k2raster {

template <bool with_time=true, bool with_value=true,
        typename value_type=int, typename size_type=size_t>
class query {

public:
    size_type xini;
    size_type xend;
    size_type yini;
    size_type yend;

    value_type valini;
    value_type valend;

    size_type tini;
    size_type tend;

public:

    //*******************************************************//
    //******************* CONSTRUCTOR ***********************//
    //*******************************************************//
    query() = default;

    query(const query& tr) {
        *this = tr;
    }

    query(query&& tr) {
        *this = std::move(tr);
    }

    query(std::istream& query_file, bool access_query) {
        if (access_query) {
            read_from_file_access(query_file);
        } else {
            read_from_file_region(query_file);
        }
    }


    //*******************************************************//
    //**************** BASIC OPERATIONS *********************//
    //*******************************************************//

    //! Move assignment operator
    query& operator=(query&& tr) {
        if (this != &tr) {
            xini = std::move(tr.xini);
            xend = std::move(tr.xend);
            yini = std::move(tr.yini);
            yend = std::move(tr.yend);

            tini = std::move(tr.tini);
            tend = std::move(tr.tend);

            valini = std::move(tr.valini);
            valend = std::move(tr.valend);
        }
        return *this;
    }

    //! Assignment operator
    query& operator=(const query& tr) {
        if (this != &tr) {
            xini = tr.xini;
            xend = tr.xend;
            yini = tr.yini;
            yend = tr.yend;

            tini = tr.tini;
            tend = tr.tend;


            valini = tr.valini;
            valend = tr.valend;
        }
        return *this;
    }


    //! Swap operator
    void swap(query& tr) {
        if (this != &tr) {
            std::swap(xini, tr.xini);
            std::swap(xend, tr.xend);
            std::swap(yini, tr.yini);
            std::swap(yend, tr.yend);

            std::swap(tini, tr.tini);
            std::swap(tend, tr.tend);

            std::swap(valini, tr.valini);
            std::swap(valend, tr.valend);
        }
    }

    //*******************************************************//
    //******************** SET VALUES ***********************//
    //*******************************************************//
    bool read_from_file_region(std::istream& query_file) {

        if (!query_file.eof() && query_file.good()){
            query_file >> xini;
            query_file >> xend;
            query_file >> yini;
            query_file >> yend;

            if (with_time) {
                query_file >> tini;
                query_file >> tend;
             }

            if (with_value) {
                query_file >> valini;
                query_file >> valend;
            }
            return query_file.good();
        }
        return false;
    }

    bool read_from_file_access(std::istream& query_file) {

        if (!query_file.eof() && query_file.good()){
            query_file >> xini;
            query_file >> yini;
            xend = xini;
            yend = yini;

            if (with_time) {
                query_file >> tini;
                tend = tini;
            }

            return query_file.good();
        }
        return false;
    }

    //*******************************************************//
    //************************ INFO *************************//
    //*******************************************************//
    void print() const {
        std::cout << "Range  [" << xini << ", " << yini << " - " << xend << ", " << yend << "] ";
        if (with_value) {
            std::cout << "-> Values [" << valini << ", " << valend << "] ";
        }
        if (with_time) {
            std::cout << "-> Time [" << tini << ", " << tend << "]" ;
        }
        std::cout << std::endl;
    }

}; // END CLASS query
} // END NAMESPACE k2raster

#endif // INCLUDED_UTILS_QUERY
