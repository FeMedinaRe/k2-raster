/*  
 * Created by Fernando Silva on 9/01/19.
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

#ifndef INCLUDED_K2_RASTER_TEMPORAL_BASE
#define INCLUDED_K2_RASTER_TEMPORAL_BASE

#include <sdsl/int_vector.hpp>
#include <utils/utils_data.hpp>

//! Namespace for k2-raster library
namespace k2raster {

    // Types
    static const ushort K2_RASTER_TEMPORAL_TYPE = 100;
    static const ushort K2_RASTER_TEMPORAL_TYPE_H = 101;
    static const ushort T_K2_RASTER_TYPE = 110;
    static const ushort AT_K2_RASTER_TYPE = 111;
    static const ushort ATH_K2_RASTER_TYPE = 112;
    static const ushort K2_RASTER_TEMPORAL_GLOBAL_TYPE = 120;


template <typename t_value=int>
class k2_raster_temporal_base
{
public:
    typedef t_value                         value_type;
    typedef sdsl::int_vector<>::size_type   size_type;

protected:
    /****** k2-raster temporal Type ******/
    ushort                      k_k2_raster_type; // k2-raster Type

    /****** Size params ******/
    size_type                   k_max_size_x;  // Max row size of all rasters
    size_type                   k_max_size_y;  // Max column size of all rasters

    /****** Values ******/
    value_type                  k_max_value;    // Max value of all rasters
    value_type                  k_min_value;    // Min value of all rasters

    /****** Size params ******/
    size_type                   k_n_times;      // Number of times (rasters)

public:
    const size_type& m_max_size_x =     k_max_size_x;
    const size_type& m_max_size_y =     k_max_size_y;

    const value_type& m_max_value =     k_max_value;
    const value_type& m_min_value =     k_min_value;
    const size_type& m_n_time =         k_n_times;

public:
    //*******************************************************//
    //******************* CONSTRUCTOR ***********************//
    //*******************************************************//
    k2_raster_temporal_base() {};

    k2_raster_temporal_base(const k2_raster_temporal_base& tr)
    {
        *this = tr;
    }

    k2_raster_temporal_base(k2_raster_temporal_base&& tr)
    {
        *this = std::move(tr);
    }

    //*******************************************************//
    //*************** BASIC OPERATIONS **********************//
    //*******************************************************//

    //! Move assignment operator
    k2_raster_temporal_base& operator=(k2_raster_temporal_base&& tr)
    {
        if (this != &tr) {
            /****** Size params ******/
            k_k2_raster_type = std::move(tr.k_k2_raster_type);
            k_max_size_x = std::move(tr.k_max_size_x);
            k_max_size_y = std::move(tr.k_max_size_y);
            k_max_value = std::move(tr.k_max_value);
            k_min_value = std::move(tr.k_min_value);
            k_n_times = std::move(tr.k_n_times);

        }
        return *this;
    }

    //! Assignment operator
    k2_raster_temporal_base& operator=(const k2_raster_temporal_base& tr)
    {
        if (this != &tr) {
            /****** Size params ******/
            k_k2_raster_type = tr.k_k2_raster_type;
            k_max_size_x = tr.k_max_size_x;
            k_max_size_y = tr.k_max_size_y;
            k_max_value = tr.k_max_value;
            k_min_value = tr.k_min_value;
            k_n_times = tr.k_n_times;
        }
        return *this;
    }

    //! Swap operator
    void swap(k2_raster_temporal_base& tr)
    {
        if (this != &tr) {
            /****** Size params ******/
            std::swap(k_k2_raster_type, tr.k_k2_raster_type);
            std::swap(k_max_size_x, tr.k_max_size_x);
            std::swap(k_max_size_y, tr.k_max_size_y);
            std::swap(k_max_value, tr.k_max_value);
            std::swap(k_min_value, tr.k_min_value);
            std::swap(k_n_times, tr.k_n_times);
        }
    }

    //! Equal operator
    bool operator==(const k2_raster_temporal_base& tr) const
    {
        if (k_k2_raster_type != tr.k_k2_raster_type) {
            return false;
        }
        if (k_max_size_x != tr.k_max_size_x && k_max_size_y != tr.k_max_size_y) {
            return false;
        }

        if (k_min_value != tr.k_min_value && k_max_value != tr.k_max_value) {
            return false;
        }

        if (k_n_times != tr.k_n_times) {
            return false;
        }
        return true;
    }

    //*******************************************************//
    //********************* GETTERS *************************//
    //*******************************************************//
    virtual inline value_type get_min_value() const { return k_min_value;}
    virtual inline value_type get_max_value() const { return k_max_value;}
    virtual inline value_type get_n_rows() const { return k_max_size_x;}
    virtual inline value_type get_n_cols() const { return k_max_size_y;}
    virtual inline value_type get_max_t() const { return k_n_times;}

    //*******************************************************//
    //******************** QUERIES **************************//
    //*******************************************************//
    virtual value_type get_cell(size_type row, size_type col, size_type time) const=0;

    virtual size_type get_cells_by_value(size_type xini, size_type xend, size_type yini, size_type yend,
                                         value_type valini, value_type valend, size_type tmin, size_type tmax,
                                         std::vector<std::vector<std::pair<size_type, size_type>>> &result)=0;

    virtual size_type get_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                        size_type tmin, size_type tmax,
                                        std::vector<std::vector<value_type>> &result)=0;

    virtual bool check_values_window(size_type xini, size_type xend, size_type yini, size_type yend,
                                     value_type valini, value_type valend,
                                     size_type tmin, size_type tmax, bool strong_check)=0;

    //*******************************************************//
    //********************** FILE ***************************//
    //*******************************************************//
    virtual size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const {

        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        /****** k2-raster Type ******/
        written_bytes += write_member(k_k2_raster_type, out, child, "k2-raster_type");

        /****** Size params ******/
        written_bytes += write_member(k_max_size_x, out, child, "max_size_x");
        written_bytes += write_member(k_max_size_y, out, child, "max_size_y");

        /****** Values ******/
        written_bytes += write_member(k_max_value, out, child, "max_value");
        written_bytes += write_member(k_min_value, out, child, "min_value");

        /****** Size params ******/
        written_bytes += write_member(k_n_times, out, child, "n_times");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    virtual void load(std::istream& in) {
        /****** k2-raster Type ******/
        sdsl::read_member(k_k2_raster_type, in);

        /****** Size params ******/
        sdsl::read_member(k_max_size_x, in);
        sdsl::read_member(k_max_size_y, in);

        /****** Values ******/
        sdsl::read_member(k_max_value, in);
        sdsl::read_member(k_min_value, in);

        /****** Size params ******/
        sdsl::read_member(k_n_times, in);
    }

    //*******************************************************//
    //********************** TEST ***************************//
    //*******************************************************//
    bool check(std::string &inputs_filename, std::string &input_path_folder, ushort scale_factor=0) {

        /**************************/
        /* Reads params from file */
        /**************************/
        std::ifstream inputs_file(inputs_filename);
        assert(inputs_file.is_open() && inputs_file.good());
        size_type time_count = 0;

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

                std::vector<value_type> values;
                read_input_data(file_path_input, n_rows, n_cols, values, scale_factor);

                /*********************/
                /* Check cell        */
                /*********************/
                value_type v1, v2;
                for (auto r = 0 ; r < n_rows; r++) {
                    for (auto c = 0; c < n_cols; c++) {
                        size_type pos = r * n_cols + c;
                        v1 = get_cell(r, c, time_count);
                        v2 = values[pos];
                        if ( v1 != v2) {
                            std::cout << "Error position (" << r << ", " << c << ") at time" << time_count;
                            std::cout << ": get " << v1 << " and expected " << v2 << std::endl;
                            return false;
                        }
                    } // END FOR c
                } // END FOR r
                time_count++;
            } // END WHiLE read_params
            inputs_file.close();
        } // END BLOCK encode k2-rasters
        return true;
    }

    virtual void print_space_by_time() const=0;

}; // END CLASS k2_raster_temporal_base
} // END NAMESPACE k2raster

#endif // INCLUDED_K2_RASTER_TEMPORAL_BASE
