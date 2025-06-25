/*
 * Created by Fernando Silva on 15/02/22.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
 *
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * Library to set and check command line arguments
 *
 * k2-raster is a compact data structure to represent raster data that
 * uses compressed space and offers indexing capabilities.
 * It uses min/max values for indexing and improving query performance.
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

#ifndef INCLUDE_UTILS_UTILS_MEM
#define INCLUDE_UTILS_UTILS_MEM

#include <malloc_count/malloc_count.h>

namespace util {

void print_memory_consumption(ushort step=0) {
    if (step > 0) {
        std::cout << "Step: " << step <<": ";
    }
    std::cout << "[[Mem. consumption]]  Current allocation: ";
    std::cout << malloc_count_current()/1024/1024 << " MB, Peak allocation: ";
    std::cout << malloc_count_peak()/1024/1024 << " MB" << std::endl;

    //if (malloc_count_peak()/1024/1024 > 163) exit(-1);
}

} // END namespace util
#endif // INCLUDE_UTILS_UTILS_MEM