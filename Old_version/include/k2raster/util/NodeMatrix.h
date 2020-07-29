/*
 * Created by Fernando Silva on 22/04/16.
 *
 * Copyright (C) 2016-current-year, Fernando Silva, all rights reserved.
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

#ifndef UTIL_NODEMATRIX_H_
#define UTIL_NODEMATRIX_H_

#include <sys/types.h>
#include <cstdlib>
#include <climits>

namespace k2raster_static {

    class NodeMatrix {

    public:
        NodeMatrix() {
            this->minValue = INT_MAX;                    // MinValue of this submatrix
            this->maxValue = -1;                        // MaxValue of this submatrix
            this->children = nullptr;
            this->plainChildren = nullptr;
        }
        ~NodeMatrix() {
            if (this->children != nullptr) {
                free(this->children);
            }
            if (this->plainChildren != nullptr) {
                free(this->plainChildren);
            }
        };

        int maxValue;
        int minValue;
        NodeMatrix **children;
        int *plainChildren;

    };

    class NodeMatrixSum {

    public:
        NodeMatrixSum() { }

        int sumValue;
        NodeMatrixSum **children;
    };
}


#endif // UTIL_NODEMATRIX_H_
