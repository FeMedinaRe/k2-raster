/*  
 * Created by Fernando Silva on 5/10/17.
 *
 * Copyright (C) 2017-current-year, Fernando Silva, all rights reserved.
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


#include <k2raster/util/InfoMinSubMatrix.h>

namespace k2raster_static {
//*******************************************************//
//************** FUNCTIONS ******************************//
//*******************************************************//

    void removeChildren(BasicSubMatrixWithChildren *matrix) {
        if (matrix->children != nullptr) {
            for (uint c = 0; c < matrix->numOfChildren; c++) {
                removeChildren(matrix->children[c]);
                delete(matrix->children[c]);
            }
            delete[](matrix->children);
            matrix->children = nullptr;
            matrix->numOfChildren = 0;
        }
        return;
    }
} // END NAMESPACE k2raster_static