/*  
 * Created by Fernando Silva on 22/09/17.
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

#include <algebra/algebra.h>

namespace k2raster_algebra_static {
//**********************************************************************//
//************************** ALGEBRA ***********************************//
//**********************************************************************//
    K2Raster* algebra(K2Raster *k2raster1, K2Raster *k2raster2, OperationRaster operation,
                      uint k1, uint k2, uint levelK1, bool deleteInputs) {
        // Check if there is a operation with that code
        printf("Operation: ");
        switch (operation) {
            case OperationRaster::OPERATION_SUM:
                printf("'Sum'");
                break;
            case OperationRaster::OPERATION_SUBT:
                printf("'Subtraction'");
                break;
            case OperationRaster::OPERATION_MULT:
                printf("'Multiplication'");
                break;
            default:
                printf("No valid operation %u\n", operation);
                return nullptr;
        }
        printf(" with k1=%u, k2=%u and levelK1=%u\n", k1, k2, levelK1);

        if (k1 == 0 && k2 == 0 && levelK1 == 0) {
            printf("Equal K's\n");
            return new K2Raster(k2raster1, k2raster2, operation, deleteInputs);
        } else {
            return new K2Raster(k1, k2, levelK1, k2raster1, k2raster2, operation, deleteInputs);
        }
    }
} // END namespace k2raster_algebra_static