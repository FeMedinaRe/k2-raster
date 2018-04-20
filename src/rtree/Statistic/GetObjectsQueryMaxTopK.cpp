/*  
 * Created by Fernando Silva on 30/10/17.
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

#include <rtree/Statistic/QueryStrategy/GetObjectsQueryMaxTopK.h>

namespace k2raster_rtree_static {

    //**********************************************************************//
    //***************************** AUX ************************************//
    //**********************************************************************//

    int GetObjectsQueryMaxTopK::getMBRValue(BasicSubMatrix *basicSubMatrix) {
        return basicSubMatrix->maxValue;
    }

    int GetObjectsQueryMaxTopK::getMBRLeafValue(uint xini, uint xend, uint yini, uint yend,
                                                BasicSubMatrix *basicSubMatrix) {
        return this->raster->getMaxValueWindow(xini, xend, yini, yend, basicSubMatrix);
    }

    ulong GetObjectsQueryMaxTopK::getCells(priorityPointerMBR *ptr) {
        return this->raster->getCellsByValue(ptr->positions->xini, ptr->positions->xend, ptr->positions->yini, ptr->positions->yend,
                                      ptr->value, ptr->value, ptr->pointer, this->positions[this->numResults], false);
    }



} // END NAMESPACE k2raster_rtree_static