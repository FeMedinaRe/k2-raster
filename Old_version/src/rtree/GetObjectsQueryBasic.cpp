/*  
 * Created by Fernando Silva on 25/10/17.
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

#include <rtree/GetObjectsQueryBasic.h>

namespace k2raster_rtree_static {

    //**********************************************************************//
    //******************* SETTERS AND GETTERS ******************************//
    //**********************************************************************//
    void GetObjectsQueryBasic::setRaster(K2Raster *raster) {
        this->raster = raster;
        this->sizeX = raster->getDimSize(1);
        this->sizeY = raster->getDimSize(2);
    }

    void GetObjectsQueryBasic::initialize() {
        this->first = true;
    }

    void GetObjectsQueryBasic::clear() {
        this->first = true;
    }

    //**********************************************************************//
    //*************************** POSITIONS ********************************//
    //**********************************************************************//
    void GetObjectsQueryBasic::calculatePositionParams(double *m_pLow, double *m_pHigh) {
        this->offsetX = m_pLow[0] * -1;
        this->offsetY = m_pLow[1] * -1;

        this->treeSizeX = m_pHigh[0] + this->offsetX;
        this->treeSizeY = m_pHigh[1] + this->offsetY;
    }

    void GetObjectsQueryBasic::calculatePositions(double *m_pLow, double *m_pHigh,
                                                     uint &xini, uint &xend, uint &yini, uint &yend) {
        xini = ((m_pLow[0] + this->offsetX) / this->treeSizeX) * (this->sizeX - 1);
        xend = ((m_pHigh[0] + this->offsetX) / this->treeSizeX) * (this->sizeX - 1);
        yini = ((m_pLow[1] + this->offsetY) / this->treeSizeY) * (this->sizeY - 1);
        yend = ((m_pHigh[1] + this->offsetY) / this->treeSizeY) * (this->sizeY - 1);

    }


} // END NAMESPACE k2raster_rtree_static