/*  
 * Created by Fernando Silva on 17/01/17.
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

#ifndef K2_RASTER_GETOBJECTSVISITORCELLS_H_
#define K2_RASTER_GETOBJECTSVISITORCELLS_H_

// System libraries
#include <map>
#include <tuple>
#include <stdexcept>      // std::out_of_range
// Own libraries
#include <include/spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

namespace k2raster_rtree_static {
    class GetObjectsVisitorCells : public IVisitor {

    public:
        /****** Constructor ******/
        GetObjectsVisitorCells(double *m_pLow, double *m_pHigh, uint rasterSizeX, uint rasterSizeY);

        /****** Visit ******/
        void visitNode(const INode &n);

        void visitData(const IData &d) {};

        void visitData(std::vector<const IData *> &v) {};

        /****** Result ******/
        std::map<id_type, std::tuple<uint *, ulong> *> getResult() const;

        ulong getAllMBRCells() const;

        /****** Getters and Setters ******/
        void clear();

    private:
        /****** Info ******/
        uint point;

        // r-tree params to transform r-tree position into raster position
        double offsetX;
        double offsetY;
        double treeSizeX;
        double treeSizeY;
        uint rasterSizeX;
        uint rasterSizeY;

        /****** Results ******/
        std::map<id_type, std::tuple<uint *, ulong> *> result;
        ulong allMBRCells;

        /****** Positions ******/
        void calculatePositionParams(double *m_pLow, double *m_pHigh, uint rasterSizeX, uint rasterSizeY);

        void calculatePositions(double *m_pLow, double *m_pHigh,
                                uint &xini, uint &xend, uint &yini, uint &yend);
    };
}

#endif // K2_RASTER_GETOBJECTSVISITORCELLS_H_
